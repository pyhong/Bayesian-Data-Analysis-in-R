
//  [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>
#include <limits.h>

using namespace RcppGSL;
//rtruncnorm
#define TAIL_LIMIT  2.983851594898812
#define EXPTLOVTL  0.00390733345803262
#define LOG_TAIL_LIMIT  1.09321494749176
#define SQ_TAIL_LIMIT  4.45168517019009
#define MAX_PTAIL  0.00390733345803262
#define basep   0.0078125
#define infinity (std::numeric_limits<double>::infinity() )
//#define NAN  (std::numeric_limits<double>::quiet_NaN())
//#define M_PI 3.141592653589793 

#define RTAIL_BIN  162
#define LTAIL_BIN  (-163)

class truncated
{
  
private:
  double *xhist;
  double *yhist;
  double *dxhist;
  int *yhistratio;
  int *whichbin;
  gsl_rng * rr;
  
  inline double rnd()
  {
    return (rand()/(RAND_MAX+1.0));
  }
  inline double itaillimit( double b)
  {
    return EXPTLOVTL - exp(SQ_TAIL_LIMIT - TAIL_LIMIT * b - LOG_TAIL_LIMIT );
  }
  
  inline double sample_tail( double a, double b)	
  {
    double sa =  1 ;
    double sb = exp(a*(a-b));
    if ( sa == sb )
      return (a*a<b*b?a:b);
    
    double u;
    do {
      u = rnd() * (sb-sa) + sa;
    } while( u == 0 );		
    u = a - log(u)/a;
    
    double v = rnd() * exp(0.5*(a-u)*(a-u));
    if ( v <= 1 )
      return u;
    else
      return NAN;
  }
  double truncate2(double a, double b)
  {
    double z;
    int y,ia,ib;
    double ah,bh;
    ah = a; bh = b;
    
    double ap = basep;
    double bp = basep;
    
    double p_tail_left = 0;
    if ( a <= -TAIL_LIMIT ) {
      ia = LTAIL_BIN;
      p_tail_left = MAX_PTAIL;
      ah = -TAIL_LIMIT;
      ap = 0;
    }
    else if ( a >= TAIL_LIMIT ) {
      do { 
        z = sample_tail( a, b );
      } while ( z != z );
      return z;
    }
    else {
      double ax = fabs(a);
      ia = whichbin[ (int)(ax*128.0) ];
      if (xhist[ia+1] <= ax )
        ++ia;
      ia = ((a>=0)?ia:-ia-1);
    }
    
    double p_tail_right = 0;
    if (b <= -TAIL_LIMIT ) {
      do {
        z = -sample_tail(-b,-a);
      } while (z != z);
      return z;
    }
    else if ( b >= TAIL_LIMIT ) {
      p_tail_right = MAX_PTAIL;
      ib = RTAIL_BIN;
      bh = TAIL_LIMIT;
      bp = 0;
    }
    else {
      double bx = fabs(b);
      ib = whichbin[ (int)(bx*128.0) ];
      if (xhist[ib+1] <= bx )
        ++ib;
      ib = ((b>=0)?ib:-ib-1);
    }
    
    if (ia == ib ) {
      do {
        z = rnd() * (bh-ah) + ah;
        y = rand();
      } while ( y >= yhistratio[ia] && y * yhist[ia] >= exp(-0.5*z*z)*(RAND_MAX+1.0) );
      return(z);
    }
    
    double middlep = (ib-ia-1)*basep;		
    do {			
      double u = rnd() * (middlep + ap + bp + p_tail_left + p_tail_right );
      if ( u < middlep ) {
        u /= basep;
        int ui = (int)u;
        u = modf(u,&z);
        int col = ui + ia  + 1;
        z = u * dxhist[col] + xhist[col];
        y = rand();
        
        if ( (y>=yhistratio[col] ) && y * yhist[col] >= exp(-0.5*z*z) * (RAND_MAX+1.0) )
          z = NAN;
        continue;
      }
      if (ap == basep )
        ap = (xhist[ia+1]-ah)*yhist[ia];
      if (u < middlep + ap ) {
        z = (u-middlep)/ap * (xhist[ia+1]-ah) + ah;
        y = rand();
        if ( y >= yhistratio[ia] && y * yhist[ia] >= exp(-0.5*z*z) * (RAND_MAX+1.0) )
          z = NAN;
        continue;
      }			
      if (bp == basep)
        bp = (bh-xhist[ib])*yhist[ib];
      if ( u < middlep + ap + bp) {
        z = (u-middlep-ap)/bp * (bh-xhist[ib])+xhist[ib];
        y = rand();
        if ( y >= yhistratio[ib] && y * yhist[ib] >= exp(-0.5*z*z) * (RAND_MAX+1.0) )
          z = NAN;
        continue;
      }			
      if ( p_tail_left == MAX_PTAIL )
        p_tail_left = itaillimit(-a);
      if ( u < middlep + ap + bp + p_tail_left ) {
        z = -sample_tail( TAIL_LIMIT, -a );
        continue;
      }
      if ( p_tail_right == MAX_PTAIL )
        p_tail_right = itaillimit(b);			
      if ( u < middlep + ap + bp + p_tail_left + p_tail_right ) {
        z = sample_tail( TAIL_LIMIT, b );
        continue;
      }			
      z = NAN;
    } while ( z != z );	
    return z;					
  }
  
public:	
  
  double draw( double mu, double sigma, double aa, double bb )
  {
    
    double a = (aa-mu)/sigma;
    double b = (bb-mu)/sigma;
    
    if ( a == -infinity && b == infinity )
      return NAN;		
    else
      return  sigma * truncate2(a,b) + mu;				
  }
  
  // creates the object by initializing the arrays we will need
  truncated()
  {
    
    rr = gsl_rng_alloc( gsl_rng_rand );
    
    double x = 0;
    int i = 0, c = 0;
    while( x < 3.5 ) {
      x += basep * exp(0.5*x*x);
      ++c;
    }
    xhist = (double*) malloc((2*c+1)*sizeof(double)) + c;
    dxhist =(double*) malloc(2*c*sizeof(double)) + c;
    yhist = (double*)malloc(2*(c+1)*sizeof(double)) + (c+1);
    
    yhistratio = (int*)malloc(2*c*sizeof(int)) + c;		
    whichbin = (int*)malloc( ((int)(TAIL_LIMIT * 128.0 )+1) * sizeof(int));
    
    x = 0; xhist[0] = 0; yhist[0] = 1; yhist[-1] = 1;
    
    for( i = 1; i <= c; ++i )  {
      x += basep * exp(0.5*x*x);	
      xhist[i] = x;
      xhist[-i] = -x;
      dxhist[i-1] = xhist[i] - xhist[i-1];
      dxhist[-i] = dxhist[i-1];			
      yhist[i] = exp(-0.5*x*x);
      yhist[-i-1] = yhist[i];
      yhistratio[i-1] = (int)ceil( yhist[i]/yhist[i-1] * (RAND_MAX+1.0) );
      yhistratio[-i] = yhistratio[i-1];
    }
    c = 0;
    for( i = 0; i <= (int)(TAIL_LIMIT / basep); ++i ) {
      while( (int)(xhist[c] / basep ) < i )
        ++c;
      whichbin[i] = c-1;
    }
  }	
};
//rtruncnorm


//  [[Rcpp::export]]
double txAx(const vector<double> &x, const matrix<double> &A)
{
  int n = x.size();
  matrix<double> vec_t(1,n);
  matrix<double> tt(1,n);
  matrix<double> vec(n,1);
  matrix<double> C(1,1);
  gsl_matrix_set_row(vec_t,0,x);
  gsl_matrix_set_col(vec,0,x);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, vec_t, A, 0.0, tt);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tt, vec, 0.0, C);
  double result = gsl_matrix_get(C,0,0);
  vec_t.free();
  tt.free();
  vec.free();
  C.free();
  return result;
}

// [[Rcpp::export]]
double VecSum(const vector<double> & x)
{
  double sum = 0;
  int n = x.size();
  for(int i = 0; i < n; ++i)
    sum+=gsl_vector_get(x, i);
  return sum;
}
//  [[Rcpp::export]]

void RowToMatrix(matrix<double> tmpX, const matrix<double> &X, int row, int k)
{
  int n = X.ncol();
  for(int i = 0; i < n/k; ++i)
    for(int j = 0; j < k; ++j)
      tmpX(j, i) = gsl_matrix_get(X,row, i * n/k + j); //use double(mat(i,j)) or gsl_matrix_get(mat,i,j)
}

void cpy(vector<double> A, Rcpp::NumericVector B)
{
  int n = A.size();
  for(int i = 0; i < n; ++i)
    A[i] = B[i];
}

void mvrnorm_cpp(const RcppGSL::Vector &mu, const RcppGSL::Matrix &cov, RcppGSL::Vector res)
{
  int times = 1; 
  gsl_set_error_handler_off();
  RcppGSL::matrix<double> chol_cov(cov.nrow(),cov.ncol());
  gsl_matrix_memcpy(chol_cov, cov);
  if(gsl_linalg_cholesky_decomp(chol_cov)) 
  {
    printf("Non def-pos!\n");
  }
  else
  {
    matrix<double> matMu(1,mu.size());
    int n = mu.size();
    RcppGSL::matrix<double> mat(times,n);
    int count = 0;
    for(int i = 0; i < times; ++i)
      for(int j = 0; j < n; ++j)
      {
        mat(i,j) = R::rnorm(0,1);
        count = count + 1;
        matMu(i,j) = mu[j];
      }
      
      gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans,CblasNonUnit,1.0, chol_cov,mat);
    gsl_matrix_add(matMu,mat);
    gsl_matrix_get_row(res, matMu, 0);
    matMu.free();
    mat.free();
    chol_cov.free();
  }
}

void GInverse(matrix<double> A)  
{  
  int n = A.ncol();  
  matrix<double> inverse(n, n);
  gsl_permutation *p = gsl_permutation_alloc(n);  
  int sign = 0;  
  gsl_linalg_LU_decomp(A, p, &sign);  
  gsl_linalg_LU_invert(A, p, inverse);  
  gsl_permutation_free(p);  
  gsl_matrix_memcpy(A, inverse);
  inverse.free();
}

void printf_matrix(const gsl_matrix * M)
{
    int nrow = M->size1;
    int ncol = M->size2;
    for(int i = 0; i < nrow; ++i)
    {
        for(int j = 0; j < ncol; ++j)
            printf("%f ", gsl_matrix_get(M, i, j));
        printf("\n");
    }
}

void printf_vector(const gsl_vector * x)
{
    int size = x->size;
    for(int i = 0; i < size; ++i)
        printf("%f ", gsl_vector_get(x, i));
    printf("\n");
}

//  [[Rcpp::export]]
void ACE_transform_cpp()
{
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  gsl_rng_default_seed = 0;
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  
  int n = 1000, n1 = 500, p = 2, MCAX = 2000, GNUM = 1800;
  matrix<double> Y(n, p);
  matrix<double> Ystar(n, p);
  matrix<double> X(n, 2 * p); 
  double tmpA[] = {1.0,0.5,0.5,1.0};
  gsl_matrix_view tmpA_view = gsl_matrix_view_array(tmpA, 2, 2);
  matrix<double> A(2, 2);
  gsl_matrix_memcpy(A, &tmpA_view.matrix);
  
  vector<double> Beta(2);
  Beta[0] = 1.0; Beta[1] = 0.5;
  double Sigma_A = 0.5;
  double Sigma_C = 0.5;
  double Sigma_E = 0.2;
  
  matrix<double> invA(2,2);
  gsl_matrix_memcpy(invA, A);
  GInverse(invA);
  matrix<double> Ua(n,2), Uc(n,2);
  
  int sp = 10;
  int MCMC_times = 500;
  Rcpp::NumericMatrix Est(sp,5), SEst(sp,5);
  gsl_matrix_set_zero(Ua);
  gsl_matrix_set_zero(Uc);
  
  matrix<double> Sigma_beta(p,p);
  vector<double> Mu_beta(p);
  
  double tmpmu[] = {0,0};
  gsl_vector_view tmpmu_view = gsl_vector_view_array(tmpmu, 2);
  vector<double> tmpMu(2);
  gsl_vector_memcpy(tmpMu, &tmpmu_view.vector);
  
  Rcpp::NumericMatrix MCMC(MCAX, 4);
  
  vector<double> Yi(2);
  matrix<double> Xi(2,2);
  vector<double> Uai(2);
  vector<double> Uci(2);
  matrix<double> tmpXi(2,2);
  vector<double> tmpx3(2);
  vector<double> tmpx4(2);
  vector<double> Ue(2);
  matrix<double> SigmaAA(2,2);
  matrix<double> Sigma_Emat(2, 2);
  matrix<double> InvADSigma_A(2, 2);  

 
  
  int kk = 25;  // the number of basis
                // thus nbreaks = 25+2-4=23
                // 
  matrix<double> D0(kk-1, kk);
  matrix<double> D1(kk-1-1, kk-1);
  double a[2] = {-1,1};
  
  matrix<double> N(kk, kk-2);
  matrix<double> M(kk, kk);
  vector<double> unit_e(kk); 
  vector<double> tmpMu_e(kk);
  gsl_vector_set_zero(tmpMu_e);
//   matrix<double> Bk(2*kk, n);
//   matrix<double> bk(2*kk, n);

  gsl_matrix * Bk[2];
  gsl_matrix * bk[2];
  for(int i = 0; i < 2; ++i)
  {
      Bk[i] = gsl_matrix_calloc(kk, n);
      Bk[i] = gsl_matrix_calloc(kk, n);
  }
  
  gsl_matrix * tmpBk = gsl_matrix_calloc(kk, kk);
  
  matrix<double> ks(2, kk+4-2);
  vector<double> Ycol(n);
  
  int bsorder = 4;
  int nbreak = kk + 2 - bsorder;
  gsl_bspline_workspace *bw = gsl_bspline_alloc(bsorder, nbreak);
  vector<double> ksRow(nbreak);
  
  vector<double> Bcoef(kk);
  int nderiv = 1;
  matrix<double> Bcoef_Driv(kk, nderiv+1);
  
  matrix<double> Gam(kk, 2);
  double tao[2];
  double sigma_r[] = {2.7, 2.7};
  vector<double> p_accept(2);
  
  matrix<double> f(2, 10); //an Aid Matrix?
  vector<double> f_row(f.ncol());
  // vector<double> y(100);
  double y;
  int BK_ncol = 100;
//   matrix<double> BK(2*kk, BK_ncol);
//   
//   matrix<double> GY(BK_ncol,2*(MCAX-GNUM));
  gsl_matrix * BK[2];
  gsl_matrix * GY[2];
  for(int i = 0; i < 2; ++i)
  {
      BK[i] = gsl_matrix_calloc(kk, BK_ncol);
      GY[i] = gsl_matrix_calloc(BK_ncol, (MCAX-GNUM));
  }


  matrix<double> gy(n, p);
  
  gsl_matrix * com_gy[sp];
  gsl_matrix * com_gy_beta[sp];
  for(int i = 0; i < sp; ++i)
  {
    com_gy[i] = gsl_matrix_calloc(n, 2);
    com_gy_beta[i] = gsl_matrix_calloc(n, 2);      
  }
  
  matrix<double> y_beta(n, p);
  double lower, upper;
  matrix<double> lo_up(2, kk);
  vector<double> lo_up_col(kk); 
  
  //我真是服了师兄r-code里面的各种神参数神变量名了...
  
  matrix<double> w1(kk, 10);
  vector<double> tmp_Unit_e(kk);
  vector<double> Gam_col(kk);
  vector<double> Bk__ij(kk);

  truncated truncnorm;
  double r1[10];
    
  for(int i = 0; i < n1; ++i)
  {
    Ua(i,0) = R::rnorm(0.0, std::sqrt(Sigma_A));
    Ua(i,1) = gsl_matrix_get(Ua,i,0);
    Uc(i,0) = R::rnorm(0.0, std::sqrt(Sigma_C));
    Uc(i,1) = gsl_matrix_get(Uc, i, 0);
    
    cpy(Ue, Rcpp::rnorm(2, 0.0, std::sqrt(Sigma_E)));
    
    for(int m = 0; m < 2; ++m) 
    {
      X(i, m) = R::rnorm(0.0, 2.0);
      X(i, m+2) = R::rnorm(0.0, 2.0);
    }
    
    RowToMatrix(Xi, X, i, 2);
    gsl_matrix_get_row(Uai, Ua, i);
    gsl_matrix_get_row(Uci, Uc, i);  
    gsl_vector_add(Uai, Uci);
    gsl_vector_add(Uai, Ue); 
    gsl_blas_dgemv(CblasNoTrans, 1.0, Xi, Beta, 1.0, Uai);
    gsl_matrix_set_row(Y, i, Uai);
  }
  
  gsl_matrix_memcpy(SigmaAA, A);
  gsl_matrix_scale(SigmaAA, Sigma_A);
  
  for(int i = n1; i < n; ++i)
  {
    mvrnorm_cpp(tmpMu, SigmaAA, tmpx3);
    gsl_matrix_set_row(Ua, i, tmpx3);
    Uc(i, 0) = R::rnorm(0.0, std::sqrt(Sigma_C));
    Uc(i, 1) = gsl_matrix_get(Uc, i, 0);
    
    cpy(Ue, Rcpp::rnorm(2, 0.0, std::sqrt(Sigma_E)));
    
    for(int m = 0; m < 2; ++m)
    {
      X(i, m) = R::rnorm(0.0, 2.0);
      X(i, m+2) = R::rnorm(0.0, 2.0);
    }
    
    RowToMatrix(Xi, X, i, 2);
    gsl_matrix_get_row(Uai, Ua, i);
    gsl_matrix_get_row(Uci, Uc, i);
    gsl_vector_add(Uai, Uci);
    gsl_vector_add(Uai, Ue);
    gsl_blas_dgemv(CblasNoTrans, 1.0, Xi, Beta, 1.0, Uai);
    gsl_matrix_set_row(Y, i, Uai);
  }
  
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < n; ++j)
      Ystar(i, j) = std::exp(gsl_matrix_get(Y,i,j));

  for(int i = 0; i < (kk-1); ++i)
  {
    D0(i,i) = a[0];
    D0(i,i+1) = a[1];
  }
  
  for(int i = 0; i < (kk-2); ++i)
  {
    D1(i,i) = a[0];
    D1(i,i+1) = a[1];
  }
  
  
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, D1, D0, 0, N);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, N, N, 0, M);
  
  for(int i = 0; i < 2; ++i)
  {
    gsl_matrix_get_col(Ycol,Y, i);
    double minYi = gsl_vector_min(Ycol);
    double maxYi = gsl_vector_max(Ycol);
    gsl_bspline_knots_uniform(minYi, maxYi, bw);
    gsl_matrix_set_row(ks, i, bw->knots);
    for(int k = 0; k < n; ++k)
    {
      double yj = gsl_vector_get(Ycol, k);
      gsl_bspline_deriv_eval(yj, nderiv, Bcoef_Driv, bw);   //something's wrong here
      gsl_matrix_get_col(Bcoef, Bcoef_Driv, 0); //Deriv of order 0
      gsl_matrix_set_row(Bk[i], k, Bcoef);
      gsl_matrix_get_col(Bcoef, Bcoef_Driv, 1); //Deriv of order 1
      gsl_matrix_set_row(bk[i], k, Bcoef);
      
    }
  }
  
  for(int i = 0; i < kk; ++i)
  {
    Gam(0, i) = (13.5-0.285)/24*i+0.285;
    Gam(1, i) = (15-0.285)/24*i+0.285;
  }
  
  for(int i = 0; i < 2; ++i)
    tao[i] = 1/(R::rgamma(1, 1/0.005)); 

  gsl_vector_set_zero(p_accept);
  gsl_matrix_set_zero(f);
  
  

  for(int i = 0; i < 2; ++i)
  {
    gsl_matrix_get_col(Ycol, Y, i);
    double minYi = gsl_vector_min(Ycol);
    double maxYi = gsl_vector_max(Ycol);
    gsl_bspline_knots_uniform(minYi, maxYi, bw);
    for(int j = 0; j < BK_ncol; ++j)
    {
      y = (maxYi - minYi)/(BK_ncol-1.0) * j + minYi;
      gsl_bspline_eval(y, Bcoef, bw);
//      printf_vector(Bcoef);
      gsl_matrix_set_row(BK[i], j, Bcoef);
    }
  }
  
//   for(int i = 0; i < 2; ++i)
//   {
//     gsl_matrix_set_zero(GY[i]);
//   }
//   gsl_matrix_set_zero(gy);
  
//   for(int i = 0; i < n1; ++i)
//   {
//     Ua(i,0) = R::rnorm(0.0, std::sqrt(Sigma_A));
//     Ua(i,1) = gsl_matrix_get(Ua,i,0);
//     Uc(i,0) = R::rnorm(0.0, std::sqrt(Sigma_C));
//     Uc(i,1) = gsl_matrix_get(Uc, i, 0);
//   }
  
//   gsl_matrix_memcpy(SigmaAA, A);
//   gsl_matrix_scale(SigmaAA, Sigma_A);
  
  
  //从这里开始，正态分布随机数生成器崩溃？
//   for(int i = n1; i < n; ++i)
//   {
//     mvrnorm_cpp(tmpMu, SigmaAA, tmpx3);
//     gsl_matrix_set_row(Ua, i, tmpx3);
//     Uc(i, 0) = R::rnorm(0.0, std::sqrt(Sigma_C));
//     Uc(i, 1) = gsl_matrix_get(Uc, i, 0);
//   }
  
//   for(int GIB = 0; GIB < MCAX; ++GIB)
//   {
//       for(int j = 0; j < 2; ++j)
//       {
//           gsl_matrix_scale(M, 1/tao[j]);
//           gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1,Bk[j], Bk[j],0,M);
//           GInverse(M);
//           if(GIB==0)
//           {
//               printf_matrix(M);
//           }
//           mvrnorm_cpp(tmpMu_e, M, unit_e);
//           gsl_vector_memcpy(tmp_Unit_e, unit_e);
//           gsl_vector_mul(tmp_Unit_e, tmp_Unit_e);
//           gsl_vector_scale(unit_e, 1/(std::sqrt(VecSum(tmp_Unit_e))));
          
//           for(int i = 0; i < kk; ++i)
//           {
//               lo_up(0, i) = std::numeric_limits<double>::min();
//               lo_up(1, i) = std::numeric_limits<double>::max();
//           }
          
//           for(int ii = 0; ii < (kk-1); ++ii)
//           {
//               if(unit_e[ii+1] > unit_e[ii])
//                 lo_up(0, ii) = (Gam(ii, j) - Gam(ii+1, j)) / (unit_e[ii+1] - unit_e[ii]);
//               else
//                 lo_up(1, ii) = (Gam(ii, j) - Gam(ii+1, j)) / (unit_e[ii+1] - unit_e[ii]);
//           }
          
//           gsl_matrix_get_col(lo_up_col, lo_up, 0);
//           lower = gsl_vector_max(lo_up_col);
//           gsl_matrix_get_col(lo_up_col, lo_up, 1);
//           upper = gsl_vector_min(lo_up_col);
          
//           for(int i = 0; i < 10; ++i)
//             r1[i] = truncnorm.draw(0, sigma_r[j], lower, upper);
          
          

          
             
          
//       }
//   }

    
  
  
  
}