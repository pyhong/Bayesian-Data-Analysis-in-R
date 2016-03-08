//  [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <limits.h>
//#include "./src/rtuncnorm.h"

using namespace RcppGSL;

//  [[Rcpp::export]]
double txAx(const Vector &x, const Matrix &A)
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
double VecSum(const Vector & x)
{
  double sum = 0;
  int n = x.size();
  for(int i = 0; i < n; ++i)
    sum+=gsl_vector_get(x, i);
  return sum;
}
//  [[Rcpp::export]]

void RowToMatrix(Matrix tmpX, const Matrix &X, int row, int k)
{
  int n = X.ncol();
  for(int i = 0; i < n/k; ++i)
    for(int j = 0; j < k; ++j)
      tmpX(j, i) = gsl_matrix_get(X,row, i * n/k + j); //use double(mat(i,j)) or gsl_matrix_get(mat,i,j)
}

void cpy(Vector A, Rcpp::NumericVector B)
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
    Rcpp::NumericVector rdvec = Rcpp::rnorm(n * times,0,1);
    RcppGSL::matrix<double> mat(times,n);
    int count = 0;
    for(int i = 0; i < times; ++i)
      for(int j = 0; j < n; ++j)
      {
        mat(i,j) = rdvec[count];
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

//  [[Rcpp::export]]
void ACE_transform_cpp()
{
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
  vector<double> unit_e(2);
 
  
  int kk = 25;  // the number of basis
                // thus nbreaks = 25+2-4=23
                // 
  matrix<double> D0(kk-1, kk);
  matrix<double> D1(kk-1-1, kk-1);
  double a[2] = {-1,1};
  
  matrix<double> N(kk, kk-2);
  matrix<double> M(kk, kk);
  
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
  matrix<double> Bcoef_Driv(nderiv, kk);
  
  matrix<double> Gam(kk, 2);
  double tao[2];
  double sigma_r[] = {2.7, 2.7};
  vector<double> p_accept(2);
  
  matrix<double> f(2, 10); //an Aid Matrix?
  
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
      Bk[i] = gsl_matrix_calloc(kk, BK_ncol);
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
  
    
  for(int i = 0; i < n1; ++i)
  {
    Ua(i,0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_A)))[0];
    Ua(i,1) = gsl_matrix_get(Ua,i,0);
    Uc(i,0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_C)))[0];
    Uc(i,1) = gsl_matrix_get(Uc,i,0);
    
    cpy(Ue, Rcpp::rnorm(2, 0.0, std::sqrt(Sigma_E)));
    
    for(int m = 0; m < 2; ++m) 
    {
      X(i,m) = (Rcpp::rnorm(1, 0.0, 2))[0];
      X(i,m+2) = (Rcpp::rnorm(1, 0.0, 2))[0];
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
    Uc(i, 0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_C)))[0];
    Uc(i, 1) = gsl_matrix_get(Uc, i, 0);
    
    cpy(Ue, Rcpp::rnorm(2, 0.0, std::sqrt(Sigma_E)));
    
    for(int m = 0; m < 2; ++m)
    {
      X(i, m) = (Rcpp::rnorm(1, 0.0, 2.0))[0];
      X(i, m+2) = (Rcpp::rnorm(1, 0.0, 2.0))[0];
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
      gsl_bspline_eval(yj, Bcoef, bw);
      gsl_matrix_set_row(Bk[i], k, Bcoef);
      
      gsl_bspline_deriv_eval(yj, nderiv, Bcoef_Driv, bw);
      gsl_matrix_get_row(Bcoef, Bcoef_Driv, 0);
      gsl_matrix_set_row(bk[i], k, Bcoef);
    }
  }
  
  for(int i = 0; i < kk; ++i)
  {
    Gam(0, i) = (13.5-0.285)/24*i+0.285;
    Gam(1, i) = (15-0.285)/24*i+0.285;
  }
  
  for(int i = 0; i < 2; ++i)
    tao[i] = 1/(Rcpp::rgamma(1, 1, 1/0.005)[0]);

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
      gsl_matrix_set_row(BK[i], j, Bcoef);
    }
  }
  
  for(int i = 0; i < 2; ++i)
  {
    gsl_matrix_set_zero(GY[i]);
  }
  gsl_matrix_set_zero(gy);
  
  for(int i = 0; i < n1; ++i)
  {
    Ua(i,0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_A)))[0];
    Ua(i,1) = gsl_matrix_get(Ua,i,0);
    Uc(i,0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_C)))[0];
    Uc(i,1) = gsl_matrix_get(Uc,i,0);
  }
  
  gsl_matrix_memcpy(SigmaAA, A);
  gsl_matrix_scale(SigmaAA, Sigma_A);
  
  for(int i = n1; i < n; ++i)
  {
    mvrnorm_cpp(tmpMu, SigmaAA, tmpx3);
    gsl_matrix_set_row(Ua, i, tmpx3);
    Uc(i, 0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_C)))[0];
    Uc(i, 1) = gsl_matrix_get(Uc, i, 0);
  }
  
  for(int GIB = 0; GIB < MCAX; ++GIB)
  {
      for(int j = 0; j < 2; ++j)
      {
          gsl_matrix_scale(M, 1/tao[j]);
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1,Bk[j], Bk[j],0,M);
          GInverse(M);
          mvrnorm_cpp(tmpMu, M, unit_e);
          gsl_vector_memcpy(tmpx4, unit_e);
          gsl_vector_mul(tmpx4, tmpx4);
          gsl_vector_scale(unit_e, 1/(std::sqrt(VecSum(tmpx4))));
          
          for(int i = 0; i < kk; ++i)
          {
              lo_up(0, i) = std::numeric_limits<double>::min();
              lo_up(1, i) = std::numeric_limits<double>::max();
          }
          
          for(int ii = 0; ii < (kk-1); ++ii)
          {
              if(unit_e[ii+1] > unit_e[ii])
                lo_up(0, ii) = (Gam(ii, j) - Gam(ii+1, j)) / (unit_e[ii+1] - unit_e[ii]);
              else
                lo_up(1, ii) = (Gam(ii, j) - Gam(ii+1, j)) / (unit_e[ii+1] - unit_e[ii]);
          }
          
          gsl_matrix_get_col(lo_up_col, lo_up, 0);
          lower = gsl_vector_max(lo_up_col);
          gsl_matrix_get_col(lo_up_col, lo_up, 1);
          upper = gsl_vector_min(lo_up_col);
          
          
          
      }
  }

    
  
  
  
}