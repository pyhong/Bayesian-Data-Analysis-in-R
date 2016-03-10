
//  [[Rcpp::depends(RcppGSL)]]

// begin of rtuncnorm.cpp
#if !defined(__SEXP_MACROS_H__)
#define __SEXP_MACROS_H__

#include <R.h>
#include <Rinternals.h>

#define CHECK_ARG_IS_REAL_MATRIX(A)					\
if (!isReal(A) || !isMatrix(A))					    \
  error("Argument '" #A "' is not a real matrix.");

#define CHECK_ARG_IS_REAL_VECTOR(A)					\
if (!isReal(A) || !isVector(A))					    \
  error("Argument '" #A "' is not a real vector.");

#define CHECK_ARG_IS_INT_VECTOR(A)					\
if (!isInteger(A) || !isVector(A))					\
  error("Argument '" #A "' is not an integer vector.");

/*
* Unpack a real matrix stored in SEXP S. 
*/
#define UNPACK_REAL_MATRIX(S, D, N, K)          \
CHECK_ARG_IS_REAL_MATRIX(S);		                  \
double *D = REAL(S);			                         \
const R_len_t N = nrows(S);			                  \
const R_len_t K = ncols(S);

/*
* Unpack a real vector stored in SEXP S.
*/
#define UNPACK_REAL_VECTOR(S, D, N)             \
CHECK_ARG_IS_REAL_VECTOR(S);		                  \
double *D = REAL(S);			                         \
const R_len_t N = length(S);                   

/*
* Unpack a single real stored in SEXP S.
*/
#define UNPACK_REAL(S, D)			  \
CHECK_ARG_IS_REAL_VECTOR(S);		\
double D = REAL(S)[0];			     \

/*
 * Unpack an integer vector stored in SEXP S.
 */
#define UNPACK_INT_VECTOR(S, I, N)             \
CHECK_ARG_IS_INT_VECTOR(S);		                  \
int *I = INTEGER(S);		                         \
const R_len_t N = length(S);                   

/*
 * Unpack a single integer stored in SEXP S.
 */
#define UNPACK_INT(S, I)			   \
CHECK_ARG_IS_INT_VECTOR(S);			\
int I = INTEGER(S)[0];			     \

#endif
#include <R.h>
#include <Rcpp.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <float.h>



#define ALLOC_REAL_VECTOR(S, D, N)		         \
SEXP S;					                                 \
PROTECT(S = allocVector(REALSXP, N));	       \
double *D = REAL(S);

#ifndef MAX
#define MAX(A, B) ((A>B)?(A):(B))
#endif

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Applic.h>

#ifdef DEBUG
#define SAMPLER_DEBUG(N, A, B) Rprintf("%8s(%f, %f)\n", N, A, B)
#else
#define SAMPLER_DEBUG(N, A, B) 
#endif

static const double t1 = 0.15;
static const double t2 = 2.18;
static const double t3 = 0.725;
static const double t4 = 0.45;

/* Exponential rejection sampling (a,inf) */
// [[Rcpp::export]]
static R_INLINE double ers_a_inf(double a) {
  SAMPLER_DEBUG("ers_a_inf", a, R_PosInf);
  const double ainv = 1.0 / a;
  double x, rho;
  do{
    x = R::rexp(ainv) + a; /* rexp works with 1/lambda */
rho = exp(-0.5 * pow((x - a), 2));
  } while (R::runif(0, 1) > rho);
  return x;
}

/* Exponential rejection sampling (a,b) */
static R_INLINE double ers_a_b(double a, double b) {
  SAMPLER_DEBUG("ers_a_b", a, b);
  const double ainv = 1.0 / a;
  double x, rho;
  do{
    x = R::rexp(ainv) + a; /* rexp works with 1/lambda */
rho = exp(-0.5 * pow((x-a), 2));
  } while (R::runif(0, 1) > rho || x > b);
  return x;
}

/* Normal rejection sampling (a,b) */
static R_INLINE double nrs_a_b(double a, double b){
  SAMPLER_DEBUG("nrs_a_b", a, b);
  double x = -DBL_MAX;
  while(x < a || x > b){
    x = R::rnorm(0, 1);
  }
  return x;
}

/* Normal rejection sampling (a,inf) */
static R_INLINE double nrs_a_inf(double a){
  SAMPLER_DEBUG("nrs_a_inf", a, R_PosInf);
  double x = -DBL_MAX;
  while(x < a){
    x = R::rnorm(0, 1);
  }
  return x;
}

/* Half-normal rejection sampling */
double hnrs_a_b(double a, double b){
  SAMPLER_DEBUG("hnrs_a_b", a, b);    
  double x = a - 1.0;
  while(x < a || x > b) {
    x = R::rnorm(0, 1);
    x = fabs(x);
  }
  return x;
}

/* Uniform rejection sampling */
static R_INLINE double urs_a_b(double a, double b){
  SAMPLER_DEBUG("urs_a_b", a, b);
  const double phi_a = R::dnorm(a, 0.0, 1.0, FALSE);
  double x = 0.0;//, u = 0.0;
  
  /* Upper bound of normal density on [a, b] */
  const double ub = a < 0 && b > 0 ? M_1_SQRT_2PI : phi_a;
  do{
    x = R::runif(a, b);
  } while (R::runif(0, 1) * ub > R::dnorm(x,0,1,0));
  return x;
}

/* Previously this was refered to as type 1 sampling: */
static inline double r_lefttruncnorm(double a, double mean, double sd) {
  const double alpha = (a - mean) / sd;
  if (alpha < t4) {
    return mean + sd * nrs_a_inf(alpha);
  } else {
    return mean + sd * ers_a_inf(alpha);
  }
}

static R_INLINE double r_righttruncnorm(double b, double mean, double sd) {
  const double beta = (b - mean) / sd;
  /* Exploit symmetry: */
  return mean - sd * r_lefttruncnorm(-beta, 0.0, 1.0);
}

static R_INLINE double r_truncnorm(double a, double b, double mean, double sd) {
  const double alpha = (a - mean) / sd;
  const double beta = (b - mean) / sd;
  const double phi_a = R::dnorm(alpha, 0.0, 1.0, FALSE);
  const double phi_b = R::dnorm(beta, 0.0, 1.0, FALSE);
  if (beta <= alpha) {
    return NA_REAL;
  } else if (alpha <= 0 && 0 <= beta) { /* 2 */
  if (phi_a <= t1 || phi_b <= t1) { /* 2 (a) */
  return mean + sd * nrs_a_b(alpha, beta);
  } else { /* 2 (b) */
  return mean + sd * urs_a_b(alpha, beta);
  }
  } else if (alpha > 0) { /* 3 */
  if (phi_a / phi_b <= t2) { /* 3 (a) */
  return mean + sd * urs_a_b(alpha, beta);
  } else {
    if (alpha < t3) { /* 3 (b) */                
  return mean + sd * hnrs_a_b(alpha, beta);
    } else { /* 3 (c) */
  return mean + sd * ers_a_b(alpha, beta);
    }
  }
  } else { /* 3s */
  if (phi_b / phi_a <= t2) { /* 3s (a) */
  return mean - sd * urs_a_b(-beta, -alpha);
  } else {
    if (beta > -t3) { /* 3s (b) */
  return mean - sd * hnrs_a_b(-beta, -alpha);
    } else { /* 3s (c) */
  return mean - sd * ers_a_b(-beta, -alpha);
    }
  }
  }
}

//[[Rcpp::export]]
SEXP do_rtruncnorm(SEXP s_n, SEXP s_a, SEXP s_b, SEXP s_mean, SEXP s_sd) {
  R_len_t i, nn;
  UNPACK_INT(s_n, n);
  if (NA_INTEGER == n)
    error("n is NA - aborting.");
  UNPACK_REAL_VECTOR(s_a   , a   , n_a);
  UNPACK_REAL_VECTOR(s_b   , b   , n_b);
  UNPACK_REAL_VECTOR(s_mean, mean, n_mean);
  UNPACK_REAL_VECTOR(s_sd  , sd  , n_sd);
  
  nn = MAX(n, MAX(MAX(n_a, n_b), MAX(n_mean, n_sd)));
  ALLOC_REAL_VECTOR(s_ret, ret, nn);
  
  GetRNGstate();
  for (i = 0; i < nn; ++i) {
    const double ca = a[i % n_a];
    const double cb = b[i % n_b];
    const double cmean = mean[i % n_mean];
    const double csd = sd[i % n_sd];
    
    if (R_FINITE(ca) && R_FINITE(cb)) {
      ret[i] = r_truncnorm(ca, cb, cmean, csd);
    } else if (R_NegInf == ca && R_FINITE(cb)) {
      ret[i] = r_righttruncnorm(cb, cmean, csd);
    } else if (R_FINITE(ca) && R_PosInf == cb) {
      ret[i] = r_lefttruncnorm(ca, cmean, csd);
    } else if (R_NegInf == ca && R_PosInf == cb) {
      ret[i] = R::rnorm(cmean, csd);
    } else {
      ret[i] = NA_REAL;
    }
    R_CheckUserInterrupt();
  }
  PutRNGstate();
  UNPROTECT(1); /* s_ret */
  return s_ret;
}

// end of rtuncnorm.cpp


#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <limits.h>

//using namespace RcppGSL;

//  [[Rcpp::export]]
double txAx(const RcppGSL::vector<double> &x, const RcppGSL::matrix<double> &A)
{
  int n = x.size();
  RcppGSL::matrix<double> vec_t(1,n);
  RcppGSL::matrix<double> tt(1,n);
  RcppGSL::matrix<double> vec(n,1);
  RcppGSL::matrix<double> C(1,1);
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
double VecSum(const RcppGSL::vector<double> & x)
{
  double sum = 0;
  int n = x.size();
  for(int i = 0; i < n; ++i)
    sum+=gsl_vector_get(x, i);
  return sum;
}
//  [[Rcpp::export]]

void RowToMatrix(RcppGSL::matrix<double> tmpX, const RcppGSL::matrix<double> &X, int row, int k)
{
  int n = X.ncol();
  for(int i = 0; i < n/k; ++i)
    for(int j = 0; j < k; ++j)
      tmpX(j, i) = gsl_matrix_get(X,row, i * n/k + j); //use double(mat(i,j)) or gsl_matrix_get(mat,i,j)
}

void cpy(RcppGSL::vector<double> A, Rcpp::NumericVector B)
{
  int n = A.size();
  for(int i = 0; i < n; ++i)
    A[i] = B[i];
}

void mvrnorm_cpp(const RcppGSL::vector<double> &mu, const RcppGSL::matrix<double> &cov, RcppGSL::vector<double> res)
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
    RcppGSL::matrix<double> matMu(1,mu.size());
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

void GInverse(RcppGSL::matrix<double> A)  
{  
  int n = A.ncol();  
  RcppGSL::matrix<double> inverse(n, n);
  gsl_permutation *p = gsl_permutation_alloc(n);  
  int sign = 0;  
  gsl_linalg_LU_decomp(A, p, &sign);  
  gsl_linalg_LU_invert(A, p, inverse);  
  gsl_permutation_free(p);  
  gsl_matrix_memcpy(A, inverse);
  inverse.free();
}

double txAy(const gsl_vector * x, const gsl_matrix *A, const gsl_matrix *y)
{
    RcppGSL::vector<double> a(x);
    RcppGSL::vector<double> c(y);
    RcppGSL::matrix<double> b(A);
    RcppGSL::matrix<double> tmp1(1,a.size());
    RcppGSL::matrix<double> tmp2(1,c.size());
    RcppGSL::matrix<double> tmp3(c.size(),1);
    RcppGSL::matrix<double> tmp4(1,1);
    gsl_matrix_set_row(tmp1, 0, a);
    gsl_matrix_set_col(tmp3, 0, c);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, tmp1, b, 0, tmp2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, tmp2, tmp3, 0, tmp4);
    double result = gsl_matrix_get(tmp4, 0, 0);
    a.free();
    c.free();
    b.free();
    tmp1.free();
    tmp2.free();
    tmp3.free();
    tmp4.free();
}

//  [[Rcpp::export]]
void ACE_transform_cpp()
{
  int n = 1000, n1 = 500, p = 2, MCAX = 2000, GNUM = 1800;
  RcppGSL::matrix<double> Y(n, p);
  RcppGSL::matrix<double> Ystar(n, p);
  RcppGSL::matrix<double> X(n, 2 * p); 
  double tmpA[] = {1.0,0.5,0.5,1.0};
  gsl_matrix_view tmpA_view = gsl_matrix_view_array(tmpA, 2, 2);
  RcppGSL::matrix<double> A(2, 2);
  gsl_matrix_memcpy(A, &tmpA_view.matrix);
  
  RcppGSL::vector<double> Beta(2);
  Beta[0] = 1.0; Beta[1] = 0.5;
  double Sigma_A = 0.5;
  double Sigma_C = 0.5;
  double Sigma_E = 0.2;
  
  RcppGSL::matrix<double> invA(2,2);
  gsl_matrix_memcpy(invA, A);
  GInverse(invA);
  RcppGSL::matrix<double> Ua(n,2), Uc(n,2);
  
  int sp = 10;
  int MCMC_times = 500;
  Rcpp::NumericMatrix Est(sp,5), SEst(sp,5);
  gsl_matrix_set_zero(Ua);
  gsl_matrix_set_zero(Uc);
  
  RcppGSL::matrix<double> Sigma_beta(p,p);
  RcppGSL::vector<double> Mu_beta(p);
  
  double tmpmu[] = {0,0};
  gsl_vector_view tmpmu_view = gsl_vector_view_array(tmpmu, 2);
  RcppGSL::vector<double> tmpMu(2);
  gsl_vector_memcpy(tmpMu, &tmpmu_view.vector);
  
  Rcpp::NumericMatrix MCMC(MCAX, 4);
  
  RcppGSL::vector<double> Yi(2);
  RcppGSL::matrix<double> Xi(2,2);
  RcppGSL::vector<double> Uai(2);
  RcppGSL::vector<double> Uci(2);
  RcppGSL::matrix<double> tmpXi(2,2);
  RcppGSL::vector<double> tmpx3(2);
  RcppGSL::vector<double> tmpx4(2);
  RcppGSL::vector<double> Ue(2);
  RcppGSL::matrix<double> SigmaAA(2,2);
  RcppGSL::matrix<double> Sigma_Emat(2, 2);
  RcppGSL::matrix<double> InvADSigma_A(2, 2);  

 
  
  int kk = 25;  // the number of basis
                // thus nbreaks = 25+2-4=23
                // 
  RcppGSL::matrix<double> D0(kk-1, kk);
  RcppGSL::matrix<double> D1(kk-1-1, kk-1);
  double a[2] = {-1,1};
  
  RcppGSL::matrix<double> N(kk, kk-2);
  RcppGSL::matrix<double> M(kk, kk);
  RcppGSL::vector<double> unit_e(kk); 
//   RcppGSL::matrix<double> Bk(2*kk, n);
//   RcppGSL::matrix<double> bk(2*kk, n);

  gsl_matrix * Bk[2];
  gsl_matrix * bk[2];
  for(int i = 0; i < 2; ++i)
  {
      Bk[i] = gsl_matrix_calloc(kk, n);
      Bk[i] = gsl_matrix_calloc(kk, n);
  }
  
  gsl_matrix * tmpBk = gsl_matrix_calloc(kk, kk);
  
  RcppGSL::matrix<double> ks(2, kk+4-2);
  RcppGSL::vector<double> Ycol(n);
  
  int bsorder = 4;
  int nbreak = kk + 2 - bsorder;
  gsl_bspline_workspace *bw = gsl_bspline_alloc(bsorder, nbreak);
  RcppGSL::vector<double> ksRow(nbreak);
  
  RcppGSL::vector<double> Bcoef(kk);
  int nderiv = 1;
  RcppGSL::matrix<double> Bcoef_Driv(nderiv, kk);
  
  RcppGSL::matrix<double> Gam(kk, 2);
  double tao[2];
  double sigma_r[] = {2.7, 2.7};
  RcppGSL::vector<double> p_accept(2);
  
  RcppGSL::matrix<double> f(2, 10); //an Aid Matrix?
  RcppGSL::vector<double> f_row(f.ncol());
  // RcppGSL::vector<double> y(100);
  double y;
  int BK_ncol = 100;
//   RcppGSL::matrix<double> BK(2*kk, BK_ncol);
//   
//   RcppGSL::matrix<double> GY(BK_ncol,2*(MCAX-GNUM));
  gsl_matrix * BK[2];
  gsl_matrix * GY[2];
  for(int i = 0; i < 2; ++i)
  {
      Bk[i] = gsl_matrix_calloc(kk, BK_ncol);
      GY[i] = gsl_matrix_calloc(BK_ncol, (MCAX-GNUM));
  }


  RcppGSL::matrix<double> gy(n, p);
  
  gsl_matrix * com_gy[sp];
  gsl_matrix * com_gy_beta[sp];
  for(int i = 0; i < sp; ++i)
  {
    com_gy[i] = gsl_matrix_calloc(n, 2);
    com_gy_beta[i] = gsl_matrix_calloc(n, 2);      
  }
  
  RcppGSL::matrix<double> y_beta(n, p);
  double lower, upper;
  RcppGSL::matrix<double> lo_up(2, kk);
  RcppGSL::vector<double> lo_up_col(kk); 
  
  //我真是服了师兄r-code里面的各种神参数神变量名了……
  
  RcppGSL::matrix<double> w1(kk, 10);
  RcppGSL::vector<double> tmp_Unit_e(kk);
  RcppGSL::vector<double> Gam_col(kk);
  RcppGSL::vector<double> Bk__ij(kk);
  //RcppGSL::vector<double> Bk__ij()
  
  //RcppGSL::vector<double> p1(f.ncol());
    
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
          
          Rcpp::NumericVector r1 = do_rtruncnorm(10, lower, upper, 0.0, sigma_r[j]);
          
          for(int ww = 0; w < 10; ++ww)
          {
              gsl_vector_memcpy(tmp_Unit_e, unit_e);
              gsl_vector_scale(r1[ww], tmp_Unit_e);
              gsl_matrix_get_col(Gam_col, Gam, j);
              gsl_vector_add(Gam_col, tmp_Unit_e;
              gsl_matrix_set_col(w1, ww, Gam_col);
              double f1 = 0, f2, f3;
              for(int i = 0; i < n; ++i)
              {
                 gsl_matrix_get_col(Bk__ij ,Bk[j], i);
                 gsl_blas_ddot(Gam_col, Bk__ij, f2);
                 f2 -= double(Ua(i, j))- double(Uc(i, j));
                 f2 *= f2;
                 gsl_matrix_get_col(Bk__ij, bk[j], i);
                 gsl_blas_ddot(Gam_col, Bk__ij, f3);
                 f3 += 2 * Sigma_E;
                 f1 = f1 - f2 / f3;              
              }
              
              f(0, ww) = f1 - txAx(Gam_col, M) / (2 * tao[j]) - double(r1[j])^2/(2*sigma_r[j]^2);
          }
          
          gsl_matrix_get_row(f_row, f, 0);
          double f1_m = gsl_vector_max(f_row);
          
          for(int ff = 0; ff < 10; ++ff)
            f(0, ff) = std::exp(double(f(0, ff)) - f1_m);
            
          gsl_vector_scale(1/VecSum(f_row), f_row);
          
             
          
      }
  }

    
  
  
  
}