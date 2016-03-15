//  [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppGSL.h>
#include <Rmath.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>

using namespace RcppGSL;

double txAx(const gsl_vector *x, const gsl_matrix *A)
{
  int n = x->size;
  gsl_matrix * vec_t = gsl_matrix_calloc(1,n);
  gsl_matrix * tt = gsl_matrix_calloc(1,n);
  gsl_matrix * vec = gsl_matrix_calloc(n,1);
  gsl_matrix * C = gsl_matrix_calloc(1,1);
  gsl_matrix_set_row(vec_t,0,x);
  gsl_matrix_set_col(vec,0,x);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, vec_t, A, 0.0, tt);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tt, vec, 0.0, C);
  double result = gsl_matrix_get(C,0,0);
  gsl_matrix_free(vec_t);
  gsl_matrix_free(tt);
  gsl_matrix_free(vec);
  gsl_matrix_free(C);
  return result;
}


double VecSum(const gsl_vector * x)
{
  double sum = 0;
  int n = x->size;
  for(int i = 0; i < n; ++i)
    sum+=gsl_vector_get(x, i);
  return sum;
}


void RowToMatrix(gsl_matrix * tmpX, const gsl_matrix * X, int row, int k)
{
  int n = X->size2;
  for(int i = 0; i < n/k; ++i)
    for(int j = 0; j < k; ++j)
    {
      double tt = gsl_matrix_get(X,row, i * n/k + j); //use double(mat(i,j)) or gsl_matrix_get(mat,i,j)
      gsl_matrix_set(tmpX, j, i, tt);
    }
}

void cpy(gsl_vector * A, Rcpp::NumericVector B)
{
  int n = A->size;
  for(int i = 0; i < n; ++i)
    gsl_vector_set(A,i,B[i]);
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
//  [[Rcpp::export]]
extern "C" SEXP ACE_standard_cpp()
{
  int n = 500, n1 = 250, p = 2;
  matrix<double> Y(n, p);
  matrix<double> X(n, 2 * p); 
  double tmpA[] = {1.0,0.5,0.5,1.0};
  gsl_matrix_view tmpA_view = gsl_matrix_view_array(tmpA, 2, 2);
  matrix<double> A(2, 2);
  gsl_matrix_memcpy(A, &tmpA_view.matrix);
  
  vector<double> Beta(2);
  double Sigma_A = 0.4;
  double Sigma_C = 0.3;
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

  Rcpp::NumericMatrix MCMC(MCMC_times, 5);
  
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
  
  for(int CIR = 0; CIR < sp; ++CIR)
  {
    Beta[0] = 1.0; Beta[1] = 0.5;
    Sigma_A = 0.4;
    Sigma_C = 0.3;
    Sigma_E = 0.2;

    
    for(int i = 0; i < n1; ++i)
    {
      Ua(i,0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_A)))[0];
      Ua(i,1) = gsl_matrix_get(Ua,i,0);
      Uc(i,0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_C)))[0];
      Uc(i,1) = gsl_matrix_get(Uc,i,0);
      
      double tmpx1 = (Rcpp::runif(1, 5, 13))[0];
      cpy(Ue, Rcpp::rnorm(2, 0.0, std::sqrt(Sigma_E)));

      for(int m = 0; m < 2; ++m) 
      {
        X(i,m) = tmpx1;
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
      
      double tmpx1 = (Rcpp::runif(1, 5, 13))[0];
      cpy(Ue, Rcpp::rnorm(2, 0.0, std::sqrt(Sigma_E)));

      for(int m = 0; m < 2; ++m)
      {
        X(i, m) = tmpx1;
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

    Beta[0] = 0.0; Beta[1] = 0.0;
    Sigma_A = Sigma_C = Sigma_E = 1.0;
    
    
    
    for(int i = 0; i < n1; ++i)
    {
      Ua(i, 0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_A)))[0];
      Ua(i, 1) = gsl_matrix_get(Ua,i,0);
      Uc(i, 0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_C)))[0];
      Uc(i, 1) = gsl_matrix_get(Uc,i,0);
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
    
    for(int GIB = 0; GIB < MCMC_times; ++GIB)
    {
      double Sigma_A2 = 1 / (2/Sigma_E + 1/Sigma_A);

      for(int i = 0; i < n1; ++i)
      {
        gsl_matrix_get_row(Yi, Y, i);
        RowToMatrix(Xi, X, i, 2);
        gsl_matrix_get_row(Uci, Uc, i);
        gsl_blas_dgemv(CblasNoTrans, 1, Xi, Beta, 1, Uci);
        gsl_vector_sub(Yi,Uci);
        double Mu_A = VecSum(Yi)*Sigma_A2/Sigma_E;
        Ua(i, 0) = (Rcpp::rnorm(1, Mu_A, std::sqrt(Sigma_A2)))[0];
        Ua(i, 1) = gsl_matrix_get(Ua,i,0);
      }

      gsl_matrix_set_identity(Sigma_Emat);
      gsl_matrix_scale(Sigma_Emat, 1/Sigma_E);
      gsl_matrix_memcpy(InvADSigma_A, invA);
      gsl_matrix_scale(InvADSigma_A, 1/Sigma_A);
      gsl_matrix_add(Sigma_Emat, InvADSigma_A);
      GInverse(Sigma_Emat);

      for(int i = n1; i < n; ++i)
      {
        gsl_matrix_get_row(Yi, Y, i);
        RowToMatrix(Xi, X, i, 2);
        gsl_matrix_get_row(Uci, Uc, i);
        gsl_blas_dgemv(CblasNoTrans, 1.0, Xi, Beta, 1.0, Uci);
        gsl_vector_sub(Yi, Uci);
        gsl_blas_dgemv(CblasNoTrans, 1/Sigma_E, Sigma_Emat, Yi, 0.0, tmpx4);
        mvrnorm_cpp(tmpx4, Sigma_Emat, tmpx3);
        gsl_matrix_set_row(Ua, i, tmpx3);
      }

      double Sigma_C2 = 1/(2/Sigma_E+1/Sigma_C);
      
      for(int i = 0; i < n; ++i)
      {
        gsl_matrix_get_row(Yi, Y, i);
        RowToMatrix(Xi, X, i, 2);
        gsl_matrix_get_row(Uai, Ua, i);
        gsl_blas_dgemv(CblasNoTrans, 1.0, Xi, Beta, 1.0, Uai);
        gsl_vector_sub(Yi, Uai);
        double Mu_C = VecSum(Yi)*Sigma_C2/Sigma_E;
        Uc(i, 0) = (Rcpp::rnorm(1, Mu_C, std::sqrt(Sigma_C2)))[0];
        Uc(i, 1) = gsl_matrix_get(Uc, i, 0);
      }
      
      double a = (n1+2*(n-n1)+1)/2, b = 0.0;
      for(int i = 0; i < n1; ++i)
      {
        double tUa = gsl_matrix_get(Ua, i, 0);
        tUa *= tUa;
        b += tUa; 
      }

      for(int i = n1; i < n; ++i)
      {
        gsl_matrix_get_row(Uai, Ua, i);
        b += txAx(Uai, invA);
      }

      b = b/2;
      Sigma_A = 1/((Rcpp::rgamma(1,a,1/b))[0]);
      MCMC(GIB,0) = Sigma_A;
      
      b = 0.0;
      a = (n+1)/2;
      vector<double> Ucj(Uc.nrow());
      gsl_matrix_get_col(Ucj, Uc, 0);
      gsl_vector_mul(Ucj, Ucj);
      b = VecSum(Ucj);
      Ucj.free();
      
      b = b/2;
      Sigma_C = 1/((Rcpp::rgamma(1,a,1/b))[0]);
      MCMC(GIB, 1) = Sigma_C;
      
      a=(2*n+1)/2; b=0;

      for(int i = 0; i < n; ++i)
      {
        gsl_matrix_get_row(Yi,Y,i);
        RowToMatrix(Xi, X, i, 2);
        gsl_matrix_get_row(Uai,Ua,i);
        gsl_matrix_get_row(Uci,Uc,i);
        gsl_vector_add(Uai, Uci);
        gsl_blas_dgemv(CblasNoTrans, 1.0, Xi, Beta, 1.0, Uai);
        gsl_vector_sub(Yi, Uai);
        gsl_vector_mul(Yi, Yi);
        b += VecSum(Yi);
      }
      
      b = b/2;
      Sigma_E = 1/((Rcpp::rgamma(1,a,1/b))[0]);
      MCMC(GIB,2) = Sigma_E;
      
      gsl_matrix_set_zero(Sigma_beta);
      gsl_vector_set_zero(Mu_beta);
      
      
      for(int i = 0; i < n; ++i)
      {
        gsl_matrix_get_row(Yi,Y,i);
        RowToMatrix(Xi, X, i, 2);
        gsl_matrix_get_row(Uai,Ua,i);
        gsl_matrix_get_row(Uci,Uc,i);
        gsl_vector_add(Uai, Uci);
        gsl_vector_sub(Yi,Uai);
        
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Xi, Xi, 0.0, tmpXi);
        gsl_matrix_add(Sigma_beta, tmpXi);
        
        gsl_blas_dgemv(CblasTrans, 1.0, Xi, Yi, 0.0, tmpx3);
        gsl_vector_add(Mu_beta, tmpx3);
      }

      
      GInverse(Sigma_beta);
      gsl_matrix_scale(Sigma_beta, Sigma_E);
      gsl_blas_dgemv(CblasNoTrans, 1/Sigma_E, Sigma_beta, Mu_beta, 0.0, tmpx3);
      mvrnorm_cpp(tmpx3, Sigma_beta, Beta);
      
      MCMC(GIB,3) = gsl_vector_get(Beta,0);
      MCMC(GIB,4) = gsl_vector_get(Beta,1);
    }
    
    
    for(int i = 0; i < 5; ++i)
    {
      Rcpp::NumericVector tmpEst(MCMC_times-300);
      for(int j = 300; j < MCMC_times; ++j)
      {
        tmpEst[j - 300] = MCMC(j,i);
      }
      Est(CIR, i) = Rcpp::mean(tmpEst);
      SEst(CIR, i) = Rcpp::sd(tmpEst);
    }
  }

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("Est")=Est,
                                      Rcpp::Named("SEst")=SEst
  );
  
  Yi.free();
  Xi.free();
  Uai.free();
  Uci.free();
  tmpXi.free();
  tmpx3.free();  
  tmpx4.free();
  Ue.free();
  InvADSigma_A.free();
  Sigma_Emat.free();
  Y.free();
  X.free();
  A.free();
  invA.free();
  Ua.free();
  Uc.free();
  Sigma_beta.free();
  Mu_beta.free();
  
  return res;
}