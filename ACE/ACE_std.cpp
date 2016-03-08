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


Matrix cMat(double *a, int nrow, int ncol)
{
  gsl_matrix_view A = gsl_matrix_view_array(a, nrow, ncol);
  matrix<double> mat(nrow,ncol);
  gsl_matrix_memcpy(mat, &A.matrix);
  return mat;
}

Vector cVec(double *a, int n)
{
  gsl_vector_view A = gsl_vector_view_array(a, n);
  vector<double> vec(n);
  gsl_vector_memcpy(vec,&A.vector);
  return vec;
}

//  [[Rcpp::export]]
Matrix AB(const double a, const matrix<double> &A, const matrix<double> &B)
{
  int Crow = A.nrow();
  int Ccol = B.ncol();
  matrix<double> C(Crow, Ccol);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, a, A, B, 0.0, C);
  return C;
}

//  [[Rcpp::export]]
double txAx(const Vector &x, const Matrix &A)
{
  int n = x.size();
  matrix<double> vec_t(1,n);
  matrix<double> vec(n,1);
  
  gsl_matrix_set_row(vec_t,0,x);
  gsl_matrix_set_col(vec,0,x);
  
  matrix<double> c1 = AB(1.0,vec_t,A);
  matrix<double> c2 = AB(1.0,c1,vec);
  
  double result = gsl_matrix_get(c2,0,0);
  
  vec_t.free();
  vec.free();
  c1.free();
  c2.free();
  
  return result;
}

// [[Rcpp::export]]
double VecSum(const Vector & x)
{
  double sum = 0;
  int n = x.size();
  vector<double> vec(n);
  gsl_vector_memcpy(vec, x);
  for(int i = 0; i < n; i++)
    sum+=gsl_vector_get(vec, i);
  return sum;
}





//  [[Rcpp::export]]
RcppGSL::Matrix mvrnorm_cpp(int times, const RcppGSL::Vector &mu, const RcppGSL::Matrix &cov)
{
  gsl_set_error_handler_off();
  RcppGSL::matrix<double> chol_cov(cov.nrow(),cov.ncol());
  gsl_matrix_memcpy(chol_cov, cov);
  if(gsl_linalg_cholesky_decomp(chol_cov)) 
  {
    printf("Non def-pos!\n");
    RcppGSL::matrix<double> errMat(1,1);
    return errMat;
  }
  else
  {
    int n = mu.size();
    Rcpp::NumericVector rdvec = Rcpp::rnorm(n * times,0,1);
    RcppGSL::matrix<double> mat(times,n);
    RcppGSL::matrix<double> matMu(times,n);
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
    mat.free();
    chol_cov.free();
    return matMu;
  }
}




//  [[Rcpp::export]]
Matrix GetInverse(matrix<double> A)  
{  
  int n = A.ncol();  
  matrix<double> tmpA(n, n);
  matrix<double> inverse(n, n);
  gsl_matrix_memcpy(tmpA, A);  
  gsl_permutation *p = gsl_permutation_alloc(n);  
  int sign = 0;  
  gsl_linalg_LU_decomp(tmpA, p, &sign);  
  gsl_linalg_LU_invert(tmpA, p, inverse);  
  gsl_permutation_free(p);  
  tmpA.free(); 
  return inverse;
}


//  [[Rcpp::export]]
Matrix RowToMatrix(matrix<double> X, int row, int k)
{
  int n = X.ncol();
  matrix<double> tmpX(n/k, k);
  for(int i = 0; i < n/k; ++i)
    for(int j = 0; j < k; ++j)
      tmpX(j, i) = gsl_matrix_get(X,row, i * n/k + j); //use double(mat(i,j)) or gsl_matrix_get(mat,i,j)
  return tmpX;  
}

//  [[Rcpp::export]] 
Vector GetRow(matrix<double> X, int row)
{
  int n = X.ncol();
  vector<double> vec(n);
  gsl_matrix_get_row(vec,X,row);
  return vec;
}

//  [[Rcpp::export]]
Vector GetCol(matrix<double> X, int col)
{
  int n = X.nrow();
  vector<double> vec(n);
  gsl_matrix_get_col(vec,X,col);
  return vec;
}

//  [[Rcpp::export]]
vector<double> add(vector<double> X, vector<double> Y)
{
  int n = Y.size();
  vector<double> vec(n);
  gsl_vector_memcpy(vec, Y);
  gsl_vector_add(vec, X);
  return vec;
}


//  [[Rcpp::export]]
vector<double> Axpy(matrix<double> A, vector<double> x, vector<double> y)
{
  int n = y.size();
  vector<double> vec(n);
  gsl_vector_memcpy(vec, y);
  gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 1.0, vec);
  return(vec);
}

vector<double> aAx(double a, matrix<double> A, vector<double> x)
{
  vector<double> vec(A.nrow());
  gsl_blas_dgemv(CblasNoTrans, a, A, x, 0.0, vec);
  return(vec);
}

//  [[Rcpp::export]]
matrix<double> aX(double a, matrix<double> X)
{
  matrix<double> mat(X.nrow(),X.ncol());
  gsl_matrix_memcpy(mat,X);
  gsl_matrix_scale(mat,a);
  return mat;
}














//  [[Rcpp::export]]
extern "C" SEXP test()
{
  //Initialization
  int n = 500, n1 = 250, p = 2;
  matrix<double> Y(n, p);
  matrix<double> X(n, 2 * p); // the array X(n,2,p) in ACE.r
  double tmpA[] = {1.0,0.5,0.5,1.0};
  matrix<double> A = cMat(tmpA,2,2);
  
  vector<double> Beta(2);
  double Sigma_A = 0.4;
  double Sigma_C = 0.3;
  double Sigma_E = 0.2;
  
  matrix<double> invA(GetInverse(A));
  matrix<double> Ua(n,2), Uc(n,2);
  
  int sp = 10;
  int MCMC_times = 500;
  Rcpp::NumericMatrix Est(sp,5), SEst(sp,5);
  gsl_matrix_set_zero(Ua);
  gsl_matrix_set_zero(Uc);

  matrix<double> Sigma_beta(p,p);
  vector<double> Mu_beta(p);
  //  matrix<double> Mu_beta(p,1);
  //  matrix<double> Mu_beta_tmp(p,1);
  
  double tmpmu[] = {0,0};
  vector<double> tmpMu(cVec(tmpmu,2));
  
  Rcpp::NumericMatrix MCMC(MCMC_times, 5);
  
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
      
      vector<double> Ue(Rcpp::rnorm(2,0.0,std::sqrt(Sigma_E)));
      
      double tmpx1 = (Rcpp::runif(1,5,13))[0];
      Rcpp::NumericVector tmpx2 = Rcpp::rnorm(2,0.0,2.0);
      for(int m = 0; m < 2; ++m) 
      {
        X(i,m) = tmpx1;
        X(i,m+2) = tmpx2[m];
      }
      
      matrix<double> Xi = RowToMatrix(X,i,2);
      
      
      vector<double> Uai = GetRow(Ua,i);
      vector<double> Uci = GetRow(Uc,i);
      vector<double> Ui = add( add( Uai, Uci), Ue);
      vector<double> Yi = Axpy(Xi,Beta,Ui);
      gsl_matrix_set_row(Y,i,Yi);
      
      //free ... why it crashed when free?
      Xi.free();
      Uai.free();
      Uci.free();
      Ui.free();
      Yi.free();
      Ue.free();
    }
    
    matrix<double> SigmaAA = aX(Sigma_A,A);
    matrix<double> UaMvrnorm = mvrnorm_cpp(n-n1,tmpMu,SigmaAA);
    Rcpp::NumericVector UcRnorm = Rcpp::rnorm(n-n1,0.0,std::sqrt(Sigma_C));
    Rcpp::NumericVector UeRnorm = Rcpp::rnorm(2*(n-n1), 0.0, std::sqrt(Sigma_E));
    Rcpp::NumericVector XRunif = Rcpp::runif(n-n1,5,13);
    Rcpp::NumericVector XRnorm = Rcpp::rnorm(2*(n-n1),0.0,2.0);
    
    SigmaAA.free();
    
    for(int i = n1; i < n; ++i)
    {
      gsl_matrix_set_row(Ua, i, GetRow(UaMvrnorm, i-n1));
      Uc(i,1) = UcRnorm[i-n1];
      Uc(i,0) = gsl_matrix_get(Uc,i,1);
      
      vector<double> Ue(2);
      Ue[0] = UeRnorm[2*(i-n1)]; Ue[1] = UeRnorm[2*(i-n1)+1];
      
      for(int m = 0; m < 2; ++m)
      {
        X(i,m) = XRunif[i-n1];
        X(i,m+2) = XRnorm[2*(i-n1)+m];
      }
      
      matrix<double> Xi = RowToMatrix(X,i,2);
      

      vector<double> Uai = GetRow(Ua,i);
      vector<double> Uci = GetRow(Uc,i);
      vector<double> Ui = add( add( Uai, Uci), Ue);
      vector<double> Yi = Axpy(Xi,Beta,Ui);
      gsl_matrix_set_row(Y,i,Yi);
      
      Xi.free();
      Uai.free();
      Uci.free();
      Ui.free();
      Yi.free();
      Ue.free();
      
    } 
    
    UaMvrnorm.free();
    

    Beta[0] = 0.0; Beta[1] = 0.0;
    Sigma_A = Sigma_C = Sigma_E = 1.0;
    
    
    
    for(int i = 0; i < n1; ++i)
    {
      Ua(i, 0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_A)))[0];
      Ua(i, 1) = gsl_matrix_get(Ua,i,0);
      Uc(i, 0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_C)))[0];
      Uc(i, 1) = gsl_matrix_get(Uc,i,0);
    }
    
    for(int i = n1; i < n; ++i)
    {
      gsl_matrix_set_row(Ua,i,GetRow(mvrnorm_cpp(1,tmpMu,aX(Sigma_A,A)),0));
      Uc(i, 0) = (Rcpp::rnorm(1, 0.0, std::sqrt(Sigma_C)))[0];
      Uc(i, 1) = gsl_matrix_get(Uc, i, 0);
    }
    
    for(int GIB = 0; GIB < MCMC_times; ++GIB)
    {
      double Sigma_A2 = 1 / (2/Sigma_E + 1/Sigma_A);
      
      for(int i = 0; i < n1; ++i)
      {
        vector<double> Yi = GetRow(Y,i);
        matrix<double> Xi = RowToMatrix(X,i,2);
        vector<double> Uci = GetRow(Uc,i);
        gsl_vector_sub(Yi,Axpy(Xi,Beta,Uci));
        double Mu_A = VecSum(Yi)*Sigma_A2/Sigma_E;
        Ua(i, 0) = (Rcpp::rnorm(1, Mu_A, std::sqrt(Sigma_A2)))[0];
        Ua(i, 1) = gsl_matrix_get(Ua,i,0);
        
        //free
        Yi.free();
        Xi.free();
        Uci.free();
        
      }
      
      matrix<double> Sigma_Emat(2, 2); ///////////
      gsl_matrix_set_identity(Sigma_Emat);
      gsl_matrix_scale(Sigma_Emat, 1/Sigma_E);
      
      matrix<double> InvADSigma_A(2, 2);
      gsl_matrix_memcpy(InvADSigma_A, invA);
      
      gsl_matrix_scale(InvADSigma_A, 1/Sigma_A);
      gsl_matrix_add(Sigma_Emat, InvADSigma_A);
      matrix<double> Sigma_A2mat = GetInverse(Sigma_Emat);
      
      InvADSigma_A.free();
      Sigma_Emat.free();
      
      
      
      for(int i = n1; i < n; ++i)
      {
        vector<double> Yi = GetRow(Y,i);
        matrix<double> Xi = RowToMatrix(X,i,2);
        vector<double> Uci = GetRow(Uc,i);
        gsl_vector_sub(Yi,Axpy(Xi,Beta,Uci));
        vector<double> Mu_Amat = aAx(1/Sigma_E, Sigma_A2mat, Yi);
        gsl_matrix_set_row(Ua,i,GetRow(mvrnorm_cpp(1, Mu_Amat, Sigma_A2mat),0));
        
        Yi.free();
        Xi.free();
        Uci.free();
        Mu_Amat.free();
      }
      
      double Sigma_C2 = 1/(2/Sigma_E+1/Sigma_C);
      for(int i = 0; i < n; ++i)
      {
        vector<double> Yi = GetRow(Y,i);
        matrix<double> Xi = RowToMatrix(X,i,2);
        vector<double> Uai = GetRow(Ua,i);
        gsl_vector_sub(Yi,Axpy(Xi,Beta,Uai));
        double Mu_C = VecSum(Yi)*Sigma_C2/Sigma_E;
        Uc(i, 0) = (Rcpp::rnorm(1, Mu_C, std::sqrt(Sigma_C2)))[0];
        Uc(i, 1) = gsl_matrix_get(Uc, i, 0);
        
        Yi.free();
        Xi.free();
        Uai.free();
      }
      
      double a = (n1+2*(n-n1)+1)/2, b = 0.0;
      //double tmpb = 0;
      for(int i = 0; i < n1; ++i)
      {
        double tUa = gsl_matrix_get(Ua, i, 0);
        tUa *= tUa;
        b += tUa; 
      }
      
      vector<double> tmpUaVec(2);
      
      for(int i = n1; i < n; ++i)
      {
        gsl_matrix_get_row(tmpUaVec, Ua, i);
        b += txAx(tmpUaVec, invA);
      }
      
      tmpUaVec.free();
      
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
      
//       //以下变量应该可以和上面那几个循环一起，
//       //都改为全局变量来节省释放内存的次数
//       //vector<double> Yi(2);
//       matrix<double> Xi(2,2);
//       vector<double> Uai(2);
//       vector<double> Uci(2);
//       vector<double> UaUc(2);
//       matrix<double> tXi(2,2);
//       //嗯，以上变量有空改，不然老是内存泄漏呀~
      
      for(int i = 0; i < n; ++i)
      {
        vector<double> Yi = GetRow(Y,i);
        //gsl_matrix_get_row(Yi,Y,i);
        matrix<double> Xi = RowToMatrix(X,i,2);
        //gsl_matrix_memcpy(Xi, RowToMatrix(X,i,2));
        vector<double> Uai = GetRow(Ua, i);
        //gsl_matrix_get_row(Uai,Ua,i);
        vector<double> Uci = GetRow(Uc, i);
        //gsl_matrix_get_row(Uci,Uc,i);
        gsl_vector_add(Uai, Uci);
        vector<double> UaUc = Axpy(Xi,Beta,Uai);
        //gsl_vector_memcpy(UaUc, Axpy(Xi,Beta,Uai));
        gsl_vector_sub(Yi,UaUc);
        gsl_vector_mul(Yi,Yi);
        b += VecSum(Yi);
        
        Yi.free();
        Xi.free();
        Uai.free();
        Uci.free();
        UaUc.free();

      }
      
      b = b/2;
      Sigma_E = 1/((Rcpp::rgamma(1,a,1/b))[0]);
      MCMC(GIB,2) = Sigma_E;
      
      gsl_matrix_set_zero(Sigma_beta);
      gsl_vector_set_zero(Mu_beta);
      
      
      for(int i = 0; i < n; ++i)
      {
        vector<double> Yi = GetRow(Y, i);
        // gsl_matrix_get_row(Yi,Y,i);
        matrix<double> Xi = RowToMatrix(X,i,2);
        //gsl_matrix_memcpy(Xi, RowToMatrix(X,i,2));
        vector<double> Uai = GetRow(Ua,i);
        //gsl_matrix_get_row(Uai,Ua,i);
        vector<double> Uci = GetRow(Uc,i);
        //gsl_matrix_get_row(Uci,Uc,i);
        gsl_vector_add(Uai, Uci);
        gsl_vector_sub(Yi,Uai);
        
        matrix<double> tXi(2,2);
        gsl_matrix_transpose_memcpy(tXi, Xi);
        
        matrix<double> ABtXiXi = AB(1,tXi,Xi);
        //gsl_matrix_memcpy(Xi, AB(1,tXi,Xi));
        gsl_matrix_add(Sigma_beta, ABtXiXi);
        
        vector<double> aAxtXiYi = aAx(1,tXi,Yi);
        //gsl_vector_memcpy(Yi, aAx(1,tXi,Yi));
        gsl_vector_add(Mu_beta, aAxtXiYi);
        
        Yi.free();
        Xi.free();
        Uai.free();
        Uci.free();
        tXi.free();
        ABtXiXi.free();
        aAxtXiYi.free();
        
      }
      
      gsl_matrix_memcpy(Sigma_beta, aX(Sigma_E, GetInverse(Sigma_beta)));
      gsl_vector_memcpy(Mu_beta, aAx(1/Sigma_E,Sigma_beta, Mu_beta));
      gsl_vector_memcpy(Beta, GetRow(mvrnorm_cpp(1,Mu_beta,Sigma_beta),0));
      
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
  
  Rcpp::NumericMatrix xx(X.nrow(),X.ncol());
  Rcpp::NumericMatrix yy(Y.nrow(),Y.ncol());

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("Est")=Est,
                                     Rcpp::Named("SEst")=SEst
                                    );
  
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