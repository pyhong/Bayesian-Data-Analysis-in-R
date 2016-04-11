//  [[Rcpp::depends(RcppGSL)]]
#include <cstring>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <limits>
#include <time.h>
#include <fstream>

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

using namespace std;

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
    return (rand() / (RAND_MAX + 1.0));
  }
  inline double itaillimit( double b)
  {
    return EXPTLOVTL - exp(SQ_TAIL_LIMIT - TAIL_LIMIT * b - LOG_TAIL_LIMIT );
  }

  inline double sample_tail( double a, double b)
  {
    double sa =  1 ;
    double sb = exp(a * (a - b));
    if ( sa == sb )
      return (a * a < b * b ? a : b);

    double u;
    do {
      u = rnd() * (sb - sa) + sa;
    } while ( u == 0 );
    u = a - log(u) / a;

    double v = rnd() * exp(0.5 * (a - u) * (a - u));
    if ( v <= 1 )
      return u;
    else
      return NAN;
  }
  double truncate2(double a, double b)
  {
    double z;
    int y, ia, ib;
    double ah, bh;
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
      ia = whichbin[ (int)(ax * 128.0) ];
      if (xhist[ia + 1] <= ax )
        ++ia;
      ia = ((a >= 0) ? ia : -ia - 1);
    }

    double p_tail_right = 0;
    if (b <= -TAIL_LIMIT ) {
      do {
        z = -sample_tail(-b, -a);
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
      ib = whichbin[ (int)(bx * 128.0) ];
      if (xhist[ib + 1] <= bx )
        ++ib;
      ib = ((b >= 0) ? ib : -ib - 1);
    }

    if (ia == ib ) {
      do {
        z = rnd() * (bh - ah) + ah;
        y = rand();
      } while ( y >= yhistratio[ia] && y * yhist[ia] >= exp(-0.5 * z * z) * (RAND_MAX + 1.0) );
      return (z);
    }

    double middlep = (ib - ia - 1) * basep;
    do {
      double u = rnd() * (middlep + ap + bp + p_tail_left + p_tail_right );
      if ( u < middlep ) {
        u /= basep;
        int ui = (int)u;
        u = modf(u, &z);
        int col = ui + ia  + 1;
        z = u * dxhist[col] + xhist[col];
        y = rand();

        if ( (y >= yhistratio[col] ) && y * yhist[col] >= exp(-0.5 * z * z) * (RAND_MAX + 1.0) )
          z = NAN;
        continue;
      }
      if (ap == basep )
        ap = (xhist[ia + 1] - ah) * yhist[ia];
      if (u < middlep + ap ) {
        z = (u - middlep) / ap * (xhist[ia + 1] - ah) + ah;
        y = rand();
        if ( y >= yhistratio[ia] && y * yhist[ia] >= exp(-0.5 * z * z) * (RAND_MAX + 1.0) )
          z = NAN;
        continue;
      }
      if (bp == basep)
        bp = (bh - xhist[ib]) * yhist[ib];
      if ( u < middlep + ap + bp) {
        z = (u - middlep - ap) / bp * (bh - xhist[ib]) + xhist[ib];
        y = rand();
        if ( y >= yhistratio[ib] && y * yhist[ib] >= exp(-0.5 * z * z) * (RAND_MAX + 1.0) )
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

    double a = (aa - mu) / sigma;
    double b = (bb - mu) / sigma;

    if ( a == -infinity && b == infinity )
      return NAN;
    else
      return  sigma * truncate2(a, b) + mu;
  }

  // creates the object by initializing the arrays we will need
  truncated()
  {

    rr = gsl_rng_alloc( gsl_rng_rand );
    gsl_rng_set(rr, time(NULL));
    double x = 0;
    int i = 0, c = 0;
    while ( x < 3.5 ) {
      x += basep * exp(0.5 * x * x);
      ++c;
    }
    xhist = (double*) malloc((2 * c + 1) * sizeof(double)) + c;
    dxhist = (double*) malloc(2 * c * sizeof(double)) + c;
    yhist = (double*)malloc(2 * (c + 1) * sizeof(double)) + (c + 1);

    yhistratio = (int*)malloc(2 * c * sizeof(int)) + c;
    whichbin = (int*)malloc( ((int)(TAIL_LIMIT * 128.0 ) + 1) * sizeof(int));

    x = 0; xhist[0] = 0; yhist[0] = 1; yhist[-1] = 1;

    for ( i = 1; i <= c; ++i )  {
      x += basep * exp(0.5 * x * x);
      xhist[i] = x;
      xhist[-i] = -x;
      dxhist[i - 1] = xhist[i] - xhist[i - 1];
      dxhist[-i] = dxhist[i - 1];
      yhist[i] = exp(-0.5 * x * x);
      yhist[-i - 1] = yhist[i];
      yhistratio[i - 1] = (int)ceil( yhist[i] / yhist[i - 1] * (RAND_MAX + 1.0) );
      yhistratio[-i] = yhistratio[i - 1];
    }
    c = 0;
    for ( i = 0; i <= (int)(TAIL_LIMIT / basep); ++i ) {
      while ( (int)(xhist[c] / basep ) < i )
        ++c;
      whichbin[i] = c - 1;
    }
  }
};
//rtruncnorm

int sp_int(gsl_rng *rr, int n, double * prob)
{
  double cdf = 0;
  double sum = 0;
  double x = gsl_rng_uniform(rr);
  for (int i = 0; i < n; ++i) sum += prob[i];
  if (sum > 0)
  {
    for (int i = 0; i < n; ++i)
    {
      cdf += prob[i] / sum;
      if (x < cdf) return i;
    }
  }
  else return -1;
}

double txAx(const gsl_vector *x, const gsl_matrix *A)
{
  int n = x->size;
  gsl_matrix * vec_t = gsl_matrix_calloc(1, n);
  gsl_matrix * tt = gsl_matrix_calloc(1, n);
  gsl_matrix * vec = gsl_matrix_calloc(n, 1);
  gsl_matrix * C = gsl_matrix_calloc(1, 1);
  gsl_matrix_set_row(vec_t, 0, x);
  gsl_matrix_set_col(vec, 0, x);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, vec_t, A, 0.0, tt);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tt, vec, 0.0, C);
  double result = gsl_matrix_get(C, 0, 0);
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
  for (int i = 0; i < n; ++i)
    sum += gsl_vector_get(x, i);
  return sum;
}


void RowToMatrix(gsl_matrix * tmpX, const gsl_matrix * X, int row, int k)
{
  int n = X->size2;
  for (int i = 0; i < n / k; ++i)
    for (int j = 0; j < k; ++j)
    {
      double tt = gsl_matrix_get(X, row, i * n / k + j); //use double(mat(i,j)) or gsl_matrix_get(mat,i,j)
      gsl_matrix_set(tmpX, j, i, tt);
    }
}

void mvrnorm_cpp(const gsl_rng *rr, const gsl_vector * mu, const gsl_matrix * cov, gsl_vector * res)
{
  int times = 1;
  gsl_set_error_handler_off();
  gsl_matrix * chol_cov = gsl_matrix_calloc(cov->size1, cov->size2);
  gsl_matrix_memcpy(chol_cov, cov);
  if (gsl_linalg_cholesky_decomp(chol_cov))
  {
    printf("Non def-pos!\n");
  }
  else
  {
    int n = mu->size;
    gsl_matrix * matMu = gsl_matrix_calloc(1, n);
    gsl_matrix * mat = gsl_matrix_calloc(times, n);
    int count = 0;
    for (int i = 0; i < times; ++i)
      for (int j = 0; j < n; ++j)
      {
        gsl_matrix_set(mat, i, j, (gsl_ran_gaussian(rr, 1)));
        count += 1;
        double muj = gsl_vector_get(mu, j);
        gsl_matrix_set(matMu, i, j, muj);
      }

    gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, chol_cov, mat);
    gsl_matrix_add(matMu, mat);
    gsl_matrix_get_row(res, matMu, 0);
    gsl_matrix_free(matMu);
    gsl_matrix_free(mat);
    gsl_matrix_free(chol_cov);
  }
}

void GInverse(gsl_matrix * A)
{
  int n = A->size2;
  gsl_matrix * inverse = gsl_matrix_calloc(n, n);
  gsl_permutation *p = gsl_permutation_alloc(n);
  int sign = 0;
  gsl_linalg_LU_decomp(A, p, &sign);
  gsl_linalg_LU_invert(A, p, inverse);
  gsl_permutation_free(p);
  gsl_matrix_memcpy(A, inverse);
  gsl_matrix_free(inverse);
}

void printf_matrix(const gsl_matrix * M)
{
  int nrow = M->size1;
  int ncol = M->size2;
  for (int i = 0; i < nrow; ++i)
  {
    for (int j = 0; j < ncol; ++j)
      printf("%f ", gsl_matrix_get(M, i, j));
    printf("\n");
  }
}

void printf_vector(const gsl_vector * x)
{
  int size = x->size;
  for (int i = 0; i < size; ++i)
    printf("%f ", gsl_vector_get(x, i));
  printf("\n");
}

void colMeanSd(const gsl_matrix *A, gsl_vector *x, gsl_vector *y)
{
  int ncol = A->size2;
  int nrow = A->size1;
  double mean = 0.0, sd = 1.0;
  gsl_vector * A_col = gsl_vector_alloc(nrow);
  for (int i = 0; i < ncol; ++i)
  {
    gsl_matrix_get_col(A_col, A, i);
    mean = gsl_stats_mean(A_col->data, 1, nrow);
    sd = gsl_stats_sd(A_col->data, 1, nrow);
    gsl_vector_set(x, i, mean);
    gsl_vector_set(y, i, sd);
  }
}


int main()
{
  const gsl_rng_type * T;
  gsl_rng * rr;
  gsl_rng_env_setup();
  gsl_rng_default_seed = 0;
  T = gsl_rng_default;
  rr = gsl_rng_alloc(T);
  gsl_rng_set(rr, time(NULL));

  int n = 500, n1 = 250, p = 2, MCAX = 1000, GNUM = 800;
  gsl_matrix * Y = gsl_matrix_calloc(n, p);
  gsl_matrix * Ystar = gsl_matrix_calloc(n, p);
  gsl_matrix * X = gsl_matrix_calloc(n, 2 * p);
  double tmpA[] = {1.0, 0.5, 0.5, 1.0};
  gsl_matrix_view tmpA_view = gsl_matrix_view_array(tmpA, 2, 2);
  gsl_matrix * A = gsl_matrix_calloc(2, 2);
  gsl_matrix_memcpy(A, &tmpA_view.matrix);

  gsl_vector * Beta = gsl_vector_calloc(2);
  gsl_vector_set(Beta, 0, 1.0);
  gsl_vector_set(Beta, 1, 0.5);
  double Sigma_A = 0.5;
  double Sigma_C = 0.3;
  double Sigma_E = 0.1;

  gsl_matrix * invA = gsl_matrix_calloc(2, 2);
  gsl_matrix_memcpy(invA, A);
  GInverse(invA);
  gsl_matrix * Ua = gsl_matrix_calloc(n, 2);
  gsl_matrix * Uc = gsl_matrix_calloc(n, 2);

  gsl_vector * Est = gsl_vector_alloc(5);
  gsl_vector * SEst = gsl_vector_alloc(5);


  gsl_matrix * Sigma_beta = gsl_matrix_calloc(p, p);
  gsl_vector * Mu_beta = gsl_vector_calloc(p);//

  double tmpmu[] = {0, 0};
  gsl_vector_view tmpmu_view = gsl_vector_view_array(tmpmu, 2);
  gsl_vector * tmpMu = gsl_vector_alloc(2);
  gsl_vector_memcpy(tmpMu, &tmpmu_view.vector);

  gsl_matrix * MCMC = gsl_matrix_calloc(MCAX, 5);

  gsl_vector * Yi = gsl_vector_calloc(2);
  gsl_matrix * Xi = gsl_matrix_calloc(2, 2);
  gsl_vector * Uai = gsl_vector_calloc(2);
  gsl_vector * Uci = gsl_vector_calloc(2);
  gsl_matrix * tmpXi = gsl_matrix_calloc(2, 2);
  gsl_vector * tmpx3 = gsl_vector_calloc(2);
  gsl_vector * tmpx4 = gsl_vector_calloc(2);
  gsl_vector * Ue = gsl_vector_calloc(2);
  gsl_matrix * SigmaAA = gsl_matrix_calloc(2, 2);
  gsl_matrix * Sigma_Emat = gsl_matrix_calloc(2, 2);
  gsl_matrix * InvADSigma_A = gsl_matrix_calloc(2, 2);



  int kk = 25;  // the number of basis
  // thus nbreaks = 25+2-4=23
  //
  gsl_matrix * D0 = gsl_matrix_calloc(kk - 1, kk);
  gsl_matrix * D1 = gsl_matrix_calloc(kk - 1 - 1, kk - 1);
  double a[2] = { -1, 1};

  gsl_matrix * N = gsl_matrix_calloc(kk - 2, kk);
  gsl_matrix * M = gsl_matrix_calloc(kk, kk);
  gsl_matrix * tmp_M = gsl_matrix_calloc(kk, kk);
  gsl_vector * unit_e = gsl_vector_calloc(kk);
  gsl_vector * tmpMu_e = gsl_vector_calloc(kk);

  gsl_matrix * Bk[2];
  gsl_matrix * bk[2];
  for (int i = 0; i < 2; ++i)
  {
    Bk[i] = gsl_matrix_calloc(kk, n);
    bk[i] = gsl_matrix_calloc(kk, n);
  }

  gsl_matrix * tmpBk = gsl_matrix_calloc(kk, kk);

  gsl_matrix * ks = gsl_matrix_calloc(2, kk + 4 - 2);
  gsl_vector * Ycol = gsl_vector_calloc(n);
  gsl_vector * Ycol_tmp = gsl_vector_calloc(n);

  int bsorder = 4;
  int nbreak = kk + 2 - bsorder;
  gsl_bspline_workspace *bw = gsl_bspline_alloc(bsorder, nbreak);
  gsl_vector * ksRow = gsl_vector_calloc(nbreak);

  gsl_vector * Bcoef = gsl_vector_calloc(kk);
  int nderiv = 1;
  gsl_matrix * Bcoef_Driv = gsl_matrix_calloc(kk, nderiv + 1);

  gsl_matrix * Gam = gsl_matrix_calloc(kk, 2);
  double tao[2];
  double sigma_r[] = {2.7, 2.7};
  gsl_vector * p_accept = gsl_vector_calloc(2);

  gsl_matrix * f = gsl_matrix_calloc(2, 10); //an Aid Matrix?
  gsl_vector * f_row = gsl_vector_calloc(f->size2);

  int BK_ncol = 100;

  gsl_matrix * BK[2];
  gsl_matrix * GY[2];
  for (int i = 0; i < 2; ++i)
  {
    BK[i] = gsl_matrix_calloc(kk, BK_ncol);
    GY[i] = gsl_matrix_calloc(BK_ncol, (MCAX - GNUM));
  }
  gsl_vector * tmpBK = gsl_vector_calloc(BK_ncol);

  gsl_matrix * gy = gsl_matrix_calloc(n, p);

  int sp = 10;
  gsl_matrix * com_gy[sp];
  gsl_matrix * com_gy_beta[sp];
  for (int i = 0; i < sp; ++i)
  {
    com_gy[i] = gsl_matrix_calloc(n, 2);
    com_gy_beta[i] = gsl_matrix_calloc(n, 2);
  }

  gsl_matrix * y_beta = gsl_matrix_calloc(n, p);
  double lower, upper;
  gsl_matrix * lo_up = gsl_matrix_calloc(2, kk);
  gsl_vector * lo_up_row = gsl_vector_calloc(kk);


  gsl_matrix * w1 = gsl_matrix_calloc(kk, 10);
  gsl_matrix * w2 = gsl_matrix_calloc(kk, 10);
  gsl_vector * wstar = gsl_vector_calloc(kk);
  gsl_vector * tmp_unit_e = gsl_vector_calloc(kk);
  gsl_vector * Gam_col = gsl_vector_calloc(kk);
  gsl_vector * w1_ff = gsl_vector_calloc(kk);
  gsl_vector * Bk__ij = gsl_vector_calloc(kk);
  gsl_vector * bk__ij = gsl_vector_calloc(kk);
  gsl_vector * f1_row = gsl_vector_calloc(10);
  gsl_vector * f2_row = gsl_vector_calloc(10);

  double Accept = 0;

  truncated truncnorm;
  double r1[10], r2[10];

  int t1_tmp[10];
  for (int ff = 0; ff < 10; ++ff)
    t1_tmp[ff] = ff;

  double Uai0, Uci0;

  for (int i = 0; i < n1; ++i)
  {
    for (int k = 0; k < 2; ++k)
      gsl_vector_set(Ue, k, gsl_ran_gaussian(rr, sqrt(Sigma_E)));

    Uai0 = gsl_ran_gaussian(rr, sqrt(Sigma_A));
    Uci0 = gsl_ran_gaussian(rr, sqrt(Sigma_C));

    for (int m = 0; m < 2; ++m)
    {
      gsl_matrix_set(X, i, m, gsl_ran_gaussian(rr, 2.0));
      gsl_matrix_set(X, i, m + 2, gsl_ran_gaussian(rr, 2.0));
      gsl_matrix_set(Ua, i, m, Uai0);
      gsl_matrix_set(Uc, i, m, Uci0);
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

  for (int i = n1; i < n; ++i)
  {
    mvrnorm_cpp(rr, tmpMu, SigmaAA, tmpx3);
    gsl_matrix_set_row(Ua, i, tmpx3);
    Uci0 = gsl_ran_gaussian(rr, sqrt(Sigma_C));

    for (int k = 0; k < 2; ++k)
      gsl_vector_set(Ue, k, gsl_ran_gaussian(rr, sqrt(Sigma_E)));

    for (int m = 0; m < 2; ++m)
    {
      gsl_matrix_set(X, i, m, gsl_ran_gaussian(rr, 2.0));
      gsl_matrix_set(X, i, m + 2, gsl_ran_gaussian(rr, 2.0));
      gsl_matrix_set(Uc, i, m, Uci0);
    }

    RowToMatrix(Xi, X, i, 2);
    gsl_matrix_get_row(Uai, Ua, i);
    gsl_matrix_get_row(Uci, Uc, i);
    gsl_vector_add(Uai, Uci);
    gsl_vector_add(Uai, Ue);
    gsl_blas_dgemv(CblasNoTrans, 1.0, Xi, Beta, 1.0, Uai);
    gsl_matrix_set_row(Y, i, Uai);
  }

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      gsl_matrix_set(Ystar, i, j, exp(gsl_matrix_get(Y, i, j)));

  for (int i = 0; i < (kk - 1); ++i)
  {
    gsl_matrix_set(D0, i, i, a[0]);
    gsl_matrix_set(D0, i, i + 1, a[1]);
  }

  for (int i = 0; i < (kk - 2); ++i)
  {
    gsl_matrix_set(D1, i, i, a[0]);
    gsl_matrix_set(D1, i, i + 1, a[1]);
  }


  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, D1, D0, 0, N);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, N, N, 0, M);

  for (int i = 0; i < 2; ++i)
  {
    gsl_matrix_get_col(Ycol, Ystar, i);
    double minYi = gsl_vector_min(Ycol);
    double maxYi = gsl_vector_max(Ycol);
    gsl_bspline_knots_uniform(minYi, maxYi, bw);
    gsl_matrix_set_row(ks, i, bw->knots);
    for (int k = 0; k < n; ++k)
    {
      double yj = gsl_vector_get(Ycol, k);
      gsl_bspline_deriv_eval(yj, nderiv, Bcoef_Driv, bw);   //something's wrong here
      gsl_matrix_get_col(Bcoef, Bcoef_Driv, 0); //Deriv of order 0
      gsl_matrix_set_col(Bk[i], k, Bcoef);
      gsl_matrix_get_col(Bcoef, Bcoef_Driv, 1); //Deriv of order 1
      gsl_matrix_set_col(bk[i], k, Bcoef);

    }
  }

  for (int i = 0; i < kk; ++i)
  {
    gsl_matrix_set(Gam, i, 0, (13.5 - 0.285) / 24 * i + 0.285);
    gsl_matrix_set(Gam, i, 1, (15 - 0.285) / 24 * i + 0.285);
  }

//  printf_matrix(Gam);

  for (int i = 0; i < 2; ++i)
  {
    tao[i] = 1 / (gsl_ran_gamma(rr, 1, 1 / 0.005));
    // cout<<tao[i]<<endl;
  }


  gsl_vector_set_zero(p_accept);
  gsl_matrix_set_zero(f);


  for (int i = 0; i < 2; ++i)
  {
    gsl_matrix_get_col(Ycol, Ystar, i);
    double minYi = gsl_vector_min(Ycol);
    double maxYi = gsl_vector_max(Ycol);
    gsl_bspline_knots_uniform(minYi, maxYi, bw);
    for (int j = 0; j < BK_ncol; ++j)
    {
      double y = (maxYi - minYi) / (BK_ncol - 1.0) * j + minYi;
      gsl_bspline_eval(y, Bcoef, bw);
      gsl_matrix_set_col(BK[i], j, Bcoef);
    }
  }


  for (int i = 0; i < 2; ++i)
  {
    gsl_matrix_set_zero(GY[i]);
  }
  gsl_matrix_set_zero(gy);

  for (int i = 0; i < n1; ++i)
  {
    Uai0 = gsl_ran_gaussian(rr, sqrt(Sigma_A));
    Uci0 = gsl_ran_gaussian(rr, sqrt(Sigma_C));
    for (int m = 0; m < 2; ++m)
    {
      gsl_matrix_set(Ua, i, m, Uai0);
      gsl_matrix_set(Uc, i, m, Uci0);
    }
  }

  gsl_matrix_memcpy(SigmaAA, A);
  gsl_matrix_scale(SigmaAA, Sigma_A);

  for (int i = n1; i < n; ++i)
  {
    mvrnorm_cpp(rr, tmpMu, SigmaAA, tmpx3);
    gsl_matrix_set_row(Ua, i, tmpx3);
    Uci0 = gsl_ran_gaussian(rr, sqrt(Sigma_C));
    for (int m = 0; m < 2; ++m)
      gsl_matrix_set(Uc, i, m, Uci0);
  }

  char filename[] = "gy.txt";
  ofstream fout(filename);

  for (int GIB = 0; GIB < MCAX; ++GIB)
  {
    for (int j = 0; j < 2; ++j)
    {
      gsl_matrix_memcpy(tmp_M, M);

      gsl_matrix_scale(tmp_M, 1 / tao[j]);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1 / Sigma_E, Bk[j], Bk[j], 1, tmp_M);

      GInverse(tmp_M);
      mvrnorm_cpp(rr, tmpMu_e, tmp_M, unit_e);
      gsl_vector_memcpy(tmp_unit_e, unit_e);
      gsl_vector_mul(tmp_unit_e, unit_e);
      gsl_vector_scale(unit_e, 1 / (sqrt(VecSum(tmp_unit_e))));

      for (int i = 0; i < kk; ++i)
      {
        gsl_matrix_set(lo_up, 0, i, -infinity);
        gsl_matrix_set(lo_up, 1, i, infinity);
      }
//      printf_matrix(lo_up);
//      cout<<min_double<<max_double<<endl;
//      printf_vector(unit_e);
      double lo_up_i, Gam_0, Gam_1, unit_e_0, unit_e_1;;
      for (int ii = 0; ii < (kk - 1); ++ii)
      {
        Gam_0 = gsl_matrix_get(Gam, ii, j);
        Gam_1 = gsl_matrix_get(Gam, ii + 1, j);
        unit_e_0 = gsl_vector_get(unit_e, ii);
        unit_e_1 = gsl_vector_get(unit_e, ii + 1);
        lo_up_i = (Gam_0 - Gam_1) / (unit_e_1 - unit_e_0);
        if (unit_e_1 > unit_e_0)
          gsl_matrix_set(lo_up, 0, ii, lo_up_i);
        else
          gsl_matrix_set(lo_up, 1, ii, lo_up_i);
      }
//      printf_vector(unit_e);

      gsl_matrix_get_row(lo_up_row, lo_up, 0);
      lower = gsl_vector_max(lo_up_row);
      gsl_matrix_get_row(lo_up_row, lo_up, 1);
      upper = gsl_vector_min(lo_up_row);
      // cout<<lower<<"   "<<upper<<endl;
//printf_matrix(lo_up);//

      gsl_matrix_get_col(Gam_col, Gam, j);
      for (int ff = 0; ff < 10; ++ff)
      {
        r1[ff] = truncnorm.draw(0.0, sigma_r[j], lower, upper);
//        cout<<lower<<upper<<endl;
//    cout<<r1[ff]<<endl;
        gsl_vector_memcpy(tmp_unit_e, unit_e);
        gsl_vector_scale(tmp_unit_e, r1[ff]);
        // gsl_matrix_get_col(Gam_col, Gam, j);
        gsl_vector_add(tmp_unit_e, Gam_col);
        gsl_matrix_set_col(w1, ff, tmp_unit_e);
        double f1 = 0;
        double f2 = 0;
        double f3 = 0;
        double r12 = 0;
        double sigma_r2 = 0;
        for (int i = 0; i < n; ++i)
        {
          Uai0 = gsl_matrix_get(Ua, i, j);
          Uci0 = gsl_matrix_get(Uc, i, j);
          gsl_matrix_get_col(Bk__ij, Bk[j], i);
          gsl_matrix_get_col(bk__ij, bk[j], i);
          gsl_vector_mul(Bk__ij, tmp_unit_e);
          gsl_vector_mul(bk__ij, tmp_unit_e);
//          printf_vector(Bk__ij);
          f2 = VecSum(Bk__ij);
          f2 -= (Uai0 + Uci0);
          f2 *= f2;
          f2 /= (2 * Sigma_E);
          f3 = log(VecSum(bk__ij));
          f1 = f1 - f2 + f3;

        }       //Debug de 到这里
        r12 = r1[ff] * r1[ff];
        sigma_r2 = sigma_r[j] * sigma_r[j];
        f3 = f1 - txAx(tmp_unit_e, M) / (2 * tao[j]) - r12 / (2 * sigma_r2);
        gsl_matrix_set(f, 0, ff, f3);
// //   cout<<f3<<" "<<" "<<f1<<endl;
      }
// //    cout<<"fuck"<<endl;

//   printf_matrix(f);  //will r1 be too small?

//
      gsl_matrix_get_row(f1_row, f, 0);
      double f1_m = gsl_vector_max(f1_row);

      for (int ff = 0; ff < 10; ++ff)
      {
        double f_1_ff = exp(gsl_matrix_get(f, 0, ff) - f1_m);
        gsl_matrix_set(f, 0, ff, f_1_ff);
      }

      gsl_matrix_get_row(f1_row, f, 0);
      gsl_vector_scale(f1_row, 1 / VecSum(f1_row));

      int m = sp_int(rr, f1_row->size, f1_row->data);
//      cout << m << endl;
      gsl_matrix_get_col(wstar, w1, m);

      for (int i = 0; i < kk; ++i)
      {
        gsl_matrix_set(lo_up, 0, i, -infinity);
        gsl_matrix_set(lo_up, 1, i, infinity);
      }

      for (int ii = 0; ii < (kk - 1); ++ii)
      {
        Gam_0 = gsl_vector_get(wstar, ii);
        Gam_1 = gsl_vector_get(wstar, ii + 1);
        unit_e_0 = gsl_vector_get(unit_e, ii);
        unit_e_1 = gsl_vector_get(unit_e, ii + 1);
        lo_up_i = (Gam_0 - Gam_1) / (unit_e_1 - unit_e_0);
        if (unit_e_1 > unit_e_0)
          gsl_matrix_set(lo_up, 0, ii, lo_up_i);
        else
          gsl_matrix_set(lo_up, 1, ii, lo_up_i);
      }

      gsl_matrix_get_row(lo_up_row, lo_up, 0);
      lower = gsl_vector_max(lo_up_row);
      gsl_matrix_get_row(lo_up_row, lo_up, 1);
      upper = gsl_vector_min(lo_up_row);
      // cout<<lower<<" "<<upper<<endl;
//      printf_matrix(lo_up);

      for (int ff = 0; ff < 10; ++ff)
      {
        if (ff < 9)
        {
          r2[ff] =  truncnorm.draw(0, sigma_r[j], lower, upper);
          gsl_vector_memcpy(tmp_unit_e, unit_e);
          gsl_vector_scale(tmp_unit_e, r2[ff]);
          gsl_vector_add(tmp_unit_e, wstar);
        }
        else
        {
          r2[ff] = r1[m];
          gsl_matrix_get_col(tmp_unit_e, Gam, j);
        }
        gsl_matrix_set_col(w2, ff, tmp_unit_e);
        double f2 = 0;
        double f3 = 0;
        double f4 = 0;

        for (int i = 0; i < n; ++i)
        {
          Uai0 = gsl_matrix_get(Ua, i, j);
          Uci0 = gsl_matrix_get(Uc, i, j);
          gsl_matrix_get_col(Bk__ij, Bk[j], i);
          gsl_matrix_get_col(bk__ij, bk[j], i);
          gsl_vector_mul(Bk__ij, tmp_unit_e);
          gsl_vector_mul(bk__ij, tmp_unit_e);
          f3 = VecSum(Bk__ij) - Uai0 - Uci0;
          f3 = f3 * f3 / (2 * Sigma_E);
          f4 = log(VecSum(bk__ij));
          f2 = f2 - f3 + f4;
        }
        f3 = exp(f2 - txAx(tmp_unit_e, M) / (2 * tao[j]) - (r2[ff] * r2[ff]) / (2 * sigma_r[j] * sigma_r[j]) - f1_m);
        gsl_matrix_set(f, 1, ff, f3);
      }
//        printf_matrix(f);

      gsl_matrix_get_row(f1_row, f, 0);
      gsl_matrix_get_row(f2_row, f, 1);

      Accept = VecSum(f1_row) / VecSum(f2_row);
//    cout<<Accept<<endl;
      double u = gsl_rng_uniform(rr);

      if (u <= Accept)
      {
        gsl_matrix_set_col(Gam, j, wstar);
        fout<<lower<<" "<<upper<<endl;
        for(int i = 0; i < wstar->size; ++i)
          fout<<gsl_vector_get(wstar, i)<<" ";
        fout<<"\n"<<endl;
        gsl_vector_set(p_accept, j, gsl_vector_get(p_accept, j) + 1);
      }

//     printf_vector(wstar);

      gsl_matrix_get_col(Gam_col, Gam, j);
      tao[j] = 1 / gsl_ran_gamma(rr, 1 + (kk - 2) / 2, 1 / (0.005 + txAx(Gam_col, M) / 2));
//     cout<<tao[j]<<endl;
//      cout<<txAx(Gam_col, M)<<" ";

      if (GIB > GNUM)
      {
        gsl_blas_dgemv(CblasTrans, 1, BK[j], Gam_col, 0, tmpBK);
        gsl_matrix_set_col(GY[j], GIB - GNUM, tmpBK);
      }
    }
    //update gy

    for (int i = 0; i < 2; ++i)
    {
      gsl_matrix_get_col(Gam_col, Gam, i);
//      printf_vector(Gam_col);
      gsl_blas_dgemv(CblasTrans, 1, Bk[i], Gam_col, 0, Ycol);
      gsl_matrix_set_col(gy, i, Ycol);
    }

    // for (int i = 0; i < 2; ++i)
    // {
    //   for(int k = 0; k < Gam->size1; ++k)
    //     fout<<gsl_matrix_get(Gam, k, i)<<" ";
    //   fout<<"\n";
    // }
    // fout<<endl;


    // double bigger = 0;
    // for (int q = 0; q < 2; ++q)
    //   for (int i = 0; i < Bk[q]->size1; ++i)
    //     for (int k = 0; k < Bk[q]->size2; ++k)
    //     {if (gsl_matrix_get(Bk[q], i, k) < 0) bigger += 1;}
    // cout << bigger << endl;
    // for (int i = 0; i < gy->size1; ++i)
    //   for (int k = 0; k < gy->size2; ++k)
    //     {if (gsl_matrix_get(gy, i, k) < 0) bigger += 1;}
    /*for(int i =0; i < Gam_col->size; ++i)
      if(gsl_vector_get(Gam_col, i)<0) bigger += 1;
    cout << bigger <<endl;*/

    // printf_vector(Gam_col);
    // cout << "gy" << endl;
    // printf_matrix(gy);

    //compare
    if (GIB > (MCAX - 10))
    {
      for (int i = 0; i < n; ++i)
      {
        RowToMatrix(Xi, X, i, 2);
        gsl_blas_dgemv(CblasNoTrans, 1, Xi, Beta, 0, tmpx3);
        gsl_matrix_set_row(y_beta, i, tmpx3);
      }
      for (int i = 0; i < 2; ++i)
      {
        gsl_matrix_get_col(Ycol, gy, i);
        gsl_matrix_set_col(com_gy[GIB - (MCAX - 10)], i, Ycol);
        gsl_matrix_get_col(Ycol_tmp, y_beta, i);
        gsl_vector_sub(Ycol, Ycol_tmp);
        gsl_matrix_set_col(com_gy_beta[GIB - (MCAX - 10)], i, Ycol);
      }
    }
    //update Ua
    double Sigma_A2 = 1 / (2 / Sigma_E + 1 / Sigma_A);
    for (int i = 0; i < n1; ++i)
    {
      RowToMatrix(Xi, X, i, 2);
      gsl_matrix_get_row(Uci, Uc, i);
      gsl_matrix_get_row(Yi, gy, i);
      gsl_blas_dgemv(CblasNoTrans, 1, Xi, Beta, 1, Uci);
      gsl_vector_sub(Yi, Uci);
      double Mu_A = VecSum(Yi) * Sigma_A2 / Sigma_E;
      Uai0 = gsl_ran_gaussian(rr, sqrt(Sigma_A2));
      Uai0 += Mu_A;
      gsl_matrix_set(Ua, i, 0, Uai0);
      gsl_matrix_set(Ua, i, 1, Uai0);
    }

    gsl_matrix_set_identity(Sigma_Emat);
    gsl_matrix_scale(Sigma_Emat, 1 / Sigma_E);
    gsl_matrix_memcpy(InvADSigma_A, invA);
    gsl_matrix_scale(InvADSigma_A, 1 / Sigma_A);
    gsl_matrix_add(Sigma_Emat, InvADSigma_A);
    GInverse(Sigma_Emat);

    for (int i = n1; i < n; ++i)
    {
      gsl_matrix_get_row(Yi, gy, i);
      RowToMatrix(Xi, X, i, 2);
      gsl_matrix_get_row(Uci, Uc, i);
      gsl_blas_dgemv(CblasNoTrans, 1.0, Xi, Beta, 1.0, Uci);
      gsl_vector_sub(Yi, Uci);
      gsl_blas_dgemv(CblasNoTrans, 1 / Sigma_E, Sigma_Emat, Yi, 0.0, tmpx4);
      mvrnorm_cpp(rr, tmpx4, Sigma_Emat, tmpx3);
      gsl_matrix_set_row(Ua, i, tmpx3);
    }

    double Sigma_C2 = 1 / (2 / Sigma_E + 1 / Sigma_C);

    for (int i = 0; i < n; ++i)
    {
      gsl_matrix_get_row(Yi, Y, i);
      RowToMatrix(Xi, X, i, 2);
      gsl_matrix_get_row(Uai, Ua, i);
      gsl_blas_dgemv(CblasNoTrans, 1.0, Xi, Beta, 1.0, Uai);
      gsl_vector_sub(Yi, Uai);
      double Mu_C = VecSum(Yi) * Sigma_C2 / Sigma_E;
      Uci0 = gsl_ran_gaussian(rr, sqrt(Sigma_C2));
      Uci0 += Mu_C;
      gsl_matrix_set(Uc, i, 0, Uci0);
      gsl_matrix_set(Uc, i, 1, Uci0);
    }

    double shape = (n1 + 2 * (n - n1) + 18) / 2, rate = 4;
    for (int i = 0; i < n1; ++i)
    {
      double tmpUa = gsl_matrix_get(Ua, i, 0);
      tmpUa *= tmpUa;
      rate += tmpUa;
    }

    for (int i = n1; i < n; ++i)
    {
      gsl_matrix_get_row(Uai, Ua, i);
      rate += txAx(Uai, invA);
    }

    rate = rate / 2;
    Sigma_A = 1 / gsl_ran_gamma(rr, shape, 1 / rate);
    // cout<<"b1:"<<rate<<endl;
    gsl_matrix_set(MCMC, GIB, 0, Sigma_A);

    shape = (n + 18) / 2;
    rate = 0;
    for (int i = 0; i < n; ++i)
    {
      double tmpUc = gsl_matrix_get(Uc, i, 0);
      tmpUc *= tmpUc;
      rate += tmpUc;
    }

    rate = rate / 2 + 4;
    Sigma_C = 1 / gsl_ran_gamma(rr, shape, 1 / rate);
    // cout<<"b2:"<<rate<<endl;
    gsl_matrix_set(MCMC, GIB, 1, Sigma_C);

//    cout<<a<<" "<<b<<endl;
//
    shape = (2 * n + 18) / 2; rate = 0;

    for (int i = 0; i < n; ++i)
    {
      gsl_matrix_get_row(Yi, gy, i);
      RowToMatrix(Xi, X, i, 2);
      gsl_matrix_get_row(Uai, Ua, i);
      gsl_matrix_get_row(Uci, Uc, i);
      gsl_vector_add(Uai, Uci);
      gsl_blas_dgemv(CblasNoTrans, 1.0, Xi, Beta, 1.0, Uai);
      gsl_vector_sub(Yi, Uai);
      gsl_vector_mul(Yi, Yi);
      rate += VecSum(Yi);
    }

    rate = rate / 2 + 4;
    Sigma_E = 1 / gsl_ran_gamma(rr, shape, 1 / rate);
    // cout<<"b3:"<<rate<<endl;
    gsl_matrix_set(MCMC, GIB, 2, Sigma_E);
//
    gsl_matrix_set_zero(Sigma_beta);
    gsl_vector_set_zero(Mu_beta);

//    cout<<a<<" "<<b<<endl;

    for (int i = 0; i < n; ++i)
    {
      gsl_matrix_get_row(Yi, gy, i);
      RowToMatrix(Xi, X, i, 2);
      gsl_matrix_get_row(Uai, Ua, i);
      gsl_matrix_get_row(Uci, Uc, i);
      gsl_vector_add(Uai, Uci);
      gsl_vector_sub(Yi, Uai);

      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Xi, Xi, 0.0, tmpXi);
      gsl_matrix_add(Sigma_beta, tmpXi);

      gsl_blas_dgemv(CblasTrans, 1.0, Xi, Yi, 0.0, tmpx3);
      gsl_vector_add(Mu_beta, tmpx3);
    }


    GInverse(Sigma_beta);
    gsl_matrix_scale(Sigma_beta, Sigma_E);
    gsl_blas_dgemv(CblasNoTrans, 1 / Sigma_E, Sigma_beta, Mu_beta, 0.0, tmpx3);
    mvrnorm_cpp(rr, tmpx3, Sigma_beta, Beta);

    gsl_matrix_set(MCMC, GIB, 3, gsl_vector_get(Beta, 0));
    gsl_matrix_set(MCMC, GIB, 4, gsl_vector_get(Beta, 1));
  }

//  printf_matrix(Ua);
  colMeanSd(MCMC, Est, SEst);
//  cout<<"MCMC:"<<endl;
//  printf_matrix(MCMC);
  cout << "Est:" << endl;
  printf_vector(Est);
  cout << "SEst:" << endl;
  printf_vector(SEst);
//  cout<<"tao"<<endl;
//  cout<<tao[0]<<tao[1]<<endl;
//  cout<<"M:"<<endl;
//  printf_matrix(M);
  // cout << "gy" << endl;
  // printf_matrix(gy);
  // cout << "Bk_0" << endl;
  // printf_matrix(Bk[1]);
  cout << "OK" << endl;

  return 0;
}
