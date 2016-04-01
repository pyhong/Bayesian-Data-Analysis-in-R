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
#include <limits>

//rtruncnorm
#define TAIL_LIMIT  2.983851594898812
#define EXPTLOVTL  0.00390733345803262
#define LOG_TAIL_LIMIT  1.09321494749176
#define SQ_TAIL_LIMIT  4.45168517019009
#define MAX_PTAIL  0.00390733345803262
#define basep   0.0078125
#define infinity (std::numeric_limits<double>::infinity() )
#define min_double (std::numeric_limits<double>::min())
#define max_double (std::numeric_limits<double>::max())
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

class sample_int
{
private:
  int * all;
  double * prob;
  double * prob_0;
  gsl_rng * rr;
  int length;
  double sum;

public:
  sample_int(int *ALL, double *PROB, int n)
  {
    rr = gsl_rng_alloc( gsl_rng_rand );
    length = n;
    all = new int(length);
    memcpy(all, ALL, length * sizeof(int));
    prob = new double(length);
    memcpy(prob, PROB, length * sizeof(double));
    prob_0 = new double(length);
    memcpy(prob_0, PROB, length * sizeof(double));
  }

  void reset()
  {
    memcpy(prob_0, prob, length * sizeof(double));
  }

  ~sample_int()
  {
    delete prob;
    delete all;
    gsl_rng_free(rr);
  }

  int between(double x)
  {
    int i = 0;
    double cdf = prob_0[0];

    while (cdf / sum < x)
    {
      i = i + 1;
      cdf += prob_0[i];
    }

    return i;
  }

  int draw_noreplace()
  {
    double x = gsl_rng_uniform(rr);
    for (int j = 0; j < length; ++j)
      sum += prob_0[j];
    int sample_i = between(x);
    prob_0[sample_i] = 0;
    return all[sample_i];
  }

  int draw()
  {
    double x = gsl_rng_uniform(rr);
    int sample_i = between(x);
    return all[sample_i];
  }

};



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

class ACE
{
private:
  const gsl_rng_type *T;
  gsl_rng *rr;
  int n, n1, p, MCAX, GNUM, MCMC_times, kk;
  gsl_matrix *Y, *Ystar, *X[2], *A, *invA, *Ua, *Uc, *Est, *SEst, *Sigma_beta,
             *Xi, *tmpXi, *SigmaAA, *Sigma_Emat, *InvADSigma_A, *D0, *D1, *N, *M;
  gsl_vector *Beta, *Mu_beta, *Yi, *Ue, *unit_e, *tmpMu_e, *Bk[2], *bk[2], *ks, *Ycol,
             *Ycol_tmp;
  double tmpA[4] = {1.0, 0.5, 0.5, 1.0}, a[2] = { -1, 1}, *MCMC;
  double Sigma_A, Sigma_C, Sigma_E;

public:
  void init_array(gsl_matrix * a[], const int array_size, const int row, const int col)
  {
    for (int i = 0; i < array_size; ++i)
    {
      a[i] = gsl_matrix_alloc(row, col);
    }
  }
  ACE()
  {
    gsl_rng_env_setup();
    gsl_rng_default = 0;
    T = gsl_rng_default;
    rr = gsl_rng_alloc(T);
    n = 1000;
    n1 = 500;
    p = 2;
    MCAX = 2000;
    GNUM = 1800;
    kk = 25;
    Y = gsl_matrix_alloc(n, p);
    Ystar = gsl_matrix_alloc(n, p);
    init_array(X, 2, n, p);
    gsl_matrix_view tmpA_view = gsl_matrix_view_array(tmpA, 2, 2);
    A = gsl_matrix_alloc(2, 2);
    invA = gsl_matrix_alloc(2, 2);
    gsl_matrix_memcpy(A, &tmpA_view.matrix);
    Beta = gsl_vector_alloc(2);
    gsl_vector_set(Beta, 0, 0.5);
    gsl_vector_set(Beta, 1, 1.0);
    Sigma_A = 0.4;
    Sigma_C = 0.3;
    Sigma_E = 0.2;
    gsl_matrix_memcpy(invA, A);
    GInverse(invA);
    Ua = gsl_matrix_alloc(n, 2);
    Uc = gsl_matrix_alloc(n, 2);
    Est = gsl_vector_alloc(5);
    SEst = gsl_vector_alloc(5);
    Sigma_beta = gsl_matrix_alloc(p, p);
    Mu_beta = gsl_vector_alloc(p);
    tmpmu = gsl_vector_calloc(2);
    MCMC =  new double[MCAX][5];
    Xi = gsl_matrix_alloc(2, 2);
    tmpXi = gsl_matrix_alloc(2, 2);
    tmpx3 = gsl_vector_alloc(2);
    tmpx4 = gsl_vector_alloc(2);
    Ue = gsl_vector_alloc(2);
    SigmaAA - gsl_matrix_alloc(2, 2);
    Sigma_Emat = gsl_matrix_alloc(2, 2);
    InvADSigma_A = gsl_matrix(2, 2);
    D0 = gsl_matrix_calloc(kk - 1, kk);
    D1 = gsl_matrix_calloc(kk - 1 - 1, kk - 1);
    unit_e = gsl_vector_calloc(kk);
    tmpMu_e = gsl_vector_calloc(kk);
    init_array(Bk, 2, kk, n);
    init_array(bk, 2, kk, n);
    ks = gsl_matrix_calloc(2, kk + 4 - 2);




  }


};

int main()
{
  ACE ace;
  return 0;
}

