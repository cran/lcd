#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>


/*
 * Matrix computation functions
 */

/*
  FUNCTION: matprod
  PURPOSE: multiply two matrices by using the LaPACK library that
           comes along with the R distribution, this code is taken from

           R-2.2.0/src/main/array.c
  RETURNS: none
*/

void
matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z) {
    char *transa = "N", *transb = "N";
    int i,  j, k;
    double one = 1.0, zero = 0.0, sum;
    Rboolean have_na = FALSE;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        /* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
         * The test is only O(n) here
         */
        for (i = 0; i < nrx*ncx; i++)
            if (ISNAN(x[i])) {have_na = TRUE; break;}
        if (!have_na)
            for (i = 0; i < nry*ncy; i++)
                if (ISNAN(y[i])) {have_na = TRUE; break;}
        if (have_na) {
            for (i = 0; i < nrx; i++)
                for (k = 0; k < ncy; k++) {
                    sum = 0.0;
                    for (j = 0; j < ncx; j++)
                        sum += x[i + j * nrx] * y[j + k * nry];
                    z[i + k * nrx] = sum;
                }
        } else
            F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
                            x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
        for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}



/*
  FUNCTION: matinv
  PURPOSE: calculates de inverse of a matrix by using the LaPACK library
           that comes along with the R distribution, this code is taken from
           the function modLa_dgesv in file

           R-2.2.0/src/modules/lapack/Lapack.c
  RETURNS: none
*/

void
matinv(double* inv, double* M, int n) {
  int     i,j;
  int     info;
  int*    ipiv;
  double* avals;
  double* work;
  double  anorm;
  double  rcond;
  double  tol = DBL_MIN;

  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      inv[i+j*n] = i == j ? 1.0 : 0.0; 

  ipiv = (int *) Calloc(n,double);
  avals = (double *) Calloc(n*n,double);
  Memcpy(avals,M,(size_t) (n*n));

  F77_CALL(dgesv)(&n,&n,avals,&n,ipiv,inv,&n,&info);
  if (info < 0)
    error("argument %d of Lapack routine %s had invalid value",-info, "dgesv");
  if (info > 0)
    error("Lapack routine dgesv: system is exactly singular");

  anorm = F77_CALL(dlange)("1", &n, &n, M, &n, (double*) NULL);

  work = (double *) Calloc(4*n,double);

  F77_CALL(dgecon)("1", &n, avals, &n, &anorm, &rcond, work, ipiv, &info);
  if (rcond < tol)
    error("system is computationally singular: reciprocal condition number = %g",rcond);

  Free(ipiv);
  Free(avals);
  Free(work);
}


/*
 * C routine that performs the conditional independence test
 * ---------------------------------------------------------
 * Note: code taken from qp_0.2-1/src/qp.c
 */

double
g_ci_test_c(double* S, int n_var, int N, int i, int j, int* Q, int q) {
  int*    subvars;
  int     subn = q + 2;
  int     k,l;
  double* Mmar;
  double  S11;
  double* S12;
  double* S21;
  double* S22;
  double* S22inv;
  double* S22inv1col;
  double* tmpmat;
  double  tmpval;
  double  betahat;
  double  sigma;
  double  se;
  double  t_value;

  subvars     = Calloc(subn,int);
  Mmar        = Calloc(subn*subn,double);
  S12         = Calloc(subn,double);
  S21         = Calloc(subn,double);
  S22         = Calloc((subn-1)*(subn-1),double);
  S22inv      = Calloc((subn-1)*(subn-1),double);
  S22inv1col  = Calloc(subn-1,double);

  subvars[0] = i; /* order here is important, first variable i */
  subvars[1] = j; /* then variable j then the conditioning set */
  for (k=2;k<subn;k++)
    subvars[k] = Q[k-2];

  for (k=0;k<subn;k++)
    for (l=0;l<subn;l++) {
      Mmar[k+l*subn] = S[subvars[k]+subvars[l]*n_var];
      if (k == 0 && l > 0)
        S12[l-1] = Mmar[k+l*subn];
      if (k > 0 && l == 0)
        S21[k-1] = Mmar[k+l*subn];
      if (k > 0 && l > 0)
        S22[k-1+(l-1)*(subn-1)] = Mmar[k+l*subn];
    }
  S11 = Mmar[0];

  matinv(S22inv,S22,subn-1);

  Memcpy(S22inv1col,S22inv,(size_t) (subn-1));
  matprod(S12,1,subn-1,S22inv1col,subn-1,1,&betahat);

  tmpmat = Calloc(subn-1,double);
  matprod(S22inv,subn-1,subn-1,S21,subn-1,1,tmpmat);
  matprod(S12,1,subn-1,tmpmat,subn-1,1,&tmpval);
  Free(tmpmat);
  sigma = sqrt( (S11 - tmpval) * (N - 1) / (N - subn) );
  se = sigma * sqrt(S22inv[0] / (N - 1));
  t_value = betahat / se;

  Free(S22inv1col);
  Free(S22inv);
  Free(S22);
  Free(S21);
  Free(S12);
  Free(Mmar);
  Free(subvars);

  return t_value;
}


/*
 * Function: dev_ci_test_c
 * ---------------------
 * Purpose: perform a conditional independence test for normally 
 * distributed data using the deviance proposed in Whittaker (1989)
 * pp180-181
 * ---------------------
 * Code written by Zongming Ma <zongming.ma@gmail.com>
 */

double
dev_ci_test_c(double * S, int n_var, int N, int i, int j, int *Q, int q) 
{
    int * subvars;
    int subn = q + 2;
    int k,l,m;
    double* Mmar;
    double* S11;
    double* S12;
    double* S21;
    double* S22;
    double* S22inv;
    double* S11_2;
    double chi_value;
    double* tmpmat1;
    double* tmpmat;
    
    subvars = Calloc(subn, int);
    Mmar    = Calloc(subn*subn, double);
    
    S11     = Calloc(2*2, double);
    S11_2   = Calloc(2*2, double);
    S12     = Calloc(2*(subn-2), double);
    S21     = Calloc((subn-2)*2, double);
    S22     = Calloc((subn-2)*(subn-2), double);
    S22inv  = Calloc((subn-2)*(subn-2), double);
    tmpmat1 = Calloc(2*(subn-2), double);
    tmpmat  = Calloc(2*2, double);
    subvars[0] = i; /* order here is important, first variable i */
    subvars[1] = j; /* then variable j then the conditioning set */
    for (k=2;k<subn;k++)
	subvars[k] = Q[k-2];
    
    for (k = 0; k < subn; k++)
	for (l = 0; l < subn; l++) {
	    Mmar[k + l*subn] = S[subvars[k] + subvars[l]*n_var];
	    if (k <= 1 && l > 1)
		S12[k + (l-2)*2] = Mmar[k + l*subn];
	    if (k > 1 && l <= 1)
		S21[k-2 + l*(subn-2)] = Mmar[k + l*subn];
	    if (k > 1 && l > 1)
		S22[k-2 + (l-2)*(subn-2)] = Mmar[k + l*subn];
	    if (k <= 1 && l <= 1)
		S11[k + l*2] = Mmar[k + l*subn];
	}
    
    if (subn - 2 > 0) {
	matinv(S22inv, S22, subn - 2); /* calculate S22inv */
	matprod(S12, 2, subn - 2, S22inv, subn - 2, subn - 2, tmpmat1);
	matprod(tmpmat1, 2, subn - 2, S21, subn - 2, 2, tmpmat);
	
	for (m = 0; m < 4; m++)
	    S11_2[m] = S11[m] - tmpmat[m];
	
	chi_value = N * log(1. - S11_2[1]*S11_2[2]/(S11_2[0]*S11_2[3]));
    } else {
	chi_value = N * log(1. - Mmar[1]*Mmar[2]/(Mmar[0]*Mmar[3]));
    }
    
    Free(tmpmat);
    Free(tmpmat1);
    Free(S22inv);
    Free(S22);
    Free(S21);
    Free(S12);
    Free(S11_2);
    Free(S11);
    Free(Mmar);
    Free(subvars);
    
    return 0. - chi_value;
}

    




/*
 * Function: g_ci_test
 * -------------------
 * Purpose: perform a conditional independence test for normally
 * distributed data
 * ---------------------------------------
 * Note: code taken from qp_0.2-1/src/qp.c; this is an R-dependent wrapper
 */
SEXP
g_ci_test(SEXP S, SEXP NR, SEXP iR, SEXP jR, SEXP C){
  int    N = INTEGER(NR)[0];
  int    n_var = INTEGER(getAttrib(S,R_DimSymbol))[0];
  int    q;
  int*   cond;
  int    i,j,k;
  double p_value;
  double chi_value;
  /*   double t_value; */
  SEXP   result;
/*   SEXP   result_names; */
/*   SEXP   result_t_val; */
/*   SEXP   result_p_val; */

  PROTECT_INDEX Spi,Cpi;

  PROTECT_WITH_INDEX(S,&Spi);
  PROTECT_WITH_INDEX(C,&Cpi);

  REPROTECT(S = coerceVector(S,REALSXP),Spi);
  REPROTECT(C = coerceVector(C,INTSXP),Cpi);

  i = INTEGER(iR)[0] - 1;
  j = INTEGER(jR)[0] - 1;
  q = length(C);

  cond = Calloc(q, int);
  for (k=0;k<q;k++)
    cond[k] = INTEGER(C)[k]-1;

/*   t_value = g_ci_test_c(REAL(S),n_var,N,i,j,cond,q); */
  chi_value = dev_ci_test_c(REAL(S),n_var,N,i,j,cond,q);
/*   p_value = 2.0 * (1.0 - pt(fabs(t_value),N-q-2,1,0)); */
  p_value = 1. - pchisq(chi_value, 1., 1, 0);

  PROTECT(result = allocVector(REALSXP,3));
  /* SET_VECTOR_ELT(result,0,result_t_val = allocVector(REALSXP,1)); */
/*   SET_VECTOR_ELT(result,1,result_p_val = allocVector(REALSXP,1)); */
/*   PROTECT(result_names = allocVector(STRSXP,2)); */
/*   SET_STRING_ELT(result_names,0,mkChar("t.value")); */
/*   SET_STRING_ELT(result_names,1,mkChar("p.value")); */
/*   setAttrib(result,R_NamesSymbol,result_names); */
/*   REAL(VECTOR_ELT(result,0))[0] = t_value; */
/*   REAL(VECTOR_ELT(result,1))[0] = p_value; */
  REAL(result)[2] = p_value;
  REAL(result)[0] = chi_value;
  /* REAL(result)[1] = N-q-2; */
  REAL(result)[1] = 1;
  UNPROTECT(3); /* S C result result_names */

  Free(cond);

  return result;
}



