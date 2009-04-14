#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP
g_ci_test(SEXP S, SEXP NR, SEXP iR, SEXP jR, SEXP C);

SEXP
is_linked(SEXP amat, SEXP iR, SEXP jR);

SEXP
cond_ind_test(SEXP tb, SEXP nvar, SEXP i0, SEXP i1);

SEXP
compress_freq_table(SEXP tb);





static R_CallMethodDef
callEntries[] = {
  {"g_ci_test", (DL_FUNC) &g_ci_test,5},
  {"is_linked", (DL_FUNC) &is_linked,3},
  {"cond_ind_test", (DL_FUNC) &cond_ind_test, 4},
  {"compress_freq_table", (DL_FUNC) &compress_freq_table, 1},
  {NULL}
};

void
R_init_lcd(DllInfo *info)
{
  /* Register routines, allocate resources. */
  R_registerRoutines(info, NULL, callEntries /*CallEntries*/,
		     NULL, NULL /*ExternEntries*/);
}
          
void
R_unload_lcd(DllInfo *info)
{
  /* Release resources. */
}
