#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
extern SEXP C_EvalGrid(SEXP, SEXP);
extern SEXP C_multiplicity(SEXP, SEXP, SEXP);
extern SEXP C_polydiv(SEXP, SEXP);
static const R_CallMethodDef CallEntries[] = {
    {"C_EvalGrid", (DL_FUNC) &C_EvalGrid, 2},
    {"C_multiplicity", (DL_FUNC) &C_multiplicity, 3},
    {"C_polydiv", (DL_FUNC) &C_polydiv, 2},
    {NULL, NULL, 0}
};
void R_init_spluti(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
