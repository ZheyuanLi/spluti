#include <stdlib.h>
#include <Rinternals.h>
void EvalGrid (double *xk, int k, int nx, double *x) {
  int np = k - 1;
  double step0 = 1.0 / (nx - 1), step, xtmp;
  double *xend = x + np * nx, *pend = x + nx, *ptrx, *ptrxk;
  for (ptrx = x, ptrxk = xk; ptrx < xend; ptrxk++, pend += nx) {
    step = (ptrxk[1] - ptrxk[0]) * step0;
    for (xtmp = *ptrxk; ptrx < pend; xtmp += step, ptrx++) *ptrx = xtmp;
    }
  }
SEXP C_EvalGrid (SEXP xk, SEXP nx) {
  int k = length(xk), c_nx = asInteger(nx);
  SEXP x = PROTECT(allocVector(REALSXP, (k - 1) * c_nx));
  EvalGrid(REAL(xk), k, c_nx, REAL(x));
  UNPROTECT(1);
  return x;
  }
