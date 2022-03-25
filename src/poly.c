#include <Rinternals.h>
#include <stdlib.h>
static inline void polydiv_inplace (double *F, size_t p, double *G, size_t q) {
  double *Gq, *Q, *ptr_G, *ptr_F, b;
  Gq = G + q;
  Q = F + q;
  F += p;
  while (F >= Q) {
    b = *F / (*Gq); *F = b;
    ptr_G = Gq; ptr_F = F;
    while (ptr_G > G) {
      ptr_F--; ptr_G--;
      *ptr_F -= b * (*ptr_G);
    }
    F--;
  }
}
void polydiv (double *F, size_t p, double *G, size_t q, double *out) {
  size_t i;
  for (i = 0; i <= p; i++) out[i] = F[i];
  polydiv_inplace(out, p, G, q);
}
SEXP C_polydiv (SEXP F, SEXP G) {
  size_t p = length(F) - 1;
  size_t q = length(G) - 1;
  SEXP out = PROTECT(allocVector(REALSXP, p + 1));
  polydiv(REAL(F), p, REAL(G), q, REAL(out));
  UNPROTECT(1);
  return out;
}
size_t mfactor (double *F, size_t p, double *G, size_t q, double tol, double *work) {
  size_t i;
  double r, tmp, *ptr;
  double *Q = work + q;
  for (i = 0; i <= p; i++) work[i] = F[i];
  i = 0;
  while (1) {
    if (p < q) break;
    polydiv_inplace(work, p, G, q);
    r = fabs(*work);
    for (ptr = work + 1; ptr < Q; ptr++) {
      tmp = fabs(*ptr);
      if (tmp > r) r = tmp;
    }
    if (r > tol) break;
    p -= q; work += q; i++;
  }
  return i;
}
void multiplicity (double *F, size_t p, double *G, size_t q, double tol, double *work,
                   size_t *m, size_t M) {
  for (size_t j = 0; j < M; j++) m[j] = mfactor(F, p, G + j * (q + 1), q, tol, work);
}
SEXP C_multiplicity (SEXP F, SEXP G, SEXP tol) {
  size_t p = length(F) - 1;
  size_t q = nrows(G) - 1;
  size_t M = ncols(G);
  void *heap = malloc((p + 1 + M) * sizeof(double));
  double *work = (double *)heap;
  size_t *m = (size_t *)(work + p + 1);
  multiplicity(REAL(F), p, REAL(G), q, asReal(tol), work, m, M);
  SEXP SEXP_m = PROTECT(allocVector(INTSXP, M));
  int *int_m = INTEGER(SEXP_m);
  for (p = 0; p < M; p++) int_m[p] = (int)m[p];
  free(heap);
  UNPROTECT(1);
  return SEXP_m;
}
