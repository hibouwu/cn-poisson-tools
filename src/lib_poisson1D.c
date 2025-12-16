/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int n = *la;
  int ldab = *lab;
  int diag = *kv + 1;
  int super = diag - 1;
  int sub = diag + 1;

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < ldab; i++) {
      AB[j * ldab + i] = 0.0; 
    }
    if (j > 0) {
      AB[j * ldab + super] = -1.0;
    }
    AB[j * ldab + diag] = 2.0;
    if (j < n - 1) {
      AB[j * ldab + sub] = -1.0;
    }
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int n = *la;
  int ldab = *lab;
  int diag = *kv + 1;

  int total = n * ldab;
  for (int idx = 0; idx < total; idx++) {
    AB[idx] = 0.0;
  }

  for (int j = 0; j < n; j++) {
    AB[j * ldab + diag] = 1.0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int n = *la;
  if (n <= 0) {
    return;
  }

  memset(RHS, 0, (size_t) n * sizeof(double));
  RHS[0] += (*BC0);
  RHS[n - 1] += (*BC1);
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int n = *la;
  double t0 = *BC0;
  double dt = *BC1 - *BC0;

  for (int i = 0; i < n; i++) {
    EX_SOL[i] = t0 + X[i] * dt;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int n = *la;
  double h = 1.0 / (double) (n + 1);

  for (int i = 0; i < n; i++) {
    x[i] = (i + 1) * h;
  }
}

double relative_forward_error(double* x, double* y, int* la){
  int n = *la;
  double *work = (double *) malloc((size_t) n * sizeof(double));
  if (work == NULL) {
    return DBL_MAX;
  }

  cblas_dcopy(n, y, 1, work, 1);
  cblas_daxpy(n, -1.0, x, 1, work, 1);

  double num = cblas_dnrm2(n, work, 1);
  double den = cblas_dnrm2(n, y, 1);
  free(work);

  if (den == 0.0) {
    return (num == 0.0) ? 0.0 : DBL_MAX;
  }
  return num / den;
}

int indexABCol(int i, int j, int *lab){
  return j * (*lab) + i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  int ncols = *n;
  int ldab = *lab;
  int diag = *kl + *ku;
  int sub  = diag + *kl;
  int super = diag - *ku;

  *info = 0;
  if (ncols <= 0) {
    return *info;
  }

  for (int i = 0; i < ncols; i++) {
    ipiv[i] = i + 1;
  }

  for (int j = 0; j < ncols - 1; j++) {
    double pivot = AB[j * ldab + diag];
    if (pivot == 0.0) {
      *info = j + 1;
      return *info;
    }
    double factor = AB[j * ldab + sub] / pivot;
    AB[j * ldab + sub] = factor;

    AB[(j + 1) * ldab + diag] -= factor * AB[(j + 1) * ldab + super];
  }

  if (AB[(ncols - 1) * ldab + diag] == 0.0) {
    *info = ncols;
  }
  return *info;
}
