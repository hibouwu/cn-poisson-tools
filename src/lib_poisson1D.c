/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>
#include <stdlib.h>
#include <float.h>

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int n = *la;      // Number of unknowns
  int ldab = *lab;  // Leading dimension of AB
  int diag = *kv + 1;   // Row index of the diagonal
  int super = diag - 1; // Row index of the super-diagonal
  int sub = diag + 1;   // Row index of the sub-diagonal
  // Initialize the whole matrix storage to zero
  memset(AB, 0, (size_t)n * ldab * sizeof(double));
  // Set up the tridiagonal matrix for 1D Poisson: -1, 2, -1
  for (int j = 1; j < n; j++) {AB[j * ldab + super] = -1.0;}
  for (int j = 0; j < n; j++) {AB[j * ldab + diag] = 2.0;}
  for (int j = 0; j < n - 1; j++) {AB[j * ldab + sub] = -1.0;}
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int n = *la;
  int ldab = *lab;
  int diag = *kv + 1; // Row index of the diagonal
  // Initialize the whole matrix storage to zero
  int total = n * ldab;
  for (int idx = 0; idx < total; idx++) {AB[idx] = 0.0;}
  // Set diagonal elements to 1
  for (int j = 0; j < n; j++) {AB[j * ldab + diag] = 1.0;}
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int n = *la;
  if (n <= 0) {return;}
  // Initialize RHS to zero
  memset(RHS, 0, (size_t) n * sizeof(double));
  RHS[0] += (*BC0);      // T0 dans le premier point (boundary T0)
  RHS[n - 1] += (*BC1);  // T1 dans le dernier point (boundary T1)
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int n = *la;
  double t0 = *BC0;
  double dt = *BC1 - *BC0;
  // Linear solution between BC0 and BC1
  for (int i = 0; i < n; i++) {EX_SOL[i] = t0 + X[i] * dt;}
}

void set_grid_points_1D(double* x, int* la){
  int n = *la;
  double h = 1.0 / (double) (n + 1); // taille du pas
  // Set grid points excluding boundaries
  // x[i] = h, 2h, ..., nh, positions between 0 and 1
  for (int i = 0; i < n; i++) {x[i] = (i + 1) * h;}
}

double relative_forward_error(double* x, double* y, int* la){
  int n = *la;
  double *work = (double *) malloc((size_t) n * sizeof(double));
  if (work == NULL) {return DBL_MAX;}
  // Compute work = x - y
  cblas_dcopy(n, y, 1, work, 1);       // work = y
  cblas_daxpy(n, -1.0, x, 1, work, 1); // work = work - x = y - x
  // Compute norms
  double num = cblas_dnrm2(n, work, 1); // ||x - y||
  double den = cblas_dnrm2(n, y, 1);    // ||y|| (reference)
  free(work);
  if (den == 0.0) {return (num == 0.0) ? 0.0 : DBL_MAX;}
  return num / den;
}

int indexABCol(int i, int j, int *lab){return j * (*lab) + i;}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  int ncols = *n;
  int ldab = *lab;
  int diag = *kl + *ku;   // Row index of diagonal
  int sub  = diag + *kl;  // Row index of sub-diagonal
  int super = diag - *ku; // Row index of super-diagonal
  *info = 0;
  if (ncols <= 0) {
    return *info;
  }
  // Initialize pivot indices (identity permutation, no pivoting implemented)
  for (int i = 0; i < ncols; i++) {ipiv[i] = i + 1;}
  // Gaussian elimination for tridiagonal matrix
  for (int j = 0; j < ncols - 1; j++) {
    double pivot = AB[j * ldab + diag];
    if (pivot == 0.0) {
      *info = j + 1; // Singular matrix
      return *info;
    }
    // Calculate multiplier (factor) for the sub-diagonal element
    double factor = AB[j * ldab + sub] / pivot;
    AB[j * ldab + sub] = factor; // Store L part in place
    // Update the next diagonal element
    AB[(j + 1) * ldab + diag] -= factor * AB[(j + 1) * ldab + super];
  }
  // Check the last diagonal element
  if (AB[(ncols - 1) * ldab + diag] == 0.0) {
    *info = ncols;
  }
  return *info;
}
