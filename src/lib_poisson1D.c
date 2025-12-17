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
  // Initialize the whole matrix storage to zero
  memset(AB, 0, (size_t)(*la) * (*lab) * sizeof(double));
  // Set up the tridiagonal matrix for 1D Poisson: -1, 2, -1
  for (int j = 1; j < *la; j++) {AB[j * (*lab) + (*kv)] = -1.0;}
  for (int j = 0; j < *la; j++) {AB[j * (*lab) + (*kv + 1)] = 2.0;}
  for (int j = 0; j < *la - 1; j++) {AB[j * (*lab) + (*kv + 2)] = -1.0;}
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // Initialize the whole matrix storage to zero
  memset(AB, 0, (size_t)(*la) * (*lab) * sizeof(double));
  // Set diagonal elements to 1
  for (int j = 0; j < *la; j++) {AB[j * (*lab) + (*kv + 1)] = 1.0;}
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  if (*la <= 0) {return;}
  // Initialize RHS to zero
  memset(RHS, 0, (size_t)(*la) * sizeof(double));
  RHS[0] += (*BC0);      // T0 dans le premier point (boundary T0)
  RHS[*la - 1] += (*BC1);  // T1 dans le dernier point (boundary T1)
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // Linear solution between BC0 and BC1
  for (int i = 0; i < *la; i++) {EX_SOL[i] = (*BC0) + X[i] * ((*BC1) - (*BC0));}
}

void set_grid_points_1D(double* x, int* la){
  double h = 1.0 / (double) (*la + 1); // taille du pas
  // Set grid points excluding boundaries
  // x[i] = h, 2h, ..., nh, positions between 0 and 1
  for (int i = 0; i < *la; i++) {x[i] = (i + 1) * h;}
}

double relative_forward_error(double* x, double* y, int* la){
  double *work = (double *) malloc((size_t) (*la) * sizeof(double));
  if (work == NULL) {return DBL_MAX;}
  // Compute work = x - y
  cblas_dcopy(*la, y, 1, work, 1);       // work = y
  cblas_daxpy(*la, -1.0, x, 1, work, 1); // work = work - x = y - x
  // Compute norms
  double num = cblas_dnrm2(*la, work, 1); // ||x - y||
  double den = cblas_dnrm2(*la, y, 1);    // ||y|| (reference)
  free(work);
  if (den == 0.0) {return (num == 0.0) ? 0.0 : DBL_MAX;}
  return num / den;
}

int indexABCol(int i, int j, int *lab){return j * (*lab) + i;}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  *info = 0;
  if (*n <= 0) {return *info;}
  // Initialize pivot indices (identity permutation, no pivoting implemented)
  for (int i = 0; i < *n; i++) {ipiv[i] = i + 1;}
  // Gaussian elimination for tridiagonal matrix
  for (int j = 0; j < *n - 1; j++) {
    double pivot = AB[j * (*lab) + (*kl + *ku)];
    if (pivot == 0.0) {
      *info = j + 1; // Singular matrix
      return *info;
    }
    // Calculate multiplier (factor) for the sub-diagonal element
    double factor = AB[j * (*lab) + (*ku + 2 * (*kl))] / pivot;
    AB[j * (*lab) + (*ku + 2 * (*kl))] = factor; // Store L part in place
    // Update the next diagonal element
    AB[(j + 1) * (*lab) + (*kl + *ku)] -= factor * AB[(j + 1) * (*lab) + (*kl)];
  }
  // Check the last diagonal element
  if (AB[(*n - 1) * (*lab) + (*kl + *ku)] == 0.0) {
    *info = *n;
  }
  return *info;
}
