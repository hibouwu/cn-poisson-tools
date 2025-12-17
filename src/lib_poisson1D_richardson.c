/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>

void eig_poisson1D(double* eigval, int *la){
  for (int k = 0; k < *la; k++) {
    eigval[k] = 2.0 - 2.0 * cos((k + 1) * M_PI / ((*la + 1)));
  }
}

double eigmax_poisson1D(int *la){
  return 2.0 - 2.0 * cos((*la) * M_PI / ((*la + 1)));
}

double eigmin_poisson1D(int *la){
  return 2.0 - 2.0 * cos(M_PI / ((*la + 1)));
}

double richardson_alpha_opt(int *la){
  double min_eig = eigmin_poisson1D(la);
  double max_eig = eigmax_poisson1D(la);
  return 2.0 / (min_eig + max_eig);
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  double *r = (double *) malloc((size_t)(*la) * sizeof(double));
  double norm_b = cblas_dnrm2(*la, RHS, 1);
  if (norm_b == 0.0) {norm_b = 1.0;}
  for (*nbite = 0; *nbite < *maxit; (*nbite)++) {
    // r = b
    cblas_dcopy(*la, RHS, 1, r, 1);
    // r = b - A * x
    // y = alpha * A * x + beta * y
    // We want r = r - A * x => alpha = -1, beta = 1
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, r, 1);
    double norm_r = cblas_dnrm2(*la, r, 1);
    resvec[*nbite] = norm_r / norm_b;
    if (resvec[*nbite] < *tol) break;
    // x = x + alpha * r
    cblas_daxpy(*la, *alpha_rich, r, 1, X, 1);
  }
  free(r);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // Initialize MB to 0 and copy the diagonal from AB.
  memset(MB, 0, (size_t)(*la) * (*lab) * sizeof(double));
  for (int j = 0; j < *la; j++) {
    // AB has Diag at row 'ku' (standard GB format)
    // MB has Diag at row 'kv + 1' (as defined by setup)
    MB[j * (*lab) + (*kv + 1)] = AB[j * (*lab) + (*ku)]; 
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  memset(MB, 0, (size_t)(*la) * (*lab) * sizeof(double));
  for (int j = 0; j < *la; j++) {
      // Diagonal
      // AB Diag at 'ku', MB Diag at 'kv + 1'
      MB[j * (*lab) + (*kv + 1)] = AB[j * (*lab) + (*ku)];
      
      // Sub-diagonal (Lower part E)
      // AB Sub at 'ku + 1'
      if (j < *la - 1) { 
         MB[j * (*lab) + (*kv + 2)] = AB[j * (*lab) + (*ku + 1)];
      }
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  double *r = (double *) malloc((size_t)(*la) * sizeof(double));
  double *z = (double *) malloc((size_t)(*la) * sizeof(double)); // Update vector M^{-1} r
  
  double norm_b = cblas_dnrm2(*la, RHS, 1);
  int lab_ab = *kl + *ku + 1; // AB stride (packed)
  if (norm_b == 0.0) {norm_b = 1.0;}
  
  for (*nbite = 0; *nbite < *maxit; (*nbite)++) {
    // r = b
    cblas_dcopy(*la, RHS, 1, r, 1);
    
    // r = b - A * x
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, lab_ab, X, 1, 1.0, r, 1);
    
    double norm_r = cblas_dnrm2(*la, r, 1);
    resvec[*nbite] = norm_r / norm_b;
    
    if (resvec[*nbite] < *tol) break;
    
    // Solve M z = r
    cblas_dcopy(*la, r, 1, z, 1);
    
    // Forward substitution
    for (int i = 0; i < *la; i++) {
        double val = r[i];
        if (i > 0) {
            val -= MB[(i-1) * (*lab) + (*ku + 1)] * z[i-1]; // M_{i, i-1} * z_{i-1}
            // Wait. Column major.
            // M_{i, i-1} is Sub-diagonal at column i-1.
            // In MB array at column i-1: MB[(i-1) * (*lab) + (*ku + 1)].
            // YES. 
        }
        z[i] = val / MB[i * (*lab) + (*ku)];
    }
    
    // x = x + z
    cblas_daxpy(*la, 1.0, z, 1, X, 1);
  }
  
  free(r);
  free(z);
}
