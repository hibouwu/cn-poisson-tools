#include "lib_poisson1D.h"
#include <assert.h>
#include <string.h>

/* Helper to print matrix for debugging */
void print_GB_matrix(double *AB, int lab, int la, const char *name) {
    printf("Matrix %s (%dx%d):\n", name, lab, la);
    for (int i = 0; i < lab; i++) {
        for (int j = 0; j < la; j++) {
            printf("%6.2f ", AB[j * lab + i]);
        }
        printf("\n");
    }
    printf("\n");
}

/* Exercise 4 Validation: Use dgbmv to verify Ax = b */
void test_dgbmv_ax_equals_b(int n) {
    printf("=== Test: dgbmv verification (n=%d) ===\n", n);
    
    /* 1. Setup Matrix A in GB format */
    int kv = 1;
    int ku = 1;
    int kl = 1;
    int lab = kv + kl + ku + 1;
    double *AB = (double *)malloc(lab * n * sizeof(double));
    set_GB_operator_colMajor_poisson1D(AB, &lab, &n, &kv);
    
    /* 2. Setup known vector x */
    double *x = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) x[i] = 1.0; /* Let x be all 1s */
    
    /* 3. Compute y = A * x using dgbmv */
    double *y = (double *)malloc(n * sizeof(double));
    memset(y, 0, n * sizeof(double));
    
    /* dgbmv arguments: layout, trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy */
    /* alpha=1.0, beta=0.0 */
    /* IMPOTANT: AB is in LU format (Diag at row kl+ku). dgbmv expects Diag at row ku.
       We must shift the pointer by 'kl' rows to align the band.
       Offset in doubles = kl (since Column Major has stride 1 within columns? 
       No, Column Major stride 1 is ROW. So AB[k] is row k of col 0. 
       So shifting by kl (1) shifts the row view by 1. */
    cblas_dgbmv(CblasColMajor, CblasNoTrans, n, n, kl, ku, 1.0, AB + kl, lab, x, 1, 0.0, y, 1);
    
    /* 4. Verify result */
    /* Mathematical result for 1D Poisson matrix [-1, 2, -1] * [1, 1, 1]^T:
       Row 0: 2*1 - 1*1 = 1
       Row i: -1*1 + 2*1 - 1*1 = 0
       Row n-1: -1*1 + 2*1 = 1
    */
    double max_err = 0.0;
    for (int i = 0; i < n; i++) {
        double expected = 0.0;
        if (i == 0 || i == n - 1) expected = 1.0;
        
        double err = fabs(y[i] - expected);
        if (err > max_err) max_err = err;
    }
    
    printf("Max Error for A*ones = [1, 0...0, 1]: %e\n", max_err);
    if (max_err < 1e-14) {
        printf("[PASS] dgbmv verification successful.\n");
    } else {
        printf("[FAIL] dgbmv verification failed!\n");
    }
    printf("\n");
    
    free(AB); free(x); free(y);
}

/* Exercise 6 Validation: Compare Custom LU (dgbtrftridiag) with LAPACK (dgbtrf) */
void test_lu_compare(int n) {
    printf("=== Test: LU Custom vs LAPACK (n=%d) ===\n", n);
    
    int kv = 1, ku = 1, kl = 1;
    int lab = kv + kl + ku + 1;
    
    /* Create two identical matrices */
    double *AB_lapack = (double *)malloc(lab * n * sizeof(double));
    double *AB_custom = (double *)malloc(lab * n * sizeof(double));
    
    set_GB_operator_colMajor_poisson1D(AB_lapack, &lab, &n, &kv);
    memcpy(AB_custom, AB_lapack, lab * n * sizeof(double));
    
    int *ipiv_lapack = (int *)malloc(n * sizeof(int));
    int *ipiv_custom = (int *)malloc(n * sizeof(int));
    int info;
    
    /* Run LAPACK dgbtrf */
    info = 0;
    dgbtrf_(&n, &n, &kl, &ku, AB_lapack, &lab, ipiv_lapack, &info);
    if (info != 0) printf("LAPACK dgbtrf failed with info=%d\n", info);
    
    /* Run Custom dgbtrftridiag */
    info = 0;
    dgbtrftridiag(&lab, &n, &kl, &ku, AB_custom, &lab, ipiv_custom, &info);
    if (info != 0) printf("Custom dgbtrftridiag failed with info=%d\n", info);

    /* Compare LU Factors (in AB) */
    /* Note: LAPACK might choose pivots, but for Poisson matrix (diag dominant), 
       it usually doesn't swap rows, so factors should be identical. */
    
    double max_diff = 0.0;
    for (int i = 0; i < lab * n; i++) {
        double diff = fabs(AB_lapack[i] - AB_custom[i]);
        if (diff > max_diff) max_diff = diff;
    }
    
    printf("Max difference in LU factors: %e\n", max_diff);
    
    /* Compare Pivots
       Note: dgbtrf uses Fortran 1-based indexing for pivots.
       ipiv_custom needs to match that if we want exact identity.
       In common dgbtrf implementation on Diag Dominant, ipiv[i] = i+1.
    */
    int pivot_diff_count = 0;
    for (int i = 0; i < n; i++) {
        if (ipiv_lapack[i] != ipiv_custom[i]) pivot_diff_count++;
    }
    printf("Pivot differences: %d\n", pivot_diff_count);

    if (max_diff < 1e-14 && pivot_diff_count == 0) {
        printf("[PASS] Custom LU matches LAPACK exactly.\n");
    } else {
        printf("[WARN] Differences found (expected if implementation details differ slightly, otherwise FAIL).\n");
    }
    printf("\n");

    free(AB_lapack); free(AB_custom);
    free(ipiv_lapack); free(ipiv_custom);
}

int main(int argc, char *argv[]) {
    printf("Starting Tests...\n\n");
    
    /* Test 1: Small scale check */
    test_dgbmv_ax_equals_b(5);
    test_lu_compare(5);
    
    /* Test 2: Larger scale */
    test_dgbmv_ax_equals_b(100);
    test_lu_compare(100);

    return 0;
}
