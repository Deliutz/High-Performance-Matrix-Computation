/*
 * Tema 2 ASC
 * 2024 Spring
 */
#include "utils.h"
#include <cblas.h>
/* 
 * Add your BLAS implementation here
 */

double* my_solver(int N, double *A, double *B, double *x) {
    double *C = (double *)malloc((size_t)N * N * sizeof(double));
    double *D = (double *)calloc((size_t)N * N, sizeof(double));
    double *y = (double *)malloc((size_t)N * sizeof(double));
    double *s = (double *)malloc((size_t)N * sizeof(double));
    double *ones = (double *)malloc((size_t)N * sizeof(double));

    int i;

    if (!C || !D || !y || !s || !ones) {
        free(C);
        free(D);
        free(y);
        free(s);
        free(ones);
        return NULL;
    }

    for (i = 0; i < N; i++) {
        ones[i] = 1.0;
    }

    cblas_dgemm(
        CblasRowMajor, 
        CblasTrans,    
        CblasNoTrans,  
        N,             
        N,             
        N,             
        1.0,           
        A,             
        N,             
        B,             
        N,             
        0.0,           
        C,             
        N              
    );

    cblas_dsyrk(
        CblasRowMajor, 
        CblasUpper,    
        CblasNoTrans,  
        N,             
        N,             
        1.0,           
        C,             
        N,             
        0.0,           
        D,             
        N             
    );

    cblas_dgemv(
        CblasRowMajor, 
        CblasTrans,    
        N,             
        N,             
        1.0,           
        C,             
        N,             
        ones,          
        1,             
        0.0,           
        s,             
        1             
    );

    cblas_dsymv(
        CblasRowMajor,
        CblasUpper,    
        N,             
        1.0,           
        D,             
        N,             
        s,             
        1,             
        0.0,          
        y,             
        1             
    );

    cblas_daxpy(
        N,   
        1.0, 
        x,   
        1,   
        y,   
        1    
    );

    free(C);
    free(D);
    free(s);
    free(ones);

    return y;
}
