/*
 * Tema 2 ASC
 * 2026 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double *B, double *x) {
	int i, j, k;

    double *C = (double *)malloc(N * N * sizeof(double));
    double *D = (double *)malloc(N * N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *tmp = (double *)malloc(N * sizeof(double));

    if (!C || !D || !y || !tmp) {
        free(C);
        free(D);
        free(y);
        free(tmp);
        return NULL;
    }

    for (i = 0; i < N; i++) {
        y[i] = 0.0;
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            double sum = 0.0;
            for (k = 0; k < N; k++) {
                sum += A[k * N + i] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }

    for (i = 0; i < N; i++) {
        for (j = i; j < N; j++) {
            double sum = 0.0;
            for (k = 0; k < N; k++) {
                sum += C[i * N + k] * C[j * N + k];
            }
            D[i * N + j] = sum;
            D[j * N + i] = sum;
        }
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            tmp[j] = 0.0;
        }

        for (j = 0; j < N; j++) {
            double vj = C[i * N + j];
            tmp[j] += D[j * N + j] * vj;

            for (k = j + 1; k < N; k++) {
                double d = D[j * N + k];
                tmp[j] += d * C[i * N + k];
                tmp[k] += d * vj;
            }
        }

        for (j = 0; j < N; j++) {
            y[j] += tmp[j];
        }
    }

    for (i = 0; i < N; i++) {
        y[i] += x[i];
    }

    free(C);
    free(D);
    free(tmp);

    return y;
}
