/*
 * Tema 2 ASC
 * 2026 Spring
 */

#include "utils.h"

#pragma GCC optimize("unroll-loops")
#pragma GCC target("avx2,fma")

#define BLOCK 40

double* my_solver(int N, double *A, double *B, double *x) {
    register int i, j, k, ii, jj, kk;

    double *AT = (double *)malloc((size_t)N * N * sizeof(double));
    double *C  = (double *)calloc((size_t)N * N, sizeof(double));
    double *CT = (double *)malloc((size_t)N * N * sizeof(double));
    double *y  = (double *)calloc((size_t)N, sizeof(double));
    double *t1 = (double *)malloc((size_t)N * sizeof(double));
    double *t2 = (double *)malloc((size_t)N * sizeof(double));

    if (!AT || !C || !CT || !y || !t1 || !t2) {
        free(AT);
        free(C);
        free(CT);
        free(y);
        free(t1);
        free(t2);
        return NULL;
    }

    /* C = A^T * B : pretranspunere A -> AT pentru acces liniar/cache-friendly */
    for (i = 0; i < N; i++) {
        register double *a_row = A + (size_t)i * N;
        for (j = 0; j + 3 < N; j += 4) {
            *(AT + (size_t)(j + 0) * N + i) = *(a_row + j + 0);
            *(AT + (size_t)(j + 1) * N + i) = *(a_row + j + 1);
            *(AT + (size_t)(j + 2) * N + i) = *(a_row + j + 2);
            *(AT + (size_t)(j + 3) * N + i) = *(a_row + j + 3);
        }
        for (; j < N; j++) {
            *(AT + (size_t)j * N + i) = *(a_row + j);
        }
    }

    /* C = A^T * B : inmultire blocked + pointeri + unrolling x4 */
    for (ii = 0; ii < N; ii += BLOCK) {
        register int i_end = (ii + BLOCK < N) ? (ii + BLOCK) : N;

        for (kk = 0; kk < N; kk += BLOCK) {
            register int k_end = (kk + BLOCK < N) ? (kk + BLOCK) : N;

            for (jj = 0; jj < N; jj += BLOCK) {
                register int j_end = (jj + BLOCK < N) ? (jj + BLOCK) : N;

                for (i = ii; i < i_end; i++) {
                    register double *c_row = C + (size_t)i * N + jj;
                    register double *at_row = AT + (size_t)i * N + kk;

                    for (k = kk; k < k_end; k++) {
                        register double aval = *at_row;
                        register const double *b_ptr = B + (size_t)k * N + jj;
                        register double *c_ptr = c_row;

                        #pragma GCC ivdep
                        for (j = jj; j + 3 < j_end; j += 4) {
                            *(c_ptr + 0) += aval * (*(b_ptr + 0));
                            *(c_ptr + 1) += aval * (*(b_ptr + 1));
                            *(c_ptr + 2) += aval * (*(b_ptr + 2));
                            *(c_ptr + 3) += aval * (*(b_ptr + 3));
                            c_ptr += 4;
                            b_ptr += 4;
                        }

                        for (; j < j_end; j++) {
                            *c_ptr += aval * (*b_ptr);
                            c_ptr++;
                            b_ptr++;
                        }

                        at_row++;
                    }
                }
            }
        }
    }

    for (i = 0; i < N; i++) {
        register double *c_row = C + (size_t)i * N;
        for (j = 0; j + 3 < N; j += 4) {
            *(CT + (size_t)(j + 0) * N + i) = *(c_row + j + 0);
            *(CT + (size_t)(j + 1) * N + i) = *(c_row + j + 1);
            *(CT + (size_t)(j + 2) * N + i) = *(c_row + j + 2);
            *(CT + (size_t)(j + 3) * N + i) = *(c_row + j + 3);
        }
        for (; j < N; j++) {
            *(CT + (size_t)j * N + i) = *(c_row + j);
        }
    }

    for (i = 0; i < N; i++) {
        register const double *row = C + (size_t)i * N;

        for (j = 0; j < N; j++) {
            register const double *m_ptr = CT + (size_t)j * N;
            register const double *v_ptr = row;
            register double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;

            #pragma GCC ivdep
            for (k = 0; k + 3 < N; k += 4) {
                s0 += (*(m_ptr + 0)) * (*(v_ptr + 0));
                s1 += (*(m_ptr + 1)) * (*(v_ptr + 1));
                s2 += (*(m_ptr + 2)) * (*(v_ptr + 2));
                s3 += (*(m_ptr + 3)) * (*(v_ptr + 3));
                m_ptr += 4;
                v_ptr += 4;
            }

            register double sum = s0 + s1 + s2 + s3;
            for (; k < N; k++) {
                sum += (*m_ptr) * (*v_ptr);
                m_ptr++;
                v_ptr++;
            }

            *(t1 + j) = sum;
        }

        for (j = 0; j < N; j++) {
            register const double *m_ptr = C + (size_t)j * N;
            register const double *v_ptr = t1;
            register double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;

            #pragma GCC ivdep
            for (k = 0; k + 3 < N; k += 4) {
                s0 += (*(m_ptr + 0)) * (*(v_ptr + 0));
                s1 += (*(m_ptr + 1)) * (*(v_ptr + 1));
                s2 += (*(m_ptr + 2)) * (*(v_ptr + 2));
                s3 += (*(m_ptr + 3)) * (*(v_ptr + 3));
                m_ptr += 4;
                v_ptr += 4;
            }

            register double sum = s0 + s1 + s2 + s3;
            for (; k < N; k++) {
                sum += (*m_ptr) * (*v_ptr);
                m_ptr++;
                v_ptr++;
            }

            *(t2 + j) = sum;
        }

        {
            register double *y_ptr = y;
            register double *t_ptr = t2;

            #pragma GCC ivdep
            for (j = 0; j + 3 < N; j += 4) {
                *(y_ptr + 0) += *(t_ptr + 0);
                *(y_ptr + 1) += *(t_ptr + 1);
                *(y_ptr + 2) += *(t_ptr + 2);
                *(y_ptr + 3) += *(t_ptr + 3);
                y_ptr += 4;
                t_ptr += 4;
            }

            for (; j < N; j++) {
                *y_ptr += *t_ptr;
                y_ptr++;
                t_ptr++;
            }
        }
    }

    {
        register double *y_ptr = y;
        register double *x_ptr = x;

        #pragma GCC ivdep
        for (i = 0; i + 3 < N; i += 4) {
            *(y_ptr + 0) += *(x_ptr + 0);
            *(y_ptr + 1) += *(x_ptr + 1);
            *(y_ptr + 2) += *(x_ptr + 2);
            *(y_ptr + 3) += *(x_ptr + 3);
            y_ptr += 4;
            x_ptr += 4;
        }

        for (; i < N; i++) {
            *y_ptr += *x_ptr;
            y_ptr++;
            x_ptr++;
        }
    }

    free(AT);
    free(C);
    free(CT);
    free(t1);
    free(t2);

    return y;
}
