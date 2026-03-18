/*
 * pzgemm_i8.c -- I8 entry point for PZGEMM.
 */
#include "pblas_i8_utils.h"
#include "PBblas.h"

extern void pzgemm_(F_CHAR_T, F_CHAR_T, Int *, Int *, Int *, double *,
                    double *, Int *, Int *, Int *,
                    double *, Int *, Int *, Int *,
                    double *,
                    double *, Int *, Int *, Int *);

void pzgemm_i8_(F_CHAR_T TRANSA, F_CHAR_T TRANSB,
                int64_t *M, int64_t *N, int64_t *K, double *ALPHA,
                double *A, int64_t *IA, int64_t *JA, int64_t *DESCA,
                double *B, int64_t *IB, int64_t *JB, int64_t *DESCB,
                double *BETA,
                double *C, int64_t *IC, int64_t *JC, int64_t *DESCC)
{
    Int m4, n4, k4, ia4, ja4, ib4, jb4, ic4, jc4;
    Int desca4[9], descb4[9], descc4[9];

    if (*M <= 0 || *N <= 0 || *K <= 0) return;

    if (!pblas_i8_narrow(*M, &m4) || !pblas_i8_narrow(*N, &n4) ||
        !pblas_i8_narrow(*K, &k4) ||
        !pblas_i8_narrow(*IA, &ia4) || !pblas_i8_narrow(*JA, &ja4) ||
        !pblas_i8_narrow(*IB, &ib4) || !pblas_i8_narrow(*JB, &jb4) ||
        !pblas_i8_narrow(*IC, &ic4) || !pblas_i8_narrow(*JC, &jc4) ||
        !pblas_i8_narrow_desc(DESCA, desca4) ||
        !pblas_i8_narrow_desc(DESCB, descb4) ||
        !pblas_i8_narrow_desc(DESCC, descc4))
    {
        pblas_i8_abort((Int)DESCA[1], "PZGEMM_I8");
        return;
    }

    pzgemm_(TRANSA, TRANSB, &m4, &n4, &k4, ALPHA,
            A, &ia4, &ja4, desca4,
            B, &ib4, &jb4, descb4,
            BETA,
            C, &ic4, &jc4, descc4);
}
