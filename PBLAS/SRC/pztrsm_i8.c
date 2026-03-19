/*
 * pztrsm_i8.c -- I8 entry point for PZTRSM.
 */
#include "pblas_i8_utils.h"
#include "PBblas.h"

extern void pztrsm_(F_CHAR_T, F_CHAR_T, F_CHAR_T, F_CHAR_T,
                     Int *, Int *, double *,
                     double *, Int *, Int *, Int *,
                     double *, Int *, Int *, Int *);

void pztrsm_i8_(F_CHAR_T SIDE, F_CHAR_T UPLO, F_CHAR_T TRANS, F_CHAR_T DIAG,
                 int64_t *M, int64_t *N, double *ALPHA,
                 double *A, int64_t *IA, int64_t *JA, int64_t *DESCA,
                 double *B, int64_t *IB, int64_t *JB, int64_t *DESCB)
{
    Int m4, n4, ia4, ja4, ib4, jb4;
    Int desca4[9], descb4[9];


    if (!pblas_i8_narrow(*M, &m4) || !pblas_i8_narrow(*N, &n4) ||
        !pblas_i8_narrow(*IA, &ia4) || !pblas_i8_narrow(*JA, &ja4) ||
        !pblas_i8_narrow(*IB, &ib4) || !pblas_i8_narrow(*JB, &jb4) ||
        !pblas_i8_narrow_desc(DESCA, desca4) ||
        !pblas_i8_narrow_desc(DESCB, descb4))
    {
        pblas_i8_abort((Int)DESCA[1], "PZTRSM_I8");
        return;
    }

    pztrsm_(SIDE, UPLO, TRANS, DIAG, &m4, &n4, ALPHA,
            A, &ia4, &ja4, desca4,
            B, &ib4, &jb4, descb4);
}
