/*
 * pdsyrk_i8.c -- I8 entry point for PDSYRK.
 */
#include "pblas_i8_utils.h"
#include "PBblas.h"

extern void pdsyrk_(F_CHAR_T, F_CHAR_T, Int *, Int *, double *,
                     double *, Int *, Int *, Int *,
                     double *,
                     double *, Int *, Int *, Int *);

void pdsyrk_i8_(F_CHAR_T UPLO, F_CHAR_T TRANS,
                 int64_t *N, int64_t *K, double *ALPHA,
                 double *A, int64_t *IA, int64_t *JA, int64_t *DESCA,
                 double *BETA,
                 double *C, int64_t *IC, int64_t *JC, int64_t *DESCC)
{
    Int n4, k4, ia4, ja4, ic4, jc4;
    Int desca4[9], descc4[9];


    if (!pblas_i8_narrow(*N, &n4) || !pblas_i8_narrow(*K, &k4) ||
        !pblas_i8_narrow(*IA, &ia4) || !pblas_i8_narrow(*JA, &ja4) ||
        !pblas_i8_narrow(*IC, &ic4) || !pblas_i8_narrow(*JC, &jc4) ||
        !pblas_i8_narrow_desc(DESCA, desca4) ||
        !pblas_i8_narrow_desc(DESCC, descc4))
    {
        pblas_i8_abort((Int)DESCA[1], "PDSYRK_I8");
        return;
    }

    pdsyrk_(UPLO, TRANS, &n4, &k4, ALPHA,
            A, &ia4, &ja4, desca4,
            BETA,
            C, &ic4, &jc4, descc4);
}
