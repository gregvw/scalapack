/*
 * pcgemv_i8.c -- I8 entry point for PCGEMV.
 */

#include "pblas_i8_utils.h"
#include "PBblas.h"

extern void pcgemv_(F_CHAR_T, Int *, Int *, float *,
                    float *, Int *, Int *, Int *,
                    float *, Int *, Int *, Int *, Int *,
                    float *,
                    float *, Int *, Int *, Int *, Int *);

void pcgemv_i8_(F_CHAR_T TRANS, int64_t *M, int64_t *N, float *ALPHA,
                float *A, int64_t *IA, int64_t *JA, int64_t *DESCA,
                float *X, int64_t *IX, int64_t *JX, int64_t *DESCX,
                int64_t *INCX,
                float *BETA,
                float *Y, int64_t *IY, int64_t *JY, int64_t *DESCY,
                int64_t *INCY)
{
    Int m4, n4, ia4, ja4, ix4, jx4, incx4, iy4, jy4, incy4;
    Int desca4[9], descx4[9], descy4[9];


    if (!pblas_i8_narrow(*M, &m4) ||
        !pblas_i8_narrow(*N, &n4) ||
        !pblas_i8_narrow(*IA, &ia4) ||
        !pblas_i8_narrow(*JA, &ja4) ||
        !pblas_i8_narrow(*IX, &ix4) ||
        !pblas_i8_narrow(*JX, &jx4) ||
        !pblas_i8_narrow(*INCX, &incx4) ||
        !pblas_i8_narrow(*IY, &iy4) ||
        !pblas_i8_narrow(*JY, &jy4) ||
        !pblas_i8_narrow(*INCY, &incy4) ||
        !pblas_i8_narrow_desc(DESCA, desca4) ||
        !pblas_i8_narrow_desc(DESCX, descx4) ||
        !pblas_i8_narrow_desc(DESCY, descy4))
    {
        pblas_i8_abort((Int)DESCA[1], "PCGEMV_I8");
        return;
    }

    pcgemv_(TRANS, &m4, &n4, ALPHA,
            A, &ia4, &ja4, desca4,
            X, &ix4, &jx4, descx4, &incx4,
            BETA,
            Y, &iy4, &jy4, descy4, &incy4);
}
