/*
 * pzgemv_i8.c -- I8 entry point for PZGEMV.
 */

#include "pblas_i8_utils.h"
#include "PBblas.h"

extern void pzgemv_(F_CHAR_T, Int *, Int *, double *,
                    double *, Int *, Int *, Int *,
                    double *, Int *, Int *, Int *, Int *,
                    double *,
                    double *, Int *, Int *, Int *, Int *);

void pzgemv_i8_(F_CHAR_T TRANS, int64_t *M, int64_t *N, double *ALPHA,
                double *A, int64_t *IA, int64_t *JA, int64_t *DESCA,
                double *X, int64_t *IX, int64_t *JX, int64_t *DESCX,
                int64_t *INCX,
                double *BETA,
                double *Y, int64_t *IY, int64_t *JY, int64_t *DESCY,
                int64_t *INCY)
{
    Int m4, n4, ia4, ja4, ix4, jx4, incx4, iy4, jy4, incy4;
    Int desca4[9], descx4[9], descy4[9];

    if (*M <= 0 || *N <= 0) return;

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
        pblas_i8_abort((Int)DESCA[1], "PZGEMV_I8");
        return;
    }

    pzgemv_(TRANS, &m4, &n4, ALPHA,
            A, &ia4, &ja4, desca4,
            X, &ix4, &jx4, descx4, &incx4,
            BETA,
            Y, &iy4, &jy4, descy4, &incy4);
}
