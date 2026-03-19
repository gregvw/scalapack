/*
 * pdger_i8.c -- I8 entry point for PDGER.
 */
#include "pblas_i8_utils.h"
#include "PBblas.h"

extern void pdger_(Int *, Int *, double *,
                   double *, Int *, Int *, Int *, Int *,
                   double *, Int *, Int *, Int *, Int *,
                   double *, Int *, Int *, Int *);

void pdger_i8_(int64_t *M, int64_t *N, double *ALPHA,
               double *X, int64_t *IX, int64_t *JX, int64_t *DESCX,
               int64_t *INCX,
               double *Y, int64_t *IY, int64_t *JY, int64_t *DESCY,
               int64_t *INCY,
               double *A, int64_t *IA, int64_t *JA, int64_t *DESCA)
{
    Int m4, n4, ix4, jx4, incx4, iy4, jy4, incy4, ia4, ja4;
    Int descx4[9], descy4[9], desca4[9];

    if (!pblas_i8_narrow(*M, &m4) || !pblas_i8_narrow(*N, &n4) ||
        !pblas_i8_narrow(*IX, &ix4) || !pblas_i8_narrow(*JX, &jx4) ||
        !pblas_i8_narrow(*INCX, &incx4) ||
        !pblas_i8_narrow(*IY, &iy4) || !pblas_i8_narrow(*JY, &jy4) ||
        !pblas_i8_narrow(*INCY, &incy4) ||
        !pblas_i8_narrow(*IA, &ia4) || !pblas_i8_narrow(*JA, &ja4) ||
        !pblas_i8_narrow_desc(DESCX, descx4) ||
        !pblas_i8_narrow_desc(DESCY, descy4) ||
        !pblas_i8_narrow_desc(DESCA, desca4))
    {
        pblas_i8_abort((Int)DESCA[1], "PDGER_I8");
        return;
    }

    pdger_(&m4, &n4, ALPHA,
           X, &ix4, &jx4, descx4, &incx4,
           Y, &iy4, &jy4, descy4, &incy4,
           A, &ia4, &ja4, desca4);
}
