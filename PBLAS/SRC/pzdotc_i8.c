/*
 * pzdotc_i8.c -- I8 entry point for PZDOTC (conjugate dot product).
 */

#include "pblas_i8_utils.h"

extern void pzdotc_(Int *, double *,
                    double *, Int *, Int *, Int *, Int *,
                    double *, Int *, Int *, Int *, Int *);

void pzdotc_i8_(int64_t *N, double *DOT,
                double *X, int64_t *IX, int64_t *JX, int64_t *DESCX,
                int64_t *INCX,
                double *Y, int64_t *IY, int64_t *JY, int64_t *DESCY,
                int64_t *INCY)
{
    Int n4, ix4, jx4, incx4, iy4, jy4, incy4;
    Int descx4[9], descy4[9];


    if (!pblas_i8_narrow(*N, &n4) ||
        !pblas_i8_narrow(*IX, &ix4) ||
        !pblas_i8_narrow(*JX, &jx4) ||
        !pblas_i8_narrow(*INCX, &incx4) ||
        !pblas_i8_narrow(*IY, &iy4) ||
        !pblas_i8_narrow(*JY, &jy4) ||
        !pblas_i8_narrow(*INCY, &incy4) ||
        !pblas_i8_narrow_desc(DESCX, descx4) ||
        !pblas_i8_narrow_desc(DESCY, descy4))
    {
        pblas_i8_abort((Int)DESCX[1], "PZDOTC_I8");
        return;
    }

    pzdotc_(&n4, DOT,
            X, &ix4, &jx4, descx4, &incx4,
            Y, &iy4, &jy4, descy4, &incy4);
}
