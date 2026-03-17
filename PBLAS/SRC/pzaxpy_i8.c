/*
 * pzaxpy_i8.c — I8 entry point for PZAXPY.
 */

#include "pblas_i8_utils.h"

typedef struct { double r, i; } cmplx8;

extern void pzaxpy_(Int *, cmplx8 *,
                    cmplx8 *, Int *, Int *, Int *, Int *,
                    cmplx8 *, Int *, Int *, Int *, Int *);

void pzaxpy_i8_(int64_t *N, cmplx8 *ALPHA,
                cmplx8 *X, int64_t *IX, int64_t *JX, int64_t *DESCX,
                int64_t *INCX,
                cmplx8 *Y, int64_t *IY, int64_t *JY, int64_t *DESCY,
                int64_t *INCY)
{
    Int n4, ix4, jx4, incx4, iy4, jy4, incy4;
    Int descx4[9], descy4[9];

    if (*N <= 0) return;

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
        pblas_i8_abort((Int)DESCX[1], "PZAXPY_I8");
        return;
    }

    pzaxpy_(&n4, ALPHA,
            X, &ix4, &jx4, descx4, &incx4,
            Y, &iy4, &jy4, descy4, &incy4);
}
