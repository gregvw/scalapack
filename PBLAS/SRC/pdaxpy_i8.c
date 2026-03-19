/*
 * pdaxpy_i8.c — I8 entry point for PDAXPY.
 *
 * Takes int64_t arguments from Fortran INTEGER*8 callers, narrows
 * to Int with overflow checks, and delegates to the legacy pdaxpy_.
 */

#include "pblas_i8_utils.h"

/* Legacy entry point */
extern void pdaxpy_(Int *, double *,
                    double *, Int *, Int *, Int *, Int *,
                    double *, Int *, Int *, Int *, Int *);

void pdaxpy_i8_(int64_t *N, double *ALPHA,
                double *X, int64_t *IX, int64_t *JX, int64_t *DESCX,
                int64_t *INCX,
                double *Y, int64_t *IY, int64_t *JY, int64_t *DESCY,
                int64_t *INCY)
{
    Int n4, ix4, jx4, incx4, iy4, jy4, incy4;
    Int descx4[9], descy4[9];
    Int ctxt;

    /* Quick return */

    /* Narrow all integer args */
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
        /* Get context for error reporting */
        ctxt = (Int)DESCX[1];  /* CTXT_ field */
        pblas_i8_abort(ctxt, "PDAXPY_I8");
        return;
    }

    pdaxpy_(&n4, ALPHA,
            X, &ix4, &jx4, descx4, &incx4,
            Y, &iy4, &jy4, descy4, &incy4);
}
