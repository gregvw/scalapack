/*
 * pdamax_i8.c -- I8 entry point for PDAMAX.
 */
#include "pblas_i8_utils.h"
#include "PBblas.h"

extern void pdamax_(Int *, double *, Int *,
                    double *, Int *, Int *, Int *, Int *);

void pdamax_i8_(int64_t *N, double *AMAX, Int *INDX,
                double *X, int64_t *IX, int64_t *JX, int64_t *DESCX,
                int64_t *INCX)
{
    Int n4, ix4, jx4, incx4;
    Int descx4[9];

    if (!pblas_i8_narrow(*N, &n4) ||
        !pblas_i8_narrow(*IX, &ix4) || !pblas_i8_narrow(*JX, &jx4) ||
        !pblas_i8_narrow(*INCX, &incx4) ||
        !pblas_i8_narrow_desc(DESCX, descx4))
    {
        pblas_i8_abort((Int)DESCX[1], "PDAMAX_I8");
        return;
    }

    pdamax_(&n4, AMAX, INDX,
            X, &ix4, &jx4, descx4, &incx4);
}
