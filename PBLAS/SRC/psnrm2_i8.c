/*
 * psnrm2_i8.c -- I8 entry point for PSNRM2.
 */

#include "pblas_i8_utils.h"

extern void psnrm2_(Int *, float *,
                    float *, Int *, Int *, Int *, Int *);

void psnrm2_i8_(int64_t *N, float *NORM2,
                float *X, int64_t *IX, int64_t *JX, int64_t *DESCX,
                int64_t *INCX)
{
    Int n4, ix4, jx4, incx4;
    Int descx4[9];


    if (!pblas_i8_narrow(*N, &n4) ||
        !pblas_i8_narrow(*IX, &ix4) ||
        !pblas_i8_narrow(*JX, &jx4) ||
        !pblas_i8_narrow(*INCX, &incx4) ||
        !pblas_i8_narrow_desc(DESCX, descx4))
    {
        pblas_i8_abort((Int)DESCX[1], "PSNRM2_I8");
        return;
    }

    psnrm2_(&n4, NORM2, X, &ix4, &jx4, descx4, &incx4);
}
