#include "scalapack-types.h"
#include <stddef.h>

#ifdef T3D
#define float double
#endif
#ifdef T3E
#define float double
#endif
#ifdef CRAY
#define float double
#endif
#ifndef Int
#define Int ScaLAPACK_ApiInt
#endif

_Static_assert(sizeof(Int) == sizeof(ScaLAPACK_ApiInt),
               "Int must match the configured ScaLAPACK API integer width.");

static inline int ScaLAPACK_RedistApiMulToSizeT(Int left, Int right,
                                                size_t *result)
{
    ScaLAPACK_Index64 product;

    if (!ScaLAPACK_Index64Mul((ScaLAPACK_Index64) left,
                              (ScaLAPACK_Index64) right, &product))
        return 0;
    return ScaLAPACK_Index64ToSizeT(product, result);
}

static inline int ScaLAPACK_RedistApiAddToSizeT(Int left, Int right,
                                                size_t *result)
{
    ScaLAPACK_Index64 sum;

    if (!ScaLAPACK_Index64Add((ScaLAPACK_Index64) left,
                              (ScaLAPACK_Index64) right, &sum))
        return 0;
    return ScaLAPACK_Index64ToSizeT(sum, result);
}

static inline int ScaLAPACK_RedistDivUpToSizeT(Int numerator, Int denominator,
                                               size_t *result)
{
    ScaLAPACK_Index64 count;

    if (numerator < 0 || denominator <= 0) return 0;
    if (numerator == 0)
    {
        *result = 0;
        return 1;
    }
    count = (((ScaLAPACK_Index64) numerator) - 1) / denominator + 1;
    return ScaLAPACK_Index64ToSizeT(count, result);
}

static inline int ScaLAPACK_RedistApiMulToApiInt(Int left, Int right,
                                                 Int *result)
{
    ScaLAPACK_Index64 product;

    if (!ScaLAPACK_Index64Mul((ScaLAPACK_Index64) left,
                              (ScaLAPACK_Index64) right, &product))
        return 0;
    return ScaLAPACK_Index64ToApiInt(product, result);
}

static inline int ScaLAPACK_RedistApiAddToApiInt(Int left, Int right,
                                                 Int *result)
{
    ScaLAPACK_Index64 sum;

    if (!ScaLAPACK_Index64Add((ScaLAPACK_Index64) left,
                              (ScaLAPACK_Index64) right, &sum))
        return 0;
    return ScaLAPACK_Index64ToApiInt(sum, result);
}

static inline int ScaLAPACK_RedistApiSubToApiInt(Int left, Int right,
                                                 Int *result)
{
    ScaLAPACK_Index64 diff;

    if (!ScaLAPACK_Index64Add((ScaLAPACK_Index64) left,
                              -((ScaLAPACK_Index64) right), &diff))
        return 0;
    return ScaLAPACK_Index64ToApiInt(diff, result);
}
