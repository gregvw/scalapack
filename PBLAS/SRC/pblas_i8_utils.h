/*
 * pblas_i8_utils.h — Shared utilities for PBLAS I8 entry points.
 *
 * Provides checked narrowing of int64_t arguments to Int (ScaLAPACK_ApiInt)
 * for delegation to legacy PBLAS internals.
 */

#ifndef PBLAS_I8_UTILS_H
#define PBLAS_I8_UTILS_H 1

#include "pblas.h"
#include "PBtools.h"
#include "PBblacs.h"
#include <stdint.h>
#include <limits.h>

/*
 * Narrow a single int64_t to Int.  Returns 0 on overflow, 1 on success.
 */
static inline int
pblas_i8_narrow(int64_t val, Int *out)
{
    if (val < (int64_t)INT_MIN || val > (int64_t)INT_MAX) return 0;
    *out = (Int)val;
    return 1;
}

/*
 * Narrow and convert an I8 ScaLAPACK descriptor (int64_t[9]) to the
 * legacy Int[9] format expected by PBLAS entry points.
 *
 * Dimension fields (M, N, MB, NB, LLD) are checked for overflow.
 * Grid-coordinate fields (DTYPE, CTXT, RSRC, CSRC) are always small.
 *
 * Returns 0 on overflow, 1 on success.
 */
static inline int
pblas_i8_narrow_desc(const int64_t desc_i8[9], Int desc_out[9])
{
    for (int k = 0; k < 9; k++) {
        if (!pblas_i8_narrow(desc_i8[k], &desc_out[k]))
            return 0;
    }
    return 1;
}

/*
 * Abort with PBLAS error reporting if narrowing fails.
 */
static inline void
pblas_i8_abort(Int ctxt, const char *routine)
{
    PB_Cabort(ctxt, routine, -2);
}

#endif /* PBLAS_I8_UTILS_H */
