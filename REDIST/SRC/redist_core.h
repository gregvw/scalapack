/*
 * redist_core.h — Internal core descriptors for PxGEMR2D / PxTRMR2D
 *
 * MDESC_CORE is the single internal representation used by the
 * redistribution core.  All dimension fields are ScaLAPACK_Index64;
 * grid-coordinate fields stay int (process counts are always small).
 *
 * Two unpackers convert ABI-facing descriptors into MDESC_CORE:
 *   unpack_desc_i4  — from ScaLAPACK_ApiInt[9]  (standard descriptors)
 *   unpack_desc_i8  — from int64_t[9]            (I8 descriptors)
 */

#ifndef REDIST_CORE_H
#define REDIST_CORE_H 1

#include "scalapack-types.h"
#include <limits.h>

/* I8 descriptor type tag — matches SL_i8_params.inc BLOCK_CYCLIC_2D_I8.
 * Uses the same semantic value (1 = 2D block-cyclic) as legacy descriptors.
 * Descriptor width is determined by the entry-point name, not by DTYPE. */
#define BLOCK_CYCLIC_2D_I8 1

/*
 * Internal core descriptor.
 *
 * Field names follow the standard ScaLAPACK descriptor convention:
 *   DTYPE_, CTXT_, M_, N_, MB_, NB_, RSRC_, CSRC_, LLD_
 *
 * Widened fields (can exceed INT_MAX for I8 callers):
 *   m, n, mb, nb, lld
 *
 * Grid-coordinate fields (always bounded by process grid size):
 *   desctype, ctxt, rsrc, csrc
 */
typedef struct {
    int                desctype;   /* descriptor type (1=2D block-cyclic) */
    int                ctxt;       /* BLACS context handle               */
    ScaLAPACK_Index64  m;          /* global number of rows              */
    ScaLAPACK_Index64  n;          /* global number of columns           */
    ScaLAPACK_Index64  mb;         /* row block size                     */
    ScaLAPACK_Index64  nb;         /* column block size                  */
    int                rsrc;       /* starting process row               */
    int                csrc;       /* starting process column            */
    ScaLAPACK_Index64  lld;        /* local leading dimension            */
} MDESC_CORE;

/*
 * Internal core interval descriptor.
 *
 * Used by scan_intervals_core to describe local block intersections.
 * Both fields are Index64 because local matrix dimensions can exceed
 * INT_MAX when global dimensions are I8-scale.
 */
typedef struct {
    ScaLAPACK_Index64  lstart;     /* starting index in local memory      */
    ScaLAPACK_Index64  len;        /* length of this interval             */
} IDESC_CORE;

/*
 * unpack_desc_i4 — unpack a standard (32-bit integer) descriptor.
 *
 * All conversions are widening (ApiInt → Index64, ApiInt → int),
 * so this always succeeds.  No return value needed.
 */
static inline void
unpack_desc_i4(const ScaLAPACK_ApiInt desc[9], MDESC_CORE *out)
{
    out->desctype = (int)desc[0];
    out->ctxt     = (int)desc[1];
    out->m        = (ScaLAPACK_Index64)desc[2];
    out->n        = (ScaLAPACK_Index64)desc[3];
    out->mb       = (ScaLAPACK_Index64)desc[4];
    out->nb       = (ScaLAPACK_Index64)desc[5];
    out->rsrc     = (int)desc[6];
    out->csrc     = (int)desc[7];
    out->lld      = (ScaLAPACK_Index64)desc[8];
}

/*
 * unpack_desc_i8 — unpack a 64-bit integer descriptor.
 *
 * Dimension fields (m, n, mb, nb, lld) copy directly.
 * Grid-coordinate fields (desctype, ctxt, rsrc, csrc) are narrowed
 * to int with bounds checking.
 *
 * Returns 1 on success, 0 if any narrowing overflows.
 */
static inline int
unpack_desc_i8(const int64_t desc[9], MDESC_CORE *out)
{
    /* desctype: must fit in int */
    if (desc[0] < INT_MIN || desc[0] > INT_MAX) return 0;
    out->desctype = (int)desc[0];

    /* ctxt: BLACS context handle, must fit in int */
    if (desc[1] < INT_MIN || desc[1] > INT_MAX) return 0;
    out->ctxt = (int)desc[1];

    /* dimension fields: direct copy */
    out->m   = desc[2];
    out->n   = desc[3];
    out->mb  = desc[4];
    out->nb  = desc[5];

    /* rsrc: process row, must fit in int */
    if (desc[6] < INT_MIN || desc[6] > INT_MAX) return 0;
    out->rsrc = (int)desc[6];

    /* csrc: process column, must fit in int */
    if (desc[7] < INT_MIN || desc[7] > INT_MAX) return 0;
    out->csrc = (int)desc[7];

    /* lld: direct copy */
    out->lld = desc[8];

    return 1;
}

/*
 * Interval descriptor for triangular redistribution (trmr).
 *
 * Uses gstart (global start within the submatrix extent) instead of
 * lstart (local memory offset) because scanD0 works in global
 * coordinates and converts to local indices via localindice_core.
 */
typedef struct {
    ScaLAPACK_Index64  gstart;     /* global start position in submatrix  */
    ScaLAPACK_Index64  len;        /* length of this interval             */
} IDESC_TR_CORE;

/* ------------------------------------------------------------------ */
/* Core helper function declarations (implemented in pgemraux_core.c) */
/* ------------------------------------------------------------------ */

/* Number of local rows (or columns) owned by process myprow in a
 * 1D block-cyclic distribution of m elements with block size nb
 * over p processes. */
extern ScaLAPACK_Index64
localsize_core(int myprow, int p, ScaLAPACK_Index64 nb, ScaLAPACK_Index64 m);

/* Total local block size (rows * cols) for the local portion of a. */
extern ScaLAPACK_Index64
memoryblocksize_core(const MDESC_CORE *a);

/* Shift the submatrix origin so that i < bs.  Updates *decal (local
 * element shift) and *newsp (new starting process coordinate).
 * Returns the residual index (i % bs). */
extern ScaLAPACK_Index64
changeorigin_core(int myp, int sp, int p,
                  ScaLAPACK_Index64 bs, ScaLAPACK_Index64 i,
                  ScaLAPACK_Index64 *decal, int *newsp);

/* Validate descriptor parameters against the process grid. */
extern void
paramcheck_core(const MDESC_CORE *a,
                ScaLAPACK_Index64 i, ScaLAPACK_Index64 j,
                ScaLAPACK_Index64 m, ScaLAPACK_Index64 n,
                int p, int q, int gcontext);

/* Local memory index of element (ig,jg) in a block-cyclic matrix. */
extern ScaLAPACK_Index64
localindice_core(ScaLAPACK_Index64 ig, ScaLAPACK_Index64 jg,
                 ScaLAPACK_Index64 templateheight,
                 ScaLAPACK_Index64 templatewidth,
                 const MDESC_CORE *a);

/* Scan two block-cyclic distributions and find intersections.
 * Returns the number of intersections (can exceed INT_MAX for
 * I8-scale matrices with small block sizes). */
extern ScaLAPACK_Index64
scan_intervals_core(char type,
                    ScaLAPACK_Index64 ja, ScaLAPACK_Index64 jb,
                    ScaLAPACK_Index64 n,
                    const MDESC_CORE *ma, const MDESC_CORE *mb,
                    int q0, int q1, int col0, int col1,
                    IDESC_CORE *result);

/* Synchronize a 64-bit parameter array across all processes in a
 * BLACS context using MPI_Allreduce with MPI_MIN.  Sentinel values
 * (MAGIC_MAX_I8) are replaced by the actual values from processes
 * that contribute them. */
extern void
redist_sync_params_i8(int gcontext,
                      ScaLAPACK_Index64 *param64, int nparam);

/* Scan two block-cyclic distributions for triangular redistribution.
 * Like scan_intervals_core but stores global start (gstart) instead
 * of local start (lstart).  Used by scanD0_core in trmr files. */
extern ScaLAPACK_Index64
scan_intervals_tr_core(char type,
                       ScaLAPACK_Index64 ja, ScaLAPACK_Index64 jb,
                       ScaLAPACK_Index64 n,
                       const MDESC_CORE *ma, const MDESC_CORE *mb,
                       int q0, int q1, int col0, int col1,
                       IDESC_TR_CORE *result);

/* Number of elements in a triangular column.
 * Returns the count of elements in column j of a trapezoid
 * (upper or lower, unit or non-unit diagonal) starting from row i.
 * *offset is set to the number of rows to skip before the first
 * triangle element (used by intersect_core). */
extern ScaLAPACK_Index64
insidemat_core(const char *uplo, const char *diag,
               ScaLAPACK_Index64 i, ScaLAPACK_Index64 j,
               ScaLAPACK_Index64 m, ScaLAPACK_Index64 n,
               ScaLAPACK_Index64 *offset);

#define NBPARAM_CORE     20
#define MAGIC_MAX_I8     INT64_C(100000000000000000) /* 10^17 */

#endif /* REDIST_CORE_H */
