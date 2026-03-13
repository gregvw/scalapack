/*
 * pgemraux_core.c — Shared core helpers for PxGEMR2D
 *
 * All functions operate on MDESC_CORE / IDESC_CORE with
 * ScaLAPACK_Index64 dimensions.  Grid-coordinate values stay int.
 *
 * These are the widened equivalents of the helpers in pgemraux.c;
 * the originals remain for the legacy Int-based entry points.
 */

#include "redist_core.h"
#include "redist.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/* MPI datatype for ScaLAPACK_Index64 (= int64_t) */
#ifdef MPI_INT64_T
#define REDIST_MPI_I64 MPI_INT64_T
#else
/* Fallback: only valid when long long is exactly int64_t */
_Static_assert(sizeof(long long) == sizeof(int64_t),
               "MPI_LONG_LONG fallback requires long long == int64_t");
#define REDIST_MPI_I64 MPI_LONG_LONG
#endif

/* BLACS C-interface declarations (narrow extern, no private headers) */
extern void Cblacs_gridinfo(Int context, Int *nprow, Int *npcol,
                            Int *myrow, Int *mycol);
extern MPI_Comm Cblacs2sys_handle(Int BlacsCtxt);
extern void Cblacs_get(Int context, Int what, Int *val);
extern void Cigesd2d(Int context, Int m, Int n, Int *A, Int lda,
                     Int rdest, Int cdest);
extern void Cigerv2d(Int context, Int m, Int n, Int *A, Int lda,
                     Int rsrc, Int csrc);

#define SHIFT(row, sprow, nbrow) \
    ((row) - (sprow) + ((row) >= (sprow) ? 0 : (nbrow)))

/* ------------------------------------------------------------------ */
/* Checked-arithmetic convenience wrappers (fatal on overflow)        */
/* ------------------------------------------------------------------ */

static ScaLAPACK_Index64
i64_mul(ScaLAPACK_Index64 a, ScaLAPACK_Index64 b)
{
    ScaLAPACK_Index64 r;
    if (!ScaLAPACK_Index64Mul(a, b, &r)) {
        fprintf(stderr, "xxGEMR2D_CORE:i64 multiply overflow\n");
        exit(1);
    }
    return r;
}

static ScaLAPACK_Index64
i64_add(ScaLAPACK_Index64 a, ScaLAPACK_Index64 b)
{
    ScaLAPACK_Index64 r;
    if (!ScaLAPACK_Index64Add(a, b, &r)) {
        fprintf(stderr, "xxGEMR2D_CORE:i64 add overflow\n");
        exit(1);
    }
    return r;
}

static ScaLAPACK_Index64
i64_sub(ScaLAPACK_Index64 a, ScaLAPACK_Index64 b)
{
    /* b is always a non-negative dimension/index, so -b is valid */
    if (b < 0) {
        fprintf(stderr, "xxGEMR2D_CORE:i64 subtract negative operand\n");
        exit(1);
    }
    return i64_add(a, -b);
}

/* ------------------------------------------------------------------ */
/* localsize_core                                                     */
/* ------------------------------------------------------------------ */

ScaLAPACK_Index64
localsize_core(int myprow, int p,
               ScaLAPACK_Index64 nb, ScaLAPACK_Index64 m)
{
    ScaLAPACK_Index64 templateheight, blockheight, remainder;
    ScaLAPACK_Index64 p64 = (ScaLAPACK_Index64)p;
    ScaLAPACK_Index64 myprow64 = (ScaLAPACK_Index64)myprow;

    templateheight = i64_mul(p64, nb);

    remainder = m % templateheight;
    if (remainder != 0) {
        if (remainder > i64_mul(nb, myprow64)) {
            if (remainder >= i64_mul(nb, myprow64 + 1)) {
                blockheight = i64_add(i64_mul(m / templateheight, nb), nb);
            } else {
                blockheight = i64_add(i64_mul(m / templateheight, nb),
                                      m % nb);
            }
        } else {
            blockheight = i64_mul(m / templateheight, nb);
        }
    } else {
        blockheight = m / p64;
    }
    return blockheight;
}

/* ------------------------------------------------------------------ */
/* memoryblocksize_core                                               */
/* ------------------------------------------------------------------ */

ScaLAPACK_Index64
memoryblocksize_core(const MDESC_CORE *a)
{
    Int nprow_b, npcol_b, myprow_b, mypcol_b;
    int myprow, mypcol, p, q;

    Cblacs_gridinfo((Int)a->ctxt, &nprow_b, &npcol_b, &myprow_b, &mypcol_b);
    p = (int)nprow_b;
    q = (int)npcol_b;
    myprow = (int)myprow_b;
    mypcol = (int)mypcol_b;

    myprow = SHIFT(myprow, a->rsrc, p);
    mypcol = SHIFT(mypcol, a->csrc, q);
    assert(myprow >= 0 && mypcol >= 0);

    return i64_mul(localsize_core(myprow, p, a->mb, a->m),
                   localsize_core(mypcol, q, a->nb, a->n));
}

/* ------------------------------------------------------------------ */
/* changeorigin_core                                                  */
/* ------------------------------------------------------------------ */

ScaLAPACK_Index64
changeorigin_core(int myp, int sp, int p,
                  ScaLAPACK_Index64 bs, ScaLAPACK_Index64 i,
                  ScaLAPACK_Index64 *decal, int *newsp)
{
    ScaLAPACK_Index64 tempheight, firsttemp;
    int firstblock;
    ScaLAPACK_Index64 p64 = (ScaLAPACK_Index64)p;

    tempheight = i64_mul(bs, p64);
    firsttemp = i / tempheight;
    firstblock = (int)((i / bs) % p64);   /* result < p, fits in int */
    *newsp = (sp + firstblock) % p;

    if (myp >= 0)
        *decal = i64_add(i64_mul(firsttemp, bs),
                         SHIFT(myp, sp, p) < firstblock ? bs : 0);
    else
        *decal = 0;

    return i % bs;
}

/* ------------------------------------------------------------------ */
/* paramcheck_core                                                    */
/* ------------------------------------------------------------------ */

void
paramcheck_core(const MDESC_CORE *a,
                ScaLAPACK_Index64 i, ScaLAPACK_Index64 j,
                ScaLAPACK_Index64 m, ScaLAPACK_Index64 n,
                int p, int q, int gcontext)
{
    Int nprow_b, npcol_b, myprow_b, mypcol_b;
    int myprow, mypcol, p2, q2;

    (void)gcontext;  /* debug-only checkequal omitted for now */

    Cblacs_gridinfo((Int)a->ctxt, &nprow_b, &npcol_b, &myprow_b, &mypcol_b);
    p2 = (int)nprow_b;
    q2 = (int)npcol_b;
    myprow = (int)myprow_b;
    mypcol = (int)mypcol_b;

    if (myprow >= p2 || mypcol >= q2)
        myprow = mypcol = -1;

    if ((myprow >= 0 || mypcol >= 0) && (p2 != p && q2 != q)) {
        fprintf(stderr, "??MR2D_CORE:incoherent p,q parameters\n");
        exit(1);
    }
    assert(myprow < p && mypcol < q);

    if (a->rsrc < 0 || a->rsrc >= p || a->csrc < 0 || a->csrc >= q) {
        fprintf(stderr, "??MR2D_CORE:Bad first processor coordinates\n");
        exit(1);
    }
    if (i < 0 || j < 0 || i + m > a->m || j + n > a->n) {
        fprintf(stderr, "??MR2D_CORE:Bad submatrix:i=%lld,j=%lld,"
                "m=%lld,n=%lld,M=%lld,N=%lld\n",
                (long long)i, (long long)j,
                (long long)m, (long long)n,
                (long long)a->m, (long long)a->n);
        exit(1);
    }
    if ((myprow >= 0 || mypcol >= 0) &&
        localsize_core(SHIFT(myprow, a->rsrc, p), p, a->mb, a->m) > a->lld) {
        fprintf(stderr, "??MR2D_CORE:bad lda arg:row=%d,m=%lld,p=%d,"
                "mb=%lld,lld=%lld,rsrc=%d\n",
                myprow, (long long)a->m, p,
                (long long)a->mb, (long long)a->lld, a->rsrc);
        exit(1);
    }
}

/* ------------------------------------------------------------------ */
/* localindice_core                                                   */
/* ------------------------------------------------------------------ */

ScaLAPACK_Index64
localindice_core(ScaLAPACK_Index64 ig, ScaLAPACK_Index64 jg,
                 ScaLAPACK_Index64 templateheight,
                 ScaLAPACK_Index64 templatewidth,
                 const MDESC_CORE *a)
{
    ScaLAPACK_Index64 vtemp, htemp, vsubtemp, hsubtemp, il, jl;

    assert(ig >= 0 && ig < a->m && jg >= 0 && jg < a->n);

    vtemp = ig / templateheight;
    htemp = jg / templatewidth;
    vsubtemp = ig % a->mb;
    hsubtemp = jg % a->nb;

    il = i64_add(i64_mul(a->mb, vtemp), vsubtemp);
    jl = i64_add(i64_mul(a->nb, htemp), hsubtemp);

    assert(il < a->lld);

    return i64_add(i64_mul(jl, a->lld), il);
}

/* ------------------------------------------------------------------ */
/* scan_intervals_core                                                */
/* ------------------------------------------------------------------ */

ScaLAPACK_Index64
scan_intervals_core(char type,
                    ScaLAPACK_Index64 ja, ScaLAPACK_Index64 jb,
                    ScaLAPACK_Index64 n,
                    const MDESC_CORE *ma, const MDESC_CORE *mb,
                    int q0, int q1, int col0, int col1,
                    IDESC_CORE *result)
{
    ScaLAPACK_Index64 offset, j0, j1, templatewidth0, templatewidth1;
    ScaLAPACK_Index64 nbcol0, nbcol1;
    ScaLAPACK_Index64 l, end0, end1;

    assert(type == 'c' || type == 'r');

    nbcol0 = (type == 'c' ? ma->nb : ma->mb);
    nbcol1 = (type == 'c' ? mb->nb : mb->mb);
    templatewidth0 = i64_mul((ScaLAPACK_Index64)q0, nbcol0);
    templatewidth1 = i64_mul((ScaLAPACK_Index64)q1, nbcol1);

    {
        int sp0 = (type == 'c' ? ma->csrc : ma->rsrc);
        int sp1 = (type == 'c' ? mb->csrc : mb->rsrc);
        ScaLAPACK_Index64 shifted0 = i64_mul(
            (ScaLAPACK_Index64)SHIFT(col0, sp0, q0), nbcol0);
        ScaLAPACK_Index64 shifted1 = i64_mul(
            (ScaLAPACK_Index64)SHIFT(col1, sp1, q1), nbcol1);
        j0 = i64_sub(shifted0, ja);
        j1 = i64_sub(shifted1, jb);
    }

    offset = 0;
    l = 0;

    assert(j0 + nbcol0 > 0);
    assert(j1 + nbcol1 > 0);

    while ((j0 < n) && (j1 < n)) {
        ScaLAPACK_Index64 start, end;
        end0 = i64_add(j0, nbcol0);
        end1 = i64_add(j1, nbcol1);

        if (end0 <= j1) {
            j0 = i64_add(j0, templatewidth0);
            l = i64_add(l, nbcol0);
            continue;
        }
        if (end1 <= j0) {
            j1 = i64_add(j1, templatewidth1);
            continue;
        }

        /* raw intersection */
        start = j0 > j1 ? j0 : j1;
        if (start < 0) start = 0;

        /* lstart = l + (start - j0) */
        result[offset].lstart = i64_add(l, i64_sub(start, j0));

        end = end0 < end1 ? end0 : end1;
        if (end0 == end) {
            j0 = i64_add(j0, templatewidth0);
            l = i64_add(l, nbcol0);
        }
        if (end1 == end)
            j1 = i64_add(j1, templatewidth1);

        /* clamp to submatrix extent */
        if (end > n) end = n;
        assert(end > start);

        result[offset].len = i64_sub(end, start);
        offset = i64_add(offset, 1);
    }
    return offset;
}

/* ------------------------------------------------------------------ */
/* redist_sync_params_i8                                              */
/* ------------------------------------------------------------------ */

void
redist_sync_params_i8(int gcontext, ScaLAPACK_Index64 *param64, int nparam)
{
    Int sys_handle;
    MPI_Comm comm;

    /* Extract MPI communicator from the BLACS context.
     * Cblacs_get(ctx, 10, &h) stores ctx's p2p comm in the system
     * context table and returns the handle; Cblacs2sys_handle maps
     * that handle back to the MPI_Comm.  Duplicate-safe: the same
     * comm is not re-registered. */
    Cblacs_get((Int)gcontext, (Int)10, &sys_handle);
    comm = Cblacs2sys_handle(sys_handle);

    MPI_Allreduce(MPI_IN_PLACE, param64, nparam,
                  REDIST_MPI_I64, MPI_MIN, comm);
}
