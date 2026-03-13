/*
 * test_pdgemr2d_i8.c — Caller-side test for the pdgemr2d_i8_ entry point
 *
 * Two test cases:
 *
 *   Test 1 ("big-block descriptor"): M=10, N=5, MB=NB=3e9 (> INT_MAX).
 *     Same grid for source and destination.  All data on process (0,0).
 *     Exercises unpack_desc_i8, redist_sync_params_i8, and 64-bit
 *     arithmetic in the core helpers.
 *
 *   Test 2 ("cross-rank redistribution"): M=12, N=8, two different
 *     grids (2×2 and 1×4), small block sizes (MB=3,NB=2 → MB=4,NB=3).
 *     Data is distributed across ranks in grid 0, then redistributed
 *     to grid 1 with a different layout, then back.  This exercises
 *     scan_intervals_core, cross-rank send/recv, proc-map extraction,
 *     and reassembly.
 *
 * Usage:  mpirun -np 4 ./xdgemr_i8
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

/* BLACS C-interface declarations */
extern void Cblacs_pinfo(int *mypnum, int *nprocs);
extern void Cblacs_get(int context, int what, int *val);
extern void Cblacs_gridinit(int *context, const char *order,
                            int nprow, int npcol);
extern void Cblacs_gridinfo(int context, int *nprow, int *npcol,
                            int *myprow, int *mypcol);
extern void Cblacs_gridexit(int context);
extern void Cblacs_exit(int status);

/* I8 entry point — Fortran-mangled name */
#if defined(Add_) || defined(f77IsF2C)
#define pdgemr2d_i8 pdgemr2d_i8_
#elif defined(UpCase)
#define pdgemr2d_i8 PDGEMR2D_I8
#else
/* default: lower-case, no underscore */
#endif

extern void pdgemr2d_i8(int64_t *m, int64_t *n,
                         double *A, int64_t *ia, int64_t *ja,
                         int64_t desc_A[9],
                         double *B, int64_t *ib, int64_t *jb,
                         int64_t desc_B[9],
                         int64_t *gcontext);

/* BLOCK_CYCLIC_2D_I8 descriptor type tag */
#define BLOCK_CYCLIC_2D_I8 501

/* Descriptor field indices */
#define DTYPE_ 0
#define CTXT_  1
#define M_     2
#define N_     3
#define MB_    4
#define NB_    5
#define RSRC_  6
#define CSRC_  7
#define LLD_   8

/* ------------------------------------------------------------------ */
/* numroc — local row/column count for block-cyclic distribution      */
/* ------------------------------------------------------------------ */
static int64_t
numroc(int64_t n, int64_t nb, int iproc, int isrcproc, int nprocs)
{
    int64_t nblocks, mydist, extra;
    mydist = (iproc + nprocs - isrcproc) % nprocs;
    nblocks = n / nb;
    extra   = n % nb;
    return (nblocks / nprocs) * nb
         + (mydist < (nblocks % nprocs) ? nb : 0)
         + (mydist == (nblocks % nprocs) ? extra : 0);
}

/* ------------------------------------------------------------------ */
/* make_desc — fill an I8 descriptor array                            */
/* ------------------------------------------------------------------ */
static void
make_desc(int64_t desc[9], int ctxt,
          int64_t m, int64_t n, int64_t mb, int64_t nb,
          int rsrc, int csrc, int64_t lld)
{
    desc[DTYPE_] = BLOCK_CYCLIC_2D_I8;
    desc[CTXT_]  = (int64_t)ctxt;
    desc[M_]     = m;
    desc[N_]     = n;
    desc[MB_]    = mb;
    desc[NB_]    = nb;
    desc[RSRC_]  = (int64_t)rsrc;
    desc[CSRC_]  = (int64_t)csrc;
    desc[LLD_]   = lld;
}

/* ------------------------------------------------------------------ */
/* global_value — A(gi,gj) = gi*1000 + gj  (0-based global indices)  */
/* ------------------------------------------------------------------ */
static double
global_value(int gi, int gj)
{
    return (double)(gi * 1000 + gj);
}

/* ------------------------------------------------------------------ */
/* init_distributed — fill local block-cyclic storage with the global */
/*                    pattern A(gi,gj) = gi*1000 + gj                */
/* ------------------------------------------------------------------ */
static void
init_distributed(double *buf, int64_t lld,
                 int64_t m, int64_t n, int64_t mb, int64_t nb,
                 int myprow, int mypcol, int nprow, int npcol,
                 int rsrc, int csrc)
{
    int64_t lr = numroc(m, mb, myprow, rsrc, nprow);
    int64_t lc = numroc(n, nb, mypcol, csrc, npcol);
    for (int64_t lj = 0; lj < lc; lj++) {
        /* global column index for local column lj */
        int64_t bcol = lj / nb;                   /* which local block */
        int64_t boff = lj % nb;                   /* offset within block */
        int64_t gj = ((bcol * npcol + ((mypcol - csrc + npcol) % npcol))
                       * nb) + boff;
        for (int64_t li = 0; li < lr; li++) {
            int64_t brow = li / mb;
            int64_t roff = li % mb;
            int64_t gi = ((brow * nprow + ((myprow - rsrc + nprow) % nprow))
                           * mb) + roff;
            buf[lj * lld + li] = global_value((int)gi, (int)gj);
        }
    }
}

/* ------------------------------------------------------------------ */
/* check_distributed — verify local storage matches the global pattern*/
/* ------------------------------------------------------------------ */
static int
check_distributed(const double *buf, int64_t lld,
                  int64_t m, int64_t n, int64_t mb, int64_t nb,
                  int myprow, int mypcol, int nprow, int npcol,
                  int rsrc, int csrc, const char *label, int mypnum)
{
    int errs = 0;
    int64_t lr = numroc(m, mb, myprow, rsrc, nprow);
    int64_t lc = numroc(n, nb, mypcol, csrc, npcol);
    for (int64_t lj = 0; lj < lc; lj++) {
        int64_t bcol = lj / nb;
        int64_t boff = lj % nb;
        int64_t gj = ((bcol * npcol + ((mypcol - csrc + npcol) % npcol))
                       * nb) + boff;
        for (int64_t li = 0; li < lr; li++) {
            int64_t brow = li / mb;
            int64_t roff = li % mb;
            int64_t gi = ((brow * nprow + ((myprow - rsrc + nprow) % nprow))
                           * mb) + roff;
            double expected = global_value((int)gi, (int)gj);
            double actual   = buf[lj * lld + li];
            if (actual != expected) {
                if (errs < 5)
                    fprintf(stderr, "%s proc %d: MISMATCH at g(%lld,%lld) "
                            "expected=%f got=%f\n",
                            label, mypnum,
                            (long long)gi, (long long)gj,
                            expected, actual);
                errs++;
            }
        }
    }
    return errs;
}

/* ------------------------------------------------------------------ */
/* collect_errors — MPI_Allreduce sum of error count                  */
/* ------------------------------------------------------------------ */
static int
collect_errors(int local_errors)
{
    int total = 0;
    MPI_Allreduce(&local_errors, &total, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    return total;
}

/* ================================================================== */
/* Test 1: big-block descriptor, single-owner round-trip              */
/* ================================================================== */
static int
test1_big_block(int mypnum, int nprocs)
{
    int gcontext;
    int nprow, npcol, myprow, mypcol;
    int errors = 0;

    int64_t M = 10, N = 5;
    int64_t BIG_MB = INT64_C(3000000000);
    int64_t BIG_NB = INT64_C(3000000000);

    double *A, *B, *C;
    int64_t desc_A[9], desc_B[9], desc_C[9];
    int64_t ia = 1, ja = 1, ib = 1, jb = 1, ic = 1, jc = 1;
    int64_t gc64, lld;

    /* 1 × nprocs grid */
    Cblacs_get(0, 0, &gcontext);
    Cblacs_gridinit(&gcontext, "R", 1, nprocs);
    Cblacs_gridinfo(gcontext, &nprow, &npcol, &myprow, &mypcol);
    if (myprow < 0) return 0;

    /* All processes share row 0, so LOCr(M)=M for everyone */
    lld = M;

    make_desc(desc_A, gcontext, M, N, BIG_MB, BIG_NB, 0, 0, lld);
    memcpy(desc_B, desc_A, sizeof(desc_A));
    memcpy(desc_C, desc_A, sizeof(desc_A));

    A = (double *)calloc((size_t)(M * N), sizeof(double));
    B = (double *)calloc((size_t)(M * N), sizeof(double));
    C = (double *)calloc((size_t)(M * N), sizeof(double));
    assert(A && B && C);

    if (mypcol == 0) {
        for (int j = 0; j < (int)N; j++)
            for (int i = 0; i < (int)M; i++)
                A[j * M + i] = global_value(i, j);
        for (int i = 0; i < (int)(M * N); i++)
            B[i] = C[i] = -999.0;
    }

    gc64 = (int64_t)gcontext;

    /* A → B → C round-trip */
    pdgemr2d_i8(&M, &N, A, &ia, &ja, desc_A,
                         B, &ib, &jb, desc_B, &gc64);
    pdgemr2d_i8(&M, &N, B, &ib, &jb, desc_B,
                         C, &ic, &jc, desc_C, &gc64);

    if (mypcol == 0) {
        for (int j = 0; j < (int)N; j++)
            for (int i = 0; i < (int)M; i++) {
                int idx = j * (int)M + i;
                if (C[idx] != A[idx]) errors++;
                if (B[idx] != A[idx]) errors++;
            }
    }

    free(A); free(B); free(C);
    Cblacs_gridexit(gcontext);
    return collect_errors(errors);
}

/* ================================================================== */
/* Test 2: cross-rank redistribution between different grids          */
/* ================================================================== */
static int
test2_cross_rank(int mypnum, int nprocs)
{
    if (nprocs < 4) {
        if (mypnum == 0)
            printf("  test2: SKIPPED (need >= 4 procs)\n");
        return 0;
    }

    int gcontext, ctx0, ctx1;
    int nprow0, npcol0, myprow0, mypcol0;
    int nprow1, npcol1, myprow1, mypcol1;
    int errors = 0;

    /* Submatrix to redistribute */
    int64_t M = 12, N = 8;

    /* Grid 0: 2×2, block sizes 3×2 */
    int p0 = 2, q0 = 2;
    int64_t mb0 = 3, nb0 = 2;

    /* Grid 1: 1×4, block sizes 4×3 */
    int p1 = 1, q1 = 4;
    int64_t mb1 = 4, nb1 = 3;

    /* Global matrix dimensions (≥ M, N; same for both descriptors) */
    int64_t GM = M, GN = N;

    double *A = NULL, *B = NULL, *C = NULL;
    int64_t desc_A[9], desc_B[9], desc_C[9];
    int64_t ia = 1, ja = 1, ib = 1, jb = 1, ic = 1, jc = 1;
    int64_t gc64;
    int64_t lld0, lld1, locr0, locr1;

    /* Global context for the redistribution call */
    Cblacs_get(0, 0, &gcontext);
    Cblacs_gridinit(&gcontext, "R", 1, nprocs);
    gc64 = (int64_t)gcontext;

    /* Grid 0: 2×2 */
    Cblacs_get(0, 0, &ctx0);
    Cblacs_gridinit(&ctx0, "R", p0, q0);
    Cblacs_gridinfo(ctx0, &nprow0, &npcol0, &myprow0, &mypcol0);
    if (myprow0 >= p0 || mypcol0 >= q0)
        myprow0 = mypcol0 = -1;

    /* Grid 1: 1×4 */
    Cblacs_get(0, 0, &ctx1);
    Cblacs_gridinit(&ctx1, "R", p1, q1);
    Cblacs_gridinfo(ctx1, &nprow1, &npcol1, &myprow1, &mypcol1);
    if (myprow1 >= p1 || mypcol1 >= q1)
        myprow1 = mypcol1 = -1;

    /* ---- Allocate and initialize A on grid 0 ---- */
    if (myprow0 >= 0 && mypcol0 >= 0) {
        locr0 = numroc(GM, mb0, myprow0, 0, p0);
        lld0 = locr0 > 0 ? locr0 : 1;
        int64_t locc0 = numroc(GN, nb0, mypcol0, 0, q0);
        A = (double *)calloc((size_t)(lld0 * (locc0 > 0 ? locc0 : 1)),
                             sizeof(double));
        C = (double *)calloc((size_t)(lld0 * (locc0 > 0 ? locc0 : 1)),
                             sizeof(double));
        assert(A && C);
        init_distributed(A, lld0, GM, GN, mb0, nb0,
                         myprow0, mypcol0, p0, q0, 0, 0);
        /* Fill C with sentinel for later round-trip check */
        for (int64_t k = 0; k < lld0 * (locc0 > 0 ? locc0 : 1); k++)
            C[k] = -999.0;
    } else {
        lld0 = 1;
        A = (double *)calloc(1, sizeof(double));
        C = (double *)calloc(1, sizeof(double));
        assert(A && C);
    }

    /* ---- Allocate B on grid 1 ---- */
    if (myprow1 >= 0 && mypcol1 >= 0) {
        locr1 = numroc(GM, mb1, myprow1, 0, p1);
        lld1 = locr1 > 0 ? locr1 : 1;
        int64_t locc1 = numroc(GN, nb1, mypcol1, 0, q1);
        B = (double *)calloc((size_t)(lld1 * (locc1 > 0 ? locc1 : 1)),
                             sizeof(double));
        assert(B);
        for (int64_t k = 0; k < lld1 * (locc1 > 0 ? locc1 : 1); k++)
            B[k] = -999.0;
    } else {
        lld1 = 1;
        B = (double *)calloc(1, sizeof(double));
        assert(B);
    }

    /* ---- Build descriptors ---- */
    make_desc(desc_A, ctx0, GM, GN, mb0, nb0, 0, 0, lld0);
    make_desc(desc_B, ctx1, GM, GN, mb1, nb1, 0, 0, lld1);
    make_desc(desc_C, ctx0, GM, GN, mb0, nb0, 0, 0, lld0);

    /* ---- Forward: A (grid 0) → B (grid 1) ---- */
    pdgemr2d_i8(&M, &N, A, &ia, &ja, desc_A,
                         B, &ib, &jb, desc_B, &gc64);

    /* Verify B on grid 1 */
    if (myprow1 >= 0 && mypcol1 >= 0) {
        errors += check_distributed(B, lld1, GM, GN, mb1, nb1,
                                    myprow1, mypcol1, p1, q1,
                                    0, 0, "test2 A->B", mypnum);
    }

    /* ---- Inverse: B (grid 1) → C (grid 0) ---- */
    pdgemr2d_i8(&M, &N, B, &ib, &jb, desc_B,
                         C, &ic, &jc, desc_C, &gc64);

    /* Verify C on grid 0 matches A */
    if (myprow0 >= 0 && mypcol0 >= 0) {
        errors += check_distributed(C, lld0, GM, GN, mb0, nb0,
                                    myprow0, mypcol0, p0, q0,
                                    0, 0, "test2 B->C", mypnum);
    }

    free(A); free(B); free(C);
    Cblacs_gridexit(ctx0);
    Cblacs_gridexit(ctx1);
    Cblacs_gridexit(gcontext);
    return collect_errors(errors);
}

/* ================================================================== */
/* main                                                               */
/* ================================================================== */
int
main(int argc, char *argv[])
{
    int mypnum, nprocs;
    int total_errors = 0;

    MPI_Init(&argc, &argv);
    Cblacs_pinfo(&mypnum, &nprocs);

    if (mypnum == 0)
        printf("test_pdgemr2d_i8: nprocs=%d\n", nprocs);

    /* Test 1: big-block descriptor */
    if (mypnum == 0) printf("  test1: big-block descriptor...\n");
    {
        int e = test1_big_block(mypnum, nprocs);
        if (mypnum == 0)
            printf("  test1: %s (%d errors)\n",
                   e == 0 ? "PASSED" : "FAILED", e);
        total_errors += e;
    }

    /* Test 2: cross-rank redistribution */
    if (mypnum == 0) printf("  test2: cross-rank redistribution...\n");
    {
        int e = test2_cross_rank(mypnum, nprocs);
        if (mypnum == 0)
            printf("  test2: %s (%d errors)\n",
                   e == 0 ? "PASSED" : "FAILED", e);
        total_errors += e;
    }

    if (mypnum == 0) {
        if (total_errors == 0)
            printf("TEST PASSED OK\n");
        else
            printf("TEST FAILED: %d total errors\n", total_errors);
    }

    Cblacs_exit(0);
    return (total_errors != 0) ? 1 : 0;
}
