/*
 * test_pdgemr2d_i8.c — Caller-side test for the pdgemr2d_i8_ entry point
 *
 * Exercises:
 *   - int64_t descriptor fields (MB_, NB_ > INT_MAX)
 *   - unpack_desc_i8 conversion
 *   - redist_sync_params_i8 (param64 MPI_Allreduce)
 *   - 64-bit arithmetic in the redistribution core
 *
 * Strategy:
 *   Global matrix is M=10, N=5 with MB=NB=3000000000 (> INT_MAX).
 *   Since MB > M and NB > N, all data lives on process (0,0) in each
 *   grid — local footprint is just 50 doubles.  We redistribute A→B
 *   on the same 1×nprocs grid, then B→C back, and verify C == A.
 *
 * Usage:  mpirun -np 4 ./xdgemr_i8
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
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

/* Descriptor field indices (matches ScaLAPACK convention) */
#define DTYPE_ 0
#define CTXT_  1
#define M_     2
#define N_     3
#define MB_    4
#define NB_    5
#define RSRC_  6
#define CSRC_  7
#define LLD_   8

int
main(int argc, char *argv[])
{
    int mypnum, nprocs;
    int gcontext;
    int nprow, npcol, myprow, mypcol;
    int errors = 0;

    /* Matrix dimensions (small) */
    int64_t M = 10, N = 5;

    /* Block sizes > INT_MAX to exercise I8 arithmetic */
    int64_t BIG_MB = INT64_C(3000000000);  /* ~3 billion */
    int64_t BIG_NB = INT64_C(3000000000);

    double *A = NULL, *B = NULL, *C = NULL;
    int64_t desc_A[9], desc_B[9], desc_C[9];
    int64_t ia, ja, ib, jb, ic, jc;
    int64_t gc64;
    int64_t lld;
    int i;

    MPI_Init(&argc, &argv);

    Cblacs_pinfo(&mypnum, &nprocs);

    /* Create a 1 × nprocs grid */
    Cblacs_get(0, 0, &gcontext);
    Cblacs_gridinit(&gcontext, "R", 1, nprocs);
    Cblacs_gridinfo(gcontext, &nprow, &npcol, &myprow, &mypcol);

    if (myprow < 0 || mypcol < 0) {
        /* Not part of grid — just finalize */
        Cblacs_exit(0);
        MPI_Finalize();
        return 0;
    }

    assert(nprow == 1 && npcol == nprocs);

    /* In a 1×nprocs grid all processes share process row 0, so
     * LOCr(M) = M = 10 for every process.  LLD must be >= LOCr(M)
     * on all processes, even those with 0 local columns. */
    lld = M;

    /* Build I8 descriptors for A, B, C — all on the same grid */
    desc_A[DTYPE_] = BLOCK_CYCLIC_2D_I8;
    desc_A[CTXT_]  = (int64_t)gcontext;
    desc_A[M_]     = M;
    desc_A[N_]     = N;
    desc_A[MB_]    = BIG_MB;
    desc_A[NB_]    = BIG_NB;
    desc_A[RSRC_]  = 0;
    desc_A[CSRC_]  = 0;
    desc_A[LLD_]   = lld;

    memcpy(desc_B, desc_A, sizeof(desc_A));
    memcpy(desc_C, desc_A, sizeof(desc_A));

    /* Allocate local storage — only (0,0) needs M*N doubles */
    if (myprow == 0 && mypcol == 0) {
        A = (double *)calloc((size_t)(M * N), sizeof(double));
        B = (double *)calloc((size_t)(M * N), sizeof(double));
        C = (double *)calloc((size_t)(M * N), sizeof(double));
        assert(A && B && C);

        /* Initialize A with known pattern: A[i,j] = i*1000 + j */
        for (int j = 0; j < (int)N; j++)
            for (int i = 0; i < (int)M; i++)
                A[j * M + i] = (double)(i * 1000 + j);

        /* Initialize B and C to sentinel */
        for (i = 0; i < (int)(M * N); i++) {
            B[i] = -999.0;
            C[i] = -999.0;
        }
    } else {
        /* Non-owner processes still need lld-sized buffers since
         * lld = M and the redistribution core may reference them. */
        A = (double *)calloc((size_t)M, sizeof(double));
        B = (double *)calloc((size_t)M, sizeof(double));
        C = (double *)calloc((size_t)M, sizeof(double));
        assert(A && B && C);
    }

    /* Redistribution parameters (1-based Fortran indices) */
    ia = 1; ja = 1;
    ib = 1; jb = 1;
    ic = 1; jc = 1;
    gc64 = (int64_t)gcontext;

    if (mypnum == 0)
        printf("test_pdgemr2d_i8: M=%lld N=%lld MB=%lld NB=%lld nprocs=%d\n",
               (long long)M, (long long)N,
               (long long)BIG_MB, (long long)BIG_NB, nprocs);

    /* Test 1: Redistribute A → B (same grid, same layout) */
    pdgemr2d_i8(&M, &N,
                A, &ia, &ja, desc_A,
                B, &ib, &jb, desc_B,
                &gc64);

    /* Test 2: Redistribute B → C (same grid, same layout) */
    pdgemr2d_i8(&M, &N,
                B, &ib, &jb, desc_B,
                C, &ic, &jc, desc_C,
                &gc64);

    /* Verify on process (0,0): C must equal A */
    if (myprow == 0 && mypcol == 0) {
        for (int j = 0; j < (int)N; j++) {
            for (int i = 0; i < (int)M; i++) {
                int idx = j * (int)M + i;
                if (C[idx] != A[idx]) {
                    if (errors < 10)
                        fprintf(stderr,
                                "MISMATCH at (%d,%d): A=%f C=%f\n",
                                i, j, A[idx], C[idx]);
                    errors++;
                }
            }
        }

        /* Also verify B == A (first redistribution) */
        for (int j = 0; j < (int)N; j++) {
            for (int i = 0; i < (int)M; i++) {
                int idx = j * (int)M + i;
                if (B[idx] != A[idx]) {
                    if (errors < 10)
                        fprintf(stderr,
                                "MISMATCH (A→B) at (%d,%d): A=%f B=%f\n",
                                i, j, A[idx], B[idx]);
                    errors++;
                }
            }
        }
    }

    free(A);
    free(B);
    free(C);

    /* Collect errors across all processes */
    {
        int total_errors = 0;
        MPI_Allreduce(&errors, &total_errors, 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD);
        errors = total_errors;
    }

    if (mypnum == 0) {
        if (errors == 0)
            printf("TEST PASSED OK\n");
        else
            printf("TEST FAILED: %d errors\n", errors);
    }

    Cblacs_gridexit(gcontext);
    Cblacs_exit(0);
    /* MPI_Finalize called by Cblacs_exit */

    return (errors != 0) ? 1 : 0;
}
