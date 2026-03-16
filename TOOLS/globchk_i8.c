/*
 * globchk_i8.c — 64-bit global consistency check for PCHK1MAT_I8 / PCHK2MAT_I8
 *
 * This is the INTEGER*8 companion of GLOBCHK (pchkxmat.f).  It checks
 * that an array of int64_t (value, position) pairs are identical across
 * all processes in a BLACS context, and returns the smallest position
 * where a mismatch is found.
 *
 * Uses MPI_Allreduce rather than BLACS integer transport, because
 * BLACS only supports default-INTEGER arrays.
 *
 * Fortran interface:
 *   CALL GLOBCHK_I8( ICTXT, N, XVAL, XPOS, IWORK8, INFO )
 *   INTEGER          ICTXT, N, INFO
 *   INTEGER*8        XVAL(N), IWORK8(N)
 *   INTEGER          XPOS(N)
 */

#include "scalapack-types.h"
#include <mpi.h>
#include <stdint.h>
#include <limits.h>

/* MPI datatype for int64_t */
#ifdef MPI_INT64_T
#define GLOBCHK_MPI_I64 MPI_INT64_T
#else
_Static_assert(sizeof(long long) == sizeof(int64_t),
               "MPI_LONG_LONG fallback requires long long == int64_t");
#define GLOBCHK_MPI_I64 MPI_LONG_LONG
#endif

/* MPI datatype matching ScaLAPACK_ApiInt (Fortran default INTEGER).
 * In an ILP64 build ScaLAPACK_ApiInt is 64-bit, so MPI_INT is wrong. */
#if SCALAPACK_FORTRAN_INT_BYTES == 8
#  define GLOBCHK_MPI_APIINT GLOBCHK_MPI_I64
#elif SCALAPACK_FORTRAN_INT_BYTES == 4
#  define GLOBCHK_MPI_APIINT MPI_INT
#else
#  error "Unsupported SCALAPACK_FORTRAN_INT_BYTES value."
#endif

/* BLACS C-interface declarations */
extern void Cblacs_gridinfo(ScaLAPACK_ApiInt context,
                            ScaLAPACK_ApiInt *nprow, ScaLAPACK_ApiInt *npcol,
                            ScaLAPACK_ApiInt *myrow, ScaLAPACK_ApiInt *mycol);
extern MPI_Comm Cblacs2sys_handle(ScaLAPACK_ApiInt BlacsCtxt);
extern void Cblacs_get(ScaLAPACK_ApiInt context, ScaLAPACK_ApiInt what,
                       ScaLAPACK_ApiInt *val);

/*
 * globchk_i8_ — Fortran-callable entry point (Add_ mangling)
 *
 * Checks that XVAL(1:N) is identical on every process in ICTXT.
 * On mismatch, INFO is set to MIN(INFO, XPOS(k)) for the first
 * disagreeing entry k.
 *
 * Strategy: process (0,0) broadcasts its values; all other processes
 * compare and report the smallest mismatch position via an
 * MPI_Allreduce with MPI_MIN.
 */

/* Fortran name mangling */
#if defined(Add_) || defined(f77IsF2C)
#define globchk_i8_fc globchk_i8_
#elif defined(UpCase)
#define globchk_i8_fc GLOBCHK_I8
#else
#define globchk_i8_fc globchk_i8
#endif

void
globchk_i8_fc(const ScaLAPACK_ApiInt *ictxt,
              const ScaLAPACK_ApiInt *n,
              const int64_t *xval,
              const ScaLAPACK_ApiInt *xpos,
              int64_t *iwork8,
              ScaLAPACK_ApiInt *info)
{
    ScaLAPACK_ApiInt nprow, npcol, myrow, mycol;
    ScaLAPACK_ApiInt sys_handle;
    MPI_Comm comm;
    int k, nn;

    nn = (int)*n;
    if (nn <= 0) return;

    Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);

    /* Get MPI communicator from BLACS context */
    Cblacs_get(*ictxt, (ScaLAPACK_ApiInt)10, &sys_handle);
    comm = Cblacs2sys_handle(sys_handle);

    /* Copy process (0,0)'s values into iwork8 via broadcast */
    for (k = 0; k < nn; k++)
        iwork8[k] = xval[k];

    /* Broadcast from root (rank 0 in the comm corresponds to grid (0,0)) */
    MPI_Bcast(iwork8, nn, GLOBCHK_MPI_I64, 0, comm);

    /* Compare: every process checks against the broadcast values */
    for (k = 0; k < nn; k++) {
        if (xval[k] != iwork8[k]) {
            if (xpos[k] < *info)
                *info = xpos[k];
        }
    }

    /* Global min of INFO across all processes */
    {
        ScaLAPACK_ApiInt local_info = *info;
        ScaLAPACK_ApiInt global_info;
        MPI_Allreduce(&local_info, &global_info, 1, GLOBCHK_MPI_APIINT,
                      MPI_MIN, comm);
        *info = global_info;
    }
}
