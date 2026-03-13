      SUBROUTINE INFOG1L_I8( GINDX, NB, NPROCS, MYROC, ISRCPROC,
     $                       LINDX, ROCSRC )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     .. Scalar Arguments ..
      INTEGER*8          GINDX, NB, LINDX
      INTEGER            NPROCS, MYROC, ISRCPROC, ROCSRC
*     ..
*
*  Purpose
*  =======
*
*  INFOG1L_I8 computes the starting local index LINDX corresponding to
*  the distributed submatrix starting globally at the entry pointed by
*  GINDX.  This routine returns the coordinate of the process in the
*  grid owning the submatrix entry of global index GINDX: ROCSRC.
*  INFOG1L_I8 is a 1-dimensional version of INFOG2L_I8.
*
*  This is the INTEGER*8 version of INFOG1L, supporting global and
*  local indices larger than 2^31-1.
*
*  Arguments
*  =========
*
*  GINDX     (global input) INTEGER*8
*            The global starting index of the submatrix.
*
*  NB        (global input) INTEGER*8
*            The block size.
*
*  NPROCS    (global input) INTEGER
*            The total number of processes over which the distributed
*            submatrix is distributed.
*
*  MYROC     (local input) INTEGER
*            The coordinate of the process calling this routine.
*
*  ISRCPROC  (global input) INTEGER
*            The coordinate of the process having the first entry of
*            the distributed submatrix.
*
*  LINDX     (local output) INTEGER*8
*            The local starting index of the distributed submatrix.
*
*  ROCSRC    (global output) INTEGER
*            The coordinate of the process that possesses the first
*            row and column of the submatrix.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER*8          GCPY, IBLK, NP8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MOD
*     ..
*     .. Executable Statements ..
*
*     Promote NPROCS to INTEGER*8 for safe mixed-kind arithmetic
*
      NP8 = NPROCS
*
      GCPY = GINDX - 1
      IBLK = GCPY / NB
      ROCSRC = INT( MOD( IBLK + ISRCPROC, NP8 ) )
*
      LINDX = ( IBLK / NP8 + 1 ) * NB + 1
*
      IF( MOD( MYROC+NP8-ISRCPROC, NP8 ) .GE.
     $    MOD( IBLK, NP8 ) ) THEN
         IF( MYROC.EQ.ROCSRC )
     $      LINDX = LINDX + MOD( GCPY, NB )
         LINDX = LINDX - NB
      END IF
*
      RETURN
*
*     End of INFOG1L_I8
*
      END
