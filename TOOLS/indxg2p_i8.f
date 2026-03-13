      INTEGER FUNCTION INDXG2P_I8( INDXGLOB, NB, IPROC, ISRCPROC,
     $                             NPROCS )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     .. Scalar Arguments ..
      INTEGER*8          INDXGLOB, NB
      INTEGER            IPROC, ISRCPROC, NPROCS
*     ..
*
*  Purpose
*  =======
*
*  INDXG2P_I8 computes the process coordinate which possesses the entry
*  of a distributed matrix specified by a global index INDXGLOB.
*
*  This is the INTEGER*8 version of INDXG2P, supporting global indices
*  larger than 2^31-1.  The return value is INTEGER since it is a
*  process coordinate (always < NPROCS).
*
*  Arguments
*  =========
*
*  INDXGLOB  (global input) INTEGER*8
*            The global index of the element.
*
*  NB        (global input) INTEGER*8
*            Block size, size of the blocks the distributed matrix is
*            split into.
*
*  IPROC     (local dummy) INTEGER
*            Dummy argument in this case in order to unify the calling
*            sequence of the tool-routines.
*
*  ISRCPROC  (global input) INTEGER
*            The coordinate of the process that possesses the first
*            row/column of the distributed matrix.
*
*  NPROCS    (global input) INTEGER
*            The total number processes over which the matrix is
*            distributed.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER*8          NP8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MOD
*     ..
*     .. Executable Statements ..
*
*     Promote NPROCS to INTEGER*8 for safe mixed-kind MOD
*
      NP8 = NPROCS
*
      INDXG2P_I8 = INT( MOD( ISRCPROC + (INDXGLOB - 1) / NB, NP8 ) )
*
      RETURN
*
*     End of INDXG2P_I8
*
      END
