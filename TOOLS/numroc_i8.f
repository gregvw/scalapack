      INTEGER*8 FUNCTION NUMROC_I8( N, NB, IPROC, ISRCPROC, NPROCS )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     .. Scalar Arguments ..
      INTEGER*8          N, NB
      INTEGER            IPROC, ISRCPROC, NPROCS
*     ..
*
*  Purpose
*  =======
*
*  NUMROC_I8 computes the NUMber of Rows Or Columns of a distributed
*  matrix owned by the process indicated by IPROC.
*
*  This is the INTEGER*8 version of NUMROC, supporting matrix dimensions
*  larger than 2^31-1.
*
*  Arguments
*  =========
*
*  N         (global input) INTEGER*8
*            The number of rows/columns in distributed matrix.
*
*  NB        (global input) INTEGER*8
*            Block size, size of the blocks the distributed matrix is
*            split into.
*
*  IPROC     (local input) INTEGER
*            The coordinate of the process whose local array row or
*            column is to be determined.
*
*  ISRCPROC  (global input) INTEGER
*            The coordinate of the process that possesses the first
*            row or column of the distributed matrix.
*
*  NPROCS    (global input) INTEGER
*            The total number processes over which the matrix is
*            distributed.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER*8          EXTRABLKS, MYDIST, NBLOCKS, NP8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     Promote NPROCS to INTEGER*8 for safe mixed-kind arithmetic
*
      NP8 = NPROCS
*
*     Figure PROC's distance from source process
*
      MYDIST = MOD( NP8+IPROC-ISRCPROC, NP8 )
*
*     Figure the total number of whole NB blocks N is split up into
*
      NBLOCKS = N / NB
*
*     Figure the minimum number of rows/cols a process can have
*
      NUMROC_I8 = (NBLOCKS/NP8) * NB
*
*     See if there are any extra blocks
*
      EXTRABLKS = MOD( NBLOCKS, NP8 )
*
*     If I have an extra block
*
      IF( MYDIST.LT.EXTRABLKS ) THEN
          NUMROC_I8 = NUMROC_I8 + NB
*
*         If I have last block, it may be a partial block
*
      ELSE IF( MYDIST.EQ.EXTRABLKS ) THEN
          NUMROC_I8 = NUMROC_I8 + MOD( N, NB )
      END IF
*
      RETURN
*
*     End of NUMROC_I8
*
      END
