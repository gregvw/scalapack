      INTEGER*8 FUNCTION INDXG2L_I8( INDXGLOB, NB, IPROC, ISRCPROC,
     $                               NPROCS )
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
*  INDXG2L_I8 computes the local index of a distributed matrix entry
*  pointed to by the global index INDXGLOB.
*
*  This is the INTEGER*8 version of INDXG2L, supporting global indices
*  larger than 2^31-1.  The NB*NPROCS product is computed in INTEGER*8
*  to avoid overflow.
*
*  Arguments
*  =========
*
*  INDXGLOB  (global input) INTEGER*8
*            The global index of the distributed matrix entry.
*
*  NB        (global input) INTEGER*8
*            Block size, size of the blocks the distributed matrix is
*            split into.
*
*  IPROC     (local dummy) INTEGER
*            Dummy argument in this case in order to unify the calling
*            sequence of the tool-routines.
*
*  ISRCPROC  (local dummy) INTEGER
*            Dummy argument in this case in order to unify the calling
*            sequence of the tool-routines.
*
*  NPROCS    (global input) INTEGER
*            The total number processes over which the distributed
*            matrix is distributed.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER*8          NP8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     Promote NPROCS to INTEGER*8 to avoid overflow in NB*NP8 product
*
      NP8 = NPROCS
*
      INDXG2L_I8 = NB*((INDXGLOB-1)/(NB*NP8))+MOD(INDXGLOB-1,NB)+1
*
      RETURN
*
*     End of INDXG2L_I8
*
      END
