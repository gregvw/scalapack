      INTEGER*8 FUNCTION INDXL2G_I8( INDXLOC, NB, IPROC, ISRCPROC,
     $                               NPROCS )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     .. Scalar Arguments ..
      INTEGER*8          INDXLOC, NB
      INTEGER            IPROC, ISRCPROC, NPROCS
*     ..
*
*  Purpose
*  =======
*
*  INDXL2G_I8 computes the global index of a distributed matrix entry
*  pointed to by the local index INDXLOC of the process indicated by
*  IPROC.
*
*  This is the INTEGER*8 version of INDXL2G, supporting indices
*  larger than 2^31-1.  All arithmetic is done in INTEGER*8 to avoid
*  overflow in the NPROCS*NB*((INDXLOC-1)/NB) product.
*
*  Arguments
*  =========
*
*  INDXLOC   (global input) INTEGER*8
*            The local index of the distributed matrix entry.
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
*            row/column of the distributed matrix.
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
*     Promote NPROCS to INTEGER*8 to avoid overflow in NP8*NB product
*
      NP8 = NPROCS
*
      INDXL2G_I8 = NP8*NB*((INDXLOC-1)/NB) + MOD(INDXLOC-1,NB) +
     $             MOD(NP8+IPROC-ISRCPROC, NP8)*NB + 1
*
      RETURN
*
*     End of INDXL2G_I8
*
      END
