      SUBROUTINE DESCSET_I8( DESC, M, N, MB, NB, IRSRC, ICSRC,
     $                       ICTXT, LLD )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     .. Scalar Arguments ..
      INTEGER*8          M, N, MB, NB, LLD
      INTEGER            IRSRC, ICSRC, ICTXT
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESC( * )
*     ..
*
*  Purpose
*  =======
*
*  DESCSET_I8 initializes an INTEGER*8 descriptor vector with the
*  8 input arguments M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD.
*
*  This is the INTEGER*8 version of DESCSET.  No validation is
*  performed.  The descriptor type is set to BLOCK_CYCLIC_2D_I8 = 501.
*
*  Arguments
*  =========
*
*  DESC    (output) INTEGER*8 array of dimension DLEN_.
*          The array descriptor of a distributed matrix to be set.
*
*  M       (global input) INTEGER*8
*          The number of rows in the distributed matrix. M >= 0.
*
*  N       (global input) INTEGER*8
*          The number of columns in the distributed matrix. N >= 0.
*
*  MB      (global input) INTEGER*8
*          The blocking factor used to distribute the rows of the
*          matrix. MB >= 1.
*
*  NB      (global input) INTEGER*8
*          The blocking factor used to distribute the columns of the
*          matrix. NB >= 1.
*
*  IRSRC   (global input) INTEGER
*          The process row over which the first row of the matrix is
*          distributed. 0 <= IRSRC < NPROW.
*
*  ICSRC   (global input) INTEGER
*          The process column over which the first column of the
*          matrix is distributed. 0 <= ICSRC < NPCOL.
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation on the matrix. The context itself is global.
*
*  LLD     (local input) INTEGER*8
*          The leading dimension of the local array storing the local
*          blocks of the distributed matrix. LLD >= MAX(1,LOCr(M)).
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
*     ..
*
*     .. Executable Statements ..
*
      DESC( DTYPE_ ) = BLOCK_CYCLIC_2D_I8
      DESC( M_ ) = M
      DESC( N_ ) = N
      DESC( MB_ ) = MB
      DESC( NB_ ) = NB
      DESC( RSRC_ ) = IRSRC
      DESC( CSRC_ ) = ICSRC
      DESC( CTXT_ ) = ICTXT
      DESC( LLD_ ) = LLD
*
      RETURN
*
*     End DESCSET_I8
*
      END
