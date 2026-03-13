      SUBROUTINE DESCINIT_I8( DESC, M, N, MB, NB, IRSRC, ICSRC,
     $                        ICTXT, LLD, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     .. Scalar Arguments ..
      INTEGER*8          M, N, MB, NB, LLD
      INTEGER            IRSRC, ICSRC, ICTXT, INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESC( * )
*     ..
*
*  Purpose
*  =======
*
*  DESCINIT_I8 initializes an INTEGER*8 descriptor vector with the
*  8 input arguments M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD.
*
*  This is the INTEGER*8 version of DESCINIT, supporting matrix
*  dimensions larger than 2^31-1.  The descriptor type is set to
*  BLOCK_CYCLIC_2D_I8 = 501 to distinguish from legacy descriptors.
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
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Note
*  ====
*
*  If the routine can recover from an erroneous input argument, it will
*  return an acceptable descriptor vector.  For example, if LLD = 0 on
*  input, DESC(LLD_) will contain the smallest leading dimension
*  required to store the specified M-by-N distributed matrix, INFO
*  will be set  -9 in that case.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
*     ..
*     .. Local Scalars ..
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( MB.LT.1 ) THEN
         INFO = -4
      ELSE IF( NB.LT.1 ) THEN
         INFO = -5
      ELSE IF( IRSRC.LT.0 .OR. IRSRC.GE.NPROW ) THEN
         INFO = -6
      ELSE IF( ICSRC.LT.0 .OR. ICSRC.GE.NPCOL ) THEN
         INFO = -7
      ELSE IF( NPROW.EQ.-1 ) THEN
         INFO = -8
      ELSE IF( LLD.LT.MAX( 1_8, NUMROC_I8( M, MB, MYROW, IRSRC,
     $                                NPROW ) ) ) THEN
         INFO = -9
      END IF
*
      IF( INFO.NE.0 )
     $   CALL PXERBLA( ICTXT, 'DESCINIT_I8', -INFO )
*
      DESC( DTYPE_ ) = BLOCK_CYCLIC_2D_I8
      DESC( M_ )  = MAX( 0_8, M )
      DESC( N_ )  = MAX( 0_8, N )
      DESC( MB_ ) = MAX( 1_8, MB )
      DESC( NB_ ) = MAX( 1_8, NB )
      DESC( RSRC_ ) = MAX( 0, MIN( IRSRC, NPROW-1 ) )
      DESC( CSRC_ ) = MAX( 0, MIN( ICSRC, NPCOL-1 ) )
      DESC( CTXT_ ) = ICTXT
      DESC( LLD_ )  = MAX( LLD, MAX( 1_8, NUMROC_I8( DESC( M_ ),
     $                     DESC( MB_ ), MYROW, INT( DESC( RSRC_ ) ),
     $                     NPROW ) ) )
*
      RETURN
*
*     End DESCINIT_I8
*
      END
