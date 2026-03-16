      SUBROUTINE DESC_CONVERT_I8( DESC_IN, DESC_OUT, INFO )
      IMPLICIT NONE
*
*
*     .. Array Arguments ..
      INTEGER*8 DESC_IN( * ), DESC_OUT( * )
      INTEGER   INFO
*     ..
*
*  Purpose
*  =======
*
*  Converts INTEGER*8 descriptors from one type to another if they
*  are compatible.
*
*  Supports *ONLY* an output descriptor type of 1D_horizontal (type
*     number 501) or 1D_vertical (number 502).
*  Supports only one-dimensional 1xP input grids if descriptor_in is 2D.
*
*  This is the INTEGER*8 companion of DESC_CONVERT.  DTYPE values use
*  the same semantic tags as legacy descriptors (1 = 2D block-cyclic,
*  501 = 1D horizontal, 502 = 1D vertical).
*
*  Arguments
*  =========
*
*  DESC_IN: (input) INTEGER*8 array, input descriptor
*
*  DESC_OUT: (input/output) INTEGER*8 array, output descriptor
*            (required to be 1D_horizontal or 1D_vertical).
*            DESC_OUT(1) must be set to the desired output type
*            (501 or 502) before calling.
*
*  INFO: (output) INTEGER return code
*          = 0: successful exit
*          = -1: grid incompatible with requested conversion
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8 DTYPE_1D_H, DTYPE_1D_V
      PARAMETER ( DTYPE_1D_H = 501, DTYPE_1D_V = 502 )
*     ..
*
*     .. Local Scalars ..
      INTEGER*8       DESC_TYPE, DESC_TYPE_IN, ICTXT8
      INTEGER*8       CSRC, RSRC, MB, NB, LLDA, M, N
      INTEGER         NPROW, NPCOL, IDUM1, IDUM2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
      DESC_TYPE_IN = DESC_IN( 1 )
*
*     .. Initialize Variables ..
*
      RSRC = 0
      NB = 0
      N = 0
      MB = 0
      M = 0
      LLDA = 0
      CSRC = 0
      NPROW = 0
      NPCOL = 0
*
      IF( DESC_TYPE_IN .EQ. BLOCK_CYCLIC_2D_I8 ) THEN
         ICTXT8 = DESC_IN( CTXT_ )
         RSRC   = DESC_IN( RSRC_ )
         CSRC   = DESC_IN( CSRC_ )
         MB     = DESC_IN( MB_ )
         NB     = DESC_IN( NB_ )
         LLDA   = DESC_IN( LLD_ )
         M      = DESC_IN( M_ )
         N      = DESC_IN( N_ )
         CALL BLACS_GRIDINFO( INT( ICTXT8 ), NPROW, NPCOL,
     $                        IDUM1, IDUM2 )
      ELSE IF( DESC_TYPE_IN .EQ. DTYPE_1D_V ) THEN
         ICTXT8 = DESC_IN( 2 )
         RSRC   = DESC_IN( 5 )
         CSRC   = 1
         MB     = DESC_IN( 4 )
         NB     = 1
         LLDA   = DESC_IN( 6 )
         M      = DESC_IN( 3 )
         N      = 1
         NPROW  = 0
         NPCOL  = 1
      ELSE IF( DESC_TYPE_IN .EQ. DTYPE_1D_H ) THEN
         ICTXT8 = DESC_IN( 2 )
         RSRC   = 1
         CSRC   = DESC_IN( 5 )
         MB     = 1
         NB     = DESC_IN( 4 )
         LLDA   = DESC_IN( 6 )
         M      = 1
         N      = DESC_IN( 3 )
         NPROW  = 1
         NPCOL  = 0
      END IF
*
*
      DESC_TYPE = DESC_OUT( 1 )
*
      IF( DESC_TYPE .EQ. DTYPE_1D_H ) THEN
         IF( NPROW .NE. 1 ) THEN
            INFO = -1
            RETURN
         END IF
         DESC_OUT( 2 ) = ICTXT8
         DESC_OUT( 5 ) = CSRC
         DESC_OUT( 4 ) = NB
         DESC_OUT( 6 ) = LLDA
         DESC_OUT( 3 ) = N
      ELSE IF( DESC_TYPE .EQ. DTYPE_1D_V ) THEN
         IF( NPCOL .NE. 1 ) THEN
            INFO = -1
            RETURN
         END IF
         DESC_OUT( 2 ) = ICTXT8
         DESC_OUT( 5 ) = RSRC
         DESC_OUT( 4 ) = MB
         DESC_OUT( 6 ) = LLDA
         DESC_OUT( 3 ) = M
      END IF
*
      RETURN
*
*     End of DESC_CONVERT_I8
*
      END
