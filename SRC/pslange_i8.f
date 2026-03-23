      REAL               FUNCTION PSLANGE_I8( NORM, M, N, A, IA, JA,
     $                                          DESCA, WORK )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PSLANGE.
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER*8          M, N, IA, JA
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      REAL               A( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
*     ..
*     .. Local Scalars ..
      INTEGER            M4, N4, IA4, JA4, ICTXT
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Functions ..
      REAL               PSLANGE
      EXTERNAL           PSLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      ICTXT = INT( DESCA( CTXT_ ) )
*
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         CALL PXERBLA( ICTXT, 'PSLANGE_I8', -1 )
         PSLANGE_I8 = -1.0
         RETURN
      END IF
*
      M4  = INT( M )
      N4  = INT( N )
      IA4 = INT( IA )
      JA4 = INT( JA )
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
*
      PSLANGE_I8 = PSLANGE( NORM, M4, N4, A, IA4, JA4, DESCA4,
     $                       WORK )
*
      RETURN
*
*     End of PSLANGE_I8
*
      END
