      SUBROUTINE PZDTTRF_I8( N, DL, D, DU, JA, DESCA, AF, LAF,
     $                       WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PZDTTRF.
*
*  Accepts I8 arguments and narrows them for the legacy routine.
*
*     .. Scalar Arguments ..
      INTEGER*8          N, JA, LAF, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      COMPLEX*16         AF( * ), D( * ), DL( * ), DU( * ), WORK( * )
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
      INTEGER            N4, JA4, LAF4, LWORK4
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PZDTTRF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
*     Quick return
*
      INFO = 0
      IF( N.LE.0 ) RETURN
*
*     Range checks
*
      IF( N.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $    LAF.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      END IF
*
*     Narrow scalars
*
      N4   = INT( N )
      JA4  = INT( JA )
      LAF4 = INT( LAF )
*
*     Handle workspace query
*
      IF( LWORK.EQ.-1 ) THEN
         LWORK4 = -1
      ELSE IF( LWORK.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      ELSE
         LWORK4 = INT( LWORK )
      END IF
*
*     Narrow descriptor
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
*
*     Call legacy routine
*
      CALL PZDTTRF( N4, DL, D, DU, JA4, DESCA4, AF, LAF4, WORK,
     $              LWORK4, INFO )
*
      RETURN
*
*     End of PZDTTRF_I8
*
      END
