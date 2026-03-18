      SUBROUTINE PZDBTRF_I8( N, BWL, BWU, A, JA, DESCA, AF, LAF,
     $                       WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     INTEGER*8 thin wrapper around PZDBTRF.
*
*  Accepts I8 arguments and narrows them for the legacy routine.
*
*     .. Scalar Arguments ..
      INTEGER*8          N, BWL, BWU, JA, LAF, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      COMPLEX*16         A( * ), AF( * ), WORK( * )
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
      INTEGER            N4, BWL4, BWU4, JA4, LAF4, LWORK4
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PZDBTRF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
*
      INFO = 0
*
*     Range checks
*
      IF( N.GT.INTMAX .OR. BWL.GT.INTMAX .OR. BWU.GT.INTMAX .OR.
     $    JA.GT.INTMAX .OR. LAF.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      END IF
*
*     Narrow scalars
*
      N4   = INT( N )
      BWL4 = INT( BWL )
      BWU4 = INT( BWU )
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
      CALL PZDBTRF( N4, BWL4, BWU4, A, JA4, DESCA4, AF, LAF4,
     $              WORK, LWORK4, INFO )
*
      RETURN
*
*     End of PZDBTRF_I8
*
      END
