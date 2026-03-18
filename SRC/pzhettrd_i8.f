      SUBROUTINE PZHETTRD_I8( UPLO, N, A, IA, JA, DESCA, D, E, TAU,
     $                       WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PZHETTRD.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          N, IA, JA, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
*     ..
*     .. Local Scalars ..
      INTEGER            N4, IA4, JA4, LWORK4
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PZHETTRD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MIN
*     ..
*
      INFO = 0
      IF( N.LE.0 ) RETURN
*
      IF( N.GT.INTMAX .OR. IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      END IF
*
      N4  = INT( N )
      IA4 = INT( IA )
      JA4 = INT( JA )
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
      CALL NARROW_DESC8( DESCA, DESCA4 )
*
      CALL PZHETTRD( UPLO, N4, A, IA4, JA4, DESCA4, D, E, TAU,
     $               WORK, LWORK4, INFO )
*
      RETURN
      END
