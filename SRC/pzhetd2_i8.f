      SUBROUTINE PZHETD2_I8( UPLO, N, A, IA, JA, DESCA, D, E, TAU,
     $                       WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PZHETD2.
*
*  Accepts I8 arguments and narrows them for the legacy routine.
*  N is always bounded by NB at the call site, so no overflow is
*  possible in practice; the range checks are purely defensive.
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
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
*     ..
*     .. Local Scalars ..
      INTEGER            N4, IA4, JA4, LWORK4
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PZHETD2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MIN
*     ..
*
*     Quick return
*
      INFO = 0
      IF( N.LE.0 ) RETURN
*
*     Range checks (defensive — N <= NB at call site)
*
      IF( N.GT.INTMAX .OR. IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      END IF
*
*     Narrow scalars
*
      N4  = INT( N )
      IA4 = INT( IA )
      JA4 = INT( JA )
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
      CALL PZHETD2( UPLO, N4, A, IA4, JA4, DESCA4, D, E, TAU,
     $              WORK, LWORK4, INFO )
*
      RETURN
*
*     End of PZHETD2_I8
*
      END
