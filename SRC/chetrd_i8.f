      SUBROUTINE CHETRD_I8( UPLO, N, A, LDA, D, E, TAU, WORK,
     $                     LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around CHETRD (LAPACK).
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          N, LDA, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * )
      COMPLEX            A( * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
*     ..
*     .. Local Scalars ..
      INTEGER            N4, LDA4, LWORK4
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHETRD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MIN
*     ..
*
      INFO = 0
*
      IF( N.GT.INTMAX .OR. LDA.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      END IF
*
      N4   = INT( N )
      LDA4 = INT( LDA )
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
      CALL CHETRD( UPLO, N4, A, LDA4, D, E, TAU, WORK,
     $             LWORK4, INFO )
*
      RETURN
      END
