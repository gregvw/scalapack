      SUBROUTINE PSSYEV_I8( JOBZ, UPLO, N, A, IA, JA, DESCA, W,
     $                      Z, IZ, JZ, DESCZ, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PSSYEV.
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER*8          N, IA, JA, IZ, JZ, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCZ( * )
      REAL               A( * ), W( * ), WORK( * ), Z( * )
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
      INTEGER            N4, IA4, JA4, IZ4, JZ4, LWORK4
      INTEGER            DESCA4( 9 ), DESCZ4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PSSYEV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      INFO = 0
      IF( N.LE.0 ) RETURN
*
      IF( N.GT.INTMAX .OR. IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $    IZ.GT.INTMAX .OR. JZ.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      END IF
*
      N4  = INT( N )
      IA4 = INT( IA )
      JA4 = INT( JA )
      IZ4 = INT( IZ )
      JZ4 = INT( JZ )
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
      CALL NARROW_DESC8( DESCZ, DESCZ4 )
*
      CALL PSSYEV( JOBZ, UPLO, N4, A, IA4, JA4, DESCA4, W,
     $             Z, IZ4, JZ4, DESCZ4, WORK, LWORK4, INFO )
*
      RETURN
*
*     End of PSSYEV_I8
*
      END
