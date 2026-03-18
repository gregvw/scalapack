      SUBROUTINE PZHEEV_I8( JOBZ, UPLO, N, A, IA, JA, DESCA, W,
     $                      Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK,
     $                      LRWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PZHEEV.
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER*8          N, IA, JA, IZ, JZ, LWORK, LRWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCZ( * )
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( * ), WORK( * ), Z( * )
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
      INTEGER            N4, IA4, JA4, IZ4, JZ4, LWORK4, LRWORK4
      INTEGER            DESCA4( 9 ), DESCZ4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PZHEEV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      INFO = 0
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
      IF( LRWORK.EQ.-1 ) THEN
         LRWORK4 = -1
      ELSE IF( LRWORK.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      ELSE
         LRWORK4 = INT( LRWORK )
      END IF
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
      CALL NARROW_DESC8( DESCZ, DESCZ4 )
*
      CALL PZHEEV( JOBZ, UPLO, N4, A, IA4, JA4, DESCA4, W,
     $             Z, IZ4, JZ4, DESCZ4, WORK, LWORK4, RWORK,
     $             LRWORK4, INFO )
*
      RETURN
*
*     End of PZHEEV_I8
*
      END
