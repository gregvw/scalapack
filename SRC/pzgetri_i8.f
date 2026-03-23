      SUBROUTINE PZGETRI_I8( N, A, IA, JA, DESCA, IPIV, WORK, LWORK,
     $                       IWORK, LIWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     INTEGER*8 thin wrapper around PZGETRI.
*
*     .. Scalar Arguments ..
      INTEGER*8          N, IA, JA, LWORK, LIWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      INTEGER            IPIV( * ), IWORK( * )
      COMPLEX*16         A( * ), WORK( * )
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
      INTEGER            N4, IA4, JA4, LWORK4, LIWORK4, ICTXT
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PZGETRI, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      INFO = 0
      ICTXT = INT( DESCA( CTXT_ ) )
*
      IF( N.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         INFO = -1
         CALL PXERBLA( ICTXT, 'PZGETRI_I8', -1 )
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
         LWORK4 = INT( INTMAX )
      ELSE
         LWORK4 = INT( LWORK )
      END IF
*
      IF( LIWORK.EQ.-1 ) THEN
         LIWORK4 = -1
      ELSE IF( LIWORK.GT.INTMAX ) THEN
         LIWORK4 = INT( INTMAX )
      ELSE
         LIWORK4 = INT( LIWORK )
      END IF
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
*
      CALL PZGETRI( N4, A, IA4, JA4, DESCA4, IPIV, WORK, LWORK4,
     $              IWORK, LIWORK4, INFO )
*
      RETURN
*
*     End of PZGETRI_I8
*
      END
