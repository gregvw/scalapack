      SUBROUTINE PDORMTR_I8( SIDE, UPLO, TRANS, M, N, A, IA, JA,
     $                       DESCA, TAU, C, IC, JC, DESCC, WORK,
     $                       LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PDORMTR.
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER*8          M, N, IA, JA, IC, JC, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCC( * )
      DOUBLE PRECISION   A( * ), C( * ), TAU( * ), WORK( * )
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
      INTEGER            M4, N4, IA4, JA4, IC4, JC4, LWORK4
      INTEGER            DESCA4( 9 ), DESCC4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PDORMTR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      INFO = 0
*
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $    IC.GT.INTMAX .OR. JC.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      END IF
*
      M4  = INT( M )
      N4  = INT( N )
      IA4 = INT( IA )
      JA4 = INT( JA )
      IC4 = INT( IC )
      JC4 = INT( JC )
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
      CALL NARROW_DESC8( DESCC, DESCC4 )
*
      CALL PDORMTR( SIDE, UPLO, TRANS, M4, N4, A, IA4, JA4, DESCA4,
     $              TAU, C, IC4, JC4, DESCC4, WORK, LWORK4, INFO )
*
      RETURN
*
*     End of PDORMTR_I8
*
      END
