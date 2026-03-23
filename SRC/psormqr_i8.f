      SUBROUTINE PSORMQR_I8( SIDE, TRANS, M, N, K, A, IA, JA, DESCA,
     $                       TAU, C, IC, JC, DESCC, WORK, LWORK,
     $                       INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     INTEGER*8 thin wrapper around PSORMQR.
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER*8          M, N, K, IA, JA, IC, JC, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCC( * )
      REAL               A( * ), C( * ), TAU( * ), WORK( * )
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
      INTEGER            M4, N4, K4, IA4, JA4, IC4, JC4, LWORK4,
     $                   ICTXT
      INTEGER            DESCA4( 9 ), DESCC4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PSORMQR, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      INFO = 0
      ICTXT = INT( DESCA( CTXT_ ) )
*
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR. K.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $    IC.GT.INTMAX .OR. JC.GT.INTMAX ) THEN
         INFO = -1
         CALL PXERBLA( ICTXT, 'PSORMQR_I8', -1 )
         RETURN
      END IF
*
      M4  = INT( M )
      N4  = INT( N )
      K4  = INT( K )
      IA4 = INT( IA )
      JA4 = INT( JA )
      IC4 = INT( IC )
      JC4 = INT( JC )
*
      IF( LWORK.EQ.-1 ) THEN
         LWORK4 = -1
      ELSE IF( LWORK.GT.INTMAX ) THEN
         LWORK4 = INT( INTMAX )
      ELSE
         LWORK4 = INT( LWORK )
      END IF
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
      CALL NARROW_DESC8( DESCC, DESCC4 )
*
      CALL PSORMQR( SIDE, TRANS, M4, N4, K4, A, IA4, JA4, DESCA4,
     $              TAU, C, IC4, JC4, DESCC4, WORK, LWORK4, INFO )
*
      RETURN
*
*     End of PSORMQR_I8
*
      END
