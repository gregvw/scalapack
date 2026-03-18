      SUBROUTINE PSLASCL_I8( TYPE, CFROM, CTO, M, N, A, IA, JA,
     $                       DESCA, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PSLASCL.
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER*8          M, N, IA, JA
      INTEGER            INFO
      REAL               CFROM, CTO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      REAL               A( * )
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
      INTEGER            M4, N4, IA4, JA4
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PSLASCL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Range-check INTEGER*8 args
*
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      END IF
*
      M4  = INT( M )
      N4  = INT( N )
      IA4 = INT( IA )
      JA4 = INT( JA )
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
*
      CALL PSLASCL( TYPE, CFROM, CTO, M4, N4, A, IA4, JA4,
     $              DESCA4, INFO )
*
      RETURN
*
*     End of PSLASCL_I8
*
      END
