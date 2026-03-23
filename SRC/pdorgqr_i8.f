      SUBROUTINE PDORGQR_I8( M, N, K, A, IA, JA, DESCA, TAU, WORK,
     $                       LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     INTEGER*8 thin wrapper around PDORGQR.
*
*  Generates an M-by-N real distributed matrix Q denoting A(IA:IA+M-1,
*  JA:JA+N-1) with orthonormal columns, which is defined as the first
*  N columns of a product of K elementary reflectors of order M.
*
*     .. Scalar Arguments ..
      INTEGER*8          M, N, K, IA, JA, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      DOUBLE PRECISION   A( * ), TAU( * ), WORK( * )
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
      INTEGER            M4, N4, K4, IA4, JA4, LWORK4, ICTXT
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PDORGQR, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      INFO = 0
      ICTXT = INT( DESCA( CTXT_ ) )
*
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR. K.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         INFO = -1
         CALL PXERBLA( ICTXT, 'PDORGQR_I8', -1 )
         RETURN
      END IF
*
      M4  = INT( M )
      N4  = INT( N )
      K4  = INT( K )
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
      CALL NARROW_DESC8( DESCA, DESCA4 )
*
      CALL PDORGQR( M4, N4, K4, A, IA4, JA4, DESCA4, TAU, WORK,
     $              LWORK4, INFO )
*
      RETURN
*
*     End of PDORGQR_I8
*
      END
