      SUBROUTINE PDGEQRF_I8( M, N, A, IA, JA, DESCA, TAU, WORK,
     $                       LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     INTEGER*8 thin wrapper around PDGEQRF.
*
*  PDGEQRF_I8 computes a QR factorization of a real distributed
*  M-by-N matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) = Q * R.
*
*  All public integer arguments except INFO are INTEGER*8.
*
*     .. Scalar Arguments ..
      INTEGER*8          M, N, IA, JA, LWORK
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
      INTEGER            M4, N4, IA4, JA4, LWORK4, ICTXT
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PDGEQRF, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      INFO = 0
      ICTXT = INT( DESCA( CTXT_ ) )
*
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         INFO = -1
         CALL PXERBLA( ICTXT, 'PDGEQRF_I8', -1 )
         RETURN
      END IF
*
      M4  = INT( M )
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
      CALL NARROW_DESC8( DESCA, DESCA4 )
*
      CALL PDGEQRF( M4, N4, A, IA4, JA4, DESCA4, TAU, WORK,
     $              LWORK4, INFO )
*
      RETURN
*
*     End of PDGEQRF_I8
*
      END
