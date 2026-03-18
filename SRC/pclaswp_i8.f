      SUBROUTINE PCLASWP_I8( DIREC, ROWCOL, N, A, IA, JA, DESCA,
     $                       K1, K2, IPIV )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PCLASWP.
*
*  Accepts I8 arguments and narrows them for the legacy routine.
*  IPIV stays INTEGER (pivot indices, bounded by matrix dimension).
*
*     .. Scalar Arguments ..
      CHARACTER          DIREC, ROWCOL
      INTEGER*8          N, IA, JA, K1, K2
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      INTEGER            IPIV( * )
      COMPLEX            A( * )
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
      INTEGER            N4, IA4, JA4, K14, K24, ICTXT
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PCLASWP, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
*     Range checks
*
      IF( N.GT.INTMAX .OR. IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $    K1.GT.INTMAX .OR. K2.GT.INTMAX ) THEN
         ICTXT = INT( DESCA( CTXT_ ) )
         CALL PXERBLA( ICTXT, 'PCLASWP_I8', 1 )
         RETURN
      END IF
*
*     Narrow scalars
*
      N4  = INT( N )
      IA4 = INT( IA )
      JA4 = INT( JA )
      K14 = INT( K1 )
      K24 = INT( K2 )
*
*     Narrow descriptor
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
*
*     Call legacy routine
*
      CALL PCLASWP( DIREC, ROWCOL, N4, A, IA4, JA4, DESCA4,
     $              K14, K24, IPIV )
*
      RETURN
*
*     End of PCLASWP_I8
*
      END
