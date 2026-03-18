      SUBROUTINE PSLAPIV_I8( DIREC, ROWCOL, PIVROC, M, N, A, IA, JA,
     $                       DESCA, IPIV, IP, JP, DESCIP, IWORK )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PSLAPIV.
*
*  Accepts I8 arguments and narrows them for the legacy routine.
*  IPIV and IWORK stay INTEGER (pivot indices, bounded).
*  Two descriptors DESCA, DESCIP are narrowed.
*
*     .. Scalar Arguments ..
      CHARACTER*1        DIREC, PIVROC, ROWCOL
      INTEGER*8          M, N, IA, JA, IP, JP
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCIP( * )
      INTEGER            IPIV( * ), IWORK( * )
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
      INTEGER            M4, N4, IA4, JA4, IP4, JP4, ICTXT
      INTEGER            DESCA4( 9 ), DESCIP4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PSLAPIV, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
*     Range checks
*
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $    IP.GT.INTMAX .OR. JP.GT.INTMAX ) THEN
         ICTXT = INT( DESCA( CTXT_ ) )
         CALL PXERBLA( ICTXT, 'PSLAPIV_I8', 1 )
         RETURN
      END IF
*
*     Narrow scalars
*
      M4  = INT( M )
      N4  = INT( N )
      IA4 = INT( IA )
      JA4 = INT( JA )
      IP4 = INT( IP )
      JP4 = INT( JP )
*
*     Narrow descriptors
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
      CALL NARROW_DESC8( DESCIP, DESCIP4 )
*
*     Call legacy routine
*
      CALL PSLAPIV( DIREC, ROWCOL, PIVROC, M4, N4, A, IA4, JA4,
     $              DESCA4, IPIV, IP4, JP4, DESCIP4, IWORK )
*
      RETURN
*
*     End of PSLAPIV_I8
*
      END
