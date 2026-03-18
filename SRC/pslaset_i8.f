      SUBROUTINE PSLASET_I8( UPLO, M, N, ALPHA, BETA, A, IA, JA,
     $                       DESCA )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PSLASET.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          M, N, IA, JA
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      REAL               A( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLASET_I8 is an INTEGER*8 wrapper around PSLASET.  It accepts
*  64-bit integer arguments, range-checks them, narrows the
*  descriptor, and delegates to the legacy 32-bit routine.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
*     ..
*     .. Local Scalars ..
      INTEGER            M4, N4, IA4, JA4, ICTXT
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PSLASET, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
*     Range-check I8 scalar arguments
*
      ICTXT = INT( DESCA( CTXT_ ) )
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         CALL PXERBLA( ICTXT, 'PSLASET_I8', -1 )
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
      CALL PSLASET( UPLO, M4, N4, ALPHA, BETA, A, IA4, JA4, DESCA4 )
*
      RETURN
*
*     End of PSLASET_I8
*
      END
