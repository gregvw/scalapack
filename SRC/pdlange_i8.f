      DOUBLE PRECISION   FUNCTION PDLANGE_I8( NORM, M, N, A, IA, JA,
     $                                          DESCA, WORK )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PDLANGE.
*
*  Returns the value of the one norm, Frobenius norm, infinity norm,
*  or element of largest absolute value of a distributed matrix.
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER*8          M, N, IA, JA
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      DOUBLE PRECISION   A( * ), WORK( * )
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
      INTEGER            M4, N4, IA4, JA4, ICTXT
      INTEGER            DESCA4( 9 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   PDLANGE
      EXTERNAL           PDLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      ICTXT = INT( DESCA( CTXT_ ) )
*
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         CALL PXERBLA( ICTXT, 'PDLANGE_I8', -1 )
         PDLANGE_I8 = -1.0D0
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
      PDLANGE_I8 = PDLANGE( NORM, M4, N4, A, IA4, JA4, DESCA4,
     $                       WORK )
*
      RETURN
*
*     End of PDLANGE_I8
*
      END
