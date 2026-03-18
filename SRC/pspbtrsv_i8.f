      SUBROUTINE PSPBTRSV_I8( UPLO, TRANS, N, BW, NRHS, A, JA,
     $                        DESCA, B, IB, DESCB, AF, LAF,
     $                        WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 thin wrapper around PSPBTRSV.
*
*  Accepts I8 arguments and narrows them for the legacy routine.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, TRANS
      INTEGER*8          N, BW, NRHS, JA, IB, LAF, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * )
      REAL               A( * ), AF( * ), B( * ), WORK( * )
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
      INTEGER            N4, BW4, NRHS4, JA4, IB4, LAF4, LWORK4
      INTEGER            DESCA4( 9 ), DESCB4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PSPBTRSV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
*     Quick return
*
      INFO = 0
      IF( N.LE.0 ) RETURN
*
*     Range checks
*
      IF( N.GT.INTMAX .OR. BW.GT.INTMAX .OR. NRHS.GT.INTMAX .OR.
     $    JA.GT.INTMAX .OR. IB.GT.INTMAX .OR. LAF.GT.INTMAX ) THEN
         INFO = -1
         RETURN
      END IF
*
*     Narrow scalars
*
      N4    = INT( N )
      BW4   = INT( BW )
      NRHS4 = INT( NRHS )
      JA4   = INT( JA )
      IB4   = INT( IB )
      LAF4  = INT( LAF )
*
*     Handle workspace query
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
*     Narrow descriptors
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
      CALL NARROW_DESC8( DESCB, DESCB4 )
*
*     Call legacy routine
*
      CALL PSPBTRSV( UPLO, TRANS, N4, BW4, NRHS4, A, JA4, DESCA4,
     $               B, IB4, DESCB4, AF, LAF4, WORK, LWORK4, INFO )
*
      RETURN
*
*     End of PSPBTRSV_I8
*
      END
