      SUBROUTINE PCGESVD_I8( JOBU, JOBVT, M, N, A, IA, JA, DESCA,
     $                       S, U, IU, JU, DESCU, VT, IVT, JVT,
     $                       DESCVT, WORK, LWORK, RWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     INTEGER*8 thin wrapper around PCGESVD.
*
*     .. Scalar Arguments ..
      CHARACTER          JOBU, JOBVT
      INTEGER*8          M, N, IA, JA, IU, JU, IVT, JVT, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCU( * ), DESCVT( * )
      COMPLEX            A( * ), U( * ), VT( * ), WORK( * )
      REAL               S( * ), RWORK( * )
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
      INTEGER            M4, N4, IA4, JA4, IU4, JU4, IVT4, JVT4,
     $                   LWORK4, ICTXT
      INTEGER            DESCA4( 9 ), DESCU4( 9 ), DESCVT4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           NARROW_DESC8, PCGESVD, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*
      INFO = 0
      ICTXT = INT( DESCA( CTXT_ ) )
*
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $    IU.GT.INTMAX .OR. JU.GT.INTMAX .OR.
     $    IVT.GT.INTMAX .OR. JVT.GT.INTMAX ) THEN
         INFO = -1
         CALL PXERBLA( ICTXT, 'PCGESVD_I8', -1 )
         RETURN
      END IF
*
      M4   = INT( M )
      N4   = INT( N )
      IA4  = INT( IA )
      JA4  = INT( JA )
      IU4  = INT( IU )
      JU4  = INT( JU )
      IVT4 = INT( IVT )
      JVT4 = INT( JVT )
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
      CALL NARROW_DESC8( DESCU, DESCU4 )
      CALL NARROW_DESC8( DESCVT, DESCVT4 )
*
      CALL PCGESVD( JOBU, JOBVT, M4, N4, A, IA4, JA4, DESCA4,
     $              S, U, IU4, JU4, DESCU4, VT, IVT4, JVT4,
     $              DESCVT4, WORK, LWORK4, RWORK, INFO )
*
      RETURN
*
*     End of PCGESVD_I8
*
      END
