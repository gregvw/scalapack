      SUBROUTINE PSLAMVE_I8( UPLO, M, N, A, IA, JA, DESCA,
     $                       B, IB, JB, DESCB, DWORK )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 version of PSLAMVE.  See PDLAMVE_I8 for full docs.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          M, N, IA, JA, IB, JB
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * )
      REAL               A( * ), B( * ), DWORK( * )
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
      LOGICAL            UPPER, LOWER, FULL
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, NPROCS
      INTEGER            M4, N4, IA4, JA4, IB4, JB4
      INTEGER            DESCA4( 9 ), DESCB4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, BLACS_ABORT, PXERBLA,
     $                   SLAMOV, PSGEMR2D_I8, PSLACPY, NARROW_DESC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*     .. Executable Statements ..
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      ICTXT = INT( DESCA( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT. UPPER ) LOWER = LSAME( UPLO, 'L' )
      FULL = (.NOT. UPPER) .AND. (.NOT. LOWER)
*
      NPROCS = NPROW * NPCOL
*
      IF( NPROCS.EQ.1 .OR. .NOT. FULL ) THEN
*
         IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $       IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $       IB.GT.INTMAX .OR. JB.GT.INTMAX .OR.
     $       DESCA( LLD_ ).GT.INTMAX .OR.
     $       DESCB( LLD_ ).GT.INTMAX ) THEN
            CALL PXERBLA( ICTXT, 'PSLAMVE_I8', -2 )
            CALL BLACS_ABORT( ICTXT, 1 )
         END IF
*
         M4  = INT( M )
         N4  = INT( N )
         IA4 = INT( IA )
         JA4 = INT( JA )
         IB4 = INT( IB )
         JB4 = INT( JB )
         CALL NARROW_DESC( DESCA, DESCA4 )
         CALL NARROW_DESC( DESCB, DESCB4 )
*
         IF( NPROCS.EQ.1 ) THEN
            CALL SLAMOV( UPLO, M4, N4,
     $           A((JA4-1)*DESCA4(9)+IA4), DESCA4(9),
     $           B((JB4-1)*DESCB4(9)+IB4), DESCB4(9) )
         ELSE
            CALL PSGEMR2D_I8( M, N, A, IA, JA, DESCA,
     $           DWORK, IB, JB, DESCB, INT( ICTXT, 8 ) )
            CALL PSLACPY( UPLO, M4, N4, DWORK, IB4, JB4, DESCB4,
     $           B, IB4, JB4, DESCB4 )
         END IF
*
      ELSE
*
         CALL PSGEMR2D_I8( M, N, A, IA, JA, DESCA,
     $        B, IB, JB, DESCB, INT( ICTXT, 8 ) )
*
      END IF
*
      RETURN
*
*     End of PSLAMVE_I8
*
      END
