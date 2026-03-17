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
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER, LOWER, FULL
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, NPROCS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, SLAMOV_I8,
     $                   PSGEMR2D_I8, PSLACPY_I8
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
      IF( NPROCS.EQ.1 ) THEN
         CALL SLAMOV_I8( UPLO, M, N,
     $        A( (JA-1)*DESCA(LLD_)+IA ), DESCA( LLD_ ),
     $        B( (JB-1)*DESCB(LLD_)+IB ), DESCB( LLD_ ) )
      ELSE IF( FULL ) THEN
         CALL PSGEMR2D_I8( M, N, A, IA, JA, DESCA,
     $        B, IB, JB, DESCB, INT( ICTXT, 8 ) )
      ELSE
         CALL PSGEMR2D_I8( M, N, A, IA, JA, DESCA,
     $        DWORK, IB, JB, DESCB, INT( ICTXT, 8 ) )
         CALL PSLACPY_I8( UPLO, M, N, DWORK, IB, JB, DESCB,
     $        B, IB, JB, DESCB )
      END IF
*
      RETURN
*
*     End of PSLAMVE_I8
*
      END
