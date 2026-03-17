      SUBROUTINE PDLAMVE_I8( UPLO, M, N, A, IA, JA, DESCA,
     $                       B, IB, JB, DESCB, DWORK )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 version of PDLAMVE.
*
*  PDLAMVE_I8 copies all or part of a distributed matrix A to another
*  distributed matrix B.  There are no alignment assumptions except
*  that A and B are of the same size.
*
*  All paths are fully I8-native:
*  - Single process: DLAMOV_I8
*  - Full copy, multi-process: PDGEMR2D_I8
*  - Triangular, multi-process: PDGEMR2D_I8 + PDLACPY_I8
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          M, N, IA, JA, IB, JB
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * ), DWORK( * )
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
      EXTERNAL           BLACS_GRIDINFO, DLAMOV_I8,
     $                   PDGEMR2D_I8, PDLACPY_I8
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
*
*        Single process: local copy via DLAMOV_I8
*
         CALL DLAMOV_I8( UPLO, M, N,
     $        A( (JA-1)*DESCA(LLD_)+IA ), DESCA( LLD_ ),
     $        B( (JB-1)*DESCB(LLD_)+IB ), DESCB( LLD_ ) )
*
      ELSE IF( FULL ) THEN
*
*        Full copy, multi-process: native I8 redistribution
*
         CALL PDGEMR2D_I8( M, N, A, IA, JA, DESCA,
     $        B, IB, JB, DESCB, INT( ICTXT, 8 ) )
*
      ELSE
*
*        Triangular, multi-process: redistribute into DWORK,
*        then extract the triangle via PDLACPY_I8
*
         CALL PDGEMR2D_I8( M, N, A, IA, JA, DESCA,
     $        DWORK, IB, JB, DESCB, INT( ICTXT, 8 ) )
         CALL PDLACPY_I8( UPLO, M, N, DWORK, IB, JB, DESCB,
     $        B, IB, JB, DESCB )
*
      END IF
*
      RETURN
*
*     End of PDLAMVE_I8
*
      END
