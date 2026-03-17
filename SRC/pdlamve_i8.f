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
*  The full-copy multi-process path uses PDGEMR2D_I8 natively.
*  The single-process and triangular paths narrow to default INTEGER
*  and delegate to legacy DLAMOV / PDLACPY; they abort if any
*  dimension exceeds INT_MAX.
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
     $                   DLAMOV, PDGEMR2D_I8, PDLACPY
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
*        Legacy path: narrow all dimensions to default INTEGER.
*        Abort if any exceed INT_MAX.
*
         IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $       IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $       IB.GT.INTMAX .OR. JB.GT.INTMAX .OR.
     $       DESCA( LLD_ ).GT.INTMAX .OR.
     $       DESCB( LLD_ ).GT.INTMAX ) THEN
            CALL PXERBLA( ICTXT, 'PDLAMVE_I8', -2 )
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
            CALL DLAMOV( UPLO, M4, N4,
     $           A((JA4-1)*DESCA4(9)+IA4), DESCA4(9),
     $           B((JB4-1)*DESCB4(9)+IB4), DESCB4(9) )
         ELSE
*           Triangular multi-process: redistribute full, then lacpy
            CALL PDGEMR2D_I8( M, N, A, IA, JA, DESCA,
     $           DWORK, IB, JB, DESCB, INT( ICTXT, 8 ) )
            CALL PDLACPY( UPLO, M4, N4, DWORK, IB4, JB4, DESCB4,
     $           B, IB4, JB4, DESCB4 )
         END IF
*
      ELSE
*
*        Full copy, multi-process: native I8 path
*
         CALL PDGEMR2D_I8( M, N, A, IA, JA, DESCA,
     $        B, IB, JB, DESCB, INT( ICTXT, 8 ) )
*
      END IF
*
      RETURN
*
*     End of PDLAMVE_I8
*
      END
*
*     ================================================================
*     NARROW_DESC — copy INTEGER*8 descriptor to INTEGER, aborting
*     on overflow.  Used only by the legacy delegation paths above.
*     ================================================================
*
      SUBROUTINE NARROW_DESC( DESC8, DESC4 )
      IMPLICIT NONE
      INTEGER*8          DESC8( 9 )
      INTEGER            DESC4( 9 )
      INTEGER            K
      INTRINSIC          INT
*
      DO 10 K = 1, 9
         DESC4( K ) = INT( DESC8( K ) )
   10 CONTINUE
*
      END
