      SUBROUTINE PCPOSV_I8( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB,
     $                     JB, DESCB, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PCPOSV.
*     All public integer arguments are INTEGER*8 except INFO (INTEGER).
*     Calls PCPOTRF_I8 and PCPOTRS_I8 internally.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          N, NRHS, IA, JA, IB, JB
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * )
      COMPLEX            A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PCPOSV_I8 computes the solution to a complex system of linear
*  equations
*
*                       sub( A ) * X = sub( B ),
*
*  where sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1) and is an N-by-N
*  hermitian distributed positive definite matrix and X and sub( B )
*  denoting B(IB:IB+N-1,JB:JB+NRHS-1) are N-by-NRHS distributed
*  matrices.
*
*  The Cholesky decomposition is used to factor sub( A ) as
*
*                    sub( A ) = U**H * U,  if UPLO = 'U', or
*
*                    sub( A ) = L * L**H,  if UPLO = 'L',
*
*  where U is an upper triangular matrix and L is a lower triangular
*  matrix.  The factored form of sub( A ) is then used to solve the
*  system of equations.
*
*  This is the INTEGER*8 version.  All integer arguments except INFO
*  are INTEGER*8.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            IAROW, IBROW, ICTXT,
     $                   MYCOL, MYROW, NPCOL, NPROW
      INTEGER*8          ICOFFA, IROFFA, IROFFB
*     ..
*     .. Local Arrays ..
      INTEGER*8          IDUM1( 1 )
      INTEGER            IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT_I8, PCHK2MAT_I8,
     $                   PCPOTRF_I8, PCPOTRS_I8, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P_I8
      EXTERNAL           INDXG2P_I8, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(700+CTXT_)
      ELSE
         UPPER = LSAME( UPLO, 'U' )
         CALL CHK1MAT_I8( N, 2, N, 2, IA, JA, DESCA, 7, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P_I8( IA, DESCA( MB_ ), MYROW,
     $                          INT(DESCA( RSRC_ )), NPROW )
            IBROW = INDXG2P_I8( IB, DESCB( MB_ ), MYROW,
     $                          INT(DESCB( RSRC_ )), NPROW )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IF ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(700+NB_)
            ELSE IF( IROFFB.NE.0 .OR. IBROW.NE.IAROW ) THEN
               INFO = -9
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(1000+NB_)
            END IF
         END IF
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
         CALL PCHK2MAT_I8( N, 2, N, 2, IA, JA, DESCA, 7, N, 2, NRHS,
     $                     3, IB, JB, DESCB, 11, 1, IDUM1, IDUM2,
     $                     INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCPOSV_I8', -INFO )
         RETURN
      END IF
*
*     Compute the Cholesky factorization sub( A ) = U'*U or L*L'.
*
      CALL PCPOTRF_I8( UPLO, N, A, IA, JA, DESCA, INFO )
*
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system sub( A ) * X = sub( B ) overwriting sub( B )
*        with X.
*
         CALL PCPOTRS_I8( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     $                    DESCB, INFO )
*
      END IF
*
      RETURN
*
*     End of PCPOSV_I8
*
      END
