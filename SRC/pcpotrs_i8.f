      SUBROUTINE PCPOTRS_I8( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB,
     $                       JB, DESCB, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PCPOTRS.
*     All public integer arguments are INTEGER*8; INFO stays INTEGER.
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
*  PCPOTRS_I8 solves a system of linear equations
*
*                      sub( A ) * X = sub( B )
*          A(IA:IA+N-1,JA:JA+N-1)*X = B(IB:IB+N-1,JB:JB+NRHS-1)
*
*  where sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1) and is a N-by-N
*  Hermitian positive definite distributed matrix using the Cholesky
*  factorization sub( A ) = U**H*U or L*L**H computed by PCPOTRF.
*  sub( B ) denotes the distributed matrix B(IB:IB+N-1,JB:JB+NRHS-1).
*
*  This is the INTEGER*8 version of PCPOTRS.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
      COMPLEX            ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            ICTXT, MYCOL, MYROW, NPCOL, NPROW
      INTEGER*8          IAROW, IBROW, IROFFA, IROFFB, ICOFFA
*     ..
*     .. Local Arrays ..
      INTEGER*8          IDUM1( 1 )
      INTEGER            IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT_I8, PCHK2MAT_I8,
     $                   PCTRSM_I8, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P_I8
      EXTERNAL           INDXG2P_I8, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, INT, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
      ICTXT = INT( DESCA( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
*
*     Range check: reject if any public integer arg exceeds INTMAX.
*
      IF( N.GT.INTMAX .OR. NRHS.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $    IB.GT.INTMAX .OR. JB.GT.INTMAX ) THEN
         INFO = -2
         CALL PXERBLA( ICTXT, 'PCPOTRS_I8', -INFO )
         RETURN
      END IF
*
*     Test the input parameters.
*
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(700+CTXT_)
      ELSE
         CALL CHK1MAT_I8( N, 2, N, 2, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT_I8( N, 2, NRHS, 3, IB, JB, DESCB, 11, INFO )
         UPPER = LSAME( UPLO, 'U' )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P_I8( IA, DESCA( MB_ ), MYROW,
     $                          DESCA( RSRC_ ), NPROW )
            IBROW = INDXG2P_I8( IB, DESCB( MB_ ), MYROW,
     $                          DESCB( RSRC_ ), NPROW )
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
               INFO = -(1100+NB_)
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
         CALL PXERBLA( ICTXT, 'PCPOTRS_I8', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Solve sub( A ) * X = sub( B ) where sub( A ) = U'*U.
*
*        Solve U'*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PCTRSM_I8( 'Left', 'Upper', 'Conjugate transpose',
     $                 'Non-unit', N, NRHS, ONE, A, IA, JA, DESCA, B,
     $                 IB, JB, DESCB )
*
*        Solve U*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PCTRSM_I8( 'Left', 'Upper', 'No transpose', 'Non-unit',
     $                 N, NRHS, ONE, A, IA, JA, DESCA, B, IB, JB,
     $                 DESCB )
      ELSE
*
*        Solve sub( A ) *X = sub( B ) where sub( A ) = L*L'.
*
*        Solve L*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PCTRSM_I8( 'Left', 'Lower', 'No transpose', 'Non-unit',
     $                 N, NRHS, ONE, A, IA, JA, DESCA, B, IB, JB,
     $                 DESCB )
*
*        Solve L'*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PCTRSM_I8( 'Left', 'Lower', 'Conjugate transpose',
     $                 'Non-unit', N, NRHS, ONE, A, IA, JA, DESCA, B,
     $                 IB, JB, DESCB )
      END IF
*
      RETURN
*
*     End of PCPOTRS_I8
*
      END
