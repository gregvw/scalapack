      SUBROUTINE PZPOTRF_I8( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PZPOTRF.
*     All public integer arguments (N, IA, JA, DESCA) are INTEGER*8;
*     INFO stays INTEGER.
*     Uses I8 calls for PZHERK and PZTRSM; PZPOTF2 stays legacy
*     (narrowed at call boundary via NARROW_DESC8).
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          N, IA, JA
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZPOTRF_I8 computes the Cholesky factorization of an N-by-N complex
*  hermitian positive definite distributed matrix sub( A ) denoting
*  A(IA:IA+N-1, JA:JA+N-1).
*
*  The factorization has the form
*
*            sub( A ) = U' * U ,  if UPLO = 'U', or
*
*            sub( A ) = L  * L',  if UPLO = 'L',
*
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This routine requires square block decomposition ( MB_A = NB_A ).
*
*  Arguments
*  =========
*
*  UPLO    (global input) CHARACTER
*          = 'U':  Upper triangle of sub( A ) is stored;
*          = 'L':  Lower triangle of sub( A ) is stored.
*
*  N       (global input) INTEGER*8
*          The order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*
*  IA      (global input) INTEGER*8
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER*8
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER*8 array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  illegal value
*          > 0:  not positive definite
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            ICTXT, JB, MYCOL, MYROW, NPCOL, NPROW
      INTEGER*8          I, ICOFF, IROFF, J, JN
*     ..
*     .. Local Arrays ..
      INTEGER*8          IDUM1( 1 )
      INTEGER            IDUM2( 1 ), DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT_I8, NARROW_DESC8,
     $                   PB_TOPGET, PB_TOPSET, PCHK1MAT_I8, PXERBLA,
     $                   PZPOTF2, PZHERK_I8, PZTRSM_I8
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, INT, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = INT( DESCA( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(600+CTXT_)
      ELSE
         CALL CHK1MAT_I8( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
         UPPER = LSAME( UPLO, 'U' )
         IF( INFO.EQ.0 ) THEN
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, DESCA( NB_ ) )
            IF ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( IROFF.NE.0 ) THEN
               INFO = -4
            ELSE IF( ICOFF.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(600+NB_)
            END IF
         END IF
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
         CALL PCHK1MAT_I8( N, 2, N, 2, IA, JA, DESCA, 6, 1, IDUM1,
     $                     IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZPOTRF_I8', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Range check: PZPOTF2 requires narrowing to INTEGER.
*     N, IA, JA must fit in INTEGER for the POTF2 calls.
*
      IF( N.GT.INTMAX .OR. IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         INFO = -1
         CALL PXERBLA( ICTXT, 'PZPOTRF_I8', 1 )
         RETURN
      END IF
*
*     Narrow descriptor for PZPOTF2 calls (legacy INTEGER interface)
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
*
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      IF( UPPER ) THEN
*
*        Split-ring topology for the communication along process
*        columns, 1-tree topology along process rows.
*
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ' ' )
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', 'S-ring' )
*
*        A is upper triangular, compute Cholesky factorization A = U'*U.
*
*        Handle the first block of columns separately
*
         JN = MIN( ( (JA + DESCA(NB_) - 1) / DESCA(NB_) )
     $        * DESCA(NB_), JA+N-1 )
         JB = INT( JN - JA + 1 )
*
*        Perform unblocked Cholesky factorization on JB block
*
         CALL PZPOTF2( UPLO, JB, A, INT( IA ), INT( JA ), DESCA4,
     $                 INFO )
         IF( INFO.NE.0 )
     $      GO TO 30
*
         IF( JB+1.LE.N ) THEN
*
*           Form the row panel of U using the triangular solver
*
            CALL PZTRSM_I8( 'Left', UPLO, 'Conjugate transpose',
     $                      'Non-Unit', INT( JB, 8 ), N-JB, CONE, A,
     $                      IA, JA, DESCA, A, IA, JA+JB, DESCA )
*
*           Update the trailing matrix, A = A - U'*U
*
            CALL PZHERK_I8( UPLO, 'Conjugate transpose', N-JB,
     $                      INT( JB, 8 ), -ONE, A, IA, JA+JB, DESCA,
     $                      ONE, A, IA+JB, JA+JB, DESCA )
         END IF
*
*        Loop over remaining block of columns
*
         DO 10 J = JN+1, JA+N-1, DESCA( NB_ )
            JB = INT( MIN( N-J+JA, DESCA( NB_ ) ) )
            I = IA + J - JA
*
*           Perform unblocked Cholesky factorization on JB block
*
            CALL PZPOTF2( UPLO, JB, A, INT( I ), INT( J ), DESCA4,
     $                    INFO )
            IF( INFO.NE.0 ) THEN
               INFO = INFO + INT( J - JA )
               GO TO 30
            END IF
*
            IF( J-JA+JB+1.LE.N ) THEN
*
*              Form the row panel of U using the triangular solver
*
               CALL PZTRSM_I8( 'Left', UPLO, 'Conjugate transpose',
     $                         'Non-Unit', INT( JB, 8 ), N-J-JB+JA,
     $                         CONE, A, I, J, DESCA, A, I, J+JB,
     $                         DESCA )
*
*              Update the trailing matrix, A = A - U'*U
*
               CALL PZHERK_I8( UPLO, 'Conjugate transpose', N-J-JB+JA,
     $                         INT( JB, 8 ), -ONE, A, I, J+JB, DESCA,
     $                         ONE, A, I+JB, J+JB, DESCA )
            END IF
   10    CONTINUE
*
      ELSE
*
*        1-tree topology for the communication along process columns,
*        Split-ring topology along process rows.
*
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'S-ring' )
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
*
*        A is lower triangular, compute Cholesky factorization A = L*L'
*        (right-looking)
*
*        Handle the first block of columns separately
*
         JN = MIN( ( (JA + DESCA(NB_) - 1) / DESCA(NB_) )
     $        * DESCA( NB_ ), JA+N-1 )
         JB = INT( JN - JA + 1 )
*
*        Perform unblocked Cholesky factorization on JB block
*
         CALL PZPOTF2( UPLO, JB, A, INT( IA ), INT( JA ), DESCA4,
     $                 INFO )
         IF( INFO.NE.0 )
     $      GO TO 30
*
         IF( JB+1.LE.N ) THEN
*
*           Form the column panel of L using the triangular solver
*
            CALL PZTRSM_I8( 'Right', UPLO, 'Conjugate transpose',
     $                      'Non-Unit', N-JB, INT( JB, 8 ), CONE, A,
     $                      IA, JA, DESCA, A, IA+JB, JA, DESCA )
*
*           Update the trailing matrix, A = A - L*L'
*
            CALL PZHERK_I8( UPLO, 'No Transpose', N-JB, INT( JB, 8 ),
     $                      -ONE, A, IA+JB, JA, DESCA, ONE, A, IA+JB,
     $                      JA+JB, DESCA )
*
         END IF
*
         DO 20 J = JN+1, JA+N-1, DESCA( NB_ )
            JB = INT( MIN( N-J+JA, DESCA( NB_ ) ) )
            I = IA + J - JA
*
*           Perform unblocked Cholesky factorization on JB block
*
            CALL PZPOTF2( UPLO, JB, A, INT( I ), INT( J ), DESCA4,
     $                    INFO )
            IF( INFO.NE.0 ) THEN
               INFO = INFO + INT( J - JA )
               GO TO 30
            END IF
*
            IF( J-JA+JB+1.LE.N ) THEN
*
*              Form the column panel of L using the triangular solver
*
               CALL PZTRSM_I8( 'Right', UPLO, 'Conjugate transpose',
     $                         'Non-Unit', N-J-JB+JA, INT( JB, 8 ),
     $                         CONE, A, I, J, DESCA, A, I+JB, J,
     $                         DESCA )
*
*              Update the trailing matrix, A = A - L*L'
*
               CALL PZHERK_I8( UPLO, 'No Transpose', N-J-JB+JA,
     $                         INT( JB, 8 ), -ONE, A, I+JB, J, DESCA,
     $                         ONE, A, I+JB, J+JB, DESCA )
*
            END IF
   20    CONTINUE
*
      END IF
*
   30 CONTINUE
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      RETURN
*
*     End of PZPOTRF_I8
*
      END
