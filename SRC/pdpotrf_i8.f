      SUBROUTINE PDPOTRF_I8( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PDPOTRF.
*     All public integer arguments are INTEGER*8; INFO stays INTEGER.
*     Uses I8 calls for PDSYRK, PDTRSM, and PDPOTF2.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          N, IA, JA
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      DOUBLE PRECISION   A( * )
*     ..
*
*  Purpose
*  =======
*
*  PDPOTRF_I8 computes the Cholesky factorization of an N-by-N real
*  symmetric positive definite distributed matrix sub( A ) denoting
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
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
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
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = K, the leading minor of order K,
*                A(IA:IA+K-1,JA:JA+K-1) is not positive definite, and
*                the factorization could not be completed.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
*
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            ICTXT, MYCOL, MYROW, NB, NPCOL, NPROW
      INTEGER*8          I, ICOFF, IROFF, J, JN, NB8
      INTEGER            JB
*     ..
*     .. Local Arrays ..
      INTEGER*8          IDUM1( 1 )
      INTEGER            IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT_I8, PCHK1MAT_I8,
     $                   PB_TOPGET, PB_TOPSET, PDPOTF2_I8, PDSYRK_I8,
     $                   PDTRSM_I8, PXERBLA
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
         NB8 = DESCA( NB_ )
         NB = INT( NB8 )
         IF( INFO.EQ.0 ) THEN
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, NB8 )
            IF ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( IROFF.NE.0 ) THEN
               INFO = -4
            ELSE IF( ICOFF.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.NB8 ) THEN
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
         CALL PXERBLA( ICTXT, 'PDPOTRF_I8', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*
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
*        ICEIL(JA, NB) = (JA + NB - 1) / NB, computed in I8
*
         JN = MIN( ((JA + NB8 - 1) / NB8) * NB8, JA+N-1 )
         JB = INT( JN - JA + 1 )
*
*        Perform unblocked Cholesky factorization on JB block
*
         CALL PDPOTF2_I8( UPLO, INT( JB, 8 ), A, IA, JA, DESCA,
     $                 INFO )
         IF( INFO.NE.0 )
     $      GO TO 30
*
         IF( JB+1.LE.N ) THEN
*
*           Form the row panel of U using the triangular solver
*
            CALL PDTRSM_I8( 'Left', UPLO, 'Transpose', 'Non-Unit',
     $                      INT( JB, 8 ), N-INT( JB, 8 ), ONE,
     $                      A, IA, JA, DESCA,
     $                      A, IA, JA+INT( JB, 8 ), DESCA )
*
*           Update the trailing matrix, A = A - U'*U
*
            CALL PDSYRK_I8( UPLO, 'Transpose', N-INT( JB, 8 ),
     $                      INT( JB, 8 ), -ONE,
     $                      A, IA, JA+INT( JB, 8 ), DESCA,
     $                      ONE, A, IA+INT( JB, 8 ),
     $                      JA+INT( JB, 8 ), DESCA )
         END IF
*
*        Loop over remaining block of columns
*
         DO 10 J = JN+1, JA+N-1, NB8
            JB = INT( MIN( N-J+JA, NB8 ) )
            I = IA + J - JA
*
*           Perform unblocked Cholesky factorization on JB block
*
            CALL PDPOTF2_I8( UPLO, INT( JB, 8 ), A, I, J, DESCA,
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
               CALL PDTRSM_I8( 'Left', UPLO, 'Transpose', 'Non-Unit',
     $                         INT( JB, 8 ), N-J-INT( JB, 8 )+JA,
     $                         ONE, A, I, J, DESCA,
     $                         A, I, J+INT( JB, 8 ), DESCA )
*
*              Update the trailing matrix, A = A - U'*U
*
               CALL PDSYRK_I8( UPLO, 'Transpose',
     $                         N-J-INT( JB, 8 )+JA, INT( JB, 8 ),
     $                         -ONE, A, I, J+INT( JB, 8 ), DESCA,
     $                         ONE, A, I+INT( JB, 8 ),
     $                         J+INT( JB, 8 ), DESCA )
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
         JN = MIN( ((JA + NB8 - 1) / NB8) * NB8, JA+N-1 )
         JB = INT( JN - JA + 1 )
*
*        Perform unblocked Cholesky factorization on JB block
*
         CALL PDPOTF2_I8( UPLO, INT( JB, 8 ), A, IA, JA, DESCA,
     $                 INFO )
         IF( INFO.NE.0 )
     $      GO TO 30
*
         IF( JB+1.LE.N ) THEN
*
*           Form the column panel of L using the triangular solver
*
            CALL PDTRSM_I8( 'Right', UPLO, 'Transpose', 'Non-Unit',
     $                      N-INT( JB, 8 ), INT( JB, 8 ), ONE,
     $                      A, IA, JA, DESCA,
     $                      A, IA+INT( JB, 8 ), JA, DESCA )
*
*           Update the trailing matrix, A = A - L*L'
*
            CALL PDSYRK_I8( UPLO, 'No Transpose', N-INT( JB, 8 ),
     $                      INT( JB, 8 ), -ONE,
     $                      A, IA+INT( JB, 8 ), JA, DESCA,
     $                      ONE, A, IA+INT( JB, 8 ),
     $                      JA+INT( JB, 8 ), DESCA )
*
         END IF
*
         DO 20 J = JN+1, JA+N-1, NB8
            JB = INT( MIN( N-J+JA, NB8 ) )
            I = IA + J - JA
*
*           Perform unblocked Cholesky factorization on JB block
*
            CALL PDPOTF2_I8( UPLO, INT( JB, 8 ), A, I, J, DESCA,
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
               CALL PDTRSM_I8( 'Right', UPLO, 'Transpose', 'Non-Unit',
     $                         N-J-INT( JB, 8 )+JA, INT( JB, 8 ),
     $                         ONE, A, I, J, DESCA,
     $                         A, I+INT( JB, 8 ), J, DESCA )
*
*              Update the trailing matrix, A = A - L*L'
*
               CALL PDSYRK_I8( UPLO, 'No Transpose',
     $                         N-J-INT( JB, 8 )+JA, INT( JB, 8 ),
     $                         -ONE, A, I+INT( JB, 8 ), J, DESCA,
     $                         ONE, A, I+INT( JB, 8 ),
     $                         J+INT( JB, 8 ), DESCA )
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
*     End of PDPOTRF_I8
*
      END
