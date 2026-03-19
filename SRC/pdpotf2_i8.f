      SUBROUTINE PDPOTF2_I8( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PDPOTF2.
*     All public integer arguments are INTEGER*8; INFO stays INTEGER.
*     Local BLAS calls (DGEMV, DSCAL, DDOT) use INT() narrowing
*     for their integer arguments.
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
*  PDPOTF2_I8 computes the Cholesky factorization of a real symmetric
*  positive definite distributed matrix sub( A )=A(IA:IA+N-1,JA:JA+N-1).
*
*  The factorization has the form
*
*            sub( A ) = U' * U ,  if UPLO = 'U', or
*
*            sub( A ) = L  * L',  if UPLO = 'L',
*
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This routine requires N <= NB_A-MOD(JA-1, NB_A) and square block
*  decomposition ( MB_A = NB_A ).
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
*  INFO    (local output) INTEGER
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
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            IACOL, IAROW, ICTXT, MYCOL, MYROW,
     $                   NPCOL, NPROW
      INTEGER*8          ICOFF, ICURR, IDIAG, IIA,
     $                   IOFFA, IROFF, J, JJA, LDA
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, CHK1MAT_I8,
     $                   DGEMV, DSCAL, IGEBR2D, IGEBS2D, INFOG2L_I8,
     $                   PB_TOPGET, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MOD, SQRT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = INT( DESCA( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters.
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(600+CTXT_)
      ELSE
         CALL CHK1MAT_I8( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
         IF( INFO.EQ.0 ) THEN
            UPPER = LSAME( UPLO, 'U' )
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, DESCA( NB_ ) )
            IF ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( N+ICOFF.GT.DESCA( NB_ ) ) THEN
               INFO = -2
            ELSE IF( IROFF.NE.0 ) THEN
               INFO = -4
            ELSE IF( ICOFF.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(600+NB_)
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDPOTF2_I8', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Compute local information
*
      CALL INFOG2L_I8( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIA, JJA, IAROW, IACOL )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      LDA = DESCA( LLD_ )
*
      IF ( UPPER ) THEN
*
*        Process (IAROW, IACOL) owns block to be factorized
*
         IF( MYROW.EQ.IAROW ) THEN
            IF( MYCOL.EQ.IACOL ) THEN
*
*              Compute the Cholesky factorization A = U'*U.
*
               IDIAG = IIA + ( JJA - 1 ) * LDA
               IOFFA = IDIAG
*
               DO 10 J = JA, JA+N-1
*
*                 Compute U(J,J) and test for non-positive-definiteness.
*
                  AJJ = A( IDIAG ) -
     $                  DDOT( INT( J-JA ), A( IOFFA ), 1,
     $                        A( IOFFA ), 1 )
                  IF( AJJ.LE.ZERO ) THEN
                     A( IDIAG ) = AJJ
                     INFO = INT( J - JA + 1 )
                     GO TO 20
                  END IF
                  AJJ = SQRT( AJJ )
                  A( IDIAG ) = AJJ
*
*                 Compute elements J+1:JA+N-1 of row J.
*
                  IF( J.LT.JA+N-1 ) THEN
                     ICURR = IDIAG + LDA
                     CALL DGEMV( 'Transpose', INT( J-JA ),
     $                           INT( JA+N-J-1 ), -ONE,
     $                           A( IOFFA+LDA ), INT( LDA ),
     $                           A( IOFFA ), 1,
     $                           ONE, A( ICURR ), INT( LDA ) )
                     CALL DSCAL( INT( N-J+JA-1 ), ONE / AJJ,
     $                           A( ICURR ), INT( LDA ) )
                  END IF
                  IDIAG = IDIAG + LDA + 1
                  IOFFA = IOFFA + LDA
   10          CONTINUE
*
   20          CONTINUE
*
*              Broadcast INFO to all processes in my IAROW.
*
               CALL IGEBS2D( ICTXT, 'Rowwise', ROWBTOP, 1, 1, INFO,
     $                       1 )
*
            ELSE
*
               CALL IGEBR2D( ICTXT, 'Rowwise', ROWBTOP, 1, 1, INFO, 1,
     $                       MYROW, IACOL )
            END IF
*
*           IAROW bcasts along columns so that everyone has INFO
*
            CALL IGEBS2D( ICTXT, 'Columnwise', COLBTOP, 1, 1, INFO,
     $                    1 )
*
         ELSE
*
            CALL IGEBR2D( ICTXT, 'Columnwise', COLBTOP, 1, 1, INFO, 1,
     $                    IAROW, MYCOL )
*
         END IF
*
      ELSE
*
*        Process (IAROW, IACOL) owns block to be factorized
*
         IF( MYCOL.EQ.IACOL ) THEN
            IF( MYROW.EQ.IAROW ) THEN
*
*              Compute the Cholesky factorization A = L*L'.
*
               IDIAG = IIA + ( JJA - 1 ) * LDA
               IOFFA = IDIAG
*
               DO 30 J = JA, JA+N-1
*
*                 Compute L(J,J) and test for non-positive-definiteness.
*
                  AJJ = A( IDIAG ) -
     $                  DDOT( INT( J-JA ), A( IOFFA ), INT( LDA ),
     $                        A( IOFFA ), INT( LDA ) )
                  IF ( AJJ.LE.ZERO ) THEN
                     A( IDIAG ) = AJJ
                     INFO = INT( J - JA + 1 )
                     GO TO 40
                  END IF
                  AJJ = SQRT( AJJ )
                  A( IDIAG ) = AJJ
*
*                 Compute elements J+1:JA+N-1 of column J.
*
                  IF( J.LT.JA+N-1 ) THEN
                     ICURR = IDIAG + 1
                     CALL DGEMV( 'No transpose', INT( JA+N-J-1 ),
     $                           INT( J-JA ), -ONE,
     $                           A( IOFFA+1 ), INT( LDA ),
     $                           A( IOFFA ), INT( LDA ),
     $                           ONE, A( ICURR ), 1 )
                     CALL DSCAL( INT( JA+N-J-1 ), ONE / AJJ,
     $                           A( ICURR ), 1 )
                  END IF
                  IDIAG = IDIAG + LDA + 1
                  IOFFA = IOFFA + 1
   30          CONTINUE
*
   40          CONTINUE
*
*              Broadcast INFO to everyone in IACOL
*
               CALL IGEBS2D( ICTXT, 'Columnwise', COLBTOP, 1, 1, INFO,
     $                       1 )
*
            ELSE
*
               CALL IGEBR2D( ICTXT, 'Columnwise', COLBTOP, 1, 1, INFO,
     $                       1, IAROW, MYCOL )
*
            END IF
*
*           IACOL bcasts INFO along rows so that everyone has it
*
            CALL IGEBS2D( ICTXT, 'Rowwise', ROWBTOP, 1, 1, INFO, 1 )
*
         ELSE
*
            CALL IGEBR2D( ICTXT, 'Rowwise', ROWBTOP, 1, 1, INFO, 1,
     $                    MYROW, IACOL )
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PDPOTF2_I8
*
      END
