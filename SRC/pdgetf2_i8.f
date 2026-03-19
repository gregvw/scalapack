      SUBROUTINE PDGETF2_I8( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PDGETF2.
*     All public integer arguments are INTEGER*8 except IPIV and INFO
*     which stay default INTEGER.
*
*     .. Scalar Arguments ..
      INTEGER*8          M, N, IA, JA
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGETF2_I8 computes an LU factorization of a general M-by-N
*  distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) using
*  partial pivoting with row interchanges.
*
*  The factorization has the form sub( A ) = P * L * U, where P is a
*  permutation matrix, L is lower triangular with unit diagonal
*  elements (lower trapezoidal if m > n), and U is upper triangular
*  (upper trapezoidal if m < n).
*
*  This is the right-looking Parallel Level 2 BLAS version of the
*  algorithm.
*
*  This routine requires N <= NB_A-MOD(JA-1, NB_A) and square block
*  decomposition ( MB_A = NB_A ).
*
*  Arguments
*  =========
*
*  M       (global input) INTEGER*8
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER*8
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ).
*          NB_A-MOD(JA-1, NB_A) >= N >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the M-by-N
*          distributed matrix sub( A ). On exit, this array contains
*          the local pieces of the factors L and U from the factoriza-
*          tion sub( A ) = P*L*U; the unit diagonal elements of L are
*          not stored.
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
*  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
*          This array contains the pivoting information.
*          IPIV(i) -> The global row local row i was swapped with.
*          This array is tied to the distributed matrix A.
*
*  INFO    (local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = K, U(IA+K-1,JA+K-1) is exactly zero.
*                The factorization has been completed, but the factor U
*                is exactly singular, and division by zero will occur if
*                it is used to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          ROWBTOP
      INTEGER            IACOL, IAROW, ICTXT, MYCOL, MYROW,
     $                   NPCOL, NPROW
      INTEGER*8          I, ICOFF, IIA, IROFF, J, JJA, MN, IPIV8
      DOUBLE PRECISION   GMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, CHK1MAT_I8,
     $                   IGEBR2D, IGEBS2D, INFOG2L_I8,
     $                   PDAMAX_I8, PDGER_I8,
     $                   PDSCAL_I8, PDSWAP_I8, PB_TOPGET, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
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
         CALL CHK1MAT_I8( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
         IF( INFO.EQ.0 ) THEN
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, DESCA( NB_ ) )
            IF( N+ICOFF.GT.DESCA( NB_ ) ) THEN
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
         CALL PXERBLA( ICTXT, 'PDGETF2_I8', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      MN = MIN( M, N )
      CALL INFOG2L_I8( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIA, JJA, IAROW, IACOL )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
*
      IF( MYCOL.EQ.IACOL ) THEN
         DO 10 J = JA, JA+MN-1
            I = IA + J - JA
*
*           Find pivot and test for singularity.
*
            CALL PDAMAX_I8( M-J+JA, GMAX, IPIV( INT(IIA+J-JA) ),
     $                      A, I, J, DESCA, 1_8 )
            IF( GMAX.NE.ZERO ) THEN
*
*              Apply the row interchanges to columns JA:JA+N-1
*
               IPIV8 = IPIV( INT(IIA+J-JA) )
               CALL PDSWAP_I8( N, A, I, JA, DESCA, DESCA( M_ ), A,
     $                         IPIV8, JA, DESCA, DESCA( M_ ) )
*
*              Compute elements I+1:IA+M-1 of J-th column.
*
               IF( J-JA+1.LT.M )
     $            CALL PDSCAL_I8( M-J+JA-1, ONE / GMAX, A, I+1, J,
     $                            DESCA, 1_8 )
            ELSE IF( INFO.EQ.0 ) THEN
               INFO = INT( J - JA + 1 )
            END IF
*
*           Update trailing submatrix
*
            IF( J-JA+1.LT.MN ) THEN
               CALL PDGER_I8( M-J+JA-1, N-J+JA-1, -ONE, A, I+1, J,
     $                        DESCA, 1_8, A, I, J+1, DESCA,
     $                        DESCA( M_ ), A, I+1, J+1, DESCA )
            END IF
   10    CONTINUE
*
         CALL IGEBS2D( ICTXT, 'Rowwise', ROWBTOP, INT(MN), 1,
     $                 IPIV( INT(IIA) ), INT(MN) )
*
      ELSE
*
         CALL IGEBR2D( ICTXT, 'Rowwise', ROWBTOP, INT(MN), 1,
     $                 IPIV( INT(IIA) ), INT(MN), MYROW, IACOL )
*
      END IF
*
      RETURN
*
*     End of PDGETF2_I8
*
      END
