      SUBROUTINE PDGETRF_I8( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PDGETRF.
*     All public integer arguments are INTEGER*8 except IPIV and INFO.
*     Uses I8 PBLAS calls (PDTRSM_I8, PDGEMM_I8, PDLASWP_I8).
*     PDGETF2 remains legacy (operates on a single block, M-J+JA <= NB).
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
*  PDGETRF_I8 computes an LU factorization of a general M-by-N
*  distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) using
*  partial pivoting with row interchanges.
*
*  The factorization has the form sub( A ) = P * L * U.
*
*  This is the INTEGER*8 version.  All integer arguments except
*  IPIV and INFO are INTEGER*8.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          COLBTOP, COLCTOP, ROWBTOP
      INTEGER            ICTXT, IINFO, MYCOL, MYROW, NPCOL, NPROW
      INTEGER            JB
      INTEGER*8          I, ICOFF, IROFF, J, JN, MN, IN, NB8
*     ..
*     .. Local Arrays ..
      INTEGER*8          IDUM1( 1 )
      INTEGER            IDUM2( 1 ), DESCA4( 9 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT_I8, IGAMN2D,
     $                   PCHK1MAT_I8, PB_TOPGET, PB_TOPSET,
     $                   PDGEMM_I8, PDGETF2, PDLASWP_I8, PDTRSM_I8,
     $                   PXERBLA, NARROW_DESC8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MIN, MOD
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
         CALL CHK1MAT_I8( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
         NB8 = DESCA( NB_ )
         IF( INFO.EQ.0 ) THEN
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, NB8 )
            IF( IROFF.NE.0 ) THEN
               INFO = -4
            ELSE IF( ICOFF.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.NB8 ) THEN
               INFO = -(600+NB_)
            END IF
         END IF
         CALL PCHK1MAT_I8( M, 1, N, 2, IA, JA, DESCA, 6, 0, IDUM1,
     $                      IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDGETRF_I8', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( DESCA( M_ ).EQ.1 ) THEN
         IPIV( 1 ) = 1
         RETURN
      ELSE IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         RETURN
      END IF
*
*     Range check: M, N, IA, JA must fit in INTEGER for PDGETF2.
*
      IF( M.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $    IA.GT.INTMAX .OR. JA.GT.INTMAX ) THEN
         INFO = -1
         CALL PXERBLA( ICTXT, 'PDGETRF_I8', 1 )
         RETURN
      END IF
*
*     Prepare narrowed descriptor for PDGETF2 (legacy INTEGER).
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
*
*     Split-ring topology for the communication along process rows
*
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'S-ring' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', ' ' )
*
*     Handle the first block of columns separately
*     ICEIL(JA, NB) = (JA + NB - 1) / NB, computed in I8
*
      MN = MIN( M, N )
      IN = MIN( ((IA + DESCA(MB_) - 1) / DESCA(MB_)) * DESCA(MB_),
     $          IA+M-1 )
      JN = MIN( ((JA + NB8 - 1) / NB8) * NB8, JA+MN-1 )
      JB = INT( JN - JA + 1 )
*
*     Factor diagonal and subdiagonal blocks and test for exact
*     singularity.
*
      CALL PDGETF2( INT( M ), JB, A, INT( IA ), INT( JA ), DESCA4,
     $              IPIV, INFO )
*
      IF( JB+1.LE.N ) THEN
*
*        Apply interchanges to columns JN+1:JA+N-1.
*
         CALL PDLASWP_I8( 'Forward', 'Rows', N-INT( JB, 8 ),
     $                    A, IA, JN+1, DESCA, IA, IN, IPIV )
*
*        Compute block row of U.
*
         CALL PDTRSM_I8( 'Left', 'Lower', 'No transpose', 'Unit',
     $                   INT( JB, 8 ), N-INT( JB, 8 ), ONE,
     $                   A, IA, JA, DESCA, A, IA, JN+1, DESCA )
*
         IF( JB+1.LE.M ) THEN
*
*           Update trailing submatrix.
*
            CALL PDGEMM_I8( 'No transpose', 'No transpose',
     $                      M-INT( JB, 8 ), N-INT( JB, 8 ),
     $                      INT( JB, 8 ), -ONE,
     $                      A, IN+1, JA, DESCA,
     $                      A, IA, JN+1, DESCA,
     $                      ONE, A, IN+1, JN+1, DESCA )
*
         END IF
      END IF
*
*     Loop over the remaining blocks of columns.
*
      DO 10 J = JN+1, JA+MN-1, NB8
         JB = INT( MIN( MN-J+JA, NB8 ) )
         I = IA + J - JA
*
*        Factor diagonal and subdiagonal blocks and test for exact
*        singularity.
*
         CALL PDGETF2( INT( M-J+JA ), JB, A, INT( I ), INT( J ),
     $                 DESCA4, IPIV, IINFO )
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO + INT( J - JA )
*
*        Apply interchanges to columns JA:J-JA.
*
         CALL PDLASWP_I8( 'Forward', 'Rowwise', J-JA,
     $                    A, IA, JA, DESCA, I, I+INT( JB, 8 )-1,
     $                    IPIV )
*
         IF( J-JA+JB+1.LE.N ) THEN
*
*           Apply interchanges to columns J+JB:JA+N-1.
*
            CALL PDLASWP_I8( 'Forward', 'Rowwise',
     $                       N-J-INT( JB, 8 )+JA,
     $                       A, IA, J+INT( JB, 8 ), DESCA,
     $                       I, I+INT( JB, 8 )-1, IPIV )
*
*           Compute block row of U.
*
            CALL PDTRSM_I8( 'Left', 'Lower', 'No transpose', 'Unit',
     $                      INT( JB, 8 ), N-J-INT( JB, 8 )+JA,
     $                      ONE, A, I, J, DESCA,
     $                      A, I, J+INT( JB, 8 ), DESCA )
*
            IF( J-JA+JB+1.LE.M ) THEN
*
*              Update trailing submatrix.
*
               CALL PDGEMM_I8( 'No transpose', 'No transpose',
     $                         M-J-INT( JB, 8 )+JA,
     $                         N-J-INT( JB, 8 )+JA,
     $                         INT( JB, 8 ), -ONE,
     $                         A, I+INT( JB, 8 ), J, DESCA,
     $                         A, I, J+INT( JB, 8 ), DESCA,
     $                         ONE, A, I+INT( JB, 8 ),
     $                         J+INT( JB, 8 ), DESCA )
*
            END IF
         END IF
*
   10 CONTINUE
*
      IF( INFO.EQ.0 )
     $   INFO = INT( MN ) + 1
      CALL IGAMN2D( ICTXT, 'Rowwise', ' ', 1, 1, INFO, 1, IDUM2,
     $              IDUM2, -1, -1, MYCOL )
      IF( INFO.EQ.INT( MN )+1 )
     $   INFO = 0
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
*
      RETURN
*
*     End of PDGETRF_I8
*
      END
