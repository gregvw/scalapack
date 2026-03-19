      SUBROUTINE PCGETRF_I8( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PCGETRF.
*     All public integer arguments are INTEGER*8 except IPIV and INFO.
*     Uses I8 PBLAS calls (PCTRSM_I8, PCGEMM_I8, PCLASWP_I8).
*     Uses PCGETF2_I8 for panel factorization.
*
*     .. Scalar Arguments ..
      INTEGER*8          M, N, IA, JA
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      INTEGER            IPIV( * )
      COMPLEX            A( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      COMPLEX            ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          COLBTOP, COLCTOP, ROWBTOP
      INTEGER            ICTXT, IINFO, MYCOL, MYROW, NPCOL, NPROW
      INTEGER            JB
      INTEGER*8          I, ICOFF, IROFF, J, JN, MN, IN, NB8
*     ..
*     .. Local Arrays ..
      INTEGER*8          IDUM1( 1 )
      INTEGER            IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT_I8, IGAMN2D,
     $                   PCHK1MAT_I8, PB_TOPGET, PB_TOPSET,
     $                   PCGEMM_I8, PCGETF2_I8, PCLASWP_I8,
     $                   PCTRSM_I8, PXERBLA
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
         CALL PXERBLA( ICTXT, 'PCGETRF_I8', -INFO )
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
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'S-ring' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', ' ' )
*
      MN = MIN( M, N )
      IN = MIN( ((IA + DESCA(MB_) - 1) / DESCA(MB_)) * DESCA(MB_),
     $          IA+M-1 )
      JN = MIN( ((JA + NB8 - 1) / NB8) * NB8, JA+MN-1 )
      JB = INT( JN - JA + 1 )
*
      CALL PCGETF2_I8( M, INT( JB, 8 ), A, IA, JA, DESCA,
     $              IPIV, INFO )
*
      IF( JB+1.LE.N ) THEN
*
         CALL PCLASWP_I8( 'Forward', 'Rows', N-INT( JB, 8 ),
     $                    A, IA, JN+1, DESCA, IA, IN, IPIV )
*
         CALL PCTRSM_I8( 'Left', 'Lower', 'No transpose', 'Unit',
     $                   INT( JB, 8 ), N-INT( JB, 8 ), ONE,
     $                   A, IA, JA, DESCA, A, IA, JN+1, DESCA )
*
         IF( JB+1.LE.M ) THEN
*
            CALL PCGEMM_I8( 'No transpose', 'No transpose',
     $                      M-INT( JB, 8 ), N-INT( JB, 8 ),
     $                      INT( JB, 8 ), -ONE,
     $                      A, IN+1, JA, DESCA,
     $                      A, IA, JN+1, DESCA,
     $                      ONE, A, IN+1, JN+1, DESCA )
*
         END IF
      END IF
*
      DO 10 J = JN+1, JA+MN-1, NB8
         JB = INT( MIN( MN-J+JA, NB8 ) )
         I = IA + J - JA
*
         CALL PCGETF2_I8( M-J+JA, INT( JB, 8 ), A, I, J,
     $                 DESCA, IPIV, IINFO )
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO + INT( J - JA )
*
         CALL PCLASWP_I8( 'Forward', 'Rowwise', J-JA,
     $                    A, IA, JA, DESCA, I, I+INT( JB, 8 )-1,
     $                    IPIV )
*
         IF( J-JA+JB+1.LE.N ) THEN
*
            CALL PCLASWP_I8( 'Forward', 'Rowwise',
     $                       N-J-INT( JB, 8 )+JA,
     $                       A, IA, J+INT( JB, 8 ), DESCA,
     $                       I, I+INT( JB, 8 )-1, IPIV )
*
            CALL PCTRSM_I8( 'Left', 'Lower', 'No transpose', 'Unit',
     $                      INT( JB, 8 ), N-J-INT( JB, 8 )+JA,
     $                      ONE, A, I, J, DESCA,
     $                      A, I, J+INT( JB, 8 ), DESCA )
*
            IF( J-JA+JB+1.LE.M ) THEN
*
               CALL PCGEMM_I8( 'No transpose', 'No transpose',
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
*     End of PCGETRF_I8
*
      END
