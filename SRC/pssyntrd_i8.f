      SUBROUTINE PSSYNTRD_I8( UPLO, N, A, IA, JA, DESCA, D, E, TAU,
     $                        WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     INTEGER*8 bridge version of PDSYNTRD.
*
*  Reduces a real symmetric distributed matrix to tridiagonal form.
*  Uses I8 descriptors and dimensions throughout, with checked
*  narrowing at PBLAS and serial LAPACK call boundaries.
*
*  BRIDGE IMPLEMENTATION NOTE:
*  PSLATRD, PSSYR2K, PSSYTD2, PSSYTTRD, and SSYTRD are called via
*  narrowed INTEGER arguments.  This routine aborts cleanly if any
*  PBLAS/LAPACK-facing quantity exceeds default INTEGER range.
*  Once native I8 PBLAS exists, these boundaries can be removed.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          N, IA, JA, LWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      REAL   A( * ), D( * ), E( * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      REAL   ONE
      PARAMETER          ( ONE = 1.0E+0 )
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      CHARACTER          COLCTOP, ROWCTOP
      INTEGER*8          IPW, INDB, INDD, INDE, INDTAU, INDW,
     $                   LWMIN, TTLWMIN, LLWORK, ONEPMIN,
     $                   NP, NQ, NPS8, I8, J8, K8, JX8
      INTEGER            ANB, CTXTB, IACOL, IAROW, ICOFFA, ICTXT,
     $                   IINFO, IROFFA, JB, KK, MINSZ,
     $                   MYCOL, MYCOLB, MYROW, MYROWB, NB, NPCOL,
     $                   NPCOLB, NPROW, NPROWB, SQNPC
*     .. Narrowed locals for PBLAS/LAPACK bridge ..
      INTEGER            N4, NPS4, DESCA4( 9 ), DESCB4( 9 ),
     $                   DESCW4( 9 ), NP4, K4, I4, J4
      DOUBLE PRECISION   DLLWORK
*     ..
*     .. Local Arrays ..
      INTEGER*8          DESCB( DLEN_ ), DESCW( DLEN_ ),
     $                   IDUM1( 2 )
      INTEGER            IDUM2( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GET, BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $                   BLACS_GRIDINIT, BLACS_ABORT,
     $                   CHK1MAT_I8, DESCSET_I8, SSYTRD,
     $                   DGAMN2D, PCHK1MAT_I8, PSELSET_I8,
     $                   PSLAMR1D_I8, PSLATRD_I8, PSSYR2K, PSSYTD2,
     $                   PSSYTTRD, PSTRMR2D_I8,
     $                   PB_TOPGET, PB_TOPSET, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER*8          NUMROC_I8
      INTEGER            INDXG2P_I8, PJLAENV
      INTEGER*8          INDXG2L_I8
      EXTERNAL           LSAME, NUMROC_I8, INDXG2P_I8, INDXG2L_I8,
     $                   PJLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, ICHAR, INT, MAX, MIN, MOD, SQRT
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
         INFO = -( 600+CTXT_ )
      ELSE
         CALL CHK1MAT_I8( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
         UPPER = LSAME( UPLO, 'U' )
         IF( INFO.EQ.0 ) THEN
            NB = INT( DESCA( NB_ ) )
            IROFFA = INT( MOD( IA-1, DESCA( MB_ ) ) )
            ICOFFA = INT( MOD( JA-1, DESCA( NB_ ) ) )
            IAROW = INDXG2P_I8( IA, DESCA( NB_ ), MYROW,
     $                          INT( DESCA( RSRC_ ) ), NPROW )
            IACOL = INDXG2P_I8( JA, DESCA( NB_ ), MYCOL,
     $                          INT( DESCA( CSRC_ ) ), NPCOL )
            NP = NUMROC_I8( N, INT( NB, 8 ), MYROW, IAROW, NPROW )
            NQ = MAX( 1_8, NUMROC_I8( N+JA-1, INT( NB, 8 ), MYCOL,
     $           INT( DESCA( CSRC_ ) ), NPCOL ) )
            LWMIN = MAX( ( NP+1 )*NB, 3_8*NB )
            ANB = PJLAENV( ICTXT, 3, 'PSSYTTRD', 'L', 0, 0, 0, 0 )
            MINSZ = PJLAENV( ICTXT, 5, 'PSSYTTRD', 'L', 0, 0, 0, 0 )
            SQNPC = INT( SQRT( REAL( NPROW*NPCOL ) ) )
            NPS8 = MAX( NUMROC_I8( N, 1_8, 0, 0, SQNPC ),
     $                  INT( 2*ANB, 8 ) )
            TTLWMIN = 2*( ANB+1 )*( 4*NPS8+2 ) + ( NPS8+4 )*NPS8
*
            WORK( 1 ) = REAL( DBLE( TTLWMIN ) )
            LQUERY = ( LWORK.EQ.-1 )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( IROFFA.NE.ICOFFA .OR. ICOFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 600+NB_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -11
            END IF
         END IF
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 11
         CALL PCHK1MAT_I8( N, 2, N, 2, IA, JA, DESCA, 6, 2,
     $                     IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSSYNTRD_I8', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      ONEPMIN = N*N + 3*N + 1
      LLWORK = LWORK
*
*     Global minimum workspace via DGAMN2D (I8-safe via DOUBLE)
*
      DLLWORK = DBLE( LLWORK )
      CALL DGAMN2D( ICTXT, 'A', ' ', 1, 1, DLLWORK, 1, 1, -1,
     $              -1, -1, -1 )
      LLWORK = INT( DLLWORK, 8 )
*
*     Use the serial or tailored path if possible
*
      NPROWB = 0
      IF( ( N.LT.MINSZ .OR. SQNPC.EQ.1 ) .AND.
     $    LLWORK.GE.ONEPMIN .AND. .NOT.UPPER ) THEN
         NPROWB = 1
         NPS8 = N
      ELSE
         IF( LLWORK.GE.TTLWMIN .AND. .NOT.UPPER ) THEN
            NPROWB = SQNPC
         END IF
      END IF
*
      IF( NPROWB.GE.1 ) THEN
*
*        Serial/tailored path: redistribute to smaller grid
*
         NPCOLB = NPROWB
         SQNPC = NPROWB
         INDB = 1
         INDD = INDB + NPS8*NPS8
         INDE = INDD + NPS8
         INDTAU = INDE + NPS8
         INDW = INDTAU + NPS8
         LLWORK = LLWORK - INDW + 1
*
         CALL BLACS_GET( ICTXT, 10, CTXTB )
         CALL BLACS_GRIDINIT( CTXTB, 'Row major', SQNPC, SQNPC )
         CALL BLACS_GRIDINFO( CTXTB, NPROWB, NPCOLB, MYROWB,
     $                        MYCOLB )
         CALL DESCSET_I8( DESCB, N, N, 1_8, 1_8, 0, 0, CTXTB,
     $                    NPS8 )
*
*        Redistribute A to work buffer using I8 redistribution
*
         CALL PSTRMR2D_I8( UPLO, 'N', N, N, A, IA, JA, DESCA,
     $                     WORK( INDB ), 1_8, 1_8, DESCB,
     $                     INT( ICTXT, 8 ) )
*
*        Only processors in context CTXTB do the reduction
*
         IF( NPROWB.GT.0 ) THEN
*
*           Bridge narrowing: SSYTRD and PSSYTTRD take INTEGER
*
            IF( NPS8.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $          LLWORK.GT.INTMAX ) THEN
               CALL PXERBLA( ICTXT, 'PSSYNTRD_I8', -2 )
               CALL BLACS_ABORT( ICTXT, 1 )
            END IF
            NPS4 = INT( NPS8 )
            N4   = INT( N )
*
*           Narrow DESCB for PSSYTTRD
*
            CALL NARROW_DESC8( DESCB, DESCB4 )
*
            IF( NPROWB.EQ.1 ) THEN
               CALL SSYTRD( UPLO, N4, WORK( INDB ), NPS4,
     $                      WORK( INDD ), WORK( INDE ),
     $                      WORK( INDTAU ), WORK( INDW ),
     $                      INT( LLWORK ), IINFO )
            ELSE
               CALL PSSYTTRD( 'L', N4, WORK( INDB ), 1, 1, DESCB4,
     $                        WORK( INDD ), WORK( INDE ),
     $                        WORK( INDTAU ), WORK( INDW ),
     $                        INT( LLWORK ), IINFO )
            END IF
         END IF
*
*        Redistribute results back using I8 redistribution
*
         CALL PSLAMR1D_I8( N-1, WORK( INDE ), 1_8, 1_8, DESCB,
     $                     E, 1_8, JA, DESCA )
         CALL PSLAMR1D_I8( N, WORK( INDD ), 1_8, 1_8, DESCB,
     $                     D, 1_8, JA, DESCA )
         CALL PSLAMR1D_I8( N, WORK( INDTAU ), 1_8, 1_8, DESCB,
     $                     TAU, 1_8, JA, DESCA )
*
         CALL PSTRMR2D_I8( UPLO, 'N', N, N, WORK( INDB ), 1_8, 1_8,
     $                     DESCB, A, IA, JA, DESCA,
     $                     INT( ICTXT, 8 ) )
*
         IF( MYROWB.GE.0 )
     $      CALL BLACS_GRIDEXIT( CTXTB )
*
      ELSE
*
*        Blocked reduction path — requires PBLAS bridge narrowing
*
*        Check that all dimensions fit in default INTEGER
*
         IF( N.GT.INTMAX .OR. IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $       NP.GT.INTMAX .OR. NQ.GT.INTMAX .OR.
     $       DESCA( LLD_ ).GT.INTMAX ) THEN
            CALL PXERBLA( ICTXT, 'PSSYNTRD_I8', -2 )
            CALL BLACS_ABORT( ICTXT, 1 )
         END IF
*
         N4  = INT( N )
         NP4 = INT( NP )
*
*        Narrow DESCA for PBLAS calls
*
         CALL NARROW_DESC8( DESCA, DESCA4 )
*
         CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
         CALL PB_TOPGET( ICTXT, 'Combine', 'Rowwise', ROWCTOP )
         CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', '1-tree' )
         CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise', '1-tree' )
*
         IPW = NP*NB + 1
*
         IF( UPPER ) THEN
*
*           Reduce the upper triangle
*
            KK = INT( MOD( JA+N-1, INT( NB, 8 ) ) )
            IF( KK.EQ.0 ) KK = NB
            CALL DESCSET_I8( DESCW, N, INT( NB, 8 ), INT( NB, 8 ),
     $           INT( NB, 8 ), IAROW,
     $           INDXG2P_I8( JA+N-KK, INT( NB, 8 ), MYCOL,
     $           INT( DESCA( CSRC_ ) ), NPCOL ),
     $           ICTXT, MAX( 1_8, NP ) )
            CALL NARROW_DESC8( DESCW, DESCW4 )
*
            DO 10 K8 = N - KK + 1, INT( NB, 8 ) + 1, -INT( NB, 8 )
               JB = INT( MIN( N-K8+1, INT( NB, 8 ) ) )
               K4 = INT( K8 )
               I4 = INT( IA + K8 - 1 )
               J4 = INT( JA + K8 - 1 )
*
               CALL PSLATRD_I8( UPLO, K8+JB-1, JB, A, IA,
     $              JA, DESCA, D, E, TAU, WORK, 1_8, 1_8,
     $              DESCW, WORK( IPW ) )
*
               CALL PSSYR2K( UPLO, 'No transpose', K4-1, JB,
     $              -ONE, A, INT( IA ), J4, DESCA4,
     $              WORK, 1, 1, DESCW4, ONE, A, INT( IA ),
     $              INT( JA ), DESCA4 )
*
               JX8 = MIN( INDXG2L_I8( JA+K8-1, INT( NB, 8 ), 0,
     $              IACOL, NPCOL ), NQ )
               CALL PSELSET_I8( A, IA+K8-2, JA+K8-1, DESCA,
     $              E( JX8 ) )
*
               DESCW( CSRC_ ) = MOD( DESCW( CSRC_ )+NPCOL-1, NPCOL )
               CALL NARROW_DESC8( DESCW, DESCW4 )
*
   10       CONTINUE
*
            CALL PSSYTD2( UPLO, MIN( N4, NB ), A, INT( IA ),
     $           INT( JA ), DESCA4, D, E, TAU, WORK,
     $           INT( MIN( LWORK, INTMAX ) ), IINFO )
*
         ELSE
*
*           Reduce the lower triangle
*
            KK = INT( MOD( JA+N-1, INT( NB, 8 ) ) )
            IF( KK.EQ.0 ) KK = NB
            CALL DESCSET_I8( DESCW, N, INT( NB, 8 ), INT( NB, 8 ),
     $           INT( NB, 8 ), IAROW, IACOL, ICTXT, MAX( 1_8, NP ) )
            CALL NARROW_DESC8( DESCW, DESCW4 )
*
            DO 20 K8 = 1, N - NB, INT( NB, 8 )
               I4 = INT( IA + K8 - 1 )
               J4 = INT( JA + K8 - 1 )
               K4 = INT( K8 )
*
               CALL PSLATRD_I8( UPLO, N-K8+1, NB, A, IA+K8-1, JA+K8-1,
     $              DESCA, D, E, TAU, WORK, K8, 1_8, DESCW,
     $              WORK( IPW ) )
*
               CALL PSSYR2K( UPLO, 'No transpose', INT( N-K8-NB+1 ),
     $              NB, -ONE, A, I4+NB, J4, DESCA4, WORK, K4+NB, 1,
     $              DESCW4, ONE, A, I4+NB, J4+NB, DESCA4 )
*
               JX8 = MIN( INDXG2L_I8( JA+K8+NB-2, INT( NB, 8 ), 0,
     $              IACOL, NPCOL ), NQ )
               CALL PSELSET_I8( A, IA+K8+NB-1, JA+K8+NB-2, DESCA,
     $              E( JX8 ) )
*
               DESCW( CSRC_ ) = MOD( DESCW( CSRC_ )+1, NPCOL )
               CALL NARROW_DESC8( DESCW, DESCW4 )
*
   20       CONTINUE
*
            CALL PSSYTD2( UPLO, KK, A, INT( IA+K8-1 ),
     $                    INT( JA+K8-1 ), DESCA4, D, E, TAU,
     $                    WORK, INT( MIN( LWORK, INTMAX ) ), IINFO )
         END IF
*
         CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
         CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise', ROWCTOP )
*
      END IF
*
      WORK( 1 ) = REAL( DBLE( TTLWMIN ) )
*
      RETURN
*
*     End of PSSYNTRD_I8
*
      END
