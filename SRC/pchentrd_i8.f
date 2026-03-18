      SUBROUTINE PCHENTRD_I8( UPLO, N, A, IA, JA, DESCA, D, E, TAU,
     $                       WORK, LWORK, RWORK, LRWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     INTEGER*8 bridge version of PCHENTRD.
*
*  Reduces a complex Hermitian distributed matrix to real tridiagonal
*  form.  Uses I8 descriptors and dimensions throughout, with checked
*  narrowing at PBLAS and serial LAPACK call boundaries.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          N, IA, JA, LWORK, LRWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      REAL               D( * ), E( * ), RWORK( * )
      COMPLEX            A( * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      CHARACTER          COLCTOP, ROWCTOP
      INTEGER*8          IPW, INDB, INDRD, INDRE, INDTAU, INDW,
     $                   LWMIN, TTLWMIN, LLWORK, ONEPMIN,
     $                   LRWMIN, TTLRWMIN, LLRWORK, ONEPRMIN,
     $                   NP, NQ, NPS8, I8, J8, K8, JX8
      INTEGER            ANB, CTXTB, IACOL, IAROW, ICOFFA, ICTXT,
     $                   IINFO, IROFFA, JB, KK, MINSZ,
     $                   MYCOL, MYCOLB, MYROW, MYROWB, NB, NPCOL,
     $                   NPCOLB, NPROW, NPROWB, SQNPC
*     .. Narrowed locals for PBLAS/LAPACK bridge ..
      INTEGER            N4, NPS4, DESCA4( 9 ), DESCB4( 9 ),
     $                   DESCW4( 9 ), NP4, K4, I4, J4
      DOUBLE PRECISION   DLLWORK, DLLRWORK
*     ..
*     .. Local Arrays ..
      INTEGER*8          DESCB( DLEN_ ), DESCW( DLEN_ ),
     $                   IDUM1( 3 )
      INTEGER            IDUM2( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GET, BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $                   BLACS_GRIDINIT, BLACS_ABORT,
     $                   CHK1MAT_I8, DESCSET_I8, CHETRD,
     $                   DGAMN2D, PCHK1MAT_I8, PCELSET_I8,
     $                   PCLAMR1D_I8, PSLAMR1D_I8,
     $                   PCLATRD_I8, PCHER2K, PCHETD2, PCHETTRD,
     $                   PCTRMR2D_I8, PSTRMR2D_I8,
     $                   PB_TOPGET, PB_TOPSET, PXERBLA,
     $                   NARROW_DESC8
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER*8          NUMROC_I8, INDXG2L_I8
      INTEGER            INDXG2P_I8, PJLAENV
      EXTERNAL           LSAME, NUMROC_I8, INDXG2P_I8, INDXG2L_I8,
     $                   PJLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, ICHAR, INT, MAX, MIN, MOD, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
      ICTXT = INT( DESCA( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
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
            ANB = PJLAENV( ICTXT, 3, 'PCHETTRD', 'L', 0, 0, 0, 0 )
            MINSZ = PJLAENV( ICTXT, 5, 'PCHETTRD', 'L', 0, 0, 0, 0 )
            SQNPC = INT( SQRT( REAL( NPROW*NPCOL ) ) )
            NPS8 = MAX( NUMROC_I8( N, 1_8, 0, 0, SQNPC ),
     $                  INT( 2*ANB, 8 ) )
            TTLWMIN = 2*( ANB+1 )*( 4*NPS8+2 ) + ( NPS8+2 )*NPS8
            LRWMIN = 1
            TTLRWMIN = 2*NPS8
*
            WORK( 1 ) = CMPLX( REAL( DBLE( TTLWMIN ) ) )
            RWORK( 1 ) = REAL( DBLE( TTLRWMIN ) )
            LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( IROFFA.NE.ICOFFA .OR. ICOFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 600+NB_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -11
            ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -13
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
         IF( LRWORK.EQ.-1 ) THEN
            IDUM1( 3 ) = -1
         ELSE
            IDUM1( 3 ) = 1
         END IF
         IDUM2( 3 ) = 13
         CALL PCHK1MAT_I8( N, 2, N, 2, IA, JA, DESCA, 6, 3,
     $                     IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCHENTRD_I8', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
*
      ONEPMIN = N*N + 3*N + 1
      LLWORK = LWORK
      DLLWORK = DBLE( LLWORK )
      CALL DGAMN2D( ICTXT, 'A', ' ', 1, 1, DLLWORK, 1, 1, -1,
     $              -1, -1, -1 )
      LLWORK = INT( DLLWORK, 8 )
*
      ONEPRMIN = 2*N
      LLRWORK = LRWORK
      DLLRWORK = DBLE( LLRWORK )
      CALL DGAMN2D( ICTXT, 'A', ' ', 1, 1, DLLRWORK, 1, 1, -1,
     $              -1, -1, -1 )
      LLRWORK = INT( DLLRWORK, 8 )
*
      NPROWB = 0
      IF( ( N.LT.MINSZ .OR. SQNPC.EQ.1 ) .AND.
     $    LLWORK.GE.ONEPMIN .AND. LLRWORK.GE.ONEPRMIN .AND.
     $    .NOT.UPPER ) THEN
         NPROWB = 1
         NPS8 = N
      ELSE
         IF( LLWORK.GE.TTLWMIN .AND. LLRWORK.GE.TTLRWMIN .AND.
     $       .NOT.UPPER ) THEN
            NPROWB = SQNPC
         END IF
      END IF
*
      IF( NPROWB.GE.1 ) THEN
*
         NPCOLB = NPROWB
         SQNPC = NPROWB
         INDB = 1
         INDRD = 1
         INDRE = INDRD + NPS8
         INDTAU = INDB + NPS8*NPS8
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
         CALL PCTRMR2D_I8( UPLO, 'N', N, N, A, IA, JA, DESCA,
     $                     WORK( INDB ), 1_8, 1_8, DESCB,
     $                     INT( ICTXT, 8 ) )
*
         IF( NPROWB.GT.0 ) THEN
*
            IF( NPS8.GT.INTMAX .OR. N.GT.INTMAX .OR.
     $          LLWORK.GT.INTMAX ) THEN
               CALL PXERBLA( ICTXT, 'PCHENTRD_I8', -2 )
               CALL BLACS_ABORT( ICTXT, 1 )
            END IF
            NPS4 = INT( NPS8 )
            N4   = INT( N )
            CALL NARROW_DESC8( DESCB, DESCB4 )
*
            IF( NPROWB.EQ.1 ) THEN
               CALL CHETRD( UPLO, N4, WORK( INDB ), NPS4,
     $                      RWORK( INDRD ), RWORK( INDRE ),
     $                      WORK( INDTAU ), WORK( INDW ),
     $                      INT( LLWORK ), IINFO )
            ELSE
               CALL PCHETTRD( 'L', N4, WORK( INDB ), 1, 1, DESCB4,
     $                        RWORK( INDRD ), RWORK( INDRE ),
     $                        WORK( INDTAU ), WORK( INDW ),
     $                        INT( LLWORK ), IINFO )
            END IF
         END IF
*
*        D, E are REAL — use PSLAMR1D_I8
*        TAU is COMPLEX — use PCLAMR1D_I8
*
         CALL PSLAMR1D_I8( N-1, RWORK( INDRE ), 1_8, 1_8, DESCB,
     $                     E, 1_8, JA, DESCA )
         CALL PSLAMR1D_I8( N, RWORK( INDRD ), 1_8, 1_8, DESCB,
     $                     D, 1_8, JA, DESCA )
         CALL PCLAMR1D_I8( N, WORK( INDTAU ), 1_8, 1_8, DESCB,
     $                     TAU, 1_8, JA, DESCA )
*
         CALL PCTRMR2D_I8( UPLO, 'N', N, N, WORK( INDB ), 1_8, 1_8,
     $                     DESCB, A, IA, JA, DESCA,
     $                     INT( ICTXT, 8 ) )
*
         IF( MYROWB.GE.0 )
     $      CALL BLACS_GRIDEXIT( CTXTB )
*
      ELSE
*
*        Blocked reduction — PBLAS bridge narrowing
*
         IF( N.GT.INTMAX .OR. IA.GT.INTMAX .OR. JA.GT.INTMAX .OR.
     $       NP.GT.INTMAX .OR. NQ.GT.INTMAX .OR.
     $       DESCA( LLD_ ).GT.INTMAX ) THEN
            CALL PXERBLA( ICTXT, 'PCHENTRD_I8', -2 )
            CALL BLACS_ABORT( ICTXT, 1 )
         END IF
*
         N4  = INT( N )
         NP4 = INT( NP )
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
               CALL PCLATRD_I8( UPLO, K8+JB-1, JB, A, IA,
     $              JA, DESCA, D, E, TAU, WORK, 1_8, 1_8,
     $              DESCW, WORK( IPW ) )
*
               CALL PCHER2K( UPLO, 'No transpose', K4-1, JB,
     $              -CONE, A, INT( IA ), J4, DESCA4,
     $              WORK, 1, 1, DESCW4, ONE, A, INT( IA ),
     $              INT( JA ), DESCA4 )
*
               JX8 = MIN( INDXG2L_I8( JA+K8-1, INT( NB, 8 ), 0,
     $              IACOL, NPCOL ), NQ )
               CALL PCELSET_I8( A, IA+K8-2, JA+K8-1, DESCA,
     $              CMPLX( E( JX8 ) ) )
*
               DESCW( CSRC_ ) = MOD( DESCW( CSRC_ )+NPCOL-1, NPCOL )
               CALL NARROW_DESC8( DESCW, DESCW4 )
*
   10       CONTINUE
*
            CALL PCHETD2( UPLO, MIN( N4, NB ), A, INT( IA ),
     $           INT( JA ), DESCA4, D, E, TAU, WORK,
     $           INT( MIN( LWORK, INTMAX ) ), IINFO )
*
         ELSE
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
               CALL PCLATRD_I8( UPLO, N-K8+1, NB, A, IA+K8-1, JA+K8-1,
     $              DESCA, D, E, TAU, WORK, K8, 1_8, DESCW,
     $              WORK( IPW ) )
*
               CALL PCHER2K( UPLO, 'No transpose', INT( N-K8-NB+1 ),
     $              NB, -CONE, A, I4+NB, J4, DESCA4, WORK, K4+NB, 1,
     $              DESCW4, ONE, A, I4+NB, J4+NB, DESCA4 )
*
               JX8 = MIN( INDXG2L_I8( JA+K8+NB-2, INT( NB, 8 ), 0,
     $              IACOL, NPCOL ), NQ )
               CALL PCELSET_I8( A, IA+K8+NB-1, JA+K8+NB-2, DESCA,
     $              CMPLX( E( JX8 ) ) )
*
               DESCW( CSRC_ ) = MOD( DESCW( CSRC_ )+1, NPCOL )
               CALL NARROW_DESC8( DESCW, DESCW4 )
*
   20       CONTINUE
*
            CALL PCHETD2( UPLO, KK, A, INT( IA+K8-1 ),
     $                    INT( JA+K8-1 ), DESCA4, D, E, TAU,
     $                    WORK, INT( MIN( LWORK, INTMAX ) ), IINFO )
         END IF
*
         CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
         CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise', ROWCTOP )
*
      END IF
*
      WORK( 1 ) = CMPLX( REAL( DBLE( TTLWMIN ) ) )
      RWORK( 1 ) = REAL( DBLE( TTLRWMIN ) )
*
      RETURN
*
*     End of PCHENTRD_I8
*
      END
