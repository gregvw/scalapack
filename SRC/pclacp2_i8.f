      SUBROUTINE PCLACP2_I8( UPLO, M, N, A, IA, JA, DESCA,
     $                       B, IB, JB, DESCB )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     INTEGER*8 version of PCLACP2.
*
*  PCLACP2_I8 copies all or part of a distributed matrix A to another
*  distributed matrix B.  No communication is performed; this is a
*  local copy.  Requires that only one dimension of the matrix
*  operands is distributed.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          IA, IB, JA, JB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * )
      COMPLEX            A( * ), B( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
*     ..
*     .. Local Scalars ..
      INTEGER*8          HEIGHT, IBASE, ICOFFA, ILEFT, IRIGHT, IROFFA,
     $                   ITOP, MBA, NBA, MP, MPAA, NQ, NQAA, WIDE,
     $                   IIA, IIAA, IIB, IIBB, IIBEGA, IIBEGB,
     $                   IIENDA, IINXTA, IINXTB,
     $                   JJA, JJAA, JJB, JJBB, JJBEGA, JJBEGB,
     $                   JJENDA, JJNXTA, JJNXTB, LDA, LDB, MYDIST8
      INTEGER            IACOL, IAROW, IBCOL, IBROW,
     $                   MYCOL, MYROW, NPCOL, NPROW, RSRC, CSRC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CLAMOV_I8, INFOG2L_I8
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER*8          NUMROC_I8
      EXTERNAL           LSAME, NUMROC_I8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD, INT
*     ..
*     .. Executable Statements ..
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( INT( DESCA( CTXT_ ) ), NPROW, NPCOL,
     $                     MYROW, MYCOL )
*
      CALL INFOG2L_I8( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIA, JJA, IAROW, IACOL )
      CALL INFOG2L_I8( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIB, JJB, IBROW, IBCOL )
*
      MBA    = DESCA( MB_ )
      NBA    = DESCA( NB_ )
      LDA    = DESCA( LLD_ )
      IROFFA = MOD( IA-1, MBA )
      ICOFFA = MOD( JA-1, NBA )
      LDB    = DESCB( LLD_ )
*
      IF( N.LE.( NBA-ICOFFA ) ) THEN
*
*        Local columns JJA:JJA+N-1 are in the same process column.
*
         IF( MYCOL.EQ.IACOL ) THEN
*
            MP = NUMROC_I8( M+IROFFA, MBA, MYROW, IAROW, NPROW )
            IF( MP.LE.0 )
     $         RETURN
            IF( MYROW.EQ.IAROW )
     $         MP = MP - IROFFA
            MYDIST8 = MOD( MYROW-IAROW+NPROW, NPROW )
            ITOP    = MYDIST8 * MBA - IROFFA
*
            IF( LSAME( UPLO, 'U' ) ) THEN
*
               ITOP   = MAX( 0_8, ITOP )
               IIBEGA = IIA
               IIENDA = IIA + MP - 1
               IINXTA = MIN( ((IIBEGA+MBA-1)/MBA)*MBA, IIENDA )
               IIBEGB = IIB
               IINXTB = IIBEGB + IINXTA - IIBEGA
*
   10          CONTINUE
               IF( ( N-ITOP ).GT.0 ) THEN
                  CALL CLAMOV_I8( UPLO, IINXTA-IIBEGA+1, N-ITOP,
     $                 A( IIBEGA+(JJA+ITOP-1)*LDA ), LDA,
     $                 B( IIBEGB+(JJB+ITOP-1)*LDB ), LDB )
                  MYDIST8 = MYDIST8 + NPROW
                  ITOP    = MYDIST8 * MBA - IROFFA
                  IIBEGA  = IINXTA + 1
                  IINXTA  = MIN( IINXTA+MBA, IIENDA )
                  IIBEGB  = IINXTB + 1
                  IINXTB  = IIBEGB + IINXTA - IIBEGA
                  GO TO 10
               END IF
*
            ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
               MPAA  = MP
               IIAA  = IIA
               JJAA  = JJA
               IIBB  = IIB
               JJBB  = JJB
               IBASE = MIN( ITOP + MBA, N )
               ITOP  = MIN( MAX( 0_8, ITOP ), N )
*
   20          CONTINUE
               IF( JJAA.LE.( JJA+N-1 ) ) THEN
                  HEIGHT = IBASE - ITOP
                  CALL CLAMOV_I8( 'All', MPAA, ITOP-JJAA+JJA,
     $                 A( IIAA+(JJAA-1)*LDA ), LDA,
     $                 B( IIBB+(JJBB-1)*LDB ), LDB )
                  CALL CLAMOV_I8( UPLO, MPAA, HEIGHT,
     $                 A( IIAA+(JJA+ITOP-1)*LDA ), LDA,
     $                 B( IIBB+(JJB+ITOP-1)*LDB ), LDB )
                  MPAA    = MAX( 0_8, MPAA - HEIGHT )
                  IIAA    = IIAA + HEIGHT
                  JJAA    = JJA  + IBASE
                  IIBB    = IIBB + HEIGHT
                  JJBB    = JJB  + IBASE
                  MYDIST8 = MYDIST8 + NPROW
                  ITOP    = MYDIST8 * MBA - IROFFA
                  IBASE   = MIN( ITOP + MBA, N )
                  ITOP    = MIN( ITOP, N )
                  GO TO 20
               END IF
*
            ELSE
*
               CALL CLAMOV_I8( 'All', MP, N, A( IIA+(JJA-1)*LDA ),
     $                         LDA, B( IIB+(JJB-1)*LDB ), LDB )
*
            END IF
*
         END IF
*
      ELSE IF( M.LE.( MBA-IROFFA ) ) THEN
*
*        Local rows IIA:IIA+M-1 are in the same process row.
*
         IF( MYROW.EQ.IAROW ) THEN
*
            NQ = NUMROC_I8( N+ICOFFA, NBA, MYCOL, IACOL, NPCOL )
            IF( NQ.LE.0 )
     $         RETURN
            IF( MYCOL.EQ.IACOL )
     $         NQ = NQ - ICOFFA
            MYDIST8 = MOD( MYCOL-IACOL+NPCOL, NPCOL )
            ILEFT   = MYDIST8 * NBA - ICOFFA
*
            IF( LSAME( UPLO, 'L' ) ) THEN
*
               ILEFT  = MAX( 0_8, ILEFT )
               JJBEGA = JJA
               JJENDA = JJA + NQ - 1
               JJNXTA = MIN( ((JJBEGA+NBA-1)/NBA)*NBA, JJENDA )
               JJBEGB = JJB
               JJNXTB = JJBEGB + JJNXTA - JJBEGA
*
   30          CONTINUE
               IF( ( M-ILEFT ).GT.0 ) THEN
                  CALL CLAMOV_I8( UPLO, M-ILEFT, JJNXTA-JJBEGA+1,
     $                 A( IIA+ILEFT+(JJBEGA-1)*LDA ), LDA,
     $                 B( IIB+ILEFT+(JJBEGB-1)*LDB ), LDB )
                  MYDIST8 = MYDIST8 + NPCOL
                  ILEFT   = MYDIST8 * NBA - ICOFFA
                  JJBEGA  = JJNXTA + 1
                  JJNXTA  = MIN( JJNXTA+NBA, JJENDA )
                  JJBEGB  = JJNXTB + 1
                  JJNXTB  = JJBEGB + JJNXTA - JJBEGA
                  GO TO 30
               END IF
*
            ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
               NQAA   = NQ
               IIAA   = IIA
               JJAA   = JJA
               IIBB   = IIB
               JJBB   = JJB
               IRIGHT = MIN( ILEFT + NBA, M )
               ILEFT  = MIN( MAX( 0_8, ILEFT ), M )
*
   40          CONTINUE
               IF( IIAA.LE.( IIA+M-1 ) ) THEN
                  WIDE = IRIGHT - ILEFT
                  CALL CLAMOV_I8( 'All', ILEFT-IIAA+IIA, NQAA,
     $                 A( IIAA+(JJAA-1)*LDA ), LDA,
     $                 B( IIBB+(JJBB-1)*LDB ), LDB )
                  CALL CLAMOV_I8( UPLO, WIDE, NQAA,
     $                 A( IIA+ILEFT+(JJAA-1)*LDA ), LDA,
     $                 B( IIB+ILEFT+(JJBB-1)*LDB ), LDB )
                  NQAA    = MAX( 0_8, NQAA - WIDE )
                  IIAA    = IIA  + IRIGHT
                  JJAA    = JJAA + WIDE
                  IIBB    = IIB  + IRIGHT
                  JJBB    = JJBB + WIDE
                  MYDIST8 = MYDIST8 + NPCOL
                  ILEFT   = MYDIST8 * NBA - ICOFFA
                  IRIGHT  = MIN( ILEFT + NBA, M )
                  ILEFT   = MIN( ILEFT, M )
                  GO TO 40
               END IF
*
            ELSE
*
               CALL CLAMOV_I8( 'All', M, NQ, A( IIA+(JJA-1)*LDA ),
     $                         LDA, B( IIB+(JJB-1)*LDB ), LDB )
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PCLACP2_I8
*
      END
