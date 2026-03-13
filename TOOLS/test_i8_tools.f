      PROGRAM TEST_I8_TOOLS
      IMPLICIT NONE
*
*  Standalone test for I8 descriptor/tool routines.
*  Runs single-process (no MPI needed) using mock grid values.
*
*  Part A: Cross-validate I8 vs legacy routines on 32-bit-sized cases.
*  Part B: I8-only tests with dimensions > 2^31-1.
*
      INCLUDE 'SL_i8_params.inc'
*     ..
*     .. Local Scalars ..
      INTEGER*8          M8, N8, NB8, LLD8
      INTEGER*8          NR8, GINDX8, LINDX8, GLOB8
      INTEGER*8          LRINDX8, LCINDX8
      INTEGER*8          DESC8(9)
      INTEGER            M4, N4, NB4, NR4, GINDX4, LINDX4, GLOB4
      INTEGER            LR4, LC4
      INTEGER            DESC4(9)
      INTEGER            NPROCS, NPROW, NPCOL, MYROW, MYCOL
      INTEGER            IPROC, ISRCPROC, IRSRC, ICSRC
      INTEGER            RSRC, CSRC, RSRC4, CSRC4
      INTEGER            ROCSRC, ROCSRC4
      INTEGER            PROC, PROC4
      INTEGER            NFAIL, NTEST
      LOGICAL            OK
*     ..
*     .. External Functions (I8) ..
      INTEGER*8          NUMROC_I8, INDXL2G_I8, INDXG2L_I8
      INTEGER            INDXG2P_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8, INDXG2L_I8
      EXTERNAL           INDXG2P_I8
*     ..
*     .. External Functions (legacy 32-bit) ..
      INTEGER            NUMROC, INDXL2G, INDXG2L, INDXG2P
      EXTERNAL           NUMROC, INDXL2G, INDXG2L, INDXG2P
*     ..
*     .. External Subroutines ..
      EXTERNAL           DESCSET_I8, INFOG2L_I8, INFOG1L_I8
      EXTERNAL           DESCSET, INFOG2L, INFOG1L
*     ..
*
      NFAIL = 0
      NTEST = 0
*
      WRITE(*,'(A)') '======================================'
      WRITE(*,'(A)') 'Part A: Cross-validation (32-bit cases)'
      WRITE(*,'(A)') '======================================'
*
*     ============================================================
*     Test A1: NUMROC_I8 matches NUMROC for 32-bit N
*     ============================================================
*
*     Try several (N, NB, NPROCS, ISRCPROC) combinations
*
      CALL CHECK_NUMROC( 1000, 64, 4, 0, NTEST, NFAIL )
      CALL CHECK_NUMROC( 1000, 64, 4, 2, NTEST, NFAIL )
      CALL CHECK_NUMROC( 1, 1, 1, 0, NTEST, NFAIL )
      CALL CHECK_NUMROC( 100000, 32, 6, 3, NTEST, NFAIL )
      CALL CHECK_NUMROC( 999, 100, 3, 1, NTEST, NFAIL )
      CALL CHECK_NUMROC( 2000000000, 128, 4, 0, NTEST, NFAIL )
      CALL CHECK_NUMROC( 2000000000, 128, 4, 3, NTEST, NFAIL )
      CALL CHECK_NUMROC( 1000000, 1, 7, 5, NTEST, NFAIL )
*
*     ============================================================
*     Test A2: INDXL2G_I8 matches INDXL2G for 32-bit values
*     ============================================================
*
      CALL CHECK_INDXL2G( 1, 64, 0, 0, 4, NTEST, NFAIL )
      CALL CHECK_INDXL2G( 65, 64, 0, 0, 4, NTEST, NFAIL )
      CALL CHECK_INDXL2G( 1, 64, 2, 0, 4, NTEST, NFAIL )
      CALL CHECK_INDXL2G( 1, 64, 1, 3, 4, NTEST, NFAIL )
      CALL CHECK_INDXL2G( 100000, 32, 3, 1, 6, NTEST, NFAIL )
      CALL CHECK_INDXL2G( 500000, 128, 0, 0, 4, NTEST, NFAIL )
*
*     ============================================================
*     Test A3: INDXG2L_I8 matches INDXG2L for 32-bit values
*     ============================================================
*
      CALL CHECK_INDXG2L( 1, 64, 0, 0, 4, NTEST, NFAIL )
      CALL CHECK_INDXG2L( 65, 64, 0, 0, 4, NTEST, NFAIL )
      CALL CHECK_INDXG2L( 129, 64, 0, 0, 4, NTEST, NFAIL )
      CALL CHECK_INDXG2L( 1000000, 32, 0, 1, 6, NTEST, NFAIL )
      CALL CHECK_INDXG2L( 2000000000, 128, 0, 0, 4, NTEST, NFAIL )
*
*     ============================================================
*     Test A4: INDXG2P_I8 matches INDXG2P for 32-bit values
*     ============================================================
*
      CALL CHECK_INDXG2P( 1, 64, 0, 0, 4, NTEST, NFAIL )
      CALL CHECK_INDXG2P( 65, 64, 0, 0, 4, NTEST, NFAIL )
      CALL CHECK_INDXG2P( 129, 64, 0, 0, 4, NTEST, NFAIL )
      CALL CHECK_INDXG2P( 1, 64, 0, 2, 4, NTEST, NFAIL )
      CALL CHECK_INDXG2P( 1000000, 32, 0, 1, 6, NTEST, NFAIL )
      CALL CHECK_INDXG2P( 2000000000, 128, 0, 0, 4, NTEST, NFAIL )
*
*     ============================================================
*     Test A5: INFOG2L_I8 matches INFOG2L for 32-bit values
*     ============================================================
*
*     Set up matching descriptors: M=N=10000, MB=NB=64, src=(0,0)
*
      M4 = 10000
      N4 = 10000
      NB4 = 64
      NPROW = 2
      NPCOL = 2
*
      CALL DESCSET( DESC4, M4, N4, NB4, NB4, 0, 0, 42,
     $              NUMROC( M4, NB4, 0, 0, NPROW ) )
      CALL DESCSET_I8( DESC8, INT(M4,8), INT(N4,8),
     $                 INT(NB4,8), INT(NB4,8), 0, 0, 42,
     $                 NUMROC_I8( INT(M4,8), INT(NB4,8), 0, 0,
     $                            NPROW ) )
*
      DO MYROW = 0, NPROW-1
         DO MYCOL = 0, NPCOL-1
*           Test several global indices
            CALL CHECK_INFOG2L( 1, 1, DESC4, DESC8, NPROW, NPCOL,
     $                          MYROW, MYCOL, NTEST, NFAIL )
            CALL CHECK_INFOG2L( 65, 65, DESC4, DESC8, NPROW, NPCOL,
     $                          MYROW, MYCOL, NTEST, NFAIL )
            CALL CHECK_INFOG2L( 1000, 5000, DESC4, DESC8, NPROW,
     $                          NPCOL, MYROW, MYCOL, NTEST, NFAIL )
            CALL CHECK_INFOG2L( 9999, 9999, DESC4, DESC8, NPROW,
     $                          NPCOL, MYROW, MYCOL, NTEST, NFAIL )
         END DO
      END DO
*
*     ============================================================
*     Test A6: INFOG1L_I8 matches INFOG1L for 32-bit values
*     ============================================================
*
      CALL CHECK_INFOG1L( 1, 64, 4, 0, 0, NTEST, NFAIL )
      CALL CHECK_INFOG1L( 65, 64, 4, 0, 0, NTEST, NFAIL )
      CALL CHECK_INFOG1L( 65, 64, 4, 1, 0, NTEST, NFAIL )
      CALL CHECK_INFOG1L( 1000, 32, 6, 3, 1, NTEST, NFAIL )
      CALL CHECK_INFOG1L( 100000, 128, 4, 2, 3, NTEST, NFAIL )
*
      WRITE(*,'(A,I3,A)') 'Part A: ', NTEST, ' cross-validation tests'
*
*     ============================================================
*
      WRITE(*,'(A)') '======================================'
      WRITE(*,'(A)') 'Part B: I8-only tests (> 2^31)'
      WRITE(*,'(A)') '======================================'
*
*     ============================================================
*     Test B1: NUMROC_I8 with N = 3 billion
*     ============================================================
*
      M8 = 3000000000_8
      NB8 = 64_8
      NPROCS = 4
      ISRCPROC = 0
*
*     Sum of NUMROC_I8 across all processes must equal M8
*
      N8 = 0
      DO IPROC = 0, NPROCS-1
         N8 = N8 + NUMROC_I8( M8, NB8, IPROC, ISRCPROC, NPROCS )
      END DO
*
      NTEST = NTEST + 1
      IF( N8 .NE. M8 ) THEN
         WRITE(*,*) 'FAIL B1: NUMROC_I8 sum =', N8, ' expected', M8
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Test B2: DESCSET_I8 with large dimensions
*     ============================================================
*
      M8 = 3000000000_8
      N8 = 2500000000_8
      NB8 = 128_8
      LLD8 = NUMROC_I8( M8, 64_8, 0, 0, 4 )
*
      CALL DESCSET_I8( DESC8, M8, N8, 64_8, NB8, 0, 0, 42, LLD8 )
*
      NTEST = NTEST + 1
      OK = ( DESC8(DTYPE_) .EQ. BLOCK_CYCLIC_2D_I8 ) .AND.
     $     ( DESC8(M_) .EQ. M8 ) .AND.
     $     ( DESC8(N_) .EQ. N8 ) .AND.
     $     ( DESC8(MB_) .EQ. 64 ) .AND.
     $     ( DESC8(NB_) .EQ. NB8 ) .AND.
     $     ( DESC8(RSRC_) .EQ. 0 ) .AND.
     $     ( DESC8(CSRC_) .EQ. 0 ) .AND.
     $     ( DESC8(CTXT_) .EQ. 42 ) .AND.
     $     ( DESC8(LLD_) .EQ. LLD8 )
      IF( .NOT. OK ) THEN
         WRITE(*,*) 'FAIL B2: DESCSET_I8 mismatch'
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Test B3: Round-trip with large local index (> 2^31 global)
*     ============================================================
*
      NB8 = 64_8
      NPROCS = 4
      ISRCPROC = 0
*
      DO IPROC = 0, NPROCS-1
         LINDX8 = 750000001_8
         GLOB8 = INDXL2G_I8(LINDX8, NB8, IPROC, ISRCPROC, NPROCS)
*
*        Global index must exceed INT_MAX
*
         NTEST = NTEST + 1
         IF( GLOB8 .LE. 2147483647_8 ) THEN
            WRITE(*,*) 'FAIL B3a: proc', IPROC,
     $                 ' global =', GLOB8, ' not > INT_MAX'
            NFAIL = NFAIL + 1
         END IF
*
*        Round-trip back to local
*
         NR8 = INDXG2L_I8(GLOB8, NB8, IPROC, ISRCPROC, NPROCS)
         PROC = INDXG2P_I8(GLOB8, NB8, IPROC, ISRCPROC, NPROCS)
*
         NTEST = NTEST + 1
         IF( NR8 .NE. LINDX8 ) THEN
            WRITE(*,*) 'FAIL B3b: proc', IPROC,
     $                 ' round-trip =', NR8, ' expected', LINDX8
            NFAIL = NFAIL + 1
         END IF
*
         NTEST = NTEST + 1
         IF( PROC .NE. IPROC ) THEN
            WRITE(*,*) 'FAIL B3c: proc', IPROC,
     $                 ' INDXG2P =', PROC
            NFAIL = NFAIL + 1
         END IF
*
*        Process coord must be in [0, NPROCS)
*
         NTEST = NTEST + 1
         IF( PROC .LT. 0 .OR. PROC .GE. NPROCS ) THEN
            WRITE(*,*) 'FAIL B3d: proc', IPROC,
     $                 ' INDXG2P out of range:', PROC
            NFAIL = NFAIL + 1
         END IF
      END DO
*
*     ============================================================
*     Test B4: INFOG2L_I8 with large global indices
*     ============================================================
*
      CALL DESCSET_I8( DESC8, 3000000000_8, 3000000000_8,
     $                 64_8, 64_8, 0, 0, 42,
     $                 NUMROC_I8( 3000000000_8, 64_8, 0, 0, 2 ) )
      NPROW = 2
      NPCOL = 2
      MYROW = 0
      MYCOL = 0
*
*     Global (1,1) -> local (1,1) on proc (0,0)
*
      CALL INFOG2L_I8( 1_8, 1_8, DESC8, NPROW, NPCOL, MYROW, MYCOL,
     $                 LRINDX8, LCINDX8, RSRC, CSRC )
      NTEST = NTEST + 1
      IF( LRINDX8.NE.1 .OR. LCINDX8.NE.1 .OR.
     $    RSRC.NE.0 .OR. CSRC.NE.0 ) THEN
         WRITE(*,*) 'FAIL B4a: INFOG2L_I8(1,1) =',
     $              LRINDX8, LCINDX8, RSRC, CSRC
         NFAIL = NFAIL + 1
      END IF
*
*     Large row index: 2000000001, row block 31250000 -> even -> proc 0
*
      CALL INFOG2L_I8( 2000000001_8, 1_8, DESC8, NPROW, NPCOL,
     $                 MYROW, MYCOL, LRINDX8, LCINDX8, RSRC, CSRC )
      NTEST = NTEST + 1
      IF( RSRC.NE.0 .OR. LRINDX8.LE.0 ) THEN
         WRITE(*,*) 'FAIL B4b: INFOG2L_I8(2e9,1) RSRC=', RSRC,
     $              ' LRINDX=', LRINDX8
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Test B5: INFOG1L_I8 with large index
*     ============================================================
*
      NB8 = 64_8
      NPROCS = 4
      GINDX8 = 3000000000_8
      CALL INFOG1L_I8( GINDX8, NB8, NPROCS, 0, 0, LINDX8, ROCSRC )
*
      NTEST = NTEST + 1
      IF( ROCSRC.LT.0 .OR. ROCSRC.GE.NPROCS ) THEN
         WRITE(*,*) 'FAIL B5a: ROCSRC out of range:', ROCSRC
         NFAIL = NFAIL + 1
      END IF
      NTEST = NTEST + 1
      IF( LINDX8 .LE. 0 ) THEN
         WRITE(*,*) 'FAIL B5b: LINDX =', LINDX8, ' expected > 0'
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Test B6: Non-zero ISRCPROC with large values
*     ============================================================
*
      NB8 = 64_8
      NPROCS = 4
      ISRCPROC = 2
      IPROC = 3
      LINDX8 = 750000001_8
*
      GLOB8 = INDXL2G_I8(LINDX8, NB8, IPROC, ISRCPROC, NPROCS)
      NR8 = INDXG2L_I8(GLOB8, NB8, IPROC, ISRCPROC, NPROCS)
      PROC = INDXG2P_I8(GLOB8, NB8, IPROC, ISRCPROC, NPROCS)
*
      NTEST = NTEST + 1
      IF( NR8 .NE. LINDX8 ) THEN
         WRITE(*,*) 'FAIL B6a: round-trip =', NR8,
     $              ' expected', LINDX8
         NFAIL = NFAIL + 1
      END IF
      NTEST = NTEST + 1
      IF( PROC .NE. IPROC ) THEN
         WRITE(*,*) 'FAIL B6b: INDXG2P =', PROC,
     $              ' expected', IPROC
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Summary
*     ============================================================
*
      WRITE(*,'(A)')    '======================================'
      WRITE(*,'(A)')    'I8 Tools Test Summary'
      WRITE(*,'(A)')    '======================================'
      WRITE(*,'(A,I4)') 'Tests run:    ', NTEST
      WRITE(*,'(A,I4)') 'Tests passed: ', NTEST - NFAIL
      WRITE(*,'(A,I4)') 'Tests failed: ', NFAIL
      IF( NFAIL .EQ. 0 ) THEN
         WRITE(*,'(A)') 'ALL TESTS PASSED'
      ELSE
         WRITE(*,'(A)') '*** FAILURES DETECTED ***'
      END IF
      WRITE(*,'(A)')    '======================================'
*
      IF( NFAIL .NE. 0 ) STOP 1
*
      END
*
*     ================================================================
*     Cross-validation helper subroutines
*     ================================================================
*
      SUBROUTINE CHECK_NUMROC( N, NB, NPROCS, ISRCPROC,
     $                         NTEST, NFAIL )
      IMPLICIT NONE
      INTEGER            N, NB, NPROCS, ISRCPROC, NTEST, NFAIL
*     ..
      INTEGER            NUMROC, IPROC, NR4
      INTEGER*8          NUMROC_I8, NR8
      EXTERNAL           NUMROC, NUMROC_I8
*
      DO IPROC = 0, NPROCS-1
         NR4 = NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
         NR8 = NUMROC_I8( INT(N,8), INT(NB,8), IPROC, ISRCPROC,
     $                     NPROCS )
         NTEST = NTEST + 1
         IF( NR8 .NE. NR4 ) THEN
            WRITE(*,*) 'FAIL NUMROC: N=', N, ' NB=', NB,
     $                 ' proc=', IPROC, ' src=', ISRCPROC,
     $                 ' legacy=', NR4, ' I8=', NR8
            NFAIL = NFAIL + 1
         END IF
      END DO
*
      END
*
      SUBROUTINE CHECK_INDXL2G( INDXLOC, NB, IPROC, ISRCPROC,
     $                          NPROCS, NTEST, NFAIL )
      IMPLICIT NONE
      INTEGER            INDXLOC, NB, IPROC, ISRCPROC, NPROCS
      INTEGER            NTEST, NFAIL
*     ..
      INTEGER            INDXL2G, G4
      INTEGER*8          INDXL2G_I8, G8
      EXTERNAL           INDXL2G, INDXL2G_I8
*
      G4 = INDXL2G( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )
      G8 = INDXL2G_I8( INT(INDXLOC,8), INT(NB,8), IPROC, ISRCPROC,
     $                  NPROCS )
      NTEST = NTEST + 1
      IF( G8 .NE. G4 ) THEN
         WRITE(*,*) 'FAIL INDXL2G: loc=', INDXLOC, ' NB=', NB,
     $              ' proc=', IPROC, ' legacy=', G4, ' I8=', G8
         NFAIL = NFAIL + 1
      END IF
*
      END
*
      SUBROUTINE CHECK_INDXG2L( INDXGLOB, NB, IPROC, ISRCPROC,
     $                          NPROCS, NTEST, NFAIL )
      IMPLICIT NONE
      INTEGER            INDXGLOB, NB, IPROC, ISRCPROC, NPROCS
      INTEGER            NTEST, NFAIL
*     ..
      INTEGER            INDXG2L, L4
      INTEGER*8          INDXG2L_I8, L8
      EXTERNAL           INDXG2L, INDXG2L_I8
*
      L4 = INDXG2L( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
      L8 = INDXG2L_I8( INT(INDXGLOB,8), INT(NB,8), IPROC, ISRCPROC,
     $                  NPROCS )
      NTEST = NTEST + 1
      IF( L8 .NE. L4 ) THEN
         WRITE(*,*) 'FAIL INDXG2L: glob=', INDXGLOB, ' NB=', NB,
     $              ' legacy=', L4, ' I8=', L8
         NFAIL = NFAIL + 1
      END IF
*
      END
*
      SUBROUTINE CHECK_INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC,
     $                          NPROCS, NTEST, NFAIL )
      IMPLICIT NONE
      INTEGER            INDXGLOB, NB, IPROC, ISRCPROC, NPROCS
      INTEGER            NTEST, NFAIL
*     ..
      INTEGER            INDXG2P, P4
      INTEGER            INDXG2P_I8, P8
      EXTERNAL           INDXG2P, INDXG2P_I8
*
      P4 = INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
      P8 = INDXG2P_I8( INT(INDXGLOB,8), INT(NB,8), IPROC, ISRCPROC,
     $                  NPROCS )
      NTEST = NTEST + 1
      IF( P8 .NE. P4 ) THEN
         WRITE(*,*) 'FAIL INDXG2P: glob=', INDXGLOB, ' NB=', NB,
     $              ' legacy=', P4, ' I8=', P8
         NFAIL = NFAIL + 1
      END IF
*
      END
*
      SUBROUTINE CHECK_INFOG2L( GR, GC, DESC4, DESC8, NPROW, NPCOL,
     $                          MYROW, MYCOL, NTEST, NFAIL )
      IMPLICIT NONE
      INTEGER            GR, GC, NPROW, NPCOL, MYROW, MYCOL
      INTEGER            DESC4( * ), NTEST, NFAIL
      INTEGER*8          DESC8( * )
*     ..
      INTEGER            LR4, LC4, RSRC4, CSRC4
      INTEGER*8          LR8, LC8
      INTEGER            RSRC8, CSRC8
      EXTERNAL           INFOG2L, INFOG2L_I8
*
      CALL INFOG2L( GR, GC, DESC4, NPROW, NPCOL, MYROW, MYCOL,
     $              LR4, LC4, RSRC4, CSRC4 )
      CALL INFOG2L_I8( INT(GR,8), INT(GC,8), DESC8, NPROW, NPCOL,
     $                 MYROW, MYCOL, LR8, LC8, RSRC8, CSRC8 )
*
      NTEST = NTEST + 1
      IF( LR8.NE.LR4 .OR. LC8.NE.LC4 .OR.
     $    RSRC8.NE.RSRC4 .OR. CSRC8.NE.CSRC4 ) THEN
         WRITE(*,*) 'FAIL INFOG2L: GR=', GR, ' GC=', GC,
     $              ' myrow=', MYROW, ' mycol=', MYCOL
         WRITE(*,*) '  legacy: LR=', LR4, ' LC=', LC4,
     $              ' RSRC=', RSRC4, ' CSRC=', CSRC4
         WRITE(*,*) '  I8:     LR=', LR8, ' LC=', LC8,
     $              ' RSRC=', RSRC8, ' CSRC=', CSRC8
         NFAIL = NFAIL + 1
      END IF
*
      END
*
      SUBROUTINE CHECK_INFOG1L( GINDX, NB, NPROCS, MYROC, ISRCPROC,
     $                          NTEST, NFAIL )
      IMPLICIT NONE
      INTEGER            GINDX, NB, NPROCS, MYROC, ISRCPROC
      INTEGER            NTEST, NFAIL
*     ..
      INTEGER            LINDX4, ROCSRC4
      INTEGER*8          LINDX8
      INTEGER            ROCSRC8
      EXTERNAL           INFOG1L, INFOG1L_I8
*
      CALL INFOG1L( GINDX, NB, NPROCS, MYROC, ISRCPROC,
     $              LINDX4, ROCSRC4 )
      CALL INFOG1L_I8( INT(GINDX,8), INT(NB,8), NPROCS, MYROC,
     $                 ISRCPROC, LINDX8, ROCSRC8 )
*
      NTEST = NTEST + 1
      IF( LINDX8.NE.LINDX4 .OR. ROCSRC8.NE.ROCSRC4 ) THEN
         WRITE(*,*) 'FAIL INFOG1L: GINDX=', GINDX, ' NB=', NB,
     $              ' myroc=', MYROC
         WRITE(*,*) '  legacy: LINDX=', LINDX4, ' ROCSRC=', ROCSRC4
         WRITE(*,*) '  I8:     LINDX=', LINDX8, ' ROCSRC=', ROCSRC8
         NFAIL = NFAIL + 1
      END IF
*
      END
