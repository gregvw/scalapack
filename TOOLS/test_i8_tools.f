      PROGRAM TEST_I8_TOOLS
      IMPLICIT NONE
*
*  Standalone test for I8 descriptor/tool routines.
*  Runs single-process (no MPI needed) using mock grid values.
*  Tests with matrix dimensions > INT_MAX (2^31-1 = 2147483647).
*
*     .. Parameters ..
      INTEGER*8          BLOCK_CYCLIC_2D_I8
      PARAMETER          ( BLOCK_CYCLIC_2D_I8 = 501 )
      INTEGER            DTYPE_, CTXT_, M_, N_, MB_, NB_,
     $                   RSRC_, CSRC_, LLD_
      PARAMETER          ( DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                     MB_ = 5, NB_ = 6, RSRC_ = 7, CSRC_ = 8,
     $                     LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER*8          M, N, MB, NB, LLD
      INTEGER*8          NR, GINDX, LINDX, GLOB_RT
      INTEGER*8          LRINDX, LCINDX
      INTEGER*8          DESC(9)
      INTEGER            NPROCS, NPROW, NPCOL, MYROW, MYCOL
      INTEGER            IPROC, ISRCPROC, IRSRC, ICSRC
      INTEGER            RSRC, CSRC, ROCSRC, PROC, INFO
      INTEGER            NFAIL, NTEST
      LOGICAL            OK
*     ..
*     .. External Functions ..
      INTEGER*8          NUMROC_I8, INDXL2G_I8, INDXG2L_I8
      INTEGER            INDXG2P_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8, INDXG2L_I8
      EXTERNAL           INDXG2P_I8
*     ..
*     .. External Subroutines ..
      EXTERNAL           DESCSET_I8, INFOG2L_I8, INFOG1L_I8
*     ..
*
      NFAIL = 0
      NTEST = 0
*
*     ============================================================
*     Test 1: NUMROC_I8 with large N
*     ============================================================
*
*     N = 3 billion, NB = 64, 4 processes, source = 0
*
      M = 3000000000_8
      NB = 64_8
      NPROCS = 4
      ISRCPROC = 0
*
*     Sum of NUMROC_I8 across all processes should equal N
*
      N = 0
      DO IPROC = 0, NPROCS-1
         N = N + NUMROC_I8( M, NB, IPROC, ISRCPROC, NPROCS )
      END DO
*
      NTEST = NTEST + 1
      IF( N .NE. M ) THEN
         WRITE(*,*) 'FAIL Test 1: NUMROC_I8 sum =', N, ' expected', M
         NFAIL = NFAIL + 1
      END IF
*
*     Check specific process: proc 0 should get ceil(N/NB/NPROCS)*NB
*     + possible partial block
*
      NR = NUMROC_I8( M, NB, 0, 0, NPROCS )
      NTEST = NTEST + 1
      IF( NR .LE. 0 .OR. NR .GT. M ) THEN
         WRITE(*,*) 'FAIL Test 1b: NUMROC_I8(proc=0) =', NR
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Test 2: DESCSET_I8 sets descriptor correctly
*     ============================================================
*
      M = 3000000000_8
      N = 2500000000_8
      MB = 64_8
      NB = 128_8
      IRSRC = 0
      ICSRC = 0
*
*     Use a mock context (just a number, not actually used)
*
      LLD = NUMROC_I8( M, MB, 0, 0, 4 )
*
      CALL DESCSET_I8( DESC, M, N, MB, NB, IRSRC, ICSRC, 42, LLD )
*
      NTEST = NTEST + 1
      OK = ( DESC(DTYPE_) .EQ. BLOCK_CYCLIC_2D_I8 ) .AND.
     $     ( DESC(M_) .EQ. M ) .AND.
     $     ( DESC(N_) .EQ. N ) .AND.
     $     ( DESC(MB_) .EQ. MB ) .AND.
     $     ( DESC(NB_) .EQ. NB ) .AND.
     $     ( DESC(RSRC_) .EQ. IRSRC ) .AND.
     $     ( DESC(CSRC_) .EQ. ICSRC ) .AND.
     $     ( DESC(CTXT_) .EQ. 42 ) .AND.
     $     ( DESC(LLD_) .EQ. LLD )
      IF( .NOT. OK ) THEN
         WRITE(*,*) 'FAIL Test 2: DESCSET_I8 descriptor mismatch'
         WRITE(*,*) '  DESC =', DESC
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Test 3: INDXL2G_I8 / INDXG2L_I8 / INDXG2P_I8 round-trip
*     ============================================================
*
*     For proc 0 on a 4-process grid with NB=64, source=0:
*     Local index 1 -> global index 1
*     Local index 65 -> global index 257 (second block on proc 0)
*
      NB = 64_8
      NPROCS = 4
      ISRCPROC = 0
      IPROC = 0
*
*     Test local index 1
*
      GLOB_RT = INDXL2G_I8( 1_8, NB, IPROC, ISRCPROC, NPROCS )
      NTEST = NTEST + 1
      IF( GLOB_RT .NE. 1 ) THEN
         WRITE(*,*) 'FAIL Test 3a: INDXL2G_I8(1) =', GLOB_RT,
     $              ' expected 1'
         NFAIL = NFAIL + 1
      END IF
*
*     Test local index 65 (first element of second local block)
*
      GLOB_RT = INDXL2G_I8( 65_8, NB, IPROC, ISRCPROC, NPROCS )
      NTEST = NTEST + 1
      IF( GLOB_RT .NE. 257 ) THEN
         WRITE(*,*) 'FAIL Test 3b: INDXL2G_I8(65) =', GLOB_RT,
     $              ' expected 257'
         NFAIL = NFAIL + 1
      END IF
*
*     Round-trip: local -> global -> local, global -> process
*     Use a large local index that forces INTEGER*8 arithmetic
*
      LINDX = 750000001_8
      GLOB_RT = INDXL2G_I8( LINDX, NB, IPROC, ISRCPROC, NPROCS )
*
*     The global index should be > INT_MAX for large local index
*
      NTEST = NTEST + 1
      IF( GLOB_RT .LE. 0 ) THEN
         WRITE(*,*) 'FAIL Test 3c: INDXL2G_I8 returned', GLOB_RT
         NFAIL = NFAIL + 1
      END IF
*
*     Convert back: global -> local
*
      NR = INDXG2L_I8( GLOB_RT, NB, IPROC, ISRCPROC, NPROCS )
      NTEST = NTEST + 1
      IF( NR .NE. LINDX ) THEN
         WRITE(*,*) 'FAIL Test 3d: round-trip local =', NR,
     $              ' expected', LINDX
         NFAIL = NFAIL + 1
      END IF
*
*     Global -> process coordinate
*
      PROC = INDXG2P_I8( GLOB_RT, NB, IPROC, ISRCPROC, NPROCS )
      NTEST = NTEST + 1
      IF( PROC .NE. IPROC ) THEN
         WRITE(*,*) 'FAIL Test 3e: INDXG2P_I8 =', PROC,
     $              ' expected', IPROC
         NFAIL = NFAIL + 1
      END IF
*
*     Verify INDXG2P_I8 returns value in [0, NPROCS)
*
      NTEST = NTEST + 1
      IF( PROC .LT. 0 .OR. PROC .GE. NPROCS ) THEN
         WRITE(*,*) 'FAIL Test 3f: INDXG2P_I8 out of range:', PROC
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Test 4: INDXG2P_I8 for different processes
*     ============================================================
*
*     Global index 65 with NB=64, src=0, nprocs=4 -> proc 1
*
      PROC = INDXG2P_I8( 65_8, 64_8, 0, 0, 4 )
      NTEST = NTEST + 1
      IF( PROC .NE. 1 ) THEN
         WRITE(*,*) 'FAIL Test 4a: INDXG2P_I8(65) =', PROC,
     $              ' expected 1'
         NFAIL = NFAIL + 1
      END IF
*
*     Global index 129 with NB=64, src=0, nprocs=4 -> proc 2
*
      PROC = INDXG2P_I8( 129_8, 64_8, 0, 0, 4 )
      NTEST = NTEST + 1
      IF( PROC .NE. 2 ) THEN
         WRITE(*,*) 'FAIL Test 4b: INDXG2P_I8(129) =', PROC,
     $              ' expected 2'
         NFAIL = NFAIL + 1
      END IF
*
*     Large global index: 3000000000 with NB=64, src=0, nprocs=4
*     Block = (3000000000-1)/64 = 46874999, mod 4 = 3 -> proc 3
*
      PROC = INDXG2P_I8( 3000000000_8, 64_8, 0, 0, 4 )
      NTEST = NTEST + 1
      IF( PROC .NE. 3 ) THEN
         WRITE(*,*) 'FAIL Test 4c: INDXG2P_I8(3e9) =', PROC,
     $              ' expected 3'
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Test 5: INFOG2L_I8 with large indices
*     ============================================================
*
*     Set up descriptor: M=3e9, N=3e9, MB=NB=64, src=(0,0)
*
      CALL DESCSET_I8( DESC, 3000000000_8, 3000000000_8,
     $                 64_8, 64_8, 0, 0, 42,
     $                 NUMROC_I8( 3000000000_8, 64_8, 0, 0, 2 ) )
*
      NPROW = 2
      NPCOL = 2
      MYROW = 0
      MYCOL = 0
*
*     Global index (1,1) should map to local (1,1) on proc (0,0)
*
      CALL INFOG2L_I8( 1_8, 1_8, DESC, NPROW, NPCOL, MYROW, MYCOL,
     $                 LRINDX, LCINDX, RSRC, CSRC )
*
      NTEST = NTEST + 1
      IF( LRINDX.NE.1 .OR. LCINDX.NE.1 .OR.
     $    RSRC.NE.0 .OR. CSRC.NE.0 ) THEN
         WRITE(*,*) 'FAIL Test 5a: INFOG2L_I8(1,1) =',
     $              LRINDX, LCINDX, RSRC, CSRC
         NFAIL = NFAIL + 1
      END IF
*
*     Global index (65,65) -> block (1,1) -> proc (1,1)
*     On proc (0,0) this should still compute a valid local offset
*
      CALL INFOG2L_I8( 65_8, 65_8, DESC, NPROW, NPCOL, MYROW, MYCOL,
     $                 LRINDX, LCINDX, RSRC, CSRC )
*
      NTEST = NTEST + 1
      IF( RSRC.NE.1 .OR. CSRC.NE.1 ) THEN
         WRITE(*,*) 'FAIL Test 5b: INFOG2L_I8(65,65) RSRC/CSRC =',
     $              RSRC, CSRC, ' expected 1,1'
         NFAIL = NFAIL + 1
      END IF
*
*     Large index: global (2000000001, 1) on (myrow=0, mycol=0)
*     Block = (2000000001-1)/64 = 31250000, proc = mod(31250000,2) = 0
*     So RSRC should be 0 for this proc
*
      CALL INFOG2L_I8( 2000000001_8, 1_8, DESC, NPROW, NPCOL,
     $                 MYROW, MYCOL, LRINDX, LCINDX, RSRC, CSRC )
*
      NTEST = NTEST + 1
      IF( RSRC.NE.0 ) THEN
         WRITE(*,*) 'FAIL Test 5c: INFOG2L_I8 RSRC =', RSRC,
     $              ' expected 0'
         NFAIL = NFAIL + 1
      END IF
*
      NTEST = NTEST + 1
      IF( LRINDX .LE. 0 ) THEN
         WRITE(*,*) 'FAIL Test 5d: INFOG2L_I8 LRINDX =', LRINDX,
     $              ' expected > 0'
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Test 6: INFOG1L_I8 with large index
*     ============================================================
*
      NB = 64_8
      NPROCS = 4
*
*     Global index 1, proc 0, src 0 -> local 1, rocsrc 0
*
      CALL INFOG1L_I8( 1_8, NB, NPROCS, 0, 0, LINDX, ROCSRC )
*
      NTEST = NTEST + 1
      IF( LINDX.NE.1 .OR. ROCSRC.NE.0 ) THEN
         WRITE(*,*) 'FAIL Test 6a: INFOG1L_I8(1) =',
     $              LINDX, ROCSRC
         NFAIL = NFAIL + 1
      END IF
*
*     Large global index
*
      GINDX = 3000000000_8
      CALL INFOG1L_I8( GINDX, NB, NPROCS, 0, 0, LINDX, ROCSRC )
*
      NTEST = NTEST + 1
      IF( ROCSRC.LT.0 .OR. ROCSRC.GE.NPROCS ) THEN
         WRITE(*,*) 'FAIL Test 6b: INFOG1L_I8 ROCSRC =', ROCSRC,
     $              ' out of range'
         NFAIL = NFAIL + 1
      END IF
*
      NTEST = NTEST + 1
      IF( LINDX .LE. 0 ) THEN
         WRITE(*,*) 'FAIL Test 6c: INFOG1L_I8 LINDX =', LINDX,
     $              ' expected > 0'
         NFAIL = NFAIL + 1
      END IF
*
*     ============================================================
*     Test 7: Comprehensive round-trip for all procs
*     ============================================================
*
*     For each proc in [0,NPROCS), pick a local index, go to global,
*     come back, verify identity.
*
      NB = 128_8
      NPROCS = 4
      ISRCPROC = 0
*
      DO IPROC = 0, NPROCS-1
         LINDX = 500000000_8
         GLOB_RT = INDXL2G_I8(LINDX, NB, IPROC, ISRCPROC, NPROCS)
         NR = INDXG2L_I8(GLOB_RT, NB, IPROC, ISRCPROC, NPROCS)
         PROC = INDXG2P_I8(GLOB_RT, NB, IPROC, ISRCPROC, NPROCS)
*
         NTEST = NTEST + 1
         IF( NR .NE. LINDX ) THEN
            WRITE(*,*) 'FAIL Test 7a: proc', IPROC,
     $                 ' round-trip local =', NR, ' expected', LINDX
            NFAIL = NFAIL + 1
         END IF
*
         NTEST = NTEST + 1
         IF( PROC .NE. IPROC ) THEN
            WRITE(*,*) 'FAIL Test 7b: proc', IPROC,
     $                 ' INDXG2P =', PROC
            NFAIL = NFAIL + 1
         END IF
      END DO
*
*     ============================================================
*     Test 8: Non-zero ISRCPROC
*     ============================================================
*
      NB = 64_8
      NPROCS = 4
      ISRCPROC = 2
      IPROC = 3
      LINDX = 750000001_8
*
      GLOB_RT = INDXL2G_I8(LINDX, NB, IPROC, ISRCPROC, NPROCS)
      NR = INDXG2L_I8(GLOB_RT, NB, IPROC, ISRCPROC, NPROCS)
      PROC = INDXG2P_I8(GLOB_RT, NB, IPROC, ISRCPROC, NPROCS)
*
      NTEST = NTEST + 1
      IF( NR .NE. LINDX ) THEN
         WRITE(*,*) 'FAIL Test 8a: non-zero src round-trip =', NR,
     $              ' expected', LINDX
         NFAIL = NFAIL + 1
      END IF
*
      NTEST = NTEST + 1
      IF( PROC .NE. IPROC ) THEN
         WRITE(*,*) 'FAIL Test 8b: non-zero src INDXG2P =', PROC,
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
      WRITE(*,'(A,I3)') 'Tests run:    ', NTEST
      WRITE(*,'(A,I3)') 'Tests passed: ', NTEST - NFAIL
      WRITE(*,'(A,I3)') 'Tests failed: ', NFAIL
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
