      PROGRAM TEST_I8_TOOLS_MPI
      IMPLICIT NONE
*
*  MPI-backed I8 tools tests — exercises paths that require a real
*  BLACS grid: DESC_CONVERT_I8 2D->1D, PCHK1MAT_I8, GLOBCHK_I8.
*
*  Usage:  mpirun -np 4 ./xi8tools_mpi
*
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL
      INTEGER            NPROCS, IAM, NFAIL, NTEST, NFAIL_G
*
      EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDEXIT, BLACS_EXIT,
     $                   IGAMX2D
*
      CALL BLACS_PINFO( IAM, NPROCS )
      IF( NPROCS .LT. 4 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'Need at least 4 processes'
         CALL BLACS_EXIT( 0 )
         STOP 1
      END IF
*
      NFAIL = 0
      NTEST = 0
*
*     ============================================================
*     Test 1: DESC_CONVERT_I8 — 2D to 1D horizontal (1xP grid)
*     ============================================================
*
      IF( IAM .EQ. 0 )
     $   WRITE(*,'(A)') '--- DESC_CONVERT_I8: 2D->1D_H ---'
      CALL TEST_DESC_CONV_2D_1DH( IAM, NPROCS, NTEST, NFAIL )
*
*     ============================================================
*     Test 2: DESC_CONVERT_I8 — 2D to 1D vertical (Px1 grid)
*     ============================================================
*
      IF( IAM .EQ. 0 )
     $   WRITE(*,'(A)') '--- DESC_CONVERT_I8: 2D->1D_V ---'
      CALL TEST_DESC_CONV_2D_1DV( IAM, NPROCS, NTEST, NFAIL )
*
*     ============================================================
*     Test 3: PCHK1MAT_I8 — consistent descriptor across 2x2 grid
*     ============================================================
*
      NPROW = 2
      NPCOL = 2
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( MYROW .GE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '--- PCHK1MAT_I8: consistent ---'
         CALL TEST_PCHK1MAT_OK( ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $                           IAM, NTEST, NFAIL )
*
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '--- PCHK1MAT_I8: inconsistent ---'
         CALL TEST_PCHK1MAT_BAD( ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $                            IAM, NTEST, NFAIL )
*
         CALL BLACS_GRIDEXIT( ICTXT )
      END IF
*
*     Summary
*
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'R', 1, NPROCS )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      IF( MYROW .GE. 0 ) THEN
         NFAIL_G = NFAIL
         CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $                 NFAIL_G, NFAIL_G, -1, -1, -1 )
         CALL BLACS_GRIDEXIT( ICTXT )
      END IF
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'I8 Tools MPI Test Summary'
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A,I4)') 'Tests run:    ', NTEST
         IF( NFAIL_G .EQ. 0 ) THEN
            WRITE(*,'(A)') 'TEST PASSED OK'
         ELSE
            WRITE(*,'(A,I6,A)') 'TEST FAILED: ',NFAIL_G,' errors'
         END IF
         WRITE(*,'(A)') '======================================'
      END IF
*
      CALL BLACS_EXIT( 0 )
      IF( NFAIL_G .NE. 0 ) STOP 1
*
      END
*
*     ================================================================
*     TEST_DESC_CONV_2D_1DH — 2D descriptor on 1xP grid -> 1D_H
*     ================================================================
*
      SUBROUTINE TEST_DESC_CONV_2D_1DH( IAM, NPROCS, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER            IAM, NPROCS, NTEST, NFAIL
*
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, INFO
      INTEGER*8          DESC2D( 9 ), DESC1D( 9 )
      INTEGER*8          N8, NB8, LLD8
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
      EXTERNAL           BLACS_GET, BLACS_GRIDINIT, BLACS_GRIDINFO,
     $                   BLACS_GRIDEXIT, DESCINIT_I8, DESC_CONVERT_I8
*
*     1xP grid (required for horizontal conversion)
      NPROW = 1
      NPCOL = NPROCS
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( MYROW .LT. 0 ) RETURN
*
      N8   = 100_8
      NB8  = 8_8
      LLD8 = MAX( NUMROC_I8( N8, NB8, MYROW, 0, NPROW ), 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESC2D, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
*
*     Convert 2D -> 1D horizontal
      DESC1D( 1 ) = 501
      CALL DESC_CONVERT_I8( DESC2D, DESC1D, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A,I4)') '  2D->1DH: FAILED INFO=', INFO
      ELSE IF( DESC1D( 3 ) .NE. N8 .OR. DESC1D( 4 ) .NE. NB8 ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  2D->1DH: FAILED wrong N or NB'
      ELSE
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  2D->1DH: PASSED'
      END IF
*
      CALL BLACS_GRIDEXIT( ICTXT )
      END
*
*     ================================================================
*     TEST_DESC_CONV_2D_1DV — 2D descriptor on Px1 grid -> 1D_V
*     ================================================================
*
      SUBROUTINE TEST_DESC_CONV_2D_1DV( IAM, NPROCS, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER            IAM, NPROCS, NTEST, NFAIL
*
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, INFO
      INTEGER*8          DESC2D( 9 ), DESC1D( 9 )
      INTEGER*8          N8, NB8, LLD8
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
      EXTERNAL           BLACS_GET, BLACS_GRIDINIT, BLACS_GRIDINFO,
     $                   BLACS_GRIDEXIT, DESCINIT_I8, DESC_CONVERT_I8
*
*     Px1 grid (required for vertical conversion)
      NPROW = NPROCS
      NPCOL = 1
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( MYROW .LT. 0 ) RETURN
*
      N8   = 100_8
      NB8  = 8_8
      LLD8 = MAX( NUMROC_I8( N8, NB8, MYROW, 0, NPROW ), 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESC2D, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
*
*     Convert 2D -> 1D vertical
      DESC1D( 1 ) = 502
      CALL DESC_CONVERT_I8( DESC2D, DESC1D, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A,I4)') '  2D->1DV: FAILED INFO=', INFO
      ELSE IF( DESC1D( 3 ) .NE. N8 .OR. DESC1D( 4 ) .NE. NB8 ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  2D->1DV: FAILED wrong M or MB'
      ELSE
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  2D->1DV: PASSED'
      END IF
*
      CALL BLACS_GRIDEXIT( ICTXT )
      END
*
*     ================================================================
*     TEST_PCHK1MAT_OK — all processes agree on descriptor
*     ================================================================
*
      SUBROUTINE TEST_PCHK1MAT_OK( ICTXT, NPROW, NPCOL, MYROW,
     $                               MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          N8, NB8, LLD8, IA8, JA8
      INTEGER*8          DESCA8( 9 ), IDUM1( 1 )
      INTEGER            IDUM2( 1 ), INFO
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
      EXTERNAL           DESCINIT_I8, PCHK1MAT_I8
*
      N8   = 50_8
      NB8  = 4_8
      IA8  = 1_8
      JA8  = 1_8
      LLD8 = MAX( NUMROC_I8( N8, NB8, MYROW, 0, NPROW ), 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A,I4)') '  PCHK1MAT consistent: DESCINIT FAILED'
         NFAIL = NFAIL + 1
         NTEST = NTEST + 1
         RETURN
      END IF
*
*     PCHK1MAT_I8 with no extra values — should return INFO=0
      INFO = 0
      CALL PCHK1MAT_I8( N8, 1, N8, 2, IA8, JA8, DESCA8, 5,
     $                   0, IDUM1, IDUM2, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A,I4)') '  PCHK1MAT consistent: FAILED INFO=',
     $                         INFO
      ELSE
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  PCHK1MAT consistent: PASSED'
      END IF
      END
*
*     ================================================================
*     TEST_PCHK1MAT_BAD — one process has wrong N, should detect
*     ================================================================
*
      SUBROUTINE TEST_PCHK1MAT_BAD( ICTXT, NPROW, NPCOL, MYROW,
     $                                MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          N8, NB8, LLD8, IA8, JA8
      INTEGER*8          DESCA8( 9 ), IDUM1( 1 )
      INTEGER            IDUM2( 1 ), INFO
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
      EXTERNAL           DESCINIT_I8, PCHK1MAT_I8
*
      N8   = 50_8
      NB8  = 4_8
      IA8  = 1_8
      JA8  = 1_8
      LLD8 = MAX( NUMROC_I8( N8, NB8, MYROW, 0, NPROW ), 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
*
*     Corrupt N on rank 1 only — PCHK1MAT_I8 should detect mismatch
      IF( MYROW .EQ. 1 .AND. MYCOL .EQ. 0 ) THEN
         N8 = 999_8
      END IF
*
      INFO = 0
      CALL PCHK1MAT_I8( N8, 1, N8, 2, IA8, JA8, DESCA8, 5,
     $                   0, IDUM1, IDUM2, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .EQ. 0 ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  PCHK1MAT inconsistent: FAILED (missed)'
      ELSE
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A,I4)') '  PCHK1MAT inconsistent: PASSED INFO=',
     $                         INFO
      END IF
      END
