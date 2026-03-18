      PROGRAM TEST_PDLAMR1D_I8
      IMPLICIT NONE
*
*  Test for PxLAMR1D_I8 -- INTEGER*8 1D row-vector redistribution.
*  Tests all four type variants (D, S, C, Z) including nonzero CSRC.
*
*  Usage:  mpirun -np 4 ./xdlamr1d_i8
*
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL
      INTEGER            NPROCS, IAM
      INTEGER            NFAIL, NTEST
*
      EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDEXIT, BLACS_EXIT,
     $                   IGAMX2D
*
      CALL BLACS_PINFO( IAM, NPROCS )
*
      IF( NPROCS .LT. 4 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'Need at least 4 processes, have', NPROCS
         CALL BLACS_EXIT( 0 )
         STOP 1
      END IF
*
      NPROW = 2
      NPCOL = 2
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( MYROW .LT. 0 ) THEN
         CALL BLACS_EXIT( 0 )
         STOP
      END IF
*
      NFAIL = 0
      NTEST = 0
*
      IF( IAM .EQ. 0 )
     $   WRITE(*,'(A)') 'test_pdlamr1d_i8: 2x2 grid'
*
*     ==== DOUBLE PRECISION (PDLAMR1D_I8) ====
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PDLAMR1D_I8 (D) ---'
      CALL RUN_D( 100_8, 3_8, 5_8, 0, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'D1:diff-nb', NTEST, NFAIL )
      CALL RUN_D( 200_8, 7_8, 7_8, 0, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'D2:same-nb', NTEST, NFAIL )
      CALL RUN_D( 1_8, 1_8, 1_8, 0, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'D3:N=1    ', NTEST, NFAIL )
      CALL RUN_D( 50_8, 50_8, 1_8, 0, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'D4:big-cyc', NTEST, NFAIL )
      CALL RUN_D( 80_8, 4_8, 6_8, 1, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'D5:csrc=1 ', NTEST, NFAIL )
*
*     ==== REAL (PSLAMR1D_I8) ====
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PSLAMR1D_I8 (S) ---'
      CALL RUN_S( 100_8, 3_8, 5_8, 0, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'S1:diff-nb', NTEST, NFAIL )
      CALL RUN_S( 80_8, 4_8, 6_8, 1, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'S2:csrc=1 ', NTEST, NFAIL )
*
*     ==== COMPLEX (PCLAMR1D_I8) ====
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PCLAMR1D_I8 (C) ---'
      CALL RUN_C( 100_8, 3_8, 5_8, 0, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'C1:diff-nb', NTEST, NFAIL )
      CALL RUN_C( 80_8, 4_8, 6_8, 1, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'C2:csrc=1 ', NTEST, NFAIL )
*
*     ==== DOUBLE COMPLEX (PZLAMR1D_I8) ====
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PZLAMR1D_I8 (Z) ---'
      CALL RUN_Z( 100_8, 3_8, 5_8, 0, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'Z1:diff-nb', NTEST, NFAIL )
      CALL RUN_Z( 80_8, 4_8, 6_8, 1, ICTXT, NPROW, NPCOL,
     $            MYROW, MYCOL, IAM, 'Z2:csrc=1 ', NTEST, NFAIL )
*
*     ==== Summary (reduce NFAIL across all processes) ====
*
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL, 1,
     $              NFAIL, NFAIL, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PxLAMR1D_I8 Test Summary'
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A,I4)') 'Tests run:    ', NTEST
         IF( NFAIL .EQ. 0 ) THEN
            WRITE(*,'(A)') 'TEST PASSED OK'
         ELSE
            WRITE(*,'(A,I6,A)') 'TEST FAILED: ', NFAIL, ' mismatches'
         END IF
         WRITE(*,'(A)') '======================================'
      END IF
*
      CALL BLACS_GRIDEXIT( ICTXT )
      CALL BLACS_EXIT( 0 )
      IF( NFAIL .NE. 0 ) STOP 1
*
      END
*
*     ================================================================
*     RUN_D — DOUBLE PRECISION test via PDLAMR1D_I8
*     ================================================================
*
      SUBROUTINE RUN_D( N8, NBA, NBB, ICSRC, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NBA, NBB
      INTEGER            ICSRC, ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          NQ_A, NQ_B, J8, GIDX
      INTEGER*8          DESCA( 9 ), DESCB( 9 )
      INTEGER            ERRS
      DOUBLE PRECISION, ALLOCATABLE :: A(:), B(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCSET_I8, PDLAMR1D_I8
      INTRINSIC          DBLE, MAX
*
      NQ_A = NUMROC_I8( N8, NBA, MYCOL, ICSRC, NPCOL )
      NQ_B = NUMROC_I8( N8, NBB, MYCOL, ICSRC, NPCOL )
*
      CALL DESCSET_I8( DESCA, 1_8, N8, 1_8, NBA, 0, ICSRC, ICTXT,
     $                 1_8 )
      CALL DESCSET_I8( DESCB, 1_8, N8, 1_8, NBB, 0, ICSRC, ICTXT,
     $                 1_8 )
*
      ALLOCATE( A( MAX( NQ_A, 1_8 ) ) )
      ALLOCATE( B( MAX( NQ_B, 1_8 ) ) )
*
      IF( MYROW .EQ. 0 ) THEN
         DO J8 = 1, NQ_A
            GIDX = INDXL2G_I8( J8, NBA, MYCOL, ICSRC, NPCOL )
            A( J8 ) = DBLE( GIDX )
         END DO
      END IF
*
      DO J8 = 1, MAX( NQ_B, 1_8 )
         B( J8 ) = -999.0D0
      END DO
*
      CALL PDLAMR1D_I8( N8, A, 1_8, 1_8, DESCA,
     $                   B, 1_8, 1_8, DESCB )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, NQ_B
         GIDX = INDXL2G_I8( J8, NBB, MYCOL, ICSRC, NPCOL )
         IF( B( J8 ) .NE. DBLE( GIDX ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ': proc(', MYROW, ',',
     $                     MYCOL, ') got=', B(J8), ' exp=', DBLE(GIDX)
            ERRS = ERRS + 1
         END IF
      END DO
      NFAIL = NFAIL + ERRS
*
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ': PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ': FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
      DEALLOCATE( A, B )
      END
*
*     ================================================================
*     RUN_S — REAL test via PSLAMR1D_I8
*     ================================================================
*
      SUBROUTINE RUN_S( N8, NBA, NBB, ICSRC, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NBA, NBB
      INTEGER            ICSRC, ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          NQ_A, NQ_B, J8, GIDX
      INTEGER*8          DESCA( 9 ), DESCB( 9 )
      INTEGER            ERRS
      REAL, ALLOCATABLE :: A(:), B(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCSET_I8, PSLAMR1D_I8
      INTRINSIC          REAL, MAX
*
      NQ_A = NUMROC_I8( N8, NBA, MYCOL, ICSRC, NPCOL )
      NQ_B = NUMROC_I8( N8, NBB, MYCOL, ICSRC, NPCOL )
*
      CALL DESCSET_I8( DESCA, 1_8, N8, 1_8, NBA, 0, ICSRC, ICTXT,
     $                 1_8 )
      CALL DESCSET_I8( DESCB, 1_8, N8, 1_8, NBB, 0, ICSRC, ICTXT,
     $                 1_8 )
*
      ALLOCATE( A( MAX( NQ_A, 1_8 ) ) )
      ALLOCATE( B( MAX( NQ_B, 1_8 ) ) )
*
      IF( MYROW .EQ. 0 ) THEN
         DO J8 = 1, NQ_A
            GIDX = INDXL2G_I8( J8, NBA, MYCOL, ICSRC, NPCOL )
            A( J8 ) = REAL( GIDX )
         END DO
      END IF
*
      DO J8 = 1, MAX( NQ_B, 1_8 )
         B( J8 ) = -999.0
      END DO
*
      CALL PSLAMR1D_I8( N8, A, 1_8, 1_8, DESCA,
     $                   B, 1_8, 1_8, DESCB )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, NQ_B
         GIDX = INDXL2G_I8( J8, NBB, MYCOL, ICSRC, NPCOL )
         IF( B( J8 ) .NE. REAL( GIDX ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ': proc(', MYROW, ',',
     $                     MYCOL, ') got=', B(J8), ' exp=', REAL(GIDX)
            ERRS = ERRS + 1
         END IF
      END DO
      NFAIL = NFAIL + ERRS
*
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ': PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ': FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
      DEALLOCATE( A, B )
      END
*
*     ================================================================
*     RUN_C — COMPLEX test via PCLAMR1D_I8
*     ================================================================
*
      SUBROUTINE RUN_C( N8, NBA, NBB, ICSRC, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NBA, NBB
      INTEGER            ICSRC, ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          NQ_A, NQ_B, J8, GIDX
      INTEGER*8          DESCA( 9 ), DESCB( 9 )
      INTEGER            ERRS
      COMPLEX, ALLOCATABLE :: A(:), B(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCSET_I8, PCLAMR1D_I8
      INTRINSIC          CMPLX, REAL, MAX
*
      NQ_A = NUMROC_I8( N8, NBA, MYCOL, ICSRC, NPCOL )
      NQ_B = NUMROC_I8( N8, NBB, MYCOL, ICSRC, NPCOL )
*
      CALL DESCSET_I8( DESCA, 1_8, N8, 1_8, NBA, 0, ICSRC, ICTXT,
     $                 1_8 )
      CALL DESCSET_I8( DESCB, 1_8, N8, 1_8, NBB, 0, ICSRC, ICTXT,
     $                 1_8 )
*
      ALLOCATE( A( MAX( NQ_A, 1_8 ) ) )
      ALLOCATE( B( MAX( NQ_B, 1_8 ) ) )
*
      IF( MYROW .EQ. 0 ) THEN
         DO J8 = 1, NQ_A
            GIDX = INDXL2G_I8( J8, NBA, MYCOL, ICSRC, NPCOL )
            A( J8 ) = CMPLX( REAL( GIDX ), -REAL( GIDX ) )
         END DO
      END IF
*
      DO J8 = 1, MAX( NQ_B, 1_8 )
         B( J8 ) = CMPLX( -999.0, -999.0 )
      END DO
*
      CALL PCLAMR1D_I8( N8, A, 1_8, 1_8, DESCA,
     $                   B, 1_8, 1_8, DESCB )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, NQ_B
         GIDX = INDXL2G_I8( J8, NBB, MYCOL, ICSRC, NPCOL )
         IF( REAL( B( J8 ) ) .NE. REAL( GIDX ) .OR.
     $       AIMAG( B( J8 ) ) .NE. -REAL( GIDX ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ': proc(', MYROW, ',',
     $                     MYCOL, ') got=', B(J8)
            ERRS = ERRS + 1
         END IF
      END DO
      NFAIL = NFAIL + ERRS
*
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ': PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ': FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
      DEALLOCATE( A, B )
      END
*
*     ================================================================
*     RUN_Z — DOUBLE COMPLEX test via PZLAMR1D_I8
*     ================================================================
*
      SUBROUTINE RUN_Z( N8, NBA, NBB, ICSRC, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NBA, NBB
      INTEGER            ICSRC, ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          NQ_A, NQ_B, J8, GIDX
      INTEGER*8          DESCA( 9 ), DESCB( 9 )
      INTEGER            ERRS
      DOUBLE COMPLEX, ALLOCATABLE :: A(:), B(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCSET_I8, PZLAMR1D_I8
      INTRINSIC          DBLE, MAX
*
      NQ_A = NUMROC_I8( N8, NBA, MYCOL, ICSRC, NPCOL )
      NQ_B = NUMROC_I8( N8, NBB, MYCOL, ICSRC, NPCOL )
*
      CALL DESCSET_I8( DESCA, 1_8, N8, 1_8, NBA, 0, ICSRC, ICTXT,
     $                 1_8 )
      CALL DESCSET_I8( DESCB, 1_8, N8, 1_8, NBB, 0, ICSRC, ICTXT,
     $                 1_8 )
*
      ALLOCATE( A( MAX( NQ_A, 1_8 ) ) )
      ALLOCATE( B( MAX( NQ_B, 1_8 ) ) )
*
      IF( MYROW .EQ. 0 ) THEN
         DO J8 = 1, NQ_A
            GIDX = INDXL2G_I8( J8, NBA, MYCOL, ICSRC, NPCOL )
            A( J8 ) = DCMPLX( DBLE( GIDX ), -DBLE( GIDX ) )
         END DO
      END IF
*
      DO J8 = 1, MAX( NQ_B, 1_8 )
         B( J8 ) = DCMPLX( -999.0D0, -999.0D0 )
      END DO
*
      CALL PZLAMR1D_I8( N8, A, 1_8, 1_8, DESCA,
     $                   B, 1_8, 1_8, DESCB )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, NQ_B
         GIDX = INDXL2G_I8( J8, NBB, MYCOL, ICSRC, NPCOL )
         IF( DBLE( B( J8 ) ) .NE. DBLE( GIDX ) .OR.
     $       DIMAG( B( J8 ) ) .NE. -DBLE( GIDX ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ': proc(', MYROW, ',',
     $                     MYCOL, ') got=', B(J8)
            ERRS = ERRS + 1
         END IF
      END DO
      NFAIL = NFAIL + ERRS
*
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ': PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ': FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
      DEALLOCATE( A, B )
      END
