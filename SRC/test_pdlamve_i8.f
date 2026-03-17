      PROGRAM TEST_PDLAMVE_I8
      IMPLICIT NONE
*
*  Test for PDLAMVE_I8 / PSLAMVE_I8 -- INTEGER*8 unaligned matrix copy.
*
*  Tests:
*    1. Full copy between different block-cyclic layouts (D)
*    2. Full copy between different block-cyclic layouts (S)
*    3. Full copy, same layout (identity redistribution, D)
*    4. Triangular upper copy (D, delegates to legacy PDLACPY)
*    5. Triangular lower copy (D)
*
*  Usage:  mpirun -np 4 ./xdlamve_i8
*
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL
      INTEGER            NPROCS, IAM
      INTEGER            NFAIL, NTEST, NFAIL_G
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
     $   WRITE(*,'(A)') 'test_pdlamve_i8: 2x2 grid'
*
*     Test 1: Full D copy, different layouts
*
      CALL RUN_D( 'A', 20_8, 15_8, 3_8, 4_8, 5_8, 3_8,
     $            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $            'D:full-diff', NTEST, NFAIL )
*
*     Test 2: Full S copy, different layouts
*
      CALL RUN_S( 'A', 20_8, 15_8, 3_8, 4_8, 5_8, 3_8,
     $            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $            'S:full-diff', NTEST, NFAIL )
*
*     Test 3: Full D copy, same layout (identity)
*
      CALL RUN_D( 'A', 30_8, 30_8, 7_8, 7_8, 7_8, 7_8,
     $            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $            'D:full-same', NTEST, NFAIL )
*
*     Test 4: Upper triangular D copy (same layout, exercises PDLACPY)
*
      CALL RUN_D( 'U', 16_8, 16_8, 4_8, 4_8, 4_8, 4_8,
     $            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $            'D:upper     ', NTEST, NFAIL )
*
*     Test 5: Lower triangular D copy (same layout)
*
      CALL RUN_D( 'L', 16_8, 16_8, 4_8, 4_8, 4_8, 4_8,
     $            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $            'D:lower     ', NTEST, NFAIL )
*
*     Test 6: Full D copy with nonzero CSRC
*
      CALL RUN_D2( 'A', 20_8, 15_8, 3_8, 4_8, 5_8, 3_8, 0, 1,
     $             ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $             'D:csrc=1    ', NTEST, NFAIL )
*
*     Test 7: Upper triangular D with nonzero RSRC
*
      CALL RUN_D2( 'U', 16_8, 16_8, 4_8, 4_8, 4_8, 4_8, 1, 0,
     $             ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $             'D:U,rsrc=1  ', NTEST, NFAIL )
*
*     Test 8: Full C copy with imaginary validation
*
      CALL RUN_C( 'A', 20_8, 15_8, 3_8, 4_8, 5_8, 3_8,
     $            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $            'C:full-diff ', NTEST, NFAIL )
*
*     Test 9: Upper triangular C copy with imaginary
*
      CALL RUN_C( 'U', 16_8, 16_8, 4_8, 4_8, 4_8, 4_8,
     $            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $            'C:upper     ', NTEST, NFAIL )
*
*     Test 10: Full Z copy with imaginary validation
*
      CALL RUN_Z( 'A', 20_8, 15_8, 3_8, 4_8, 5_8, 3_8,
     $            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $            'Z:full-diff ', NTEST, NFAIL )
*
*     Test 11: Lower triangular Z copy with imaginary
*
      CALL RUN_Z( 'L', 16_8, 16_8, 4_8, 4_8, 4_8, 4_8,
     $            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $            'Z:lower     ', NTEST, NFAIL )
*
*     Reduce errors globally
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NTEST, NTEST, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PxLAMVE_I8 Test Summary'
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A,I4)') 'Tests run:    ', NTEST
         IF( NFAIL_G .EQ. 0 ) THEN
            WRITE(*,'(A)') 'TEST PASSED OK'
         ELSE
            WRITE(*,'(A,I6,A)') 'TEST FAILED: ',NFAIL_G,' mismatches'
         END IF
         WRITE(*,'(A)') '======================================'
      END IF
*
      CALL BLACS_GRIDEXIT( ICTXT )
      CALL BLACS_EXIT( 0 )
      IF( NFAIL_G .NE. 0 ) STOP 1
*
      END
*
*     ================================================================
*     RUN_D — DOUBLE PRECISION test via PDLAMVE_I8
*     ================================================================
*
      SUBROUTINE RUN_D( UPLO, M8, N8, MBA, NBA, MBB, NBB,
     $                   ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $                   LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      CHARACTER*(*)      UPLO, LABEL
      INTEGER*8          M8, N8, MBA, NBA, MBB, NBB
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR_A, LC_A, LR_B, LC_B
      INTEGER*8          DESCA( 9 ), DESCB( 9 )
      INTEGER*8          I8, J8, GI, GJ, LLD_A, LLD_B
      INTEGER            ERRS
      DOUBLE PRECISION, ALLOCATABLE :: A(:), B(:), DWORK(:)
      LOGICAL            IN_TRI
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCINIT_I8, PDLAMVE_I8
      LOGICAL            LSAME
      EXTERNAL           LSAME
      INTRINSIC          DBLE, MAX, INT
*
      LR_A = NUMROC_I8( M8, MBA, MYROW, 0, NPROW )
      LC_A = NUMROC_I8( N8, NBA, MYCOL, 0, NPCOL )
      LLD_A = MAX( LR_A, 1_8 )
      LR_B = NUMROC_I8( M8, MBB, MYROW, 0, NPROW )
      LC_B = NUMROC_I8( N8, NBB, MYCOL, 0, NPCOL )
      LLD_B = MAX( LR_B, 1_8 )
*
      ERRS = 0
      CALL DESCINIT_I8( DESCA, M8, N8, MBA, NBA, 0, 0, ICTXT,
     $                  LLD_A, ERRS )
      CALL DESCINIT_I8( DESCB, M8, N8, MBB, NBB, 0, 0, ICTXT,
     $                  LLD_B, ERRS )
*
      ALLOCATE( A( MAX( LLD_A * LC_A, 1_8 ) ) )
      ALLOCATE( B( MAX( LLD_B * LC_B, 1_8 ) ) )
      ALLOCATE( DWORK( MAX( LLD_B * LC_B, 1_8 ) ) )
*
*     Fill A: A(gi,gj) = gi*1000 + gj (1-based global indices)
*
      DO J8 = 1, LC_A
         GJ = INDXL2G_I8( J8, NBA, MYCOL, 0, NPCOL )
         DO I8 = 1, LR_A
            GI = INDXL2G_I8( I8, MBA, MYROW, 0, NPROW )
            A( (J8-1)*LLD_A + I8 ) = DBLE( GI*1000 + GJ )
         END DO
      END DO
*
*     Sentinel fill for B
*
      DO I8 = 1, MAX( LLD_B * LC_B, 1_8 )
         B( I8 ) = -999.0D0
      END DO
*
      CALL PDLAMVE_I8( UPLO, M8, N8, A, 1_8, 1_8, DESCA,
     $                 B, 1_8, 1_8, DESCB, DWORK )
*
*     Verify B
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC_B
         GJ = INDXL2G_I8( J8, NBB, MYCOL, 0, NPCOL )
         DO I8 = 1, LR_B
            GI = INDXL2G_I8( I8, MBB, MYROW, 0, NPROW )
*
*           Determine expected value based on UPLO
*
            IN_TRI = .TRUE.
            IF( LSAME( UPLO, 'U' ) ) THEN
               IN_TRI = ( GI .LE. GJ )
            ELSE IF( LSAME( UPLO, 'L' ) ) THEN
               IN_TRI = ( GI .GE. GJ )
            END IF
*
            IF( IN_TRI ) THEN
               IF( B( (J8-1)*LLD_B + I8 ) .NE.
     $             DBLE( GI*1000 + GJ ) ) THEN
                  IF( ERRS .LT. 3 )
     $               WRITE(*,*) 'FAIL ', LABEL, ': g(', GI, ',',
     $                  GJ, ') got=', B((J8-1)*LLD_B+I8)
                  ERRS = ERRS + 1
               END IF
            ELSE
               IF( B( (J8-1)*LLD_B + I8 ) .NE. -999.0D0 ) THEN
                  IF( ERRS .LT. 3 )
     $               WRITE(*,*) 'FAIL ', LABEL, ': non-tri g(',
     $                  GI, ',', GJ, ') modified'
                  ERRS = ERRS + 1
               END IF
            END IF
*
         END DO
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
      DEALLOCATE( A, B, DWORK )
      END
*
*     ================================================================
*     RUN_D2 — DOUBLE PRECISION with explicit RSRC/CSRC
*     ================================================================
*
      SUBROUTINE RUN_D2( UPLO, M8, N8, MBA, NBA, MBB, NBB,
     $                    IRSRC, ICSRC,
     $                    ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $                    LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      CHARACTER*(*)      UPLO, LABEL
      INTEGER*8          M8, N8, MBA, NBA, MBB, NBB
      INTEGER            IRSRC, ICSRC
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR_A, LC_A, LR_B, LC_B
      INTEGER*8          DESCA( 9 ), DESCB( 9 )
      INTEGER*8          I8, J8, GI, GJ, LLD_A, LLD_B
      INTEGER            ERRS
      LOGICAL            IN_TRI
      DOUBLE PRECISION, ALLOCATABLE :: A(:), B(:), DWORK(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCINIT_I8, PDLAMVE_I8
      LOGICAL            LSAME
      EXTERNAL           LSAME
      INTRINSIC          DBLE, MAX, INT
*
      LR_A = NUMROC_I8( M8, MBA, MYROW, IRSRC, NPROW )
      LC_A = NUMROC_I8( N8, NBA, MYCOL, ICSRC, NPCOL )
      LLD_A = MAX( LR_A, 1_8 )
      LR_B = NUMROC_I8( M8, MBB, MYROW, IRSRC, NPROW )
      LC_B = NUMROC_I8( N8, NBB, MYCOL, ICSRC, NPCOL )
      LLD_B = MAX( LR_B, 1_8 )
*
      ERRS = 0
      CALL DESCINIT_I8( DESCA, M8, N8, MBA, NBA, IRSRC, ICSRC,
     $                  ICTXT, LLD_A, ERRS )
      CALL DESCINIT_I8( DESCB, M8, N8, MBB, NBB, IRSRC, ICSRC,
     $                  ICTXT, LLD_B, ERRS )
*
      ALLOCATE( A( MAX( LLD_A * LC_A, 1_8 ) ) )
      ALLOCATE( B( MAX( LLD_B * LC_B, 1_8 ) ) )
      ALLOCATE( DWORK( MAX( LLD_B * LC_B, 1_8 ) ) )
*
      DO J8 = 1, LC_A
         GJ = INDXL2G_I8( J8, NBA, MYCOL, ICSRC, NPCOL )
         DO I8 = 1, LR_A
            GI = INDXL2G_I8( I8, MBA, MYROW, IRSRC, NPROW )
            A( (J8-1)*LLD_A + I8 ) = DBLE( GI*1000 + GJ )
         END DO
      END DO
*
      DO I8 = 1, MAX( LLD_B * LC_B, 1_8 )
         B( I8 ) = -999.0D0
      END DO
*
      CALL PDLAMVE_I8( UPLO, M8, N8, A, 1_8, 1_8, DESCA,
     $                 B, 1_8, 1_8, DESCB, DWORK )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC_B
         GJ = INDXL2G_I8( J8, NBB, MYCOL, ICSRC, NPCOL )
         DO I8 = 1, LR_B
            GI = INDXL2G_I8( I8, MBB, MYROW, IRSRC, NPROW )
            IN_TRI = .TRUE.
            IF( LSAME( UPLO, 'U' ) ) THEN
               IN_TRI = ( GI .LE. GJ )
            ELSE IF( LSAME( UPLO, 'L' ) ) THEN
               IN_TRI = ( GI .GE. GJ )
            END IF
            IF( IN_TRI ) THEN
               IF( B( (J8-1)*LLD_B + I8 ) .NE.
     $             DBLE( GI*1000 + GJ ) ) THEN
                  IF( ERRS .LT. 3 )
     $               WRITE(*,*) 'FAIL ', LABEL, ': g(', GI, ',',
     $                  GJ, ') got=', B((J8-1)*LLD_B+I8)
                  ERRS = ERRS + 1
               END IF
            END IF
         END DO
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
      DEALLOCATE( A, B, DWORK )
      END
*
*     ================================================================
*     RUN_C — COMPLEX test via PCGEMR2D_I8 + PCLACPY_I8
*     ================================================================
*
      SUBROUTINE RUN_C( UPLO, M8, N8, MBA, NBA, MBB, NBB,
     $                   ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $                   LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      CHARACTER*(*)      UPLO, LABEL
      INTEGER*8          M8, N8, MBA, NBA, MBB, NBB
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR_A, LC_A, LR_B, LC_B
      INTEGER*8          DESCA( 9 ), DESCB( 9 )
      INTEGER*8          I8, J8, GI, GJ, LLD_A, LLD_B
      INTEGER            ERRS
      LOGICAL            IN_TRI
      COMPLEX, ALLOCATABLE :: A(:), B(:), DWORK(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCINIT_I8, PCGEMR2D_I8, PCLACPY_I8
      LOGICAL            LSAME
      EXTERNAL           LSAME
      INTRINSIC          REAL, AIMAG, CMPLX, MAX, INT
*
      LR_A = NUMROC_I8( M8, MBA, MYROW, 0, NPROW )
      LC_A = NUMROC_I8( N8, NBA, MYCOL, 0, NPCOL )
      LLD_A = MAX( LR_A, 1_8 )
      LR_B = NUMROC_I8( M8, MBB, MYROW, 0, NPROW )
      LC_B = NUMROC_I8( N8, NBB, MYCOL, 0, NPCOL )
      LLD_B = MAX( LR_B, 1_8 )
*
      ERRS = 0
      CALL DESCINIT_I8( DESCA, M8, N8, MBA, NBA, 0, 0, ICTXT,
     $                  LLD_A, ERRS )
      CALL DESCINIT_I8( DESCB, M8, N8, MBB, NBB, 0, 0, ICTXT,
     $                  LLD_B, ERRS )
*
      ALLOCATE( A( MAX( LLD_A * LC_A, 1_8 ) ) )
      ALLOCATE( B( MAX( LLD_B * LC_B, 1_8 ) ) )
      ALLOCATE( DWORK( MAX( LLD_B * LC_B, 1_8 ) ) )
*
      DO J8 = 1, LC_A
         GJ = INDXL2G_I8( J8, NBA, MYCOL, 0, NPCOL )
         DO I8 = 1, LR_A
            GI = INDXL2G_I8( I8, MBA, MYROW, 0, NPROW )
            A( (J8-1)*LLD_A + I8 ) = CMPLX( REAL( GI*1000 + GJ ),
     $                                       -REAL( GI*1000 + GJ ) )
         END DO
      END DO
*
      DO I8 = 1, MAX( LLD_B * LC_B, 1_8 )
         B( I8 ) = CMPLX( -999.0, -999.0 )
      END DO
*
      CALL PCGEMR2D_I8( M8, N8, A, 1_8, 1_8, DESCA,
     $     DWORK, 1_8, 1_8, DESCB, INT( ICTXT, 8 ) )
      IF( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) THEN
         CALL PCLACPY_I8( UPLO, M8, N8, DWORK, 1_8, 1_8, DESCB,
     $        B, 1_8, 1_8, DESCB )
      ELSE
         DO I8 = 1, MAX( LLD_B * LC_B, 1_8 )
            B( I8 ) = DWORK( I8 )
         END DO
      END IF
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC_B
         GJ = INDXL2G_I8( J8, NBB, MYCOL, 0, NPCOL )
         DO I8 = 1, LR_B
            GI = INDXL2G_I8( I8, MBB, MYROW, 0, NPROW )
            IN_TRI = .TRUE.
            IF( LSAME( UPLO, 'U' ) ) THEN
               IN_TRI = ( GI .LE. GJ )
            ELSE IF( LSAME( UPLO, 'L' ) ) THEN
               IN_TRI = ( GI .GE. GJ )
            END IF
            IF( IN_TRI ) THEN
               IF( REAL( B( (J8-1)*LLD_B + I8 ) ) .NE.
     $             REAL( GI*1000 + GJ ) .OR.
     $             AIMAG( B( (J8-1)*LLD_B + I8 ) ) .NE.
     $             -REAL( GI*1000 + GJ ) ) THEN
                  IF( ERRS .LT. 3 )
     $               WRITE(*,*) 'FAIL ', LABEL, ': g(', GI, ',',
     $                  GJ, ') got=', B((J8-1)*LLD_B+I8)
                  ERRS = ERRS + 1
               END IF
            ELSE
               IF( B( (J8-1)*LLD_B + I8 ) .NE.
     $             CMPLX( -999.0, -999.0 ) ) THEN
                  IF( ERRS .LT. 3 )
     $               WRITE(*,*) 'FAIL ', LABEL, ': non-tri modified'
                  ERRS = ERRS + 1
               END IF
            END IF
         END DO
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
      DEALLOCATE( A, B, DWORK )
      END
*
*     ================================================================
*     RUN_Z — DOUBLE COMPLEX test via PZGEMR2D_I8 + PZLACPY_I8
*     ================================================================
*
      SUBROUTINE RUN_Z( UPLO, M8, N8, MBA, NBA, MBB, NBB,
     $                   ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $                   LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      CHARACTER*(*)      UPLO, LABEL
      INTEGER*8          M8, N8, MBA, NBA, MBB, NBB
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR_A, LC_A, LR_B, LC_B
      INTEGER*8          DESCA( 9 ), DESCB( 9 )
      INTEGER*8          I8, J8, GI, GJ, LLD_A, LLD_B
      INTEGER            ERRS
      LOGICAL            IN_TRI
      DOUBLE COMPLEX, ALLOCATABLE :: A(:), B(:), DWORK(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCINIT_I8, PZGEMR2D_I8, PZLACPY_I8
      LOGICAL            LSAME
      EXTERNAL           LSAME
      INTRINSIC          DBLE, MAX, INT
*
      LR_A = NUMROC_I8( M8, MBA, MYROW, 0, NPROW )
      LC_A = NUMROC_I8( N8, NBA, MYCOL, 0, NPCOL )
      LLD_A = MAX( LR_A, 1_8 )
      LR_B = NUMROC_I8( M8, MBB, MYROW, 0, NPROW )
      LC_B = NUMROC_I8( N8, NBB, MYCOL, 0, NPCOL )
      LLD_B = MAX( LR_B, 1_8 )
*
      ERRS = 0
      CALL DESCINIT_I8( DESCA, M8, N8, MBA, NBA, 0, 0, ICTXT,
     $                  LLD_A, ERRS )
      CALL DESCINIT_I8( DESCB, M8, N8, MBB, NBB, 0, 0, ICTXT,
     $                  LLD_B, ERRS )
*
      ALLOCATE( A( MAX( LLD_A * LC_A, 1_8 ) ) )
      ALLOCATE( B( MAX( LLD_B * LC_B, 1_8 ) ) )
      ALLOCATE( DWORK( MAX( LLD_B * LC_B, 1_8 ) ) )
*
      DO J8 = 1, LC_A
         GJ = INDXL2G_I8( J8, NBA, MYCOL, 0, NPCOL )
         DO I8 = 1, LR_A
            GI = INDXL2G_I8( I8, MBA, MYROW, 0, NPROW )
            A( (J8-1)*LLD_A + I8 ) = DCMPLX( DBLE( GI*1000 + GJ ),
     $                                        -DBLE( GI*1000 + GJ ) )
         END DO
      END DO
*
      DO I8 = 1, MAX( LLD_B * LC_B, 1_8 )
         B( I8 ) = DCMPLX( -999.0D0, -999.0D0 )
      END DO
*
      CALL PZGEMR2D_I8( M8, N8, A, 1_8, 1_8, DESCA,
     $     DWORK, 1_8, 1_8, DESCB, INT( ICTXT, 8 ) )
      IF( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) THEN
         CALL PZLACPY_I8( UPLO, M8, N8, DWORK, 1_8, 1_8, DESCB,
     $        B, 1_8, 1_8, DESCB )
      ELSE
         DO I8 = 1, MAX( LLD_B * LC_B, 1_8 )
            B( I8 ) = DWORK( I8 )
         END DO
      END IF
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC_B
         GJ = INDXL2G_I8( J8, NBB, MYCOL, 0, NPCOL )
         DO I8 = 1, LR_B
            GI = INDXL2G_I8( I8, MBB, MYROW, 0, NPROW )
            IN_TRI = .TRUE.
            IF( LSAME( UPLO, 'U' ) ) THEN
               IN_TRI = ( GI .LE. GJ )
            ELSE IF( LSAME( UPLO, 'L' ) ) THEN
               IN_TRI = ( GI .GE. GJ )
            END IF
            IF( IN_TRI ) THEN
               IF( DBLE( B( (J8-1)*LLD_B + I8 ) ) .NE.
     $             DBLE( GI*1000 + GJ ) .OR.
     $             DIMAG( B( (J8-1)*LLD_B + I8 ) ) .NE.
     $             -DBLE( GI*1000 + GJ ) ) THEN
                  IF( ERRS .LT. 3 )
     $               WRITE(*,*) 'FAIL ', LABEL, ': g(', GI, ',',
     $                  GJ, ') got=', B((J8-1)*LLD_B+I8)
                  ERRS = ERRS + 1
               END IF
            ELSE
               IF( B( (J8-1)*LLD_B + I8 ) .NE.
     $             DCMPLX( -999.0D0, -999.0D0 ) ) THEN
                  IF( ERRS .LT. 3 )
     $               WRITE(*,*) 'FAIL ', LABEL, ': non-tri modified'
                  ERRS = ERRS + 1
               END IF
            END IF
         END DO
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
      DEALLOCATE( A, B, DWORK )
      END
*
*     ================================================================
*     RUN_S — REAL test via PSLAMVE_I8
*     ================================================================
*
      SUBROUTINE RUN_S( UPLO, M8, N8, MBA, NBA, MBB, NBB,
     $                   ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM,
     $                   LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      CHARACTER*(*)      UPLO, LABEL
      INTEGER*8          M8, N8, MBA, NBA, MBB, NBB
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR_A, LC_A, LR_B, LC_B
      INTEGER*8          DESCA( 9 ), DESCB( 9 )
      INTEGER*8          I8, J8, GI, GJ, LLD_A, LLD_B
      INTEGER            ERRS
      REAL, ALLOCATABLE :: A(:), B(:), DWORK(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCINIT_I8, PSLAMVE_I8
      INTRINSIC          REAL, MAX, INT
*
      LR_A = NUMROC_I8( M8, MBA, MYROW, 0, NPROW )
      LC_A = NUMROC_I8( N8, NBA, MYCOL, 0, NPCOL )
      LLD_A = MAX( LR_A, 1_8 )
      LR_B = NUMROC_I8( M8, MBB, MYROW, 0, NPROW )
      LC_B = NUMROC_I8( N8, NBB, MYCOL, 0, NPCOL )
      LLD_B = MAX( LR_B, 1_8 )
*
      ERRS = 0
      CALL DESCINIT_I8( DESCA, M8, N8, MBA, NBA, 0, 0, ICTXT,
     $                  LLD_A, ERRS )
      CALL DESCINIT_I8( DESCB, M8, N8, MBB, NBB, 0, 0, ICTXT,
     $                  LLD_B, ERRS )
*
      ALLOCATE( A( MAX( LLD_A * LC_A, 1_8 ) ) )
      ALLOCATE( B( MAX( LLD_B * LC_B, 1_8 ) ) )
      ALLOCATE( DWORK( MAX( LLD_B * LC_B, 1_8 ) ) )
*
      DO J8 = 1, LC_A
         GJ = INDXL2G_I8( J8, NBA, MYCOL, 0, NPCOL )
         DO I8 = 1, LR_A
            GI = INDXL2G_I8( I8, MBA, MYROW, 0, NPROW )
            A( (J8-1)*LLD_A + I8 ) = REAL( GI*1000 + GJ )
         END DO
      END DO
*
      DO I8 = 1, MAX( LLD_B * LC_B, 1_8 )
         B( I8 ) = -999.0
      END DO
*
      CALL PSLAMVE_I8( UPLO, M8, N8, A, 1_8, 1_8, DESCA,
     $                 B, 1_8, 1_8, DESCB, DWORK )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC_B
         GJ = INDXL2G_I8( J8, NBB, MYCOL, 0, NPCOL )
         DO I8 = 1, LR_B
            GI = INDXL2G_I8( I8, MBB, MYROW, 0, NPROW )
            IF( B( (J8-1)*LLD_B + I8 ) .NE.
     $          REAL( GI*1000 + GJ ) ) THEN
               IF( ERRS .LT. 3 )
     $            WRITE(*,*) 'FAIL ', LABEL, ': g(', GI, ',',
     $               GJ, ') got=', B((J8-1)*LLD_B+I8)
               ERRS = ERRS + 1
            END IF
         END DO
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
      DEALLOCATE( A, B, DWORK )
      END
