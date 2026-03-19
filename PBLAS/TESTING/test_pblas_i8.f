      PROGRAM TEST_PBLAS_I8
      IMPLICIT NONE
*
*  Per-kernel tests for PBLAS I8 entry points.
*  Compares each I8 wrapper against its legacy counterpart.
*
*  Usage: mpirun -np 4 ./xpblas_i8
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
      IF( NPROCS .LT. 4 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'Need at least 4 processes'
         CALL BLACS_EXIT( 0 )
         STOP 1
      END IF
*
      NPROW = 2
      NPCOL = 2
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      IF( MYROW .LT. 0 ) THEN
         CALL BLACS_EXIT( 0 )
         STOP
      END IF
*
      NFAIL = 0
      NTEST = 0
*
      IF( IAM .EQ. 0 )
     $   WRITE(*,'(A)') 'test_pblas_i8: 2x2 grid'
*
*     === PDAXPY / PDSCAL ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 1: AXPY/SCAL ---'
      CALL TEST_DAXPY_DSCAL( 100_8, 4_8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PDDOT ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 1: DOT ---'
      CALL TEST_DDOT( 100_8, 4_8, ICTXT, NPROW, NPCOL,
     $                MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PDGEMV ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 2: GEMV ---'
      CALL TEST_DGEMV( 30_8, 20_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PDSYMV ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 2: SYMV ---'
      CALL TEST_DSYMV( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PCHEMV (complex Hermitian) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 2: CHEMV ---'
      CALL TEST_CHEMV( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PDNRM2 ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 1: NRM2 ---'
      CALL TEST_DNRM2( 100_8, 4_8, ICTXT, NPROW, NPCOL,
     $                  MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PCDOTC (complex conjugate dot) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 1: CDOTC ---'
      CALL TEST_CDOTC( 80_8, 4_8, ICTXT, NPROW, NPCOL,
     $                  MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PDSYR2K (double symmetric rank-2k) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 3: SYR2K ---'
      CALL TEST_DSYR2K( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PCHER2K (complex Hermitian rank-2k) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 3: HER2K ---'
      CALL TEST_CHER2K( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PDSYRK (double symmetric rank-k) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 3: SYRK ---'
      CALL TEST_DSYRK( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                  MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PDTRSM (double triangular solve) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 3: TRSM ---'
      CALL TEST_DTRSM( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                  MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PDGEMM (double general multiply) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 3: GEMM ---'
      CALL TEST_DGEMM( 30_8, 20_8, 4_8, ICTXT, NPROW, NPCOL,
     $                  MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PSTRSM (single triangular solve) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 3: STRSM ---'
      CALL TEST_STRSM( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                  MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PCTRSM (complex triangular solve) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 3: CTRSM ---'
      CALL TEST_CTRSM( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                  MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PCHERK (complex Hermitian rank-k) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 3: HERK ---'
      CALL TEST_CHERK( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                  MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     === PCGEMM (complex general multiply) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Level 3: CGEMM ---'
      CALL TEST_CGEMM( 30_8, 20_8, 4_8, ICTXT, NPROW, NPCOL,
     $                  MYROW, MYCOL, IAM, NTEST, NFAIL )
*
*     Summary
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NFAIL_G, NFAIL_G, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PBLAS I8 Kernel Test Summary'
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
      CALL BLACS_GRIDEXIT( ICTXT )
      CALL BLACS_EXIT( 0 )
      IF( NFAIL_G .NE. 0 ) STOP 1
*
      END
*
*     ================================================================
*     TEST_DAXPY_DSCAL — PDAXPY_I8 + PDSCAL_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_DAXPY_DSCAL( N8, NB8, ICTXT, NPROW, NPCOL,
     $                              MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LLD8, DESCA8( 9 ), I8
      INTEGER            DESCA4( 9 ), INFO, ERRS, N4, NB4
      DOUBLE PRECISION, ALLOCATABLE :: X(:), Y(:), XREF(:), YREF(:)
      DOUBLE PRECISION   ALPHA
*
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PDAXPY_I8, PDAXPY, PDSCAL_I8, PDSCAL
      INTRINSIC          DBLE, MAX, INT
*
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, 1_8, NB8, 1_8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, 1, NB4, 1, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( X( MAX( LR, 1_8 ) ) )
      ALLOCATE( Y( MAX( LR, 1_8 ) ) )
      ALLOCATE( XREF( MAX( LR, 1_8 ) ) )
      ALLOCATE( YREF( MAX( LR, 1_8 ) ) )
*
*     Fill X, Y with pattern
      DO I8 = 1, LR
         X( I8 ) = DBLE( I8 )
         Y( I8 ) = DBLE( I8 * 2 )
         XREF( I8 ) = X( I8 )
         YREF( I8 ) = Y( I8 )
      END DO
*
*     Test PDAXPY: Y = Y + alpha*X
      ALPHA = 3.0D0
      CALL PDAXPY_I8( N8, ALPHA, X, 1_8, 1_8, DESCA8, 1_8,
     $                Y, 1_8, 1_8, DESCA8, 1_8 )
      CALL PDAXPY( N4, ALPHA, XREF, 1, 1, DESCA4, 1,
     $             YREF, 1, 1, DESCA4, 1 )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO I8 = 1, LR
         IF( Y( I8 ) .NE. YREF( I8 ) ) ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PDAXPY_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PDAXPY_I8: FAILED ', ERRS
         END IF
      END IF
*
*     Test PDSCAL: Y = alpha*Y
      ALPHA = 0.5D0
      CALL PDSCAL_I8( N8, ALPHA, Y, 1_8, 1_8, DESCA8, 1_8 )
      CALL PDSCAL( N4, ALPHA, YREF, 1, 1, DESCA4, 1 )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO I8 = 1, LR
         IF( Y( I8 ) .NE. YREF( I8 ) ) ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PDSCAL_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PDSCAL_I8: FAILED ', ERRS
         END IF
      END IF
*
      DEALLOCATE( X, Y, XREF, YREF )
      END
*
*     ================================================================
*     TEST_DDOT — PDDOT_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_DDOT( N8, NB8, ICTXT, NPROW, NPCOL,
     $                       MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LLD8, DESCA8( 9 ), I8
      INTEGER            DESCA4( 9 ), INFO, N4, NB4
      DOUBLE PRECISION, ALLOCATABLE :: X(:), Y(:)
      DOUBLE PRECISION   DOT8, DOT4
*
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT, PDDOT_I8, PDDOT
      INTRINSIC          DBLE, MAX, INT
*
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, 1_8, NB8, 1_8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, 1, NB4, 1, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( X( MAX( LR, 1_8 ) ) )
      ALLOCATE( Y( MAX( LR, 1_8 ) ) )
      DO I8 = 1, LR
         X( I8 ) = DBLE( I8 )
         Y( I8 ) = 1.0D0
      END DO
*
      CALL PDDOT_I8( N8, DOT8, X, 1_8, 1_8, DESCA8, 1_8,
     $               Y, 1_8, 1_8, DESCA8, 1_8 )
      CALL PDDOT( N4, DOT4, X, 1, 1, DESCA4, 1,
     $            Y, 1, 1, DESCA4, 1 )
*
      NTEST = NTEST + 1
      IF( DOT8 .NE. DOT4 ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A,F12.1,A,F12.1)') '  PDDOT_I8: FAILED i8=',
     $            DOT8, ' ref=', DOT4
      ELSE
         IF( IAM .EQ. 0 ) WRITE(*,'(A)') '  PDDOT_I8: PASSED'
      END IF
*
      DEALLOCATE( X, Y )
      END
*
*     ================================================================
*     TEST_DGEMV — PDGEMV_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_DGEMV( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LRX, LLDA8, LLDX8
      INTEGER*8          DESCA8( 9 ), DESCX8( 9 ), DESCY8( 9 )
      INTEGER*8          I8, J8, GI, GJ
      INTEGER            DESCA4( 9 ), DESCX4( 9 ), DESCY4( 9 )
      INTEGER            INFO, ERRS, M4, N4, NB4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), X(:), Y(:), YREF(:)
      DOUBLE PRECISION   ALPHA, BETA
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PDGEMV_I8, PDGEMV
      INTRINSIC          DBLE, MAX, INT
*
      M4 = INT( M8 )
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      LRA = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA8 = MAX( LRA, 1_8 )
      LRX = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LLDX8 = MAX( LRX, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA8, INFO )
      CALL DESCINIT_I8( DESCX8, N8, 1_8, NB8, 1_8, 0, 0, ICTXT,
     $                  LLDX8, INFO )
      CALL DESCINIT_I8( DESCY8, M8, 1_8, NB8, 1_8, 0, 0, ICTXT,
     $                  LLDA8, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA8 ), INFO )
      CALL DESCINIT( DESCX4, N4, 1, NB4, 1, 0, 0, ICTXT,
     $               INT( LLDX8 ), INFO )
      CALL DESCINIT( DESCY4, M4, 1, NB4, 1, 0, 0, ICTXT,
     $               INT( LLDA8 ), INFO )
*
      ALLOCATE( A( MAX( LLDA8 * LCA, 1_8 ) ) )
      ALLOCATE( X( MAX( LRX, 1_8 ) ) )
      ALLOCATE( Y( MAX( LRA, 1_8 ) ) )
      ALLOCATE( YREF( MAX( LRA, 1_8 ) ) )
*
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            A( (J8-1)*LLDA8 + I8 ) = 1.0D0
         END DO
      END DO
      DO I8 = 1, MAX( LRX, 1_8 )
         X( I8 ) = 2.0D0
      END DO
      DO I8 = 1, MAX( LRA, 1_8 )
         Y( I8 ) = 0.0D0
         YREF( I8 ) = 0.0D0
      END DO
*
      ALPHA = 1.0D0
      BETA = 0.0D0
      CALL PDGEMV_I8( 'N', M8, N8, ALPHA, A, 1_8, 1_8, DESCA8,
     $                X, 1_8, 1_8, DESCX8, 1_8, BETA,
     $                Y, 1_8, 1_8, DESCY8, 1_8 )
      CALL PDGEMV( 'N', M4, N4, ALPHA, A, 1, 1, DESCA4,
     $             X, 1, 1, DESCX4, 1, BETA,
     $             YREF, 1, 1, DESCY4, 1 )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO I8 = 1, LRA
         IF( Y( I8 ) .NE. YREF( I8 ) ) ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PDGEMV_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PDGEMV_I8: FAILED ', ERRS
         END IF
      END IF
*
      DEALLOCATE( A, X, Y, YREF )
      END
*
*     ================================================================
*     TEST_DSYMV — PDSYMV_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_DSYMV( N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LLD8, DESCA8( 9 ), DESCX8( 9 )
      INTEGER*8          I8, J8, GI, GJ
      INTEGER            DESCA4( 9 ), DESCX4( 9 )
      INTEGER            INFO, ERRS, N4, NB4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), X(:), Y(:), YREF(:)
      DOUBLE PRECISION   ALPHA, BETA
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PDSYMV_I8, PDSYMV
      INTRINSIC          DBLE, MAX, INT, ABS
*
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCX8, N8, 1_8, NB8, 1_8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCX4, N4, 1, NB4, 1, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8 * LC, 1_8 ) ) )
      ALLOCATE( X( MAX( LR, 1_8 ) ) )
      ALLOCATE( Y( MAX( LR, 1_8 ) ) )
      ALLOCATE( YREF( MAX( LR, 1_8 ) ) )
*
*     Symmetric matrix: A(i,j) = 1/(1+|i-j|)
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            A( (J8-1)*LLD8 + I8 ) = 1.0D0 / DBLE(1+ABS(GI-GJ))
         END DO
      END DO
      DO I8 = 1, LR
         X( I8 ) = 1.0D0
         Y( I8 ) = 0.0D0
         YREF( I8 ) = 0.0D0
      END DO
*
      ALPHA = 1.0D0
      BETA = 0.0D0
      CALL PDSYMV_I8( 'L', N8, ALPHA, A, 1_8, 1_8, DESCA8,
     $                X, 1_8, 1_8, DESCX8, 1_8, BETA,
     $                Y, 1_8, 1_8, DESCX8, 1_8 )
      CALL PDSYMV( 'L', N4, ALPHA, A, 1, 1, DESCA4,
     $             X, 1, 1, DESCX4, 1, BETA,
     $             YREF, 1, 1, DESCX4, 1 )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO I8 = 1, LR
         IF( Y( I8 ) .NE. YREF( I8 ) ) ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PDSYMV_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PDSYMV_I8: FAILED ', ERRS
         END IF
      END IF
*
      DEALLOCATE( A, X, Y, YREF )
      END
*
*     ================================================================
*     TEST_CHEMV — PCHEMV_I8 vs legacy (complex, imaginary-sensitive)
*     ================================================================
*
      SUBROUTINE TEST_CHEMV( N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LLD8, DESCA8( 9 ), DESCX8( 9 )
      INTEGER*8          I8, J8, GI, GJ
      INTEGER            DESCA4( 9 ), DESCX4( 9 )
      INTEGER            INFO, ERRS, N4, NB4
      COMPLEX, ALLOCATABLE :: A(:), X(:), Y(:), YREF(:)
      COMPLEX            ALPHA, BETA
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PCHEMV_I8, PCHEMV
      INTRINSIC          CMPLX, REAL, AIMAG, MAX, INT, ABS
*
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCX8, N8, 1_8, NB8, 1_8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCX4, N4, 1, NB4, 1, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8 * LC, 1_8 ) ) )
      ALLOCATE( X( MAX( LR, 1_8 ) ) )
      ALLOCATE( Y( MAX( LR, 1_8 ) ) )
      ALLOCATE( YREF( MAX( LR, 1_8 ) ) )
*
*     Hermitian matrix: diag real, off-diag has imaginary
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) = CMPLX( REAL( N8 ), 0.0 )
            ELSE
               A( (J8-1)*LLD8 + I8 ) = CMPLX(
     $            1.0 / REAL( 1+ABS(GI-GJ) ),
     $            REAL( GI-GJ ) * 0.1 )
            END IF
         END DO
      END DO
      DO I8 = 1, LR
         X( I8 ) = CMPLX( 1.0, 0.5 )
         Y( I8 ) = CMPLX( 0.0, 0.0 )
         YREF( I8 ) = CMPLX( 0.0, 0.0 )
      END DO
*
      ALPHA = CMPLX( 1.0, 0.0 )
      BETA = CMPLX( 0.0, 0.0 )
      CALL PCHEMV_I8( 'L', N8, ALPHA, A, 1_8, 1_8, DESCA8,
     $                X, 1_8, 1_8, DESCX8, 1_8, BETA,
     $                Y, 1_8, 1_8, DESCX8, 1_8 )
      CALL PCHEMV( 'L', N4, ALPHA, A, 1, 1, DESCA4,
     $             X, 1, 1, DESCX4, 1, BETA,
     $             YREF, 1, 1, DESCX4, 1 )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO I8 = 1, LR
         IF( REAL( Y( I8 ) ) .NE. REAL( YREF( I8 ) ) .OR.
     $       AIMAG( Y( I8 ) ) .NE. AIMAG( YREF( I8 ) ) )
     $      ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PCHEMV_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PCHEMV_I8: FAILED ', ERRS
         END IF
      END IF
*
      DEALLOCATE( A, X, Y, YREF )
      END
*
*     ================================================================
*     TEST_DNRM2 — PDNRM2_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_DNRM2( N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
      INTEGER*8          LR, LLD8, DESCA8( 9 ), I8
      INTEGER            DESCA4( 9 ), INFO, N4, NB4
      DOUBLE PRECISION, ALLOCATABLE :: X(:)
      DOUBLE PRECISION   NRM8, NRM4
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT, PDNRM2_I8, PDNRM2
      INTRINSIC          DBLE, MAX, INT
*
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LLD8 = MAX( LR, 1_8 )
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, 1_8, NB8, 1_8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, 1, NB4, 1, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      ALLOCATE( X( MAX( LR, 1_8 ) ) )
      DO I8 = 1, LR
         X( I8 ) = 1.0D0
      END DO
*
      CALL PDNRM2_I8( N8, NRM8, X, 1_8, 1_8, DESCA8, 1_8 )
      CALL PDNRM2( N4, NRM4, X, 1, 1, DESCA4, 1 )
*
      NTEST = NTEST + 1
      IF( NRM8 .NE. NRM4 ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  PDNRM2_I8: FAILED'
      ELSE
         IF( IAM .EQ. 0 ) WRITE(*,'(A)') '  PDNRM2_I8: PASSED'
      END IF
      DEALLOCATE( X )
      END
*
*     ================================================================
*     TEST_CDOTC — PCDOTC_I8 vs legacy (complex conjugate dot)
*     ================================================================
*
      SUBROUTINE TEST_CDOTC( N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
      INTEGER*8          LR, LLD8, DESCA8( 9 ), I8
      INTEGER            DESCA4( 9 ), INFO, N4, NB4
      COMPLEX, ALLOCATABLE :: X(:), Y(:)
      COMPLEX            DOT8, DOT4
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT, PCDOTC_I8, PCDOTC
      INTRINSIC          CMPLX, REAL, AIMAG, MAX, INT
*
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LLD8 = MAX( LR, 1_8 )
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, 1_8, NB8, 1_8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, 1, NB4, 1, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      ALLOCATE( X( MAX( LR, 1_8 ) ) )
      ALLOCATE( Y( MAX( LR, 1_8 ) ) )
      DO I8 = 1, LR
         X( I8 ) = CMPLX( 1.0, 0.5 )
         Y( I8 ) = CMPLX( 2.0, -1.0 )
      END DO
*
      CALL PCDOTC_I8( N8, DOT8, X, 1_8, 1_8, DESCA8, 1_8,
     $                Y, 1_8, 1_8, DESCA8, 1_8 )
      CALL PCDOTC( N4, DOT4, X, 1, 1, DESCA4, 1,
     $             Y, 1, 1, DESCA4, 1 )
*
      NTEST = NTEST + 1
      IF( REAL( DOT8 ) .NE. REAL( DOT4 ) .OR.
     $    AIMAG( DOT8 ) .NE. AIMAG( DOT4 ) ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  PCDOTC_I8: FAILED'
      ELSE
         IF( IAM .EQ. 0 ) WRITE(*,'(A)') '  PCDOTC_I8: PASSED'
      END IF
      DEALLOCATE( X, Y )
      END
*
*     ================================================================
*     TEST_DSYR2K — PDSYR2K_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_DSYR2K( N8, NB8, ICTXT, NPROW, NPCOL,
     $                         MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LLD8, K8
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 )
      INTEGER*8          I8, J8, GI, GJ
      INTEGER            DESCA4( 9 ), DESCB4( 9 )
      INTEGER            INFO, ERRS, N4, NB4, K4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), B(:), C(:), CREF(:)
      DOUBLE PRECISION   ALPHA, BETA
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PDSYR2K_I8, PDSYR2K
      INTRINSIC          DBLE, MAX, INT, ABS
*
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      K8 = NB8
      K4 = NB4
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCB8, N8, K8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCB4, N4, K4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( B( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( C( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( CREF( MAX( LLD8*LC, 1_8 ) ) )
*
*     Fill A (N x K), B (N x K), C (symmetric N x N)
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            A( (J8-1)*LLD8 + I8 ) = DBLE( GI + GJ )
            B( (J8-1)*LLD8 + I8 ) = DBLE( GI - GJ + N8 )
            C( (J8-1)*LLD8 + I8 ) =
     $           1.0D0 / DBLE( 1 + ABS(GI-GJ) )
            CREF( (J8-1)*LLD8 + I8 ) = C( (J8-1)*LLD8 + I8 )
         END DO
      END DO
*
      ALPHA = -1.0D0
      BETA = 1.0D0
      CALL PDSYR2K_I8( 'L', 'N', N8, K8, ALPHA,
     $     A, 1_8, 1_8, DESCB8,
     $     B, 1_8, 1_8, DESCB8, BETA,
     $     C, 1_8, 1_8, DESCA8 )
      CALL PDSYR2K( 'L', 'N', N4, K4, ALPHA,
     $     A, 1, 1, DESCB4,
     $     B, 1, 1, DESCB4, BETA,
     $     CREF, 1, 1, DESCA4 )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC
         DO I8 = 1, LR
            IF( C( (J8-1)*LLD8+I8 ) .NE.
     $          CREF( (J8-1)*LLD8+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PDSYR2K_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PDSYR2K_I8: FAILED ', ERRS
         END IF
      END IF
*
      DEALLOCATE( A, B, C, CREF )
      END
*
*     ================================================================
*     TEST_CHER2K — PCHER2K_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_CHER2K( N8, NB8, ICTXT, NPROW, NPCOL,
     $                         MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LLD8, K8
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 )
      INTEGER*8          I8, J8, GI, GJ
      INTEGER            DESCA4( 9 ), DESCB4( 9 )
      INTEGER            INFO, ERRS, N4, NB4, K4
      COMPLEX, ALLOCATABLE :: A(:), B(:), C(:), CREF(:)
      COMPLEX            ALPHA
      REAL               BETA
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PCHER2K_I8, PCHER2K
      INTRINSIC          CMPLX, REAL, MAX, INT, ABS
*
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      K8 = NB8
      K4 = NB4
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCB8, N8, K8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCB4, N4, K4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( B( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( C( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( CREF( MAX( LLD8*LC, 1_8 ) ) )
*
*     Fill A (N x K), B (N x K), C (Hermitian N x N)
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            A( (J8-1)*LLD8 + I8 ) = CMPLX( REAL(GI+GJ),
     $                                       REAL(GI-GJ) )
            B( (J8-1)*LLD8 + I8 ) = CMPLX( REAL(GI),
     $                                       REAL(GJ) )
            IF( GI .EQ. GJ ) THEN
               C( (J8-1)*LLD8 + I8 ) = CMPLX( REAL(2*N8), 0.0 )
            ELSE
               C( (J8-1)*LLD8 + I8 ) = CMPLX(
     $              1.0 / REAL(1+ABS(GI-GJ)),
     $              REAL(GI-GJ) / REAL(N8) )
            END IF
            CREF( (J8-1)*LLD8 + I8 ) = C( (J8-1)*LLD8 + I8 )
         END DO
      END DO
*
      ALPHA = CMPLX( -1.0, 0.0 )
      BETA = 1.0
      CALL PCHER2K_I8( 'L', 'N', N8, K8, ALPHA,
     $     A, 1_8, 1_8, DESCB8,
     $     B, 1_8, 1_8, DESCB8, BETA,
     $     C, 1_8, 1_8, DESCA8 )
      CALL PCHER2K( 'L', 'N', N4, K4, ALPHA,
     $     A, 1, 1, DESCB4,
     $     B, 1, 1, DESCB4, BETA,
     $     CREF, 1, 1, DESCA4 )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC
         DO I8 = 1, LR
            IF( C( (J8-1)*LLD8+I8 ) .NE.
     $          CREF( (J8-1)*LLD8+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PCHER2K_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PCHER2K_I8: FAILED ', ERRS
         END IF
      END IF
*
      DEALLOCATE( A, B, C, CREF )
      END
*
*     ================================================================
*     TEST_DSYRK — PDSYRK_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_DSYRK( N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LLD8, K8, I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCC8( 9 )
      INTEGER            DESCA4( 9 ), DESCC4( 9 )
      INTEGER            INFO, ERRS, N4, NB4, K4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), C(:), CREF(:)
      DOUBLE PRECISION   ALPHA, BETA
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT, PDSYRK_I8, PDSYRK
      INTRINSIC          DBLE, MAX, INT, ABS
*
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      K8 = NB8
      K4 = NB4
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCC8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCA8, N8, K8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCC4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCA4, N4, K4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( C( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( CREF( MAX( LLD8*LC, 1_8 ) ) )
*
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            A( (J8-1)*LLD8 + I8 ) = DBLE( GI + GJ )
            C( (J8-1)*LLD8 + I8 ) =
     $           1.0D0 / DBLE( 1 + ABS(GI-GJ) )
            CREF( (J8-1)*LLD8 + I8 ) = C( (J8-1)*LLD8 + I8 )
         END DO
      END DO
*
      ALPHA = 1.0D0
      BETA = 1.0D0
      CALL PDSYRK_I8( 'L', 'N', N8, K8, ALPHA,
     $     A, 1_8, 1_8, DESCA8, BETA,
     $     C, 1_8, 1_8, DESCC8 )
      CALL PDSYRK( 'L', 'N', N4, K4, ALPHA,
     $     A, 1, 1, DESCA4, BETA,
     $     CREF, 1, 1, DESCC4 )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC
         DO I8 = 1, LR
            IF( C( (J8-1)*LLD8+I8 ) .NE.
     $          CREF( (J8-1)*LLD8+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PDSYRK_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PDSYRK_I8: FAILED ', ERRS
         END IF
      END IF
      DEALLOCATE( A, C, CREF )
      END
*
*     ================================================================
*     TEST_DTRSM — PDTRSM_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_DTRSM( N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LLD8, I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 )
      INTEGER            DESCA4( 9 ), DESCB4( 9 )
      INTEGER            INFO, ERRS, N4, NB4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), B(:), BREF(:)
      DOUBLE PRECISION   ALPHA
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT, PDTRSM_I8, PDTRSM
      INTRINSIC          DBLE, MAX, INT, ABS
*
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCB8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCB4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( B( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( BREF( MAX( LLD8*LC, 1_8 ) ) )
*
*     Upper triangular with strong diagonal
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) = DBLE( N8 )
            ELSE IF( GI .LT. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) =
     $              1.0D0 / DBLE( 1 + ABS(GI-GJ) )
            ELSE
               A( (J8-1)*LLD8 + I8 ) = 0.0D0
            END IF
            B( (J8-1)*LLD8 + I8 ) = 1.0D0
            BREF( (J8-1)*LLD8 + I8 ) = 1.0D0
         END DO
      END DO
*
      ALPHA = 1.0D0
      CALL PDTRSM_I8( 'L', 'U', 'N', 'N', N8, N8, ALPHA,
     $     A, 1_8, 1_8, DESCA8,
     $     B, 1_8, 1_8, DESCB8 )
      CALL PDTRSM( 'L', 'U', 'N', 'N', N4, N4, ALPHA,
     $     A, 1, 1, DESCA4,
     $     BREF, 1, 1, DESCB4 )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC
         DO I8 = 1, LR
            IF( B( (J8-1)*LLD8+I8 ) .NE.
     $          BREF( (J8-1)*LLD8+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PDTRSM_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PDTRSM_I8: FAILED ', ERRS
         END IF
      END IF
      DEALLOCATE( A, B, BREF )
      END
*
*     ================================================================
*     TEST_DGEMM — PDGEMM_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_DGEMM( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRM, LCM, LRN, LCN, LLDM, LLDN, K8
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 ), DESCC8( 9 )
      INTEGER            DESCA4( 9 ), DESCB4( 9 ), DESCC4( 9 )
      INTEGER            INFO, ERRS, M4, N4, NB4, K4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), B(:), C(:), CREF(:)
      DOUBLE PRECISION   ALPHA, BETA
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT, PDGEMM_I8, PDGEMM
      INTRINSIC          DBLE, MAX, INT
*
      M4 = INT( M8 )
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      K8 = NB8
      K4 = NB4
*
      LRM = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCM = NUMROC_I8( M8, NB8, MYCOL, 0, NPCOL )
      LRN = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LCN = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDM = MAX( LRM, 1_8 )
      LLDN = MAX( LRN, 1_8 )
*
      INFO = 0
*     A is M x K
      CALL DESCINIT_I8( DESCA8, M8, K8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDM, INFO )
*     B is K x N
      CALL DESCINIT_I8( DESCB8, K8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  MAX( NUMROC_I8(K8,NB8,MYROW,0,NPROW), 1_8 ),
     $                  INFO )
*     C is M x N
      CALL DESCINIT_I8( DESCC8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDM, INFO )
      CALL DESCINIT( DESCA4, M4, K4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDM ), INFO )
      CALL DESCINIT( DESCB4, K4, N4, NB4, NB4, 0, 0, ICTXT,
     $               MAX(NUMROC(K4,NB4,MYROW,0,NPROW),1), INFO )
      CALL DESCINIT( DESCC4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDM ), INFO )
*
      ALLOCATE( A( MAX( LLDM*LCM, 1_8 ) ) )
      ALLOCATE( B( MAX( LLDN*LCN, 1_8 ) ) )
      ALLOCATE( C( MAX( LLDM*LCN, 1_8 ) ) )
      ALLOCATE( CREF( MAX( LLDM*LCN, 1_8 ) ) )
*
      DO I8 = 1, MAX( LLDM*LCM, 1_8 )
         A( I8 ) = 1.0D0
      END DO
      DO I8 = 1, MAX( LLDN*LCN, 1_8 )
         B( I8 ) = 2.0D0
      END DO
      DO I8 = 1, MAX( LLDM*LCN, 1_8 )
         C( I8 ) = 0.0D0
         CREF( I8 ) = 0.0D0
      END DO
*
      ALPHA = 1.0D0
      BETA = 0.0D0
      CALL PDGEMM_I8( 'N', 'N', M8, N8, K8, ALPHA,
     $     A, 1_8, 1_8, DESCA8,
     $     B, 1_8, 1_8, DESCB8, BETA,
     $     C, 1_8, 1_8, DESCC8 )
      CALL PDGEMM( 'N', 'N', M4, N4, K4, ALPHA,
     $     A, 1, 1, DESCA4,
     $     B, 1, 1, DESCB4, BETA,
     $     CREF, 1, 1, DESCC4 )
*
      NTEST = NTEST + 1
      ERRS = 0
      DO I8 = 1, MAX( LLDM*LCN, 1_8 )
         IF( C( I8 ) .NE. CREF( I8 ) ) ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PDGEMM_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PDGEMM_I8: FAILED ', ERRS
         END IF
      END IF
      DEALLOCATE( A, B, C, CREF )
      END
*
*     ================================================================
*     TEST_STRSM — PSTRSM_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_STRSM( N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
      INTEGER*8          LR, LC, LLD8, I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 )
      INTEGER            DESCA4( 9 ), DESCB4( 9 ), INFO, ERRS, N4, NB4
      REAL, ALLOCATABLE :: A(:), B(:), BREF(:)
      REAL               ALPHA
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC, DESCINIT_I8, DESCINIT,
     $                   PSTRSM_I8, PSTRSM
      INTRINSIC          REAL, MAX, INT, ABS
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCB8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCB4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      ALLOCATE( A( MAX(LLD8*LC,1_8) ), B( MAX(LLD8*LC,1_8) ),
     $          BREF( MAX(LLD8*LC,1_8) ) )
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) = REAL( N8 )
            ELSE IF( GI .LT. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) =
     $              1.0 / REAL( 1 + ABS(GI-GJ) )
            ELSE
               A( (J8-1)*LLD8 + I8 ) = 0.0
            END IF
            B( (J8-1)*LLD8 + I8 ) = 1.0
            BREF( (J8-1)*LLD8 + I8 ) = 1.0
         END DO
      END DO
      ALPHA = 1.0
      CALL PSTRSM_I8( 'L', 'U', 'N', 'N', N8, N8, ALPHA,
     $     A, 1_8, 1_8, DESCA8, B, 1_8, 1_8, DESCB8 )
      CALL PSTRSM( 'L', 'U', 'N', 'N', N4, N4, ALPHA,
     $     A, 1, 1, DESCA4, BREF, 1, 1, DESCB4 )
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC
         DO I8 = 1, LR
            IF( B((J8-1)*LLD8+I8) .NE.
     $          BREF((J8-1)*LLD8+I8) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PSTRSM_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PSTRSM_I8: FAILED ', ERRS
         END IF
      END IF
      DEALLOCATE( A, B, BREF )
      END
*
*     ================================================================
*     TEST_CTRSM — PCTRSM_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_CTRSM( N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
      INTEGER*8          LR, LC, LLD8, I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 )
      INTEGER            DESCA4( 9 ), DESCB4( 9 ), INFO, ERRS, N4, NB4
      COMPLEX, ALLOCATABLE :: A(:), B(:), BREF(:)
      COMPLEX            ALPHA
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC, DESCINIT_I8, DESCINIT,
     $                   PCTRSM_I8, PCTRSM
      INTRINSIC          CMPLX, REAL, MAX, INT, ABS
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCB8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCB4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      ALLOCATE( A( MAX(LLD8*LC,1_8) ), B( MAX(LLD8*LC,1_8) ),
     $          BREF( MAX(LLD8*LC,1_8) ) )
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) = CMPLX( REAL(N8), 0.0 )
            ELSE IF( GI .LT. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) = CMPLX(
     $              1.0/REAL(1+ABS(GI-GJ)), 0.1 )
            ELSE
               A( (J8-1)*LLD8 + I8 ) = CMPLX( 0.0, 0.0 )
            END IF
            B( (J8-1)*LLD8 + I8 ) = CMPLX( 1.0, 0.0 )
            BREF( (J8-1)*LLD8 + I8 ) = CMPLX( 1.0, 0.0 )
         END DO
      END DO
      ALPHA = CMPLX( 1.0, 0.0 )
      CALL PCTRSM_I8( 'L', 'U', 'N', 'N', N8, N8, ALPHA,
     $     A, 1_8, 1_8, DESCA8, B, 1_8, 1_8, DESCB8 )
      CALL PCTRSM( 'L', 'U', 'N', 'N', N4, N4, ALPHA,
     $     A, 1, 1, DESCA4, BREF, 1, 1, DESCB4 )
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC
         DO I8 = 1, LR
            IF( B((J8-1)*LLD8+I8) .NE.
     $          BREF((J8-1)*LLD8+I8) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PCTRSM_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PCTRSM_I8: FAILED ', ERRS
         END IF
      END IF
      DEALLOCATE( A, B, BREF )
      END
*
*     ================================================================
*     TEST_CHERK — PCHERK_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_CHERK( N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
      INTEGER*8          LR, LC, LLD8, K8, I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCC8( 9 )
      INTEGER            DESCA4( 9 ), DESCC4( 9 )
      INTEGER            INFO, ERRS, N4, NB4, K4
      COMPLEX, ALLOCATABLE :: A(:), C(:), CREF(:)
      REAL               ALPHA, BETA
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC, DESCINIT_I8, DESCINIT,
     $                   PCHERK_I8, PCHERK
      INTRINSIC          CMPLX, REAL, MAX, INT, ABS
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      K8 = NB8
      K4 = NB4
      LR = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
      INFO = 0
      CALL DESCINIT_I8( DESCC8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCA8, N8, K8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCC4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCA4, N4, K4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      ALLOCATE( A( MAX(LLD8*LC,1_8) ), C( MAX(LLD8*LC,1_8) ),
     $          CREF( MAX(LLD8*LC,1_8) ) )
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            A( (J8-1)*LLD8 + I8 ) = CMPLX( REAL(GI+GJ),
     $                                       REAL(GI-GJ) )
            IF( GI .EQ. GJ ) THEN
               C( (J8-1)*LLD8 + I8 ) = CMPLX( REAL(2*N8), 0.0 )
            ELSE
               C( (J8-1)*LLD8 + I8 ) = CMPLX(
     $              1.0/REAL(1+ABS(GI-GJ)),
     $              REAL(GI-GJ)/REAL(N8) )
            END IF
            CREF( (J8-1)*LLD8 + I8 ) = C( (J8-1)*LLD8 + I8 )
         END DO
      END DO
      ALPHA = 1.0
      BETA = 1.0
      CALL PCHERK_I8( 'L', 'N', N8, K8, ALPHA,
     $     A, 1_8, 1_8, DESCA8, BETA,
     $     C, 1_8, 1_8, DESCC8 )
      CALL PCHERK( 'L', 'N', N4, K4, ALPHA,
     $     A, 1, 1, DESCA4, BETA,
     $     CREF, 1, 1, DESCC4 )
      NTEST = NTEST + 1
      ERRS = 0
      DO J8 = 1, LC
         DO I8 = 1, LR
            IF( C((J8-1)*LLD8+I8) .NE.
     $          CREF((J8-1)*LLD8+I8) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PCHERK_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PCHERK_I8: FAILED ', ERRS
         END IF
      END IF
      DEALLOCATE( A, C, CREF )
      END
*
*     ================================================================
*     TEST_CGEMM — PCGEMM_I8 vs legacy
*     ================================================================
*
      SUBROUTINE TEST_CGEMM( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
      INTEGER*8          LRM, LCM, LRN, LCN, LLDM, LLDN, K8, I8
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 ), DESCC8( 9 )
      INTEGER            DESCA4( 9 ), DESCB4( 9 ), DESCC4( 9 )
      INTEGER            INFO, ERRS, M4, N4, NB4, K4
      COMPLEX, ALLOCATABLE :: A(:), B(:), C(:), CREF(:)
      COMPLEX            ALPHA, BETA
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC, DESCINIT_I8, DESCINIT,
     $                   PCGEMM_I8, PCGEMM
      INTRINSIC          CMPLX, MAX, INT
      M4 = INT( M8 )
      N4 = INT( N8 )
      NB4 = INT( NB8 )
      K8 = NB8
      K4 = NB4
      LRM = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCM = NUMROC_I8( M8, NB8, MYCOL, 0, NPCOL )
      LRN = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LCN = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDM = MAX( LRM, 1_8 )
      LLDN = MAX( LRN, 1_8 )
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, K8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDM, INFO )
      CALL DESCINIT_I8( DESCB8, K8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  MAX(NUMROC_I8(K8,NB8,MYROW,0,NPROW),1_8),
     $                  INFO )
      CALL DESCINIT_I8( DESCC8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDM, INFO )
      CALL DESCINIT( DESCA4, M4, K4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDM ), INFO )
      CALL DESCINIT( DESCB4, K4, N4, NB4, NB4, 0, 0, ICTXT,
     $               MAX(NUMROC(K4,NB4,MYROW,0,NPROW),1), INFO )
      CALL DESCINIT( DESCC4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDM ), INFO )
      ALLOCATE( A( MAX(LLDM*LCM,1_8) ), B( MAX(LLDN*LCN,1_8) ),
     $          C( MAX(LLDM*LCN,1_8) ), CREF( MAX(LLDM*LCN,1_8) ) )
      DO I8 = 1, MAX( LLDM*LCM, 1_8 )
         A( I8 ) = CMPLX( 1.0, 0.0 )
      END DO
      DO I8 = 1, MAX( LLDN*LCN, 1_8 )
         B( I8 ) = CMPLX( 2.0, 0.0 )
      END DO
      DO I8 = 1, MAX( LLDM*LCN, 1_8 )
         C( I8 ) = CMPLX( 0.0, 0.0 )
         CREF( I8 ) = CMPLX( 0.0, 0.0 )
      END DO
      ALPHA = CMPLX( 1.0, 0.0 )
      BETA = CMPLX( 0.0, 0.0 )
      CALL PCGEMM_I8( 'N', 'N', M8, N8, K8, ALPHA,
     $     A, 1_8, 1_8, DESCA8,
     $     B, 1_8, 1_8, DESCB8, BETA,
     $     C, 1_8, 1_8, DESCC8 )
      CALL PCGEMM( 'N', 'N', M4, N4, K4, ALPHA,
     $     A, 1, 1, DESCA4,
     $     B, 1, 1, DESCB4, BETA,
     $     CREF, 1, 1, DESCC4 )
      NTEST = NTEST + 1
      ERRS = 0
      DO I8 = 1, MAX( LLDM*LCN, 1_8 )
         IF( C( I8 ) .NE. CREF( I8 ) ) ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A)') '  PCGEMM_I8: PASSED'
         ELSE
            WRITE(*,'(A,I6)') '  PCGEMM_I8: FAILED ', ERRS
         END IF
      END IF
      DEALLOCATE( A, B, C, CREF )
      END
