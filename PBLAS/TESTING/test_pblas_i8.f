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
