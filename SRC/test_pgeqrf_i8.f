      PROGRAM TEST_PGEQRF_I8
      IMPLICIT NONE
*
*  Test PxGEQRF_I8, PxORGQR_I8/PxUNGQR_I8, PxORMQR_I8/PxUNMQR_I8
*  Compare against legacy routines — bit-identical.
*
*  Usage:  mpirun -np 4 ./xgeqrf_i8
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
     $   WRITE(*,'(A)') 'test_pgeqrf_i8: 2x2 grid'
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PDGEQRF_I8 (D) ---'
      CALL RUN_DGEQRF( 40_8, 30_8, 4_8, 5_8, ICTXT, NPROW, NPCOL,
     $              MYROW, MYCOL, IAM, 'D:M=40,N=30   ', NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PSGEQRF_I8 (S) ---'
      CALL RUN_SGEQRF( 40_8, 30_8, 4_8, 5_8, ICTXT, NPROW, NPCOL,
     $              MYROW, MYCOL, IAM, 'S:M=40,N=30   ', NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PZGEQRF_I8 (Z) ---'
      CALL RUN_ZGEQRF( 40_8, 30_8, 4_8, 5_8, ICTXT, NPROW, NPCOL,
     $              MYROW, MYCOL, IAM, 'Z:M=40,N=30   ', NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PCGEQRF_I8 (C) ---'
      CALL RUN_CGEQRF( 40_8, 30_8, 4_8, 5_8, ICTXT, NPROW, NPCOL,
     $              MYROW, MYCOL, IAM, 'C:M=40,N=30   ', NTEST, NFAIL )
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NFAIL_G, NFAIL_G, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PxGEQRF_I8 Test Summary'
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
*     RUN_DGEQRF
*     ================================================================
*
      SUBROUTINE RUN_DGEQRF( M8, N8, NB8, NRHS8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8, NRHS8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LRC, LCC, LLDA, LLDC
      INTEGER*8          I8, J8, GI, GJ, K8, LTAU
      INTEGER*8          DESCA8( 9 ), DESCC8( 9 )
      INTEGER            DESCA4( 9 ), DESCC4( 9 ), INFO, ERRS
      INTEGER            M4, N4, NB4, K4, NRHS4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), ACOPY(:), TAU(:),
     $                   TAUREF(:), Q(:), QREF(:), C(:), CREF(:),
     $                   WORK(:), WREF(:)
      INTEGER            LWORK4, LWREF4
      DOUBLE PRECISION   WQUERY( 1 ), WQREF( 1 )
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PDGEQRF_I8, PDGEQRF,
     $                   PDORGQR_I8, PDORGQR,
     $                   PDORMQR_I8, PDORMQR
      INTRINSIC          DBLE, MAX, MIN, INT, ABS
*
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
      K8  = MIN( M8, N8 )
      K4  = INT( K8 )
      NRHS4 = INT( NRHS8 )
      LRA  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LTAU = NUMROC_I8( K8, NB8, MYCOL, 0, NPCOL )
      LRC  = LRA
      LCC  = NUMROC_I8( NRHS8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
      LLDC = MAX( LRC, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCC8, M8, NRHS8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDC, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCC4, M4, NRHS4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDC ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( TAU( MAX( LTAU, 1_8 ) ) )
      ALLOCATE( TAUREF( MAX( LTAU, 1_8 ) ) )
      ALLOCATE( Q( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( QREF( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( C( MAX( LLDC*LCC, 1_8 ) ) )
      ALLOCATE( CREF( MAX( LLDC*LCC, 1_8 ) ) )
*
*     Fill A(i,j) = 1/(1+|i-j|) + N*delta(i,j)
*
      DO J8 = 1, LCA
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LRA
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLDA + I8 ) = DBLE( N8 )
            ELSE
               A( (J8-1)*LLDA + I8 ) =
     $              1.0D0 / DBLE( 1 + ABS( GI - GJ ) )
            END IF
         END DO
      END DO
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
*     --- Test 1: GEQRF ---
*
*     Workspace query for I8
      CALL PDGEQRF_I8( M8, N8, A, 1_8, 1_8, DESCA8, TAU,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( WQUERY( 1 ) )
*     Workspace query for legacy
      CALL PDGEQRF( M4, N4, ACOPY, 1, 1, DESCA4, TAUREF,
     $              WQREF, -1, INFO )
      LWREF4 = INT( WQREF( 1 ) )
*
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PDGEQRF_I8( M8, N8, A, 1_8, 1_8, DESCA8, TAU,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' GEQRF I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PDGEQRF( M4, N4, ACOPY, 1, 1, DESCA4, TAUREF,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( A( (J8-1)*LLDA+I8 ) .NE.
     $          ACOPY( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      DO I8 = 1, MAX( LTAU, 1_8 )
         IF( TAU( I8 ) .NE. TAUREF( I8 ) ) ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' GEQRF: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' GEQRF: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
*     --- Test 2: ORGQR ---
*
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         Q( I8 )    = A( I8 )
         QREF( I8 ) = ACOPY( I8 )
      END DO
*
      DEALLOCATE( WORK, WREF )
*     Workspace query
      CALL PDORGQR_I8( M8, K8, K8, Q, 1_8, 1_8, DESCA8, TAU,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( WQUERY( 1 ) )
      CALL PDORGQR( M4, K4, K4, QREF, 1, 1, DESCA4, TAUREF,
     $              WQREF, -1, INFO )
      LWREF4 = INT( WQREF( 1 ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PDORGQR_I8( M8, K8, K8, Q, 1_8, 1_8, DESCA8, TAU,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' ORGQR I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PDORGQR( M4, K4, K4, QREF, 1, 1, DESCA4, TAUREF,
     $              WQREF, -1, INFO )
      LWREF4 = INT( WQREF( 1 ) )
      CALL PDORGQR( M4, K4, K4, QREF, 1, 1, DESCA4, TAUREF,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( Q( (J8-1)*LLDA+I8 ) .NE.
     $          QREF( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' ORGQR: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' ORGQR: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
*     --- Test 3: ORMQR ---
*
      DO J8 = 1, LCC
         DO I8 = 1, LRC
            C( (J8-1)*LLDC + I8 )    = 1.0D0
            CREF( (J8-1)*LLDC + I8 ) = 1.0D0
         END DO
      END DO
*
      DEALLOCATE( WORK, WREF )
*     Workspace query
      CALL PDORMQR_I8( 'L', 'T', M8, NRHS8, K8, A, 1_8, 1_8,
     $                 DESCA8, TAU, C, 1_8, 1_8, DESCC8,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( WQUERY( 1 ) )
      CALL PDORMQR( 'L', 'T', M4, NRHS4, K4, ACOPY, 1, 1,
     $              DESCA4, TAUREF, CREF, 1, 1, DESCC4,
     $              WQREF, -1, INFO )
      LWREF4 = INT( WQREF( 1 ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PDORMQR_I8( 'L', 'T', M8, NRHS8, K8, A, 1_8, 1_8,
     $                 DESCA8, TAU, C, 1_8, 1_8, DESCC8,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' ORMQR I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PDORMQR( 'L', 'T', M4, NRHS4, K4, ACOPY, 1, 1,
     $              DESCA4, TAUREF, CREF, 1, 1, DESCC4,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCC
         DO I8 = 1, LRC
            IF( C( (J8-1)*LLDC+I8 ) .NE.
     $          CREF( (J8-1)*LLDC+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' ORMQR: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' ORMQR: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
      DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $            WORK, WREF )
      END
*
*     ================================================================
*     RUN_SGEQRF
*     ================================================================
*
      SUBROUTINE RUN_SGEQRF( M8, N8, NB8, NRHS8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8, NRHS8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LRC, LCC, LLDA, LLDC
      INTEGER*8          I8, J8, GI, GJ, K8, LTAU
      INTEGER*8          DESCA8( 9 ), DESCC8( 9 )
      INTEGER            DESCA4( 9 ), DESCC4( 9 ), INFO, ERRS
      INTEGER            M4, N4, NB4, K4, NRHS4
      REAL, ALLOCATABLE :: A(:), ACOPY(:), TAU(:),
     $                   TAUREF(:), Q(:), QREF(:), C(:), CREF(:),
     $                   WORK(:), WREF(:)
      INTEGER            LWORK4, LWREF4
      REAL               WQUERY( 1 ), WQREF( 1 )
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PSGEQRF_I8, PSGEQRF,
     $                   PSORGQR_I8, PSORGQR,
     $                   PSORMQR_I8, PSORMQR
      INTRINSIC          REAL, MAX, MIN, INT, ABS
*
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
      K8  = MIN( M8, N8 )
      K4  = INT( K8 )
      NRHS4 = INT( NRHS8 )
      LRA  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LTAU = NUMROC_I8( K8, NB8, MYCOL, 0, NPCOL )
      LRC  = LRA
      LCC  = NUMROC_I8( NRHS8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
      LLDC = MAX( LRC, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCC8, M8, NRHS8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDC, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCC4, M4, NRHS4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDC ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( TAU( MAX( LTAU, 1_8 ) ) )
      ALLOCATE( TAUREF( MAX( LTAU, 1_8 ) ) )
      ALLOCATE( Q( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( QREF( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( C( MAX( LLDC*LCC, 1_8 ) ) )
      ALLOCATE( CREF( MAX( LLDC*LCC, 1_8 ) ) )
*
      DO J8 = 1, LCA
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LRA
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLDA + I8 ) = REAL( N8 )
            ELSE
               A( (J8-1)*LLDA + I8 ) =
     $              1.0 / REAL( 1 + ABS( GI - GJ ) )
            END IF
         END DO
      END DO
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
*     --- Test 1: GEQRF ---
*
      CALL PSGEQRF_I8( M8, N8, A, 1_8, 1_8, DESCA8, TAU,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( WQUERY( 1 ) )
      CALL PSGEQRF( M4, N4, ACOPY, 1, 1, DESCA4, TAUREF,
     $              WQREF, -1, INFO )
      LWREF4 = INT( WQREF( 1 ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PSGEQRF_I8( M8, N8, A, 1_8, 1_8, DESCA8, TAU,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' GEQRF I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PSGEQRF( M4, N4, ACOPY, 1, 1, DESCA4, TAUREF,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( A( (J8-1)*LLDA+I8 ) .NE.
     $          ACOPY( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      DO I8 = 1, MAX( LTAU, 1_8 )
         IF( TAU( I8 ) .NE. TAUREF( I8 ) ) ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' GEQRF: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' GEQRF: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
*     --- Test 2: ORGQR ---
*
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         Q( I8 )    = A( I8 )
         QREF( I8 ) = ACOPY( I8 )
      END DO
*
      DEALLOCATE( WORK, WREF )
      CALL PSORGQR_I8( M8, K8, K8, Q, 1_8, 1_8, DESCA8, TAU,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( WQUERY( 1 ) )
      CALL PSORGQR( M4, K4, K4, QREF, 1, 1, DESCA4, TAUREF,
     $              WQREF, -1, INFO )
      LWREF4 = INT( WQREF( 1 ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PSORGQR_I8( M8, K8, K8, Q, 1_8, 1_8, DESCA8, TAU,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' ORGQR I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PSORGQR( M4, K4, K4, QREF, 1, 1, DESCA4, TAUREF,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( Q( (J8-1)*LLDA+I8 ) .NE.
     $          QREF( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' ORGQR: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' ORGQR: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
*     --- Test 3: ORMQR ---
*
      DO J8 = 1, LCC
         DO I8 = 1, LRC
            C( (J8-1)*LLDC + I8 )    = 1.0
            CREF( (J8-1)*LLDC + I8 ) = 1.0
         END DO
      END DO
*
      DEALLOCATE( WORK, WREF )
      CALL PSORMQR_I8( 'L', 'T', M8, NRHS8, K8, A, 1_8, 1_8,
     $                 DESCA8, TAU, C, 1_8, 1_8, DESCC8,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( WQUERY( 1 ) )
      CALL PSORMQR( 'L', 'T', M4, NRHS4, K4, ACOPY, 1, 1,
     $              DESCA4, TAUREF, CREF, 1, 1, DESCC4,
     $              WQREF, -1, INFO )
      LWREF4 = INT( WQREF( 1 ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PSORMQR_I8( 'L', 'T', M8, NRHS8, K8, A, 1_8, 1_8,
     $                 DESCA8, TAU, C, 1_8, 1_8, DESCC8,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' ORMQR I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PSORMQR( 'L', 'T', M4, NRHS4, K4, ACOPY, 1, 1,
     $              DESCA4, TAUREF, CREF, 1, 1, DESCC4,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCC
         DO I8 = 1, LRC
            IF( C( (J8-1)*LLDC+I8 ) .NE.
     $          CREF( (J8-1)*LLDC+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' ORMQR: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' ORMQR: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
      DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $            WORK, WREF )
      END
*
*     ================================================================
*     RUN_ZGEQRF
*     ================================================================
*
      SUBROUTINE RUN_ZGEQRF( M8, N8, NB8, NRHS8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8, NRHS8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LRC, LCC, LLDA, LLDC
      INTEGER*8          I8, J8, GI, GJ, K8, LTAU
      INTEGER*8          DESCA8( 9 ), DESCC8( 9 )
      INTEGER            DESCA4( 9 ), DESCC4( 9 ), INFO, ERRS
      INTEGER            M4, N4, NB4, K4, NRHS4
      DOUBLE COMPLEX, ALLOCATABLE :: A(:), ACOPY(:), TAU(:),
     $                   TAUREF(:), Q(:), QREF(:), C(:), CREF(:),
     $                   WORK(:), WREF(:)
      INTEGER            LWORK4, LWREF4
      DOUBLE COMPLEX     WQUERY( 1 ), WQREF( 1 )
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PZGEQRF_I8, PZGEQRF,
     $                   PZUNGQR_I8, PZUNGQR,
     $                   PZUNMQR_I8, PZUNMQR
      INTRINSIC          DBLE, MAX, MIN, INT, ABS
*
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
      K8  = MIN( M8, N8 )
      K4  = INT( K8 )
      NRHS4 = INT( NRHS8 )
      LRA  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LTAU = NUMROC_I8( K8, NB8, MYCOL, 0, NPCOL )
      LRC  = LRA
      LCC  = NUMROC_I8( NRHS8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
      LLDC = MAX( LRC, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCC8, M8, NRHS8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDC, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCC4, M4, NRHS4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDC ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( TAU( MAX( LTAU, 1_8 ) ) )
      ALLOCATE( TAUREF( MAX( LTAU, 1_8 ) ) )
      ALLOCATE( Q( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( QREF( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( C( MAX( LLDC*LCC, 1_8 ) ) )
      ALLOCATE( CREF( MAX( LLDC*LCC, 1_8 ) ) )
*
      DO J8 = 1, LCA
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LRA
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLDA + I8 ) = DCMPLX( DBLE(N8), 0.0D0 )
            ELSE
               A( (J8-1)*LLDA + I8 ) = DCMPLX(
     $              1.0D0 / DBLE( 1+ABS(GI-GJ) ),
     $              DBLE( GI-GJ ) / DBLE( N8 ) )
            END IF
         END DO
      END DO
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
*     --- Test 1: GEQRF ---
*
      CALL PZGEQRF_I8( M8, N8, A, 1_8, 1_8, DESCA8, TAU,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( DBLE( WQUERY( 1 ) ) )
      CALL PZGEQRF( M4, N4, ACOPY, 1, 1, DESCA4, TAUREF,
     $              WQREF, -1, INFO )
      LWREF4 = INT( DBLE( WQREF( 1 ) ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PZGEQRF_I8( M8, N8, A, 1_8, 1_8, DESCA8, TAU,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' GEQRF I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PZGEQRF( M4, N4, ACOPY, 1, 1, DESCA4, TAUREF,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( A( (J8-1)*LLDA+I8 ) .NE.
     $          ACOPY( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      DO I8 = 1, MAX( LTAU, 1_8 )
         IF( TAU( I8 ) .NE. TAUREF( I8 ) ) ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' GEQRF: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' GEQRF: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
*     --- Test 2: UNGQR ---
*
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         Q( I8 )    = A( I8 )
         QREF( I8 ) = ACOPY( I8 )
      END DO
*
      DEALLOCATE( WORK, WREF )
      CALL PZUNGQR_I8( M8, K8, K8, Q, 1_8, 1_8, DESCA8, TAU,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( DBLE( WQUERY( 1 ) ) )
      CALL PZUNGQR( M4, K4, K4, QREF, 1, 1, DESCA4, TAUREF,
     $              WQREF, -1, INFO )
      LWREF4 = INT( DBLE( WQREF( 1 ) ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PZUNGQR_I8( M8, K8, K8, Q, 1_8, 1_8, DESCA8, TAU,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' UNGQR I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PZUNGQR( M4, K4, K4, QREF, 1, 1, DESCA4, TAUREF,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( Q( (J8-1)*LLDA+I8 ) .NE.
     $          QREF( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' UNGQR: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' UNGQR: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
*     --- Test 3: UNMQR ---
*
      DO J8 = 1, LCC
         DO I8 = 1, LRC
            C( (J8-1)*LLDC + I8 )    = DCMPLX( 1.0D0, 0.0D0 )
            CREF( (J8-1)*LLDC + I8 ) = DCMPLX( 1.0D0, 0.0D0 )
         END DO
      END DO
*
      DEALLOCATE( WORK, WREF )
      CALL PZUNMQR_I8( 'L', 'C', M8, NRHS8, K8, A, 1_8, 1_8,
     $                 DESCA8, TAU, C, 1_8, 1_8, DESCC8,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( DBLE( WQUERY( 1 ) ) )
      CALL PZUNMQR( 'L', 'C', M4, NRHS4, K4, ACOPY, 1, 1,
     $              DESCA4, TAUREF, CREF, 1, 1, DESCC4,
     $              WQREF, -1, INFO )
      LWREF4 = INT( DBLE( WQREF( 1 ) ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PZUNMQR_I8( 'L', 'C', M8, NRHS8, K8, A, 1_8, 1_8,
     $                 DESCA8, TAU, C, 1_8, 1_8, DESCC8,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' UNMQR I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PZUNMQR( 'L', 'C', M4, NRHS4, K4, ACOPY, 1, 1,
     $              DESCA4, TAUREF, CREF, 1, 1, DESCC4,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCC
         DO I8 = 1, LRC
            IF( C( (J8-1)*LLDC+I8 ) .NE.
     $          CREF( (J8-1)*LLDC+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' UNMQR: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' UNMQR: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
      DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $            WORK, WREF )
      END
*
*     ================================================================
*     RUN_CGEQRF
*     ================================================================
*
      SUBROUTINE RUN_CGEQRF( M8, N8, NB8, NRHS8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8, NRHS8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LRC, LCC, LLDA, LLDC
      INTEGER*8          I8, J8, GI, GJ, K8, LTAU
      INTEGER*8          DESCA8( 9 ), DESCC8( 9 )
      INTEGER            DESCA4( 9 ), DESCC4( 9 ), INFO, ERRS
      INTEGER            M4, N4, NB4, K4, NRHS4
      COMPLEX, ALLOCATABLE :: A(:), ACOPY(:), TAU(:),
     $                   TAUREF(:), Q(:), QREF(:), C(:), CREF(:),
     $                   WORK(:), WREF(:)
      INTEGER            LWORK4, LWREF4
      COMPLEX            WQUERY( 1 ), WQREF( 1 )
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PCGEQRF_I8, PCGEQRF,
     $                   PCUNGQR_I8, PCUNGQR,
     $                   PCUNMQR_I8, PCUNMQR
      INTRINSIC          CMPLX, REAL, MAX, MIN, INT, ABS
*
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
      K8  = MIN( M8, N8 )
      K4  = INT( K8 )
      NRHS4 = INT( NRHS8 )
      LRA  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LTAU = NUMROC_I8( K8, NB8, MYCOL, 0, NPCOL )
      LRC  = LRA
      LCC  = NUMROC_I8( NRHS8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
      LLDC = MAX( LRC, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCC8, M8, NRHS8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDC, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCC4, M4, NRHS4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDC ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( TAU( MAX( LTAU, 1_8 ) ) )
      ALLOCATE( TAUREF( MAX( LTAU, 1_8 ) ) )
      ALLOCATE( Q( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( QREF( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( C( MAX( LLDC*LCC, 1_8 ) ) )
      ALLOCATE( CREF( MAX( LLDC*LCC, 1_8 ) ) )
*
      DO J8 = 1, LCA
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LRA
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLDA + I8 ) = CMPLX( REAL(N8), 0.0 )
            ELSE
               A( (J8-1)*LLDA + I8 ) = CMPLX(
     $              1.0 / REAL( 1+ABS(GI-GJ) ),
     $              REAL( GI-GJ ) / REAL( N8 ) )
            END IF
         END DO
      END DO
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
*     --- Test 1: GEQRF ---
*
      CALL PCGEQRF_I8( M8, N8, A, 1_8, 1_8, DESCA8, TAU,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( REAL( WQUERY( 1 ) ) )
      CALL PCGEQRF( M4, N4, ACOPY, 1, 1, DESCA4, TAUREF,
     $              WQREF, -1, INFO )
      LWREF4 = INT( REAL( WQREF( 1 ) ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PCGEQRF_I8( M8, N8, A, 1_8, 1_8, DESCA8, TAU,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' GEQRF I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PCGEQRF( M4, N4, ACOPY, 1, 1, DESCA4, TAUREF,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( A( (J8-1)*LLDA+I8 ) .NE.
     $          ACOPY( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      DO I8 = 1, MAX( LTAU, 1_8 )
         IF( TAU( I8 ) .NE. TAUREF( I8 ) ) ERRS = ERRS + 1
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' GEQRF: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' GEQRF: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
*     --- Test 2: UNGQR ---
*
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         Q( I8 )    = A( I8 )
         QREF( I8 ) = ACOPY( I8 )
      END DO
*
      DEALLOCATE( WORK, WREF )
      CALL PCUNGQR_I8( M8, K8, K8, Q, 1_8, 1_8, DESCA8, TAU,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( REAL( WQUERY( 1 ) ) )
      CALL PCUNGQR( M4, K4, K4, QREF, 1, 1, DESCA4, TAUREF,
     $              WQREF, -1, INFO )
      LWREF4 = INT( REAL( WQREF( 1 ) ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PCUNGQR_I8( M8, K8, K8, Q, 1_8, 1_8, DESCA8, TAU,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' UNGQR I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PCUNGQR( M4, K4, K4, QREF, 1, 1, DESCA4, TAUREF,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( Q( (J8-1)*LLDA+I8 ) .NE.
     $          QREF( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' UNGQR: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' UNGQR: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
*     --- Test 3: UNMQR ---
*
      DO J8 = 1, LCC
         DO I8 = 1, LRC
            C( (J8-1)*LLDC + I8 )    = CMPLX( 1.0, 0.0 )
            CREF( (J8-1)*LLDC + I8 ) = CMPLX( 1.0, 0.0 )
         END DO
      END DO
*
      DEALLOCATE( WORK, WREF )
      CALL PCUNMQR_I8( 'L', 'C', M8, NRHS8, K8, A, 1_8, 1_8,
     $                 DESCA8, TAU, C, 1_8, 1_8, DESCC8,
     $                 WQUERY, -1_8, INFO )
      LWORK4 = INT( REAL( WQUERY( 1 ) ) )
      CALL PCUNMQR( 'L', 'C', M4, NRHS4, K4, ACOPY, 1, 1,
     $              DESCA4, TAUREF, CREF, 1, 1, DESCC4,
     $              WQREF, -1, INFO )
      LWREF4 = INT( REAL( WQREF( 1 ) ) )
      ALLOCATE( WORK( MAX( LWORK4, 1 ) ) )
      ALLOCATE( WREF( MAX( LWREF4, 1 ) ) )
*
      INFO = 0
      CALL PCUNMQR_I8( 'L', 'C', M8, NRHS8, K8, A, 1_8, 1_8,
     $                 DESCA8, TAU, C, 1_8, 1_8, DESCC8,
     $                 WORK, INT( LWORK4, 8 ), INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ' UNMQR I8 INFO=', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $               WORK, WREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PCUNMQR( 'L', 'C', M4, NRHS4, K4, ACOPY, 1, 1,
     $              DESCA4, TAUREF, CREF, 1, 1, DESCC4,
     $              WREF, LWREF4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCC
         DO I8 = 1, LRC
            IF( C( (J8-1)*LLDC+I8 ) .NE.
     $          CREF( (J8-1)*LLDC+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      NFAIL = NFAIL + ERRS
      IF( IAM .EQ. 0 ) THEN
         IF( ERRS .EQ. 0 ) THEN
            WRITE(*,'(A,A,A)') '  ', LABEL, ' UNMQR: PASSED'
         ELSE
            WRITE(*,'(A,A,A,I6,A)') '  ', LABEL, ' UNMQR: FAILED (',
     $                               ERRS, ' errors)'
         END IF
      END IF
*
      DEALLOCATE( A, ACOPY, TAU, TAUREF, Q, QREF, C, CREF,
     $            WORK, WREF )
      END
