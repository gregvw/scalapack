      PROGRAM TEST_PCHENTRD_I8
      IMPLICIT NONE
*
*  Test for PCHENTRD_I8 / PZHENTRD_I8 — compare against legacy.
*
*  Usage:  mpirun -np 4 ./xchentrd_i8
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
     $   WRITE(*,'(A)') 'test_pchentrd_i8: 2x2 grid'
*
*     === COMPLEX (PCHENTRD_I8) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PCHENTRD_I8 (C) ---'
      CALL RUN_C( 30_8, 4_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $            IAM, 'C:N=30,NB=4 ', NTEST, NFAIL )
      CALL RUN_C( 50_8, 8_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $            IAM, 'C:N=50,NB=8 ', NTEST, NFAIL )
*
*     === DOUBLE COMPLEX (PZHENTRD_I8) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PZHENTRD_I8 (Z) ---'
      CALL RUN_Z( 30_8, 4_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $            IAM, 'Z:N=30,NB=4 ', NTEST, NFAIL )
      CALL RUN_Z( 50_8, 8_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $            IAM, 'Z:N=50,NB=8 ', NTEST, NFAIL )
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NTEST, NTEST, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PxHENTRD_I8 Test Summary'
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
*     RUN_C — compare PCHENTRD_I8 vs legacy PCHENTRD
*     ================================================================
*
      SUBROUTINE RUN_C( N8, NB8, ICTXT, NPROW, NPCOL, MYROW,
     $                   MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LCE, LLD8, I8, J8, GI, GJ
      INTEGER*8          LWORK8, LRWORK8
      INTEGER*8          DESCA8( 9 )
      INTEGER            DESCA4( 9 ), INFO, ERRS, NB4, N4
      INTEGER            LWORK4, LRWORK4
      COMPLEX, ALLOCATABLE :: A(:), ACOPY(:), TAU(:), TAUREF(:),
     $                        WORK(:), WREF(:)
      REAL, ALLOCATABLE :: D(:), E(:), DREF(:), EREF(:),
     $                     RWORK(:), RWREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PCHENTRD_I8, PCHENTRD
      INTRINSIC          CMPLX, REAL, MAX, INT
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      LR  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LCE = NUMROC_I8( N8-1, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( D( N8 ), E( N8 ), TAU( N8 ) )
      ALLOCATE( DREF( N8 ), EREF( N8 ), TAUREF( N8 ) )
*
*     Hermitian matrix: A(i,j) = (N-|i-j|) + i*(i-j), diag real
*
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) = CMPLX( REAL( 2*N8 ), 0.0 )
            ELSE
               A( (J8-1)*LLD8 + I8 ) = CMPLX(
     $              REAL( N8 - ABS(GI-GJ) ),
     $              REAL( GI - GJ ) )
            END IF
         END DO
      END DO
      DO I8 = 1, MAX( LLD8*LC, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
*     === PCHENTRD_I8 ===
*
      LWORK8 = -1
      LRWORK8 = -1
      ALLOCATE( WORK( 1 ), RWORK( 1 ) )
      CALL PCHENTRD_I8( 'L', N8, A, 1_8, 1_8, DESCA8, D, E, TAU,
     $                   WORK, LWORK8, RWORK, LRWORK8, INFO )
      LWORK8 = INT( REAL( WORK( 1 ) ), 8 )
      LRWORK8 = INT( RWORK( 1 ), 8 )
      DEALLOCATE( WORK, RWORK )
      ALLOCATE( WORK( LWORK8 ), RWORK( MAX( LRWORK8, 1_8 ) ) )
*
      INFO = 0
      CALL PCHENTRD_I8( 'L', N8, A, 1_8, 1_8, DESCA8, D, E, TAU,
     $                   WORK, LWORK8, RWORK, LRWORK8, INFO )
      DEALLOCATE( WORK, RWORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
      END IF
*
*     === Legacy PCHENTRD ===
*
      LWORK4 = -1
      LRWORK4 = -1
      ALLOCATE( WREF( 1 ), RWREF( 1 ) )
      CALL PCHENTRD( 'L', N4, ACOPY, 1, 1, DESCA4, DREF, EREF,
     $               TAUREF, WREF, LWORK4, RWREF, LRWORK4, INFO )
      LWORK4 = INT( REAL( WREF( 1 ) ) )
      LRWORK4 = INT( RWREF( 1 ) )
      DEALLOCATE( WREF, RWREF )
      ALLOCATE( WREF( LWORK4 ), RWREF( MAX( LRWORK4, 1 ) ) )
*
      INFO = 0
      CALL PCHENTRD( 'L', N4, ACOPY, 1, 1, DESCA4, DREF, EREF,
     $               TAUREF, WREF, LWORK4, RWREF, LRWORK4, INFO )
      DEALLOCATE( WREF, RWREF )
*
*     === Compare D, E (locally-owned portions) ===
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, LC
         IF( D( I8 ) .NE. DREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' D: i=', I8
            ERRS = ERRS + 1
         END IF
      END DO
      NTEST = NTEST + 1
      DO I8 = 1, LCE
         IF( E( I8 ) .NE. EREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' E: i=', I8
            ERRS = ERRS + 1
         END IF
      END DO
      NTEST = NTEST + 1
      DO I8 = 1, LCE
         IF( TAU( I8 ) .NE. TAUREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' TAU: i=', I8
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
      DEALLOCATE( A, ACOPY, D, E, TAU, DREF, EREF, TAUREF )
      END
*
*     ================================================================
*     RUN_Z — compare PZHENTRD_I8 vs legacy PZHENTRD
*     ================================================================
*
      SUBROUTINE RUN_Z( N8, NB8, ICTXT, NPROW, NPCOL, MYROW,
     $                   MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LCE, LLD8, I8, J8, GI, GJ
      INTEGER*8          LWORK8, LRWORK8
      INTEGER*8          DESCA8( 9 )
      INTEGER            DESCA4( 9 ), INFO, ERRS, NB4, N4
      INTEGER            LWORK4, LRWORK4
      DOUBLE COMPLEX, ALLOCATABLE :: A(:), ACOPY(:), TAU(:),
     $                               TAUREF(:), WORK(:), WREF(:)
      DOUBLE PRECISION, ALLOCATABLE :: D(:), E(:), DREF(:), EREF(:),
     $                                 RWORK(:), RWREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PZHENTRD_I8, PZHENTRD
      INTRINSIC          DBLE, MAX, INT
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      LR  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LCE = NUMROC_I8( N8-1, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( D( N8 ), E( N8 ), TAU( N8 ) )
      ALLOCATE( DREF( N8 ), EREF( N8 ), TAUREF( N8 ) )
*
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) = DCMPLX( DBLE( 2*N8 ), 0.0D0 )
            ELSE
               A( (J8-1)*LLD8 + I8 ) = DCMPLX(
     $              DBLE( N8 - ABS(GI-GJ) ),
     $              DBLE( GI - GJ ) )
            END IF
         END DO
      END DO
      DO I8 = 1, MAX( LLD8*LC, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
*     === PZHENTRD_I8 ===
*
      LWORK8 = -1
      LRWORK8 = -1
      ALLOCATE( WORK( 1 ), RWORK( 1 ) )
      CALL PZHENTRD_I8( 'L', N8, A, 1_8, 1_8, DESCA8, D, E, TAU,
     $                   WORK, LWORK8, RWORK, LRWORK8, INFO )
      LWORK8 = INT( DBLE( WORK( 1 ) ), 8 )
      LRWORK8 = INT( RWORK( 1 ), 8 )
      DEALLOCATE( WORK, RWORK )
      ALLOCATE( WORK( LWORK8 ), RWORK( MAX( LRWORK8, 1_8 ) ) )
*
      INFO = 0
      CALL PZHENTRD_I8( 'L', N8, A, 1_8, 1_8, DESCA8, D, E, TAU,
     $                   WORK, LWORK8, RWORK, LRWORK8, INFO )
      DEALLOCATE( WORK, RWORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
      END IF
*
*     === Legacy PZHENTRD ===
*
      LWORK4 = -1
      LRWORK4 = -1
      ALLOCATE( WREF( 1 ), RWREF( 1 ) )
      CALL PZHENTRD( 'L', N4, ACOPY, 1, 1, DESCA4, DREF, EREF,
     $               TAUREF, WREF, LWORK4, RWREF, LRWORK4, INFO )
      LWORK4 = INT( DBLE( WREF( 1 ) ) )
      LRWORK4 = INT( RWREF( 1 ) )
      DEALLOCATE( WREF, RWREF )
      ALLOCATE( WREF( LWORK4 ), RWREF( MAX( LRWORK4, 1 ) ) )
*
      INFO = 0
      CALL PZHENTRD( 'L', N4, ACOPY, 1, 1, DESCA4, DREF, EREF,
     $               TAUREF, WREF, LWORK4, RWREF, LRWORK4, INFO )
      DEALLOCATE( WREF, RWREF )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, LC
         IF( D( I8 ) .NE. DREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' D: i=', I8
            ERRS = ERRS + 1
         END IF
      END DO
      NTEST = NTEST + 1
      DO I8 = 1, LCE
         IF( E( I8 ) .NE. EREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' E: i=', I8
            ERRS = ERRS + 1
         END IF
      END DO
      NTEST = NTEST + 1
      DO I8 = 1, LCE
         IF( TAU( I8 ) .NE. TAUREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' TAU: i=', I8
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
      DEALLOCATE( A, ACOPY, D, E, TAU, DREF, EREF, TAUREF )
      END
