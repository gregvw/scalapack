      PROGRAM TEST_PSYEV_I8
      IMPLICIT NONE
*
*  Test for PDSYEV_I8 / PSSYEV_I8 / PCHEEV_I8 / PZHEEV_I8
*  — compare eigenvalues against legacy counterparts.
*
*  Usage:  mpirun -np 4 ./xsyev_i8
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
     $   WRITE(*,'(A)') 'test_psyev_i8: 2x2 grid'
*
*     === DOUBLE PRECISION (PDSYEV_I8) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PDSYEV_I8 (D) ---'
      CALL RUN_DSYEV( 30_8, 4_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $                IAM, 'D:N=30,NB=4 ', NTEST, NFAIL )
      CALL RUN_DSYEV( 50_8, 8_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $                IAM, 'D:N=50,NB=8 ', NTEST, NFAIL )
*
*     === REAL (PSSYEV_I8) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PSSYEV_I8 (S) ---'
      CALL RUN_SSYEV( 30_8, 4_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $                IAM, 'S:N=30,NB=4 ', NTEST, NFAIL )
*
*     === COMPLEX (PCHEEV_I8) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PCHEEV_I8 (C) ---'
      CALL RUN_CHEEV( 30_8, 4_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $                IAM, 'C:N=30,NB=4 ', NTEST, NFAIL )
*
*     === DOUBLE COMPLEX (PZHEEV_I8) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PZHEEV_I8 (Z) ---'
      CALL RUN_ZHEEV( 30_8, 4_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $                IAM, 'Z:N=30,NB=4 ', NTEST, NFAIL )
*
*     Summary
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NFAIL_G, NFAIL_G, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PxSYEV/HEEV_I8 Test Summary'
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
*     RUN_DSYEV — compare PDSYEV_I8 eigenvalues vs legacy PDSYEV
*     ================================================================
*
      SUBROUTINE RUN_DSYEV( N8, NB8, ICTXT, NPROW, NPCOL, MYROW,
     $                       MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LLD8, I8, J8, GI, GJ, LWORK8
      INTEGER*8          DESCA8( 9 ), DESCZ8( 9 )
      INTEGER            DESCA4( 9 ), DESCZ4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, LWORK4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), ACOPY(:), Z(:), ZREF(:),
     $                   W(:), WREF(:), WORK(:), WKREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PDSYEV_I8, PDSYEV
      INTRINSIC          DBLE, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      LR  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCZ8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCZ4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( Z( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ZREF( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( W( N8 ), WREF( N8 ) )
*
*     Fill symmetric matrix: A(i,j) = N-|i-j|, diag=2*N
*
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) = 2.0D0 * DBLE( N8 )
            ELSE
               A( (J8-1)*LLD8 + I8 ) = DBLE( N8 - ABS( GI - GJ ) )
            END IF
         END DO
      END DO
      DO I8 = 1, MAX( LLD8*LC, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
*     === PDSYEV_I8 ===
*
      LWORK8 = -1
      ALLOCATE( WORK( 1 ) )
      CALL PDSYEV_I8( 'V', 'L', N8, A, 1_8, 1_8, DESCA8, W,
     $                Z, 1_8, 1_8, DESCZ8, WORK, LWORK8, INFO )
      LWORK8 = INT( WORK( 1 ), 8 )
      DEALLOCATE( WORK )
      ALLOCATE( WORK( LWORK8 ) )
      INFO = 0
      CALL PDSYEV_I8( 'V', 'L', N8, A, 1_8, 1_8, DESCA8, W,
     $                Z, 1_8, 1_8, DESCZ8, WORK, LWORK8, INFO )
      DEALLOCATE( WORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
      END IF
*
*     === Legacy PDSYEV ===
*
      LWORK4 = -1
      ALLOCATE( WKREF( 1 ) )
      CALL PDSYEV( 'V', 'L', N4, ACOPY, 1, 1, DESCA4, WREF,
     $             ZREF, 1, 1, DESCZ4, WKREF, LWORK4, INFO )
      LWORK4 = INT( WKREF( 1 ) )
      DEALLOCATE( WKREF )
      ALLOCATE( WKREF( LWORK4 ) )
      INFO = 0
      CALL PDSYEV( 'V', 'L', N4, ACOPY, 1, 1, DESCA4, WREF,
     $             ZREF, 1, 1, DESCZ4, WKREF, LWORK4, INFO )
      DEALLOCATE( WKREF )
*
*     === Compare eigenvalues W (global, all processes have same copy) ===
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, N8
         IF( W( I8 ) .NE. WREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' W: i=', I8,
     $                     ' I8=', W(I8), ' ref=', WREF(I8)
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
      DEALLOCATE( A, ACOPY, Z, ZREF, W, WREF )
      END
*
*     ================================================================
*     RUN_CHEEV — compare PCHEEV_I8 eigenvalues vs legacy PCHEEV
*     ================================================================
*
      SUBROUTINE RUN_CHEEV( N8, NB8, ICTXT, NPROW, NPCOL, MYROW,
     $                       MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LLD8, I8, J8, GI, GJ
      INTEGER*8          LWORK8, LRWORK8
      INTEGER*8          DESCA8( 9 ), DESCZ8( 9 )
      INTEGER            DESCA4( 9 ), DESCZ4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, LWORK4, LRWORK4
      COMPLEX, ALLOCATABLE :: A(:), ACOPY(:), Z(:), ZREF(:),
     $                        WORK(:), WKREF(:)
      REAL, ALLOCATABLE :: W(:), WREF(:), RWORK(:), RWREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PCHEEV_I8, PCHEEV
      INTRINSIC          CMPLX, REAL, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      LR  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCZ8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCZ4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( Z( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ZREF( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( W( N8 ), WREF( N8 ) )
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
*     === PCHEEV_I8 ===
*
      LWORK8 = -1
      LRWORK8 = -1
      ALLOCATE( WORK( 1 ), RWORK( 1 ) )
      CALL PCHEEV_I8( 'V', 'L', N8, A, 1_8, 1_8, DESCA8, W,
     $                Z, 1_8, 1_8, DESCZ8, WORK, LWORK8,
     $                RWORK, LRWORK8, INFO )
      LWORK8 = INT( REAL( WORK( 1 ) ), 8 )
      LRWORK8 = INT( RWORK( 1 ), 8 )
      DEALLOCATE( WORK, RWORK )
      ALLOCATE( WORK( LWORK8 ), RWORK( MAX( LRWORK8, 1_8 ) ) )
      INFO = 0
      CALL PCHEEV_I8( 'V', 'L', N8, A, 1_8, 1_8, DESCA8, W,
     $                Z, 1_8, 1_8, DESCZ8, WORK, LWORK8,
     $                RWORK, LRWORK8, INFO )
      DEALLOCATE( WORK, RWORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
      END IF
*
*     === Legacy PCHEEV ===
*
      LWORK4 = -1
      LRWORK4 = -1
      ALLOCATE( WKREF( 1 ), RWREF( 1 ) )
      CALL PCHEEV( 'V', 'L', N4, ACOPY, 1, 1, DESCA4, WREF,
     $             ZREF, 1, 1, DESCZ4, WKREF, LWORK4,
     $             RWREF, LRWORK4, INFO )
      LWORK4 = INT( REAL( WKREF( 1 ) ) )
      LRWORK4 = INT( RWREF( 1 ) )
      DEALLOCATE( WKREF, RWREF )
      ALLOCATE( WKREF( LWORK4 ), RWREF( MAX( LRWORK4, 1 ) ) )
      INFO = 0
      CALL PCHEEV( 'V', 'L', N4, ACOPY, 1, 1, DESCA4, WREF,
     $             ZREF, 1, 1, DESCZ4, WKREF, LWORK4,
     $             RWREF, LRWORK4, INFO )
      DEALLOCATE( WKREF, RWREF )
*
*     === Compare eigenvalues W ===
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, N8
         IF( W( I8 ) .NE. WREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' W: i=', I8
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
      DEALLOCATE( A, ACOPY, Z, ZREF, W, WREF )
      END
*
*     ================================================================
*     RUN_SSYEV — compare PSSYEV_I8 eigenvalues vs legacy PSSYEV
*     ================================================================
*
      SUBROUTINE RUN_SSYEV( N8, NB8, ICTXT, NPROW, NPCOL, MYROW,
     $                       MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LLD8, I8, J8, GI, GJ, LWORK8
      INTEGER*8          DESCA8( 9 ), DESCZ8( 9 )
      INTEGER            DESCA4( 9 ), DESCZ4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, LWORK4
      REAL, ALLOCATABLE :: A(:), ACOPY(:), Z(:), ZREF(:),
     $                   W(:), WREF(:), WORK(:), WKREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PSSYEV_I8, PSSYEV
      INTRINSIC          REAL, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      LR  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCZ8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCZ4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( Z( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ZREF( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( W( N8 ), WREF( N8 ) )
*
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) = 2.0 * REAL( N8 )
            ELSE
               A( (J8-1)*LLD8 + I8 ) = REAL( N8 - ABS( GI - GJ ) )
            END IF
         END DO
      END DO
      DO I8 = 1, MAX( LLD8*LC, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
      LWORK8 = -1
      ALLOCATE( WORK( 1 ) )
      CALL PSSYEV_I8( 'V', 'L', N8, A, 1_8, 1_8, DESCA8, W,
     $                Z, 1_8, 1_8, DESCZ8, WORK, LWORK8, INFO )
      LWORK8 = INT( WORK( 1 ), 8 )
      DEALLOCATE( WORK )
      ALLOCATE( WORK( LWORK8 ) )
      INFO = 0
      CALL PSSYEV_I8( 'V', 'L', N8, A, 1_8, 1_8, DESCA8, W,
     $                Z, 1_8, 1_8, DESCZ8, WORK, LWORK8, INFO )
      DEALLOCATE( WORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
      END IF
*
      LWORK4 = -1
      ALLOCATE( WKREF( 1 ) )
      CALL PSSYEV( 'V', 'L', N4, ACOPY, 1, 1, DESCA4, WREF,
     $             ZREF, 1, 1, DESCZ4, WKREF, LWORK4, INFO )
      LWORK4 = INT( WKREF( 1 ) )
      DEALLOCATE( WKREF )
      ALLOCATE( WKREF( LWORK4 ) )
      INFO = 0
      CALL PSSYEV( 'V', 'L', N4, ACOPY, 1, 1, DESCA4, WREF,
     $             ZREF, 1, 1, DESCZ4, WKREF, LWORK4, INFO )
      DEALLOCATE( WKREF )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, N8
         IF( W( I8 ) .NE. WREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' W: i=', I8,
     $                     ' I8=', W(I8), ' ref=', WREF(I8)
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
      DEALLOCATE( A, ACOPY, Z, ZREF, W, WREF )
      END
*
*     ================================================================
*     RUN_ZHEEV — compare PZHEEV_I8 eigenvalues vs legacy PZHEEV
*     ================================================================
*
      SUBROUTINE RUN_ZHEEV( N8, NB8, ICTXT, NPROW, NPCOL, MYROW,
     $                       MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LLD8, I8, J8, GI, GJ
      INTEGER*8          LWORK8, LRWORK8
      INTEGER*8          DESCA8( 9 ), DESCZ8( 9 )
      INTEGER            DESCA4( 9 ), DESCZ4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, LWORK4, LRWORK4
      DOUBLE COMPLEX, ALLOCATABLE :: A(:), ACOPY(:), Z(:), ZREF(:),
     $                               WORK(:), WKREF(:)
      DOUBLE PRECISION, ALLOCATABLE :: W(:), WREF(:),
     $                                 RWORK(:), RWREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PZHEEV_I8, PZHEEV
      INTRINSIC          DBLE, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      LR  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT_I8( DESCZ8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
      CALL DESCINIT( DESCZ4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
      ALLOCATE( A( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( Z( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ZREF( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( W( N8 ), WREF( N8 ) )
*
      DO J8 = 1, LC
         GJ = INDXL2G_I8( J8, NB8, MYCOL, 0, NPCOL )
         DO I8 = 1, LR
            GI = INDXL2G_I8( I8, NB8, MYROW, 0, NPROW )
            IF( GI .EQ. GJ ) THEN
               A( (J8-1)*LLD8 + I8 ) = DCMPLX( DBLE(2*N8), 0.0D0 )
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
      LWORK8 = -1
      LRWORK8 = -1
      ALLOCATE( WORK( 1 ), RWORK( 1 ) )
      CALL PZHEEV_I8( 'V', 'L', N8, A, 1_8, 1_8, DESCA8, W,
     $                Z, 1_8, 1_8, DESCZ8, WORK, LWORK8,
     $                RWORK, LRWORK8, INFO )
      LWORK8 = INT( DBLE( WORK( 1 ) ), 8 )
      LRWORK8 = INT( RWORK( 1 ), 8 )
      DEALLOCATE( WORK, RWORK )
      ALLOCATE( WORK( LWORK8 ), RWORK( MAX( LRWORK8, 1_8 ) ) )
      INFO = 0
      CALL PZHEEV_I8( 'V', 'L', N8, A, 1_8, 1_8, DESCA8, W,
     $                Z, 1_8, 1_8, DESCZ8, WORK, LWORK8,
     $                RWORK, LRWORK8, INFO )
      DEALLOCATE( WORK, RWORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
      END IF
*
      LWORK4 = -1
      LRWORK4 = -1
      ALLOCATE( WKREF( 1 ), RWREF( 1 ) )
      CALL PZHEEV( 'V', 'L', N4, ACOPY, 1, 1, DESCA4, WREF,
     $             ZREF, 1, 1, DESCZ4, WKREF, LWORK4,
     $             RWREF, LRWORK4, INFO )
      LWORK4 = INT( DBLE( WKREF( 1 ) ) )
      LRWORK4 = INT( RWREF( 1 ) )
      DEALLOCATE( WKREF, RWREF )
      ALLOCATE( WKREF( LWORK4 ), RWREF( MAX( LRWORK4, 1 ) ) )
      INFO = 0
      CALL PZHEEV( 'V', 'L', N4, ACOPY, 1, 1, DESCA4, WREF,
     $             ZREF, 1, 1, DESCZ4, WKREF, LWORK4,
     $             RWREF, LRWORK4, INFO )
      DEALLOCATE( WKREF, RWREF )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, N8
         IF( W( I8 ) .NE. WREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' W: i=', I8,
     $                     ' I8=', W(I8), ' ref=', WREF(I8)
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
      DEALLOCATE( A, ACOPY, Z, ZREF, W, WREF )
      END
