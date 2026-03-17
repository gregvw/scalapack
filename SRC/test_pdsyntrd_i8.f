      PROGRAM TEST_PDSYNTRD_I8
      IMPLICIT NONE
*
*  Test for PDSYNTRD_I8 — compare against legacy PDSYNTRD.
*  Both should produce bit-identical D, E, TAU since PDSYNTRD_I8
*  is a bridge that calls the same underlying PBLAS/LAPACK code.
*
*  Usage:  mpirun -np 4 ./xdsyntrd_i8
*
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL
      INTEGER            NPROCS, IAM, INFO
      INTEGER*8          N8, NB8, LLD8, LR, LC
      INTEGER*8          DESCA8( 9 )
      INTEGER            DESCA4( 9 )
      INTEGER*8          I8, J8, GI, GJ, LWORK8
      INTEGER            NFAIL, NTEST, NFAIL_G, LWORK4
*
      DOUBLE PRECISION, ALLOCATABLE :: A(:), ACOPY(:)
      DOUBLE PRECISION, ALLOCATABLE :: D(:), E(:), TAU(:), WORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: DREF(:), EREF(:), TAUREF(:),
     $                                 WREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDEXIT, BLACS_EXIT,
     $                   DESCINIT_I8, DESCINIT, PDSYNTRD_I8,
     $                   PDSYNTRD, IGAMX2D
      INTRINSIC          DBLE, MAX, INT, ABS
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
     $   WRITE(*,'(A)') 'test_pdsyntrd_i8: 2x2 grid'
*
*     === DOUBLE PRECISION (PDSYNTRD_I8) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PDSYNTRD_I8 (D) ---'
      CALL RUN_D( 30_8, 4_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $            IAM, 'D:N=30,NB=4 ', NTEST, NFAIL )
      CALL RUN_D( 50_8, 8_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $            IAM, 'D:N=50,NB=8 ', NTEST, NFAIL )
      CALL RUN_D( 16_8, 16_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $            IAM, 'D:N=16,NB=16', NTEST, NFAIL )
*
*     === REAL (PSSYNTRD_I8) ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PSSYNTRD_I8 (S) ---'
      CALL RUN_S( 30_8, 4_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $            IAM, 'S:N=30,NB=4 ', NTEST, NFAIL )
      CALL RUN_S( 50_8, 8_8, ICTXT, NPROW, NPCOL, MYROW, MYCOL,
     $            IAM, 'S:N=50,NB=8 ', NTEST, NFAIL )
*
*     Summary
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NTEST, NTEST, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PDSYNTRD_I8 Test Summary'
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
*     RUN_D — compare PDSYNTRD_I8 vs legacy PDSYNTRD
*     ================================================================
*
      SUBROUTINE RUN_D( N8, NB8, ICTXT, NPROW, NPCOL, MYROW,
     $                   MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LCE, LLD8, I8, J8, GI, GJ, LWORK8
      INTEGER*8          DESCA8( 9 )
      INTEGER            DESCA4( 9 ), INFO, ERRS, NB4, N4, LWORK4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), ACOPY(:),
     $                   D(:), E(:), TAU(:), WORK(:),
     $                   DREF(:), EREF(:), TAUREF(:), WREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PDSYNTRD_I8, PDSYNTRD
      INTRINSIC          DBLE, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
*
      LR  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LC  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LCE = NUMROC_I8( N8-1, NB8, MYCOL, 0, NPCOL )
      LLD8 = MAX( LR, 1_8 )
*
*     Build I8 descriptor
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLD8, INFO )
*
*     Build matching legacy descriptor
*
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLD8 ), INFO )
*
*     Allocate
*
      ALLOCATE( A(     MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLD8*LC, 1_8 ) ) )
      ALLOCATE( D( N8 ), E( N8 ), TAU( N8 ) )
      ALLOCATE( DREF( N8 ), EREF( N8 ), TAUREF( N8 ) )
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
*
*     Save copy for the reference run
*
      DO I8 = 1, MAX( LLD8*LC, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
*     === Run PDSYNTRD_I8 ===
*
*     Workspace query
      LWORK8 = -1
      ALLOCATE( WORK( 1 ) )
      CALL PDSYNTRD_I8( 'L', N8, A, 1_8, 1_8, DESCA8, D, E, TAU,
     $                   WORK, LWORK8, INFO )
      LWORK8 = INT( WORK( 1 ), 8 )
      DEALLOCATE( WORK )
      ALLOCATE( WORK( LWORK8 ) )
*
*     Actual reduction
      INFO = 0
      CALL PDSYNTRD_I8( 'L', N8, A, 1_8, 1_8, DESCA8, D, E, TAU,
     $                   WORK, LWORK8, INFO )
      DEALLOCATE( WORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
      END IF
*
*     === Run legacy PDSYNTRD on ACOPY ===
*
      LWORK4 = -1
      ALLOCATE( WREF( 1 ) )
      CALL PDSYNTRD( 'L', N4, ACOPY, 1, 1, DESCA4, DREF, EREF,
     $               TAUREF, WREF, LWORK4, INFO )
      LWORK4 = INT( WREF( 1 ) )
      DEALLOCATE( WREF )
      ALLOCATE( WREF( LWORK4 ) )
*
      INFO = 0
      CALL PDSYNTRD( 'L', N4, ACOPY, 1, 1, DESCA4, DREF, EREF,
     $               TAUREF, WREF, LWORK4, INFO )
      DEALLOCATE( WREF )
*
*     === Compare D, E over the locally-owned portion only ===
*     D, E, TAU are distributed 1D vectors with the same column
*     distribution as A.  Only the first LC elements are valid on
*     this process; the rest are uninitialized.
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, LC
         IF( D( I8 ) .NE. DREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' D: i=', I8,
     $                     ' I8=', D(I8), ' ref=', DREF(I8)
            ERRS = ERRS + 1
         END IF
      END DO
*
      NTEST = NTEST + 1
      DO I8 = 1, LCE
         IF( E( I8 ) .NE. EREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' E: i=', I8,
     $                     ' I8=', E(I8), ' ref=', EREF(I8)
            ERRS = ERRS + 1
         END IF
      END DO
*
      NTEST = NTEST + 1
      DO I8 = 1, LCE
         IF( TAU( I8 ) .NE. TAUREF( I8 ) ) THEN
            IF( ERRS .LT. 3 )
     $         WRITE(*,*) 'FAIL ', LABEL, ' TAU: i=', I8,
     $                     ' I8=', TAU(I8), ' ref=', TAUREF(I8)
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
*     RUN_S — compare PSSYNTRD_I8 vs legacy PSSYNTRD
*     ================================================================
*
      SUBROUTINE RUN_S( N8, NB8, ICTXT, NPROW, NPCOL, MYROW,
     $                   MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LR, LC, LCE, LLD8, I8, J8, GI, GJ, LWORK8
      INTEGER*8          DESCA8( 9 )
      INTEGER            DESCA4( 9 ), INFO, ERRS, NB4, N4, LWORK4
      REAL, ALLOCATABLE :: A(:), ACOPY(:),
     $                   D(:), E(:), TAU(:), WORK(:),
     $                   DREF(:), EREF(:), TAUREF(:), WREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PSSYNTRD_I8, PSSYNTRD
      INTRINSIC          REAL, MAX, INT, ABS
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
*     I8 run
      LWORK8 = -1
      ALLOCATE( WORK( 1 ) )
      CALL PSSYNTRD_I8( 'L', N8, A, 1_8, 1_8, DESCA8, D, E, TAU,
     $                   WORK, LWORK8, INFO )
      LWORK8 = INT( WORK( 1 ), 8 )
      DEALLOCATE( WORK )
      ALLOCATE( WORK( LWORK8 ) )
      INFO = 0
      CALL PSSYNTRD_I8( 'L', N8, A, 1_8, 1_8, DESCA8, D, E, TAU,
     $                   WORK, LWORK8, INFO )
      DEALLOCATE( WORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
      END IF
*
*     Legacy run
      LWORK4 = -1
      ALLOCATE( WREF( 1 ) )
      CALL PSSYNTRD( 'L', N4, ACOPY, 1, 1, DESCA4, DREF, EREF,
     $               TAUREF, WREF, LWORK4, INFO )
      LWORK4 = INT( WREF( 1 ) )
      DEALLOCATE( WREF )
      ALLOCATE( WREF( LWORK4 ) )
      INFO = 0
      CALL PSSYNTRD( 'L', N4, ACOPY, 1, 1, DESCA4, DREF, EREF,
     $               TAUREF, WREF, LWORK4, INFO )
      DEALLOCATE( WREF )
*
*     Compare D, E, TAU (locally-owned portions)
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
