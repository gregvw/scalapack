      PROGRAM TEST_PGETRI_I8
      IMPLICIT NONE
*
*  Test PxGETRI_I8 — compare against legacy PxGETRI.
*  Verifies GETRF_I8 + GETRI_I8 produce bit-identical inverses.
*
*  Usage:  mpirun -np 4 ./xgetri_i8
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
     $   WRITE(*,'(A)') 'test_pgetri_i8: 2x2 grid'
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PDGETRI_I8 (D) ---'
      CALL RUN_DGETRI( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, 'D:N=30,NB=4 ',
     $                 NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PSGETRI_I8 (S) ---'
      CALL RUN_SGETRI( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, 'S:N=30,NB=4 ',
     $                 NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PZGETRI_I8 (Z) ---'
      CALL RUN_ZGETRI( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, 'Z:N=30,NB=4 ',
     $                 NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PCGETRI_I8 (C) ---'
      CALL RUN_CGETRI( 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, 'C:N=30,NB=4 ',
     $                 NTEST, NFAIL )
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NFAIL_G, NFAIL_G, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PxGETRI_I8 Test Summary'
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
*     RUN_DGETRI
*     ================================================================
*
      SUBROUTINE RUN_DGETRI( N8, NB8, ICTXT, NPROW, NPCOL,
     $                       MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LLDA
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 )
      INTEGER*8          LWORK8, LIWORK8
      INTEGER            DESCA4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, LWORK4, LIWORK4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), ACOPY(:),
     $                   WORK(:), WREF(:)
      INTEGER, ALLOCATABLE :: IPIV(:), IPREF(:),
     $                        IWORK(:), IWREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PDGETRF_I8, PDGETRI_I8,
     $                   PDGETRF, PDGETRI
      INTRINSIC          DBLE, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      LRA  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( IPIV( INT( LLDA ) + NB4 ) )
      ALLOCATE( IPREF( INT( LLDA ) + NB4 ) )
*
*     General matrix: A(i,j) = 1/(1+|i-j|), diag boosted
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
*
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
*     --- I8 path: factor then invert ---
*
      INFO = 0
      CALL PDGETRF_I8( N8, N8, A, 1_8, 1_8, DESCA8, IPIV, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 GETRF INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, IPIV, IPREF )
         RETURN
      END IF
*
*     Workspace query
*
      LWORK8 = -1_8
      LIWORK8 = -1_8
      ALLOCATE( WORK( 1 ) )
      ALLOCATE( IWORK( 1 ) )
      CALL PDGETRI_I8( N8, A, 1_8, 1_8, DESCA8, IPIV,
     $                 WORK, LWORK8, IWORK, LIWORK8, INFO )
      LWORK8 = INT( WORK( 1 ), 8 )
      LIWORK8 = INT( IWORK( 1 ), 8 )
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( LWORK8 ) )
      ALLOCATE( IWORK( LIWORK8 ) )
*
      INFO = 0
      CALL PDGETRI_I8( N8, A, 1_8, 1_8, DESCA8, IPIV,
     $                 WORK, LWORK8, IWORK, LIWORK8, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 GETRI INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, IPIV, IPREF, WORK, IWORK )
         RETURN
      END IF
*
*     --- Legacy path: factor then invert ---
*
      INFO = 0
      CALL PDGETRF( N4, N4, ACOPY, 1, 1, DESCA4, IPREF, INFO )
*
      LWORK4 = -1
      LIWORK4 = -1
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( 1 ) )
      ALLOCATE( IWORK( 1 ) )
      CALL PDGETRI( N4, ACOPY, 1, 1, DESCA4, IPREF,
     $              WORK, LWORK4, IWORK, LIWORK4, INFO )
      LWORK4 = INT( WORK( 1 ) )
      LIWORK4 = IWORK( 1 )
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( LWORK4 ) )
      ALLOCATE( IWORK( LIWORK4 ) )
*
      INFO = 0
      CALL PDGETRI( N4, ACOPY, 1, 1, DESCA4, IPREF,
     $              WORK, LWORK4, IWORK, LIWORK4, INFO )
*
*     --- Compare ---
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( A( (J8-1)*LLDA+I8 ) .NE.
     $          ACOPY( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, IPIV, IPREF, WORK, IWORK )
      END
*
*     ================================================================
*     RUN_SGETRI
*     ================================================================
*
      SUBROUTINE RUN_SGETRI( N8, NB8, ICTXT, NPROW, NPCOL,
     $                       MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LLDA
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 )
      INTEGER*8          LWORK8, LIWORK8
      INTEGER            DESCA4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, LWORK4, LIWORK4
      REAL, ALLOCATABLE :: A(:), ACOPY(:), WORK(:), WREF(:)
      INTEGER, ALLOCATABLE :: IPIV(:), IPREF(:),
     $                        IWORK(:), IWREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PSGETRF_I8, PSGETRI_I8,
     $                   PSGETRF, PSGETRI
      INTRINSIC          REAL, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      LRA  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( IPIV( INT( LLDA ) + NB4 ) )
      ALLOCATE( IPREF( INT( LLDA ) + NB4 ) )
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
*
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
      INFO = 0
      CALL PSGETRF_I8( N8, N8, A, 1_8, 1_8, DESCA8, IPIV, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 GETRF INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, IPIV, IPREF )
         RETURN
      END IF
*
      LWORK8 = -1_8
      LIWORK8 = -1_8
      ALLOCATE( WORK( 1 ) )
      ALLOCATE( IWORK( 1 ) )
      CALL PSGETRI_I8( N8, A, 1_8, 1_8, DESCA8, IPIV,
     $                 WORK, LWORK8, IWORK, LIWORK8, INFO )
      LWORK8 = INT( WORK( 1 ), 8 )
      LIWORK8 = INT( IWORK( 1 ), 8 )
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( LWORK8 ) )
      ALLOCATE( IWORK( LIWORK8 ) )
*
      INFO = 0
      CALL PSGETRI_I8( N8, A, 1_8, 1_8, DESCA8, IPIV,
     $                 WORK, LWORK8, IWORK, LIWORK8, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 GETRI INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, IPIV, IPREF, WORK, IWORK )
         RETURN
      END IF
*
      INFO = 0
      CALL PSGETRF( N4, N4, ACOPY, 1, 1, DESCA4, IPREF, INFO )
*
      LWORK4 = -1
      LIWORK4 = -1
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( 1 ) )
      ALLOCATE( IWORK( 1 ) )
      CALL PSGETRI( N4, ACOPY, 1, 1, DESCA4, IPREF,
     $              WORK, LWORK4, IWORK, LIWORK4, INFO )
      LWORK4 = INT( WORK( 1 ) )
      LIWORK4 = IWORK( 1 )
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( LWORK4 ) )
      ALLOCATE( IWORK( LIWORK4 ) )
*
      INFO = 0
      CALL PSGETRI( N4, ACOPY, 1, 1, DESCA4, IPREF,
     $              WORK, LWORK4, IWORK, LIWORK4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( A( (J8-1)*LLDA+I8 ) .NE.
     $          ACOPY( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, IPIV, IPREF, WORK, IWORK )
      END
*
*     ================================================================
*     RUN_ZGETRI
*     ================================================================
*
      SUBROUTINE RUN_ZGETRI( N8, NB8, ICTXT, NPROW, NPCOL,
     $                       MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LLDA
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 )
      INTEGER*8          LWORK8, LIWORK8
      INTEGER            DESCA4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, LWORK4, LIWORK4
      DOUBLE COMPLEX, ALLOCATABLE :: A(:), ACOPY(:),
     $                               WORK(:), WREF(:)
      INTEGER, ALLOCATABLE :: IPIV(:), IPREF(:),
     $                        IWORK(:), IWREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PZGETRF_I8, PZGETRI_I8,
     $                   PZGETRF, PZGETRI
      INTRINSIC          DBLE, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      LRA  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( IPIV( INT( LLDA ) + NB4 ) )
      ALLOCATE( IPREF( INT( LLDA ) + NB4 ) )
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
*
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
      INFO = 0
      CALL PZGETRF_I8( N8, N8, A, 1_8, 1_8, DESCA8, IPIV, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 GETRF INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, IPIV, IPREF )
         RETURN
      END IF
*
      LWORK8 = -1_8
      LIWORK8 = -1_8
      ALLOCATE( WORK( 1 ) )
      ALLOCATE( IWORK( 1 ) )
      CALL PZGETRI_I8( N8, A, 1_8, 1_8, DESCA8, IPIV,
     $                 WORK, LWORK8, IWORK, LIWORK8, INFO )
      LWORK8 = INT( WORK( 1 ), 8 )
      LIWORK8 = INT( IWORK( 1 ), 8 )
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( LWORK8 ) )
      ALLOCATE( IWORK( LIWORK8 ) )
*
      INFO = 0
      CALL PZGETRI_I8( N8, A, 1_8, 1_8, DESCA8, IPIV,
     $                 WORK, LWORK8, IWORK, LIWORK8, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 GETRI INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, IPIV, IPREF, WORK, IWORK )
         RETURN
      END IF
*
      INFO = 0
      CALL PZGETRF( N4, N4, ACOPY, 1, 1, DESCA4, IPREF, INFO )
*
      LWORK4 = -1
      LIWORK4 = -1
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( 1 ) )
      ALLOCATE( IWORK( 1 ) )
      CALL PZGETRI( N4, ACOPY, 1, 1, DESCA4, IPREF,
     $              WORK, LWORK4, IWORK, LIWORK4, INFO )
      LWORK4 = INT( DBLE( WORK( 1 ) ) )
      LIWORK4 = IWORK( 1 )
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( LWORK4 ) )
      ALLOCATE( IWORK( LIWORK4 ) )
*
      INFO = 0
      CALL PZGETRI( N4, ACOPY, 1, 1, DESCA4, IPREF,
     $              WORK, LWORK4, IWORK, LIWORK4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( A( (J8-1)*LLDA+I8 ) .NE.
     $          ACOPY( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, IPIV, IPREF, WORK, IWORK )
      END
*
*     ================================================================
*     RUN_CGETRI
*     ================================================================
*
      SUBROUTINE RUN_CGETRI( N8, NB8, ICTXT, NPROW, NPCOL,
     $                       MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LLDA
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 )
      INTEGER*8          LWORK8, LIWORK8
      INTEGER            DESCA4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, LWORK4, LIWORK4
      COMPLEX, ALLOCATABLE :: A(:), ACOPY(:), WORK(:), WREF(:)
      INTEGER, ALLOCATABLE :: IPIV(:), IPREF(:),
     $                        IWORK(:), IWREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PCGETRF_I8, PCGETRI_I8,
     $                   PCGETRF, PCGETRI
      INTRINSIC          CMPLX, REAL, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      LRA  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( IPIV( INT( LLDA ) + NB4 ) )
      ALLOCATE( IPREF( INT( LLDA ) + NB4 ) )
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
*
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
*
      INFO = 0
      CALL PCGETRF_I8( N8, N8, A, 1_8, 1_8, DESCA8, IPIV, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 GETRF INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, IPIV, IPREF )
         RETURN
      END IF
*
      LWORK8 = -1_8
      LIWORK8 = -1_8
      ALLOCATE( WORK( 1 ) )
      ALLOCATE( IWORK( 1 ) )
      CALL PCGETRI_I8( N8, A, 1_8, 1_8, DESCA8, IPIV,
     $                 WORK, LWORK8, IWORK, LIWORK8, INFO )
      LWORK8 = INT( WORK( 1 ), 8 )
      LIWORK8 = INT( IWORK( 1 ), 8 )
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( LWORK8 ) )
      ALLOCATE( IWORK( LIWORK8 ) )
*
      INFO = 0
      CALL PCGETRI_I8( N8, A, 1_8, 1_8, DESCA8, IPIV,
     $                 WORK, LWORK8, IWORK, LIWORK8, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 GETRI INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, IPIV, IPREF, WORK, IWORK )
         RETURN
      END IF
*
      INFO = 0
      CALL PCGETRF( N4, N4, ACOPY, 1, 1, DESCA4, IPREF, INFO )
*
      LWORK4 = -1
      LIWORK4 = -1
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( 1 ) )
      ALLOCATE( IWORK( 1 ) )
      CALL PCGETRI( N4, ACOPY, 1, 1, DESCA4, IPREF,
     $              WORK, LWORK4, IWORK, LIWORK4, INFO )
      LWORK4 = INT( REAL( WORK( 1 ) ) )
      LIWORK4 = IWORK( 1 )
      DEALLOCATE( WORK, IWORK )
      ALLOCATE( WORK( LWORK4 ) )
      ALLOCATE( IWORK( LIWORK4 ) )
*
      INFO = 0
      CALL PCGETRI( N4, ACOPY, 1, 1, DESCA4, IPREF,
     $              WORK, LWORK4, IWORK, LIWORK4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCA
         DO I8 = 1, LRA
            IF( A( (J8-1)*LLDA+I8 ) .NE.
     $          ACOPY( (J8-1)*LLDA+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, IPIV, IPREF, WORK, IWORK )
      END
