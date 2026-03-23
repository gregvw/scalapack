      PROGRAM TEST_PGESVD_I8
      IMPLICIT NONE
*
*  Test PxGESVD_I8 — compare against legacy PxGESVD.
*  Verifies thin SVD produces bit-identical singular values and vectors.
*
*  Usage:  mpirun -np 4 ./xgesvd_i8
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
     $   WRITE(*,'(A)') 'test_pgesvd_i8: 2x2 grid'
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PDGESVD_I8 (D) ---'
      CALL RUN_DGESVD( 40_8, 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, 'D:M=40,N=30   ',
     $                 NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PSGESVD_I8 (S) ---'
      CALL RUN_SGESVD( 40_8, 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, 'S:M=40,N=30   ',
     $                 NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PZGESVD_I8 (Z) ---'
      CALL RUN_ZGESVD( 40_8, 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, 'Z:M=40,N=30   ',
     $                 NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PCGESVD_I8 (C) ---'
      CALL RUN_CGESVD( 40_8, 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                 MYROW, MYCOL, IAM, 'C:M=40,N=30   ',
     $                 NTEST, NFAIL )
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NFAIL_G, NFAIL_G, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PxGESVD_I8 Test Summary'
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
*     RUN_DGESVD
*     ================================================================
*
      SUBROUTINE RUN_DGESVD( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          K8, LRA, LCA, LRU, LCU, LRVT, LCVT
      INTEGER*8          LLDA, LLDU, LLDVT, LWORK8
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCU8( 9 ), DESCVT8( 9 )
      INTEGER            DESCA4( 9 ), DESCU4( 9 ), DESCVT4( 9 )
      INTEGER            INFO, ERRS, K4, LWORK4
      INTEGER            M4, N4, NB4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), ACOPY(:),
     $                   U(:), UREF(:), VT(:), VTREF(:),
     $                   WORK(:), WREF(:)
      DOUBLE PRECISION, ALLOCATABLE :: S(:), SREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PDGESVD_I8, PDGESVD
      INTRINSIC          DBLE, MAX, MIN, INT, ABS
*
      K8  = MIN( M8, N8 )
      K4  = INT( K8 )
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
*
      LRA  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LRU  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCU  = NUMROC_I8( K8, NB8, MYCOL, 0, NPCOL )
      LRVT = NUMROC_I8( K8, NB8, MYROW, 0, NPROW )
      LCVT = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA  = MAX( LRA, 1_8 )
      LLDU  = MAX( LRU, 1_8 )
      LLDVT = MAX( LRVT, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCU8, M8, K8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDU, INFO )
      CALL DESCINIT_I8( DESCVT8, K8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDVT, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCU4, M4, K4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDU ), INFO )
      CALL DESCINIT( DESCVT4, K4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDVT ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( U( MAX( LLDU*LCU, 1_8 ) ) )
      ALLOCATE( UREF( MAX( LLDU*LCU, 1_8 ) ) )
      ALLOCATE( VT( MAX( LLDVT*LCVT, 1_8 ) ) )
      ALLOCATE( VTREF( MAX( LLDVT*LCVT, 1_8 ) ) )
      ALLOCATE( S( K4 ) )
      ALLOCATE( SREF( K4 ) )
*
*     Well-conditioned matrix: A(i,j) = 1/(1+|i-j|) + N*delta(i,j)
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
*     Workspace query
*
      ALLOCATE( WORK( 1 ) )
      LWORK8 = -1_8
      INFO = 0
      CALL PDGESVD_I8( 'V', 'V', M8, N8, A, 1_8, 1_8, DESCA8,
     $                 S, U, 1_8, 1_8, DESCU8,
     $                 VT, 1_8, 1_8, DESCVT8, WORK, LWORK8, INFO )
      LWORK8 = INT( WORK( 1 ), 8 )
      DEALLOCATE( WORK )
      ALLOCATE( WORK( LWORK8 ) )
*
*     I8 SVD
*
      INFO = 0
      CALL PDGESVD_I8( 'V', 'V', M8, N8, A, 1_8, 1_8, DESCA8,
     $                 S, U, 1_8, 1_8, DESCU8,
     $                 VT, 1_8, 1_8, DESCVT8, WORK, LWORK8, INFO )
      DEALLOCATE( WORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, U, UREF, VT, VTREF, S, SREF )
         RETURN
      END IF
*
*     Legacy SVD workspace query
*
      ALLOCATE( WREF( 1 ) )
      LWORK4 = -1
      INFO = 0
      CALL PDGESVD( 'V', 'V', M4, N4, ACOPY, 1, 1, DESCA4,
     $              SREF, UREF, 1, 1, DESCU4,
     $              VTREF, 1, 1, DESCVT4, WREF, LWORK4, INFO )
      LWORK4 = INT( WREF( 1 ) )
      DEALLOCATE( WREF )
      ALLOCATE( WREF( LWORK4 ) )
*
      INFO = 0
      CALL PDGESVD( 'V', 'V', M4, N4, ACOPY, 1, 1, DESCA4,
     $              SREF, UREF, 1, 1, DESCU4,
     $              VTREF, 1, 1, DESCVT4, WREF, LWORK4, INFO )
      DEALLOCATE( WREF )
*
*     Compare singular values
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, K8
         IF( S( I8 ) .NE. SREF( I8 ) ) ERRS = ERRS + 1
      END DO
*
*     Compare U
*
      DO J8 = 1, LCU
         DO I8 = 1, LRU
            IF( U( (J8-1)*LLDU+I8 ) .NE.
     $          UREF( (J8-1)*LLDU+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
*
*     Compare VT
*
      DO J8 = 1, LCVT
         DO I8 = 1, LRVT
            IF( VT( (J8-1)*LLDVT+I8 ) .NE.
     $          VTREF( (J8-1)*LLDVT+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, U, UREF, VT, VTREF, S, SREF )
      END
*
*     ================================================================
*     RUN_SGESVD
*     ================================================================
*
      SUBROUTINE RUN_SGESVD( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          K8, LRA, LCA, LRU, LCU, LRVT, LCVT
      INTEGER*8          LLDA, LLDU, LLDVT, LWORK8
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCU8( 9 ), DESCVT8( 9 )
      INTEGER            DESCA4( 9 ), DESCU4( 9 ), DESCVT4( 9 )
      INTEGER            INFO, ERRS, K4, LWORK4
      INTEGER            M4, N4, NB4
      REAL, ALLOCATABLE :: A(:), ACOPY(:),
     $                   U(:), UREF(:), VT(:), VTREF(:),
     $                   WORK(:), WREF(:)
      REAL, ALLOCATABLE :: S(:), SREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PSGESVD_I8, PSGESVD
      INTRINSIC          REAL, MAX, MIN, INT, ABS
*
      K8  = MIN( M8, N8 )
      K4  = INT( K8 )
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
*
      LRA  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LRU  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCU  = NUMROC_I8( K8, NB8, MYCOL, 0, NPCOL )
      LRVT = NUMROC_I8( K8, NB8, MYROW, 0, NPROW )
      LCVT = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA  = MAX( LRA, 1_8 )
      LLDU  = MAX( LRU, 1_8 )
      LLDVT = MAX( LRVT, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCU8, M8, K8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDU, INFO )
      CALL DESCINIT_I8( DESCVT8, K8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDVT, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCU4, M4, K4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDU ), INFO )
      CALL DESCINIT( DESCVT4, K4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDVT ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( U( MAX( LLDU*LCU, 1_8 ) ) )
      ALLOCATE( UREF( MAX( LLDU*LCU, 1_8 ) ) )
      ALLOCATE( VT( MAX( LLDVT*LCVT, 1_8 ) ) )
      ALLOCATE( VTREF( MAX( LLDVT*LCVT, 1_8 ) ) )
      ALLOCATE( S( K4 ) )
      ALLOCATE( SREF( K4 ) )
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
*     Workspace query
*
      ALLOCATE( WORK( 1 ) )
      LWORK8 = -1_8
      INFO = 0
      CALL PSGESVD_I8( 'V', 'V', M8, N8, A, 1_8, 1_8, DESCA8,
     $                 S, U, 1_8, 1_8, DESCU8,
     $                 VT, 1_8, 1_8, DESCVT8, WORK, LWORK8, INFO )
      LWORK8 = INT( WORK( 1 ), 8 )
      DEALLOCATE( WORK )
      ALLOCATE( WORK( LWORK8 ) )
*
      INFO = 0
      CALL PSGESVD_I8( 'V', 'V', M8, N8, A, 1_8, 1_8, DESCA8,
     $                 S, U, 1_8, 1_8, DESCU8,
     $                 VT, 1_8, 1_8, DESCVT8, WORK, LWORK8, INFO )
      DEALLOCATE( WORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, U, UREF, VT, VTREF, S, SREF )
         RETURN
      END IF
*
*     Legacy SVD
*
      ALLOCATE( WREF( 1 ) )
      LWORK4 = -1
      INFO = 0
      CALL PSGESVD( 'V', 'V', M4, N4, ACOPY, 1, 1, DESCA4,
     $              SREF, UREF, 1, 1, DESCU4,
     $              VTREF, 1, 1, DESCVT4, WREF, LWORK4, INFO )
      LWORK4 = INT( WREF( 1 ) )
      DEALLOCATE( WREF )
      ALLOCATE( WREF( LWORK4 ) )
*
      INFO = 0
      CALL PSGESVD( 'V', 'V', M4, N4, ACOPY, 1, 1, DESCA4,
     $              SREF, UREF, 1, 1, DESCU4,
     $              VTREF, 1, 1, DESCVT4, WREF, LWORK4, INFO )
      DEALLOCATE( WREF )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, K8
         IF( S( I8 ) .NE. SREF( I8 ) ) ERRS = ERRS + 1
      END DO
      DO J8 = 1, LCU
         DO I8 = 1, LRU
            IF( U( (J8-1)*LLDU+I8 ) .NE.
     $          UREF( (J8-1)*LLDU+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      DO J8 = 1, LCVT
         DO I8 = 1, LRVT
            IF( VT( (J8-1)*LLDVT+I8 ) .NE.
     $          VTREF( (J8-1)*LLDVT+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, U, UREF, VT, VTREF, S, SREF )
      END
*
*     ================================================================
*     RUN_ZGESVD
*     ================================================================
*
      SUBROUTINE RUN_ZGESVD( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          K8, LRA, LCA, LRU, LCU, LRVT, LCVT
      INTEGER*8          LLDA, LLDU, LLDVT, LWORK8
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCU8( 9 ), DESCVT8( 9 )
      INTEGER            DESCA4( 9 ), DESCU4( 9 ), DESCVT4( 9 )
      INTEGER            INFO, ERRS, K4, LWORK4
      INTEGER            M4, N4, NB4, LRWORK
      DOUBLE COMPLEX, ALLOCATABLE :: A(:), ACOPY(:),
     $                   U(:), UREF(:), VT(:), VTREF(:),
     $                   WORK(:), WREF(:)
      DOUBLE PRECISION, ALLOCATABLE :: S(:), SREF(:),
     $                   RWORK(:), RWREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PZGESVD_I8, PZGESVD
      INTRINSIC          DBLE, MAX, MIN, INT, ABS
*
      K8  = MIN( M8, N8 )
      K4  = INT( K8 )
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
      LRWORK = 1 + 4*MAX( M4, N4 ) + 10
*
      LRA  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LRU  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCU  = NUMROC_I8( K8, NB8, MYCOL, 0, NPCOL )
      LRVT = NUMROC_I8( K8, NB8, MYROW, 0, NPROW )
      LCVT = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA  = MAX( LRA, 1_8 )
      LLDU  = MAX( LRU, 1_8 )
      LLDVT = MAX( LRVT, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCU8, M8, K8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDU, INFO )
      CALL DESCINIT_I8( DESCVT8, K8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDVT, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCU4, M4, K4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDU ), INFO )
      CALL DESCINIT( DESCVT4, K4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDVT ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( U( MAX( LLDU*LCU, 1_8 ) ) )
      ALLOCATE( UREF( MAX( LLDU*LCU, 1_8 ) ) )
      ALLOCATE( VT( MAX( LLDVT*LCVT, 1_8 ) ) )
      ALLOCATE( VTREF( MAX( LLDVT*LCVT, 1_8 ) ) )
      ALLOCATE( S( K4 ) )
      ALLOCATE( SREF( K4 ) )
      ALLOCATE( RWORK( LRWORK ) )
      ALLOCATE( RWREF( LRWORK ) )
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
*     Workspace query
*
      ALLOCATE( WORK( 1 ) )
      LWORK8 = -1_8
      INFO = 0
      CALL PZGESVD_I8( 'V', 'V', M8, N8, A, 1_8, 1_8, DESCA8,
     $                 S, U, 1_8, 1_8, DESCU8,
     $                 VT, 1_8, 1_8, DESCVT8,
     $                 WORK, LWORK8, RWORK, INFO )
      LWORK8 = INT( DBLE( WORK( 1 ) ), 8 )
      DEALLOCATE( WORK )
      ALLOCATE( WORK( LWORK8 ) )
*
      INFO = 0
      CALL PZGESVD_I8( 'V', 'V', M8, N8, A, 1_8, 1_8, DESCA8,
     $                 S, U, 1_8, 1_8, DESCU8,
     $                 VT, 1_8, 1_8, DESCVT8,
     $                 WORK, LWORK8, RWORK, INFO )
      DEALLOCATE( WORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, U, UREF, VT, VTREF, S, SREF,
     $               RWORK, RWREF )
         RETURN
      END IF
*
*     Legacy SVD
*
      ALLOCATE( WREF( 1 ) )
      LWORK4 = -1
      INFO = 0
      CALL PZGESVD( 'V', 'V', M4, N4, ACOPY, 1, 1, DESCA4,
     $              SREF, UREF, 1, 1, DESCU4,
     $              VTREF, 1, 1, DESCVT4,
     $              WREF, LWORK4, RWREF, INFO )
      LWORK4 = INT( DBLE( WREF( 1 ) ) )
      DEALLOCATE( WREF )
      ALLOCATE( WREF( LWORK4 ) )
*
      INFO = 0
      CALL PZGESVD( 'V', 'V', M4, N4, ACOPY, 1, 1, DESCA4,
     $              SREF, UREF, 1, 1, DESCU4,
     $              VTREF, 1, 1, DESCVT4,
     $              WREF, LWORK4, RWREF, INFO )
      DEALLOCATE( WREF )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, K8
         IF( S( I8 ) .NE. SREF( I8 ) ) ERRS = ERRS + 1
      END DO
      DO J8 = 1, LCU
         DO I8 = 1, LRU
            IF( U( (J8-1)*LLDU+I8 ) .NE.
     $          UREF( (J8-1)*LLDU+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      DO J8 = 1, LCVT
         DO I8 = 1, LRVT
            IF( VT( (J8-1)*LLDVT+I8 ) .NE.
     $          VTREF( (J8-1)*LLDVT+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, U, UREF, VT, VTREF, S, SREF,
     $            RWORK, RWREF )
      END
*
*     ================================================================
*     RUN_CGESVD
*     ================================================================
*
      SUBROUTINE RUN_CGESVD( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                        MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          K8, LRA, LCA, LRU, LCU, LRVT, LCVT
      INTEGER*8          LLDA, LLDU, LLDVT, LWORK8
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCU8( 9 ), DESCVT8( 9 )
      INTEGER            DESCA4( 9 ), DESCU4( 9 ), DESCVT4( 9 )
      INTEGER            INFO, ERRS, K4, LWORK4
      INTEGER            M4, N4, NB4, LRWORK
      COMPLEX, ALLOCATABLE :: A(:), ACOPY(:),
     $                   U(:), UREF(:), VT(:), VTREF(:),
     $                   WORK(:), WREF(:)
      REAL, ALLOCATABLE :: S(:), SREF(:),
     $                   RWORK(:), RWREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PCGESVD_I8, PCGESVD
      INTRINSIC          REAL, MAX, MIN, INT, ABS
*
      K8  = MIN( M8, N8 )
      K4  = INT( K8 )
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
      LRWORK = 1 + 4*MAX( M4, N4 ) + 10
*
      LRA  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LRU  = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCU  = NUMROC_I8( K8, NB8, MYCOL, 0, NPCOL )
      LRVT = NUMROC_I8( K8, NB8, MYROW, 0, NPROW )
      LCVT = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA  = MAX( LRA, 1_8 )
      LLDU  = MAX( LRU, 1_8 )
      LLDVT = MAX( LRVT, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCU8, M8, K8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDU, INFO )
      CALL DESCINIT_I8( DESCVT8, K8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDVT, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCU4, M4, K4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDU ), INFO )
      CALL DESCINIT( DESCVT4, K4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDVT ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( U( MAX( LLDU*LCU, 1_8 ) ) )
      ALLOCATE( UREF( MAX( LLDU*LCU, 1_8 ) ) )
      ALLOCATE( VT( MAX( LLDVT*LCVT, 1_8 ) ) )
      ALLOCATE( VTREF( MAX( LLDVT*LCVT, 1_8 ) ) )
      ALLOCATE( S( K4 ) )
      ALLOCATE( SREF( K4 ) )
      ALLOCATE( RWORK( LRWORK ) )
      ALLOCATE( RWREF( LRWORK ) )
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
*     Workspace query
*
      ALLOCATE( WORK( 1 ) )
      LWORK8 = -1_8
      INFO = 0
      CALL PCGESVD_I8( 'V', 'V', M8, N8, A, 1_8, 1_8, DESCA8,
     $                 S, U, 1_8, 1_8, DESCU8,
     $                 VT, 1_8, 1_8, DESCVT8,
     $                 WORK, LWORK8, RWORK, INFO )
      LWORK8 = INT( REAL( WORK( 1 ) ), 8 )
      DEALLOCATE( WORK )
      ALLOCATE( WORK( LWORK8 ) )
*
      INFO = 0
      CALL PCGESVD_I8( 'V', 'V', M8, N8, A, 1_8, 1_8, DESCA8,
     $                 S, U, 1_8, 1_8, DESCU8,
     $                 VT, 1_8, 1_8, DESCVT8,
     $                 WORK, LWORK8, RWORK, INFO )
      DEALLOCATE( WORK )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, U, UREF, VT, VTREF, S, SREF,
     $               RWORK, RWREF )
         RETURN
      END IF
*
*     Legacy SVD
*
      ALLOCATE( WREF( 1 ) )
      LWORK4 = -1
      INFO = 0
      CALL PCGESVD( 'V', 'V', M4, N4, ACOPY, 1, 1, DESCA4,
     $              SREF, UREF, 1, 1, DESCU4,
     $              VTREF, 1, 1, DESCVT4,
     $              WREF, LWORK4, RWREF, INFO )
      LWORK4 = INT( REAL( WREF( 1 ) ) )
      DEALLOCATE( WREF )
      ALLOCATE( WREF( LWORK4 ) )
*
      INFO = 0
      CALL PCGESVD( 'V', 'V', M4, N4, ACOPY, 1, 1, DESCA4,
     $              SREF, UREF, 1, 1, DESCU4,
     $              VTREF, 1, 1, DESCVT4,
     $              WREF, LWORK4, RWREF, INFO )
      DEALLOCATE( WREF )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO I8 = 1, K8
         IF( S( I8 ) .NE. SREF( I8 ) ) ERRS = ERRS + 1
      END DO
      DO J8 = 1, LCU
         DO I8 = 1, LRU
            IF( U( (J8-1)*LLDU+I8 ) .NE.
     $          UREF( (J8-1)*LLDU+I8 ) ) ERRS = ERRS + 1
         END DO
      END DO
      DO J8 = 1, LCVT
         DO I8 = 1, LRVT
            IF( VT( (J8-1)*LLDVT+I8 ) .NE.
     $          VTREF( (J8-1)*LLDVT+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, U, UREF, VT, VTREF, S, SREF,
     $            RWORK, RWREF )
      END
