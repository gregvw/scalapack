      PROGRAM TEST_PGESV_I8
      IMPLICIT NONE
*
*  Test PxGESV_I8 — compare against legacy PxGESV.
*  Verifies GETRF_I8 + GETRS_I8 produce bit-identical solutions.
*
*  Usage:  mpirun -np 4 ./xgesv_i8
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
     $   WRITE(*,'(A)') 'test_pgesv_i8: 2x2 grid'
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PDGESV_I8 (D) ---'
      CALL RUN_DGESV( 30_8, 4_8, 5_8, ICTXT, NPROW, NPCOL,
     $                MYROW, MYCOL, IAM, 'D:N=30,NB=4 ',
     $                NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PSGESV_I8 (S) ---'
      CALL RUN_SGESV( 30_8, 4_8, 5_8, ICTXT, NPROW, NPCOL,
     $                MYROW, MYCOL, IAM, 'S:N=30,NB=4 ',
     $                NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PCGESV_I8 (C) ---'
      CALL RUN_CGESV( 30_8, 4_8, 5_8, ICTXT, NPROW, NPCOL,
     $                MYROW, MYCOL, IAM, 'C:N=30,NB=4 ',
     $                NTEST, NFAIL )
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- PZGESV_I8 (Z) ---'
      CALL RUN_ZGESV( 30_8, 4_8, 5_8, ICTXT, NPROW, NPCOL,
     $                MYROW, MYCOL, IAM, 'Z:N=30,NB=4 ',
     $                NTEST, NFAIL )
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NFAIL_G, NFAIL_G, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PxGESV_I8 Test Summary'
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
*     RUN_DGESV
*     ================================================================
*
      SUBROUTINE RUN_DGESV( N8, NB8, NRHS8, ICTXT, NPROW, NPCOL,
     $                       MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8, NRHS8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LRB, LCB, LLDA, LLDB
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 )
      INTEGER            DESCA4( 9 ), DESCB4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, NRHS4
      DOUBLE PRECISION, ALLOCATABLE :: A(:), ACOPY(:),
     $                   B(:), BCOPY(:), BREF(:)
      INTEGER, ALLOCATABLE :: IPIV(:), IPREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PDGESV_I8, PDGESV
      INTRINSIC          DBLE, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      NRHS4 = INT( NRHS8 )
      LRA  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LRB  = LRA
      LCB  = NUMROC_I8( NRHS8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
      LLDB = MAX( LRB, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCB8, N8, NRHS8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDB, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCB4, N4, NRHS4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDB ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( B( MAX( LLDB*LCB, 1_8 ) ) )
      ALLOCATE( BCOPY( MAX( LLDB*LCB, 1_8 ) ) )
      ALLOCATE( BREF( MAX( LLDB*LCB, 1_8 ) ) )
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
      DO J8 = 1, LCB
         DO I8 = 1, LRB
            B( (J8-1)*LLDB + I8 ) = 1.0D0
         END DO
      END DO
*
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
      DO I8 = 1, MAX( LLDB*LCB, 1_8 )
         BCOPY( I8 ) = B( I8 )
         BREF( I8 ) = B( I8 )
      END DO
*
      INFO = 0
      CALL PDGESV_I8( N8, NRHS8, A, 1_8, 1_8, DESCA8, IPIV,
     $                B, 1_8, 1_8, DESCB8, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, B, BCOPY, BREF, IPIV, IPREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PDGESV( N4, NRHS4, ACOPY, 1, 1, DESCA4, IPREF,
     $             BREF, 1, 1, DESCB4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCB
         DO I8 = 1, LRB
            IF( B( (J8-1)*LLDB+I8 ) .NE.
     $          BREF( (J8-1)*LLDB+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, B, BCOPY, BREF, IPIV, IPREF )
      END
*
*     ================================================================
*     RUN_ZGESV
*     ================================================================
*
      SUBROUTINE RUN_ZGESV( N8, NB8, NRHS8, ICTXT, NPROW, NPCOL,
     $                       MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8, NRHS8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LRB, LCB, LLDA, LLDB
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 )
      INTEGER            DESCA4( 9 ), DESCB4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, NRHS4
      DOUBLE COMPLEX, ALLOCATABLE :: A(:), ACOPY(:),
     $                               B(:), BCOPY(:), BREF(:)
      INTEGER, ALLOCATABLE :: IPIV(:), IPREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT,
     $                   PZGESV_I8, PZGESV
      INTRINSIC          DBLE, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      NRHS4 = INT( NRHS8 )
      LRA  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LRB  = LRA
      LCB  = NUMROC_I8( NRHS8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
      LLDB = MAX( LRB, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCB8, N8, NRHS8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDB, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCB4, N4, NRHS4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDB ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( B( MAX( LLDB*LCB, 1_8 ) ) )
      ALLOCATE( BCOPY( MAX( LLDB*LCB, 1_8 ) ) )
      ALLOCATE( BREF( MAX( LLDB*LCB, 1_8 ) ) )
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
      DO J8 = 1, LCB
         DO I8 = 1, LRB
            B( (J8-1)*LLDB + I8 ) = DCMPLX( 1.0D0, 0.0D0 )
         END DO
      END DO
*
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
      DO I8 = 1, MAX( LLDB*LCB, 1_8 )
         BCOPY( I8 ) = B( I8 )
         BREF( I8 ) = B( I8 )
      END DO
*
      INFO = 0
      CALL PZGESV_I8( N8, NRHS8, A, 1_8, 1_8, DESCA8, IPIV,
     $                B, 1_8, 1_8, DESCB8, INFO )
*
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, B, BCOPY, BREF, IPIV, IPREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PZGESV( N4, NRHS4, ACOPY, 1, 1, DESCA4, IPREF,
     $             BREF, 1, 1, DESCB4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCB
         DO I8 = 1, LRB
            IF( B( (J8-1)*LLDB+I8 ) .NE.
     $          BREF( (J8-1)*LLDB+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, B, BCOPY, BREF, IPIV, IPREF )
      END
*
*     ================================================================
*     RUN_SGESV
*     ================================================================
*
      SUBROUTINE RUN_SGESV( N8, NB8, NRHS8, ICTXT, NPROW, NPCOL,
     $                       MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8, NRHS8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LRB, LCB, LLDA, LLDB
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 )
      INTEGER            DESCA4( 9 ), DESCB4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, NRHS4
      REAL, ALLOCATABLE :: A(:), ACOPY(:),
     $                   B(:), BCOPY(:), BREF(:)
      INTEGER, ALLOCATABLE :: IPIV(:), IPREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT, PSGESV_I8, PSGESV
      INTRINSIC          REAL, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      NRHS4 = INT( NRHS8 )
      LRA  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LRB  = LRA
      LCB  = NUMROC_I8( NRHS8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
      LLDB = MAX( LRB, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCB8, N8, NRHS8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDB, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCB4, N4, NRHS4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDB ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( B( MAX( LLDB*LCB, 1_8 ) ) )
      ALLOCATE( BCOPY( MAX( LLDB*LCB, 1_8 ) ) )
      ALLOCATE( BREF( MAX( LLDB*LCB, 1_8 ) ) )
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
      DO J8 = 1, LCB
         DO I8 = 1, LRB
            B( (J8-1)*LLDB + I8 ) = 1.0
         END DO
      END DO
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
      DO I8 = 1, MAX( LLDB*LCB, 1_8 )
         BCOPY( I8 ) = B( I8 )
         BREF( I8 ) = B( I8 )
      END DO
*
      INFO = 0
      CALL PSGESV_I8( N8, NRHS8, A, 1_8, 1_8, DESCA8, IPIV,
     $                B, 1_8, 1_8, DESCB8, INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, B, BCOPY, BREF, IPIV, IPREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PSGESV( N4, NRHS4, ACOPY, 1, 1, DESCA4, IPREF,
     $             BREF, 1, 1, DESCB4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCB
         DO I8 = 1, LRB
            IF( B( (J8-1)*LLDB+I8 ) .NE.
     $          BREF( (J8-1)*LLDB+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, B, BCOPY, BREF, IPIV, IPREF )
      END
*
*     ================================================================
*     RUN_CGESV
*     ================================================================
*
      SUBROUTINE RUN_CGESV( N8, NB8, NRHS8, ICTXT, NPROW, NPCOL,
     $                       MYROW, MYCOL, IAM, LABEL, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          N8, NB8, NRHS8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      CHARACTER*(*)      LABEL
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LRB, LCB, LLDA, LLDB
      INTEGER*8          I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 ), DESCB8( 9 )
      INTEGER            DESCA4( 9 ), DESCB4( 9 ), INFO, ERRS
      INTEGER            NB4, N4, NRHS4
      COMPLEX, ALLOCATABLE :: A(:), ACOPY(:),
     $                        B(:), BCOPY(:), BREF(:)
      INTEGER, ALLOCATABLE :: IPIV(:), IPREF(:)
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      EXTERNAL           DESCINIT_I8, DESCINIT, PCGESV_I8, PCGESV
      INTRINSIC          CMPLX, REAL, MAX, INT, ABS
*
      NB4 = INT( NB8 )
      N4  = INT( N8 )
      NRHS4 = INT( NRHS8 )
      LRA  = NUMROC_I8( N8, NB8, MYROW, 0, NPROW )
      LCA  = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LRB  = LRA
      LCB  = NUMROC_I8( NRHS8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
      LLDB = MAX( LRB, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, N8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT_I8( DESCB8, N8, NRHS8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDB, INFO )
      CALL DESCINIT( DESCA4, N4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
      CALL DESCINIT( DESCB4, N4, NRHS4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDB ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( ACOPY( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( B( MAX( LLDB*LCB, 1_8 ) ) )
      ALLOCATE( BCOPY( MAX( LLDB*LCB, 1_8 ) ) )
      ALLOCATE( BREF( MAX( LLDB*LCB, 1_8 ) ) )
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
      DO J8 = 1, LCB
         DO I8 = 1, LRB
            B( (J8-1)*LLDB + I8 ) = CMPLX( 1.0, 0.0 )
         END DO
      END DO
      DO I8 = 1, MAX( LLDA*LCA, 1_8 )
         ACOPY( I8 ) = A( I8 )
      END DO
      DO I8 = 1, MAX( LLDB*LCB, 1_8 )
         BCOPY( I8 ) = B( I8 )
         BREF( I8 ) = B( I8 )
      END DO
*
      INFO = 0
      CALL PCGESV_I8( N8, NRHS8, A, 1_8, 1_8, DESCA8, IPIV,
     $                B, 1_8, 1_8, DESCB8, INFO )
      NTEST = NTEST + 1
      IF( INFO .NE. 0 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,*) 'FAIL ', LABEL, ': I8 INFO =', INFO
         NFAIL = NFAIL + 1
         DEALLOCATE( A, ACOPY, B, BCOPY, BREF, IPIV, IPREF )
         RETURN
      END IF
*
      INFO = 0
      CALL PCGESV( N4, NRHS4, ACOPY, 1, 1, DESCA4, IPREF,
     $             BREF, 1, 1, DESCB4, INFO )
*
      ERRS = 0
      NTEST = NTEST + 1
      DO J8 = 1, LCB
         DO I8 = 1, LRB
            IF( B( (J8-1)*LLDB+I8 ) .NE.
     $          BREF( (J8-1)*LLDB+I8 ) ) ERRS = ERRS + 1
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
      DEALLOCATE( A, ACOPY, B, BCOPY, BREF, IPIV, IPREF )
      END
