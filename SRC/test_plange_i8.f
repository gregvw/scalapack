      PROGRAM TEST_PLANGE_I8
      IMPLICIT NONE
*
*  Test PxLANGE_I8 — compare against legacy PxLANGE.
*  Verifies bit-identical norm values for all 4 types and norms.
*
*  Usage:  mpirun -np 4 ./xlange_i8
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
     $   WRITE(*,'(A)') 'test_plange_i8: 2x2 grid'
*
      CALL TEST_DLANGE( 40_8, 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, IAM, NTEST, NFAIL )
      CALL TEST_SLANGE( 40_8, 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, IAM, NTEST, NFAIL )
      CALL TEST_ZLANGE( 40_8, 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, IAM, NTEST, NFAIL )
      CALL TEST_CLANGE( 40_8, 30_8, 4_8, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, IAM, NTEST, NFAIL )
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NFAIL_G, NFAIL_G, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'PxLANGE_I8 Test Summary'
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
*     TEST_DLANGE
*     ================================================================
*
      SUBROUTINE TEST_DLANGE( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                         MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LLDA, I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 )
      INTEGER            DESCA4( 9 ), M4, N4, NB4, INFO
      DOUBLE PRECISION, ALLOCATABLE :: A(:), WORK(:)
      DOUBLE PRECISION   VAL8, VAL4
      CHARACTER          NORMS( 4 )
      INTEGER            K
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      DOUBLE PRECISION   PDLANGE_I8, PDLANGE
      EXTERNAL           PDLANGE_I8, PDLANGE
      EXTERNAL           DESCINIT_I8, DESCINIT
      INTRINSIC          DBLE, MAX, INT, ABS
*
      DATA NORMS / 'M', '1', 'I', 'F' /
*
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
      LRA = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( WORK( MAX( LRA, LCA, 1_8 ) ) )
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
      DO K = 1, 4
         VAL8 = PDLANGE_I8( NORMS( K ), M8, N8, A, 1_8, 1_8,
     $                       DESCA8, WORK )
         VAL4 = PDLANGE( NORMS( K ), M4, N4, A, 1, 1,
     $                    DESCA4, WORK )
         NTEST = NTEST + 1
         IF( VAL8 .NE. VAL4 ) THEN
            NFAIL = NFAIL + 1
            IF( IAM .EQ. 0 )
     $         WRITE(*,'(A,A,A)') '  PDLANGE_I8 norm=',
     $                            NORMS( K ), ': FAILED'
         ELSE
            IF( IAM .EQ. 0 )
     $         WRITE(*,'(A,A,A)') '  PDLANGE_I8 norm=',
     $                            NORMS( K ), ': PASSED'
         END IF
      END DO
*
      DEALLOCATE( A, WORK )
      END
*
*     ================================================================
*     TEST_SLANGE
*     ================================================================
*
      SUBROUTINE TEST_SLANGE( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                         MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LLDA, I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 )
      INTEGER            DESCA4( 9 ), M4, N4, NB4, INFO
      REAL, ALLOCATABLE :: A(:), WORK(:)
      REAL               VAL8, VAL4
      CHARACTER          NORMS( 4 )
      INTEGER            K
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      REAL               PSLANGE_I8, PSLANGE
      EXTERNAL           PSLANGE_I8, PSLANGE
      EXTERNAL           DESCINIT_I8, DESCINIT
      INTRINSIC          REAL, MAX, INT, ABS
*
      DATA NORMS / 'M', '1', 'I', 'F' /
*
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
      LRA = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( WORK( MAX( LRA, LCA, 1_8 ) ) )
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
      DO K = 1, 4
         VAL8 = PSLANGE_I8( NORMS( K ), M8, N8, A, 1_8, 1_8,
     $                       DESCA8, WORK )
         VAL4 = PSLANGE( NORMS( K ), M4, N4, A, 1, 1,
     $                    DESCA4, WORK )
         NTEST = NTEST + 1
         IF( VAL8 .NE. VAL4 ) THEN
            NFAIL = NFAIL + 1
            IF( IAM .EQ. 0 )
     $         WRITE(*,'(A,A,A)') '  PSLANGE_I8 norm=',
     $                            NORMS( K ), ': FAILED'
         ELSE
            IF( IAM .EQ. 0 )
     $         WRITE(*,'(A,A,A)') '  PSLANGE_I8 norm=',
     $                            NORMS( K ), ': PASSED'
         END IF
      END DO
*
      DEALLOCATE( A, WORK )
      END
*
*     ================================================================
*     TEST_ZLANGE
*     ================================================================
*
      SUBROUTINE TEST_ZLANGE( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                         MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LLDA, I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 )
      INTEGER            DESCA4( 9 ), M4, N4, NB4, INFO
      COMPLEX*16, ALLOCATABLE :: A(:)
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
      DOUBLE PRECISION   VAL8, VAL4
      CHARACTER          NORMS( 4 )
      INTEGER            K
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      DOUBLE PRECISION   PZLANGE_I8, PZLANGE
      EXTERNAL           PZLANGE_I8, PZLANGE
      EXTERNAL           DESCINIT_I8, DESCINIT
      INTRINSIC          DBLE, MAX, INT, ABS
*
      DATA NORMS / 'M', '1', 'I', 'F' /
*
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
      LRA = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( WORK( MAX( LRA, LCA, 1_8 ) ) )
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
      DO K = 1, 4
         VAL8 = PZLANGE_I8( NORMS( K ), M8, N8, A, 1_8, 1_8,
     $                       DESCA8, WORK )
         VAL4 = PZLANGE( NORMS( K ), M4, N4, A, 1, 1,
     $                    DESCA4, WORK )
         NTEST = NTEST + 1
         IF( VAL8 .NE. VAL4 ) THEN
            NFAIL = NFAIL + 1
            IF( IAM .EQ. 0 )
     $         WRITE(*,'(A,A,A)') '  PZLANGE_I8 norm=',
     $                            NORMS( K ), ': FAILED'
         ELSE
            IF( IAM .EQ. 0 )
     $         WRITE(*,'(A,A,A)') '  PZLANGE_I8 norm=',
     $                            NORMS( K ), ': PASSED'
         END IF
      END DO
*
      DEALLOCATE( A, WORK )
      END
*
*     ================================================================
*     TEST_CLANGE
*     ================================================================
*
      SUBROUTINE TEST_CLANGE( M8, N8, NB8, ICTXT, NPROW, NPCOL,
     $                         MYROW, MYCOL, IAM, NTEST, NFAIL )
      IMPLICIT NONE
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER*8          M8, N8, NB8
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, IAM
      INTEGER            NTEST, NFAIL
*
      INTEGER*8          LRA, LCA, LLDA, I8, J8, GI, GJ
      INTEGER*8          DESCA8( 9 )
      INTEGER            DESCA4( 9 ), M4, N4, NB4, INFO
      COMPLEX, ALLOCATABLE :: A(:)
      REAL, ALLOCATABLE :: WORK(:)
      REAL               VAL8, VAL4
      CHARACTER          NORMS( 4 )
      INTEGER            K
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      REAL               PCLANGE_I8, PCLANGE
      EXTERNAL           PCLANGE_I8, PCLANGE
      EXTERNAL           DESCINIT_I8, DESCINIT
      INTRINSIC          CMPLX, REAL, MAX, INT, ABS
*
      DATA NORMS / 'M', '1', 'I', 'F' /
*
      M4  = INT( M8 )
      N4  = INT( N8 )
      NB4 = INT( NB8 )
      LRA = NUMROC_I8( M8, NB8, MYROW, 0, NPROW )
      LCA = NUMROC_I8( N8, NB8, MYCOL, 0, NPCOL )
      LLDA = MAX( LRA, 1_8 )
*
      INFO = 0
      CALL DESCINIT_I8( DESCA8, M8, N8, NB8, NB8, 0, 0, ICTXT,
     $                  LLDA, INFO )
      CALL DESCINIT( DESCA4, M4, N4, NB4, NB4, 0, 0, ICTXT,
     $               INT( LLDA ), INFO )
*
      ALLOCATE( A( MAX( LLDA*LCA, 1_8 ) ) )
      ALLOCATE( WORK( MAX( LRA, LCA, 1_8 ) ) )
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
      DO K = 1, 4
         VAL8 = PCLANGE_I8( NORMS( K ), M8, N8, A, 1_8, 1_8,
     $                       DESCA8, WORK )
         VAL4 = PCLANGE( NORMS( K ), M4, N4, A, 1, 1,
     $                    DESCA4, WORK )
         NTEST = NTEST + 1
         IF( VAL8 .NE. VAL4 ) THEN
            NFAIL = NFAIL + 1
            IF( IAM .EQ. 0 )
     $         WRITE(*,'(A,A,A)') '  PCLANGE_I8 norm=',
     $                            NORMS( K ), ': FAILED'
         ELSE
            IF( IAM .EQ. 0 )
     $         WRITE(*,'(A,A,A)') '  PCLANGE_I8 norm=',
     $                            NORMS( K ), ': PASSED'
         END IF
      END DO
*
      DEALLOCATE( A, WORK )
      END
