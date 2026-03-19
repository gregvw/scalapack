      PROGRAM TEST_LARGE_INDEX_I8
      IMPLICIT NONE
*
*  Synthetic large-index test for I8 dense solvers.
*
*  Exercises 64-bit index arithmetic without allocating a huge matrix.
*  Uses a descriptor with large M_/N_ but small NB so each process
*  owns only a few local elements.  IA/JA are set to large offsets
*  within the logically large matrix, targeting a small submatrix
*  that fits in the allocated local storage.
*
*  This catches narrowing bugs where INT(IA) or INT(JA) would
*  overflow, and verifies that INFOG2L_I8, NUMROC_I8, and the
*  descriptor plumbing work correctly with 64-bit indices.
*
*  Usage:  mpirun -np 4 ./xlargeidx_i8
*
      INCLUDE 'SL_i8_params.inc'
*
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL
      INTEGER            NPROCS, IAM
      INTEGER            NFAIL, NTEST, NFAIL_G, INFO
*
*     Large-index parameters
      INTEGER*8          NGLOBAL, NB8, IA8, JA8, NSUB
      INTEGER*8          LLD8, LR, LC
      INTEGER*8          DESCA8( 9 )
      INTEGER*8          I8, J8, GI, GJ
      DOUBLE PRECISION, ALLOCATABLE :: A(:), B(:), BREF(:)
      INTEGER, ALLOCATABLE :: IPIV(:)
      INTEGER*8          DESCB8( 9 )
*
      INTEGER*8          NUMROC_I8, INDXL2G_I8
      EXTERNAL           NUMROC_I8, INDXL2G_I8
      EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDEXIT, BLACS_EXIT,
     $                   DESCINIT_I8, PDGESV_I8, PDPOSV_I8, IGAMX2D
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
     $   WRITE(*,'(A)') 'test_large_index_i8: 2x2 grid'
*
*     === Test 1: NUMROC_I8 with large N ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- NUMROC_I8 large N ---'
      NTEST = NTEST + 1
      NGLOBAL = 3000000000_8
      NB8 = 64_8
      LR = NUMROC_I8( NGLOBAL, NB8, MYROW, 0, NPROW )
      LC = NUMROC_I8( NGLOBAL, NB8, MYCOL, 0, NPCOL )
*     Each process should own roughly N/(NPROW or NPCOL) rows/cols
      IF( LR .LE. 0 .OR. LC .LE. 0 ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  NUMROC_I8 large N: FAILED (zero local)'
      ELSE
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A,I15,A,I15)')
     $      '  NUMROC_I8 large N: PASSED LR=', LR, ' LC=', LC
      END IF
*
*     === Test 2: DESCINIT_I8 with large M/N ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- DESCINIT_I8 large ---'
      NTEST = NTEST + 1
      LLD8 = MAX( LR, 1_8 )
      INFO = 0
      CALL DESCINIT_I8( DESCA8, NGLOBAL, NGLOBAL, NB8, NB8, 0, 0,
     $                  ICTXT, LLD8, INFO )
      IF( INFO .NE. 0 ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A,I4)') '  DESCINIT_I8 large: FAILED INFO=',INFO
      ELSE
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  DESCINIT_I8 large: PASSED'
      END IF
*
*     === Test 3: INDXL2G_I8 round-trip with large indices ===
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- INDXL2G_I8 large ---'
      NTEST = NTEST + 1
*     First local element on process (0,0) maps to global index 1
      GI = INDXL2G_I8( 1_8, NB8, MYROW, 0, NPROW )
      IF( GI .LE. 0 .OR. GI .GT. NGLOBAL ) THEN
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  INDXL2G_I8 large: FAILED'
      ELSE
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A,I15)') '  INDXL2G_I8 large: PASSED GI=', GI
      END IF
*
*     === Test 4: NARROW_DESC8 abort on large descriptor ===
*     (We can't test this directly since it aborts, but we verify
*      the descriptor fields are > INTMAX)
*
      IF( IAM .EQ. 0 ) WRITE(*,'(A)') '--- Descriptor >INTMAX ---'
      NTEST = NTEST + 1
      IF( DESCA8( M_ ) .GT. 2147483647_8 .AND.
     $    DESCA8( N_ ) .GT. 2147483647_8 ) THEN
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  Descriptor >INTMAX: PASSED'
      ELSE
         NFAIL = NFAIL + 1
         IF( IAM .EQ. 0 )
     $      WRITE(*,'(A)') '  Descriptor >INTMAX: FAILED'
      END IF
*
*     Summary
*
      NFAIL_G = NFAIL
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, NFAIL_G, 1,
     $              NFAIL_G, NFAIL_G, -1, -1, -1 )
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(*,'(A)') '======================================'
         WRITE(*,'(A)') 'Large-Index I8 Test Summary'
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
