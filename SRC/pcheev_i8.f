      SUBROUTINE PCHEEV_I8( JOBZ, UPLO, N, A, IA, JA, DESCA, W,
     $                      Z, IZ, JZ, DESCZ, WORK, LWORK, RWORK,
     $                      LRWORK, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PCHEEV.
*     All public integer arguments are INTEGER*8; INFO stays INTEGER.
*     Uses I8 calls internally; PCLANHE is the one remaining narrowing
*     boundary (requires legacy INTEGER descriptor via NARROW_DESC8).
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER*8          N, IA, JA, IZ, JZ, LWORK, LRWORK
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCZ( * )
      REAL               RWORK( * ), W( * )
      COMPLEX            A( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PCHEEV_I8 computes all eigenvalues and, optionally, eigenvectors
*  of a complex Hermitian matrix A by calling the recommended sequence
*  of ScaLAPACK routines.  This is the INTEGER*8 version of PCHEEV.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
      INTEGER            ITHVAL
      PARAMETER          ( ITHVAL = 10 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, WANTZ
      INTEGER            CONTEXTC, ICTXT, IACOL, IAROW,
     $                   IINFO, ISCALE,
     $                   MYCOL, MYPCOLC, MYPROWC, MYROW, NB,
     $                   NPCOL, NPCOLC, NPROCS, NPROW, NPROWC
      INTEGER*8          CSRC_A8, I8, ICOFFA, INDD, INDE,
     $                   INDRD, INDRE, INDRWORK,
     $                   INDTAU, INDWORK, INDWORK2,
     $                   IROFFA, IROFFZ, IZROW8,
     $                   J8, K8, LDC, LLWORK, LLRWORK,
     $                   LRWMIN, LWMIN,
     $                   MB_A8, MB_Z8, NB_A8, NB_Z8,
     $                   NP, NQ, NRC,
     $                   RSIZECSTEQR2, RSRC_A8, RSRC_Z8,
     $                   SIZECSTEQR2, SIZEPCHETRD, SIZEPCUNMTR
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER*8          DESCQR( 9 ), IDUM1( 3 )
      INTEGER            IDUM2( 3 ), DESCA4( 9 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P_I8, SL_GRIDRESHAPE
      INTEGER*8          NUMROC_I8
      REAL               PCLANHE, PSLAMCH
      EXTERNAL           LSAME, NUMROC_I8, PSLAMCH, PCLANHE,
     $                   SL_GRIDRESHAPE, INDXG2P_I8
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $                   CHK1MAT_I8, CSTEQR2,
     $                   DESCINIT_I8,
     $                   NARROW_DESC8,
     $                   PCHK1MAT_I8, PCHK2MAT_I8,
     $                   PCELGET_I8, PCGEMR2D_I8,
     $                   PCLASCL_I8, PCLASET_I8, PCUNMTR_I8,
     $                   PCHENTRD_I8, PXERBLA,
     $                   SCOPY, SGAMN2D, SGAMX2D, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, ICHAR, INT, MAX, MIN, MOD, REAL,
     $                   SQRT
*     ..
*     .. Executable Statements ..
*
*     Quick return
*
      IF( N.EQ.0 ) RETURN
*
*     Get grid parameters.
*
      ICTXT = INT( DESCA( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
*
*     Initialize pointers to some safe value
*
      INDTAU = 1
      INDD = 1
      INDE = 1
      INDWORK = 1
      INDWORK2 = 1
*
      INDRE = 1
      INDRD = 1
      INDRWORK = 1
*
      WANTZ = LSAME( JOBZ, 'V' )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 700+CTXT_ )
      ELSE IF( WANTZ ) THEN
         IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
            INFO = -( 1200+CTXT_ )
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT_I8( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         IF( WANTZ )
     $      CALL CHK1MAT_I8( N, 3, N, 3, IZ, JZ, DESCZ, 12, INFO )
*
         IF( INFO.EQ.0 ) THEN
*
*           Get machine constants.
*
            SAFMIN = PSLAMCH( ICTXT, 'Safe minimum' )
            EPS = PSLAMCH( ICTXT, 'Precision' )
            SMLNUM = SAFMIN / EPS
            BIGNUM = ONE / SMLNUM
            RMIN = SQRT( SMLNUM )
            RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
            NPROCS = NPROW*NPCOL
            NB_A8 = DESCA( NB_ )
            MB_A8 = DESCA( MB_ )
            NB = INT( NB_A8 )
            LOWER = LSAME( UPLO, 'L' )
*
            RSRC_A8 = DESCA( RSRC_ )
            CSRC_A8 = DESCA( CSRC_ )
            IROFFA = MOD( IA-1, MB_A8 )
            ICOFFA = MOD( JA-1, NB_A8 )
            IAROW = INDXG2P_I8( 1_8, NB_A8, MYROW,
     $                          INT( RSRC_A8 ), NPROW )
            IACOL = INDXG2P_I8( 1_8, MB_A8, MYCOL,
     $                          INT( CSRC_A8 ), NPCOL )
            NP = NUMROC_I8( N+IROFFA, INT( NB, 8 ), MYROW,
     $                      IAROW, NPROW )
            NQ = NUMROC_I8( N+ICOFFA, INT( NB, 8 ), MYCOL,
     $                      IACOL, NPCOL )
*
            IF( WANTZ ) THEN
               NB_Z8 = DESCZ( NB_ )
               MB_Z8 = DESCZ( MB_ )
               RSRC_Z8 = DESCZ( RSRC_ )
               IROFFZ = MOD( IZ-1, MB_A8 )
               IZROW8 = INDXG2P_I8( 1_8, NB_A8, MYROW,
     $                               INT( RSRC_Z8 ), NPROW )
            ELSE
               IROFFZ = 0
               IZROW8 = 0
            END IF
*
*           COMPLEX work space for PCHENTRD_I8
*
            CALL PCHENTRD_I8( UPLO, N, A, IA, JA, DESCA,
     $                        RWORK( INDD ),
     $                        RWORK( INDE ), WORK( INDTAU ),
     $                        WORK( INDWORK ), -1_8, IINFO )
            SIZEPCHETRD = INT( ABS( WORK( 1 ) ), 8 )
*
*           COMPLEX work space for PCUNMTR_I8
*
            IF( WANTZ ) THEN
               CALL PCUNMTR_I8( 'L', UPLO, 'N', N, N, A, IA, JA,
     $                          DESCA,
     $                          WORK( INDTAU ), Z, IZ, JZ, DESCZ,
     $                          WORK( INDWORK ), -1_8, IINFO )
               SIZEPCUNMTR = INT( ABS( WORK( 1 ) ), 8 )
            ELSE
               SIZEPCUNMTR = 0
            END IF
*
*           REAL work space for CSTEQR2
*
            IF( WANTZ ) THEN
               RSIZECSTEQR2 = MAX( 1_8, 2*N-2 )
            ELSE
               RSIZECSTEQR2 = 0
            END IF
*
*           Initialize the context of the single column distributed
*           matrix required by CSTEQR2. This specific distribution
*           allows each process to do 1/pth of the work updating matrix
*           Q during CSTEQR2 and achieve some parallelization to an
*           otherwise serial subroutine.
*
            LDC = 0
            IF( WANTZ ) THEN
               CONTEXTC = SL_GRIDRESHAPE( ICTXT, 0, 1, 1,
     $                    NPROCS, 1 )
               CALL BLACS_GRIDINFO( CONTEXTC, NPROWC, NPCOLC,
     $                              MYPROWC, MYPCOLC )
               NRC = NUMROC_I8( N, NB_A8, MYPROWC, 0, NPROCS )
               LDC = MAX( 1_8, NRC )
               CALL DESCINIT_I8( DESCQR, N, N, INT( NB, 8 ),
     $                           INT( NB, 8 ), 0, 0, CONTEXTC,
     $                           LDC, INFO )
            END IF
*
*           COMPLEX work space for CSTEQR2
*
            IF( WANTZ ) THEN
               SIZECSTEQR2 = N*LDC
            ELSE
               SIZECSTEQR2 = 0
            END IF
*
*           Set up pointers into the WORK array
*
            INDTAU = 1
            INDD = INDTAU + N
            INDE = INDD + N
            INDWORK = INDE + N
            INDWORK2 = INDWORK + N*LDC
            LLWORK = LWORK - INDWORK + 1
*
*           Set up pointers into the RWORK array
*
            INDRE = 1
            INDRD = INDRE + N
            INDRWORK = INDRD + N
            LLRWORK = LRWORK - INDRWORK + 1
*
*           Compute the total amount of space needed
*
            LRWMIN = 2*N + RSIZECSTEQR2
            LWMIN = 3*N + MAX( SIZEPCHETRD, SIZEPCUNMTR,
     $                         SIZECSTEQR2 )
*
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
               INFO = -1
            ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
               INFO = -2
            ELSE IF( LWORK.LT.LWMIN .AND. LWORK.NE.-1 ) THEN
               INFO = -14
            ELSE IF( LRWORK.LT.LRWMIN .AND. LRWORK.NE.-1 ) THEN
               INFO = -16
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            END IF
            IF( WANTZ ) THEN
               IF( IROFFA.NE.IROFFZ ) THEN
                  INFO = -10
               ELSE IF( IAROW.NE.INT( IZROW8 ) ) THEN
                  INFO = -10
               ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
                  INFO = -( 1200+M_ )
               ELSE IF( DESCA( N_ ).NE.DESCZ( N_ ) ) THEN
                  INFO = -( 1200+N_ )
               ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
                  INFO = -( 1200+MB_ )
               ELSE IF( DESCA( NB_ ).NE.DESCZ( NB_ ) ) THEN
                  INFO = -( 1200+NB_ )
               ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
                  INFO = -( 1200+RSRC_ )
               ELSE IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
                  INFO = -( 1200+CTXT_ )
               END IF
            END IF
         END IF
         IF( WANTZ ) THEN
            IDUM1( 1 ) = ICHAR( 'V' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'N' )
         END IF
         IDUM2( 1 ) = 1
         IF( LOWER ) THEN
            IDUM1( 2 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'U' )
         END IF
         IDUM2( 2 ) = 2
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 3 ) = -1
         ELSE
            IDUM1( 3 ) = 1
         END IF
         IDUM2( 3 ) = 3
         IF( WANTZ ) THEN
            CALL PCHK2MAT_I8( N, 3, N, 3, IA, JA, DESCA, 7,
     $                        N, 3, N, 3,
     $                        IZ, JZ, DESCZ, 12, 3, IDUM1, IDUM2,
     $                        INFO )
         ELSE
            CALL PCHK1MAT_I8( N, 3, N, 3, IA, JA, DESCA, 7, 3,
     $                        IDUM1, IDUM2, INFO )
         END IF
*
*     Write the required workspace for lwork/lrwork queries.
*
         WORK( 1 ) = CMPLX( REAL( LWMIN ) )
         RWORK( 1 ) = REAL( LRWMIN )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCHEEV_I8', -INFO )
         IF( WANTZ ) CALL BLACS_GRIDEXIT( CONTEXTC )
         RETURN
      ELSE IF( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 ) THEN
         IF( WANTZ ) CALL BLACS_GRIDEXIT( CONTEXTC )
         RETURN
      END IF
*
*     Scale matrix to allowable range, if necessary.
*     PCLANHE has no _I8 version -- narrow DESCA for this call only.
*
      ISCALE = 0
*
      CALL NARROW_DESC8( DESCA, DESCA4 )
      ANRM = PCLANHE( 'M', UPLO, INT( N ), A, INT( IA ), INT( JA ),
     $                 DESCA4, RWORK( INDRWORK ) )
*
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
*
      IF( ISCALE.EQ.1 ) THEN
         CALL PCLASCL_I8( UPLO, ONE, SIGMA, N, N, A, IA, JA,
     $                     DESCA, IINFO )
      END IF
*
*     Reduce Hermitian matrix to tridiagonal form.
*
      CALL PCHENTRD_I8( UPLO, N, A, IA, JA, DESCA, RWORK( INDRD ),
     $                   RWORK( INDRE ), WORK( INDTAU ),
     $                   WORK( INDWORK ), LLWORK, IINFO )
*
*     Copy the values of D, E to all processes.
*
      DO 10 I8 = 1, N
         CALL PCELGET_I8( 'A', ' ', WORK( INDD+I8-1 ), A,
     $                    I8+IA-1, I8+JA-1, DESCA )
         RWORK( INDRD+I8-1 ) = REAL( WORK( INDD+I8-1 ) )
 10   CONTINUE
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 I8 = 1, N-1
            CALL PCELGET_I8( 'A', ' ', WORK( INDE+I8-1 ), A,
     $                       I8+IA-1, I8+JA, DESCA )
            RWORK( INDRE+I8-1 ) = REAL( WORK( INDE+I8-1 ) )
 20      CONTINUE
      ELSE
         DO 30 I8 = 1, N-1
            CALL PCELGET_I8( 'A', ' ', WORK( INDE+I8-1 ), A,
     $                       I8+IA, I8+JA-1, DESCA )
            RWORK( INDRE+I8-1 ) = REAL( WORK( INDE+I8-1 ) )
 30      CONTINUE
      END IF
*
      IF( WANTZ ) THEN
*
         CALL PCLASET_I8( 'Full', N, N, CZERO, CONE,
     $                    WORK( INDWORK ), 1_8, 1_8, DESCQR )
*
*        CSTEQR2 is a modified version of LAPACK's CSTEQR.  The
*        modifications allow each process to perform partial updates
*        to matrix Q.
*
         CALL CSTEQR2( 'I', INT( N ), RWORK( INDRD ), RWORK( INDRE ),
     $                 WORK( INDWORK ), INT( LDC ), INT( NRC ),
     $                 RWORK( INDRWORK ), INFO )
*
         CALL PCGEMR2D_I8( N, N, WORK( INDWORK ), 1_8, 1_8, DESCQR,
     $                     Z, IA, JA, DESCZ, INT( CONTEXTC, 8 ) )
*
         CALL PCUNMTR_I8( 'L', UPLO, 'N', N, N, A, IA, JA, DESCA,
     $                    WORK( INDTAU ), Z, IZ, JZ, DESCZ,
     $                    WORK( INDWORK ), LLWORK, IINFO )
*
      ELSE
*
         CALL CSTEQR2( 'N', INT( N ), RWORK( INDRD ), RWORK( INDRE ),
     $                 WORK( INDWORK ), 1, 1, RWORK( INDRWORK ),
     $                 INFO )
      END IF
*
*     Copy eigenvalues from workspace to output array
*
      CALL SCOPY( INT( N ), RWORK( INDRD ), 1, W, 1 )
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         CALL SSCAL( INT( N ), ONE / SIGMA, W, 1 )
      END IF
*
*     Free up resources
*
      IF( WANTZ ) THEN
         CALL BLACS_GRIDEXIT( CONTEXTC )
      END IF
*
*     Compare every ith eigenvalue, or all if there are only a few,
*     across the process grid to check for heterogeneity.
*
      IF( N.LE.ITHVAL ) THEN
         J8 = N
         K8 = 1
      ELSE
         J8 = N / ITHVAL
         K8 = ITHVAL
      END IF
*
      DO 40 I8 = 1, J8
         RWORK( I8 ) = W( (I8-1)*K8+1 )
         RWORK( I8+J8 ) = W( (I8-1)*K8+1 )
 40   CONTINUE
*
      CALL SGAMN2D( ICTXT, 'All', ' ', INT( J8 ), 1,
     $              RWORK( 1 ), INT( J8 ), 1, 1, -1, -1, 0 )
      CALL SGAMX2D( ICTXT, 'All', ' ', INT( J8 ), 1,
     $              RWORK( 1+J8 ), INT( J8 ), 1, 1, -1, -1, 0 )
*
      DO 50 I8 = 1, J8
         IF( INFO.EQ.0 .AND. ( RWORK( I8 )-RWORK( I8+J8 )
     $        .NE. ZERO ) ) THEN
            INFO = INT( N )+1
         END IF
 50   CONTINUE
*
      RETURN
*
*     End of PCHEEV_I8
*
      END
