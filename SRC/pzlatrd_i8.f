      SUBROUTINE PZLATRD_I8( UPLO, N, NB, A, IA, JA, DESCA, D, E,
     $                       TAU, W, IW, JW, DESCW, WORK )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 version of PCLATRD (complex Hermitian).
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          N, IA, JA, IW, JW
      INTEGER            NB
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCW( * )
      DOUBLE PRECISION   D( * ), E( * )
      DOUBLE COMPLEX     A( * ), TAU( * ), W( * ), WORK( * )
*     ..
*
*  =====================================================================
*
      INCLUDE 'SL_i8_params.inc'
      DOUBLE COMPLEX     HALF, ONE, ZERO
      PARAMETER          ( HALF = ( 0.5D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*
      INTEGER            IACOL, IAROW, ICTXT, MYCOL, MYROW, NPCOL,
     $                   NPROW
      INTEGER*8          I8, J8, K8, KW8, II8, JJ8, JP8, JWK8, NQ8
      DOUBLE COMPLEX     AII, ALPHA, BETA
*
      INTEGER*8          DESCD( DLEN_ ), DESCE( DLEN_ ),
     $                   DESCWK( DLEN_ )
*
      EXTERNAL           BLACS_GRIDINFO, DESCSET_I8, DGEBR2D,
     $                   DGEBS2D, INFOG2L_I8,
     $                   PZAXPY_I8, PZDOTC_I8, PZELGET_I8,
     $                   PZELSET_I8, PZGEMV_I8, PZLARFG_I8,
     $                   PZSCAL_I8, PZHEMV_I8, PZLACGV_I8,
     $                   PDELSET_I8
      LOGICAL            LSAME
      INTEGER*8          NUMROC_I8
      EXTERNAL           LSAME, NUMROC_I8
      INTRINSIC          DCMPLX, MIN, DBLE, MOD, INT
*
      IF( N.LE.0 ) RETURN
*
      ICTXT = INT( DESCA( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NQ8 = MAX( 1_8, NUMROC_I8( JA+N-1, DESCA( NB_ ), MYCOL,
     $           INT( DESCA( CSRC_ ) ), NPCOL ) )
      CALL DESCSET_I8( DESCD, 1_8, JA+N-1, 1_8, DESCA( NB_ ),
     $                 MYROW, INT( DESCA( CSRC_ ) ), ICTXT, 1_8 )
      AII = ZERO
      BETA = ZERO
      JP8 = 1
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
         CALL INFOG2L_I8( N+IA-NB, N+JA-NB, DESCA, NPROW, NPCOL,
     $                    MYROW, MYCOL, II8, JJ8, IAROW, IACOL )
         CALL DESCSET_I8( DESCWK, 1_8, INT( DESCW( NB_ ), 8 ),
     $                    1_8, INT( DESCW( NB_ ), 8 ), IAROW, IACOL,
     $                    ICTXT, 1_8 )
         CALL DESCSET_I8( DESCE, 1_8, JA+N-1, 1_8, DESCA( NB_ ),
     $                    MYROW, INT( DESCA( CSRC_ ) ), ICTXT, 1_8 )
*
         DO 10 J8 = JA+N-1, JA+N-NB, -1
            I8 = IA + J8 - JA
            K8 = J8 - JA + 1
            KW8 = MOD( K8-1, DESCA( MB_ ) ) + 1
*
            CALL PZELGET_I8( 'E', ' ', AII, A, I8, J8, DESCA )
            CALL PZELSET_I8( A, I8, J8, DESCA,
     $                       DCMPLX( DBLE( AII ) ) )
            CALL PZLACGV_I8( N-K8, W, IW+K8-1, JW+KW8, DESCW,
     $                       DESCW( M_ ) )
            CALL PZGEMV_I8( 'No transpose', K8, N-K8, -ONE, A, IA,
     $                      J8+1, DESCA, W, IW+K8-1, JW+KW8, DESCW,
     $                      DESCW( M_ ), ONE, A, IA, J8, DESCA,
     $                      1_8 )
            CALL PZLACGV_I8( N-K8, W, IW+K8-1, JW+KW8, DESCW,
     $                       DESCW( M_ ) )
            CALL PZLACGV_I8( N-K8, A, I8, J8+1, DESCA,
     $                       DESCA( M_ ) )
            CALL PZGEMV_I8( 'No transpose', K8, N-K8, -ONE, W, IW,
     $                      JW+KW8, DESCW, A, I8, J8+1, DESCA,
     $                      DESCA( M_ ), ONE, A, IA, J8, DESCA,
     $                      1_8 )
            CALL PZLACGV_I8( N-K8, A, I8, J8+1, DESCA,
     $                       DESCA( M_ ) )
            CALL PZELGET_I8( 'E', ' ', AII, A, I8, J8, DESCA )
            CALL PZELSET_I8( A, I8, J8, DESCA,
     $                       DCMPLX( DBLE( AII ) ) )
            IF( N-K8.GT.0 )
     $         CALL PZELSET_I8( A, I8, J8+1, DESCA,
     $                          DCMPLX( E( JP8 ) ) )
*
            JP8 = MIN( JJ8+KW8-1, NQ8 )
            CALL PZLARFG_I8( K8-1, BETA, I8-1, J8, A, IA, J8,
     $                       DESCA, 1_8, TAU )
            CALL PDELSET_I8( E, 1_8, J8, DESCE, DBLE( BETA ) )
            CALL PZELSET_I8( A, I8-1, J8, DESCA, ONE )
*
            CALL PZHEMV_I8( 'Upper', K8-1, ONE, A, IA, JA, DESCA,
     $                      A, IA, J8, DESCA, 1_8, ZERO, W, IW,
     $                      JW+KW8-1, DESCW, 1_8 )
*
            JWK8 = MOD( K8-1, DESCWK( NB_ ) ) + 2
            CALL PZGEMV_I8( 'Conjugate transpose', K8-1, N-K8, ONE,
     $                      W, IW, JW+KW8, DESCW, A, IA, J8, DESCA,
     $                      1_8, ZERO, WORK, 1_8, JWK8, DESCWK,
     $                      DESCWK( M_ ) )
            CALL PZGEMV_I8( 'No transpose', K8-1, N-K8, -ONE, A,
     $                      IA, J8+1, DESCA, WORK, 1_8, JWK8,
     $                      DESCWK, DESCWK( M_ ), ONE, W, IW,
     $                      JW+KW8-1, DESCW, 1_8 )
            CALL PZGEMV_I8( 'Conjugate transpose', K8-1, N-K8, ONE,
     $                      A, IA, J8+1, DESCA, A, IA, J8, DESCA,
     $                      1_8, ZERO, WORK, 1_8, JWK8, DESCWK,
     $                      DESCWK( M_ ) )
            CALL PZGEMV_I8( 'No transpose', K8-1, N-K8, -ONE, W,
     $                      IW, JW+KW8, DESCW, WORK, 1_8, JWK8,
     $                      DESCWK, DESCWK( M_ ), ONE, W, IW,
     $                      JW+KW8-1, DESCW, 1_8 )
            CALL PZSCAL_I8( K8-1, TAU( JP8 ), W, IW, JW+KW8-1,
     $                      DESCW, 1_8 )
*
            CALL PZDOTC_I8( K8-1, ALPHA, W, IW, JW+KW8-1, DESCW,
     $                      1_8, A, IA, J8, DESCA, 1_8 )
            IF( MYCOL.EQ.IACOL )
     $         ALPHA = -HALF*TAU( JP8 )*ALPHA
            CALL PZAXPY_I8( K8-1, ALPHA, A, IA, J8, DESCA, 1_8,
     $                      W, IW, JW+KW8-1, DESCW, 1_8 )
            CALL PZELGET_I8( 'E', ' ', BETA, A, I8, J8, DESCA )
            CALL PDELSET_I8( D, 1_8, J8, DESCD, DBLE( BETA ) )
*
   10    CONTINUE
*
      ELSE
*
         CALL INFOG2L_I8( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    II8, JJ8, IAROW, IACOL )
         CALL DESCSET_I8( DESCWK, 1_8, INT( DESCW( NB_ ), 8 ),
     $                    1_8, INT( DESCW( NB_ ), 8 ), IAROW, IACOL,
     $                    ICTXT, 1_8 )
         CALL DESCSET_I8( DESCE, 1_8, JA+N-2, 1_8, DESCA( NB_ ),
     $                    MYROW, INT( DESCA( CSRC_ ) ), ICTXT, 1_8 )
*
         DO 20 J8 = JA, JA+NB-1
            I8 = IA + J8 - JA
            K8 = J8 - JA + 1
*
            CALL PZELGET_I8( 'E', ' ', AII, A, I8, J8, DESCA )
            CALL PZELSET_I8( A, I8, J8, DESCA,
     $                       DCMPLX( DBLE( AII ) ) )
            CALL PZLACGV_I8( K8-1, W, IW+K8-1, JW, DESCW,
     $                       DESCW( M_ ) )
            CALL PZGEMV_I8( 'No transpose', N-K8+1, K8-1, -ONE, A,
     $                      I8, JA, DESCA, W, IW+K8-1, JW, DESCW,
     $                      DESCW( M_ ), ONE, A, I8, J8, DESCA,
     $                      1_8 )
            CALL PZLACGV_I8( K8-1, W, IW+K8-1, JW, DESCW,
     $                       DESCW( M_ ) )
            CALL PZLACGV_I8( K8-1, A, I8, JA, DESCA,
     $                       DESCA( M_ ) )
            CALL PZGEMV_I8( 'No transpose', N-K8+1, K8-1, -ONE, W,
     $                      IW+K8-1, JW, DESCW, A, I8, JA, DESCA,
     $                      DESCA( M_ ), ONE, A, I8, J8, DESCA,
     $                      1_8 )
            CALL PZLACGV_I8( K8-1, A, I8, JA, DESCA,
     $                       DESCA( M_ ) )
            CALL PZELGET_I8( 'E', ' ', AII, A, I8, J8, DESCA )
            CALL PZELSET_I8( A, I8, J8, DESCA,
     $                       DCMPLX( DBLE( AII ) ) )
            IF( K8.GT.1 )
     $         CALL PZELSET_I8( A, I8, J8-1, DESCA,
     $                          DCMPLX( E( JP8 ) ) )
*
            JP8 = MIN( JJ8+K8-1, NQ8 )
            CALL PZLARFG_I8( N-K8, BETA, I8+1, J8, A, I8+2, J8,
     $                       DESCA, 1_8, TAU )
            CALL PDELSET_I8( E, 1_8, J8, DESCE, DBLE( BETA ) )
            CALL PZELSET_I8( A, I8+1, J8, DESCA, ONE )
*
            CALL PZHEMV_I8( 'Lower', N-K8, ONE, A, I8+1, J8+1,
     $                      DESCA, A, I8+1, J8, DESCA, 1_8, ZERO,
     $                      W, IW+K8, JW+K8-1, DESCW, 1_8 )
*
            CALL PZGEMV_I8( 'Conjugate Transpose', N-K8, K8-1, ONE,
     $                      W, IW+K8, JW, DESCW, A, I8+1, J8,
     $                      DESCA, 1_8, ZERO, WORK, 1_8, 1_8,
     $                      DESCWK, DESCWK( M_ ) )
            CALL PZGEMV_I8( 'No transpose', N-K8, K8-1, -ONE, A,
     $                      I8+1, JA, DESCA, WORK, 1_8, 1_8,
     $                      DESCWK, DESCWK( M_ ), ONE, W, IW+K8,
     $                      JW+K8-1, DESCW, 1_8 )
            CALL PZGEMV_I8( 'Conjugate transpose', N-K8, K8-1, ONE,
     $                      A, I8+1, JA, DESCA, A, I8+1, J8, DESCA,
     $                      1_8, ZERO, WORK, 1_8, 1_8, DESCWK,
     $                      DESCWK( M_ ) )
            CALL PZGEMV_I8( 'No transpose', N-K8, K8-1, -ONE, W,
     $                      IW+K8, JW, DESCW, WORK, 1_8, 1_8,
     $                      DESCWK, DESCWK( M_ ), ONE, W, IW+K8,
     $                      JW+K8-1, DESCW, 1_8 )
            CALL PZSCAL_I8( N-K8, TAU( JP8 ), W, IW+K8, JW+K8-1,
     $                      DESCW, 1_8 )
            CALL PZDOTC_I8( N-K8, ALPHA, W, IW+K8, JW+K8-1, DESCW,
     $                      1_8, A, I8+1, J8, DESCA, 1_8 )
            IF( MYCOL.EQ.IACOL )
     $         ALPHA = -HALF*TAU( JP8 )*ALPHA
            CALL PZAXPY_I8( N-K8, ALPHA, A, I8+1, J8, DESCA, 1_8,
     $                      W, IW+K8, JW+K8-1, DESCW, 1_8 )
            CALL PZELGET_I8( 'E', ' ', BETA, A, I8, J8, DESCA )
            CALL PDELSET_I8( D, 1_8, J8, DESCD, DBLE( BETA ) )
*
   20    CONTINUE
*
      END IF
*
      IF( MYCOL.EQ.IACOL ) THEN
         IF( MYROW.EQ.IAROW ) THEN
            CALL DGEBS2D( ICTXT, 'Columnwise', ' ', 1, NB,
     $                    D( JJ8 ), 1 )
         ELSE
            CALL DGEBR2D( ICTXT, 'Columnwise', ' ', 1, NB,
     $                    D( JJ8 ), 1, IAROW, MYCOL )
         END IF
      END IF
*
      RETURN
*
*     End of PZLATRD_I8
*
      END
