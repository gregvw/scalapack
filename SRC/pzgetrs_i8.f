      SUBROUTINE PZGETRS_I8( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV,
     $                       B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PZGETRS.
*     All public integer arguments are INTEGER*8 except IPIV and INFO.
*     Calls PZLAPIV_I8 and PZTRSM_I8 internally.
*     For complex types, TRANS='C' means conjugate transpose.
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER*8          N, NRHS, IA, JA, IB, JB
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * )
      INTEGER            IPIV( * )
      COMPLEX*16         A( * ), B( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
      COMPLEX*16         ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            ICTXT, MYCOL, MYROW, NPCOL, NPROW
      INTEGER            IAROW, IBROW
      INTEGER*8          ICOFFA, IROFFA, IROFFB
*     ..
*     .. Local Arrays ..
      INTEGER*8          DESCIP( DLEN_ ), IDUM1( 1 )
      INTEGER            IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT_I8, DESCSET_I8,
     $                   PCHK2MAT_I8, PXERBLA, PZLAPIV_I8, PZTRSM_I8
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P_I8
      INTEGER*8          NUMROC_I8
      EXTERNAL           INDXG2P_I8, LSAME, NUMROC_I8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, INT, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = INT( DESCA( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(700+CTXT_)
      ELSE
         NOTRAN = LSAME( TRANS, 'N' )
         CALL CHK1MAT_I8( N, 2, N, 2, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT_I8( N, 2, NRHS, 3, IB, JB, DESCB, 12, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P_I8( IA, DESCA( MB_ ), MYROW,
     $                          INT(DESCA( RSRC_ )), NPROW )
            IBROW = INDXG2P_I8( IB, DESCB( MB_ ), MYROW,
     $                          INT(DESCB( RSRC_ )), NPROW )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND.
     $          .NOT.LSAME( TRANS, 'C' ) ) THEN
               INFO = -1
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(700+NB_)
            ELSE IF( IROFFB.NE.0 .OR. IBROW.NE.IAROW ) THEN
               INFO = -10
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(1200+NB_)
            ELSE IF( ICTXT.NE.INT( DESCB( CTXT_ ) ) ) THEN
               INFO = -(1200+CTXT_)
            END IF
         END IF
         IF( NOTRAN ) THEN
            IDUM1( 1 ) = ICHAR( 'N' )
         ELSE IF( LSAME( TRANS, 'T' ) ) THEN
            IDUM1( 1 ) = ICHAR( 'T' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'C' )
         END IF
         IDUM2( 1 ) = 1
         CALL PCHK2MAT_I8( N, 2, N, 2, IA, JA, DESCA, 7, N, 2, NRHS,
     $                     3, IB, JB, DESCB, 12, 1, IDUM1, IDUM2,
     $                     INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZGETRS_I8', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*     Build I8 descriptor for IPIV
*
      CALL DESCSET_I8( DESCIP,
     $                 DESCA( M_ ) + DESCA( MB_ )*INT( NPROW, 8 ),
     $                 INT( 1, 8 ),
     $                 DESCA( MB_ ), INT( 1, 8 ),
     $                 INT( DESCA( RSRC_ ) ), MYCOL, ICTXT,
     $                 DESCA( MB_ ) + NUMROC_I8( DESCA( M_ ),
     $                 DESCA( MB_ ), MYROW,
     $                 INT( DESCA( RSRC_ ) ), NPROW ) )
*
      IF( NOTRAN ) THEN
*
         CALL PZLAPIV_I8( 'Forward', 'Row', 'Col', N, NRHS, B, IB, JB,
     $                    DESCB, IPIV, IA, INT( 1, 8 ), DESCIP, IDUM2 )
*
         CALL PZTRSM_I8( 'Left', 'Lower', 'No transpose', 'Unit', N,
     $                   NRHS, ONE, A, IA, JA, DESCA, B, IB, JB,
     $                   DESCB )
*
         CALL PZTRSM_I8( 'Left', 'Upper', 'No transpose', 'Non-unit',
     $                   N, NRHS, ONE, A, IA, JA, DESCA, B, IB, JB,
     $                   DESCB )
      ELSE
*
*        Solve sub( A )**T or sub( A )**H * X = sub( B ).
*        Pass TRANS through to handle 'T' vs 'C' correctly.
*
         CALL PZTRSM_I8( 'Left', 'Upper', TRANS, 'Non-unit', N,
     $                   NRHS, ONE, A, IA, JA, DESCA, B, IB, JB,
     $                   DESCB )
*
         CALL PZTRSM_I8( 'Left', 'Lower', TRANS, 'Unit', N,
     $                   NRHS, ONE, A, IA, JA, DESCA, B, IB, JB,
     $                   DESCB )
*
         CALL PZLAPIV_I8( 'Backward', 'Row', 'Col', N, NRHS, B, IB,
     $                    JB, DESCB, IPIV, IA, INT( 1, 8 ), DESCIP,
     $                    IDUM2 )
*
      END IF
*
      RETURN
*
*     End of PZGETRS_I8
*
      END
