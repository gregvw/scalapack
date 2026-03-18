      SUBROUTINE PCGESV_I8( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB,
     $                      JB, DESCB, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 native version of PCGESV.
*     All public integer arguments are INTEGER*8 except IPIV and INFO.
*     Calls PCGETRF_I8 and PCGETRS_I8 internally.
*
*     .. Scalar Arguments ..
      INTEGER*8          N, NRHS, IA, JA, IB, JB
      INTEGER            INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * )
      INTEGER            IPIV( * )
      COMPLEX            A( * ), B( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, MYCOL, MYROW, NPCOL, NPROW
      INTEGER            IAROW, IBROW
      INTEGER*8          ICOFFA, IROFFA, IROFFB
*     ..
*     .. Local Arrays ..
      INTEGER*8          IDUM1( 1 )
      INTEGER            IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT_I8, PCHK2MAT_I8,
     $                   PCGETRF_I8, PCGETRS_I8, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P_I8
      EXTERNAL           INDXG2P_I8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MOD
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
         INFO = -(600+CTXT_)
      ELSE
         CALL CHK1MAT_I8( N, 1, N, 1, IA, JA, DESCA, 6, INFO )
         CALL CHK1MAT_I8( N, 1, NRHS, 2, IB, JB, DESCB, 11, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P_I8( IA, DESCA( MB_ ), MYROW,
     $                          INT(DESCA( RSRC_ )), NPROW )
            IBROW = INDXG2P_I8( IB, DESCB( MB_ ), MYROW,
     $                          INT(DESCB( RSRC_ )), NPROW )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            IF( IROFFA.NE.0 ) THEN
               INFO = -4
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(600+NB_)
            ELSE IF( IBROW.NE.IAROW .OR. ICOFFA.NE.IROFFB ) THEN
               INFO = -9
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(1100+NB_)
            ELSE IF( ICTXT.NE.INT( DESCB( CTXT_ ) ) ) THEN
               INFO = -(1100+CTXT_)
            END IF
         END IF
         CALL PCHK2MAT_I8( N, 1, N, 1, IA, JA, DESCA, 6, N, 1, NRHS,
     $                     2, IB, JB, DESCB, 11, 0, IDUM1, IDUM2,
     $                     INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCGESV_I8', -INFO )
         RETURN
      END IF
*
*     Compute the LU factorization of sub( A ).
*
      CALL PCGETRF_I8( N, N, A, IA, JA, DESCA, IPIV, INFO )
*
      IF( INFO.EQ.0 ) THEN
*
         CALL PCGETRS_I8( 'No transpose', N, NRHS, A, IA, JA, DESCA,
     $                    IPIV, B, IB, JB, DESCB, INFO )
*
      END IF
*
      RETURN
*
*     End of PCGESV_I8
*
      END
