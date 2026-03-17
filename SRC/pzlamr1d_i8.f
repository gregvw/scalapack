      SUBROUTINE PZLAMR1D_I8( N, A, IA, JA, DESCA,
     $                        B, IB, JB, DESCB )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 version of PZLAMR1D.  See PDLAMR1D_I8 for full docs.
*
*     .. Scalar Arguments ..
      INTEGER*8          N, IA, JA, IB, JB
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * )
      DOUBLE COMPLEX     A( * ), B( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICTXT, MYCOL, MYROW, NPCOL, NPROW, NQ4
      INTEGER*8          NQ8
*     ..
*     .. Local Arrays ..
      INTEGER*8          DESCAA( DLEN_ ), DESCBB( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, BLACS_ABORT, PXERBLA,
     $                   ZGEBR2D, ZGEBS2D, PZGEMR2D_I8
*     ..
*     .. External Functions ..
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 )
     $   RETURN
*
      DO 10 I = 1, DLEN_
         DESCAA( I ) = DESCA( I )
         DESCBB( I ) = DESCB( I )
   10 CONTINUE
*
      DESCAA( M_ ) = 1
      DESCBB( M_ ) = 1
      DESCAA( LLD_ ) = 1
      DESCBB( LLD_ ) = 1
*
      ICTXT = INT( DESCB( CTXT_ ) )
      CALL PZGEMR2D_I8( 1_8, N, A, IA, JA, DESCAA,
     $                   B, IB, JB, DESCBB, INT( ICTXT, 8 ) )
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NQ8 = NUMROC_I8( N, DESCB( NB_ ), MYCOL,
     $                  INT( DESCB( CSRC_ ) ), NPCOL )
*
      IF( NQ8 .GT. 2147483647_8 ) THEN
         CALL PXERBLA( ICTXT, 'PZLAMR1D_I8', -1 )
         CALL BLACS_ABORT( ICTXT, 1 )
      END IF
      NQ4 = INT( NQ8 )
*
      IF( MYROW.EQ.0 ) THEN
         CALL ZGEBS2D( ICTXT, 'C', ' ', NQ4, 1, B, NQ4 )
      ELSE
         CALL ZGEBR2D( ICTXT, 'C', ' ', NQ4, 1, B, NQ4, 0, MYCOL )
      END IF
*
      RETURN
*
*     End of PZLAMR1D_I8
*
      END
