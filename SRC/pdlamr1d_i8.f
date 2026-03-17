      SUBROUTINE PDLAMR1D_I8( N, A, IA, JA, DESCA,
     $                        B, IB, JB, DESCB )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 version of PDLAMR1D.  Supports matrix dimensions
*     larger than 2^31-1.  Descriptor width is determined by the
*     entry-point name (_I8 suffix); DTYPE remains semantic (1 = 2D).
*
*     .. Scalar Arguments ..
      INTEGER*8          N, IA, JA, IB, JB
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAMR1D_I8 redistributes a one-dimensional row vector from one
*  data decomposition to another.
*
*  This is an auxiliary routine called by I8-width reduction routines
*  to redistribute D, E and TAU vectors.
*
*  Notes
*  =====
*
*  Although all processes call PDGEMR2D_I8, only the processes that
*  own the first column of A send data and only processes that own
*  the first column of B receive data.  The calls to DGEBS2D/DGEBR2D
*  spread the data down.
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER*8
*          The size of the matrix to be transposed.
*
*  A       (local output) DOUBLE PRECISION pointer into the local
*          memory to an array of dimension (LOCc(JA+N-1)).
*
*  IA      (global input) INTEGER*8
*          A's global row index.
*
*  JA      (global input) INTEGER*8
*          A's global column index.
*
*  DESCA   (global and local input) INTEGER*8 array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  B       (local input/local output) DOUBLE PRECISION pointer into
*          the local memory to an array of dimension (LOCc(JB+N-1)).
*
*  IB      (global input) INTEGER*8
*          B's global row index, NOT USED.
*
*  JB      (global input) INTEGER*8
*          B's global column index.
*
*  DESCB   (global and local input) INTEGER*8 array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
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
     $                   DGEBR2D, DGEBS2D, PDGEMR2D_I8
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
*     Quick return if possible
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
      CALL PDGEMR2D_I8( 1_8, N, A, IA, JA, DESCAA,
     $                   B, IB, JB, DESCBB, INT( ICTXT, 8 ) )
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NQ8 = NUMROC_I8( N, DESCB( NB_ ), MYCOL,
     $                  INT( DESCB( CSRC_ ) ), NPCOL )
*
*     Narrow NQ for BLACS broadcast (which takes default INTEGER).
*     For a 1-row redistribution NQ is the local column count and
*     can exceed INT_MAX only for extreme N with tiny NB.
*     If it overflows, abort collectively — a partial result after
*     PDGEMR2D_I8 without the column broadcast is unsafe.
*
      IF( NQ8 .GT. 2147483647_8 ) THEN
         CALL PXERBLA( ICTXT, 'PDLAMR1D_I8', -1 )
         CALL BLACS_ABORT( ICTXT, 1 )
      END IF
      NQ4 = INT( NQ8 )
*
      IF( MYROW.EQ.0 ) THEN
         CALL DGEBS2D( ICTXT, 'C', ' ', NQ4, 1, B, NQ4 )
      ELSE
         CALL DGEBR2D( ICTXT, 'C', ' ', NQ4, 1, B, NQ4, 0, MYCOL )
      END IF
*
      RETURN
*
*     End of PDLAMR1D_I8
*
      END
