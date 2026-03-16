      SUBROUTINE PCHK1MAT_I8( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     $                        DESCAPOS0, NEXTRA, EX, EXPOS, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 version of PCHK1MAT.  Checks that the values associated
*     with one distributed matrix are consistent across the entire
*     process grid.  Uses GLOBCHK_I8 (C helper) for 64-bit transport.
*
*     .. Scalar Arguments ..
      INTEGER*8          MA, NA, IA, JA
      INTEGER            MAPOS0, NAPOS0, DESCAPOS0, NEXTRA, INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), EX( NEXTRA )
      INTEGER            EXPOS( NEXTRA )
*     ..
*
*  Purpose
*  =======
*
*  PCHK1MAT_I8 checks that INTEGER*8 values associated with one
*  distributed matrix are consistent across the entire process grid.
*
*  Notes
*  =====
*
*  This routine checks that all values are the same across the grid.
*  It does no local checking; it is therefore legal to abuse the
*  definitions of the non-descriptor arguments.
*
*  Arguments
*  =========
*
*  MA      (global input) INTEGER*8
*          The global number of matrix rows of A being operated on.
*
*  MAPOS0  (global input) INTEGER
*          Where in the calling routine's parameter list MA appears.
*
*  NA      (global input) INTEGER*8
*          The global number of matrix columns of A being operated on.
*
*  NAPOS0  (global input) INTEGER
*          Where in the calling routine's parameter list NA appears.
*
*  IA      (global input) INTEGER*8
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER*8
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER*8 array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  DESCAPOS0 (global input) INTEGER
*          Where in the calling routine's parameter list DESCA
*          appears.  Note that we assume IA and JA are respectively 2
*          and 1 entries behind DESCA.
*
*  NEXTRA  (global input) INTEGER
*          The number of extra INTEGER*8 parameters to check.
*          NEXTRA <= LDW - 11.
*
*  EX      (local input) INTEGER*8 array of dimension (NEXTRA)
*          The values of these extra parameters.
*
*  EXPOS   (local input) INTEGER array of dimension (NEXTRA)
*          The parameter list positions of these extra values.
*
*  INFO    (local input/global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER            BIGNUM, DESCMULT, LDW
      PARAMETER          ( DESCMULT = 100, BIGNUM = DESCMULT * DESCMULT,
     $                     LDW = 25 )
*     ..
*     .. Local Scalars ..
      INTEGER            DESCPOS, K
*     ..
*     .. Local Arrays ..
      INTEGER*8          IVAL( LDW ), IWORK8( LDW )
      INTEGER            IPOS( LDW )
*     ..
*     .. External Subroutines ..
      EXTERNAL           GLOBCHK_I8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*     .. Executable Statements ..
*
*     Want to find errors with MIN( ), so if no error, set it to a big
*     number. If there already is an error, multiply by the the
*     descriptor multiplier.
*
      IF( INFO.GE.0 ) THEN
         INFO = BIGNUM
      ELSE IF( INFO.LT.-DESCMULT ) THEN
         INFO = -INFO
      ELSE
         INFO = -INFO * DESCMULT
      END IF
*
*     Pack values and their positions in the parameter list, factoring
*     in the descriptor multiplier
*
      IVAL( 1 ) = MA
      IPOS( 1 ) = MAPOS0 * DESCMULT
      IVAL( 2 ) = NA
      IPOS( 2 ) = NAPOS0 * DESCMULT
      IVAL( 3 ) = IA
      IPOS( 3 ) = (DESCAPOS0-2) * DESCMULT
      IVAL( 4 ) = JA
      IPOS( 4 ) = (DESCAPOS0-1) * DESCMULT
      DESCPOS = DESCAPOS0 * DESCMULT
*
      IVAL(  5 ) = DESCA( DTYPE_ )
      IPOS(  5 ) = DESCPOS + DTYPE_
      IVAL(  6 ) = DESCA( M_ )
      IPOS(  6 ) = DESCPOS + M_
      IVAL(  7 ) = DESCA( N_ )
      IPOS(  7 ) = DESCPOS + N_
      IVAL(  8 ) = DESCA( MB_ )
      IPOS(  8 ) = DESCPOS + MB_
      IVAL(  9 ) = DESCA( NB_ )
      IPOS(  9 ) = DESCPOS + NB_
      IVAL( 10 ) = DESCA( RSRC_ )
      IPOS( 10 ) = DESCPOS + RSRC_
      IVAL( 11 ) = DESCA( CSRC_ )
      IPOS( 11 ) = DESCPOS + CSRC_
*
      IF( NEXTRA.GT.0 ) THEN
         DO 10 K = 1, NEXTRA
            IVAL( 11+K ) = EX( K )
            IPOS( 11+K ) = EXPOS( K )
   10    CONTINUE
      END IF
      K = 11 + NEXTRA
*
*     Get the smallest error detected anywhere (BIGNUM if no error)
*
      CALL GLOBCHK_I8( INT( DESCA( CTXT_ ) ), K, IVAL, IPOS,
     $                 IWORK8, INFO )
*
*     Prepare output: set info = 0 if no error, and divide by DESCMULT if
*     error is not in a descriptor entry
*
      IF( INFO .EQ. BIGNUM ) THEN
         INFO = 0
      ELSE IF( MOD( INFO, DESCMULT ) .EQ. 0 ) THEN
         INFO = -INFO / DESCMULT
      ELSE
         INFO = -INFO
      END IF
*
      RETURN
*
*     End of PCHK1MAT_I8
*
      END
*
      SUBROUTINE PCHK2MAT_I8( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     $                        DESCAPOS0, MB, MBPOS0, NB, NBPOS0,
     $                        IB, JB, DESCB, DESCBPOS0,
     $                        NEXTRA, EX, EXPOS, INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     INTEGER*8 version of PCHK2MAT.  Checks that the values associated
*     with two distributed matrices are consistent across the entire
*     process grid.
*
*     .. Scalar Arguments ..
      INTEGER*8          MA, NA, IA, JA, MB, NB, IB, JB
      INTEGER            MAPOS0, NAPOS0, DESCAPOS0
      INTEGER            MBPOS0, NBPOS0, DESCBPOS0
      INTEGER            NEXTRA, INFO
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * ), EX( NEXTRA )
      INTEGER            EXPOS( NEXTRA )
*     ..
*
*  Purpose
*  =======
*
*  PCHK2MAT_I8 checks that INTEGER*8 values associated with two
*  distributed matrices are consistent across the entire process grid.
*
*  Arguments — same as PCHK2MAT but with INTEGER*8 matrix dimensions
*  and descriptor arrays.  See PCHK1MAT_I8 for detailed descriptions.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      INTEGER            DESCMULT, BIGNUM, LDW
      PARAMETER          ( DESCMULT = 100, BIGNUM = DESCMULT * DESCMULT,
     $                     LDW = 35 )
*     ..
*     .. Local Scalars ..
      INTEGER            K, DESCPOS
*     ..
*     .. Local Arrays ..
      INTEGER*8          IVAL( LDW ), IWORK8( LDW )
      INTEGER            IPOS( LDW )
*     ..
*     .. External Subroutines ..
      EXTERNAL           GLOBCHK_I8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD, INT
*     ..
*     .. Executable Statements ..
*
*     Want to find errors with MIN( ), so if no error, set it to a big
*     number. If there already is an error, multiply by the the
*     descriptor multiplier.
*
      IF( INFO.GE.0 ) THEN
         INFO = BIGNUM
      ELSE IF( INFO.LT.-DESCMULT ) THEN
         INFO = -INFO
      ELSE
         INFO = -INFO * DESCMULT
      END IF
*
*     Pack values and their positions in the parameter list, factoring
*     in the descriptor multiplier
*
      IVAL( 1 ) = MA
      IPOS( 1 ) = MAPOS0 * DESCMULT
      IVAL( 2 ) = NA
      IPOS( 2 ) = NAPOS0 * DESCMULT
      IVAL( 3 ) = IA
      IPOS( 3 ) = (DESCAPOS0-2) * DESCMULT
      IVAL( 4 ) = JA
      IPOS( 4 ) = (DESCAPOS0-1) * DESCMULT
      DESCPOS = DESCAPOS0 * DESCMULT
*
      IVAL(  5 ) = DESCA( DTYPE_ )
      IPOS(  5 ) = DESCPOS + DTYPE_
      IVAL(  6 ) = DESCA( M_ )
      IPOS(  6 ) = DESCPOS + M_
      IVAL(  7 ) = DESCA( N_ )
      IPOS(  7 ) = DESCPOS + N_
      IVAL(  8 ) = DESCA( MB_ )
      IPOS(  8 ) = DESCPOS + MB_
      IVAL(  9 ) = DESCA( NB_ )
      IPOS(  9 ) = DESCPOS + NB_
      IVAL( 10 ) = DESCA( RSRC_ )
      IPOS( 10 ) = DESCPOS + RSRC_
      IVAL( 11 ) = DESCA( CSRC_ )
      IPOS( 11 ) = DESCPOS + CSRC_
*
      IVAL( 12 ) = MB
      IPOS( 12 ) = MBPOS0 * DESCMULT
      IVAL( 13 ) = NB
      IPOS( 13 ) = NBPOS0 * DESCMULT
      IVAL( 14 ) = IB
      IPOS( 14 ) = (DESCBPOS0-2) * DESCMULT
      IVAL( 15 ) = JB
      IPOS( 15 ) = (DESCBPOS0-1) * DESCMULT
      DESCPOS = DESCBPOS0 * DESCMULT
*
      IVAL( 16 ) = DESCB( DTYPE_ )
      IPOS( 16 ) = DESCPOS + DTYPE_
      IVAL( 17 ) = DESCB( M_ )
      IPOS( 17 ) = DESCPOS + M_
      IVAL( 18 ) = DESCB( N_ )
      IPOS( 18 ) = DESCPOS + N_
      IVAL( 19 ) = DESCB( MB_ )
      IPOS( 19 ) = DESCPOS + MB_
      IVAL( 20 ) = DESCB( NB_ )
      IPOS( 20 ) = DESCPOS + NB_
      IVAL( 21 ) = DESCB( RSRC_ )
      IPOS( 21 ) = DESCPOS + RSRC_
      IVAL( 22 ) = DESCB( CSRC_ )
      IPOS( 22 ) = DESCPOS + CSRC_
*
      IF( NEXTRA.GT.0 ) THEN
         DO 10 K = 1, NEXTRA
            IVAL( 22+K ) = EX( K )
            IPOS( 22+K ) = EXPOS( K )
   10    CONTINUE
      END IF
      K = 22 + NEXTRA
*
*     Get the smallest error detected anywhere (BIGNUM if no error)
*
      CALL GLOBCHK_I8( INT( DESCA( CTXT_ ) ), K, IVAL, IPOS,
     $                 IWORK8, INFO )
*
*     Prepare output: set info = 0 if no error, and divide by DESCMULT
*     if error is not in a descriptor entry.
*
      IF( INFO.EQ.BIGNUM ) THEN
         INFO = 0
      ELSE IF( MOD( INFO, DESCMULT ) .EQ. 0 ) THEN
         INFO = -INFO / DESCMULT
      ELSE
         INFO = -INFO
      END IF
*
      RETURN
*
*     End of PCHK2MAT_I8
*
      END
