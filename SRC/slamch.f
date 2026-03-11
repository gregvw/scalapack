      REAL FUNCTION SLAMCH( CMACH )
      IMPLICIT NONE
*
*  -- Local LAPACK compatibility routine --
*
*  Provide SLAMCH from this tree so single-precision code does not rely on
*  vendor LAPACK implementations with incompatible Fortran ABIs.
*
*  .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      REAL               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      EPS = EPSILON( ONE )
      BASE = REAL( RADIX( ONE ) )
      T = REAL( DIGITS( ONE ) )
      RND = ONE
      EMIN = REAL( MINEXPONENT( ONE ) )
      RMIN = TINY( ONE )
      EMAX = REAL( MAXEXPONENT( ONE ) )
      RMAX = HUGE( ONE )
      PREC = EPS * BASE
      SFMIN = RMIN
      SMALL = ONE / RMAX
      IF( SMALL.GE.SFMIN ) THEN
         SFMIN = SMALL * ( ONE + EPS )
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      ELSE
         RMACH = ZERO
      END IF
*
      SLAMCH = RMACH
      RETURN
*
*     End of SLAMCH
*
      END
