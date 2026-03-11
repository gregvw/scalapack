

/* ---------------------------------------------------------------------
*
*  -- ScaLAPACK routine (version 1.5) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*  ---------------------------------------------------------------------
*/

/*
 * Include Files
 */
#include "pxsyevx.h"
#include "pblas.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define  proto(x)	()

static void scalapack_double_words(double x, ScaLAPACK_UWord32 words[2])
{
   memcpy(words, &x, sizeof(ScaLAPACK_UWord32) * 2);
}

static Int scalapack_double_signbit(double x, Int ieflag)
{
   ScaLAPACK_UWord32 words[2];
   scalapack_double_words(x, words);
   if( ieflag == 1 ) return (Int) ((words[0] >> 31) & 1U);
   return (Int) ((words[1] >> 31) & 1U);
}


void pdlasnbt_( Int *ieflag )
{
/* 
*
*  Purpose
*  ======= 
*
*  pdalsnbt finds the position of the signbit of a double
*  double precision floating point number. This routine assumes IEEE
*  arithmetic, and hence, tests only the 32nd and 64th bits as
*  possibilities for the sign bit.
*
*  Note : For this release, we assume that sizeof(int) is 4 bytes.
*
*  Note : If a compile time flag (NO_IEEE) indicates that the
*  machine does not have IEEE arithmetic, IEFLAG = 0 is returned.
*
*  Arguments
*  =========
*
*  IEFLAG   (output) INTEGER
*           This indicates the position of the signbit of any double
*           precision floating point number.
*           IEFLAG = 0 if the compile time flag, NO_IEEE, indicates
*           that the machine does not have IEEE  arithmetic, or if
*           sizeof(int) is different from 4 bytes.
*           IEFLAG = 1 indicates that the sign bit is the 32nd
*           bit ( Big Endian ).
*           IEFLAG = 2 indicates that the sign bit is the 64th
*           bit ( Little Endian ).
*
*  =====================================================================
*
*  .. Local Scalars ..
*/
   double x;
   ScaLAPACK_UWord32 words[2];
/* ..
*  .. Executable Statements ..
*/

#ifdef NO_IEEE
   *ieflag = 0;
#else
   x = (double) -1.0;
   scalapack_double_words(x, words);
   if(( words[0] == 0xbff00000U ) && ( words[1] == 0x0U ) )
   {
      *ieflag = 1;
   } else if(( words[1] == 0xbff00000U ) && ( words[0] == 0x0U ) ) {
      *ieflag = 2;
   } else {
      *ieflag = 0; 
   }
#endif
}

void pdlaiectb_( double *sigma, Int *n, double *d, Int *count )
{
/* 
*
*  Purpose
*  ======= 
*
*  pdlaiectb computes the number of negative eigenvalues of (A- SIGMA I).
*  This implementation of the Sturm Sequence loop exploits IEEE Arithmetic
*  and has no conditionals in the innermost loop. To extract the signbit,
*  this routine assumes that the double precision word is stored in
*  "Big Endian" word order, i.e, the signbit is assumed to be bit 32.
*
*  Note that all arguments are call-by-reference so that this routine
*  can be directly called from Fortran code.
*
*  This is a ScaLAPACK internal subroutine and arguments are not
*  checked for unreasonable values.
*
*  Arguments
*  =========
*
*  SIGMA    (input) DOUBLE PRECISION
*           The shift. pdlaiectb finds the number of eigenvalues
*           less than equal to SIGMA.
*
*  N        (input) INTEGER
*           The order of the tridiagonal matrix T. N >= 1.
*
*  D        (input) DOUBLE PRECISION array, dimension (2*N - 1)
*           Contains the diagonals and the squares of the off-diagonal
*           elements of the tridiagonal matrix T. These elements are
*           assumed to be interleaved in memory for better cache
*           performance. The diagonal entries of T are in the entries
*           D(1),D(3),...,D(2*N-1), while the squares of the off-diagonal
*           entries are D(2),D(4),...,D(2*N-2). To avoid overflow, the
*           matrix must be scaled so that its largest entry is no greater
*           than overflow**(1/2) * underflow**(1/4) in absolute value,
*           and for greatest accuracy, it should not be much smaller
*           than that.
*
*  COUNT    (output) INTEGER
*           The count of the number of eigenvalues of T less than or
*           equal to SIGMA.
*
*  =====================================================================
*
*  .. Local Scalars ..
*/
   double      lsigma, tmp, *pd, *pe2;
   Int         i;
/* ..
*  .. Executable Statements ..
*/

   lsigma = *sigma;
   pd = d; pe2 = d+1;
   tmp = *pd - lsigma; pd += 2;
   *count = scalapack_double_signbit(tmp, 1);
   for(i = 1;i < *n;i++){
      tmp = *pd - *pe2/tmp - lsigma;
      pd += 2; pe2 += 2;
      *count += scalapack_double_signbit(tmp, 1);
   }
}

void pdlaiectl_( double *sigma, Int *n, double *d, Int *count )
{
/* 
*
*  Purpose
*  ======= 
*
*  pdlaiectl computes the number of negative eigenvalues of (A- SIGMA I).
*  This implementation of the Sturm Sequence loop exploits IEEE Arithmetic
*  and has no conditionals in the innermost loop. To extract the signbit,
*  this routine assumes that the double precision word is stored in
*  "Little Endian" word order, i.e, the signbit is assumed to be bit 64.
*
*  Note that all arguments are call-by-reference so that this routine
*  can be directly called from Fortran code.
*
*  This is a ScaLAPACK internal subroutine and arguments are not
*  checked for unreasonable values.
*
*  Arguments
*  =========
*
*  SIGMA    (input) DOUBLE PRECISION
*           The shift. pdlaiectl finds the number of eigenvalues
*           less than equal to SIGMA.
*
*  N        (input) INTEGER
*           The order of the tridiagonal matrix T. N >= 1.
*
*  D        (input) DOUBLE PRECISION array, dimension (2*N - 1)
*           Contains the diagonals and the squares of the off-diagonal
*           elements of the tridiagonal matrix T. These elements are
*           assumed to be interleaved in memory for better cache
*           performance. The diagonal entries of T are in the entries
*           D(1),D(3),...,D(2*N-1), while the squares of the off-diagonal
*           entries are D(2),D(4),...,D(2*N-2). To avoid overflow, the
*           matrix must be scaled so that its largest entry is no greater
*           than overflow**(1/2) * underflow**(1/4) in absolute value,
*           and for greatest accuracy, it should not be much smaller
*           than that.
*
*  COUNT    (output) INTEGER
*           The count of the number of eigenvalues of T less than or
*           equal to SIGMA.
*
*  =====================================================================
*
*  .. Local Scalars ..
*/
   double      lsigma, tmp, *pd, *pe2;
   Int         i;
/* ..
*  .. Executable Statements ..
*/

   lsigma = *sigma;
   pd = d; pe2 = d+1;
   tmp = *pd - lsigma; pd += 2;
   *count = scalapack_double_signbit(tmp, 2);
   for(i = 1;i < *n;i++){
      tmp = *pd - *pe2/tmp - lsigma;
      pd += 2; pe2 += 2;
      *count += scalapack_double_signbit(tmp, 2);
   }
}

void pdlachkieee_( Int *isieee, double *rmax, double *rmin )
{
/* 
*
*  Purpose
*  ======= 
*
*  pdlachkieee performs a simple check to make sure that the features
*  of the IEEE standard that we rely on are implemented.  In some
*  implementations, pdlachkieee may not return.
*
*  Note that all arguments are call-by-reference so that this routine
*  can be directly called from Fortran code.
*
*  This is a ScaLAPACK internal subroutine and arguments are not
*  checked for unreasonable values.
*
*  Arguments
*  =========
*
*  ISIEEE   (local output) INTEGER
*           On exit, ISIEEE = 1 implies that all the features of the
*           IEEE standard that we rely on are implemented.
*           On exit, ISIEEE = 0 implies that some the features of the
*           IEEE standard that we rely on are missing.
*
*  RMAX     (local input) DOUBLE PRECISION
*           The overflow threshold ( = DLAMCH('O') ).
*
*  RMIN     (local input) DOUBLE PRECISION
*           The underflow threshold ( = DLAMCH('U') ).
*
*  =====================================================================
*
*  .. Local Scalars ..
*/
   double pinf, pzero, ninf, nzero;
   Int         ieflag, sbit1, sbit2;
/* ..
*  .. Executable Statements ..
*/

   *isieee = 1;
   pdlasnbt_( &ieflag );
   if( ( ieflag != 1 ) && ( ieflag != 2 ) ){
      *isieee = 0;
      return;
   }

   pinf = *rmax / *rmin;
   pzero = 1.0 / pinf;
   pinf = 1.0 / pzero;

   if( pzero != 0.0 ){
      printf("pzero = %g should be zero\n",pzero);
      *isieee = 0; 
      return ;
   }
   sbit1 = scalapack_double_signbit(pzero, ieflag);
   sbit2 = scalapack_double_signbit(pinf, ieflag);
   if( sbit1 == 1 ){
      printf("Sign of positive infinity is incorrect\n");
      *isieee = 0;
   }
   if( sbit2 == 1 ){
      printf("Sign of positive zero is incorrect\n");
      *isieee = 0;
   }

   ninf = -pinf;
   nzero = 1.0 / ninf;
   ninf = 1.0 / nzero;

   if( nzero != 0.0 ){
      printf("nzero = %g should be zero\n",nzero);
      *isieee = 0;
   }
   sbit1 = scalapack_double_signbit(nzero, ieflag);
   sbit2 = scalapack_double_signbit(ninf, ieflag);
   if( sbit1 == 0 ){
      printf("Sign of negative infinity is incorrect\n");
      *isieee = 0;
   }
   if( sbit2 == 0 ){
      printf("Sign of negative zero is incorrect\n");
      *isieee = 0;
   }
}
