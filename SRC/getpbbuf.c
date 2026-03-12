#include "tools.h"

void blacs_abort_( Int *ictxt, Int *errornum );

static char * pblasbuf = NULL;
static ScaLAPACK_ByteCount pbbuflen = 0;
static Int  mone = -1;

char * getpbbuf64( char *mess, ScaLAPACK_ByteCount length )
{
/*
*  Purpose
*  =======
*
*  getpbbuf returns a pointer to a working buffer of size length alloca-
*  ted for the PBLAS routines.
*
* ======================================================================
*
*  .. Local Scalars ..
*/
/* ..
*  .. Executable Statements ..
*/
   if( length > pbbuflen )
   {
      if( pblasbuf )
         free( pblasbuf );
      pblasbuf = (char *) malloc( length );
      if( !pblasbuf )
      {
         fprintf( stderr,
                  "PBLAS %s ERROR: Memory allocation failed\n",
                  mess );
         blacs_abort_( &mone, &mone );
      }
      pbbuflen = length;
   }
   return( pblasbuf );
}

char * getpbbuf( char *mess, Int length )
{
   ScaLAPACK_ByteCount length_bytes;

   if( length < 0 )
   {
      if( pblasbuf )
      {
         free( pblasbuf );
         pblasbuf = NULL;
         pbbuflen = 0;
      }
      return( pblasbuf );
   }
   if( !ScaLAPACK_Index64ToSizeT( (ScaLAPACK_Index64) length, &length_bytes ) )
   {
      fprintf( stderr, "PBLAS %s ERROR: Invalid buffer length\n", mess );
      blacs_abort_( &mone, &mone );
   }
   return( getpbbuf64( mess, length_bytes ) );
}
