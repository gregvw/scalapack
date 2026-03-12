/* ---------------------------------------------------------------------
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*  ---------------------------------------------------------------------
*/
/*
*  Include files
*/
#include "../pblas.h"
#include "../PBpblas.h"
#include "../PBtools.h"
#include "../PBblacs.h"
#include "../PBblas.h"

static char    * pblasbuf = NULL;
static ScaLAPACK_ByteCount pbbuflen = 0;

static char * PB_Cgetbuf_impl( char * MESS, ScaLAPACK_ByteCount LENGTH )
{
   if( LENGTH > pbbuflen )
   {
      if( pblasbuf ) free( pblasbuf );
      pblasbuf = (char *) malloc( LENGTH );
      if( !pblasbuf )
      {
         (void) fprintf( stderr, "ERROR: Memory allocation failed\n%s\n",
                         MESS );
         Cblacs_abort( -1, -1 );
      }
      pbbuflen = LENGTH;
   }
   return( pblasbuf );
}

#ifdef __STDC__
char * PB_Cgetbuf64( char * MESS, ScaLAPACK_ByteCount LENGTH )
#else
char * PB_Cgetbuf64( MESS, LENGTH )
   ScaLAPACK_ByteCount LENGTH;
   char           * MESS;
#endif
{
   return( PB_Cgetbuf_impl( MESS, LENGTH ) );
}

#ifdef __STDC__
char * PB_Cgetbuf( char * MESS, Int LENGTH )
#else
char * PB_Cgetbuf( MESS, LENGTH )
   Int            LENGTH;
   char           * MESS;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cgetbuf  allocates a dynamic memory buffer. The routine checks the
*  size of the already allocated buffer against the value of the  formal
*  parameter LENGTH. If the current buffer is large enough, this a poin-
*  ter to it is returned. Otherwise, this function tries to allocate it.
*  In case of failure,  the program is stopped by calling  Cblacs_abort.
*  When LENGTH is zero, this function returns a NULL pointer. If the va-
*  lue of LENGTH is strictly less than zero, the buffer is released.
*
*  Arguments
*  =========
*
*  MESS    (local input) pointer to CHAR
*          On entry, MESS is a string containing a message to be printed
*          in case of allocation failure.
*
*  LENGTH  (local input) INTEGER
*          On entry, LENGTH  specifies the length in bytes of the buffer
*          to be allocated.  If LENGTH is less or equal than zero,  this
*          function returns NULL.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   ScaLAPACK_ByteCount requested;
/* ..
*  .. Executable Statements ..
*
*/
   if( LENGTH >= 0 )
   {
      if( LENGTH > pbbuflen )
      {
         if( !PB_CSizeFromInt( LENGTH, &requested ) ) Cblacs_abort( -1, -1 );
         pblasbuf = PB_Cgetbuf_impl( MESS, requested );
      }
   }
   else if( pblasbuf )
   {
      free( pblasbuf );
      pblasbuf = NULL;
      pbbuflen = 0;
   }
   return( pblasbuf );
/*
*  End of PB_Cgetbuf
*/
}
