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

static char * PB_Cmalloc_impl( ScaLAPACK_ByteCount LENGTH )
{
   char           * bufptr = NULL;

   if( LENGTH > 0 )
   {
      if( !( bufptr = (char *) malloc( LENGTH ) ) )
      {
         (void) fprintf( stderr, "Not enough memory on line %d of file %s!!\n",
                         __LINE__, __FILE__ );
         Cblacs_abort( -1, -1 );
      }
   }
   return( bufptr );
}

#ifdef __STDC__
char * PB_Cmalloc64( ScaLAPACK_ByteCount LENGTH )
#else
char * PB_Cmalloc64( LENGTH )
   ScaLAPACK_ByteCount LENGTH;
#endif
{
   return( PB_Cmalloc_impl( LENGTH ) );
}

#ifdef __STDC__
char * PB_Cmalloc( Int LENGTH )
#else
char * PB_Cmalloc( LENGTH )
   Int            LENGTH;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cmalloc allocates a dynamic memory buffer. In case of failure, the
*  program is stopped by calling Cblacs_abort.
*
*  Arguments
*  =========
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
   ScaLAPACK_ByteCount length_bytes;

   if( LENGTH <= 0 ) return( NULL );
   if( !PB_CSizeFromInt( LENGTH, &length_bytes ) ) Cblacs_abort( -1, -1 );
   return( PB_Cmalloc_impl( length_bytes ) );
/*
*  End of PB_Cmalloc
*/
}
