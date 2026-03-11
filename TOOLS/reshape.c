#include <stdlib.h>
#include "scalapack-types.h"

#ifndef Int
#define Int int
#endif

_Static_assert(sizeof(Int) == sizeof(ScaLAPACK_ApiInt),
               "Int must match the configured ScaLAPACK API integer width.");

void Creshape( Int context_in, Int major_in, Int* context_out, Int major_out,
                    Int first_proc, Int nprow_new, Int npcol_new )
/* major in, major out represent whether processors go row major (1) or
column major (2) in the input and output grids */
{

   /** called subprograms **/
   void proc_inc( Int* myrow, Int* mycol, Int nprow, Int npcol, Int major );
   void Cblacs_gridinfo( Int context, Int* nprow, Int* npcol, Int* myrow, Int* mycol );
   Int Cblacs_pnum( Int context, Int prow, Int pcol );
   void Cblacs_abort( Int context, Int errornum );
   void Cblacs_get( Int context, Int what, Int* val );
   void Cblacs_gridmap( Int* context, Int* usermap, Int ldumap, Int nprow, Int npcol );

   /** variables **/
   ScaLAPACK_Index64 grid_offset, nprocs_new64, proc_index;
   size_t grid_elems;
   Int j;
   Int nprow_in, npcol_in, myrow_in, mycol_in;
   Int myrow_old, mycol_old, myrow_new, mycol_new;
   Int pnum;
   Int *grid_new;

/********** executable statements ************/

   if( !ScaLAPACK_Index64Mul( (ScaLAPACK_Index64) nprow_new,
                              (ScaLAPACK_Index64) npcol_new,
                              &nprocs_new64 ) ||
       !ScaLAPACK_Index64ToSizeT( nprocs_new64, &grid_elems ) ||
       grid_elems > SIZE_MAX / sizeof( Int ) )
   {
      Cblacs_abort( context_in, -24 );
      return;
   }

   Cblacs_gridinfo( context_in, &nprow_in, &npcol_in, &myrow_in, &mycol_in );

   /* Quick return if possible */
   if( ( nprow_in == nprow_new ) && ( npcol_in == npcol_new ) &&
       ( first_proc == 0 ) && ( major_in == major_out ) )
   {
      *context_out = context_in;
      return;
   }

   /* allocate space for new process mapping */
   grid_new = (Int *) malloc( grid_elems * sizeof( Int ) );
   if( grid_new == NULL )
   {
      Cblacs_abort( context_in, -25 );
      return;
   }

   /* set place in old grid to start grabbing processors for new grid */
   myrow_old = 0; mycol_old = 0;
   if ( major_in == 1 ) /* row major */
   {
      myrow_old = first_proc / nprow_in;
      mycol_old = first_proc % nprow_in;
   }
   else                  /* col major */
   {
      myrow_old = first_proc % nprow_in;
      mycol_old = first_proc / nprow_in;
   }

   myrow_new = 0; mycol_new = 0;

   /* Set up array of process numbers for new grid */
   for( proc_index = 0; proc_index < nprocs_new64; ++proc_index )
   {
      if( !ScaLAPACK_Index64Mul( (ScaLAPACK_Index64) mycol_new,
                                 (ScaLAPACK_Index64) nprow_new,
                                 &grid_offset ) ||
          !ScaLAPACK_Index64Add( grid_offset,
                                 (ScaLAPACK_Index64) myrow_new,
                                 &grid_offset ) ||
          !ScaLAPACK_Index64ToSizeT( grid_offset, &grid_elems ) )
      {
         free( grid_new );
         Cblacs_abort( context_in, -26 );
         return;
      }
      pnum = Cblacs_pnum( context_in, myrow_old, mycol_old );
      grid_new[grid_elems] = pnum;
      proc_inc( &myrow_old, &mycol_old, nprow_in, npcol_in, major_in );
      proc_inc( &myrow_new, &mycol_new, nprow_new, npcol_new, major_out );
   }

   /* get context */
   Cblacs_get( context_in, 10, context_out );

   /* allocate grid */
   Cblacs_gridmap( context_out, grid_new, nprow_new, nprow_new, npcol_new );

   /* free malloced space */
   free( grid_new );
}

/*************************************************************************/
void reshape( Int* context_in, Int* major_in, Int* context_out, Int* major_out,
                    Int* first_proc, Int* nprow_new, Int* npcol_new )
{
   Creshape( *context_in, *major_in, context_out, *major_out,
                    *first_proc, *nprow_new, *npcol_new );
}
/*************************************************************************/
void RESHAPE( Int* context_in, Int* major_in, Int* context_out, Int* major_out,
                    Int* first_proc, Int* nprow_new, Int* npcol_new )
{
   Creshape( *context_in, *major_in, context_out, *major_out,
                    *first_proc, *nprow_new, *npcol_new );
}
/*************************************************************************/
void reshape_( Int* context_in, Int* major_in, Int* context_out, Int* major_out,
                    Int* first_proc, Int* nprow_new, Int* npcol_new )
{
   Creshape( *context_in, *major_in, context_out, *major_out,
                    *first_proc, *nprow_new, *npcol_new );
}
/*************************************************************************/
void proc_inc( Int* myrow, Int* mycol, Int nprow, Int npcol, Int major )
{
   if( major == 1) /* row major */
   {
      if( *mycol == npcol-1 )
      {
         *mycol = 0;
         if( *myrow == nprow-1 )
         {
            *myrow = 0;
         }
         else
         {
            *myrow = *myrow + 1;
         }
      }
      else
      {
         *mycol = *mycol + 1;
      }
   }
   else            /* col major */
   {
      if( *myrow == nprow-1 )
      {
         *myrow = 0;
         if( *mycol == npcol-1 )
         {
            *mycol = 0;
         }
         else
         {
            *mycol = *mycol + 1;
         }
      }
      else
      {
         *myrow = *myrow + 1;
      }
   }
}
