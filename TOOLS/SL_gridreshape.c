#include <stdlib.h>
#include <stdio.h>
#include "scalapack-types.h"

#ifndef Int
#define Int int
#endif

_Static_assert(sizeof(Int) == sizeof(ScaLAPACK_ApiInt),
               "Int must match the configured ScaLAPACK API integer width.");

Int SL_Cgridreshape(Int ctxt, Int pstart, Int row_major_in, Int row_major_out, Int P, Int Q)
{
   void Cblacs_gridinfo( Int context, Int *nprow, Int *npcol, Int *myrow, Int *mycol );
   void Cblacs_abort( Int context, Int errornum );
   void Cblacs_get( Int context, Int what, Int *val );
   void Cblacs_gridmap( Int *context, Int *usermap, Int ldumap, Int nprow, Int npcol );
   Int Cblacs_pnum( Int context, Int prow, Int pcol );
   ScaLAPACK_ApiInt prow, pcol;
   ScaLAPACK_Index64 Np, grid_index, limit, source_rank, user_rank;
   size_t grid_elems, grid_slot;
   Int nctxt, P0, Q0, mycol, myrow, *g;

   Cblacs_gridinfo(ctxt, &P0, &Q0, &myrow, &mycol);
   if (!ScaLAPACK_Index64Mul((ScaLAPACK_Index64) P, (ScaLAPACK_Index64) Q, &Np) ||
       !ScaLAPACK_Index64Mul((ScaLAPACK_Index64) P0, (ScaLAPACK_Index64) Q0, &limit) ||
       !ScaLAPACK_Index64Add((ScaLAPACK_Index64) pstart, Np, &source_rank) ||
       source_rank > limit)
   {
      fprintf(stderr, "Illegal reshape command in %s\n",__FILE__);
      Cblacs_abort(ctxt, -22);
      return -1;
   }
   if (!ScaLAPACK_Index64ToSizeT(Np, &grid_elems) ||
       grid_elems > SIZE_MAX / sizeof(Int))
   {
      fprintf(stderr, "Illegal reshape command in %s\n",__FILE__);
      Cblacs_abort(ctxt, -22);
      return -1;
   }
   g = (Int *) malloc(grid_elems * sizeof(Int));
   if (!g)
   {
      fprintf(stderr, "Cannot allocate memory in %s\n",__FILE__);
      Cblacs_abort(ctxt, -23);
      return -1;
   }
   if (row_major_in)  /* Read in in row-major order */
   {
      if (row_major_out)
         for (user_rank = 0; user_rank != Np; ++user_rank)
         {
            if (!ScaLAPACK_Index64Mul(user_rank % Q, P, &grid_index) ||
                !ScaLAPACK_Index64Add(grid_index, user_rank / Q, &grid_index) ||
                !ScaLAPACK_Index64ToSizeT(grid_index, &grid_slot) ||
                !ScaLAPACK_Index64Add((ScaLAPACK_Index64) pstart, user_rank, &source_rank) ||
                !ScaLAPACK_Index64ToApiInt(source_rank / Q0, &prow) ||
                !ScaLAPACK_Index64ToApiInt(source_rank % Q0, &pcol))
            {
               free(g);
               fprintf(stderr, "Illegal reshape command in %s\n", __FILE__);
               Cblacs_abort(ctxt, -22);
               return -1;
            }
            g[grid_slot] = Cblacs_pnum(ctxt, prow, pcol);
         }
      else
         for (user_rank = 0; user_rank != Np; ++user_rank)
         {
            if (!ScaLAPACK_Index64ToSizeT(user_rank, &grid_slot) ||
                !ScaLAPACK_Index64Add((ScaLAPACK_Index64) pstart, user_rank, &source_rank) ||
                !ScaLAPACK_Index64ToApiInt(source_rank / Q0, &prow) ||
                !ScaLAPACK_Index64ToApiInt(source_rank % Q0, &pcol))
            {
               free(g);
               fprintf(stderr, "Illegal reshape command in %s\n", __FILE__);
               Cblacs_abort(ctxt, -22);
               return -1;
            }
            g[grid_slot] = Cblacs_pnum(ctxt, prow, pcol);
         }
   }
   else /* read in in column-major order */
   {
      if (row_major_out)
         for (user_rank = 0; user_rank != Np; ++user_rank)
         {
            if (!ScaLAPACK_Index64Mul(user_rank % Q, P, &grid_index) ||
                !ScaLAPACK_Index64Add(grid_index, user_rank / Q, &grid_index) ||
                !ScaLAPACK_Index64ToSizeT(grid_index, &grid_slot) ||
                !ScaLAPACK_Index64Add((ScaLAPACK_Index64) pstart, user_rank, &source_rank) ||
                !ScaLAPACK_Index64ToApiInt(source_rank % P0, &prow) ||
                !ScaLAPACK_Index64ToApiInt(source_rank / P0, &pcol))
            {
               free(g);
               fprintf(stderr, "Illegal reshape command in %s\n", __FILE__);
               Cblacs_abort(ctxt, -22);
               return -1;
            }
            g[grid_slot] = Cblacs_pnum(ctxt, prow, pcol);
         }
      else
         for (user_rank = 0; user_rank != Np; ++user_rank)
         {
            if (!ScaLAPACK_Index64ToSizeT(user_rank, &grid_slot) ||
                !ScaLAPACK_Index64Add((ScaLAPACK_Index64) pstart, user_rank, &source_rank) ||
                !ScaLAPACK_Index64ToApiInt(source_rank % P0, &prow) ||
                !ScaLAPACK_Index64ToApiInt(source_rank / P0, &pcol))
            {
               free(g);
               fprintf(stderr, "Illegal reshape command in %s\n", __FILE__);
               Cblacs_abort(ctxt, -22);
               return -1;
            }
            g[grid_slot] = Cblacs_pnum(ctxt, prow, pcol);
         }
   }
   Cblacs_get(ctxt, 10, &nctxt);
   Cblacs_gridmap(&nctxt, g, P, P, Q);
   free(g);

   return(nctxt);
}

Int sl_gridreshape_(Int *ctxt, Int *pstart, Int *row_major_in, Int *row_major_out, Int *P, Int *Q)
{
   return( SL_Cgridreshape(*ctxt, *pstart, *row_major_in, *row_major_out,
                           *P, *Q) );
}

Int SL_GRIDRESHAPE(Int *ctxt, Int *pstart, Int *row_major_in, Int *row_major_out, Int *P, Int *Q)
{
   return( SL_Cgridreshape(*ctxt, *pstart, *row_major_in, *row_major_out,
                           *P, *Q) );
}

Int sl_gridreshape__(Int *ctxt, Int *pstart, Int *row_major_in, Int *row_major_out, Int *P, Int *Q)
{
   return( SL_Cgridreshape(*ctxt, *pstart, *row_major_in, *row_major_out,
                           *P, *Q) );
}

Int sl_gridreshape(Int *ctxt, Int *pstart, Int *row_major_in, Int *row_major_out, Int *P, Int *Q)
{
   return( SL_Cgridreshape(*ctxt, *pstart, *row_major_in, *row_major_out,
                           *P, *Q) );
}
