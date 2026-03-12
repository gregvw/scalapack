#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_gridinit(Int *ConTxt, char *order, Int nprow, Int npcol)
#else
F_VOID_FUNC blacs_gridinit_(Int *ConTxt, F_CHAR order, Int *nprow, Int *npcol)
#endif
{
#if (INTFACE == C_CALL)
   void Cblacs_gridmap(Int *, Int *, Int, Int, Int);
#else
   F_VOID_FUNC blacs_gridmap_(Int *ConTxt, Int *usermap, Int *ldup, Int *nprow0,
                              Int *npcol0);
#endif
   ScaLAPACK_Index64 grid_elems64;
   size_t grid_bytes;
   Int *tmpgrid, *iptr;
   Int i, j;

/*
 * Grid can be row- or column-major natural ordering when blacs_gridinit is
 * called.  Define a tmpgrid to reflect this, and call blacs_gridmap to
 * set it up
 */
   if (!ScaLAPACK_Index64Mul((ScaLAPACK_Index64) Mpval(nprow),
                             (ScaLAPACK_Index64) Mpval(npcol),
                             &grid_elems64) ||
       !ScaLAPACK_Index64ToSizeT(grid_elems64, &grid_bytes) ||
       !ScaLAPACK_SizeTMul(grid_bytes, sizeof(*tmpgrid), &grid_bytes))
      BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDINIT",
                  "Temporary grid allocation overflow");

   iptr = tmpgrid = (Int*) malloc(grid_bytes);
   if (tmpgrid == NULL)
      BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDINIT",
                  "Cannot allocate temporary grid");
   if (Mlowcase(F2C_CharTrans(order)) == 'c')
   {
      i = Mpval(npcol) * Mpval(nprow);
      for (j=0; j < i; j++) iptr[j] = j;
   }
   else
   {
      for (j=0; j < Mpval(npcol); j++)
      {
         for (i=0; i < Mpval(nprow); i++) iptr[i] = i * Mpval(npcol) + j;
         iptr += Mpval(nprow);
      }
   }
#if (INTFACE == C_CALL)
   Cblacs_gridmap(ConTxt, tmpgrid, nprow, nprow, npcol);
#else
   blacs_gridmap_(ConTxt, tmpgrid, nprow, nprow, npcol);
#endif
   free(tmpgrid);
}
