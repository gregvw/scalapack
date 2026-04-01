#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_abort(Int ConTxt, Int ErrNo)
#else
F_VOID_FUNC blacs_abort_(Int *ConTxt, Int *ErrNo)
#endif
{
   void Cblacs_gridinfo(Int, Int *, Int *, Int *, Int *);
   void BI_BlacsAbort(Int ErrNo);
   Int nprow, npcol, myrow, mycol;
   extern Int BI_Iam;

   Cblacs_gridinfo(Mpval(ConTxt), &nprow, &npcol, &myrow, &mycol);
   fprintf(stderr,
"{" SLINT_FMT "," SLINT_FMT "}, pnum=" SLINT_FMT ", Contxt=" SLINT_FMT ", killed other procs, exiting with error #" SLINT_FMT ".\n\n",
           myrow, mycol, BI_Iam, Mpval(ConTxt), Mpval(ErrNo));

   BI_BlacsAbort(Mpval(ErrNo));
}
