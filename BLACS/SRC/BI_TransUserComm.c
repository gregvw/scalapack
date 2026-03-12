#include "Bdef.h"

MPI_Comm BI_TransUserComm(Int Ucomm, int Np, Int *pmap)
{
   MPI_Comm bcomm, ucomm;
   MPI_Group bgrp, ugrp;
   Int i;
   size_t map_bytes;

   if (!ScaLAPACK_Index64ToSizeT((ScaLAPACK_Index64) Np, &map_bytes) ||
       !ScaLAPACK_SizeTMul(map_bytes, sizeof(int), &map_bytes))
      BI_BlacsErr(-1, __LINE__, __FILE__, "MPI rank map allocation overflow");

   int *mpmap = (int *)malloc(map_bytes);
   if (mpmap == NULL)
      BI_BlacsErr(-1, __LINE__, __FILE__, "Cannot allocate MPI rank map");
   for (i=0; i<Np; i++)
   {
      if (!ScaLAPACK_ApiIntToCInt((ScaLAPACK_ApiInt) pmap[i], &mpmap[i]))
      {
         free(mpmap);
         BI_BlacsErr(-1, __LINE__, __FILE__, "MPI rank map value out of range");
      }
   }

   ucomm = MPI_Comm_f2c(Ucomm);
   i=MPI_Comm_group(ucomm, &ugrp);
   i=MPI_Group_incl(ugrp, Np, mpmap, &bgrp);
   i=MPI_Comm_create(ucomm, bgrp, &bcomm);
   i=MPI_Group_free(&ugrp);
   i=MPI_Group_free(&bgrp);

   free(mpmap);

   return(bcomm);
}
