#include "Bdef.h"

/***************************************************************************
 *  If there is insufficient space to allocate a needed buffer, this       *
 *  routine is called.  It moniters active buffers for the time defined by *
 *  the user-changeable macro value BUFWAIT.  If in that time no active    *
 *  buffer becomes inactive, a hang is assumed, and the grid is killed.    *
 ***************************************************************************/
void BI_EmergencyBuff(Int length)
{
   BI_EmergencyBuffS((ScaLAPACK_BufLen) length);
}

void BI_EmergencyBuffS(ScaLAPACK_BufLen length)
{
   void BI_UpdateBuffs(BLACBUFF *);

   char *cptr;
   size_t active_ops, i, j, padding;
   double Mwalltime(void);
   double t1;
   extern Int BI_Np;
   extern BLACBUFF *BI_ReadyB, *BI_ActiveQ;

   if (!ScaLAPACK_Index64ToSizeT((ScaLAPACK_Index64) BI_Np, &active_ops) ||
       active_ops > (SIZE_MAX - sizeof(BLACBUFF)) / sizeof(MPI_Request))
      BI_BlacsErr(-1, __LINE__, __FILE__, "BLACS buffer metadata overflow");

   j = sizeof(BLACBUFF);
   if (j % sizeof(MPI_Request))
      if (!ScaLAPACK_SizeTAdd(j, sizeof(MPI_Request) - j % sizeof(MPI_Request), &j))
         BI_BlacsErr(-1, __LINE__, __FILE__, "BLACS buffer metadata overflow");
   if (!ScaLAPACK_SizeTMul(active_ops, sizeof(MPI_Request), &padding) ||
       !ScaLAPACK_SizeTAdd(j, padding, &i))
      BI_BlacsErr(-1, __LINE__, __FILE__, "BLACS buffer metadata overflow");
   if (i % BUFFALIGN)
   {
      padding = BUFFALIGN - i % BUFFALIGN;
      if (!ScaLAPACK_SizeTAdd(i, padding, &i))
         BI_BlacsErr(-1, __LINE__, __FILE__, "BLACS buffer metadata overflow");
   }
   t1 =  Mwalltime();
   while ( (BI_ActiveQ) && (Mwalltime() - t1 < BUFWAIT) && !(BI_ReadyB) )
   {
      BI_UpdateBuffs(NULL);
      if (BI_ReadyB)
      {
         if (BI_ReadyB->BufLen < length)
         {
            free(BI_ReadyB);
            if (!ScaLAPACK_SizeTAdd(i, (size_t) length, &padding))
               BI_BlacsErr(-1, __LINE__, __FILE__,
                           "BLACS buffer allocation overflow");
            cptr = malloc(padding);
            BI_ReadyB = (BLACBUFF *) cptr;
            if (BI_ReadyB)
            {
               BI_ReadyB->BufLen = length;
               BI_ReadyB->Len = 0;
               BI_ReadyB->nAops = 0;
               BI_ReadyB->Aops = (MPI_Request *) &cptr[j];
               BI_ReadyB->Buff = &cptr[i];
            }
         }
      }
   }
   if (BI_ReadyB == NULL)
   {
      BI_BlacsErr(-1, __LINE__, __FILE__, "BLACS out of buffer space");
   }
}
