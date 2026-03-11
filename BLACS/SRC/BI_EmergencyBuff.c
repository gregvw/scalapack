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
   size_t active_ops, i, j;
   double Mwalltime(void);
   double t1;
   extern Int BI_Np;
   extern BLACBUFF *BI_ReadyB, *BI_ActiveQ;

   if (!ScaLAPACK_Index64ToSizeT((ScaLAPACK_Index64) BI_Np, &active_ops) ||
       active_ops > (SIZE_MAX - sizeof(BLACBUFF)) / sizeof(MPI_Request))
      BI_BlacsErr(-1, __LINE__, __FILE__, "BLACS buffer metadata overflow");

   j = sizeof(BLACBUFF);
   if (j % sizeof(MPI_Request))
      j += sizeof(MPI_Request) - j % sizeof(MPI_Request);
   i = j + active_ops * sizeof(MPI_Request);
   if (i % BUFFALIGN) i += BUFFALIGN - i % BUFFALIGN;
   t1 =  Mwalltime();
   while ( (BI_ActiveQ) && (Mwalltime() - t1 < BUFWAIT) && !(BI_ReadyB) )
   {
      BI_UpdateBuffs(NULL);
      if (BI_ReadyB)
      {
         if (BI_ReadyB->BufLen < length)
         {
	    free(BI_ReadyB);
            if ((size_t) length > SIZE_MAX - i)
               BI_BlacsErr(-1, __LINE__, __FILE__,
                           "BLACS buffer allocation overflow");
            cptr = malloc((size_t) length + i);
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
