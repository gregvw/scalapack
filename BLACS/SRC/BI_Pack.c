#include "Bdef.h"
BLACBUFF *BI_Pack(BLACSCONTEXT *ctxt,BVOID *A,BLACBUFF *bp,MPI_Datatype Dtype)
{
   ScaLAPACK_BufLen i = 0;
   int info;
   MpiInt one=1;
#ifdef ZeroByteTypeBug
   char *cptr;
   extern BLACBUFF BI_AuxBuff;
   extern Int BI_Np;
#endif

/*
 * Some versions of mpich and its derivitives cannot handle 0 byte typedefs,
 * so we have set MPI_BYTE as a flag for a 0 byte message
 */
#ifdef ZeroByteTypeBug
   if (Dtype == MPI_BYTE)
   {
      info = sizeof(BLACBUFF);
      if (info % sizeof(MPI_Request))
         info += sizeof(MPI_Request) - info % sizeof(MPI_Request);
      i = info + BI_Np*sizeof(MPI_Request);
      if (i % BUFFALIGN) i += BUFFALIGN - i % BUFFALIGN;
      cptr = malloc(i);
      if (cptr)
      {
         bp = (BLACBUFF *) cptr;
         bp->BufLen = 0;
         bp->Len = bp->nAops = 0;
         bp->N = 0;
         bp->Aops = (MPI_Request *) &cptr[info];
         bp->Buff = (char *) &bp->BufLen;
         bp->dtype = MPI_BYTE;
         return(bp);
      }
      else BI_BlacsErr(BI_ContxtNum(ctxt), __LINE__, __FILE__, 
                       "Not enough memory to allocate 0 byte buffer\n");
   }
#endif
   if (bp == NULL)
   {
      info = _MPI_Pack_size(one, Dtype, ctxt->scp->comm, &i);
      bp = BI_GetBuffS(i);
   }

   i = 0;
   info = _MPI_Pack(A, one, Dtype, bp->Buff, bp->BufLen, &i, ctxt->scp->comm);
   bp->dtype = MPI_PACKED;
   bp->N = (MpiInt) i;

   return(bp);
}
