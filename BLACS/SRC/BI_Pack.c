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
   size_t header_bytes, aops_count, aops_bytes, total_bytes;
#endif

/*
 * Some versions of mpich and its derivitives cannot handle 0 byte typedefs,
 * so we have set MPI_BYTE as a flag for a 0 byte message
 */
#ifdef ZeroByteTypeBug
   if (Dtype == MPI_BYTE)
   {
      header_bytes = sizeof(BLACBUFF);
      if (!ScaLAPACK_SizeTAlignUp(header_bytes, sizeof(MPI_Request), &header_bytes) ||
          !ScaLAPACK_Index64ToSizeT((ScaLAPACK_Index64) BI_Np, &aops_count) ||
          !ScaLAPACK_SizeTMul(aops_count, sizeof(MPI_Request), &aops_bytes) ||
          !ScaLAPACK_SizeTAdd(header_bytes, aops_bytes, &total_bytes) ||
          !ScaLAPACK_SizeTAlignUp(total_bytes, BUFFALIGN, &total_bytes))
         BI_BlacsErr(BI_ContxtNum(ctxt), __LINE__, __FILE__,
                     "0 byte buffer workspace overflow");
      cptr = malloc(total_bytes);
      if (cptr)
      {
         bp = (BLACBUFF *) cptr;
         bp->BufLen = 0;
         bp->Len = bp->nAops = 0;
         bp->N = 0;
         bp->Aops = (MPI_Request *) &cptr[header_bytes];
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
