#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cfree_blacs_system_handle(Int ISysCtxt)
#else
void free_blacs_system_handle_(Int *ISysCxt)
#endif
{
#if (INTFACE == C_CALL)
   Int i, j, DEF_WORLD;
   size_t alloc_elems, alloc_bytes;
   MPI_Comm *tSysCtxt;
   extern Int BI_MaxNSysCtxt;
   extern MPI_Comm *BI_SysContxts;


   if ( (ISysCtxt < BI_MaxNSysCtxt) && (ISysCtxt > 0) )
   {
      if (BI_SysContxts[ISysCtxt] != MPI_COMM_NULL)
         BI_SysContxts[ISysCtxt] = MPI_COMM_NULL;
      else BI_BlacsWarn(-1, __LINE__, __FILE__,
          "Trying to free non-existent system context handle %d", ISysCtxt);
   }
   else if (ISysCtxt == 0) return;  /* never free MPI_COMM_WORLD */
   else BI_BlacsWarn(-1, __LINE__, __FILE__,
        "Trying to free non-existent system context handle %d", ISysCtxt);

/*
 * See if we have freed enough space to decrease the size of our table
 */
   for (i=j=0; i < BI_MaxNSysCtxt; i++)
      if (BI_SysContxts[i] == MPI_COMM_NULL) j++;
/*
 * If needed, get a smaller system context array
 */
   if (j > 2*MAXNSYSCTXT)
   {
      j = BI_MaxNSysCtxt - MAXNSYSCTXT;
      if (!ScaLAPACK_Index64ToSizeT((ScaLAPACK_Index64) j, &alloc_elems) ||
          !ScaLAPACK_SizeTMul(alloc_elems, sizeof(MPI_Comm), &alloc_bytes))
         BI_BlacsWarn(-1, __LINE__, __FILE__,
             "System context table resize overflow");
      tSysCtxt = (MPI_Comm *) malloc(alloc_bytes);
      if (!tSysCtxt)
         BI_BlacsWarn(-1, __LINE__, __FILE__,
             "Cannot shrink system context table");
      for (i=j=0; i < BI_MaxNSysCtxt; i++)
      {
         if (BI_SysContxts[i] != MPI_COMM_NULL)
            tSysCtxt[j++] = BI_SysContxts[i];
      }
      BI_MaxNSysCtxt -= MAXNSYSCTXT;
      for(; j < BI_MaxNSysCtxt; j++) tSysCtxt[j] = MPI_COMM_NULL;
      free(BI_SysContxts);
      BI_SysContxts = tSysCtxt;
   }
#endif
}
