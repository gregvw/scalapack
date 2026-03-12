#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_gridmap(Int *ConTxt, Int *usermap, Int ldup, Int nprow0, Int npcol0)
#else
F_VOID_FUNC blacs_gridmap_(Int *ConTxt, Int *usermap, Int *ldup, Int *nprow0,
                           Int *npcol0)
#endif
{
   void Cblacs_pinfo(Int *, Int *);
   void Cblacs_get(Int, Int, Int *);

   MPI_Comm Cblacs2sys_handle(Int BlacsCtxt);
   MPI_Comm BI_TransUserComm(Int, Int, Int *);

   ScaLAPACK_ApiInt ng_api;
   int Iam;
   ScaLAPACK_Index64 ng64;
   size_t alloc_elems, alloc_bytes;
   Int info, i, j, *iptr;
   Int myrow, mycol, nprow, npcol, Ng;
   BLACSCONTEXT *ctxt, **tCTxts;
   MPI_Comm comm, tcomm;
   MPI_Group grp, tgrp;

   extern BLACSCONTEXT **BI_MyContxts;
   extern BLACBUFF BI_AuxBuff;
   extern Int BI_Iam, BI_Np, BI_MaxNCtxt;
   extern MPI_Status *BI_Stats;

/*
 * If first call to blacs_gridmap
 */
   if (BI_MaxNCtxt == 0)
   {
      Cblacs_pinfo(&BI_Iam, &BI_Np);
      BI_AuxBuff.nAops = 0;
      if (!ScaLAPACK_Index64ToSizeT((ScaLAPACK_Index64) BI_Np, &alloc_elems) ||
          !ScaLAPACK_SizeTMul(alloc_elems, sizeof(*BI_AuxBuff.Aops), &alloc_bytes))
         BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                     "Auxiliary request allocation overflow");
      BI_AuxBuff.Aops = (MPI_Request*)malloc(alloc_bytes);
      if (!ScaLAPACK_SizeTMul(alloc_elems, sizeof(MPI_Status), &alloc_bytes))
         BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                     "Auxiliary status allocation overflow");
      BI_Stats = (MPI_Status *) malloc(alloc_bytes);
      if (BI_AuxBuff.Aops == NULL || BI_Stats == NULL)
         BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                     "Cannot allocate auxiliary BLACS buffers");
   }

   nprow = Mpval(nprow0);
   npcol = Mpval(npcol0);
   if (!ScaLAPACK_Index64Mul((ScaLAPACK_Index64) nprow,
                             (ScaLAPACK_Index64) npcol, &ng64) ||
       !ScaLAPACK_Index64ToApiInt(ng64, &ng_api))
      BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDINIT/BLACS_GRIDMAP",
                  "Illegal grid (%d x %d), #procs=%d", nprow, npcol, BI_Np);
   Ng = (Int) ng_api;
   if ( (Ng > BI_Np) || (nprow < 1) || (npcol < 1) )
      BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDINIT/BLACS_GRIDMAP",
                  "Illegal grid (%d x %d), #procs=%d", nprow, npcol, BI_Np);
/*
 * Form MPI communicator for scope = 'all'
 */
   if (Ng > 2) i = Ng;
   else i = 2;
   if (!ScaLAPACK_Index64ToSizeT((ScaLAPACK_Index64) i, &alloc_elems) ||
       !ScaLAPACK_SizeTMul(alloc_elems, sizeof(Int), &alloc_bytes))
      BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                  "User rank map allocation overflow");
   iptr = (Int *) malloc(alloc_bytes);
   if (iptr == NULL)
      BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                  "Cannot allocate user rank map");
   for (j=0; j < npcol; j++)
   {
      for (i=0; i < nprow; i++) iptr[i*npcol+j] = usermap[j*Mpval(ldup)+i];
   }
#if (INTFACE == C_CALL)
   if (!ScaLAPACK_Index64ToSizeT((ScaLAPACK_Index64) Ng, &alloc_elems) ||
       !ScaLAPACK_SizeTMul(alloc_elems, sizeof(int), &alloc_bytes))
      BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                  "MPI rank map allocation overflow");
   int *miptr = (int *) malloc(alloc_bytes);
   if (miptr == NULL)
      BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                  "Cannot allocate MPI rank map");
   for (j=0; j < Ng; j++)
   {
      if (!ScaLAPACK_ApiIntToCInt((ScaLAPACK_ApiInt) iptr[j], &miptr[j]))
      {
         free(miptr);
         free(iptr);
         BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                     "MPI rank map value out of range");
      }
   }
   tcomm = Cblacs2sys_handle(*ConTxt);
   MPI_Comm_group(tcomm, &grp);           /* find input comm's group */
   MPI_Group_incl(grp, Ng, miptr, &tgrp);  /* form new group */
   MPI_Comm_create(tcomm, tgrp, &comm);   /* create new comm */
   MPI_Group_free(&tgrp);
   MPI_Group_free(&grp);
   free(miptr);
#else  /* gridmap called from fortran */
   comm = BI_TransUserComm(*ConTxt, Ng, iptr);
#endif

/*
 * Weed out callers who are not participating in present grid
 */
   if (comm == MPI_COMM_NULL)
   {
      *ConTxt = NOTINCONTEXT;
      free(iptr);
      return;
   }

/*
 * ==================================================
 * Get new context and add it to my array of contexts
 * ==================================================
 */
   ctxt = (BLACSCONTEXT *) malloc(sizeof(BLACSCONTEXT));
   if (ctxt == NULL)
      BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                  "Cannot allocate BLACS context");
/*
 * Find free slot in my context array
 */
   for (i=0; i < BI_MaxNCtxt; i++) if (BI_MyContxts[i] == NULL) break;
/*
 * Get bigger context pointer array, if needed
 */
   if (i == BI_MaxNCtxt)
   {
      j = BI_MaxNCtxt + MAXNCTXT;
      if (!ScaLAPACK_Index64ToSizeT((ScaLAPACK_Index64) j, &alloc_elems) ||
          !ScaLAPACK_SizeTMul(alloc_elems, sizeof(*tCTxts), &alloc_bytes))
         BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                     "Context table allocation overflow");
      tCTxts = (BLACSCONTEXT **) malloc(alloc_bytes);
      if (tCTxts == NULL)
         BI_BlacsErr((Int)-1, (Int)-1, "BLACS_GRIDMAP",
                     "Cannot allocate context table");
      for (i=0; i < BI_MaxNCtxt; i++) tCTxts[i] = BI_MyContxts[i];
      BI_MaxNCtxt = j;
      for(j=i; j < BI_MaxNCtxt; j++) tCTxts[j] = NULL;
      if (BI_MyContxts) free(BI_MyContxts);
      BI_MyContxts = tCTxts;
   }
   BI_MyContxts[i] = ctxt;
   *ConTxt = i;

   ctxt->ascp.comm = comm;
   MPI_Comm_dup(comm, &ctxt->pscp.comm); /* copy acomm for pcomm */
   MPI_Comm_rank(comm, &Iam);            /* find my rank in new comm */
   myrow = Iam / npcol;
   mycol = Iam % npcol;

/*
 * Form MPI communicators for scope = 'row'
 */
   MPI_Comm_split(comm, myrow, mycol, &ctxt->rscp.comm);
/*
 * Form MPI communicators for scope = 'Column'
 */
   MPI_Comm_split(comm, mycol, myrow, &ctxt->cscp.comm);

   ctxt->rscp.Np = npcol;
   ctxt->rscp.Iam = mycol;
   ctxt->cscp.Np = nprow;
   ctxt->cscp.Iam = myrow;
   ctxt->pscp.Np = ctxt->ascp.Np = Ng;
   ctxt->pscp.Iam = ctxt->ascp.Iam = Iam;
   ctxt->Nr_bs = ctxt->Nr_co = 1;
   ctxt->Nb_bs = ctxt->Nb_co = 2;
   ctxt->TopsRepeat = ctxt->TopsCohrnt = 0;

/*
 * ===========================
 * Set up the message id stuff
 * ===========================
 */
   Cblacs_get(-1, 1, iptr);
   ctxt->pscp.MinId = ctxt->rscp.MinId = ctxt->cscp.MinId =
   ctxt->ascp.MinId = ctxt->pscp.ScpId = ctxt->rscp.ScpId =
   ctxt->cscp.ScpId = ctxt->ascp.ScpId = iptr[0];
   ctxt->pscp.MaxId = ctxt->rscp.MaxId = ctxt->cscp.MaxId =
   ctxt->ascp.MaxId = iptr[1];
   free(iptr);

}
