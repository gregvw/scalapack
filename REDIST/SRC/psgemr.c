#include "redist.h"
#include "redist_core.h"
#include <stddef.h>
/** $Id: psgemr.c,v 1.1.1.1 2000/02/15 18:04:09 susan Exp $
  ------------------------------------------------------------------------

    -- ScaLAPACK routine (version 1.7) --
       Oak Ridge National Laboratory, Univ. of Tennessee, and Univ. of
       California, Berkeley.
       October 31, 1994.

      SUBROUTINE PSGEMR2D( M, N,
     $                     A, IA, JA, ADESC,
     $                     B, IB, JB, BDESC,
     $                     CTXT)
  ------------------------------------------------------------------------
    Purpose
    =======

    PSGEMR2D copies a submatrix of A on a submatrix of B.
    A and B can have different distributions: they can be on different
    processor grids, they can have different blocksizes, the beginning
    of the area to be copied can be at a different places on A and B.

    The parameters can be confusing when the grids of A and B are
    partially or completly disjoint, in the case a processor calls
    this routines but is either not in the A context or B context, the
    ADESC[CTXT] or BDESC[CTXT] must be equal to -1, to ensure the
    routine recognise this situation.
    To summarize the rule:
    - If a processor is in A context, all parameters related to A must be valid.
    - If a processor is in B context, all parameters related to B must be valid.
    -  ADESC[CTXT] and BDESC[CTXT] must be either valid contexts or equal to -1.
    - M and N must be valid for everyone.
    - other parameters are not examined.


    Notes
    =====

    A description vector is associated with each 2D block-cyclicly dis-
    tributed matrix.  This vector stores the information required to
    establish the mapping between a matrix entry and its corresponding
    process and memory location.

    In the following comments, the character _ should be read as
    "of the distributed matrix".  Let A be a generic term for any 2D
    block cyclicly distributed matrix.  Its description vector is DESC_A:

   NOTATION        STORED IN      EXPLANATION
   --------------- -------------- --------------------------------------
   DT_A   (global) DESCA( DT_ )   The descriptor type.
   CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
                                  the BLACS process grid A is distribu-
                                  ted over. The context itself is glo-
                                  bal, but the handle (the integer
                                  value) may vary.
   M_A    (global) DESCA( M_ )    The number of rows in the distributed
                                  matrix A.
   N_A    (global) DESCA( N_ )    The number of columns in the distri-
                                  buted matrix A.
   MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
                                  the rows of A.
   NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
                                  the columns of A.
   RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
                                  row of the matrix A is distributed.
   CSRC_A (global) DESCA( CSRC_ ) The process column over which the
                                  first column of A is distributed.
   LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
                                  array storing the local blocks of the
                                  distributed matrix A.
                                  LLD_A >= MAX(1,LOCp(M_A)).



    Important notice
    ================
     The parameters of the routine have changed in April 1996
     There is a new last argument. It must be a context englobing
     all processors involved in the initial and final distribution.

     Be aware that all processors  included in this
      context must call the redistribution routine.

    Parameters
    ==========


    M        (input) INTEGER.
             On entry, M specifies the number of rows of the
             submatrix to be copied.  M must be at least zero.
             Unchanged on exit.

    N        (input) INTEGER.
             On entry, N specifies the number of cols of the submatrix
             to be redistributed.rows of B.  M must be at least zero.
             Unchanged on exit.

    A        (input) REAL
             On entry, the source matrix.
             Unchanged on exit.

    IA,JA    (input) INTEGER
             On entry,the coordinates of the beginning of the submatrix
             of A to copy.
             1 <= IA <= M_A - M + 1,1 <= JA <= N_A - N + 1,
             Unchanged on exit.

    ADESC    (input) A description vector (see Notes above)
             If the current processor is not part of the context of A
             the ADESC[CTXT] must be equal to -1.


    B        (output) REAL
             On entry, the destination matrix.
             The portion corresponding to the defined submatrix are updated.

    IB,JB    (input) INTEGER
             On entry,the coordinates of the beginning of the submatrix
             of B that will be updated.
             1 <= IB <= M_B - M + 1,1 <= JB <= N_B - N + 1,
             Unchanged on exit.

    BDESC    (input) B description vector (see Notes above)
             For processors not part of the context of B
             BDESC[CTXT] must be equal to -1.

    CTXT     (input) a context englobing at least all processors included
                in either A context or B context



   Memory requirement :
   ====================

   for the processors belonging to grid 0, one buffer of size block 0
   and for the processors belonging to grid 1, also one buffer of size
   block 1.

   ============================================================
   Created March 1993 by B. Tourancheau (See sccs for modifications).
   Modifications by Loic PRYLLI 1995
   ============================================================ */
#define static2 static
#if defined(Add_) || defined(f77IsF2C)
#define fortran_mr2d psgemr2do_
#define fortran_mr2dnew psgemr2d_
#elif defined(UpCase)
#define fortran_mr2dnew PSGEMR2D
#define fortran_mr2d PSGEMR2DO
#define scopy_ SCOPY
#define slacpy_ SLACPY
#else
#define fortran_mr2d psgemr2do
#define fortran_mr2dnew psgemr2d
#define scopy_ scopy
#define slacpy_ slacpy
#endif
/* I8 entry point name mangling */
#if defined(Add_) || defined(f77IsF2C)
#define fortran_mr2dnew_i8 psgemr2d_i8_
#elif defined(UpCase)
#define fortran_mr2dnew_i8 PSGEMR2D_I8
#else
#define fortran_mr2dnew_i8 psgemr2d_i8
#endif
#define Clacpy Csgelacpy
void  Clacpy( Int m, Int n, float *a, Int lda, float *b, Int ldb );
typedef struct {
  Int   desctype;
  Int   ctxt;
  Int   m;
  Int   n;
  Int   nbrow;
  Int   nbcol;
  Int   sprow;
  Int   spcol;
  Int   lda;
}     MDESC;
#define BLOCK_CYCLIC_2D 1
typedef struct {
  Int   lstart;
  Int   len;
}     IDESC;
#define SHIFT(row,sprow,nbrow) ((row)-(sprow)+ ((row) >= (sprow) ? 0 : (nbrow)))
#define max(A,B) ((A)>(B)?(A):(B))
#define min(A,B) ((A)>(B)?(B):(A))
#define DIVUP(a,b) ( ((a)-1) /(b)+1)
#define ROUNDUP(a,b) (DIVUP(a,b)*(b))
#ifdef MALLOCDEBUG
#define malloc mymalloc
#define free myfree
#define realloc myrealloc
#endif
/* Cblacs */
extern void Cblacs_pcoord( Int context, Int pnum, Int* prow, Int* pcol );
extern Int Cblacs_pnum( Int context, Int prow, Int pcol );
extern void Csetpvmtids();
extern void Cblacs_get( Int context, Int what, Int* val );
extern void Cblacs_pinfo( Int* mypnum, Int* nprocs );
extern void Cblacs_gridinfo( Int context, Int* nprow, Int* npcol, Int* myrow, Int* mycol );
extern void Cblacs_gridinit( Int* context, char* order, Int nprow, Int npcol );
extern void Cblacs_gridmap( Int* context, Int* usermap, Int ldumap, Int nprow, Int npcol );
extern void Cblacs_exit( Int continue_blacs );
extern void Cblacs_gridexit( Int context );
extern void Cblacs_setup( Int* mypnum, Int* nprocs );
extern void Cigebs2d( Int context, char* scope, char* top, Int m, Int n, Int* A, Int lda );
extern void Cigebr2d( Int context, char* scope, char* top, Int m, Int n, Int* A, Int lda, Int rsrc, Int csrc );
extern void Cigesd2d( Int context, Int m, Int n, Int* A, Int lda, Int rdest, Int cdest );
extern void Cigerv2d( Int context, Int m, Int n, Int* A, Int lda, Int rsrc, Int csrc );
extern void Cigsum2d( Int context, char* scope, char* top, Int m, Int n, Int* A, Int lda, Int rdest, Int cdest );
extern void Cigamn2d( Int context, char* scope, char* top, Int m, Int n, Int* A, Int lda, Int* RA, Int* CA, Int rcflag, Int rdest, Int cdest );
extern void Cigamx2d( Int context, char* scope, char* top, Int m, Int n, Int* A, Int lda, Int* RA, Int* CA, Int rcflag, Int rdest, Int cdest );
extern void Csgesd2d( Int context, Int m, Int n, float* A, Int lda, Int rdest, Int cdest );
extern void Csgerv2d( Int context, Int m, Int n, float* A, Int lda, Int rsrc, Int csrc );
/* lapack */
void  slacpy_();
/* aux fonctions */
extern Int localindice( Int ig, Int jg, Int templateheight, Int templatewidth, MDESC *a );
extern void *mr2d_malloc( size_t n );
extern Int ppcm( Int a, Int b );
extern Int localsize( Int myprow, Int p, Int nbrow, Int m );
extern Int memoryblocksize( MDESC *a );
extern Int changeorigin( Int myp, Int sp, Int p, Int bs, Int i, Int *decal, Int *newsp );
extern void paramcheck( MDESC *a, Int i, Int j, Int m, Int n, Int p, Int q, Int gcontext );
/* tools and others function */
#define scanD0 sgescanD0
#define dispmat sgedispmat
#define setmemory sgesetmemory
#define freememory sgefreememory
#define scan_intervals sgescan_intervals
extern void scanD0();
extern void dispmat();
extern void setmemory( float** ptr, Int size );
extern void freememory( char* ptr );
extern Int scan_intervals( char type, Int ja, Int jb, Int n, MDESC *ma, MDESC *mb, Int q0, Int q1, Int col0, Int col1, IDESC *result );
extern void Cpsgemr2do( Int m, Int n, float *ptrmyblock, Int ia, Int ja, MDESC *ma, float *ptrmynewblock, Int ib, Int jb, MDESC *mb );
extern void Cpsgemr2d( Int m, Int n, float *ptrmyblock, Int ia, Int ja, MDESC *ma, float *ptrmynewblock, Int ib, Int jb, MDESC *mb, Int globcontext );
extern void Cpsgemr2d_core( ScaLAPACK_Index64 m, ScaLAPACK_Index64 n, float *ptrmyblock, ScaLAPACK_Index64 ia, ScaLAPACK_Index64 ja, MDESC_CORE *ma, float *ptrmynewblock, ScaLAPACK_Index64 ib, ScaLAPACK_Index64 jb, MDESC_CORE *mb, int globcontext );
/* some defines for Cpsgemr2do */
#define SENDBUFF 0
#define RECVBUFF 1
#define SIZEBUFF 2
#if 0
#define DEBUG
#endif
#ifndef DEBUG
#define NDEBUG
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define DESCLEN 9
void
fortran_mr2d(Int *m, Int *n, float *A, Int *ia, Int *ja, Int desc_A[DESCLEN],
	     float *B, Int *ib, Int *jb, Int desc_B[DESCLEN])
{
  Cpsgemr2do(*m, *n, A, *ia, *ja, (MDESC *) desc_A,
	     B, *ib, *jb, (MDESC *) desc_B);
  return;
}
void
fortran_mr2dnew(Int *m, Int *n, float *A, Int *ia, Int *ja, Int desc_A[DESCLEN],
		float *B, Int *ib, Int *jb, Int desc_B[DESCLEN], Int *gcontext)
{
  Cpsgemr2d(*m, *n, A, *ia, *ja, (MDESC *) desc_A,
	    B, *ib, *jb, (MDESC *) desc_B, *gcontext);
  return;
}
static2 void init_chenille( Int mypnum, Int nprocs, Int n0, Int *proc0, Int n1, Int *proc1, Int **psend, Int **precv, Int *myrang );
static2 Int inter_len( Int hinb, IDESC *hi, Int vinb, IDESC *vi );
static2 Int block2buff( IDESC *vi, Int vinb, IDESC *hi, Int hinb, float *ptra, MDESC *ma, float *buff );
static2 void buff2block( IDESC *vi, Int vinb, IDESC *hi, Int hinb, float *buff, float *ptrb, MDESC *mb );
static2 void gridreshape( Int *ctxtp );
void
Cpsgemr2do(m, n,
	   ptrmyblock, ia, ja, ma,
	   ptrmynewblock, ib, jb, mb)
  float *ptrmyblock, *ptrmynewblock;
/* pointers to the memory location of the matrix and the redistributed matrix */
  MDESC *ma;
  MDESC *mb;
  Int   ia, ja, ib, jb, m, n;
{
  Int   dummy, nprocs;
  Int   gcontext;
  /* first we initialize a global grid which serve as a reference to
   * communicate from grid a to grid b */
  Cblacs_pinfo(&dummy, &nprocs);
  Cblacs_get((Int)0, (Int)0, &gcontext);
  Cblacs_gridinit(&gcontext, "R", (Int)1, nprocs);
  Cpsgemr2d(m, n, ptrmyblock, ia, ja, ma,
	    ptrmynewblock, ib, jb, mb, gcontext);
  Cblacs_gridexit(gcontext);
}
/* Cpsgemr2d: thin wrapper — unpacks MDESC to MDESC_CORE and delegates */
void
Cpsgemr2d(Int m, Int n,
	  float *ptrmyblock, Int ia, Int ja, MDESC *ma,
	  float *ptrmynewblock, Int ib, Int jb, MDESC *mb,
	  Int globcontext)
{
  MDESC_CORE core_a, core_b;
  core_a.desctype = (int)ma->desctype;
  core_a.ctxt     = (int)ma->ctxt;
  core_a.m        = (ScaLAPACK_Index64)ma->m;
  core_a.n        = (ScaLAPACK_Index64)ma->n;
  core_a.mb       = (ScaLAPACK_Index64)ma->nbrow;
  core_a.nb       = (ScaLAPACK_Index64)ma->nbcol;
  core_a.rsrc     = (int)ma->sprow;
  core_a.csrc     = (int)ma->spcol;
  core_a.lld      = (ScaLAPACK_Index64)ma->lda;
  core_b.desctype = (int)mb->desctype;
  core_b.ctxt     = (int)mb->ctxt;
  core_b.m        = (ScaLAPACK_Index64)mb->m;
  core_b.n        = (ScaLAPACK_Index64)mb->n;
  core_b.mb       = (ScaLAPACK_Index64)mb->nbrow;
  core_b.nb       = (ScaLAPACK_Index64)mb->nbcol;
  core_b.rsrc     = (int)mb->sprow;
  core_b.csrc     = (int)mb->spcol;
  core_b.lld      = (ScaLAPACK_Index64)mb->lda;
  Cpsgemr2d_core((ScaLAPACK_Index64)m, (ScaLAPACK_Index64)n,
		 ptrmyblock,
		 (ScaLAPACK_Index64)ia, (ScaLAPACK_Index64)ja, &core_a,
		 ptrmynewblock,
		 (ScaLAPACK_Index64)ib, (ScaLAPACK_Index64)jb, &core_b,
		 (int)globcontext);
}
static2 void
init_chenille(Int mypnum, Int nprocs, Int n0, Int *proc0, Int n1, Int *proc1, Int **psend, Int **precv, Int *myrang)
{
  Int   ns, nr, i, tot;
  size_t alloc_count, alloc_bytes;
  Int  *sender, *recver, *g0, *g1;
  tot = max(n0, n1);
  if (!ScaLAPACK_RedistApiAddToSizeT(nprocs, tot, &alloc_count) ||
      !ScaLAPACK_SizeTMul(alloc_count, sizeof(Int), &alloc_bytes) ||
      !ScaLAPACK_SizeTMul(alloc_bytes, 2, &alloc_bytes)) {
    fprintf(stderr, "xxGEMR2D:sender workspace overflow\n");
    exit(1);
  }
  sender = (Int *) mr2d_malloc(alloc_bytes);
  recver = sender + tot;
  *psend = sender;
  *precv = recver;
  g0 = recver + tot;
  g1 = g0 + nprocs;
  for (i = 0; i < nprocs; i++) {
    g0[i] = -1;
    g1[i] = -1;
  }
  for (i = 0; i < tot; i++) {
    sender[i] = -1;
    recver[i] = -1;
  }
  for (i = 0; i < n0; i++)
    g0[proc0[i]] = i;
  for (i = 0; i < n1; i++)
    g1[proc1[i]] = i;
  ns = 0;
  nr = 0;
  *myrang = -1;
  for (i = 0; i < nprocs; i++)
    if (g0[i] >= 0 && g1[i] >= 0) {
      if (i == mypnum)
	*myrang = nr;
      sender[ns] = g0[i];
      ns += 1;
      recver[nr] = g1[i];
      nr += 1;
      assert(ns <= n0 && nr <= n1 && nr == ns);
    }
  for (i = 0; i < nprocs; i++)
    if (g0[i] >= 0 && g1[i] < 0) {
      if (i == mypnum)
	*myrang = ns;
      sender[ns] = g0[i];
      ns += 1;
      assert(ns <= n0);
    }
  for (i = 0; i < nprocs; i++)
    if (g1[i] >= 0 && g0[i] < 0) {
      if (i == mypnum)
	*myrang = nr;
      recver[nr] = g1[i];
      nr += 1;
      assert(nr <= n1);
    }
}
#define Mlacpy(mo,no,ao,ldao,bo,ldbo) \
{ \
float *_a,*_b; \
Int _m,_n,_lda,_ldb; \
    Int _i,_j; \
    _m = (mo);_n = (no); \
    _a = (ao);_b = (bo); \
    _lda = (ldao) - _m; \
    _ldb = (ldbo) - _m; \
    assert(_lda >= 0 && _ldb >= 0); \
    for (_j=0;_j<_n;_j++) { \
      for (_i=0;_i<_m;_i++) \
        *_b++ = *_a++; \
      _b += _ldb; \
      _a += _lda; \
    } \
} (void)0
static2 Int
block2buff(IDESC *vi, Int vinb, IDESC *hi, Int hinb, float *ptra, MDESC *ma, float *buff)
{
  Int   h, v, sizebuff, block_elems;
  size_t row_offset_elems;
  float *ptr2;
  sizebuff = 0;
  for (h = 0; h < hinb; h++) {
    if (!ScaLAPACK_RedistApiMulToSizeT(hi[h].lstart, ma->lda,
                                       &row_offset_elems)) {
      fprintf(stderr, "xxGEMR2D:block row offset overflow\n");
      exit(1);
    }
    ptr2 = ptra + row_offset_elems;
    for (v = 0; v < vinb; v++) {
      Mlacpy(vi[v].len, hi[h].len,
	     ptr2 + vi[v].lstart,
	     ma->lda,
	     buff + sizebuff, vi[v].len);
      if (!ScaLAPACK_RedistApiMulToApiInt(hi[h].len, vi[v].len, &block_elems) ||
          !ScaLAPACK_RedistApiAddToApiInt(sizebuff, block_elems, &sizebuff)) {
        fprintf(stderr, "xxGEMR2D:block size overflow\n");
        exit(1);
      }
    }
  }
  return sizebuff;
}
static2 void
buff2block(IDESC *vi, Int vinb, IDESC *hi, Int hinb, float *buff, float *ptrb, MDESC *mb)
{
  Int   h, v, sizebuff, block_elems;
  size_t row_offset_elems;
  float *ptr2;
  sizebuff = 0;
  for (h = 0; h < hinb; h++) {
    if (!ScaLAPACK_RedistApiMulToSizeT(hi[h].lstart, mb->lda,
                                       &row_offset_elems)) {
      fprintf(stderr, "xxGEMR2D:block row offset overflow\n");
      exit(1);
    }
    ptr2 = ptrb + row_offset_elems;
    for (v = 0; v < vinb; v++) {
      Mlacpy(vi[v].len, hi[h].len,
	     buff + sizebuff, vi[v].len,
	     ptr2 + vi[v].lstart,
	     mb->lda);
      if (!ScaLAPACK_RedistApiMulToApiInt(hi[h].len, vi[v].len, &block_elems) ||
          !ScaLAPACK_RedistApiAddToApiInt(sizebuff, block_elems, &sizebuff)) {
        fprintf(stderr, "xxGEMR2D:block size overflow\n");
        exit(1);
      }
    }
  }
}
static2 Int
inter_len(Int hinb, IDESC *hi, Int vinb, IDESC *vi)
{
  Int   hlen, vlen, h, v, total;
  hlen = 0;
  for (h = 0; h < hinb; h++)
    if (!ScaLAPACK_RedistApiAddToApiInt(hlen, hi[h].len, &hlen)) {
      fprintf(stderr, "xxGEMR2D:horizontal span overflow\n");
      exit(1);
    }
  vlen = 0;
  for (v = 0; v < vinb; v++)
    if (!ScaLAPACK_RedistApiAddToApiInt(vlen, vi[v].len, &vlen)) {
      fprintf(stderr, "xxGEMR2D:vertical span overflow\n");
      exit(1);
    }
  if (!ScaLAPACK_RedistApiMulToApiInt(hlen, vlen, &total)) {
    fprintf(stderr, "xxGEMR2D:intersection size overflow\n");
    exit(1);
  }
  return total;
}
void
Clacpy(Int m, Int n, float *a, Int lda, float *b, Int ldb)
{
  Int   i, j;
  lda -= m;
  ldb -= m;
  assert(lda >= 0 && ldb >= 0);
  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++)
      *b++ = *a++;
    b += ldb;
    a += lda;
  }
}
static2 void
gridreshape(Int *ctxtp)
{
  Int   ori, final;	/* original context, and new context created, with
			 * line form */
  ScaLAPACK_ApiInt line_np_api;
  size_t alloc_count, alloc_bytes;
  Int   nprow, npcol, myrow, mycol;
  Int  *usermap;
  Int   i, j;
  ori = *ctxtp;
  Cblacs_gridinfo(ori, &nprow, &npcol, &myrow, &mycol);
  if (!ScaLAPACK_RedistApiMulToSizeT(nprow, npcol, &alloc_count) ||
      !ScaLAPACK_Index64ToApiInt((ScaLAPACK_Index64)alloc_count, &line_np_api) ||
      !ScaLAPACK_SizeTMul(alloc_count, sizeof(Int), &alloc_bytes)) {
    fprintf(stderr, "xxGEMR2D:grid reshape workspace overflow\n");
    exit(1);
  }
  usermap = mr2d_malloc(alloc_bytes);
  for (i = 0; i < nprow; i++)
    for (j = 0; j < npcol; j++) {
      usermap[i + j * nprow] = Cblacs_pnum(ori, i, j);
    }
  /* Cblacs_get(0, 0, &final); */
  Cblacs_get(ori, (Int)10, &final);
  Cblacs_gridmap(&final, usermap, (Int)1, (Int)1, (Int)line_np_api);
  *ctxtp = final;
  free(usermap);
}

/* ================================================================== */
/* Core implementation — all dimensions are ScaLAPACK_Index64         */
/* ================================================================== */

#define Mlacpy_core(mo,no,ao,ldao,bo,ldbo) \
{ \
float *_a,*_b; \
ScaLAPACK_Index64 _m,_n,_lda,_ldb,_i,_j; \
    _m = (mo);_n = (no); \
    _a = (ao);_b = (bo); \
    _lda = (ldao) - _m; \
    _ldb = (ldbo) - _m; \
    assert(_lda >= 0 && _ldb >= 0); \
    for (_j=0;_j<_n;_j++) { \
      for (_i=0;_i<_m;_i++) \
        *_b++ = *_a++; \
      _b += _ldb; \
      _a += _lda; \
    } \
} (void)0

static void
setmemory_core(float **adpointer, ScaLAPACK_Index64 blocksize)
{
  size_t alloc_bytes;
  assert(blocksize >= 0);
  if (blocksize == 0) {
    *adpointer = NULL;
    return;
  }
  if (!ScaLAPACK_Index64ToSizeT(blocksize, &alloc_bytes) ||
      !ScaLAPACK_SizeTMul(alloc_bytes, sizeof(float), &alloc_bytes)) {
    fprintf(stderr, "xxGEMR2D_CORE:buffer workspace overflow\n");
    exit(1);
  }
  *adpointer = (float *) mr2d_malloc(alloc_bytes);
}

static ScaLAPACK_Index64
block2buff_core(IDESC_CORE *vi, ScaLAPACK_Index64 vinb,
		IDESC_CORE *hi, ScaLAPACK_Index64 hinb,
		float *ptra, const MDESC_CORE *ma, float *buff)
{
  ScaLAPACK_Index64 h, v, sizebuff;
  sizebuff = 0;
  for (h = 0; h < hinb; h++) {
    ScaLAPACK_Index64 row_off64;
    size_t row_off;
    float *ptr2;
    if (!ScaLAPACK_Index64Mul(hi[h].lstart, ma->lld, &row_off64) ||
	!ScaLAPACK_Index64ToSizeT(row_off64, &row_off)) {
      fprintf(stderr, "xxGEMR2D_CORE:block row offset overflow\n");
      exit(1);
    }
    ptr2 = ptra + row_off;
    for (v = 0; v < vinb; v++) {
      size_t vs;
      if (!ScaLAPACK_Index64ToSizeT(vi[v].lstart, &vs)) {
	fprintf(stderr, "xxGEMR2D_CORE:block lstart overflow\n");
	exit(1);
      }
      Mlacpy_core(vi[v].len, hi[h].len,
		  ptr2 + vs, ma->lld,
		  buff + sizebuff, vi[v].len);
      {
	ScaLAPACK_Index64 block_elems;
	if (!ScaLAPACK_Index64Mul(hi[h].len, vi[v].len, &block_elems) ||
	    !ScaLAPACK_Index64Add(sizebuff, block_elems, &sizebuff)) {
	  fprintf(stderr, "xxGEMR2D_CORE:block size overflow\n");
	  exit(1);
	}
      }
    }
  }
  return sizebuff;
}

static void
buff2block_core(IDESC_CORE *vi, ScaLAPACK_Index64 vinb,
		IDESC_CORE *hi, ScaLAPACK_Index64 hinb,
		float *buff, float *ptrb, const MDESC_CORE *mb)
{
  ScaLAPACK_Index64 h, v, sizebuff;
  sizebuff = 0;
  for (h = 0; h < hinb; h++) {
    ScaLAPACK_Index64 row_off64;
    size_t row_off;
    float *ptr2;
    if (!ScaLAPACK_Index64Mul(hi[h].lstart, mb->lld, &row_off64) ||
	!ScaLAPACK_Index64ToSizeT(row_off64, &row_off)) {
      fprintf(stderr, "xxGEMR2D_CORE:block row offset overflow\n");
      exit(1);
    }
    ptr2 = ptrb + row_off;
    for (v = 0; v < vinb; v++) {
      size_t vs;
      if (!ScaLAPACK_Index64ToSizeT(vi[v].lstart, &vs)) {
	fprintf(stderr, "xxGEMR2D_CORE:block lstart overflow\n");
	exit(1);
      }
      Mlacpy_core(vi[v].len, hi[h].len,
		  buff + sizebuff, vi[v].len,
		  ptr2 + vs, mb->lld);
      {
	ScaLAPACK_Index64 block_elems;
	if (!ScaLAPACK_Index64Mul(hi[h].len, vi[v].len, &block_elems) ||
	    !ScaLAPACK_Index64Add(sizebuff, block_elems, &sizebuff)) {
	  fprintf(stderr, "xxGEMR2D_CORE:block size overflow\n");
	  exit(1);
	}
      }
    }
  }
}

static ScaLAPACK_Index64
inter_len_core(ScaLAPACK_Index64 hinb, IDESC_CORE *hi,
	       ScaLAPACK_Index64 vinb, IDESC_CORE *vi)
{
  ScaLAPACK_Index64 hlen, vlen, h, v, total;
  hlen = 0;
  for (h = 0; h < hinb; h++)
    if (!ScaLAPACK_Index64Add(hlen, hi[h].len, &hlen)) {
      fprintf(stderr, "xxGEMR2D_CORE:horizontal span overflow\n");
      exit(1);
    }
  vlen = 0;
  for (v = 0; v < vinb; v++)
    if (!ScaLAPACK_Index64Add(vlen, vi[v].len, &vlen)) {
      fprintf(stderr, "xxGEMR2D_CORE:vertical span overflow\n");
      exit(1);
    }
  if (!ScaLAPACK_Index64Mul(hlen, vlen, &total)) {
    fprintf(stderr, "xxGEMR2D_CORE:intersection size overflow\n");
    exit(1);
  }
  return total;
}

static void
Clacpy_core(ScaLAPACK_Index64 m, ScaLAPACK_Index64 n,
	    float *a, ScaLAPACK_Index64 lda,
	    float *b, ScaLAPACK_Index64 ldb)
{
  ScaLAPACK_Index64 i, j;
  lda -= m;
  ldb -= m;
  assert(lda >= 0 && ldb >= 0);
  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++)
      *b++ = *a++;
    b += ldb;
    a += lda;
  }
}

/* Checked narrowing from Index64 to Int for BLACS message counts */
static Int
i64_to_blacs_count(ScaLAPACK_Index64 v)
{
  ScaLAPACK_ApiInt r;
  if (!ScaLAPACK_Index64ToApiInt(v, &r)) {
    fprintf(stderr, "xxGEMR2D_CORE:BLACS message count overflow "
	    "(count=%lld exceeds Int range)\n", (long long)v);
    exit(1);
  }
  return (Int)r;
}

/* ------------------------------------------------------------------ */
/* Cpsgemr2d_core — full redistribution with Index64 dimensions       */
/* ------------------------------------------------------------------ */

void
Cpsgemr2d_core(ScaLAPACK_Index64 m, ScaLAPACK_Index64 n,
	       float *ptrmyblock,
	       ScaLAPACK_Index64 ia, ScaLAPACK_Index64 ja,
	       MDESC_CORE *ma,
	       float *ptrmynewblock,
	       ScaLAPACK_Index64 ib, ScaLAPACK_Index64 jb,
	       MDESC_CORE *mb,
	       int globcontext)
{
  float *ptrsendbuff = NULL, *ptrrecvbuff = NULL;
  float *recvptr;
  MDESC_CORE newa, newb;
  size_t alloc_count, alloc_bytes;
  Int   mypnum_b, nprow_b, npcol_b, dummy_b;
  Int   p0_b, q0_b, myprow0_b, mypcol0_b;
  Int   p1_b, q1_b, myprow1_b, mypcol1_b;
  int   mypnum, myprow0, mypcol0, myprow1, mypcol1, nprocs;
  int   p0, q0, p1, q1;
  Int   gcontext;
  ScaLAPACK_Index64 *param64;
  ScaLAPACK_Index64 param_span64;
  Int  *proc0, *proc1;
  int   proc0_span, proc1_span;
  IDESC_CORE *h_inter, *v_inter;
  ScaLAPACK_Index64 hinter_nb, vinter_nb;
  ScaLAPACK_Index64 sendsize = 0, recvsize;
  int   i;

  if (m == 0 || n == 0)
    return;
  /* Convert from 1-based Fortran to 0-based */
  ia -= 1;
  ja -= 1;
  ib -= 1;
  jb -= 1;

  /* Get global grid info */
  Cblacs_gridinfo((Int)globcontext, &nprow_b, &npcol_b, &dummy_b, &mypnum_b);
  gcontext = (Int)globcontext;
  nprocs = (int)nprow_b * (int)npcol_b;
  mypnum = (int)mypnum_b;

  if ((int)nprow_b != 1) {
    gridreshape(&gcontext);
    Cblacs_gridinfo(gcontext, &dummy_b, &dummy_b, &dummy_b, &mypnum_b);
    mypnum = (int)mypnum_b;
  }

  /* Get source grid info */
  Cblacs_gridinfo((Int)ma->ctxt, &p0_b, &q0_b, &myprow0_b, &mypcol0_b);
  p0 = (int)p0_b; q0 = (int)q0_b;
  myprow0 = (int)myprow0_b; mypcol0 = (int)mypcol0_b;
  if (myprow0 >= p0 || mypcol0 >= q0)
    myprow0 = mypcol0 = -1;
  assert((myprow0 < p0 && mypcol0 < q0) || (myprow0 == -1 && mypcol0 == -1));

  /* Get destination grid info */
  Cblacs_gridinfo((Int)mb->ctxt, &p1_b, &q1_b, &myprow1_b, &mypcol1_b);
  p1 = (int)p1_b; q1 = (int)q1_b;
  myprow1 = (int)myprow1_b; mypcol1 = (int)mypcol1_b;
  if (myprow1 >= p1 || mypcol1 >= q1)
    myprow1 = mypcol1 = -1;
  assert((myprow1 < p1 && mypcol1 < q1) || (myprow1 == -1 && mypcol1 == -1));

  /* ----- Allocate and fill param64 sync array ----- */
  param_span64 = (ScaLAPACK_Index64)NBPARAM_CORE + 2 * (ScaLAPACK_Index64)nprocs;
  if (!ScaLAPACK_Index64ToSizeT(param_span64, &alloc_count) ||
      !ScaLAPACK_SizeTMul(alloc_count, sizeof(ScaLAPACK_Index64), &alloc_bytes)) {
    fprintf(stderr, "xxGEMR2D_CORE:parameter workspace overflow\n");
    exit(1);
  }
  param64 = (ScaLAPACK_Index64 *) mr2d_malloc(alloc_bytes);
  for (i = 0; i < (int)param_span64; i++)
    param64[i] = MAGIC_MAX_I8;

  if (myprow0 >= 0) {
    param64[myprow0 * q0 + mypcol0 + NBPARAM_CORE] = (ScaLAPACK_Index64)mypnum;
    param64[0]  = (ScaLAPACK_Index64)p0;
    param64[1]  = (ScaLAPACK_Index64)q0;
    param64[4]  = ma->m;
    param64[5]  = ma->n;
    param64[6]  = ma->mb;
    param64[7]  = ma->nb;
    param64[8]  = (ScaLAPACK_Index64)ma->rsrc;
    param64[9]  = (ScaLAPACK_Index64)ma->csrc;
    param64[10] = ia;
    param64[11] = ja;
  }
  if (myprow1 >= 0) {
    param64[myprow1 * q1 + mypcol1 + NBPARAM_CORE + nprocs] = (ScaLAPACK_Index64)mypnum;
    param64[2]  = (ScaLAPACK_Index64)p1;
    param64[3]  = (ScaLAPACK_Index64)q1;
    param64[12] = mb->m;
    param64[13] = mb->n;
    param64[14] = mb->mb;
    param64[15] = mb->nb;
    param64[16] = (ScaLAPACK_Index64)mb->rsrc;
    param64[17] = (ScaLAPACK_Index64)mb->csrc;
    param64[18] = ib;
    param64[19] = jb;
  }

  redist_sync_params_i8((int)gcontext, param64, (int)param_span64);

  /* ----- Extract synced parameters ----- */
  newa = *ma;
  newb = *mb;
  ma = &newa;
  mb = &newb;

  if (myprow0 == -1) {
    p0 = (int)param64[0];  q0 = (int)param64[1];
    ma->m    = param64[4];  ma->n    = param64[5];
    ma->mb   = param64[6];  ma->nb   = param64[7];
    ma->rsrc = (int)param64[8];  ma->csrc = (int)param64[9];
    ia = param64[10]; ja = param64[11];
  }
  if (myprow1 == -1) {
    p1 = (int)param64[2];  q1 = (int)param64[3];
    mb->m    = param64[12]; mb->n    = param64[13];
    mb->mb   = param64[14]; mb->nb   = param64[15];
    mb->rsrc = (int)param64[16]; mb->csrc = (int)param64[17];
    ib = param64[18]; jb = param64[19];
  }

  /* Extract proc0/proc1 as Int arrays (process numbers are small) */
  proc0_span = p0 * q0;
  proc1_span = p1 * q1;
  if (!ScaLAPACK_SizeTMul((size_t)(proc0_span + proc1_span), sizeof(Int),
			  &alloc_bytes)) {
    fprintf(stderr, "xxGEMR2D_CORE:proc workspace overflow\n");
    exit(1);
  }
  proc0 = (Int *) mr2d_malloc(alloc_bytes);
  proc1 = proc0 + proc0_span;
  for (i = 0; i < proc0_span; i++)
    proc0[i] = (Int)param64[NBPARAM_CORE + i];
  for (i = 0; i < proc1_span; i++)
    proc1[i] = (Int)param64[NBPARAM_CORE + nprocs + i];

  /* Verify all params were synced */
  for (i = 0; i < NBPARAM_CORE; i++) {
    if (param64[i] == MAGIC_MAX_I8) {
      fprintf(stderr, "xxGEMR2D_CORE:something wrong in the parameters\n");
      exit(1);
    }
  }
#ifndef NDEBUG
  for (i = 0; i < proc0_span; i++)
    assert(proc0[i] >= 0 && proc0[i] < nprocs);
  for (i = 0; i < proc1_span; i++)
    assert(proc1[i] >= 0 && proc1[i] < nprocs);
#endif

  /* Validate parameters */
  paramcheck_core(ma, ia, ja, m, n, p0, q0, (int)gcontext);
  paramcheck_core(mb, ib, jb, m, n, p1, q1, (int)gcontext);

  /* Change origin so that ia < mb, ja < nb, etc. */
  {
    ScaLAPACK_Index64 decal;
    size_t shift_bytes;
    ScaLAPACK_Index64 shift64;

    ia = changeorigin_core(myprow0, ma->rsrc, p0,
			   ma->mb, ia, &decal, &ma->rsrc);
    if (!ScaLAPACK_Index64ToSizeT(decal, &shift_bytes) ||
	!ScaLAPACK_SizeTMul(shift_bytes, sizeof(*ptrmyblock), &shift_bytes)) {
      fprintf(stderr, "xxGEMR2D_CORE:local block pointer overflow\n");
      exit(1);
    }
    ptrmyblock = (float *) ((char *) ptrmyblock + shift_bytes);

    ja = changeorigin_core(mypcol0, ma->csrc, q0,
			   ma->nb, ja, &decal, &ma->csrc);
    if (!ScaLAPACK_Index64Mul(decal, ma->lld, &shift64) ||
	!ScaLAPACK_Index64ToSizeT(shift64, &shift_bytes) ||
	!ScaLAPACK_SizeTMul(shift_bytes, sizeof(*ptrmyblock), &shift_bytes)) {
      fprintf(stderr, "xxGEMR2D_CORE:local block pointer overflow\n");
      exit(1);
    }
    ptrmyblock = (float *) ((char *) ptrmyblock + shift_bytes);

    ma->m = ia + m;
    ma->n = ja + n;

    ib = changeorigin_core(myprow1, mb->rsrc, p1,
			   mb->mb, ib, &decal, &mb->rsrc);
    if (!ScaLAPACK_Index64ToSizeT(decal, &shift_bytes) ||
	!ScaLAPACK_SizeTMul(shift_bytes, sizeof(*ptrmynewblock), &shift_bytes)) {
      fprintf(stderr, "xxGEMR2D_CORE:local destination pointer overflow\n");
      exit(1);
    }
    ptrmynewblock = (float *) ((char *) ptrmynewblock + shift_bytes);

    jb = changeorigin_core(mypcol1, mb->csrc, q1,
			   mb->nb, jb, &decal, &mb->csrc);
    if (!ScaLAPACK_Index64Mul(decal, mb->lld, &shift64) ||
	!ScaLAPACK_Index64ToSizeT(shift64, &shift_bytes) ||
	!ScaLAPACK_SizeTMul(shift_bytes, sizeof(*ptrmynewblock), &shift_bytes)) {
      fprintf(stderr, "xxGEMR2D_CORE:local destination pointer overflow\n");
      exit(1);
    }
    ptrmynewblock = (float *) ((char *) ptrmynewblock + shift_bytes);

    mb->m = ib + m;
    mb->n = jb + n;

    if (p0 == 1) ma->mb = ma->m;
    if (q0 == 1) ma->nb = ma->n;
    if (p1 == 1) mb->mb = mb->m;
    if (q1 == 1) mb->nb = mb->n;

#ifndef NDEBUG
    paramcheck_core(ma, ia, ja, m, n, p0, q0, (int)gcontext);
    paramcheck_core(mb, ib, jb, m, n, p1, q1, (int)gcontext);
#endif
  }

  /* Allocate send/recv buffers */
  if (myprow0 >= 0 && mypcol0 >= 0)
    setmemory_core(&ptrsendbuff, memoryblocksize_core(ma));
  if (myprow1 >= 0 && mypcol1 >= 0)
    setmemory_core(&ptrrecvbuff, memoryblocksize_core(mb));

  /* Allocate IDESC_CORE interval arrays (worst-case size) */
  {
    ScaLAPACK_Index64 tw0, tw1, hint_count, vint_count;
    tw0 = (ScaLAPACK_Index64)q0 * ma->nb;
    if (tw0 > 0)
      hint_count = ((ma->n - 1) / tw0 + 1) * ma->nb;
    else
      hint_count = 0;
    if (!ScaLAPACK_Index64ToSizeT(hint_count, &alloc_count) ||
	!ScaLAPACK_SizeTMul(alloc_count, sizeof(IDESC_CORE), &alloc_bytes)) {
      fprintf(stderr, "xxGEMR2D_CORE:horizontal interval workspace overflow\n");
      exit(1);
    }
    h_inter = (IDESC_CORE *) mr2d_malloc(alloc_bytes);

    tw1 = (ScaLAPACK_Index64)p0 * ma->mb;
    if (tw1 > 0)
      vint_count = ((ma->m - 1) / tw1 + 1) * ma->mb;
    else
      vint_count = 0;
    if (!ScaLAPACK_Index64ToSizeT(vint_count, &alloc_count) ||
	!ScaLAPACK_SizeTMul(alloc_count, sizeof(IDESC_CORE), &alloc_bytes)) {
      fprintf(stderr, "xxGEMR2D_CORE:vertical interval workspace overflow\n");
      exit(1);
    }
    v_inter = (IDESC_CORE *) mr2d_malloc(alloc_bytes);
  }

  /* Communication loop */
  recvptr = ptrrecvbuff;
  {
    Int   tot, myrang, step, sens;
    Int  *sender, *recver;
    Int   mesending, merecving;
    Int   ii, jj;
    tot = max(proc0_span, proc1_span);
    init_chenille((Int)mypnum, (Int)nprocs,
		  (Int)proc0_span, proc0, (Int)proc1_span, proc1,
		  &sender, &recver, &myrang);
    if (myrang == -1)
      goto after_comm;
    mesending = myprow0 >= 0;
    assert(sender[myrang] >= 0 || !mesending);
    assert(!mesending || proc0[sender[myrang]] == (Int)mypnum);
    merecving = myprow1 >= 0;
    assert(recver[myrang] >= 0 || !merecving);
    assert(!merecving || proc1[recver[myrang]] == (Int)mypnum);
    step = tot - 1 - myrang;
    do {
      for (sens = 0; sens < 2; sens++) {
	if (mesending && recver[step] >= 0 && (sens == 0)) {
	  ii = recver[step] / q1;
	  jj = recver[step] % q1;
	  vinter_nb = scan_intervals_core('r', ia, ib, m, ma, mb,
					  p0, p1, myprow0, (int)ii, v_inter);
	  hinter_nb = scan_intervals_core('c', ja, jb, n, ma, mb,
					  q0, q1, mypcol0, (int)jj, h_inter);
	  sendsize = block2buff_core(v_inter, vinter_nb, h_inter, hinter_nb,
				     ptrmyblock, ma, ptrsendbuff);
	}
	if (mesending && recver[step] >= 0 && (sens == myrang > step)) {
	  ii = recver[step] / q1;
	  jj = recver[step] % q1;
	  if (sendsize > 0 && (step != myrang || !merecving)) {
	    Int ss = i64_to_blacs_count(sendsize);
	    Csgesd2d(gcontext, ss, (Int)1, ptrsendbuff, ss,
		     (Int)0, proc1[ii * q1 + jj]);
	  }
	}
	if (merecving && sender[step] >= 0 && (sens == myrang <= step)) {
	  ii = sender[step] / q0;
	  jj = sender[step] % q0;
	  vinter_nb = scan_intervals_core('r', ib, ia, m, mb, ma,
					  p1, p0, myprow1, (int)ii, v_inter);
	  hinter_nb = scan_intervals_core('c', jb, ja, n, mb, ma,
					  q1, q0, mypcol1, (int)jj, h_inter);
	  recvsize = inter_len_core(hinter_nb, h_inter, vinter_nb, v_inter);
	  if (recvsize > 0) {
	    if (step == myrang && mesending) {
	      Clacpy_core(recvsize, 1,
			  ptrsendbuff, recvsize,
			  ptrrecvbuff, recvsize);
	    } else {
	      Int rs = i64_to_blacs_count(recvsize);
	      Csgerv2d(gcontext, rs, (Int)1, ptrrecvbuff, rs,
		       (Int)0, proc0[ii * q0 + jj]);
	    }
	  }
	}
	if (merecving && sender[step] >= 0 && sens == 1) {
	  buff2block_core(v_inter, vinter_nb, h_inter, hinter_nb,
			  recvptr, ptrmynewblock, mb);
	}
      }
      step -= 1;
      if (step < 0)
	step = tot - 1;
    } while (step != tot - 1 - myrang);
after_comm:
    free(sender);
  }

  /* Cleanup */
  if (myprow1 >= 0 && mypcol1 >= 0)
    free(ptrrecvbuff);
  if (myprow0 >= 0 && mypcol0 >= 0)
    free(ptrsendbuff);
  if ((int)nprow_b != 1)
    Cblacs_gridexit(gcontext);
  free(v_inter);
  free(h_inter);
  free(proc0);
  free(param64);
}

/* ------------------------------------------------------------------ */
/* I8 Fortran entry point                                             */
/* ------------------------------------------------------------------ */

void
fortran_mr2dnew_i8(int64_t *m, int64_t *n,
		   float *A, int64_t *ia, int64_t *ja, int64_t desc_A[9],
		   float *B, int64_t *ib, int64_t *jb, int64_t desc_B[9],
		   int64_t *gcontext)
{
  MDESC_CORE core_a, core_b;
  int gc;

  if (!unpack_desc_i8(desc_A, &core_a)) {
    fprintf(stderr, "PSGEMR2D_I8:bad source descriptor (narrowing overflow)\n");
    exit(1);
  }
  if (!unpack_desc_i8(desc_B, &core_b)) {
    fprintf(stderr, "PSGEMR2D_I8:bad destination descriptor (narrowing overflow)\n");
    exit(1);
  }

  /* gcontext is a BLACS handle — must fit in int */
  if (*gcontext < INT_MIN || *gcontext > INT_MAX) {
    fprintf(stderr, "PSGEMR2D_I8:bad global context (narrowing overflow)\n");
    exit(1);
  }
  gc = (int)*gcontext;

  Cpsgemr2d_core(*m, *n, A, *ia, *ja, &core_a,
		 B, *ib, *jb, &core_b, gc);
}
