#include "redist.h"
#include <stddef.h>
/* $Id: psgemr2.c,v 1.1.1.1 2000/02/15 18:04:09 susan Exp $
 * 
 * some functions used by the psgemr2d routine see file psgemr.c for more
 * documentation.
 * 
 * Created March 1993 by B. Tourancheau (See sccs for modifications). */
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
extern void freememory( float* ptr );
extern Int scan_intervals( char type, Int ja, Int jb, Int n, MDESC *ma, MDESC *mb, Int q0, Int q1, Int col0, Int col1, IDESC *result );
extern void Cpsgemr2do();
extern void Cpsgemr2d();
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
#include <string.h>
#include <assert.h>
#include <ctype.h>
/* Created March 1993 by B. Tourancheau (See sccs for modifications). */
/************************************************************************/
/* Set the memory space with the malloc function */
void
setmemory(float **adpointer, Int blocksize)
{
  size_t alloc_count, alloc_bytes;
  assert(blocksize >= 0);
  if (blocksize == 0) {
    *adpointer = NULL;
    return;
  }
  if (!ScaLAPACK_Index64ToSizeT((ScaLAPACK_Index64)blocksize, &alloc_count) ||
      !ScaLAPACK_SizeTMul(alloc_count, sizeof(float), &alloc_bytes)) {
    fprintf(stderr, "xxGEMR2D:buffer workspace overflow\n");
    exit(1);
  }
  *adpointer = (float *) mr2d_malloc(alloc_bytes);
}
/******************************************************************/
/* Free the memory space after the malloc */
void
freememory(float *ptrtobefreed)
{
  if (ptrtobefreed == NULL)
    return;
  free((char *) ptrtobefreed);
}
/* extern functions for intersect() extern scopy_(); */
/* scan_intervals: scans two distributions in one dimension, and compute the
 * intersections on the local processor. result must be long enough to
 * contains the result that are stocked in IDESC structure, the function
 * returns the number of intersections found */
Int 
scan_intervals(type, ja, jb, n, ma, mb, q0, q1, col0, col1,
	       result)
  char  type;
  Int   ja, jb, n, q0, q1, col0, col1;
  MDESC *ma, *mb;
  IDESC *result;
{
  Int   offset, j0, j1, templatewidth0, templatewidth1, nbcol0, nbcol1;
  Int   l, shifted0, shifted1, end0, end1;	/* local indice on the beginning of the interval */
  assert(type == 'c' || type == 'r');
  nbcol0 = (type == 'c' ? ma->nbcol : ma->nbrow);
  nbcol1 = (type == 'c' ? mb->nbcol : mb->nbrow);
  if (!ScaLAPACK_RedistApiMulToApiInt(q0, nbcol0, &templatewidth0) ||
      !ScaLAPACK_RedistApiMulToApiInt(q1, nbcol1, &templatewidth1)) {
    fprintf(stderr, "xxGEMR2D:template width overflow\n");
    exit(1);
  }
  {
    Int   sp0 = (type == 'c' ? ma->spcol : ma->sprow);
    Int   sp1 = (type == 'c' ? mb->spcol : mb->sprow);
    if (!ScaLAPACK_RedistApiMulToApiInt(SHIFT(col0, sp0, q0), nbcol0, &shifted0) ||
        !ScaLAPACK_RedistApiMulToApiInt(SHIFT(col1, sp1, q1), nbcol1, &shifted1) ||
        !ScaLAPACK_RedistApiSubToApiInt(shifted0, ja, &j0) ||
        !ScaLAPACK_RedistApiSubToApiInt(shifted1, jb, &j1)) {
      fprintf(stderr, "xxGEMR2D:template offset overflow\n");
      exit(1);
    }
  }
  offset = 0;
  l = 0;
  /* a small check to verify that the submatrix begin inside the first block
   * of the original matrix, this done by a sort of coordinate change at the
   * beginning of the Cpsgemr2d */
  assert(j0 + nbcol0 > 0);
  assert(j1 + nbcol1 > 0);
  while ((j0 < n) && (j1 < n)) {
    Int   start, end;
    if (!ScaLAPACK_RedistApiAddToApiInt(j0, nbcol0, &end0) ||
        !ScaLAPACK_RedistApiAddToApiInt(j1, nbcol1, &end1)) {
      fprintf(stderr, "xxGEMR2D:interval bound overflow\n");
      exit(1);
    }
    if (end0 <= j1) {
      if (!ScaLAPACK_RedistApiAddToApiInt(j0, templatewidth0, &j0) ||
          !ScaLAPACK_RedistApiAddToApiInt(l, nbcol0, &l)) {
        fprintf(stderr, "xxGEMR2D:interval advance overflow\n");
        exit(1);
      }
      continue;
    }
    if (end1 <= j0) {
      if (!ScaLAPACK_RedistApiAddToApiInt(j1, templatewidth1, &j1)) {
        fprintf(stderr, "xxGEMR2D:interval advance overflow\n");
        exit(1);
      }
      continue;
    }
    /* compute the raw intersection */
    start = max(j0, j1);
    start = max(start, 0);
    /* the start is correct now, update the corresponding fields */
    if (!ScaLAPACK_RedistApiAddToApiInt(l, start, &result[offset].lstart) ||
        !ScaLAPACK_RedistApiSubToApiInt(result[offset].lstart, j0,
                                        &result[offset].lstart)) {
      fprintf(stderr, "xxGEMR2D:local start overflow\n");
      exit(1);
    }
    end = min(end0, end1);
    if (end0 == end) {
      if (!ScaLAPACK_RedistApiAddToApiInt(j0, templatewidth0, &j0) ||
          !ScaLAPACK_RedistApiAddToApiInt(l, nbcol0, &l)) {
        fprintf(stderr, "xxGEMR2D:interval advance overflow\n");
        exit(1);
      }
    }
    if (end1 == end)
      if (!ScaLAPACK_RedistApiAddToApiInt(j1, templatewidth1, &j1)) {
        fprintf(stderr, "xxGEMR2D:interval advance overflow\n");
        exit(1);
      }
    /* throw the limit if they go out of the matrix */
    end = min(end, n);
    assert(end > start);
    /* it is a bit tricky to see why the length is always positive after all
     * this min and max, first we have the property that every interval
     * considered is at least partly into the submatrix, second we arrive
     * here only if the raw intersection is non-void, if we remove a limit
     * that means the corresponding frontier is in both intervals which
     * proove the final interval is non-void, clear ?? */
    if (!ScaLAPACK_RedistApiSubToApiInt(end, start, &result[offset].len) ||
        !ScaLAPACK_RedistApiAddToApiInt(offset, (Int)1, &offset)) {
      fprintf(stderr, "xxGEMR2D:interval descriptor overflow\n");
      exit(1);
    }
  }	/* while */
  return offset;
}
