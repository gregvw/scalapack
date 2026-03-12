//
//  lamov.h
//
//  Written by Lee Killough 04/19/2012
//  

#include "pblas.h"
#include <ctype.h>

extern void xerbla_(const char *, const F_INTG_FCT *, size_t);

#define LAMOV_ELEM(base_, row_, col_, ld_) \
   ((base_)[ScaLAPACK_Index64MatrixOffset((ScaLAPACK_Index64)(row_), \
                                          (ScaLAPACK_Index64)(col_), \
                                          (ScaLAPACK_Index64)(ld_), 1)])

void LACPY(const char *UPLO,
           const F_INTG_FCT *M,
           const F_INTG_FCT *N,
           const TYPE *A,
           const F_INTG_FCT *LDA,
           TYPE *B,
           const F_INTG_FCT *LDB);

void LAMOV(const char *UPLO,
           const F_INTG_FCT *M,
           const F_INTG_FCT *N,
           const TYPE *A,
           const F_INTG_FCT *LDA,
           TYPE *B,
           const F_INTG_FCT *LDB)
{
   const F_INTG_FCT m = *M;
   const F_INTG_FCT n = *N;
   const F_INTG_FCT lda = *LDA;
   const F_INTG_FCT ldb = *LDB;
   ScaLAPACK_Index64 a_last_offset, b_last_offset;
   const TYPE *a_last;
   TYPE *b_last;

   if (m <= 0 || n <= 0)
     return;

   if (!ScaLAPACK_Index64MatrixOffsetChecked((ScaLAPACK_Index64)(m - 1),
                                             (ScaLAPACK_Index64)(n - 1),
                                             (ScaLAPACK_Index64)lda, 1,
                                             &a_last_offset) ||
       !ScaLAPACK_Index64MatrixOffsetChecked((ScaLAPACK_Index64)(m - 1),
                                             (ScaLAPACK_Index64)(n - 1),
                                             (ScaLAPACK_Index64)ldb, 1,
                                             &b_last_offset))
     {
       F_INTG_FCT info = -1;
       const char func[] = FUNC;
       xerbla_(func, &info, sizeof func);
       return;
     }

   a_last = A + a_last_offset;
   b_last = B + b_last_offset;

   if (b_last < A || a_last < B)
     {
       LACPY(UPLO, M, N, A, LDA, B, LDB);
     }
   else if (lda != ldb)
     {
       TYPE *tmp;
       ScaLAPACK_Index64 elem_count;
       size_t alloc_elems, alloc_bytes;
       if (!ScaLAPACK_Index64Mul((ScaLAPACK_Index64) m,
                                 (ScaLAPACK_Index64) n,
                                 &elem_count) ||
           !ScaLAPACK_Index64ToSizeT(elem_count, &alloc_elems) ||
           !ScaLAPACK_SizeTMul(alloc_elems, sizeof(*A), &alloc_bytes))
         tmp = NULL;
       else
         tmp = malloc(alloc_bytes);
       if (!tmp)
         {
           F_INTG_FCT info = -1;
           const char func[] = FUNC;
           xerbla_(func, &info, sizeof func);
         }
       else
         {
           LACPY(UPLO, M, N,   A, LDA, tmp,  &m);
           LACPY(UPLO, M, N, tmp,  &m,   B, LDB);
           free(tmp);
         }
     }
   else
     {
       F_INTG_FCT i, j;
       switch (toupper(*UPLO))
         {
         case 'U':
           if (A > B)
             {
               for (j=0; j<n; j++)
                 for (i=0; i<j && i<m; i++)
                   LAMOV_ELEM(B, i, j, ldb) = LAMOV_ELEM(A, i, j, lda);
             }
           else
             {
               for (j=n; --j>=0;)
                 for (i=j<m ? j : m; --i>=0;)
                   LAMOV_ELEM(B, i, j, ldb) = LAMOV_ELEM(A, i, j, lda);
             }
           break;
         
         case 'L':
           if (A > B)
             {
               for (j=0; j<n; j++)
                 for (i=j; i<m; i++)
                   LAMOV_ELEM(B, i, j, ldb) = LAMOV_ELEM(A, i, j, lda);
             }
           else
             {
               for (j=m<n ? m : n; --j>=0;)
                 for (i=m; --i>=j;)
                   LAMOV_ELEM(B, i, j, ldb) = LAMOV_ELEM(A, i, j, lda);
             }
           break;
         
         default:
           if (A > B)
             {
               for (j=0; j<n; j++)
                 for (i=0; i<m; i++)
                   LAMOV_ELEM(B, i, j, ldb) = LAMOV_ELEM(A, i, j, lda);
             }
           else
             {
               for (j=n; --j>=0;)
                 for (i=m; --i>=0;)
                   LAMOV_ELEM(B, i, j, ldb) = LAMOV_ELEM(A, i, j, lda);
             }
           break;
         }
     }
}

#undef LAMOV_ELEM
