#include "Bdef.h"
void BI_zvvamx(Int N, char *vec1, char *vec2)
{
   DCOMPLEX *v1=(DCOMPLEX*)vec1, *v2=(DCOMPLEX*)vec2;
   double diff;
   BI_DistType *dist1, *dist2;
   size_t count, dist_offset;
   Int k;

   count = (size_t) N;
   ScaLAPACK_SizeTMul(count, sizeof(DCOMPLEX), &dist_offset);
   ScaLAPACK_SizeTAlignUp(dist_offset, sizeof(BI_DistType), &dist_offset);
   dist1 = (BI_DistType *) &vec1[dist_offset];
   dist2 = (BI_DistType *) &vec2[dist_offset];

   for (k=0; k < N; k++)
   {
      diff = Cabs(v1[k]) - Cabs(v2[k]);
      if (diff < 0)
      {
         v1[k].r = v2[k].r;
         v1[k].i = v2[k].i;
         dist1[k] = dist2[k];
      }
      else if (diff == 0)
      {
         if (dist1[k] > dist2[k])
         {
            v1[k].r = v2[k].r;
            v1[k].i = v2[k].i;
            dist1[k] = dist2[k];
         }
      }
   }
}
