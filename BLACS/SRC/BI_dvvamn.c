#include "Bdef.h"
void BI_dvvamn(Int N, char *vec1, char *vec2)
{
   double *v1=(double*)vec1, *v2=(double*)vec2;
   double diff;
   BI_DistType *dist1, *dist2;
   size_t count, dist_offset;
   Int k;

   count = (size_t) N;
   ScaLAPACK_SizeTMul(count, sizeof(double), &dist_offset);
   ScaLAPACK_SizeTAlignUp(dist_offset, sizeof(BI_DistType), &dist_offset);
   dist1 = (BI_DistType *) &vec1[dist_offset];
   dist2 = (BI_DistType *) &vec2[dist_offset];

   for (k=0; k < N; k++)
   {
      diff = Rabs(v1[k]) - Rabs(v2[k]);
      if (diff > 0)
      {
         v1[k] = v2[k];
         dist1[k] = dist2[k];
      }
      else if (diff == 0)
      {
         if (dist1[k] > dist2[k])
         {
            v1[k] = v2[k];
            dist1[k] = dist2[k];
         }
      }
   }
}
