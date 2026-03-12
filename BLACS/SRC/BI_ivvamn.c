#include "Bdef.h"
void BI_ivvamn(MpiInt N, char *vec1, char *vec2)
{
   Int *v1=(Int*)vec1, *v2=(Int*)vec2;
   Int diff;
   BI_DistType *dist1, *dist2;
   size_t count, dist_offset;
   MpiInt k;

   count = (size_t) N;
   ScaLAPACK_SizeTMul(count, sizeof(Int), &dist_offset);
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
