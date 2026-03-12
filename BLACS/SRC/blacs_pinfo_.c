#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_pinfo(Int *mypnum, Int *nprocs)
#else
F_VOID_FUNC blacs_pinfo_(Int *mypnum, Int *nprocs)
#endif
{
   int ierr;
   size_t alloc_bytes;
   extern Int BI_Iam, BI_Np;
   int flag, Iam = BI_Iam, Np = BI_Np;
   int argc=0;
   char **argv=NULL;
   if (BI_COMM_WORLD == NULL)
   {
      MPI_Initialized(&flag);

      if (!flag) 
         ierr = MPI_Init(&argc,&argv);  // call Init and ignore argc and argv

      if (!ScaLAPACK_SizeTMul((size_t)1, sizeof(Int), &alloc_bytes))
         BI_BlacsErr(-1, __LINE__, __FILE__, "Cannot size BI_COMM_WORLD");
      BI_COMM_WORLD = (Int *) malloc(alloc_bytes);
      if (!BI_COMM_WORLD)
         BI_BlacsErr(-1, __LINE__, __FILE__, "Cannot allocate BI_COMM_WORLD");
      *BI_COMM_WORLD = MPI_Comm_c2f(MPI_COMM_WORLD);
   }
   MPI_Comm_size(MPI_COMM_WORLD, &Np);
   MPI_Comm_rank(MPI_COMM_WORLD, &Iam);
   *mypnum = BI_Iam = Iam;
   *nprocs = BI_Np  = Np;
}
