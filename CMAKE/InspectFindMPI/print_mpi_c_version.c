#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{
  int ierr;
  int version_len = 0;
  char library_version[MPI_MAX_LIBRARY_VERSION_STRING];

  printf("C header MPI_VERSION=%d\n", MPI_VERSION);
  printf("C header MPI_SUBVERSION=%d\n", MPI_SUBVERSION);

  ierr = MPI_Init(&argc, &argv);
  printf("MPI_Init ierr=%d\n", ierr);
  if (ierr != MPI_SUCCESS) {
    return ierr;
  }

  ierr = MPI_Get_library_version(library_version, &version_len);
  printf("MPI_Get_library_version ierr=%d\n", ierr);
  if (ierr == MPI_SUCCESS) {
    printf("MPI library version=%.*s\n", version_len, library_version);
  }

  ierr = MPI_Finalize();
  printf("MPI_Finalize ierr=%d\n", ierr);
  return ierr;
}

