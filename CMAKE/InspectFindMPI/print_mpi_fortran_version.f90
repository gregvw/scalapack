program scalapack_probe_mpi_fortran
  use mpi
  implicit none

  integer :: ierr
  integer :: version_len
  character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: library_version

  print *, 'Fortran module MPI_VERSION=', MPI_VERSION
  print *, 'Fortran module MPI_SUBVERSION=', MPI_SUBVERSION

  call MPI_Init(ierr)
  print *, 'MPI_Init ierr=', ierr
  if (ierr .ne. MPI_SUCCESS) stop 1

  call MPI_Get_library_version(library_version, version_len, ierr)
  print *, 'MPI_Get_library_version ierr=', ierr
  if (ierr .eq. MPI_SUCCESS) then
     print *, 'MPI library version=', trim(library_version(:version_len))
  end if

  call MPI_Finalize(ierr)
  print *, 'MPI_Finalize ierr=', ierr
end program scalapack_probe_mpi_fortran
