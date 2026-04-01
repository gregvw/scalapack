include_guard(GLOBAL)

macro(scalapack_configure_mpi)
  find_package(MPI)
  if (MPI_FOUND)
     message(STATUS "Found MPI_LIBRARY : ${MPI_FOUND} ")

     find_program(MPI_C_COMPILER
        NAMES mpicc
        HINTS "${MPI_BASE_DIR}"
        PATH_SUFFIXES bin
        DOC "MPI C compiler.")
     mark_as_advanced(MPI_C_COMPILER)
     if ("${MPI_C_COMPILER}" STREQUAL "MPI_C_COMPILER-NOTFOUND")
        message(ERROR "--> MPI C Compiler NOT FOUND (please set MPI_BASE_DIR accordingly")
     else()
        message(STATUS "--> MPI C Compiler : ${MPI_C_COMPILER}")
        message(STATUS "--> C Compiler : ${CMAKE_C_COMPILER}")
     endif()

     find_program(MPI_Fortran_COMPILER
        NAMES mpif77
        HINTS "${MPI_BASE_DIR}"
        PATH_SUFFIXES bin
        DOC "MPI Fortran compiler.")
     mark_as_advanced(MPI_Fortran_COMPILER)
     if ("${MPI_Fortran_COMPILER}" STREQUAL "MPI_Fortran_COMPILER-NOTFOUND")
        message(ERROR "--> MPI Fortran Compiler NOT FOUND (please set MPI_BASE_DIR accordingly")
     else()
        message(STATUS "--> MPI Fortran Compiler : ${MPI_Fortran_COMPILER}")
        set(Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")
        message(STATUS "--> Fortran Compiler : ${CMAKE_Fortran_COMPILER}")
     endif()
  else()
     message(STATUS "Found MPI_LIBRARY : ${MPI_FOUND} ")
     set(MPI_BASE_DIR ${MPI_BASE_DIR} CACHE PATH "MPI Path")
     unset(MPIEXEC CACHE)
     unset(MPIEXEC_POSTFLAGS CACHE)
     unset(MPIEXEC_PREFLAGS CACHE)
     unset(MPIEXEC_MAX_NUMPROCS CACHE)
     unset(MPIEXEC_NUMPROC_FLAG CACHE)
     unset(MPI_COMPILE_FLAGS CACHE)
     unset(MPI_LINK_FLAGS CACHE)
     unset(MPI_INCLUDE_PATH CACHE)
     message(FATAL_ERROR "--> MPI Library NOT FOUND -- please set MPI_BASE_DIR accordingly --")
  endif()

  message(STATUS "--> MPI C Compiler : ${MPI_C_COMPILER}")
  message(STATUS "--> C Compiler : ${CMAKE_C_COMPILER}")
  message(STATUS "--> MPI Fortran Compiler : ${MPI_Fortran_COMPILER}")
  message(STATUS "--> Fortran Compiler : ${CMAKE_Fortran_COMPILER}")
endmacro()
