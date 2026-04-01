cmake_minimum_required(VERSION 3.26)

project(ScaLAPACKInspectFindMPI LANGUAGES C Fortran)

include(CheckCSourceCompiles)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/..")
include(CheckMpiLargeCount)

option(MPI_PROBE_ENABLE_FIND_DEBUG_MODE
       "Enable CMake's CMAKE_FIND_DEBUG_MODE while running FindMPI"
       OFF)
if(MPI_PROBE_ENABLE_FIND_DEBUG_MODE)
  set(CMAKE_FIND_DEBUG_MODE TRUE)
endif()

# Ask FindMPI to interrogate the implementation as thoroughly as it can.
set(MPI_DETERMINE_LIBRARY_VERSION TRUE CACHE BOOL "" FORCE)
set(MPI_DETERMINE_Fortran_CAPABILITIES TRUE CACHE BOOL "" FORCE)

function(_probe_print_var name)
  if(DEFINED ${name})
    message(STATUS "${name}=[${${name}}]")
  else()
    message(STATUS "${name}=<UNDEFINED>")
  endif()
endfunction()

message(STATUS "===== Probe Environment =====")
_probe_print_var(CMAKE_VERSION)
_probe_print_var(CMAKE_ROOT)
_probe_print_var(CMAKE_GENERATOR)
_probe_print_var(CMAKE_MODULE_PATH)
_probe_print_var(CMAKE_C_COMPILER)
_probe_print_var(CMAKE_Fortran_COMPILER)
_probe_print_var(CMAKE_C_COMPILER_ID)
_probe_print_var(CMAKE_Fortran_COMPILER_ID)

set(_probe_findmpi_candidates "")
foreach(_probe_dir IN LISTS CMAKE_MODULE_PATH)
  if(EXISTS "${_probe_dir}/FindMPI.cmake")
    list(APPEND _probe_findmpi_candidates "${_probe_dir}/FindMPI.cmake")
  endif()
endforeach()
if(EXISTS "${CMAKE_ROOT}/Modules/FindMPI.cmake")
  list(APPEND _probe_findmpi_candidates "${CMAKE_ROOT}/Modules/FindMPI.cmake")
endif()
list(REMOVE_DUPLICATES _probe_findmpi_candidates)
message(STATUS "FindMPI candidates searched first-to-last: ${_probe_findmpi_candidates}")

message(STATUS "===== Running find_package(MPI) =====")
find_package(MPI REQUIRED COMPONENTS C Fortran)

message(STATUS "===== FindMPI Results =====")
foreach(_probe_var IN ITEMS
    MPI_FOUND
    MPI_VERSION
    MPI_VERSION_MAJOR
    MPI_VERSION_MINOR
    MPI_C_FOUND
    MPI_C_VERSION
    MPI_C_VERSION_MAJOR
    MPI_C_VERSION_MINOR
    MPI_C_COMPILER
    MPI_C_INCLUDE_DIRS
    MPI_C_LIBRARIES
    MPI_C_COMPILE_OPTIONS
    MPI_C_COMPILE_DEFINITIONS
    MPI_C_LINK_FLAGS
    MPI_C_HEADER_DIR
    MPI_C_WORKS
    MPI_C_LIBRARY_VERSION_STRING
    MPI_Fortran_FOUND
    MPI_Fortran_VERSION
    MPI_Fortran_VERSION_MAJOR
    MPI_Fortran_VERSION_MINOR
    MPI_Fortran_COMPILER
    MPI_Fortran_INCLUDE_DIRS
    MPI_Fortran_LIBRARIES
    MPI_Fortran_COMPILE_OPTIONS
    MPI_Fortran_COMPILE_DEFINITIONS
    MPI_Fortran_LINK_FLAGS
    MPI_Fortran_HAVE_F77_HEADER
    MPI_Fortran_HAVE_F90_MODULE
    MPI_Fortran_HAVE_F08_MODULE
    MPI_Fortran_F77_HEADER_SUBARRAYS
    MPI_Fortran_F90_MODULE_SUBARRAYS
    MPI_Fortran_F08_MODULE_SUBARRAYS
    MPI_Fortran_F77_HEADER_ASYNCPROT
    MPI_Fortran_F90_MODULE_ASYNCPROT
    MPI_Fortran_F08_MODULE_ASYNCPROT
    MPI_Fortran_LIBRARY_VERSION_STRING
    MPIEXEC_EXECUTABLE
    MPIEXEC_NUMPROC_FLAG
    MPIEXEC_PREFLAGS
    MPIEXEC_POSTFLAGS
    MPIEXEC_MAX_NUMPROCS)
  _probe_print_var(${_probe_var})
endforeach()

CheckMpiLargeCount(SCALAPACK_PROBE_HAVE_MPI_LARGE_COUNT_APIS
                   OUTPUT_LOG_VAR SCALAPACK_PROBE_MPI_LARGE_COUNT_LOG)

message(STATUS "===== Large-Count API Probe =====")
_probe_print_var(SCALAPACK_PROBE_HAVE_MPI_LARGE_COUNT_APIS)
message(STATUS "SCALAPACK_PROBE_MPI_LARGE_COUNT_LOG=[${SCALAPACK_PROBE_MPI_LARGE_COUNT_LOG}]")

message(STATUS "===== C Header Probe =====")
try_run(
  SCALAPACK_PROBE_C_RUN_RESULT
  SCALAPACK_PROBE_C_COMPILE_RESULT
  "${CMAKE_CURRENT_BINARY_DIR}"
  "${CMAKE_CURRENT_LIST_DIR}/print_mpi_c_version.c"
  LINK_LIBRARIES MPI::MPI_C
  COMPILE_OUTPUT_VARIABLE SCALAPACK_PROBE_C_COMPILE_OUTPUT
  RUN_OUTPUT_VARIABLE SCALAPACK_PROBE_C_RUN_OUTPUT
)
_probe_print_var(SCALAPACK_PROBE_C_COMPILE_RESULT)
_probe_print_var(SCALAPACK_PROBE_C_RUN_RESULT)
message(STATUS "SCALAPACK_PROBE_C_COMPILE_OUTPUT=[${SCALAPACK_PROBE_C_COMPILE_OUTPUT}]")
message(STATUS "SCALAPACK_PROBE_C_RUN_OUTPUT=[${SCALAPACK_PROBE_C_RUN_OUTPUT}]")

message(STATUS "===== Fortran Module Probe =====")
try_run(
  SCALAPACK_PROBE_Fortran_RUN_RESULT
  SCALAPACK_PROBE_Fortran_COMPILE_RESULT
  "${CMAKE_CURRENT_BINARY_DIR}"
  "${CMAKE_CURRENT_LIST_DIR}/print_mpi_fortran_version.f90"
  LINK_LIBRARIES MPI::MPI_Fortran
  COMPILE_OUTPUT_VARIABLE SCALAPACK_PROBE_Fortran_COMPILE_OUTPUT
  RUN_OUTPUT_VARIABLE SCALAPACK_PROBE_Fortran_RUN_OUTPUT
)
_probe_print_var(SCALAPACK_PROBE_Fortran_COMPILE_RESULT)
_probe_print_var(SCALAPACK_PROBE_Fortran_RUN_RESULT)
message(STATUS "SCALAPACK_PROBE_Fortran_COMPILE_OUTPUT=[${SCALAPACK_PROBE_Fortran_COMPILE_OUTPUT}]")
message(STATUS "SCALAPACK_PROBE_Fortran_RUN_OUTPUT=[${SCALAPACK_PROBE_Fortran_RUN_OUTPUT}]")

message(STATUS "===== Summary =====")
if(SCALAPACK_PROBE_HAVE_MPI_LARGE_COUNT_APIS)
  message(STATUS "BLACS large-count _c APIs: AVAILABLE")
else()
  message(STATUS "BLACS large-count _c APIs: NOT AVAILABLE")
endif()
