include_guard(GLOBAL)

macro(scalapack_configure_integer_abi)
  option(SCALAPACK_ENABLE_ILP64
         "Build ScaLAPACK with 8-byte default Fortran INTEGER and ILP64 BLAS/LAPACK"
         OFF)
  set(SCALAPACK_FORTRAN_INT_BYTES "" CACHE STRING
      "Detected default Fortran INTEGER width in bytes for the current compiler configuration (set by CMake; supported values: 4 or 8)")
  set(SCALAPACK_EXPECT_BLAS_INT_BYTES "" CACHE STRING
      "Advanced override for the expected BLAS integer width in bytes (auto-derived from SCALAPACK_ENABLE_ILP64 when empty; supported values: 4 or 8)")
  set(SCALAPACK_BLACS_BUFFER_SIZE_TYPE "size_t" CACHE STRING
      "C type used for BLACS internal packed-buffer byte lengths")

  if(SCALAPACK_ENABLE_ILP64)
     set(_scalapack_requested_int_bytes "8")
  else()
     set(_scalapack_requested_int_bytes "4")
  endif()

  if(DEFINED SCALAPACK_FORTRAN_INT_BYTES AND NOT "${SCALAPACK_FORTRAN_INT_BYTES}" STREQUAL "")
     if(NOT SCALAPACK_FORTRAN_INT_BYTES MATCHES "^(4|8)$")
        message(FATAL_ERROR
                "SCALAPACK_FORTRAN_INT_BYTES must be 4 or 8 when present in the cache.")
     endif()
     if(NOT "${SCALAPACK_FORTRAN_INT_BYTES}" STREQUAL "${_scalapack_requested_int_bytes}")
        message(STATUS
                "Ignoring cached SCALAPACK_FORTRAN_INT_BYTES=${SCALAPACK_FORTRAN_INT_BYTES}; "
                "the default Fortran INTEGER width is auto-detected from the current "
                "compiler configuration on each configure.")
     endif()
  endif()

  if(DEFINED SCALAPACK_EXPECT_BLAS_INT_BYTES AND NOT "${SCALAPACK_EXPECT_BLAS_INT_BYTES}" STREQUAL "")
     if(NOT SCALAPACK_EXPECT_BLAS_INT_BYTES MATCHES "^(4|8)$")
        message(FATAL_ERROR
                "SCALAPACK_EXPECT_BLAS_INT_BYTES must be 4 or 8 when set explicitly.")
     endif()
     if(NOT "${SCALAPACK_EXPECT_BLAS_INT_BYTES}" STREQUAL "${_scalapack_requested_int_bytes}")
        message(FATAL_ERROR
                "SCALAPACK_ENABLE_ILP64=${SCALAPACK_ENABLE_ILP64} conflicts with "
                "SCALAPACK_EXPECT_BLAS_INT_BYTES=${SCALAPACK_EXPECT_BLAS_INT_BYTES}. "
                "Prefer using SCALAPACK_ENABLE_ILP64 as the top-level ILP64 switch.")
     endif()
  else()
     set(SCALAPACK_EXPECT_BLAS_INT_BYTES "${_scalapack_requested_int_bytes}" CACHE STRING
         "Advanced override for the expected BLAS integer width in bytes (auto-derived from SCALAPACK_ENABLE_ILP64 when empty; supported values: 4 or 8)" FORCE)
  endif()

  if(SCALAPACK_ENABLE_ILP64)
     set(_scalapack_ilp64_fortran_flag "")
     if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|LLVMFlang|Flang|NAG")
        set(_scalapack_ilp64_fortran_flag "-fdefault-integer-8")
     elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM|NVHPC|PGI")
        set(_scalapack_ilp64_fortran_flag "-i8")
     elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
        set(_scalapack_ilp64_fortran_flag "-qintsize=8")
     else()
        message(FATAL_ERROR
                "SCALAPACK_ENABLE_ILP64=ON, but no default-8-byte-INTEGER compiler flag "
                "is known for Fortran compiler ID ${CMAKE_Fortran_COMPILER_ID}.")
     endif()

     check_fortran_compiler_flag("${_scalapack_ilp64_fortran_flag}"
                                 SCALAPACK_HAVE_ILP64_FORTRAN_FLAG)
     if(NOT SCALAPACK_HAVE_ILP64_FORTRAN_FLAG)
        message(FATAL_ERROR
                "SCALAPACK_ENABLE_ILP64=ON, but ${CMAKE_Fortran_COMPILER} does not accept "
                "the required flag ${_scalapack_ilp64_fortran_flag}.")
     endif()

     string(FIND " ${CMAKE_Fortran_FLAGS} " " ${_scalapack_ilp64_fortran_flag} "
            _scalapack_ilp64_flag_index)
     if(_scalapack_ilp64_flag_index EQUAL -1)
        string(STRIP "${CMAKE_Fortran_FLAGS} ${_scalapack_ilp64_fortran_flag}"
               CMAKE_Fortran_FLAGS)
     endif()
     message(STATUS
             "SCALAPACK_ENABLE_ILP64=ON: using Fortran flag ${_scalapack_ilp64_fortran_flag}")
  endif()

  CheckFortranIntegerABI(_scalapack_fortran_int_bytes)
  set(SCALAPACK_FORTRAN_INT_BYTES "${_scalapack_fortran_int_bytes}" CACHE STRING
      "Detected default Fortran INTEGER width in bytes for the current compiler configuration (set by CMake; supported values: 4 or 8)" FORCE)
  message(STATUS "Default Fortran INTEGER width: ${SCALAPACK_FORTRAN_INT_BYTES} bytes")

  set(_scalapack_mpi_version "")
  if(DEFINED MPI_VERSION AND NOT "${MPI_VERSION}" STREQUAL "")
     set(_scalapack_mpi_version "${MPI_VERSION}")
  elseif(DEFINED MPI_C_VERSION AND NOT "${MPI_C_VERSION}" STREQUAL "")
     set(_scalapack_mpi_version "${MPI_C_VERSION}")
  elseif(DEFINED MPI_Fortran_VERSION AND NOT "${MPI_Fortran_VERSION}" STREQUAL "")
     set(_scalapack_mpi_version "${MPI_Fortran_VERSION}")
  endif()

  set(_scalapack_mpi_library_version "")
  if(DEFINED MPI_C_LIBRARY_VERSION_STRING AND NOT "${MPI_C_LIBRARY_VERSION_STRING}" STREQUAL "")
     set(_scalapack_mpi_library_version "${MPI_C_LIBRARY_VERSION_STRING}")
  elseif(DEFINED MPI_Fortran_LIBRARY_VERSION_STRING AND NOT "${MPI_Fortran_LIBRARY_VERSION_STRING}" STREQUAL "")
     set(_scalapack_mpi_library_version "${MPI_Fortran_LIBRARY_VERSION_STRING}")
  endif()

  if(SCALAPACK_ENABLE_ILP64)
     CheckMpiLargeCount(SCALAPACK_HAVE_MPI_LARGE_COUNT_APIS
                        OUTPUT_LOG_VAR SCALAPACK_MPI_LARGE_COUNT_LOG)
     if(NOT SCALAPACK_HAVE_MPI_LARGE_COUNT_APIS)
        message(FATAL_ERROR
                "ILP64 builds require the MPI large-count _c APIs used by BLACS "
                "(for example MPI_Send_c, MPI_Type_vector_c, MPI_Pack_c, and "
                "MPI_Op_create_c), but the current MPI C interface does not provide "
                "them. Detected MPI standard version from CMake/headers: "
                "${_scalapack_mpi_version}. MPI library version string: "
                "${_scalapack_mpi_library_version}\n"
                "MPI large-count probe build log:\n${SCALAPACK_MPI_LARGE_COUNT_LOG}")
     endif()
     message(STATUS
             "ILP64 build enabled with MPI large-count _c support "
             "(reported MPI standard version: ${_scalapack_mpi_version})")
  endif()

  if(SCALAPACK_FORTRAN_INT_BYTES STREQUAL "4")
     set(SCALAPACK_API_INT_C_TYPE "int32_t")
  elseif(SCALAPACK_FORTRAN_INT_BYTES STREQUAL "8")
     set(SCALAPACK_API_INT_C_TYPE "int64_t")
  else()
     message(FATAL_ERROR "SCALAPACK_FORTRAN_INT_BYTES must be 4 or 8.")
  endif()

  if(SCALAPACK_EXPECT_BLAS_INT_BYTES STREQUAL "4")
     set(SCALAPACK_BLAS_INT_C_TYPE "int32_t")
  elseif(SCALAPACK_EXPECT_BLAS_INT_BYTES STREQUAL "8")
     set(SCALAPACK_BLAS_INT_C_TYPE "int64_t")
  else()
     message(FATAL_ERROR "SCALAPACK_EXPECT_BLAS_INT_BYTES must be 4 or 8.")
  endif()

  if(DEFINED BLA_SIZEOF_INTEGER AND NOT "${BLA_SIZEOF_INTEGER}" STREQUAL "")
     if(NOT "${BLA_SIZEOF_INTEGER}" STREQUAL "${SCALAPACK_EXPECT_BLAS_INT_BYTES}")
        message(FATAL_ERROR
                "BLA_SIZEOF_INTEGER=${BLA_SIZEOF_INTEGER} conflicts with "
                "SCALAPACK_EXPECT_BLAS_INT_BYTES=${SCALAPACK_EXPECT_BLAS_INT_BYTES}. "
                "Use matching values so CMake searches for a BLAS/LAPACK with the "
                "expected integer ABI.")
     endif()
  else()
     set(BLA_SIZEOF_INTEGER "${SCALAPACK_EXPECT_BLAS_INT_BYTES}" CACHE STRING
         "Requested BLAS/LAPACK integer size for CMake FindBLAS/FindLAPACK." FORCE)
  endif()
  message(STATUS "Requesting BLAS/LAPACK integer ABI: ${BLA_SIZEOF_INTEGER} bytes")

  if(NOT SCALAPACK_FORTRAN_INT_BYTES STREQUAL SCALAPACK_EXPECT_BLAS_INT_BYTES)
     set(_scalapack_fortran_i8_hint
         "Reconfigure with a compiler setup whose default INTEGER is 8 bytes, or set SCALAPACK_FORTRAN_INT_BYTES explicitly if auto-detection is not possible.")
     if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|LLVMFlang|Flang|NAG")
        set(_scalapack_fortran_i8_hint
            "GNU/Flang-style compilers default to 4-byte INTEGER unless you add -fdefault-integer-8. Reconfigure with -DCMAKE_Fortran_FLAGS=-fdefault-integer-8 (or pass that flag through your mpif90 wrapper) so ScaLAPACK and the ILP64 BLAS use the same INTEGER ABI.")
     elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM|NVHPC|PGI")
        set(_scalapack_fortran_i8_hint
            "Intel/NVHPC-style compilers usually need -i8 to make default INTEGER 8 bytes. Reconfigure with -DCMAKE_Fortran_FLAGS=-i8 (or pass that flag through your mpif90 wrapper) so ScaLAPACK and the ILP64 BLAS use the same INTEGER ABI.")
     elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
        set(_scalapack_fortran_i8_hint
            "XL Fortran usually needs -qintsize=8 to make default INTEGER 8 bytes. Reconfigure with -DCMAKE_Fortran_FLAGS=-qintsize=8 so ScaLAPACK and the ILP64 BLAS use the same INTEGER ABI.")
     endif()
     message(FATAL_ERROR
             "Mixed-width ScaLAPACK/BLAS builds are not supported yet: "
             "Fortran INTEGER is ${SCALAPACK_FORTRAN_INT_BYTES} bytes, "
             "but SCALAPACK_EXPECT_BLAS_INT_BYTES is ${SCALAPACK_EXPECT_BLAS_INT_BYTES}. "
             "The Fortran INTEGER width was auto-detected from the current compiler "
             "configuration for `${CMAKE_Fortran_COMPILER}` (${CMAKE_Fortran_COMPILER_ID}). "
             "${_scalapack_fortran_i8_hint}")
  endif()

  if(SCALAPACK_ENABLE_ILP64 AND NOT SCALAPACK_FORTRAN_INT_BYTES STREQUAL "8")
     message(FATAL_ERROR
             "SCALAPACK_ENABLE_ILP64=ON, but the current Fortran compiler configuration "
             "still defaults to ${SCALAPACK_FORTRAN_INT_BYTES}-byte INTEGER. "
             "Check the ILP64 compiler flag handling in CMake or your mpif90 wrapper.")
  endif()

  check_c_source_compiles(
"#include <stddef.h>
typedef ${SCALAPACK_BLACS_BUFFER_SIZE_TYPE} scalapack_buflen_t;
typedef char scalapack_buflen_must_fit_size_t[(sizeof(scalapack_buflen_t) <= sizeof(size_t)) ? 1 : -1];
int main(void) {
  scalapack_buflen_t x = 0;
  return (int) x;
}
" SCALAPACK_BLACS_BUFFER_SIZE_TYPE_VALID)

  if(NOT SCALAPACK_BLACS_BUFFER_SIZE_TYPE_VALID)
     message(FATAL_ERROR
             "SCALAPACK_BLACS_BUFFER_SIZE_TYPE='${SCALAPACK_BLACS_BUFFER_SIZE_TYPE}' is not a valid C type or is wider than size_t.")
  endif()

  configure_file(${SCALAPACK_SOURCE_DIR}/CMAKE/scalapack-types.h.in
                 ${SCALAPACK_BINARY_DIR}/scalapack-types.h
                 @ONLY)
  include_directories(${SCALAPACK_BINARY_DIR})
  include_directories(${SCALAPACK_SOURCE_DIR}/TOOLS)
endmacro()
