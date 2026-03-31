cmake_minimum_required(VERSION 3.26...4.0)

project(SCALAPACK VERSION 2.2.2 LANGUAGES C Fortran)
set(CMAKE_C_STANDARD 17)
set(CMAKE_C_EXTENSIONS OFF)
include(CheckCSourceCompiles)
include(CheckCCompilerFlag)
include(CheckFortranCompilerFlag)

# Add the CMake directory for custon CMake modules
set(CMAKE_MODULE_PATH "${SCALAPACK_SOURCE_DIR}/CMAKE" ${CMAKE_MODULE_PATH})
include(CheckFortranIntegerABI)

if (UNIX)
   if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
      set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fltconsistency -fp_port" )
   endif ()
endif ()

#
#   MPI
#

find_package(MPI)
if (MPI_FOUND)
   message(STATUS "Found MPI_LIBRARY : ${MPI_FOUND} ")

   find_program(MPI_C_COMPILER
      NAMES mpicc
      HINTS "${MPI_BASE_DIR}"
      PATH_SUFFIXES bin
      DOC "MPI C compiler.")
   MARK_AS_ADVANCED(MPI_C_COMPILER)
   if ("${MPI_C_COMPILER}" STREQUAL "MPI_C_COMPILER-NOTFOUND")
      message(ERROR "--> MPI C Compiler NOT FOUND (please set MPI_BASE_DIR accordingly")
   else()
      message(STATUS "--> MPI C Compiler : ${MPI_C_COMPILER}")
      #SET(CMAKE_C_COMPILER "${MPI_C_COMPILER}")
      message(STATUS "--> C Compiler : ${CMAKE_C_COMPILER}")
   endif()
   find_program(MPI_Fortran_COMPILER
      NAMES mpif77
      HINTS "${MPI_BASE_DIR}"
      PATH_SUFFIXES bin
      DOC "MPI Fortran compiler.")
   MARK_AS_ADVANCED(MPI_Fortran_COMPILER)

   if ("${MPI_Fortran_COMPILER}" STREQUAL "MPI_Fortran_COMPILER-NOTFOUND")
      message(ERROR "--> MPI Fortran Compiler NOT FOUND (please set MPI_BASE_DIR accordingly")
   else()
      message(STATUS "--> MPI Fortran Compiler : ${MPI_Fortran_COMPILER}")
      SET(Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")
      #SET(CMAKE_Fortran_COMPILER "${MPI_Fortran_COMPILER}")
      message(STATUS "--> Fortran Compiler : ${CMAKE_Fortran_COMPILER}")
   endif()

else()
   message(STATUS "Found MPI_LIBRARY : ${MPI_FOUND} ")
   set(MPI_BASE_DIR ${MPI_BASE_DIR} CACHE PATH "MPI Path")
   UNSET(MPIEXEC CACHE)
   UNSET(MPIEXEC_POSTFLAGS CACHE)
   UNSET(MPIEXEC_PREFLAGS CACHE)
   UNSET(MPIEXEC_MAX_NUMPROCS CACHE)
   UNSET(MPIEXEC_NUMPROC_FLAG CACHE)
   UNSET(MPI_COMPILE_FLAGS CACHE)
   UNSET(MPI_LINK_FLAGS CACHE)
   UNSET(MPI_INCLUDE_PATH CACHE)
   message(FATAL_ERROR "--> MPI Library NOT FOUND -- please set MPI_BASE_DIR accordingly --")
endif()

message(STATUS "--> MPI C Compiler : ${MPI_C_COMPILER}")
message(STATUS "--> C Compiler : ${CMAKE_C_COMPILER}")
message(STATUS "--> MPI Fortran Compiler : ${MPI_Fortran_COMPILER}")
message(STATUS "--> Fortran Compiler : ${CMAKE_Fortran_COMPILER}")

# --------------------------------------------------
# ILP64 Support via MPI 4.0+
#
# ScaLAPACK will automatically use MPI 4.0+ large count support if available.
# MPI 4.0+ provides native support for counts > INT_MAX via _c variants.
#
# Note: BigMPI is NOT supported because it lacks MPI_Op_create_c functionality
# which is required for ScaLAPACK's custom reduction operations.
message(STATUS "ILP64 support: Will use MPI 4.0+ _c variants if MPI_VERSION >= 4")

if (UNIX)
   if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
      set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fltconsistency -fp_port" )
   endif ()
endif ()

macro(SCALAPACK_install_library lib)
  install(TARGETS ${lib} EXPORT scalapack-targets
    RUNTIME DESTINATION Testing
  )
endmacro()

# --------------------------------------------------
# Testing
SET(DART_TESTING_TIMEOUT 600)

enable_testing()
# --------------------------------------------------

# Organize output files.  On Windows this also keeps .dll files next
# to the .exe files that need them, making tests easy to run.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCALAPACK_BINARY_DIR}/TESTING)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${SCALAPACK_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SCALAPACK_BINARY_DIR}/lib)

# --------------------------------------------------
# By default static library
function (invertBoolean varName varValue)
  if (${varValue})
    set(${varName} false PARENT_SCOPE)
  else ()
    set(${varName} true PARENT_SCOPE)
  endif ()
endfunction ()
option(BUILD_SHARED_LIBS "Build shared libraries" OFF )
invertBoolean("BUILD_STATIC_LIBS" ${BUILD_SHARED_LIBS})
if ((${BUILD_SHARED_LIBS} EQUAL ON) AND NOT CMAKE_POSITION_INDEPENDENT_CODE)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()
option(SCALAPACK_ENABLE_SANITIZERS
       "Enable ASan/UBSan for C code in Debug builds"
       OFF)
option(SCALAPACK_ENABLE_FORTRAN_RUNTIME_CHECKS
       "Enable extra GNU Fortran runtime checks in Debug builds"
       OFF)
set(SCALAPACK_SANITIZERS "address;undefined" CACHE STRING
    "Semicolon-separated sanitizer list for Debug builds")
option(SCALAPACK_ENABLE_STRICT_C_WARNINGS
       "Enable a curated Clang warning bundle for C sources"
       OFF)
set(SCALAPACK_STRICT_C_WARNING_FLAGS
    "-Wall;-Wextra;-Wconversion;-Wsign-conversion;-Wshadow;-Wundef;-Wformat=2;-Wcast-qual;-Wwrite-strings;-Wstrict-prototypes;-Wmissing-prototypes;-Wimplicit-fallthrough;-Wnull-dereference;-Wdouble-promotion"
    CACHE STRING
    "Semicolon-separated C warning flags to enable when SCALAPACK_ENABLE_STRICT_C_WARNINGS=ON")

if(SCALAPACK_ENABLE_STRICT_C_WARNINGS)
  if(CMAKE_C_COMPILER_ID MATCHES "^(AppleClang|Clang)$")
    set(_scalapack_strict_c_warning_flags)
    foreach(_scalapack_warning_flag IN LISTS SCALAPACK_STRICT_C_WARNING_FLAGS)
      string(MAKE_C_IDENTIFIER "${_scalapack_warning_flag}" _scalapack_warning_ident)
      check_c_compiler_flag("${_scalapack_warning_flag}"
                            "SCALAPACK_HAVE_C_WARNING_${_scalapack_warning_ident}")
      if(SCALAPACK_HAVE_C_WARNING_${_scalapack_warning_ident})
        list(APPEND _scalapack_strict_c_warning_flags "${_scalapack_warning_flag}")
      else()
        message(STATUS
                "Skipping unsupported C warning flag: ${_scalapack_warning_flag}")
      endif()
    endforeach()

    if(_scalapack_strict_c_warning_flags)
      add_compile_options(
        $<$<COMPILE_LANGUAGE:C>:${_scalapack_strict_c_warning_flags}>
      )
    else()
      message(WARNING
              "SCALAPACK_ENABLE_STRICT_C_WARNINGS is ON, but no requested warning flags are supported.")
    endif()
  else()
    message(STATUS
            "Ignoring SCALAPACK_ENABLE_STRICT_C_WARNINGS because the C compiler is ${CMAKE_C_COMPILER_ID}, not Clang.")
  endif()
endif()

if(SCALAPACK_ENABLE_SANITIZERS)
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING
            "SCALAPACK_ENABLE_SANITIZERS is most useful with Debug builds.")
  endif()

  string(REPLACE ";" "," _scalapack_sanitizer_list "${SCALAPACK_SANITIZERS}")
  foreach(_scalapack_sanitizer IN LISTS SCALAPACK_SANITIZERS)
    set(CMAKE_REQUIRED_FLAGS "-fsanitize=${_scalapack_sanitizer}")
    set(CMAKE_REQUIRED_LINK_OPTIONS "-fsanitize=${_scalapack_sanitizer}")
    check_c_source_compiles(
"int main(void) { return 0; }"
      "SCALAPACK_HAVE_C_SANITIZER_${_scalapack_sanitizer}")
    unset(CMAKE_REQUIRED_FLAGS)
    unset(CMAKE_REQUIRED_LINK_OPTIONS)
    if(NOT SCALAPACK_HAVE_C_SANITIZER_${_scalapack_sanitizer})
      message(FATAL_ERROR
              "Requested C sanitizer flag is not supported: -fsanitize=${_scalapack_sanitizer}")
    endif()
  endforeach()

  add_compile_options(
    $<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Debug>>:-O1>
    $<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Debug>>:-g>
    $<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Debug>>:-fno-omit-frame-pointer>
    $<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Debug>>:-fsanitize=${_scalapack_sanitizer_list}>
  )
  add_link_options(
    $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
    $<$<CONFIG:Debug>:-fsanitize=${_scalapack_sanitizer_list}>
  )
endif()

if(SCALAPACK_ENABLE_FORTRAN_RUNTIME_CHECKS)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    check_fortran_compiler_flag("-fcheck=all" SCALAPACK_HAVE_FORTRAN_FCHECK_ALL)
    check_fortran_compiler_flag("-fbacktrace" SCALAPACK_HAVE_FORTRAN_FBACKTRACE)
    check_fortran_compiler_flag("-ffpe-trap=invalid,zero,overflow"
                                SCALAPACK_HAVE_FORTRAN_FFPE_TRAP)
    check_fortran_compiler_flag("-finit-real=snan" SCALAPACK_HAVE_FORTRAN_FINIT_REAL_SNAN)

    add_compile_options(
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O1>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-g>
    )

    if(SCALAPACK_HAVE_FORTRAN_FCHECK_ALL)
      add_compile_options(
        $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fcheck=all>
      )
    endif()
    if(SCALAPACK_HAVE_FORTRAN_FBACKTRACE)
      add_compile_options(
        $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace>
      )
    endif()
    if(SCALAPACK_HAVE_FORTRAN_FFPE_TRAP)
      add_compile_options(
        $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-ffpe-trap=invalid,zero,overflow>
      )
    endif()
    if(SCALAPACK_HAVE_FORTRAN_FINIT_REAL_SNAN)
      add_compile_options(
        $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-finit-real=snan>
      )
    endif()
  else()
    message(WARNING
            "SCALAPACK_ENABLE_FORTRAN_RUNTIME_CHECKS currently only adds GNU Fortran flags.")
  endif()
endif()

# --------------------------------------------------
# Check for any necessary platform specific compiler flags
include( CheckBLACSCompilerFlags )
CheckBLACSCompilerFlags()

include(GNUInstallDirs)
set(PKG_CONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig CACHE PATH "pkg-config install path")

# --------------------------------------------------
# BLACS Internal variables
#
#   Fortran Mangling, MPI Tests and BLACS settings
#
include(FortranMangling)

#MESSAGE(STATUS "Setting CDEFS = ${CDEFS}")
#set(CDEFS ${CDEFS} CACHE STRING "Fortran Mangling" FORCE)
MESSAGE(STATUS "=========")

# --------------------------------------------------
# Compiler Flags
ADD_DEFINITIONS( "-D${CDEFS}")

# --------------------------------------------------
# ScaLAPACK needs BLAS and LAPACK

option(USE_OPTIMIZED_LAPACK_BLAS "Whether or not to search for optimized LAPACK and BLAS libraries on your machine (if not found, Reference LAPACK and BLAS will be downloaded and installed)" ON)
option(SCALAPACK_ENABLE_ILP64
       "Build ScaLAPACK with 8-byte default Fortran INTEGER and ILP64 BLAS/LAPACK"
       OFF)
option(SCALAPACK_PREFER_OPENBLAS_ON_APPLE_SILICON
       "Prefer OpenBLAS over Accelerate on macOS arm64 unless the user explicitly selects another BLAS/LAPACK backend"
       ON)
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

if(SCALAPACK_ENABLE_ILP64)
   if(_scalapack_mpi_version STREQUAL "" OR _scalapack_mpi_version VERSION_LESS 4)
      message(FATAL_ERROR
              "ILP64 builds require MPI 4.0 or newer because BLACS relies on MPI "
              "large-count _c APIs. Detected MPI version: ${_scalapack_mpi_version}.")
   endif()
   message(STATUS "ILP64 build enabled with MPI ${_scalapack_mpi_version} large-count support")
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

message(STATUS "CHECKING BLAS AND LAPACK LIBRARIES")
include(CheckFortranFunctionExists)

set(_scalapack_force_openblas FALSE)
set(_scalapack_openblas_formula "openblas")
set(_scalapack_openblas_pkgconfig_module "openblas")
set(_scalapack_openblas_example_package "openblas")
set(_scalapack_openblas_direct_lib "")
if(SCALAPACK_EXPECT_BLAS_INT_BYTES STREQUAL "8")
   set(_scalapack_openblas_formula "openblas64")
   set(_scalapack_openblas_pkgconfig_module "openblas64")
   set(_scalapack_openblas_example_package "openblas64")
   set(_scalapack_openblas_direct_lib "/opt/homebrew/opt/openblas64/lib/libopenblas64_.dylib")
else()
   set(_scalapack_openblas_direct_lib "/opt/homebrew/opt/openblas/lib/libopenblas.dylib")
endif()
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin"
   AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$"
   AND SCALAPACK_PREFER_OPENBLAS_ON_APPLE_SILICON
   AND NOT LAPACK_LIBRARIES
   AND NOT DEFINED BLA_VENDOR)
   set(_scalapack_force_openblas TRUE)
   set(BLA_VENDOR OpenBLAS)
   set(BLA_PREFER_PKGCONFIG ON)
   set(BLA_PKGCONFIG_BLAS "${_scalapack_openblas_pkgconfig_module}")
   set(BLA_PKGCONFIG_LAPACK "${_scalapack_openblas_pkgconfig_module}")
   foreach(_scalapack_openblas_hint
           /opt/homebrew/opt/${_scalapack_openblas_formula}
           /opt/homebrew/opt/${_scalapack_openblas_formula}/lib
           /opt/homebrew/opt/${_scalapack_openblas_formula}/lib/pkgconfig
           /opt/local/lib/openblas
           /opt/local/lib)
      if(EXISTS "${_scalapack_openblas_hint}")
         list(PREPEND CMAKE_PREFIX_PATH "${_scalapack_openblas_hint}")
      endif()
   endforeach()
   if(EXISTS "/opt/homebrew/opt/${_scalapack_openblas_formula}/lib/pkgconfig")
      if(DEFINED ENV{PKG_CONFIG_PATH} AND NOT "$ENV{PKG_CONFIG_PATH}" STREQUAL "")
         set(ENV{PKG_CONFIG_PATH}
             "/opt/homebrew/opt/${_scalapack_openblas_formula}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
      else()
         set(ENV{PKG_CONFIG_PATH}
             "/opt/homebrew/opt/${_scalapack_openblas_formula}/lib/pkgconfig")
      endif()
   endif()
   message(STATUS
           "macOS arm64 detected: preferring ${_scalapack_openblas_formula} over Accelerate. "
           "Install ${_scalapack_openblas_example_package} or set "
           "-DSCALAPACK_PREFER_OPENBLAS_ON_APPLE_SILICON=OFF to allow Accelerate.")
endif()

IF(LAPACK_LIBRARIES)
  message(STATUS "--> LAPACK supplied by user is ${LAPACK_LIBRARIES}.")
  set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
  CHECK_FORTRAN_FUNCTION_EXISTS("dgesv" LAPACK_FOUND)
  unset( CMAKE_REQUIRED_LIBRARIES )
  message(STATUS "--> LAPACK routine dgesv is found: ${LAPACK_FOUND}.")
ENDIF()

if(LAPACK_FOUND)
      message(STATUS "--> LAPACK supplied by user is WORKING, will use ${LAPACK_LIBRARIES}.")
else(LAPACK_FOUND)
   if(USE_OPTIMIZED_LAPACK_BLAS)
      message(STATUS "--> Searching for optimized LAPACK and BLAS libraries on your machine.")
      find_package(LAPACK)
      if(_scalapack_force_openblas
         AND NOT LAPACK_FOUND
         AND EXISTS "${_scalapack_openblas_direct_lib}")
         set(LAPACK_LIBRARIES "${_scalapack_openblas_direct_lib}")
         set(BLAS_LIBRARIES "${_scalapack_openblas_direct_lib}")
         set(CMAKE_REQUIRED_LIBRARIES "${LAPACK_LIBRARIES}")
         CHECK_FORTRAN_FUNCTION_EXISTS("dgesv" LAPACK_FOUND)
         unset(CMAKE_REQUIRED_LIBRARIES)
         if(LAPACK_FOUND)
            message(STATUS
                    "--> Using Homebrew ${_scalapack_openblas_formula} fallback: ${LAPACK_LIBRARIES}")
         endif()
      endif()
      if(_scalapack_force_openblas AND NOT LAPACK_FOUND)
         message(FATAL_ERROR
                 "macOS arm64 builds default to ${_scalapack_openblas_formula} because Accelerate can cause "
                 "pathological eigensolver test behavior on Apple Silicon. "
                 "Install ${_scalapack_openblas_example_package} (for example: brew install ${_scalapack_openblas_example_package}) or configure "
                 "with -DSCALAPACK_PREFER_OPENBLAS_ON_APPLE_SILICON=OFF to allow "
                 "CMake to fall back to Accelerate.")
      endif()
   ENDIF(USE_OPTIMIZED_LAPACK_BLAS)
   if(NOT LAPACK_FOUND)
      message(STATUS "--> LAPACK and BLAS were not found. Reference LAPACK and BLAS will be downloaded and installed")
      include(ExternalProject)
      ExternalProject_Add(
        lapack
        URL http://www.netlib.org/lapack/lapack.tgz
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${SCALAPACK_BINARY_DIR}
        PREFIX ${SCALAPACK_BINARY_DIR}/dependencies
      )
      if (NOT MSVC)
         SET(LAPACK_LIBRARIES ${SCALAPACK_BINARY_DIR}/lib/liblapack.a CACHE STRING "LAPACK library" FORCE)
         SET(BLAS_LIBRARIES ${SCALAPACK_BINARY_DIR}/lib/libblas.a CACHE STRING "BLAS library" FORCE)
      else () # MSVC
         SET(LAPACK_LIBRARIES ${SCALAPACK_BINARY_DIR}/lib/liblapack.lib CACHE STRING "LAPACK library" FORCE)
         SET(BLAS_LIBRARIES ${SCALAPACK_BINARY_DIR}/lib/libblas.lib CACHE STRING "BLAS library" FORCE)
      endif ()
   ENDIF()
ENDIF()

message(STATUS "BLAS library: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK library: ${LAPACK_LIBRARIES}")

if(LAPACK_FOUND)
   include(CheckBlasIntegerABI)
   CheckBlasIntegerABI(${SCALAPACK_EXPECT_BLAS_INT_BYTES})
endif()

set(SCALAPACK_BLAS_SYMBOL_ALIAS_FILE "")
if(APPLE
   AND LAPACK_FOUND
   AND SCALAPACK_BLAS_SYMBOL_STYLE STREQUAL "suffix64")
   set(SCALAPACK_BLAS_SYMBOL_ALIAS_FILE
       "${SCALAPACK_BINARY_DIR}/scalapack-blas-symbol-aliases.txt")
   set(_scalapack_alias_inputs "${BLAS_LIBRARIES};${LAPACK_LIBRARIES}")
   file(WRITE "${SCALAPACK_BLAS_SYMBOL_ALIAS_FILE}"
        "# Generated during the build after libscalapack.a is available.\n")
   set(_scalapack_blas_alias_link_option_build
       "-Wl,-alias_list,${SCALAPACK_BLAS_SYMBOL_ALIAS_FILE}")
   set(_scalapack_blas_alias_link_option_install
       "-Wl,-alias_list,${CMAKE_INSTALL_FULL_LIBDIR}/cmake/scalapack-${SCALAPACK_VERSION}/scalapack-blas-symbol-aliases.txt")
   message(STATUS
           "Darwin BLAS/LAPACK alias list will be generated at build time: ${SCALAPACK_BLAS_SYMBOL_ALIAS_FILE}")
endif()

MESSAGE(STATUS "=========")

# --------------------------------------------------
# Subdirectories that need to be processed

macro(append_subdir_files variable dirname)
get_directory_property(holder DIRECTORY ${dirname} DEFINITION ${variable})
foreach(depfile ${holder})
  list(APPEND ${variable} "${dirname}/${depfile}")
endforeach()
endmacro()

function(scalapack_apply_blas_symbol_aliasing target)
  if(NOT APPLE OR "${SCALAPACK_BLAS_SYMBOL_ALIAS_FILE}" STREQUAL "")
    return()
  endif()

  target_link_options(${target} PUBLIC
    $<BUILD_INTERFACE:${_scalapack_blas_alias_link_option_build}>
    $<INSTALL_INTERFACE:${_scalapack_blas_alias_link_option_install}>
  )
endfunction()


#
# BLACS
#
add_subdirectory(BLACS)
append_subdir_files(blacs "BLACS/SRC")

#
# TOOLS
#
add_subdirectory(TOOLS)
append_subdir_files(tools TOOLS)
append_subdir_files(tools-C TOOLS)
append_subdir_files(extra_lapack "TOOLS/LAPACK")

#
# PBLAS
#
add_subdirectory(PBLAS)
append_subdir_files(pblas "PBLAS/SRC")
append_subdir_files(pblas-F "PBLAS/SRC")
append_subdir_files(pbblas "PBLAS/SRC/PBBLAS")
append_subdir_files(ptzblas "PBLAS/SRC/PTZBLAS")
append_subdir_files(ptools "PBLAS/SRC/PTOOLS")

#
# REDIST
#
add_subdirectory(REDIST)
append_subdir_files(redist "REDIST/SRC")

#
# SRC
#
add_subdirectory(SRC)
append_subdir_files(src "SRC")
append_subdir_files(src-C "SRC")

if (NOT MSVC)
   set_source_files_properties(TOOLS/SL_gridreshape.c PROPERTIES COMPILE_OPTIONS "-Wno-error=implicit-function-declaration")
   add_library(scalapack ${blacs} ${tools} ${tools-C} ${extra_lapack} ${pblas} ${pblas-F} ${ptzblas} ${ptools} ${pbblas} ${redist} ${src} ${src-C})
   set_target_properties(scalapack PROPERTIES
    VERSION ${SCALAPACK_VERSION}
    SOVERSION ${SCALAPACK_VERSION_MAJOR}.${SCALAPACK_VERSION_MINOR})
   target_link_libraries( scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_C)
   if(APPLE AND NOT "${SCALAPACK_BLAS_SYMBOL_ALIAS_FILE}" STREQUAL "")
      add_custom_command(TARGET scalapack POST_BUILD
        COMMAND ${CMAKE_COMMAND}
                "-DINPUTS=${_scalapack_alias_inputs}"
                "-DREQUIRED_INPUTS=$<TARGET_FILE:scalapack>"
                "-DOUTPUT=${SCALAPACK_BLAS_SYMBOL_ALIAS_FILE}"
                "-DNM=${CMAKE_NM}"
                -P "${SCALAPACK_SOURCE_DIR}/CMAKE/GenerateDarwinBlasAliasList.cmake"
        BYPRODUCTS "${SCALAPACK_BLAS_SYMBOL_ALIAS_FILE}"
        VERBATIM)
   endif()
   scalapack_apply_blas_symbol_aliasing(scalapack)
   scalapack_install_library(scalapack)
else () # Need to separate Fortran and C Code
   OPTION(BUILD_SHARED_LIBS "Build shared libraries" ON )
   add_library(scalapack ${blacs} ${tools-C} ${pblas} ${ptools} ${redist} ${src-C})
   set_target_properties(scalapack PROPERTIES
    VERSION ${SCALAPACK_VERSION}
    SOVERSION ${SCALAPACK_VERSION_MAJOR}.${SCALAPACK_VERSION_MINOR})
   target_link_libraries( scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_C)
   scalapack_apply_blas_symbol_aliasing(scalapack)
   add_library(scalapack-F ${pblas-F} ${pbblas} ${ptzblas} ${tools} ${src} ${extra_lapack} )
   set_target_properties(scalapack-F PROPERTIES
    VERSION ${SCALAPACK_VERSION}
    SOVERSION ${SCALAPACK_VERSION_MAJOR}.${SCALAPACK_VERSION_MINOR})
   target_link_libraries( scalapack-F ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
   scalapack_apply_blas_symbol_aliasing(scalapack-F)
   scalapack_install_library(scalapack)
   scalapack_install_library(scalapack-F)
endif ()

if(TARGET lapack)
   add_dependencies(scalapack lapack)
endif()

option(SCALAPACK_BUILD_TESTS "Build all tests of the ScaLAPACK library" ON)
if(${SCALAPACK_BUILD_TESTS})
  add_subdirectory(TESTING)
endif()

# --------------------------------------------------
# CPACK Packaging

SET(CPACK_PACKAGE_NAME "ScaLAPACK")
SET(CPACK_PACKAGE_VENDOR "University of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ScaLAPACK- Linear Algebra Package")
set(CPACK_PACKAGE_VERSION_MAJOR ${SCALAPACK_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${SCALAPACK_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${SCALAPACK_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "SCALAPACK")
IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\http://icl.cs.utk.edu/lapack-forum")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.netlib.org/scalapack")
  SET(CPACK_NSIS_CONTACT "scalapack@eecs.utk.edu")
  SET(CPACK_NSIS_MODIFY_PATH ON)
  SET(CPACK_NSIS_DISPLAY_NAME "SCALAPACK-${SCALAPACK_VERSION}")
  set(CPACK_PACKAGE_RELOCATABLE "true")
ELSE(WIN32 AND NOT UNIX)
  SET(CPACK_GENERATOR "TGZ")
  SET(CPACK_SOURCE_GENERATOR TGZ)
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "scalapack-${SCALAPACK_VERSION}" )
  SET(CPACK_SOURCE_IGNORE_FILES ~$ .svn ${CPACK_SOURCE_IGNORE_FILES} )
ENDIF(WIN32 AND NOT UNIX)
INCLUDE(CPack)


# --------------------------------------------------


export(TARGETS scalapack FILE scalapack-targets.cmake)

if( NOT LAPACK_FOUND )
 install(FILES
  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
   )
endif( NOT LAPACK_FOUND )

configure_file(${SCALAPACK_SOURCE_DIR}/CMAKE/scalapack-config-version.cmake.in
  ${SCALAPACK_BINARY_DIR}/scalapack-config-version.cmake @ONLY)
configure_file(${SCALAPACK_SOURCE_DIR}/CMAKE/scalapack-config-build.cmake.in
  ${SCALAPACK_BINARY_DIR}/scalapack-config.cmake @ONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scalapack.pc.in ${CMAKE_CURRENT_BINARY_DIR}/scalapack.pc)
 install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/scalapack.pc
  DESTINATION ${PKG_CONFIG_DIR}
   )

configure_file(${SCALAPACK_SOURCE_DIR}/CMAKE/scalapack-config-install.cmake.in
  ${SCALAPACK_BINARY_DIR}/CMakeFiles/scalapack-config.cmake @ONLY)
install(FILES
  ${SCALAPACK_BINARY_DIR}/CMakeFiles/scalapack-config.cmake
  ${SCALAPACK_BINARY_DIR}/scalapack-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scalapack-${SCALAPACK_VERSION}
  )

if("${SCALAPACK_BLAS_SYMBOL_ALIAS_FILE}" STREQUAL "")
else()
 install(FILES
  ${SCALAPACK_BLAS_SYMBOL_ALIAS_FILE}
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scalapack-${SCALAPACK_VERSION}
   )
endif()

install(EXPORT scalapack-targets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scalapack-${SCALAPACK_VERSION})
