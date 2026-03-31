include(CheckCSourceCompiles)
include(CheckCSourceRuns)

function(CheckBlasIntegerABI expected_bytes)
  if(NOT expected_bytes MATCHES "^(4|8)$")
    message(FATAL_ERROR "SCALAPACK_EXPECT_BLAS_INT_BYTES must be 4 or 8.")
  endif()

  if(CMAKE_CROSSCOMPILING)
    message(STATUS "Skipping OpenBLAS integer ABI check while cross-compiling.")
    return()
  endif()

  set(_saved_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
  set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  set(_scalapack_blas_symbol_style "plain")

  unset(SCALAPACK_HAVE_BLAS_UNDERSCORE_SYMBOL CACHE)
  unset(SCALAPACK_HAVE_BLAS_64_UNDERSCORE_SYMBOL CACHE)

  check_c_source_compiles(
"extern void sgemm_(void);
int main(void) {
  sgemm_();
  return 0;
}
" SCALAPACK_HAVE_BLAS_UNDERSCORE_SYMBOL)

  check_c_source_compiles(
"extern void sgemm_64_(void);
int main(void) {
  sgemm_64_();
  return 0;
}
" SCALAPACK_HAVE_BLAS_64_UNDERSCORE_SYMBOL)

  if(expected_bytes STREQUAL "8"
     AND NOT SCALAPACK_HAVE_BLAS_UNDERSCORE_SYMBOL
     AND SCALAPACK_HAVE_BLAS_64_UNDERSCORE_SYMBOL)
    set(_scalapack_blas_symbol_style "suffix64")
    if(APPLE)
      message(STATUS
        "Detected an ILP64 BLAS/LAPACK that exports suffixed Fortran symbols like "
        "sgemm_64_. ScaLAPACK will enable Darwin linker aliases so unsuffixed calls "
        "such as sgemm_ resolve correctly.")
    else()
      message(FATAL_ERROR
        "Detected an ILP64 BLAS/LAPACK that exports suffixed Fortran symbols like "
        "sgemm_64_ instead of the unsuffixed names ScaLAPACK currently calls "
        "(for example sgemm_). This build is not link-compatible with that BLAS "
        "yet; a BLAS/LAPACK symbol translation layer is still needed on this "
        "platform.")
    endif()
  endif()

  unset(SCALAPACK_HAVE_OPENBLAS_GET_CONFIG CACHE)
  check_c_source_compiles(
"extern const char *openblas_get_config(void);
int main(void) {
  return openblas_get_config() == 0;
}
" SCALAPACK_HAVE_OPENBLAS_GET_CONFIG)

  if(SCALAPACK_HAVE_OPENBLAS_GET_CONFIG)
    set(_scalapack_openblas_probe_dir
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CheckBlasIntegerABI")
    set(_scalapack_openblas_probe_src
        "${_scalapack_openblas_probe_dir}/openblas_config_probe.c")
    file(MAKE_DIRECTORY "${_scalapack_openblas_probe_dir}")
    file(WRITE "${_scalapack_openblas_probe_src}"
"#include <stdio.h>
#include <stdlib.h>
#include <string.h>
extern const char *openblas_get_config(void);
int main(void) {
  const char *config = openblas_get_config();
  const int detected_bytes = (config != NULL && strstr(config, \"USE64BITINT\") != NULL) ? 8 : 4;
  if (config == NULL) {
    puts(\"OPENBLAS_CONFIG=<null>\");
  } else {
    printf(\"OPENBLAS_CONFIG=%s\\n\", config);
  }
  printf(\"OPENBLAS_DETECTED_INT_BYTES=%d\\n\", detected_bytes);
  return detected_bytes == ${expected_bytes} ? EXIT_SUCCESS : EXIT_FAILURE;
}
")

    try_run(_scalapack_openblas_probe_run_result
            _scalapack_openblas_probe_compile_result
            SOURCES "${_scalapack_openblas_probe_src}"
            NO_CACHE
            LINK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
            COMPILE_OUTPUT_VARIABLE _scalapack_openblas_probe_compile_output
            RUN_OUTPUT_VARIABLE _scalapack_openblas_probe_run_output)

    string(STRIP "${_scalapack_openblas_probe_compile_output}"
           _scalapack_openblas_probe_compile_output)
    string(STRIP "${_scalapack_openblas_probe_run_output}"
           _scalapack_openblas_probe_run_output)

    if(NOT _scalapack_openblas_probe_compile_result)
      message(FATAL_ERROR
        "Failed to compile the OpenBLAS integer ABI probe even though openblas_get_config() "
        "linked earlier. This usually indicates an inconsistent BLAS/LAPACK link line.\n"
        "Compile output:\n${_scalapack_openblas_probe_compile_output}")
    endif()

    if(_scalapack_openblas_probe_run_result MATCHES "^[0-9-]+$"
       AND _scalapack_openblas_probe_run_result EQUAL 0)
      if(NOT "${_scalapack_openblas_probe_run_output}" STREQUAL "")
        message(STATUS
          "OpenBLAS integer ABI matches SCALAPACK_EXPECT_BLAS_INT_BYTES=${expected_bytes}. "
          "${_scalapack_openblas_probe_run_output}")
      else()
        message(STATUS
          "OpenBLAS integer ABI matches SCALAPACK_EXPECT_BLAS_INT_BYTES=${expected_bytes}.")
      endif()
    elseif(_scalapack_openblas_probe_run_result MATCHES "^[0-9-]+$")
      message(FATAL_ERROR
        "Detected an OpenBLAS integer ABI mismatch. "
        "SCALAPACK_EXPECT_BLAS_INT_BYTES=${expected_bytes}, but the linked OpenBLAS probe "
        "reported a different integer width.\n"
        "Probe output:\n${_scalapack_openblas_probe_run_output}\n"
        "Point CMake at a matching BLAS/LAPACK, or change SCALAPACK_EXPECT_BLAS_INT_BYTES "
        "explicitly if you really intend to use a different ABI.")
    else()
      message(FATAL_ERROR
        "Failed to execute the OpenBLAS integer ABI probe. "
        "This is distinct from an ABI mismatch: the probe executable built successfully but "
        "could not be run. On clusters this often means a runtime loader issue such as a "
        "missing RPATH or `LD_LIBRARY_PATH` entry for the BLAS shared library.\n"
        "Run result: ${_scalapack_openblas_probe_run_result}\n"
        "Probe output:\n${_scalapack_openblas_probe_run_output}")
    endif()
  else()
    message(STATUS "OpenBLAS integer ABI check skipped: resolved BLAS does not expose openblas_get_config().")
  endif()

  set(CMAKE_REQUIRED_LIBRARIES "${_saved_required_libraries}")
  set(SCALAPACK_BLAS_SYMBOL_STYLE "${_scalapack_blas_symbol_style}" PARENT_SCOPE)
endfunction()
