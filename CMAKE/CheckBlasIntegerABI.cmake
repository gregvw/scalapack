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

  check_c_source_compiles(
"extern const char *openblas_get_config(void);
int main(void) {
  return openblas_get_config() == 0;
}
" SCALAPACK_HAVE_OPENBLAS_GET_CONFIG)

  if(SCALAPACK_HAVE_OPENBLAS_GET_CONFIG)
    check_c_source_runs(
"#include <stdlib.h>
#include <string.h>
extern const char *openblas_get_config(void);
int main(void) {
  const char *config = openblas_get_config();
  const int detected_bytes = (config != NULL && strstr(config, \"USE64BITINT\") != NULL) ? 8 : 4;
  return detected_bytes == ${expected_bytes} ? EXIT_SUCCESS : EXIT_FAILURE;
}
" SCALAPACK_OPENBLAS_INT_ABI_MATCHES)

    if(NOT SCALAPACK_OPENBLAS_INT_ABI_MATCHES)
      message(FATAL_ERROR
        "Detected an OpenBLAS integer ABI mismatch. "
        "SCALAPACK_EXPECT_BLAS_INT_BYTES=${expected_bytes}, but the linked OpenBLAS reports "
        "${expected_bytes}-byte integers are not in use. Point CMake at a matching BLAS "
        "or change SCALAPACK_EXPECT_BLAS_INT_BYTES explicitly.")
    endif()

    message(STATUS "OpenBLAS integer ABI matches SCALAPACK_EXPECT_BLAS_INT_BYTES=${expected_bytes}.")
  else()
    message(STATUS "OpenBLAS integer ABI check skipped: resolved BLAS does not expose openblas_get_config().")
  endif()

  set(CMAKE_REQUIRED_LIBRARIES "${_saved_required_libraries}")
  set(SCALAPACK_BLAS_SYMBOL_STYLE "${_scalapack_blas_symbol_style}" PARENT_SCOPE)
endfunction()
