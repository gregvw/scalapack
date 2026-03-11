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
endfunction()
