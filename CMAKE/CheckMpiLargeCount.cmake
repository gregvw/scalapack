include_guard(GLOBAL)

function(CheckMpiLargeCount output_var)
  set(options)
  set(oneValueArgs OUTPUT_LOG_VAR)
  set(multiValueArgs)
  cmake_parse_arguments(PARSE_ARGV 1 _scalapack_check_mpi
                        "${options}" "${oneValueArgs}" "${multiValueArgs}")

  if(NOT TARGET MPI::MPI_C)
    message(FATAL_ERROR "CheckMpiLargeCount requires the imported target MPI::MPI_C.")
  endif()

  try_compile(
    _scalapack_have_mpi_large_count_apis
    SOURCE_FROM_FILE scalapack_check_mpi_large_count.c
                     "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/check_mpi_large_count.c"
    LINK_LIBRARIES MPI::MPI_C
    LOG_DESCRIPTION "Checking MPI large-count _c APIs required by ScaLAPACK BLACS"
    NO_CACHE
    OUTPUT_VARIABLE _scalapack_check_mpi_log
  )

  set(${output_var} "${_scalapack_have_mpi_large_count_apis}" PARENT_SCOPE)
  if(DEFINED _scalapack_check_mpi_OUTPUT_LOG_VAR
     AND NOT "${_scalapack_check_mpi_OUTPUT_LOG_VAR}" STREQUAL "")
    set(${_scalapack_check_mpi_OUTPUT_LOG_VAR}
        "${_scalapack_check_mpi_log}" PARENT_SCOPE)
  endif()
endfunction()
