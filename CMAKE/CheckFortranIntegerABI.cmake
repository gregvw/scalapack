include(CheckFortranSourceRuns)

function(CheckFortranIntegerABI out_var)
  if(CMAKE_CROSSCOMPILING)
    if(DEFINED SCALAPACK_FORTRAN_INT_BYTES AND
       NOT "${SCALAPACK_FORTRAN_INT_BYTES}" STREQUAL "")
      if(NOT SCALAPACK_FORTRAN_INT_BYTES MATCHES "^(4|8)$")
        message(FATAL_ERROR
          "SCALAPACK_FORTRAN_INT_BYTES must be 4 or 8 when set explicitly.")
      endif()
      set(${out_var} "${SCALAPACK_FORTRAN_INT_BYTES}" PARENT_SCOPE)
      return()
    endif()
    message(FATAL_ERROR
      "Cannot auto-detect the default Fortran INTEGER width while cross-compiling. "
      "Set SCALAPACK_FORTRAN_INT_BYTES to 4 or 8.")
  endif()

  unset(SCALAPACK_FORTRAN_INT_IS_4 CACHE)
  unset(SCALAPACK_FORTRAN_INT_IS_8 CACHE)

  check_fortran_source_runs(
"program main
  integer :: i
  if (storage_size(i) / 8 /= 4) stop 1
end program main
" SCALAPACK_FORTRAN_INT_IS_4)

  check_fortran_source_runs(
"program main
  integer :: i
  if (storage_size(i) / 8 /= 8) stop 1
end program main
" SCALAPACK_FORTRAN_INT_IS_8)

  if(SCALAPACK_FORTRAN_INT_IS_4)
    set(${out_var} "4" PARENT_SCOPE)
  elseif(SCALAPACK_FORTRAN_INT_IS_8)
    set(${out_var} "8" PARENT_SCOPE)
  else()
    message(FATAL_ERROR
      "Unable to determine the default Fortran INTEGER width. "
      "Set SCALAPACK_FORTRAN_INT_BYTES to 4 or 8 explicitly.")
  endif()
endfunction()
