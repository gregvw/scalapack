include_guard(GLOBAL)

macro(SCALAPACK_install_library lib)
  install(TARGETS ${lib} EXPORT scalapack-targets
    RUNTIME DESTINATION Testing
  )
endmacro()

function(invertBoolean varName varValue)
  if (${varValue})
    set(${varName} false PARENT_SCOPE)
  else ()
    set(${varName} true PARENT_SCOPE)
  endif ()
endfunction ()

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
