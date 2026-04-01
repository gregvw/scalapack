include_guard(GLOBAL)

macro(scalapack_configure_build_options)
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

  include(CheckBLACSCompilerFlags)
  CheckBLACSCompilerFlags()

  include(GNUInstallDirs)
  set(PKG_CONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig CACHE PATH "pkg-config install path")
endmacro()
