include_guard(GLOBAL)

macro(scalapack_configure_lapack_blas)
  option(USE_OPTIMIZED_LAPACK_BLAS "Whether or not to search for optimized LAPACK and BLAS libraries on your machine (if not found, Reference LAPACK and BLAS will be downloaded and installed)" ON)
  option(SCALAPACK_PREFER_OPENBLAS_ON_APPLE_SILICON
         "Prefer OpenBLAS over Accelerate on macOS arm64 unless the user explicitly selects another BLAS/LAPACK backend"
         ON)

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

  if(LAPACK_LIBRARIES)
    message(STATUS "--> LAPACK supplied by user is ${LAPACK_LIBRARIES}.")
    set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
    check_fortran_function_exists("dgesv" LAPACK_FOUND)
    unset(CMAKE_REQUIRED_LIBRARIES)
    message(STATUS "--> LAPACK routine dgesv is found: ${LAPACK_FOUND}.")
  endif()

  if(LAPACK_FOUND)
    message(STATUS "--> LAPACK supplied by user is WORKING, will use ${LAPACK_LIBRARIES}.")
  else()
    if(USE_OPTIMIZED_LAPACK_BLAS)
      message(STATUS "--> Searching for optimized LAPACK and BLAS libraries on your machine.")
      find_package(LAPACK)
      if(_scalapack_force_openblas
         AND NOT LAPACK_FOUND
         AND EXISTS "${_scalapack_openblas_direct_lib}")
         set(LAPACK_LIBRARIES "${_scalapack_openblas_direct_lib}")
         set(BLAS_LIBRARIES "${_scalapack_openblas_direct_lib}")
         set(CMAKE_REQUIRED_LIBRARIES "${LAPACK_LIBRARIES}")
         check_fortran_function_exists("dgesv" LAPACK_FOUND)
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
    endif()
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
         set(LAPACK_LIBRARIES ${SCALAPACK_BINARY_DIR}/lib/liblapack.a CACHE STRING "LAPACK library" FORCE)
         set(BLAS_LIBRARIES ${SCALAPACK_BINARY_DIR}/lib/libblas.a CACHE STRING "BLAS library" FORCE)
      else()
         set(LAPACK_LIBRARIES ${SCALAPACK_BINARY_DIR}/lib/liblapack.lib CACHE STRING "LAPACK library" FORCE)
         set(BLAS_LIBRARIES ${SCALAPACK_BINARY_DIR}/lib/libblas.lib CACHE STRING "BLAS library" FORCE)
      endif ()
    endif()
  endif()

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
endmacro()
