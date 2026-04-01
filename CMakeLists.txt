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
include(CheckMpiLargeCount)
include(ScalapackHelpers)
include(ScalapackMPI)
include(ScalapackBuildOptions)
include(ScalapackIntegerAbi)
include(ScalapackLapackBlas)

if (UNIX)
   if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
      set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fltconsistency -fp_port" )
   endif ()
endif ()

scalapack_configure_mpi()
message(STATUS "ILP64 support: probing MPI large-count _c APIs required by BLACS")

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

scalapack_configure_build_options()

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

scalapack_configure_integer_abi()
scalapack_configure_lapack_blas()

MESSAGE(STATUS "=========")

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
