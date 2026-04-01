set(FTestObj
   blacstest.f btprim.f tools.f)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set_source_files_properties(blacstest.f PROPERTIES COMPILE_FLAGS "-std=legacy")
endif()

# The upstream BLACS test harness (blacstest.f) communicates process
# coordinates via MPI_INTEGER, which is 4 bytes regardless of the
# Fortran default INTEGER width.  In ILP64 builds the default INTEGER
# is 8 bytes, so the communicated values are corrupted.  Until the
# test harness is ported to ILP64, skip these tests.
if(NOT SCALAPACK_ENABLE_ILP64)

add_executable(xFbtest ${FTestObj})
target_link_libraries(xFbtest scalapack MPI::MPI_Fortran)

set(CTestObj
   Cbt.c)

set_property(
   SOURCE Cbt.c
   APPEND PROPERTY COMPILE_DEFINITIONS BTCINTFACE
   )

add_executable(xCbtest ${CTestObj} ${FTestObj})
target_link_libraries(xCbtest scalapack MPI::MPI_Fortran)

file(COPY bsbr.dat DESTINATION ${SCALAPACK_BINARY_DIR}/BLACS/TESTING)
file(COPY bt.dat DESTINATION ${SCALAPACK_BINARY_DIR}/BLACS/TESTING)
file(COPY comb.dat DESTINATION ${SCALAPACK_BINARY_DIR}/BLACS/TESTING)
file(COPY sdrv.dat DESTINATION ${SCALAPACK_BINARY_DIR}/BLACS/TESTING)

# The BLACS tests exit abnormally by design, so a direct add_test
# would always report failure.  The runtest.cmake wrapper checks
# that the output contains "DONE BLACS_GRIDEXIT" instead.
add_test(xCbtest
     ${CMAKE_COMMAND}
    -DMPIEXEC=${MPIEXEC}
    -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
    -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
    -DTEST_PROG=xCbtest
    -DOUTPUTDIR=${SCALAPACK_BINARY_DIR}/BLACS/TESTING
    -DRUNTIMEDIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake
    )

add_test(xFbtest
     ${CMAKE_COMMAND}
    -DMPIEXEC=${MPIEXEC}
    -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
    -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
    -DTEST_PROG=xFbtest
    -DOUTPUTDIR=${SCALAPACK_BINARY_DIR}/BLACS/TESTING
    -DRUNTIMEDIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake
    )

endif() # NOT SCALAPACK_ENABLE_ILP64
