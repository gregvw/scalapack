# Phase 4 Plan

  1. Fix the descriptor contract first.
      - Update SL_i8_params.inc:1 so I8 descriptors use the same semantic DTYPE values as legacy descriptors.
      - Update descinit_i8.f:1 and descset_i8.f:1 to emit DTYPE=1 for 2D block-cyclic descriptors, not 501.
      - Update redist_core.h:1 so unpack_desc_i8 accepts semantic descriptor tags, not a width-encoded tag.
      - Update the phase notes in PHASE3-NOTES.md:1 so the design record matches the corrected ABI strategy.
  2. Build the missing I8 validation layer in TOOLS.
      - Add chk1mat_i8.f beside chk1mat.f:1.
      - Add pchkxmat_i8.f beside pchkxmat.f:1, containing PCHK1MAT_I8 and PCHK2MAT_I8.
      - Add desc_convert_i8.f beside desc_convert.f:1.
      - Register them in TOOLS/CMakeLists.txt:3.
  3. Solve the global-check transport cleanly.
      - Current PCHK1MAT packs values into default INTEGER work arrays in pchkxmat.f:1, which is wrong for I8-scale values.
      - Best option: add a small C helper that does an MPI_Allreduce over int64_t (value, position) pairs and call it from PCHK1MAT_I8.
      - Alternative: add a true Fortran GLOBCHK_I8, but that is more awkward and lower leverage.
  4. Add automated tests before migrating real routines.
      - Wire test_i8_tools.f:1 into CMake as a real test target.
      - Extend it with negative tests:
          - invalid DTYPE
          - invalid LLD
          - out-of-range RSRC/CSRC
          - descriptor conversion cases for 2D-to-1D I8
      - Keep the existing redistribution tests in REDIST/TESTING/CMakeLists.txt:30, but update them to use the corrected descriptor tag.
  5. Migrate the first Fortran _I8 routines in a narrow band.
      - Start with routines that already call redistribution and have relatively contained argument checking.
      - Good discovery set is the 26 SRC files that reference PxGEMR2D/PxTRMR2D; examples include pdgesvd.f:1, pcgesvd.f:1, pzgesvd.f:1, pdsyev.f:1,
        pssyev.f:1, pcheev.f:1, and pzheev.f:1.
      - I would not start with gesvd or eigensolver drivers despite the redistribution dependency; they are high-risk. Prefer a smaller
        redistribution-adjacent routine first if you can find one with limited workspace and simpler validation.
  6. Defer the 1D-descriptor-heavy banded/tridiagonal family until later.
      - There are about 56 SRC files using DESC_CONVERT or legacy 501/502 1D logic.
      - Those routines are the most sensitive to descriptor-tag confusion, so they become much safer once DESC_CONVERT_I8 exists and DTYPE semantics
        are normalized.
      - This family should be a separate phase after the first successful _I8 routine wave.

#  Suggested implementation order

  1. Descriptor tag correction.
  2. CHK1MAT_I8.
  3. DESC_CONVERT_I8.
  4. PCHK1MAT_I8 / PCHK2MAT_I8 with a 64-bit global-check path.
  5. CMake + test_i8_tools.
  6. One small Fortran _I8 routine end-to-end.
  7. Broader routine rollout.
  8. Banded/tridiagonal 1D-descriptor family.

#  Exit criteria for Phase 4

  - I8 2D descriptors use semantic DTYPE=1.
  - CHK1MAT_I8, PCHK1MAT_I8, PCHK2MAT_I8, and DESC_CONVERT_I8 exist and are tested.
  - test_i8_tools runs under ctest.
  - At least one non-redistribution Fortran _I8 entry point is implemented and tested end-to-end.
  - No code path infers descriptor width from DTYPE.
