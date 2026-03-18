# ScaLAPACK I8 Refactor

Starting with an initial commit (b935167ca4d244735abc04a3cd4f6d56699702a0) after forking
ScaLAPACK, we have been refactoring the codebase to be compatible with 64-bit integers
for MPI 4+. Planning and progress notes are in `docs/i8-refactor/`.

## Current state (Phase 6 complete)

The blocked reduction loop is fully I8-native across all four driver families.

### I8 surface

- **Redistribution core:** PxGEMR2D_I8, PxTRMR2D_I8 (C, all types)
- **Validation layer:** CHK1MAT_I8, PCHK1MAT_I8, PCHK2MAT_I8, DESC_CONVERT_I8, GLOBCHK_I8
- **Tool routines:** NUMROC_I8, INDXL2G_I8, INDXG2L_I8, INDXG2P_I8, INFOG1L_I8, INFOG2L_I8, DESCINIT_I8, DESCSET_I8, PxELSET_I8, PxELGET_I8, NARROW_DESC8
- **Copy utilities:** xLAMOV_I8, PxLACP2_I8, PxLACPY_I8 (all 4 types)
- **Redistribution wrappers:** PxLAMR1D_I8 (all 4 types), PxLAMVE_I8 (D, S)
- **PBLAS I8 entry points:** 30 C wrappers via pblas_i8_utils.h
  - Level 1: PxAXPY, PxSCAL, PCSSCAL, PZDSCAL, PxNRM2, PxDOT/DOTC (18)
  - Level 2: PxGEMV, PxSYMV/HEMV (8)
  - Level 3: PxSYR2K/HER2K (4)
- **Fortran I8 auxiliaries:** PxLARFG (4), PxLACGV (2), PxLATRD (4)
- **Reduction drivers:** PDSYNTRD_I8, PSSYNTRD_I8, PCHENTRD_I8, PZHENTRD_I8
  - Blocked loop: fully I8-native (zero narrowing)
  - Bit-identical to legacy counterparts

8 ctest targets, 104 tests passing on x86_64 Linux, all passing on macOS arm64.
Zero ASan/UBSan errors in I8 code (clang-20 sanitizer build).

### Remaining narrowed boundaries

These execute once per reduction call, not per iteration:

| Call | When | Frequency |
|------|------|-----------|
| PxSYTD2 / PxHETD2 | Last block of blocked path | Once per reduction |
| xSYTRD / xHETRD | Serial LAPACK path | Once, small-N only |
| PxSYTTRD / PxHETTRD | Tailored parallel path | Once, large-workspace only |

## Design rules

- DTYPE encodes layout (1=2D, 501=1D-H, 502=1D-V). Width lives in symbol names (_I8 suffix).
- All public _I8 integer arguments are INTEGER*8. INFO stays default INTEGER.
- BLACS context narrowed to INTEGER at call boundary.
- PBLAS I8 wrappers use pblas_i8_utils.h with SCALAPACK_FORTRAN_INT_BYTES for ILP64-safe range checks.
- Workspace reductions use DBLE/DGAMN2D (never REAL/SGAMN2D) for I8-safe global min.
