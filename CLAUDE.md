# ScaLAPACK I8 Refactor

Starting with an initial commit (b935167ca4d244735abc04a3cd4f6d56699702a0) after forking
ScaLAPACK, we have been refactoring the codebase to be compatible with 64-bit integers
for MPI 4+. Planning and progress notes are in `docs/i8-refactor/`.

## Current state (Phase 8 complete)

The tridiagonal reduction driver cone is fully I8-native.  All narrowing is
encapsulated in dedicated _I8 wrapper files — zero inline INT() narrowing
remains in the four reduction drivers.

### I8 surface

- **Redistribution core:** PxGEMR2D_I8, PxTRMR2D_I8 (C, all types)
- **Validation layer:** CHK1MAT_I8, PCHK1MAT_I8, PCHK2MAT_I8, DESC_CONVERT_I8, GLOBCHK_I8
- **Tool routines:** NUMROC_I8, INDXL2G_I8, INDXG2L_I8, INDXG2P_I8, INFOG1L_I8, INFOG2L_I8, DESCINIT_I8, DESCSET_I8, PxELSET_I8, PxELGET_I8, NARROW_DESC8
- **Copy utilities:** xLAMOV_I8, PxLACP2_I8, PxLACPY_I8 (all 4 types)
- **Redistribution wrappers:** PxLAMR1D_I8 (all 4 types), PxLAMVE_I8 (D, S)
- **PBLAS I8 entry points:** 42 C wrappers via pblas_i8_utils.h
  - Level 1: PxAXPY, PxSCAL, PCSSCAL, PZDSCAL, PxNRM2, PxDOT/DOTC (18)
  - Level 2: PxGEMV, PxSYMV/HEMV (8)
  - Level 3: PxSYR2K/HER2K, PxSYRK/HERK, PxTRSM, PxGEMM (16)
- **Fortran I8 auxiliaries:** PxLARFG (4), PxLACGV (2), PxLATRD (4), PxSYTD2_I8 (2), PxHETD2_I8 (2)
- **Serial LAPACK I8 wrappers:** xSYTRD_I8 (2), xHETRD_I8 (2)
- **Tailored parallel I8 wrappers:** PxSYTTRD_I8 (2), PxHETTRD_I8 (2)
- **Banded/tridiagonal I8 wrappers:** 56 routines (DB, DT, GB, PB, PT families, all 4 types)
- **Eigenvalue solvers:** PxSYEV_I8 (2), PxHEEV_I8 (2)
- **Eigensolver support:** PxLASCL_I8 (4), PxLASET_I8 (4), PxORMTR_I8 (2), PxUNMTR_I8 (2)
- **Cholesky solvers:** PxPOTRF_I8 (4), PxPOTRS_I8 (4), PxPOSV_I8 (4)
- **LU solvers:** PxGETRF_I8 (4), PxGETRS_I8 (4), PxGESV_I8 (4)
- **LU support:** PxLASWP_I8 (4), PxLAPIV_I8 (4)
- **Reduction drivers:** PDSYNTRD_I8, PSSYNTRD_I8, PCHENTRD_I8, PZHENTRD_I8
  - All paths (blocked, serial, tailored): fully I8-native (zero inline narrowing)
  - Bit-identical to legacy counterparts

12 ctest targets, all passing on macOS arm64 and x86_64 Linux.

### Test coverage

- Driver tests exercise both UPLO='L' (serial/tailored path) and UPLO='U' (blocked path)
- PBLAS I8 kernel tests cover Level 1 (AXPY, SCAL, DOT, NRM2, DOTC), Level 2 (GEMV, SYMV, HEMV), Level 3 (SYR2K, HER2K)
- All comparisons are bit-identical against legacy routines

## Design rules

- DTYPE encodes layout (1=2D, 501=1D-H, 502=1D-V). Width lives in symbol names (_I8 suffix).
- All public _I8 integer arguments are INTEGER*8. INFO stays default INTEGER.
- BLACS context narrowed to INTEGER at call boundary.
- PBLAS I8 wrappers use pblas_i8_utils.h with SCALAPACK_FORTRAN_INT_BYTES for ILP64-safe range checks.
- Workspace reductions use DBLE/DGAMN2D (never REAL/SGAMN2D) for I8-safe global min.
