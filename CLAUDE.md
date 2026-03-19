# ScaLAPACK I8 Refactor

Starting with an initial commit (b935167ca4d244735abc04a3cd4f6d56699702a0) after forking
ScaLAPACK, we have been refactoring the codebase to be compatible with 64-bit integers
for MPI 4+. Planning and progress notes are in `docs/i8-refactor/`.

## Current state (Phase 11 complete — large-N solver closure)

Dense direct-solve cones (LU, Cholesky) are fully large-N capable with no
INTMAX entry guards.  All narrowing in the LU/Cholesky paths is encapsulated
in native I8 panel routines (PxGETF2_I8, PxPOTF2_I8).

### I8 surface

- **Redistribution core:** PxGEMR2D_I8, PxTRMR2D_I8 (C, all types)
- **Validation layer:** CHK1MAT_I8, PCHK1MAT_I8, PCHK2MAT_I8, DESC_CONVERT_I8, GLOBCHK_I8
- **Tool routines:** NUMROC_I8, INDXL2G_I8, INDXG2L_I8, INDXG2P_I8, INFOG1L_I8, INFOG2L_I8, DESCINIT_I8, DESCSET_I8, PxELSET_I8, PxELGET_I8, NARROW_DESC8
- **Copy utilities:** xLAMOV_I8, PxLACP2_I8, PxLACPY_I8 (all 4 types)
- **Redistribution wrappers:** PxLAMR1D_I8 (all 4 types), PxLAMVE_I8 (D, S)
- **PBLAS I8 entry points:** 54 C wrappers via pblas_i8_utils.h
  - Level 1: PxAXPY, PxSCAL, PCSSCAL, PZDSCAL, PxNRM2, PxDOT/DOTC, PxAMAX, PxSWAP (26)
  - Level 2: PxGEMV, PxSYMV/HEMV, PxGER/GERU (16)
  - Level 3: PxSYR2K/HER2K, PxSYRK/HERK, PxTRSM, PxGEMM (16)
- **Tridiagonal reduction cone:** PxSYNTRD_I8, PxHENTRD_I8 (4 native drivers), PxLATRD_I8 (4), PxLARFG_I8 (4), PxLACGV_I8 (2), PxSYTD2_I8 (2), PxHETD2_I8 (2), xSYTRD_I8 (2), xHETRD_I8 (2), PxSYTTRD_I8 (2), PxHETTRD_I8 (2)
- **Banded/tridiagonal I8 wrappers:** 56 routines (DB, DT, GB, PB, PT families, all 4 types)
- **Eigenvalue solvers:** PxSYEV_I8 (2), PxHEEV_I8 (2) + support: PxLASCL_I8 (4), PxLASET_I8 (4), PxORMTR_I8 (2), PxUNMTR_I8 (2)
- **Cholesky solvers:** PxPOTRF_I8 (4), PxPOTRS_I8 (4), PxPOSV_I8 (4)
- **LU solvers:** PxGETRF_I8 (4), PxGETRS_I8 (4), PxGESV_I8 (4) + support: PxLASWP_I8 (4), PxLAPIV_I8 (4)
- **Unblocked panel I8:** PxPOTF2_I8 (4), PxGETF2_I8 (4)

13 ctest targets, all passing on macOS arm64 and x86_64 Linux.

### Large-N status

The dense solver cones (Cholesky, LU) are fully large-N capable: PxGETF2_I8 and
PxPOTF2_I8 are native I8 panel routines, and the blocked drivers have no INTMAX
entry guards.  54 PBLAS I8 entry points support the full call tree.

The eigenvalue cone still has an INTMAX guard because PxLANSY/PxLANHE (matrix
norm) lacks an _I8 version.

The banded/tridiagonal wrappers are thin wrappers that narrow all arguments —
they do not support N > INTMAX.

### Test coverage

- **Solver tests:** xgesv_i8 (D, S, C, Z), xposv_i8 (D, S, C, Z + bad-UPLO N=0 edge case)
- **Eigensolver tests:** xsyev_i8 (D, S, C, Z)
- **Reduction tests:** xdsyntrd_i8 (D, S, UPLO=L+U), xchentrd_i8 (C, Z, UPLO=L+U)
- **PBLAS kernel tests:** xpblas_i8 (18 tests: L1/L2/L3 + K=0 beta edge case, D and C types)
- **Tools tests:** xi8tools (standalone), xi8tools_mpi (DESC_CONVERT 2D→1D, PCHK1MAT)
- **Large-index test:** xlargeidx_i8 (NUMROC/DESCINIT/INDXL2G with N=3B, descriptor > INTMAX)
- **Redistribution tests:** xdgemr_i8, xdtrmr_i8, xdlamr1d_i8, xdlamve_i8
- All comparisons are bit-identical against legacy routines

## Design rules

- DTYPE encodes layout (1=2D, 501=1D-H, 502=1D-V). Width lives in symbol names (_I8 suffix).
- All public _I8 integer arguments are INTEGER*8. INFO stays default INTEGER.
- BLACS context narrowed to INTEGER at call boundary.
- PBLAS I8 wrappers use pblas_i8_utils.h. No early returns — legacy handles all edge cases.
- NARROW_DESC8 aborts on overflow (fail-fast, not clamp).
- Workspace queries pass through to legacy (no premature quick returns in thin wrappers).
- Workspace reductions use DBLE/DGAMN2D (never REAL/SGAMN2D) for I8-safe global min.
