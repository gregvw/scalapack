# ScaLAPACK I8 Refactor

Starting with an initial commit (b935167ca4d244735abc04a3cd4f6d56699702a0) after forking
ScaLAPACK, we have been refactoring the codebase to be compatible with 64-bit integers
for MPI 4+. Planning and progress notes are in `docs/i8-refactor/`.

## Current state (Phase 4 complete)

The I8 surface includes:
- **Redistribution core:** PxGEMR2D_I8, PxTRMR2D_I8 (C, all types)
- **Validation layer:** CHK1MAT_I8, PCHK1MAT_I8, PCHK2MAT_I8, DESC_CONVERT_I8, GLOBCHK_I8
- **Tool routines:** NUMROC_I8, INDXL2G_I8, INDXG2L_I8, INDXG2P_I8, INFOG1L_I8, INFOG2L_I8, DESCINIT_I8, DESCSET_I8
- **Copy utilities:** xLAMOV_I8, PxLACP2_I8, PxLACPY_I8 (all 4 types)
- **Element setter:** PxELSET_I8 (all 4 types)
- **Redistribution wrappers:** PxLAMR1D_I8 (all 4 types), PxLAMVE_I8 (D, S)
- **Bridge drivers:** PDSYNTRD_I8, PSSYNTRD_I8, PCHENTRD_I8, PZHENTRD_I8
  - Fully I8-native internal arithmetic
  - Checked narrowing at PBLAS/LAPACK boundaries
  - Bit-identical to legacy counterparts

7 ctest targets, 103+ tests passing on x86_64 Linux and macOS arm64.

## What work remains (Phase 5)

See `docs/i8-refactor/PHASE5-PLAN.md` for the native I8 PBLAS slice plan.

Key decision pending: duplicate _I8 C files vs parameterize shared PBLAS cores.

First-wave PBLAS kernels (driven by PxLATRD dependency cone):
PxAXPY, PxSCAL, PxNRM2, PxDOT/DOTC, PxGEMV, PxSYMV/HEMV

## Design rules

- DTYPE encodes layout (1=2D, 501=1D-H, 502=1D-V). Width lives in symbol names (_I8 suffix).
- All public _I8 integer arguments are INTEGER*8. INFO stays default INTEGER.
- BLACS context narrowed to INTEGER at call boundary.
- Bridge drivers narrow at PBLAS/LAPACK boundaries with overflow checks + BLACS_ABORT.
- Workspace reductions use DGAMN2D/SGAMN2D (via DOUBLE/REAL) for I8-safe global min.
