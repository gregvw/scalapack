# Phase 9-11 Notes: Banded Solvers, Eigensolvers, Dense Solvers, Large-N Closure (Complete)

## Overview

Phases 9-11 broadened the I8 surface from the tridiagonal reduction cone
to complete solver families, culminating in truly large-N capable dense
direct solvers.

## Phase 9: Banded/Tridiagonal I8 Wrappers

56 thin wrappers for the 1D-descriptor-heavy banded family:

- DB: PxDBTRF_I8, PxDBTRS_I8, PxDBTRSV_I8 (12 files)
- DT: PxDTTRF_I8, PxDTTRS_I8, PxDTTRSV_I8 (12 files)
- GB: PxGBTRF_I8, PxGBTRS_I8 (8 files)
- PB: PxPBTRF_I8, PxPBTRS_I8, PxPBTRSV_I8 (12 files)
- PT: PxPTTRF_I8, PxPTTRS_I8, PxPTTRSV_I8 (12 files)

All are thin wrappers that narrow I8 arguments to call legacy routines.
They do not support N > INTMAX.

## Phase 10: Eigensolver Cone

### Native I8 eigensolver drivers (4)

PxSYEV_I8 (D, S), PxHEEV_I8 (C, Z) — forked from legacy with all
internal calls replaced by I8 variants:

- PxSYNTRD_I8 / PxHENTRD_I8 (tridiagonal reduction)
- PxELGET_I8 (element extraction)
- PxLASCL_I8 (scaling)
- PxLASET_I8 (initialization)
- PxORMTR_I8 / PxUNMTR_I8 (transformation apply)
- PxGEMR2D_I8 (redistribution)
- DESCINIT_I8, NUMROC_I8, INDXG2P_I8 (tools)
- CHK1MAT_I8, PCHK1MAT_I8, PCHK2MAT_I8 (validation)

One remaining narrowing boundary: PxLANSY / PxLANHE (matrix norm).
INTMAX guard rejects N > INTMAX because of this.

### Support wrappers (12)

- PxLASCL_I8 (4) — matrix scaling
- PxLASET_I8 (4) — matrix initialization
- PxORMTR_I8 (2) / PxUNMTR_I8 (2) — transformation apply

### Bugs found and fixed during review

- PZHEEV_I8: INT(ABS(WORK(1))) → INT(ABS(WORK(1)),8) for workspace sizes
- PCHEEV_I8/PZHEEV_I8: missing RWORK/LRWORK args to PCHENTRD_I8/PZHENTRD_I8

## Phase 11: Dense Direct Solver Cones + Large-N Closure

### PBLAS Level 3 I8 (16 C wrappers, Phase 11a)

- PxTRSM_I8 (4) — triangular solve
- PxSYRK_I8 (2) / PxHERK_I8 (2) — rank-k updates
- PxGEMM_I8 (4) — general matrix multiply

### Cholesky cone (12 native Fortran drivers, Phase 11b)

- PxPOTRF_I8 (4) — blocked factorization
- PxPOTRS_I8 (4) — triangular solves
- PxPOSV_I8 (4) — composition

### LU cone (20 native + thin Fortran files, Phase 11c)

- PxLASWP_I8 (4) — row swaps (thin wrappers)
- PxLAPIV_I8 (4) — pivot apply (thin wrappers)
- PxGETRF_I8 (4) — blocked factorization
- PxGETRS_I8 (4) — solve with pivots
- PxGESV_I8 (4) — composition

### Large-N closure (20 new files, Phase 11d)

Removed the INTMAX entry guards from POTRF_I8 and GETRF_I8 by
introducing native I8 unblocked panel routines:

New PBLAS I8 wrappers (12 C files):
- PxAMAX_I8 (4), PxSWAP_I8 (4), PxGER/GERU_I8 (4)

Native I8 panel routines (8 Fortran files):
- PxPOTF2_I8 (4) — uses INFOG2L_I8 + local BLAS
- PxGETF2_I8 (4) — uses PxAMAX/SWAP/SCAL/GER_I8

Blocked driver updates: removed INTMAX guards, NARROW_DESC8 calls,
DESCA4 arrays. Dense solver cones now handle N > 2^31-1.

Bug fix: IPIV INTEGER→INTEGER*8 mismatch at PDSWAP_I8 call in
GETF2_I8 (all 4 types).

### Cross-cutting fixes applied during review

- PBLAS I8 early returns removed (all 42→54 C wrappers)
- PxPOTRS_I8 validation order: quick return moved after arg checking
- NARROW_DESC8: now aborts on overflow instead of silent truncation
- Thin wrapper quick returns removed (80 files) to preserve workspace query
- Edge-case regression tests: K=0 beta scaling, bad-UPLO on N=0

## Test coverage

12 I8 ctest targets:
- xi8tools, xi8tools_mpi
- xpblas_i8 (18 kernel tests: L1/L2/L3 + K=0 edge case)
- xdgemr_i8, xdtrmr_i8
- xdsyntrd_i8, xchentrd_i8
- xgesv_i8 (D, S, C, Z), xposv_i8 (D, S, C, Z + bad-UPLO edge)
- xsyev_i8 (D, S, C, Z)
- xdlamr1d_i8, xdlamve_i8

108 total tests on Linux, all passing.

## PBLAS I8 surface: 54 entry points

- Level 1: AXPY, SCAL, CSSCAL, ZDSCAL, NRM2, DOT/DOTC, AMAX, SWAP (26)
- Level 2: GEMV, SYMV/HEMV, GER/GERU (16)
- Level 3: SYR2K/HER2K, SYRK/HERK, TRSM, GEMM (16)
Note: some sub-families have fewer than 4 types (e.g., CSSCAL/ZDSCAL are 2)

## Current state

The dense direct-solve cones (LU, Cholesky) are fully large-N capable.
The eigensolver cone has one remaining internal narrowing boundary
(PxLANSY/PxLANHE). The banded/tridiagonal wrappers are thin wrappers
with N bounded by INTMAX.
