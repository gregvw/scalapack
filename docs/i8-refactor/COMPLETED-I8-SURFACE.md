# Completed I8 Surface

Inventory of all public `_I8` routines implemented so far.

## Foundation layer

- Redistribution core: `PxGEMR2D_I8`, `PxTRMR2D_I8`
- Validation/tools: `CHK1MAT_I8`, `PCHK1MAT_I8`, `PCHK2MAT_I8`,
  `DESC_CONVERT_I8`, `GLOBCHK_I8`, `DESCINIT_I8`, `DESCSET_I8`,
  `NUMROC_I8`, `INDX*_I8`, `INFOG*L_I8`, `NARROW_DESC8`
- Copy/move: `xLAMOV_I8`, `PxLACP2_I8`, `PxLACPY_I8` (all 4 types)
- Redistribution wrappers: `PxLAMR1D_I8` (4), `PxLAMVE_I8` (D, S)
- Element access: `PxELSET_I8`, `PxELGET_I8` (all 4 types)

## PBLAS I8 entry points (58 C wrappers)

- Level 1: `PxAXPY`, `PxSCAL`, `PCSSCAL`, `PZDSCAL`, `PxNRM2`,
  `PxDOT/DOTC`, `PxAMAX`, `PxSWAP`, `PxCOPY` (30)
- Level 2: `PxGEMV`, `PxSYMV/HEMV`, `PxGER/GERU` (16)
- Level 3: `PxSYR2K/HER2K`, `PxSYRK/HERK`, `PxTRSM`, `PxGEMM` (16)

## Tridiagonal reduction cone

- Reduction drivers: `PxSYNTRD_I8`, `PxHENTRD_I8` (4 native drivers)
- Support: `PxLATRD_I8` (4), `PxLARFG_I8` (4), `PxLACGV_I8` (2)
- Unblocked: `PxSYTD2_I8` (2), `PxHETD2_I8` (2)
- Serial LAPACK: `xSYTRD_I8` (2), `xHETRD_I8` (2)
- Tailored parallel: `PxSYTTRD_I8` (2), `PxHETTRD_I8` (2)

## Banded/tridiagonal solvers (56 thin wrappers)

- DB family: `PxDBTRF_I8`, `PxDBTRS_I8`, `PxDBTRSV_I8` (12)
- DT family: `PxDTTRF_I8`, `PxDTTRS_I8`, `PxDTTRSV_I8` (12)
- GB family: `PxGBTRF_I8`, `PxGBTRS_I8` (8)
- PB family: `PxPBTRF_I8`, `PxPBTRS_I8`, `PxPBTRSV_I8` (12)
- PT family: `PxPTTRF_I8`, `PxPTTRS_I8`, `PxPTTRSV_I8` (12)

## Eigenvalue solvers (4 native drivers + 12 support)

- Drivers: `PxSYEV_I8` (2), `PxHEEV_I8` (2)
- Support: `PxLASCL_I8` (4), `PxLASET_I8` (4), `PxORMTR_I8` (2),
  `PxUNMTR_I8` (2)

## Dense direct solvers

- Cholesky: `PxPOTRF_I8` (4), `PxPOTRS_I8` (4), `PxPOSV_I8` (4)
- LU: `PxGETRF_I8` (4), `PxGETRS_I8` (4), `PxGESV_I8` (4)
- LU support: `PxLASWP_I8` (4), `PxLAPIV_I8` (4)
- Unblocked panels: `PxPOTF2_I8` (4), `PxGETF2_I8` (4)
- Matrix inverse: `PxGETRI_I8` (4)

## Matrix norms (4 thin wrappers)

- `PxLANGE_I8` (4) — general matrix norm (M, 1, I, F)

## QR factorization cone (16 thin wrappers)

- QR factorization: `PxGEQRF_I8` (4)
- Q generation: `PxORGQR_I8` (2), `PxUNGQR_I8` (2)
- Q multiply: `PxORMQR_I8` (2), `PxUNMQR_I8` (2)

## Large-N status

Dense solver cones (LU, Cholesky) are fully large-N capable: `PxGETF2_I8`
and `PxPOTF2_I8` are native I8 panel routines with no INTMAX entry guards.

The eigenvalue cone has an INTMAX guard because `PxLANSY`/`PxLANHE`
(matrix norm) lacks an `_I8` version.

The banded/tridiagonal wrappers are thin wrappers that narrow all
arguments — they do not support N > INTMAX.

## Notes

- Internal `PB_*` PBLAS helpers are intentionally excluded.
  This manifest covers public or phase-planning surface only.
- This file uses exact file-stem matching.  If a future `_I8` routine
  deliberately reuses a legacy stem through internal dispatch rather than
  adding a new `*_i8` file, update this manifest manually.
