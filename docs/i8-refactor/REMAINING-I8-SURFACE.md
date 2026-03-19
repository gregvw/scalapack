# Remaining I8 Surface

Snapshot of the public routine surface, organized by I8 completion status.

This manifest is meant for planning.  It focuses on public or
planning-relevant routine families rather than internal C helpers.

## Already complete

### Foundation layer

- Redistribution core: `PxGEMR2D_I8`, `PxTRMR2D_I8`
- Validation/tools: `CHK1MAT_I8`, `PCHK1MAT_I8`, `PCHK2MAT_I8`,
  `DESC_CONVERT_I8`, `GLOBCHK_I8`, `DESCINIT_I8`, `DESCSET_I8`,
  `NUMROC_I8`, `INDX*_I8`, `INFOG*L_I8`, `NARROW_DESC8`
- Copy/move: `xLAMOV_I8`, `PxLACP2_I8`, `PxLACPY_I8` (all 4 types)
- Redistribution wrappers: `PxLAMR1D_I8` (4), `PxLAMVE_I8` (D, S)
- Element access: `PxELSET_I8`, `PxELGET_I8` (all 4 types)

### PBLAS I8 entry points (54 C wrappers)

- Level 1: `PxAXPY`, `PxSCAL`, `PCSSCAL`, `PZDSCAL`, `PxNRM2`,
  `PxDOT/DOTC`, `PxAMAX`, `PxSWAP` (26)
- Level 2: `PxGEMV`, `PxSYMV/HEMV`, `PxGER/GERU` (16)
- Level 3: `PxSYR2K/HER2K`, `PxSYRK/HERK`, `PxTRSM`, `PxGEMM` (16)

### Tridiagonal reduction cone

- Reduction drivers: `PxSYNTRD_I8`, `PxHENTRD_I8` (4 native drivers)
- Support: `PxLATRD_I8` (4), `PxLARFG_I8` (4), `PxLACGV_I8` (2)
- Unblocked: `PxSYTD2_I8` (2), `PxHETD2_I8` (2)
- Serial LAPACK: `xSYTRD_I8` (2), `xHETRD_I8` (2)
- Tailored parallel: `PxSYTTRD_I8` (2), `PxHETTRD_I8` (2)

### Banded/tridiagonal solvers (56 thin wrappers)

- DB family: `PxDBTRF_I8`, `PxDBTRS_I8`, `PxDBTRSV_I8` (12)
- DT family: `PxDTTRF_I8`, `PxDTTRS_I8`, `PxDTTRSV_I8` (12)
- GB family: `PxGBTRF_I8`, `PxGBTRS_I8` (8)
- PB family: `PxPBTRF_I8`, `PxPBTRS_I8`, `PxPBTRSV_I8` (12)
- PT family: `PxPTTRF_I8`, `PxPTTRS_I8`, `PxPTTRSV_I8` (12)

### Eigenvalue solvers (4 native drivers + 12 support)

- Drivers: `PxSYEV_I8` (2), `PxHEEV_I8` (2)
- Support: `PxLASCL_I8` (4), `PxLASET_I8` (4), `PxORMTR_I8` (2),
  `PxUNMTR_I8` (2)

### Dense direct solvers

- Cholesky: `PxPOTRF_I8` (4), `PxPOTRS_I8` (4), `PxPOSV_I8` (4)
- LU: `PxGETRF_I8` (4), `PxGETRS_I8` (4), `PxGESV_I8` (4)
- LU support: `PxLASWP_I8` (4), `PxLAPIV_I8` (4)
- Unblocked panels: `PxPOTF2_I8` (4), `PxGETF2_I8` (4)

Dense solver cones are fully large-N capable: PxGETF2_I8 and
PxPOTF2_I8 are native I8 panel routines with no INTMAX entry guards.

## Remaining: SVD cone

These require `PxGEBRD_I8` and `PxORMBR_I8` as foundational pieces
before the driver can be made native I8.

- `pcgesvd`, `pdgesvd`, `psgesvd`, `pzgesvd`

## Remaining: Schur / eigenvalue variants

Deferred — recursive/index-heavy families requiring deep refactoring.

- `pdlaqr0`, `pdlaqr2`, `pdlaqr3`, `pdlaqr4`
- `pslaqr0`, `pslaqr2`, `pslaqr3`, `pslaqr4`
- `pdsyevd`, `pssyevd`, `pcheevd`, `pzheevd`
- `pdsyevx`, `pssyevx`, `pcheevx`, `pzheevx`
- `pdsyevr`, `pssyevr`, `pcheevr`, `pzheevr`
- `pdsygvx`, `pssygvx`, `pchegvx`, `pzhegvx`

## Remaining: public PBLAS families without `_I8`

Current missing public PBLAS surface: 69 entry points.

### Level 3 / matrix-matrix

- `pctrmm`, `pdtrmm`, `pstrmm`, `pztrmm`
- `pcsymm`, `pdsymm`, `pssymm`, `pzsymm`
- `pcsyr2k`, `pzsyr2k`
- `pchemm`, `pzhemm`

### Level 2 / matrix-vector and rank updates

- `pcagemv`, `pdagemv`, `psagemv`, `pzagemv`
- `pcahemv`, `pzahemv`
- `pdasymv`, `psasymv`
- `pcatrmv`, `pdatrmv`, `psatrmv`, `pzatrmv`
- `pcgerc`, `pzgerc`
- `pcher`, `pzher`
- `pcher2`, `pzher2`
- `pdsyr`, `pssyr`
- `pdsyr2`, `pssyr2`
- `pctranc`, `pztranc`
- `pdtran`, `pstran`
- `pctranu`, `pztranu`
- `pctradd`, `pdtradd`, `pstradd`, `pztradd`

### Level 1 / vector

- `pdasum`, `psasum`
- `pscasum`, `pdzasum`
- `pccopy`, `pdcopy`, `picopy`, `pscopy`, `pzcopy`
- `pcdotu`, `pzdotu`
- `pctrmv`, `pdtrmv`, `pstrmv`, `pztrmv`
- `pctrsv`, `pdtrsv`, `pstrsv`, `pztrsv`
- `pcgeadd`, `pdgeadd`, `psgeadd`, `pzgeadd`

## Remaining: TOOLS surface without `_I8`

Mostly utility routines, not blockers for current solver cones.

### Scalar / local helpers

- `iceil`, `ilacpy`, `ilcm`, `npreroc`
- `dsnrm2`, `dscnrm2`, `dsasum`, `dscasum`
- `ssdot`, `dddot`, `ccdotu`, `ccdotc`, `zzdotu`, `zzdotc`
- `slatcpy`, `dlatcpy`, `clatcpy`, `zlatcpy`
- `smatadd`, `dmatadd`, `cmatadd`, `zmatadd`
- `sltimer`

### Distributed utility families

- `PxCHEKPAD`, `PxFILLPAD`, `PxLAPRNT` (4 types each)
- `PxLAREAD`, `PxLAWRITE` (4 types each)
- `PxMATADD`, `PxCOL2ROW`, `PxROW2COL` (4 types each)
- `PxTREECOMB` (4 types), `PxELSET2` (4 types)

### Integer tools / reshape

- `PICHEKPAD`, `PIFILLPAD`, `PICOL2ROW`, `PIROW2COL`
- `PILAPRNT`, `PITREECOMB`, `PIELGET`, `PIELSET`, `PIELSET2`
- `reshape`, `SL_gridreshape`, `SL_init`

## Notes

- This file uses exact file-stem matching.  If a future `_I8` routine
  deliberately reuses a legacy stem through internal dispatch rather than
  adding a new `*_i8` file, update this manifest manually.
- Internal `PB_*` PBLAS helpers are intentionally excluded from the PBLAS
  list above.  This manifest is about public or phase-planning surface.
- `PxGETF2_I8` and `PxPOTF2_I8` are now native I8 — the dense solver
  cones have no remaining internal narrowing boundaries.
- `PxLANSY` / `PxLANHE` remain internal narrowing boundaries in the
  eigensolver cone only.
