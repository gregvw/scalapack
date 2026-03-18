# Remaining I8 Surface

Snapshot of the remaining routine surface without a stem-matched `_I8`
entry point after Phase 8.

This manifest is meant for planning.  It focuses on public or
planning-relevant routine families rather than internal C helpers.

## Already complete

The tridiagonal reduction cone is complete:

- Redistribution core: `PxGEMR2D_I8`, `PxTRMR2D_I8`
- Validation/tools layer used by the cone: `CHK1MAT_I8`,
  `PCHK1MAT_I8`, `PCHK2MAT_I8`, `DESC_CONVERT_I8`, `GLOBCHK_I8`,
  `DESCINIT_I8`, `DESCSET_I8`, `NUMROC_I8`, `INDX*_*I8`, `INFOG*L_I8`
- Reduction-driver cone: `PDSYNTRD_I8`, `PSSYNTRD_I8`, `PCHENTRD_I8`,
  `PZHENTRD_I8`, `PxLATRD_I8`, `PxSYTD2_I8`, `PxHETD2_I8`,
  `xSYTRD_I8`, `xHETRD_I8`, `PxSYTTRD_I8`, `PxHETTRD_I8`
- Native PBLAS slice used by the cone:
  `PxAXPY_I8`, `PxSCAL_I8`, `PxNRM2_I8`, `PxDOT/DOTC_I8`,
  `PxGEMV_I8`, `PxSYMV/HEMV_I8`, `PxSYR2K/HER2K_I8`

## Priority 1: 1D-descriptor-heavy banded / tridiagonal family

These 56 routines still depend on `DESC_CONVERT` and are the most
natural next phase if broader `_I8` coverage is desired.

- `pcdbtrf`
- `pcdbtrs`
- `pcdbtrsv`
- `pcdttrf`
- `pcdttrs`
- `pcdttrsv`
- `pcgbtrf`
- `pcgbtrs`
- `pcpbtrf`
- `pcpbtrs`
- `pcpbtrsv`
- `pcpttrf`
- `pcpttrs`
- `pcpttrsv`
- `pddbtrf`
- `pddbtrs`
- `pddbtrsv`
- `pddttrf`
- `pddttrs`
- `pddttrsv`
- `pdgbtrf`
- `pdgbtrs`
- `pdpbtrf`
- `pdpbtrs`
- `pdpbtrsv`
- `pdpttrf`
- `pdpttrs`
- `pdpttrsv`
- `psdbtrf`
- `psdbtrs`
- `psdbtrsv`
- `psdttrf`
- `psdttrs`
- `psdttrsv`
- `psgbtrf`
- `psgbtrs`
- `pspbtrf`
- `pspbtrs`
- `pspbtrsv`
- `pspttrf`
- `pspttrs`
- `pspttrsv`
- `pzdbtrf`
- `pzdbtrs`
- `pzdbtrsv`
- `pzdttrf`
- `pzdttrs`
- `pzdttrsv`
- `pzgbtrf`
- `pzgbtrs`
- `pzpbtrf`
- `pzpbtrs`
- `pzpbtrsv`
- `pzpttrf`
- `pzpttrs`
- `pzpttrsv`

## Priority 2: redistribution-adjacent ScaLAPACK routines

These routines are close to the work already completed and are the
shortest path to broadening public `_I8` coverage beyond the reduction
drivers.

- `pcgesvd`
- `pdgesvd`
- `psgesvd`
- `pzgesvd`
- `pdsyev`
- `pssyev`
- `pcheev`
- `pzheev`
- `pdlaqr0`
- `pdlaqr2`
- `pdlaqr3`
- `pdlaqr4`
- `pslaqr0`
- `pslaqr2`
- `pslaqr3`
- `pslaqr4`

## Remaining public PBLAS families without `_I8`

Current missing public PBLAS surface: 93 entry points.

### Level 3 / matrix-matrix

- `pcgemm`, `pdgemm`, `psgemm`, `pzgemm`
- `pctrmm`, `pdtrmm`, `pstrmm`, `pztrmm`
- `pctrsm`, `pdtrsm`, `pstrsm`, `pztrsm`
- `pcsymm`, `pdsymm`, `pssymm`, `pzsymm`
- `pcsyrk`, `pdsyrk`, `pssyrk`, `pzsyrk`
- `pcsyr2k`, `pzsyr2k`
- `pchemm`, `pzhemm`
- `pcherk`, `pzherk`

### Level 2 / matrix-vector and rank updates

- `pcagemv`, `pdagemv`, `psagemv`, `pzagemv`
- `pcahemv`, `pzahemv`
- `pdasymv`, `psasymv`
- `pcatrmv`, `pdatrmv`, `psatrmv`, `pzatrmv`
- `pcgerc`, `pzgerc`
- `pcgeru`, `pzgeru`
- `pdger`, `psger`
- `pcher`, `pzher`
- `pcher2`, `pzher2`
- `pdsyr`, `pssyr`
- `pdsyr2`, `pssyr2`
- `pctranc`, `pztranc`
- `pdtran`, `pstran`
- `pctranu`, `pztranu`
- `pctradd`, `pdtradd`, `pstradd`, `pztradd`

### Level 1 / vector

- `pcamax`, `pdamax`, `psamax`, `pzamax`
- `pdasum`, `psasum`
- `pscasum`, `pdzasum`
- `pccopy`, `pdcopy`, `picopy`, `pscopy`, `pzcopy`
- `pcdotu`, `pzdotu`
- `pcswap`, `pdswap`, `psswap`, `pzswap`
- `pctrmv`, `pdtrmv`, `pstrmv`, `pztrmv`
- `pctrsv`, `pdtrsv`, `pstrsv`, `pztrsv`
- `pcgeadd`, `pdgeadd`, `psgeadd`, `pzgeadd`

## Remaining TOOLS surface without `_I8`

These are mostly utility routines rather than blockers for the current
cone.  They matter if the goal is broad, systematic `_I8` coverage.

### Scalar / local helpers

- `iceil`
- `ilacpy`
- `ilcm`
- `npreroc`
- `dsnrm2`
- `dscnrm2`
- `dsasum`
- `dscasum`
- `ssdot`
- `dddot`
- `ccdotu`
- `ccdotc`
- `zzdotu`
- `zzdotc`
- `slatcpy`
- `dlatcpy`
- `clatcpy`
- `zlatcpy`
- `smatadd`
- `dmatadd`
- `cmatadd`
- `zmatadd`
- `sltimer`

### Distributed utility families

- `pcchekpad`, `pdchekpad`, `pschekpad`, `pzchekpad`
- `pcfillpad`, `pdfillpad`, `psfillpad`, `pzfillpad`
- `pclaprnt`, `pdlaprnt`, `pslaprnt`, `pzlaprnt`
- `pclaread`, `pdlaread`, `pslaread`, `pzlaread`
- `pclawrite`, `pdlawrite`, `pslawrite`, `pzlawrite`
- `pcmatadd`, `pdmatadd`, `psmatadd`, `pzmatadd`
- `pccol2row`, `pdcol2row`, `pscol2row`, `pzcol2row`
- `pcrow2col`, `pdrow2col`, `psrow2col`, `pzrow2col`
- `pctreecomb`, `pdtreecomb`, `pstreecomb`, `pztreecomb`
- `pcelset2`, `pdelset2`, `pselset2`, `pzelset2`

### Integer tools / reshape

- `pichekpad`
- `pifillpad`
- `picol2row`
- `pirow2col`
- `pilaprnt`
- `pitreecomb`
- `pielget`
- `pielset`
- `pielset2`
- `reshape`
- `SL_gridreshape`
- `SL_init`

## Notes

- This file uses exact file-stem matching.  If a future `_I8` routine
  deliberately reuses a legacy stem through internal dispatch rather than
  adding a new `*_i8` file, update this manifest manually.
- Internal `PB_*` PBLAS helpers are intentionally excluded from the PBLAS
  list above.  This manifest is about public or phase-planning surface.
