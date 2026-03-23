# Remaining I8 Surface

Public routine families that do not yet have `_I8` variants.

## Recently completed (Phase 12)

These families now have `_I8` thin wrappers:

- `pdgeqrf`, `psgeqrf`, `pcgeqrf`, `pzgeqrf` — QR factorization
- `pdorgqr`, `psorgqr`, `pcungqr`, `pzungqr` — explicit Q generation
- `pdormqr`, `psormqr`, `pcunmqr`, `pzunmqr` — Q multiply
- `pdgetri`, `psgetri`, `pcgetri`, `pzgetri` — matrix inverse

## SVD cone

Requires `PxGEBRD_I8` and `PxORMBR_I8` as foundational pieces
before the driver can be made native I8.

- `pcgesvd`, `pdgesvd`, `psgesvd`, `pzgesvd`

## Schur / eigenvalue variants

Deferred — recursive/index-heavy families requiring deep refactoring.

- `pdlaqr0`, `pdlaqr2`, `pdlaqr3`, `pdlaqr4`
- `pslaqr0`, `pslaqr2`, `pslaqr3`, `pslaqr4`
- `pdsyevd`, `pssyevd`, `pcheevd`, `pzheevd`
- `pdsyevx`, `pssyevx`, `pcheevx`, `pzheevx`
- `pdsyevr`, `pssyevr`, `pcheevr`, `pzheevr`
- `pdsygvx`, `pssygvx`, `pchegvx`, `pzhegvx`

## Public PBLAS families without `_I8`

69 entry points.

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

## TOOLS surface without `_I8`

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

## Internal narrowing boundaries

- `PxLANSY` / `PxLANHE` — blocks the eigensolver cone from supporting
  N > INTMAX.  This is the only remaining narrowing boundary in a
  solver-level call tree.

## Notes

- Internal `PB_*` PBLAS helpers are intentionally excluded.
  This manifest covers public or phase-planning surface only.
- This file uses exact file-stem matching.  If a future `_I8` routine
  deliberately reuses a legacy stem through internal dispatch rather than
  adding a new `*_i8` file, update this manifest manually.
