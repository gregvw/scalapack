# Phase 7 Notes: PxSYTD2_I8 / PxHETD2_I8 Thin Wrappers (Complete)

## Outcome

The last-block PxSYTD2/PxHETD2 call sites in all four reduction drivers are now
routed through thin I8 wrappers.  The blocked reduction path is fully I8-native
from entry to exit — zero inline narrowing remains on that code path.

## What was delivered

4 new thin wrapper files:

| File | Wraps | Array types |
|------|-------|-------------|
| pdsytd2_i8.f | PDSYTD2 | DOUBLE PRECISION |
| pssytd2_i8.f | PSSYTD2 | REAL |
| pchetd2_i8.f | PCHETD2 | COMPLEX (D,E are REAL) |
| pzhetd2_i8.f | PZHETD2 | COMPLEX*16 (D,E are DOUBLE PRECISION) |

Each wrapper:
1. Accepts INTEGER*8 for N, IA, JA, DESCA, LWORK; INFO stays INTEGER
2. Range-checks N, IA, JA vs INTMAX (defensive; N <= NB at call site)
3. Handles workspace query (LWORK = -1)
4. Calls NARROW_DESC8 + legacy routine

4 driver updates (upper + lower call sites in each):
  PDSYNTRD_I8, PSSYNTRD_I8, PCHENTRD_I8, PZHENTRD_I8

At upper path (last block after blocked loop):
  Before: CALL PxSYTD2(UPLO, MIN(N4,NB), A, INT(IA), INT(JA), DESCA4, ...)
  After:  CALL PxSYTD2_I8(UPLO, MIN(N,INT(NB,8)), A, IA, JA, DESCA, ...)

At lower path (trailing block):
  Before: CALL PxSYTD2(UPLO, KK, A, INT(IA+K8-1), INT(JA+K8-1), DESCA4, ...)
  After:  CALL PxSYTD2_I8(UPLO, INT(KK,8), A, IA+K8-1, JA+K8-1, DESCA, ...)

## Design rationale

N is bounded by NB (block size, always default INTEGER) at these call sites,
so a full I8 reimplementation of PxSYTD2 would add complexity with no practical
benefit.  The thin-wrapper pattern encapsulates the narrowing, keeping the driver
code clean and consistent with the PBLAS I8 wrapper approach.

## Remaining narrowed boundaries (outside the blocked path)

These are on the serial/tailored path, which operates on a redistributed copy
of the matrix in a square process sub-grid:

| Call site | When | Frequency |
|-----------|------|-----------|
| xSYTRD / xHETRD | Serial LAPACK (NPROWB = 1) | Once, serial path only |
| PxSYTTRD / PxHETTRD | Tailored parallel (NPROWB > 1) | Once, tailored path only |

## Validated platforms

- macOS arm64 (Apple clang): 108 tests, all passing
- x86_64 Linux (GCC): 104 tests, all passing

## Known coverage gap

Post-phase review identified that the driver tests (test_pdsyntrd_i8.f,
test_pchentrd_i8.f) use UPLO='L' with full queried workspace, which routes
through the serial/tailored path — not the blocked path.  The blocked loop
(PxLATRD_I8 → PxSYR2K_I8/PxHER2K_I8 → PxSYTD2_I8/PxHETD2_I8) is therefore
validated only indirectly.  PxSYR2K_I8/PxHER2K_I8 also lack dedicated tests
in test_pblas_i8.f.  Addressed in Phase 8 plan (Part A).

## I8 surface additions

Fortran I8 auxiliaries: PxSYTD2_I8 (2), PxHETD2_I8 (2)
Total I8 Fortran wrappers/routines: 16 (PxLARFG 4, PxLACGV 2, PxLATRD 4,
PxSYTD2_I8 2, PxHETD2_I8 2)
