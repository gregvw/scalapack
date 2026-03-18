# Phase 8 Notes: Serial/Tailored Path Wrappers + Test Coverage (Complete)

## Outcome

The tridiagonal reduction driver cone is fully I8-native.  All narrowing is
encapsulated in dedicated _I8 wrapper files.  Zero inline INT() narrowing
remains in the four reduction drivers.  Test coverage now exercises both the
blocked and serial/tailored paths.

## Part A: Test coverage improvements

### A1. Blocked-path driver tests (UPLO='U')

Added UPLO parameter to all driver test subroutines (RUN_D, RUN_S, RUN_C, RUN_Z)
and added UPLO='U' test cases alongside the existing 'L' cases.

UPLO='U' forces the blocked path in PDSYNTRD/PSSYNTRD/PCHENTRD/PZHENTRD because
the serial/tailored path only activates for UPLO='L'.  This directly exercises:
- PxLATRD_I8 (upper triangle variant)
- PxSYR2K_I8 / PxHER2K_I8 (blocked loop)
- PxSYTD2_I8 / PxHETD2_I8 (last-block wrapper, upper path)

### A2. Direct PxSYR2K/HER2K_I8 tests

Added TEST_DSYR2K and TEST_CHER2K to PBLAS/TESTING/test_pblas_i8.f.
Both compare I8 wrapper output against legacy, bit-identical.

PBLAS I8 kernel test suite now covers 10 operations across Level 1/2/3.

## Part B: Serial/tailored path wrappers

### 8 new thin wrapper files

Serial LAPACK wrappers (no descriptor, just N/LDA/LWORK):

| File | Wraps | Array types |
|------|-------|-------------|
| dsytrd_i8.f | DSYTRD | DOUBLE PRECISION |
| ssytrd_i8.f | SSYTRD | REAL |
| chetrd_i8.f | CHETRD | COMPLEX (D,E are REAL) |
| zhetrd_i8.f | ZHETRD | COMPLEX*16 (D,E are DOUBLE PRECISION) |

Tailored parallel wrappers (descriptor + N/IA/JA/LWORK):

| File | Wraps | Array types |
|------|-------|-------------|
| pdsyttrd_i8.f | PDSYTTRD | DOUBLE PRECISION |
| pssyttrd_i8.f | PSSYTTRD | REAL |
| pchettrd_i8.f | PCHETTRD | COMPLEX (D,E are REAL) |
| pzhettrd_i8.f | PZHETTRD | COMPLEX*16 (D,E are DOUBLE PRECISION) |

### 4 driver updates

Replaced inline narrowing blocks at serial/tailored call sites with wrapper calls.

Before (each driver):
```
IF( NPS8.GT.INTMAX .OR. ... ) THEN BLACS_ABORT END IF
NPS4 = INT(NPS8) ; N4 = INT(N)
CALL NARROW_DESC8(DESCB, DESCB4)
IF( NPROWB.EQ.1 ) THEN
   CALL xSYTRD( UPLO, N4, ..., NPS4, ..., INT(LLWORK), ... )
ELSE
   CALL PxSYTTRD( 'L', N4, ..., 1, 1, DESCB4, ..., INT(LLWORK), ... )
END IF
```

After:
```
IF( NPROWB.EQ.1 ) THEN
   CALL xSYTRD_I8( UPLO, N, ..., NPS8, ..., LLWORK, ... )
ELSE
   CALL PxSYTTRD_I8( 'L', N, ..., 1_8, 1_8, DESCB, ..., LLWORK, ... )
END IF
```

### Dead variable cleanup

Removed from all 4 drivers:
- N4, NPS4, NP4 (narrowed scalars — no longer used anywhere)
- DESCA4(9), DESCB4(9) (narrowed descriptors)
- NARROW_DESC8 (from EXTERNAL lists in complex drivers)
- Inline overflow-check + BLACS_ABORT blocks at serial/tailored call sites

The blocked-path overflow check is retained (PBLAS I8 wrappers still narrow
internally via pblas_i8_utils.h), but the call boundary is now inside the
wrappers, not inline in the driver.

## Validated platforms

- macOS arm64 (Apple clang): 8 I8 tests pass, full suite available
- x86_64 Linux (GCC): to be validated

## I8 surface additions

- Serial LAPACK wrappers: xSYTRD_I8 (2), xHETRD_I8 (2)
- Tailored parallel wrappers: PxSYTTRD_I8 (2), PxHETTRD_I8 (2)
- Total new files: 8 wrappers

## Hardening pass

### LWORK fail-fast fix

All 12 thin wrappers (4 Phase 7 + 8 Phase 8) originally clamped LWORK with
INT(MIN(LWORK, INTMAX)), silently truncating oversized workspace.  Changed to
reject LWORK > INTMAX with INFO = -1, matching the fail-fast semantics of the
removed driver-side overflow guards.

### MPI-backed I8 tools test (xi8tools_mpi)

New test target exercises paths that require a real BLACS grid:
- DESC_CONVERT_I8 2D→1D_H (on 1xP grid)
- DESC_CONVERT_I8 2D→1D_V (on Px1 grid)
- PCHK1MAT_I8 consistent case (all processes agree)
- PCHK1MAT_I8 inconsistent case (one process has wrong N, detected)

GLOBCHK_I8 is exercised implicitly via PCHK1MAT_I8.

### Repo hygiene

- Added docs/i8-refactor/scripts/.gitignore to exclude __pycache__/*.pyc

## Remaining work (beyond this cone)

The tridiagonal reduction cone is complete.  Future expansion paths:
- 1D-descriptor-heavy banded/tridiagonal family
- Broader public _I8 routine rollout beyond reduction drivers
- CI integration: run 9 I8 tests + i8_boundary_audit.py on I8-touching PRs
