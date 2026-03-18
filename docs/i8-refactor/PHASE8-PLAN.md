# Phase 8 Plan: Serial/Tailored Path Wrappers + Test Coverage

## Context

After Phase 7, the blocked reduction path is fully I8-native.  The two remaining
narrowed boundaries are on the serial/tailored path — a separate code path that
redistributes the matrix into a square sub-grid and reduces it there.

These calls execute at most once per reduction, operate on a redistributed
workspace copy (not the original distributed matrix), and the sub-grid dimensions
are bounded by SQNPC (square root of process count).  The practical overflow risk
is zero, but wrapping them completes the I8 encapsulation of all call sites.

Review also identified that the existing tests do not actually exercise the blocked
path or the Phase 6 SYR2K/HER2K wrappers directly, and that the I8 tools layer
lacks MPI-backed tests for DESC_CONVERT_I8 2D->1D and PCHK*MAT_I8/GLOBCHK_I8.
This phase addresses both the wrapper work and the coverage gaps.

## Part A: Test coverage improvements

### A1. Force blocked-path driver tests

The current driver tests (test_pdsyntrd_i8.f, test_pchentrd_i8.f) use UPLO='L'
with queried workspace, which always takes the serial/tailored path.  The blocked
path is never exercised directly.

The driver branch logic (pdsyntrd_i8.f:163-171):
- Serial: (N < MINSZ or SQNPC == 1) AND LLWORK >= ONEPMIN AND UPLO='L'
- Tailored: LLWORK >= TTLWMIN AND UPLO='L'
- Blocked: everything else (UPLO='U', or LWORK < TTLWMIN)

Two approaches to force the blocked path:
1. **UPLO='U'**: always takes the blocked path (simplest)
2. **Restrict LWORK**: query, then pass only LWMIN (blocked minimum) instead
   of TTLWMIN.  More surgical but fragile if workspace accounting changes.

Plan: Add UPLO='U' test cases to both test files.  This exercises:
- PxLATRD_I8 (upper triangle variant)
- PxSYR2K_I8 / PxHER2K_I8 (called from blocked loop)
- PxSYTD2_I8 / PxHETD2_I8 (last-block wrapper, upper path)

### A2. Direct PxSYR2K/HER2K_I8 tests in test_pblas_i8.f

The PBLAS I8 test covers Level 1 (AXPY, SCAL, DOT, NRM2, DOTC) and Level 2
(GEMV, SYMV, HEMV) but not Level 3.  Add TEST_DSYR2K and TEST_CHER2K
subroutines comparing PxSYR2K_I8/PxHER2K_I8 against legacy.

### A3. MPI-backed I8 tools test (deferred — optional)

DESC_CONVERT_I8 2D->1D, PCHK1MAT_I8, PCHK2MAT_I8, and GLOBCHK_I8 need a real
BLACS context.  This requires a new MPI test target (xi8tools_mpi).  Lower
priority than A1/A2 since these routines are well-exercised indirectly by the
driver tests.  Defer unless time permits.

## Part B: Serial/tailored path wrappers

### Current state of the serial/tailored path

In each driver (PDSYNTRD_I8 etc.), the path is:

```
PxTRMR2D_I8  →  redistribute A into WORK(INDB)
                 ↓
           NPROWB == 1?
          /              \
    xSYTRD/xHETRD    PxSYTTRD/PxHETTRD
    (serial LAPACK)   (tailored parallel)
                 ↓
PxLAMR1D_I8  →  redistribute D, E, TAU back
PxTRMR2D_I8  →  redistribute A back
```

The narrowing occurs at the xSYTRD/PxSYTTRD call site:
- DESCB narrowed via NARROW_DESC8
- N, NPS, LLWORK narrowed to N4, NPS4, INT(LLWORK)
- Overflow check + BLACS_ABORT before narrowing

### 8 thin wrapper files (same pattern as Phase 7)

Serial LAPACK wrappers:

| File | Wraps | Notes |
|------|-------|-------|
| SRC/dsytrd_i8.f | DSYTRD | No descriptor — just N, LDA, LWORK |
| SRC/ssytrd_i8.f | SSYTRD | |
| SRC/chetrd_i8.f | CHETRD | |
| SRC/zhetrd_i8.f | ZHETRD | |

These are trivially thin: narrow N, LDA, LWORK and call the LAPACK routine.
No descriptors involved.

Tailored parallel wrappers:

| File | Wraps | Notes |
|------|-------|-------|
| SRC/pdsyttrd_i8.f | PDSYTTRD | Same signature as PxSYTD2 wrappers |
| SRC/pssyttrd_i8.f | PSSYTTRD | |
| SRC/pchettrd_i8.f | PCHETTRD | |
| SRC/pzhettrd_i8.f | PZHETTRD | |

These narrow N, IA, JA, DESCA, LWORK like the PxSYTD2_I8 wrappers.
Note: the drivers always pass IA=1, JA=1 for the sub-grid copy, so
the I8 values are always small.

### 4 driver updates

Replace inline narrowing at the serial/tailored call sites with wrapper calls.

In PDSYNTRD_I8 (and analogously for the other 3), the block:
```fortran
      IF( NPS8.GT.INTMAX .OR. ... ) THEN
         CALL PXERBLA / BLACS_ABORT
      END IF
      NPS4 = INT( NPS8 )
      N4   = INT( N )
      CALL NARROW_DESC8( DESCB, DESCB4 )
      IF( NPROWB.EQ.1 ) THEN
         CALL DSYTRD( UPLO, N4, WORK(INDB), NPS4, ... INT(LLWORK), IINFO )
      ELSE
         CALL PDSYTTRD( 'L', N4, WORK(INDB), 1, 1, DESCB4, ... INT(LLWORK), IINFO )
      END IF
```
becomes:
```fortran
      IF( NPROWB.EQ.1 ) THEN
         CALL DSYTRD_I8( UPLO, N, WORK(INDB), NPS8, ... LLWORK, IINFO )
      ELSE
         CALL PDSYTTRD_I8( 'L', N, WORK(INDB), 1_8, 1_8, DESCB, ... LLWORK, IINFO )
      END IF
```

After this, the N4, NPS4, DESCB4 variables may become dead in the drivers
and can be removed.  Also verify whether DESCA4 (used only by PxSYTD2 before
Phase 7) is now fully dead.

### Build and doc updates

- SRC/CMakeLists.txt: add 8 new files
- CLAUDE.md: add to I8 surface, remove remaining narrowed boundaries table
- PHASE8-NOTES.md: record completion

## Execution order

1. A1 + A2 first (test coverage) — validates the pre-existing blocked path
2. B (wrappers + driver updates) — then re-run to confirm no regression
3. Dead variable cleanup in drivers
4. Doc updates (CLAUDE.md, PHASE8-NOTES.md)

## What this does NOT touch

- The PJLAENV calls that query ANB/MINSZ use the tailored routine's name
  string.  These return INTEGER and are bounded by block size — no wrapping
  needed.
- BLACS_GET, BLACS_GRIDINIT, BLACS_GRIDEXIT, BLACS_GRIDINFO — these are
  context management calls with inherently small INTEGER arguments.
- The blocked path (already fully I8-native after Phase 7).

## Verification

1. `cmake --build . --target scalapack` (build succeeds)
2. `ctest -R i8` (all I8 tests pass, including new blocked-path and SYR2K tests)
3. `ctest` (full suite, no regressions)
4. `python3 docs/i8-refactor/scripts/i8_boundary_audit.py SRC/pdsyntrd_i8.f`
   (zero INT() narrowing findings at DSYTRD/PDSYTTRD calls)
5. Confirm DESCA4/N4/NPS4/DESCB4 are dead in drivers and removed

## Risk assessment

Low.  The serial/tailored path arguments are all bounded:
- N <= matrix dimension (already validated at driver entry)
- NPS <= NUMROC(N, 1, ..., SQNPC) — bounded by N/SQNPC
- IA, JA are always 1 on the sub-grid copy
- LLWORK is workspace length, validated at entry

The wrappers add defensive range checks but the real value is
encapsulation: all narrowing lives in dedicated _I8 wrapper files
rather than scattered inline in the drivers.
