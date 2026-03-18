# Phase 6 Notes: Level-3 PBLAS I8 and Blocked-Loop Completion (Complete)

## Outcome

The main blocked-iteration loop in all four reduction driver families is now
fully I8-native. Zero narrowing occurs on the hot path.

## What was delivered

PBLAS I8 entry points added: 4 C wrappers
  Level 3: PDSYR2K_I8, PSSYR2K_I8, PCHER2K_I8, PZHER2K_I8

All four bridge drivers updated:
  PDSYNTRD_I8, PSSYNTRD_I8, PCHENTRD_I8, PZHENTRD_I8

Dead narrowed variables removed from all drivers:
  DESCW4, I4, J4, K4 (no longer referenced after LATRD + SYR2K/HER2K migration)

PBLAS I8 surface total: 30 C entry points.

## Hot-path call chain (fully I8-native)

Per blocked iteration, the driver calls:

    PxLATRD_I8
      → PxGEMV_I8 (8 calls)
      → PxSYMV_I8 / PxHEMV_I8
      → PxDOT_I8 / PxDOTC_I8
      → PxSCAL_I8
      → PxAXPY_I8
      → PxLARFG_I8 → PxNRM2_I8 + PxSCAL_I8
      → PxELSET_I8, PxELGET_I8
      → PxLACGV_I8 (complex only)

    PxSYR2K_I8 / PxHER2K_I8

No narrowing occurs in this path. All integer arguments, descriptors,
and local index arithmetic are INTEGER*8 through PxLATRD_I8. The
SYR2K/HER2K wrappers narrow internally to legacy PBLAS internals
via pblas_i8_utils.h, but the caller-facing interface is fully I8.

## Remaining narrowed boundaries (outside the blocked loop)

These execute once per reduction call, not once per iteration:

| Call site | When | Frequency |
|-----------|------|-----------|
| PxSYTD2 / PxHETD2 | Last block of blocked path | Once per reduction |
| xSYTRD / xHETRD | Serial LAPACK (small N or large workspace) | Once, serial path only |
| PxSYTTRD / PxHETTRD | Tailored parallel (large workspace, lower triangle) | Once, tailored path only |

All narrowed at call boundary with overflow check + BLACS_ABORT.

## Validated platforms

- macOS arm64 (Apple clang): 8 ctest I8 targets, all passing
- x86_64 Linux (GCC): 104 total tests, all passing
- x86_64 Linux (clang-20 + ASan/UBSan): zero I8-code sanitizer errors

## Test coverage

8 ctest targets:
  xi8tools, xdgemr_i8, xdtrmr_i8, xdsyntrd_i8, xchentrd_i8,
  xdlamr1d_i8, xdlamve_i8, xpblas_i8

xpblas_i8 covers: PDAXPY, PDSCAL, PDDOT, PDNRM2, PCDOTC, PDGEMV,
PDSYMV, PCHEMV (8 kernel-level bit-identical comparisons).

Driver tests verify D, E, TAU bit-identical to legacy for all four types.
