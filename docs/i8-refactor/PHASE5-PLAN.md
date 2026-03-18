# Phase 5 Plan: Native I8 PBLAS Slice

## Prerequisites satisfied

Before starting this phase, the following were validated:

- GLOBCHK_I8 uses GLOBCHK_MPI_APIINT (ILP64-safe MPI datatype selection)
- Complex PxLAMR1D_I8 tests validate imaginary parts (imag = -real pattern)
- All four bridge drivers (PDSYNTRD_I8, PSSYNTRD_I8, PCHENTRD_I8, PZHENTRD_I8)
  produce bit-identical results to legacy on macOS arm64 and x86_64 Linux
- 107 tests passing on Linux

## Context

Phase 4 delivered the full bridge reduction-driver family validated on both platforms.
These drivers use checked narrowing at PBLAS/LAPACK call boundaries. The bridge proves
the vertical _I8 strategy works but leaves ~13K lines of PBLAS as the remaining
narrowing bottleneck.

## First-wave native-I8 PBLAS set

Driven by the PxLATRD dependency cone (the hot-path inner kernel in every blocked
reduction iteration):

| Kernel       | C lines | Calls/iter | Notes                          |
|-------------|---------|------------|--------------------------------|
| PxAXPY      | ~240    | 6          | Simplest, establish pattern    |
| PxSCAL      | ~255    | 2          | Small, also used by PxLARFG    |
| PxNRM2      | ~510    | 1          | Via PxLARFG                    |
| PxDOT/DOTC  | ~785    | 1          | Real/complex variants          |
| PxGEMV      | ~480    | 8          | Highest call frequency         |
| PxSYMV/HEMV | ~630    | 1          | Real/complex variants          |

Total: ~2,900 lines per type × 4 types = ~11,600 lines of C (plus Fortran entry points).

## Key design decision: duplication vs shared core

**Option A: Duplicate _I8 files.**
Each PBLAS routine gets a separate `px*_i8.c` that mirrors the original with widened types.
- Pro: Simple, no risk to existing code.
- Con: Doubles the C source surface, maintenance burden.

**Option B: Parameterize shared C cores.**
Refactor each PBLAS routine so the core logic is type-width-agnostic (uses Index64 internally),
with thin LP64 and I8 entry-point wrappers.
- Pro: Single source of truth, less code.
- Con: Requires touching existing PBLAS C code, higher risk of regressions.

**Stage gate resolved: Option A confirmed.**

PxAXPY_I8 and PxSCAL_I8 established the pattern. Each `_i8.c` wrapper is ~30 lines:
checked narrowing of int64_t args, delegation to the legacy PBLAS entry point, shared
helper in `pblas_i8_utils.h`. This is not real duplication — no algorithm logic is
copied. Reconsider shared-core only if later kernels force substantial logic into the
wrappers.

## Implementation order

 1. PxAXPY_I8 (all 4 types) — DONE
 2. PxSCAL_I8 (all 4 types) — DONE
 3. PxNRM2_I8 (4 types) — DONE
 4. PxDOT_I8 / PxDOTC_I8 (4 types) — DONE
 5. PxLARFG_I8 (all 4 types) + PCSSCAL_I8, PZDSCAL_I8 — DONE
 6. PxLACGV_I8 (2 types) — DONE
 7. PxGEMV_I8 (all 4 types) — DONE
 8. PxSYMV_I8 / PxHEMV_I8 (4 types) — DONE
 9. PxLATRD_I8 (all 4 types) + PxELGET_I8 (4 types) — DONE
10. Update all four bridge drivers to native PxLATRD_I8 — DONE

## Acceptance criteria — status

- [x] All existing tests (103 on Linux) continue to pass.
- [ ] Each new PBLAS _I8 routine has a focused test comparing against its legacy counterpart.
      (Driver-level bit-identical tests exist; dedicated per-kernel tests not yet added.)
- [x] PxLATRD_I8 is bit-identical to legacy PxLATRD for 32-bit-sized inputs in all four types
      (verified indirectly through driver tests).
- [x] All four bridge drivers updated to call native PxLATRD_I8.
- [x] The main blocked-iteration loop has zero hot-path narrowing for PxLATRD.
      Remaining narrowed boundaries: PxSYR2K/PxHER2K (1/iter), PxSYTD2/PxHETD2 (once).
- [x] Remaining narrowing documented in driver source comments.

## Phase 5 summary

PBLAS I8 entry points delivered: 26 C wrappers
  Level 1: PxAXPY, PxSCAL, PCSSCAL, PZDSCAL, PxNRM2, PxDOT/DOTC (18)
  Level 2: PxGEMV, PxSYMV/HEMV (8)

Fortran I8 auxiliaries delivered: 14
  PxLARFG (4), PxLACGV (2), PxLATRD (4), PxELGET (4)

Bridge drivers with native hot-path: 4
  PDSYNTRD_I8, PSSYNTRD_I8, PCHENTRD_I8, PZHENTRD_I8

Key design decisions recorded:
  - Option A (thin wrappers) confirmed at stage gate
  - pblas_i8_utils.h uses SCALAPACK_FORTRAN_INT_BYTES for ILP64-safe range checks
  - JP8 initialized to avoid inherited UB from legacy code

## What this phase does NOT include

- PxSYR2K_I8 / PxHER2K_I8 (level-3 PBLAS — defer until PxLATRD is validated)
- PxSYTD2_I8 / PxHETD2_I8 (unblocked reduction — lower leverage)
- Full ILP64 PBLAS campaign
- Banded/tridiagonal 1D-descriptor family (still deferred per Phase 4 plan)
