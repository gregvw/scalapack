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

**Recommendation:** Start with Option A for PxAXPY and PxSCAL to establish the pattern
and validate the I8 PBLAS entry-point convention. The stage gate after step 2 determines
whether to continue with Option A or switch to Option B before scaling to larger kernels.

## Implementation order

 1. PxAXPY_I8 (all 4 types) — establish PBLAS I8 C entry-point pattern
 2. PxSCAL_I8 (all 4 types)

**--- Stage gate: decide duplication vs shared-core before proceeding ---**

 3. PxNRM2_I8 (4 types: PDNRM2, PSNRM2, PSCNRM2, PDZNRM2)
 4. PxLARFG_I8 (all 4 types) — first Fortran auxiliary using native I8 PBLAS
 5. PxLACGV_I8 (2 types: PCLACGV, PZLACGV) — complex vector conjugation
 6. PxDOT_I8 / PxDOTC_I8 (4 types)
 7. PxGEMV_I8 (all 4 types)
 8. PxSYMV_I8 / PxHEMV_I8 (4 types)
 9. PxLATRD_I8 (all 4 types) — first fully native I8 inner kernel
10. Update bridge drivers to use native PxLATRD_I8

## Acceptance criteria

- All existing tests (107+ on Linux) continue to pass.
- Each new PBLAS _I8 routine has a focused test comparing against its legacy counterpart.
- PxLATRD_I8 is bit-identical to legacy PxLATRD for 32-bit-sized inputs in all four types.
- At least one bridge driver family (all four PxSYNTRD_I8 / PxHENTRD_I8) is updated to
  call PxLATRD_I8 natively.
- The main blocked-iteration loop in each updated driver has zero hot-path narrowing for
  PxLATRD; only PxSYR2K/PxHER2K and the final PxSYTD2/PxHETD2 call remain as narrowed
  boundaries.
- The remaining narrowing points are documented explicitly in driver source comments.

## What this phase does NOT include

- PxSYR2K_I8 / PxHER2K_I8 (level-3 PBLAS — defer until PxLATRD is validated)
- PxSYTD2_I8 / PxHETD2_I8 (unblocked reduction — lower leverage)
- Full ILP64 PBLAS campaign
- Banded/tridiagonal 1D-descriptor family (still deferred per Phase 4 plan)
