# I8 Refactor Phase History

Consolidated design rationale and key decisions from Phases 1–13.
Each section preserves the reasoning behind architectural choices that
would not be obvious from reading the code alone.

---

## Phase 1: Semantic Types and MPI Large-Count Foundation

Established the type foundation for all later work.

### Semantic types introduced

| Type | Purpose |
|------|---------|
| `ScaLAPACK_ApiInt` | Public/API/descriptor integers (Fortran `INTEGER` ABI) |
| `ScaLAPACK_BlasInt` | BLAS/LAPACK integer ABI type |
| `ScaLAPACK_BufLen` / `ScaLAPACK_ByteCount` | Internal byte counts and packed-buffer sizing |
| `ScaLAPACK_Index64` | Widened arithmetic for products/offsets |
| `ScaLAPACK_UWord32` | Exact 32-bit IEEE bit-pattern inspection |

### Key changes

- MPI-4-aware wrapper paths for large-count MPI operations in BLACS.
- Separate internal packed-buffer byte-length path (buffers can exceed 2^31
  without widening the public API integer type).
- Configure-time Fortran and BLAS integer-width checks tied to C typedefs.
- IEEE helper files refactored from `Int`-width assumptions to fixed-width
  32-bit word handling.

### Platform note

On Apple Silicon macOS, Accelerate caused pathological SEP test behavior.
CMake now prefers OpenBLAS over Accelerate on `Darwin arm64`.

---

## Phase 2: ABI-Neutral C-Side Hardening

Propagated semantic types through internal allocation, sizing, and offset
paths across BLACS, REDIST, TOOLS, SRC, and PBLAS.

### Design principle

A raw grep for `Int` is not a bug signal — many uses are correct (API args,
descriptor fields, process coordinates, block sizes). The audit focused on:
- `Int` products feeding allocation sizes
- `Int` products feeding pointer arithmetic or flattened local offsets
- `Int` used as byte count, `memcpy` length, or MPI size/count

### Key changes

- Allocation and temporary-buffer sizing hardened across all subsystems.
- REDIST interval sizing, local index arithmetic, and origin-shift pointer
  updates converted to fail-fast on overflow.
- Shared checked helpers added for `size_t`, byte-count, and `Index64`
  arithmetic.

---

## Phase 3: I8 Redistribution Layer

Built PxGEMR2D_I8 and PxTRMR2D_I8 — the first I8 entry points.

### Descriptor model

DTYPE encodes layout (1=2D, 501=1D-H, 502=1D-V). Width lives in symbol
names (`_I8` suffix), not in DTYPE. I8 callers pass `int64_t[9]` descriptors
with the same 9-field layout as legacy.

### Core types (`redist_core.h`)

| Type | Purpose |
|------|---------|
| `MDESC_CORE` | Internal descriptor; dimension fields are `Index64`, grid fields are `int` |
| `IDESC_CORE` | Interval descriptor for gemr (`lstart` = local memory offset) |
| `IDESC_TR_CORE` | Interval descriptor for trmr (`gstart` = global position) |

### Thin wrapper pattern (established here, reused throughout)

Each type variant has three entry points:
1. `Cpxtrmr2d` (legacy C API) — widens to `MDESC_CORE`, calls `_core`
2. `Cpxtrmr2d_core` (internal) — all arithmetic uses `Index64`
3. `pxtrmr2d_i8_` (Fortran I8 entry) — unpacks `int64_t[9]`, calls `_core`

### Key difference: gemr vs trmr intervals

- **gemr** uses `lstart` (local memory offset) because pack/unpack works
  directly between local storage and message buffers.
- **trmr** uses `gstart` (global position) because triangle geometry is
  computed in global coordinates.

---

## Phase 4: Validation Layer and Descriptor Contracts

Delivered CHK1MAT_I8, PCHK1MAT_I8, PCHK2MAT_I8, DESC_CONVERT_I8,
GLOBCHK_I8, and the xi8tools/xi8tools_mpi test targets.

### Key decision: DTYPE semantics

I8 2D descriptors use semantic DTYPE=1 (same as legacy). Width is a
compile-time property of the routine name, not a runtime descriptor tag.

### GLOBCHK_I8

Uses a small C helper (`globchk_mpi_i8.c`) for MPI_Allreduce over int64_t
(value, position) pairs, avoiding Fortran-side MPI type complexity.
Uses `GLOBCHK_MPI_APIINT` for ILP64-safe MPI datatype selection.

---

## Phase 5: First-Wave Native I8 PBLAS

Delivered 26 PBLAS I8 C wrappers + Fortran auxiliaries, driven by the
PxLATRD dependency cone.

### Stage gate: duplication vs shared core

**Option A (thin wrappers) confirmed.** Each `_i8.c` is ~30 lines: checked
narrowing of int64_t args via `pblas_i8_utils.h`, delegation to legacy.
No algorithm logic is copied. Shared-core (Option B) was rejected as
higher risk for zero practical benefit.

### `pblas_i8_utils.h` design

- `pblas_i8_narrow()`: narrows single int64_t, returns 0 on overflow
- `pblas_i8_narrow_desc()`: narrows all 9 descriptor entries
- `pblas_i8_abort()`: PBLAS error reporting on overflow
- Uses `SCALAPACK_FORTRAN_INT_BYTES` for ILP64-safe range checks

### Workspace reduction convention

All drivers use `DBLE`/`DGAMN2D` (never `REAL`/`SGAMN2D`) for I8-safe
global workspace minimum computation.

---

## Phase 6: Blocked-Loop Completion

Added PxSYR2K_I8/PxHER2K_I8. The main blocked-iteration loop in all four
reduction drivers became fully I8-native with zero hot-path narrowing.

---

## Phase 7: Unblocked Panel Wrappers (PxSYTD2_I8/PxHETD2_I8)

### Design rationale

N is bounded by NB (block size, always default INTEGER) at the last-block
call sites. Full I8 reimplementation would add complexity with no practical
benefit — thin wrappers encapsulate the narrowing cleanly.

---

## Phase 8: Serial/Tailored Path + Test Hardening

Completed the tridiagonal reduction cone. Zero inline INT() narrowing
remains in the four reduction drivers.

### Key additions

- xSYTRD_I8/xHETRD_I8 (4) — serial LAPACK wrappers
- PxSYTTRD_I8/PxHETTRD_I8 (4) — tailored parallel wrappers
- UPLO='U' test cases added (exercises blocked path directly)
- PxSYR2K/HER2K direct kernel tests added to xpblas_i8

### LWORK fail-fast fix

All thin wrappers changed from clamping LWORK to rejecting LWORK > INTMAX
with INFO=-1 (fail-fast semantics). Note: Phase 12 later adopted LWORK
capping instead (see Phase 12 notes).

---

## Phases 9–11: Banded Solvers, Eigensolvers, Dense Solvers, Large-N Closure

### Phase 9: Banded/Tridiagonal I8 (56 thin wrappers)

DB, DT, GB, PB, PT families — all 4 types × factor/solve/solve-vector.
Thin wrappers that narrow all arguments; do not support N > INTMAX.

### Phase 10: Eigensolver cone

Native I8 drivers PxSYEV_I8 (2), PxHEEV_I8 (2) with full internal I8
call tree. One remaining narrowing boundary: PxLANSY/PxLANHE (matrix norm)
causes an INTMAX entry guard.

Support wrappers: PxLASCL_I8 (4), PxLASET_I8 (4), PxORMTR_I8 (2),
PxUNMTR_I8 (2).

Bugs found: PZHEEV_I8 workspace sizing used wrong INT kind; PCHEEV_I8 and
PZHEEV_I8 were missing RWORK/LRWORK args to reduction drivers.

### Phase 11: Dense direct solvers + large-N closure

**Cholesky:** PxPOTRF_I8 (4), PxPOTRS_I8 (4), PxPOSV_I8 (4) — native I8.

**LU:** PxGETRF_I8 (4), PxGETRS_I8 (4), PxGESV_I8 (4) + support
PxLASWP_I8 (4), PxLAPIV_I8 (4) — native I8.

**Large-N closure:** Native I8 unblocked panel routines PxPOTF2_I8 (4) and
PxGETF2_I8 (4) eliminated INTMAX entry guards from blocked drivers.
Additional PBLAS I8: PxAMAX_I8 (4), PxSWAP_I8 (4), PxGER/GERU_I8 (4).

PBLAS total reached 54 C entry points. Dense solver cones became fully
large-N capable (N > 2^31-1).

### Cross-cutting fixes

- PBLAS I8 early returns removed (all wrappers — legacy handles edge cases)
- NARROW_DESC8 changed to abort on overflow (fail-fast, not clamp)
- Thin wrapper quick returns removed (80 files) to preserve workspace queries
- Edge-case regression tests: K=0 beta scaling, bad-UPLO on N=0

---

## Phase 12: QR Cone, Matrix Inverse, PBLAS COPY, Matrix Norms

Thin wrappers enabling ButterflyPACK's core compression and factorization
paths.

### QR factorization cone (12 wrappers)

- PxGEQRF_I8 (4) — QR factorization
- PxORGQR_I8 (2) / PxUNGQR_I8 (2) — explicit Q generation
- PxORMQR_I8 (2) / PxUNMQR_I8 (2) — Q multiply

### Supporting additions

- PxGETRI_I8 (4) — matrix inverse
- PxLANGE_I8 (4) — matrix norms (first FUNCTION I8 wrappers)
- PxCOPY_I8 (4) — PBLAS L1 vector copy (PBLAS total: 58)

### LWORK policy change

Phase 12 wrappers cap LWORK > INTMAX at INTMAX (instead of rejecting).
Rationale: if the caller has 3B of workspace, we shouldn't fail — legacy
will only use what it needs. Pre-existing wrappers (e.g., PDORMTR_I8)
still reject; this inconsistency is noted for future harmonization.

### INFO on overflow

All Phase 12+ wrappers set INFO=-1 before calling PXERBLA on INTMAX
overflow (code review finding — earlier wrappers left INFO=0).

---

## Phase 13: SVD Cone

Completed the full SVD dependency chain as thin wrappers:

- PxORMLQ_I8 (2) / PxUNMLQ_I8 (2) — LQ multiply (needed by PxORMBR)
- PxORMBR_I8 (2) / PxUNMBR_I8 (2) — bidiag orthogonal multiply
- PxGEBRD_I8 (4) — bidiagonal reduction
- PxGESVD_I8 (4) — SVD driver (complex variants include RWORK)

Note: PDGESVD only accepts JOBU/JOBVT = 'V' or 'N' (no 'S' thin-SVD
option exists in ScaLAPACK).

This completed all ScaLAPACK routines needed by ButterflyPACK's
compression (`PComputeRange`), SVD truncation (`PSVD_Truncate`), and
block inversion paths.
