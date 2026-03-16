# Phase 3: I8 Redistribution Layer

Notes on the I8 (64-bit descriptor) support added to PxGEMR2D and PxTRMR2D.

## Design

### Descriptor model

I8 callers pass a `BLOCK_CYCLIC_2D_I8 = 501` descriptor as `int64_t[9]`.
The 9 fields mirror the standard descriptor but dimension fields
(M, N, MB, NB, LLD) are 64-bit.  Grid-coordinate fields (DTYPE, CTXT,
RSRC, CSRC) stay semantically 32-bit; they are narrowed with bounds
checking on entry.

### Internal core types

All redistribution logic operates on widened internal types defined in
`REDIST/SRC/redist_core.h`:

| Type            | Purpose                                         |
|-----------------|-------------------------------------------------|
| `MDESC_CORE`    | Internal descriptor; dimension fields are `ScaLAPACK_Index64`, grid fields are `int` |
| `IDESC_CORE`    | Interval descriptor for gemr (`lstart` = local memory offset) |
| `IDESC_TR_CORE` | Interval descriptor for trmr (`gstart` = global position in submatrix) |

Two unpackers convert ABI-facing descriptors into `MDESC_CORE`:
- `unpack_desc_i4` -- from `ScaLAPACK_ApiInt[9]` (always succeeds, widening)
- `unpack_desc_i8` -- from `int64_t[9]` (bounds-checks grid fields)

### Thin wrapper pattern

Each type variant (d/s/c/z/i) has three entry points:

1. **`Cpxtrmr2d`** (legacy C API) -- unpacks `MDESC` to `MDESC_CORE` via
   field-by-field widening, delegates to `_core`.
2. **`Cpxtrmr2d_core`** (internal) -- all arithmetic uses `ScaLAPACK_Index64`;
   BLACS message counts are narrowed via `i64_to_blacs_count`.
3. **`pxtrmr2d_i8_`** (Fortran I8 entry) -- unpacks `int64_t[9]` descriptors
   via `unpack_desc_i8`, calls `_core`.

### Shared helpers (`pgemraux_core.c`)

| Function                   | Purpose |
|----------------------------|---------|
| `localsize_core`           | Local row/column count for block-cyclic distribution |
| `memoryblocksize_core`     | Total local block size |
| `changeorigin_core`        | Shift submatrix origin |
| `paramcheck_core`          | Validate descriptor against process grid |
| `localindice_core`         | Local memory index of global element |
| `scan_intervals_core`      | Block intersection scan for gemr (stores `lstart`) |
| `scan_intervals_tr_core`   | Block intersection scan for trmr (stores `gstart`) |
| `insidemat_core`           | Triangle element geometry (shared by all trmr types) |
| `redist_sync_params_i8`    | `MPI_Allreduce(MPI_MIN)` sync of 64-bit parameter array |

### Key difference: gemr vs trmr intervals

- **gemr** uses `IDESC_CORE` with `lstart` (local memory offset) because
  `block2buff_core`/`buff2block_core` copy directly between local storage
  and message buffers.
- **trmr** uses `IDESC_TR_CORE` with `gstart` (global position in submatrix)
  because `scanD0_core` works in global coordinates and calls
  `localindice_core` to convert per-column.

## Files modified

### Shared infrastructure
- `REDIST/SRC/redist_core.h` -- core types, unpacker inlines, function declarations
- `REDIST/SRC/pgemraux_core.c` -- shared helper implementations

### PxGEMR2D_I8 (general redistribution)
- `REDIST/SRC/p{d,s,c,z,i}gemr.c` -- each file gets `_core` extraction + `_i8_` entry
- `REDIST/TESTING/test_pdgemr2d_i8.c` -- caller-side test

### PxTRMR2D_I8 (triangular redistribution)
- `REDIST/SRC/p{d,s,c,z,i}trmr.c` -- each file gets `_core` extraction + `_i8_` entry
- `REDIST/TESTING/test_pdtrmr2d_i8.c` -- caller-side test

## Tests

### test_pdgemr2d_i8 (`xdgemr_i8`)

Two test cases, run with `mpirun -np 4`:

1. **Big-block descriptor**: M=10, N=5, MB=NB=3e9 (> INT_MAX).  Single-owner
   round-trip on a 1xN grid.  Exercises `unpack_desc_i8`,
   `redist_sync_params_i8`, and 64-bit arithmetic in core helpers.

2. **Cross-rank redistribution**: M=12, N=8.  Two different grids (2x2 and 1x4)
   with different block sizes (3x2 -> 4x3).  Forward redistribution, verify
   element-level correctness on destination grid, then inverse and verify
   round-trip.

### test_pdtrmr2d_i8 (`xdtrmr_i8`)

Two test cases, each testing all 4 uplo/diag combinations
(upper/lower x unit/non-unit), run with `mpirun -np 4`:

1. **Big-block descriptor**: M=N=8, MB=NB=3e9.  Single-owner round-trip.
   Verifies that only triangular elements are redistributed and
   non-triangular elements remain unchanged (sentinel values).

2. **Cross-rank redistribution**: M=N=12, grids 2x2 and 1x4, block sizes
   3x2 -> 4x3.  Round-trip with element-level triangle-aware verification.

## Verification state

- Linux (Open MPI 4.1.6, GCC): 98/98 pass
- macOS (MPICH 4.3.2, AppleClang): 102/102 pass

## What stays next

The Fortran/API boundary work is deferred to a subsequent session:
- `CHK1MAT_I8` and descriptor-validation helpers
- First PBLAS/ScaLAPACK I8 routine surface beyond redistribution
- Broader descriptor rollout strategy
