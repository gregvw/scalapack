# Phase 2 Notes

## Scope

Phase 2 is an ABI-neutral C-side hardening pass.

- Public `Int` remains the Fortran-facing API integer width.
- The work in this phase is about making internal byte counts, offsets,
  workspace sizing, and selected internal products explicit and checked.
- This phase does not introduce a new public 64-bit ScaLAPACK interface.

Recent phase-2 commits:

- `78da0ff` Harden internal byte-count and offset arithmetic across BLACS/PBLAS/REDIST
- `84e135a` Harden REDIST interval and index arithmetic
- `84d046a` Check REDIST core origin-shift pointer arithmetic

## Semantic Types

- `ScaLAPACK_ApiInt`
  Public/API/descriptor integers. This is the semantic replacement for raw
  API-width `Int`.
- `ScaLAPACK_BlasInt`
  BLAS/LAPACK integer ABI type. This is configured independently, but mixed
  BLAS/API widths are still treated as unsupported for now.
- `ScaLAPACK_ByteCount`
  Byte sizes, allocation lengths, and `memcpy`/buffer sizing.
- `ScaLAPACK_Index64`
  Selected internal products and flattened offsets that should not rely on
  unchecked API-width integer arithmetic.
- `ScaLAPACK_UWord32`
  Exact 32-bit IEEE bit-pattern inspection in the `p{s,d}laiect` helpers.

## What Changed

- BLACS packed-buffer handling and MPI count/datatype wrappers were updated
  to support large internal buffer sizes without widening the public API.
- Shared checked helpers were added for `size_t`, byte-count, and `Index64`
  arithmetic.
- Allocation and temporary-buffer sizing was hardened across `BLACS`,
  `REDIST`, `TOOLS`, `SRC`, and `PBLAS`.
- Selected flattened pointer/offset arithmetic was converted to checked
  helper-based forms rather than raw `Int` products.
- `REDIST` interval sizing, local index arithmetic, and origin-shift pointer
  updates were hardened so they fail fast on overflow instead of silently
  wrapping.
- The IEEE helper files were refactored to use fixed-width 32-bit word logic
  instead of relying on `sizeof(Int) == 4`.

## What Was Intentionally Not Changed

- Public BLACS/PBLAS/ScaLAPACK interfaces still use API-width integers.
- Descriptor integer ABI remains unchanged.
- Broad Fortran-side integer widths were not changed in this phase.
- Mixed-width ScaLAPACK-vs-BLAS support is still not supported.
- This phase did not attempt a blanket removal of `Int` from the C sources.
  Many `Int` uses are still correct because they describe public API values,
  descriptor fields, process coordinates, leading dimensions, and related
  control data.

## How To Read Remaining `Int` Uses

A raw grep for `Int` is no longer a useful bug signal by itself.

Usually fine:

- public entrypoint arguments
- descriptor and BLACS context fields
- process-grid coordinates and ranks
- block sizes, leading dimensions, and loop/control variables tied to API
  parameters

Still worth auditing if seen in future work:

- `Int` products feeding allocation sizes
- `Int` products feeding pointer arithmetic or flattened local offsets
- `Int` used as a byte count, `memcpy` length, or MPI size/count surrogate

## Platform Note

On Apple Silicon macOS, `Accelerate` caused pathological SEP test behavior.
The CMake configuration now prefers `OpenBLAS` over `Accelerate` on
`Darwin arm64`, with an explicit opt-out for users who want the previous
behavior.

## Remaining Work

- Fortran/API-width audit and any future public 64-bit interface strategy
- Any remaining diminishing-return C audit items that are discovered later
- Optional mixed-width BLAS/API support, if that ever becomes a goal
