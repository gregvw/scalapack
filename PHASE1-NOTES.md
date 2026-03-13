# Phase 1 Notes

## Scope

Phase 1 established the large-count and semantic-type foundation for the
later C-side hardening work.

- Public `Int` remained the Fortran-facing API integer width.
- The goal in this phase was to make large internal MPI buffer/count paths
  safe on newer MPI implementations and to introduce the semantic typing and
  ABI checks needed for later work.
- This phase also captured the Apple Silicon/OpenBLAS platform requirement
  discovered during verification.

Key phase-1 commits:

- `2ee2f60` Add ILP64 support for MPI operations via MPI 4.0+ wrappers
- `5726d57` Fixed platform-dependent pointer conversion warnings
- `b687f01` Prefer OpenBLAS over Accelerate on Apple Silicon and fix large-count MPI paths
- `51955af` removing Int and MpiInt from BI_TransUserComm where not appropriate for MPI
- `3f885e8` wrapper for MPI_Type_indexed and MpiInt for BI_GetMpiTrType
- `88e10a6` removing Int and MpiInt from additional places where not appropriate for MPI

## Semantic Types

Phase 1 introduced the main semantic type layer used by later work:

- `ScaLAPACK_ApiInt`
  Public/API/descriptor integers matching the default Fortran `INTEGER` ABI.
- `ScaLAPACK_BlasInt`
  BLAS/LAPACK integer ABI type.
- `ScaLAPACK_BufLen` / `ScaLAPACK_ByteCount`
  Internal byte-count and packed-buffer sizing.
- `ScaLAPACK_Index64`
  Internal widened arithmetic helper type for selected products/offsets.
- `ScaLAPACK_UWord32`
  Exact 32-bit IEEE bit-pattern inspection.

## What Changed

- Added MPI-4-aware wrapper paths for large-count MPI operations used by
  BLACS, including pack/unpack/count/datatype-related calls.
- Introduced a separate internal packed-buffer byte-length path so BLACS
  packed buffers can exceed legacy 32-bit size assumptions without widening
  the public API integer type.
- Added configure-time Fortran and BLAS integer-width checks and tied the C
  semantic typedefs to those checks.
- Added compile-time assertions for the semantic integer contracts.
- Refactored the IEEE helper files away from raw `Int`-width assumptions and
  into fixed-width 32-bit word handling.
- Tightened several MPI boundary sites where raw `Int` or `MpiInt` usage was
  not semantically correct, including communicator/group helper paths and MPI
  datatype helper plumbing.
- Added the Apple Silicon/OpenBLAS CMake policy after Accelerate showed
  pathological SEP test behavior on `Darwin arm64`.

## What Was Intentionally Not Changed

- Public BLACS/PBLAS/ScaLAPACK routine signatures still use API-width
  integers.
- Descriptor integer ABI remains unchanged.
- Mixed-width ScaLAPACK-vs-BLAS support is still treated as unsupported.
- Phase 1 did not attempt the broad C-side allocator/offset cleanup that was
  handled later in phase 2.
- Phase 1 also did not attempt a Fortran-wide integer-width refactor.

## Platform Note

On Apple Silicon macOS, `Accelerate` caused pathological SEP test behavior.
The CMake configuration now prefers `OpenBLAS` over `Accelerate` on
`Darwin arm64`, with an explicit opt-out for users who want the previous
behavior.

## Relationship To Phase 2

Phase 1 created the semantic types, ABI checks, MPI wrapper layer, and large
packed-buffer foundation. Phase 2 then propagated those semantics through the
main C allocation, sizing, and offset paths across `BLACS`, `REDIST`,
`TOOLS`, `SRC`, and `PBLAS`.

## Remaining Work After Phase 1

- Broader ABI-neutral C audit and hardening work, completed later in phase 2
- Fortran/API-width audit and any future public 64-bit interface strategy
- Optional mixed-width BLAS/API support, if that ever becomes a goal
