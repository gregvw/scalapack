# MPI Large-Count Scope for ILP64

This document scopes the MPI side of the `_I8` / ILP64 work.

The central distinction is:

- ScaLAPACK integer-width cleanup is about 64-bit-safe descriptor, index,
  workspace, and local-size arithmetic.
- MPI large-count support is about whether a single communication call can
  legally move more than `INT_MAX` elements.

These are related, but they are not the same problem.

## Current branch contract

The current branch does **not** implement a chunked fallback transport for
legacy MPI count APIs.

Instead, ILP64 builds require the MPI large-count `_c` entry points used by
BLACS, and configuration must fail hard if they are not available.

That is the current intended contract:

- `SCALAPACK_ENABLE_ILP64=ON`
  - requires 8-byte default Fortran `INTEGER`
  - requires an ILP64 BLAS/LAPACK
  - requires MPI large-count `_c` APIs in the C MPI interface

This is implemented by:

- `CMAKE/check_mpi_large_count.c`
- `CMAKE/CheckMpiLargeCount.cmake`

and enforced from `CMakeLists.txt`.

## Why this branch requires `_c` APIs

BLACS currently routes ILP64 communication through wrappers in
`BLACS/SRC/Bconfig.h`.

In ILP64 mode those wrappers use large-count MPI entry points such as:

- `MPI_Get_count_c`
- `MPI_Isend_c`
- `MPI_Irecv_c`
- `MPI_Send_c`
- `MPI_Rsend_c`
- `MPI_Recv_c`
- `MPI_Sendrecv_c`
- `MPI_Op_create_c`
- `MPI_Type_create_struct_c`
- `MPI_Type_indexed_c`
- `MPI_Type_vector_c`
- `MPI_Bcast_c`
- `MPI_Reduce_c`
- `MPI_Allreduce_c`
- `MPI_Pack_size_c`
- `MPI_Pack_c`
- `MPI_Unpack_c`

If any of these are missing, the current branch is not transport-safe for
general ILP64 use.

## What is in scope now

The current MPI scope is narrow and deliberate.

### In scope

- Detect the actual MPI large-count `_c` capability at configure time.
- Hard-fail ILP64 builds when those APIs are unavailable.
- Keep the BLACS transport on a single code path that assumes verified
  large-count support.
- Document which MPI implementations are currently usable for ILP64.

### Out of scope

- Chunking large messages into multiple legacy `int`-count MPI calls.
- Supporting MPI-3.1-era transports for ILP64 through alternate code paths.
- Adding a local replacement for CMake's `FindMPI.cmake`.
- Accepting an MPI implementation merely because its package release name is
  new enough.

## Why the probe is capability-based

MPI implementation release numbers are not sufficient for this decision.

Examples observed during development:

- Open MPI `4.1.6`
  - reports MPI standard `3.1`
  - does not provide the required `_c` APIs
- Open MPI `5.0.8`
  - still reports MPI standard `3.1`
  - does not provide the required `_c` APIs
- MPICH `4.3.x`
  - reports MPI standard `4.1`
  - does provide the required `_c` APIs

So the branch must test the concrete entry points it uses, not a launcher
version string and not an implementation package version.

## Current acceptance rule

An MPI stack is accepted for ILP64 on this branch only if all of the following
are true:

1. `find_package(MPI)` succeeds for C and Fortran.
2. The Fortran compiler configuration is ILP64-compatible for the build.
3. `CMAKE/check_mpi_large_count.c` compiles and links against `MPI::MPI_C`.
4. BLAS/LAPACK integer ABI is compatible with the requested ILP64 build.

If any of these fail, configure should stop immediately.

## Practical interpretation

This means the branch currently supports:

- large global indices and descriptor arithmetic in `_I8` surfaces
- large MPI counts only on MPI implementations with validated `_c` support

It does **not** mean:

- every `_I8` routine is automatically safe on any MPI-3.1 implementation
- ILP64 can be enabled merely because the MPI package release is recent

## Follow-on work for broader MPI portability

If later we want ILP64 support on MPI implementations without `_c`
large-count entry points, that is a separate phase of work.

That future phase would require:

- auditing every BLACS/MPI boundary for count overflow
- proving each call is bounded by `INT_MAX`, or
- introducing chunked transport helpers, or
- introducing an alternate large-count backend

The highest-priority future audit buckets would be:

- redistribution paths (`PxGEMR2D_I8`, `PxTRMR2D_I8`, `PxLAMR1D_I8`,
  `PxLAMVE_I8`)
- panel-factorization traffic (`PxGETF2_I8`, `PxPOTF2_I8`, support swaps)
- PBLAS reductions and packed collective payloads
- validation/global-check transport helpers

That work is explicitly **not** part of the current branch contract.

## Recommended project wording

For this branch, the supported-policy statement should be:

> ILP64 ScaLAPACK currently requires an MPI implementation whose C interface
> provides the large-count `_c` routines used by BLACS. Builds fail at
> configure time when those APIs are unavailable.

If a later branch adds chunked legacy transport, that statement can be widened.

## Tooling

Use the standalone inspector in:

- `CMAKE/InspectFindMPI/CMakeLists.txt`

to diagnose an MPI installation. It reports:

- CMake's `FindMPI` results
- header/module MPI standard version
- MPI library version string
- the large-count `_c` capability probe log

This inspector is the recommended first step when a new platform is being
qualified for ILP64.
