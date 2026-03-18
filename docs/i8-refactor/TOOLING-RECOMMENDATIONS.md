# Tooling Recommendations for the Native I8 PBLAS Refactor

## Purpose

Phase 5 is the first stage where the remaining work is dominated by substantial C code
in `PBLAS/SRC`. The goal of this note is to recommend tooling that reduces manual porting,
improves review quality, and keeps the native-I8 PBLAS effort scoped and testable.

This document is not a replacement for the phase plan. It describes supporting tools and
validation infrastructure that should make the Phase 5 work more disciplined.

## Summary

The recommended tooling stack is:

1. Python scripts for inventory, dependency mapping, and report generation
2. `clang-tidy` for continuous detection of narrowing, cast, and boundary issues
3. Clang AST/Transformer tooling for mechanical C rewrites after the target pattern is settled
4. ILP64 BLAS/LAPACK backends for validation:
   - OpenBLAS with `INTERFACE64=1`
   - Intel MKL ILP64

Recommendation: start with Python plus `clang-tidy`, then add Clang-based rewriting only
after the first native-I8 PBLAS pattern is proven on a small kernel set.

## What the backends do and do not solve

Both OpenBLAS `INTERFACE64=1` and MKL ILP64 are valuable, but only as backend and ABI
validation aids.

They help with:

- 64-bit local BLAS/LAPACK call interfaces
- exercising local dense-kernel call boundaries in ILP64 builds
- catching accidental narrowing in bridge and PBLAS call paths
- providing backend diversity during testing

They do not solve:

- ScaLAPACK descriptor width
- BLACS/MPI count and datatype width
- distributed ownership/index calculations
- native-I8 PBLAS algorithms
- redistribution and distributed-copy semantics

Conclusion: ILP64 BLAS backends are useful for testing the refactor, but they are not a
substitute for native-I8 PBLAS work.

## Recommended validation matrix

Use the same `_I8` tests across multiple backend/build configurations:

- LP64 baseline: current OpenBLAS/reference BLAS build
- ILP64 backend lane 1: OpenBLAS with `INTERFACE64=1`
- ILP64 backend lane 2: MKL ILP64, when available

This gives coverage across two independent ILP64 dense-BLAS implementations and makes it
easier to separate backend ABI issues from ScaLAPACK/PBLAS logic bugs.

## Python tools

Python is the fastest way to get high-leverage visibility before touching the PBLAS C code.
The first scripts should produce reports, not code edits.

### 1. `pblas_callgraph.py`

Purpose:

- build the direct and transitive call graph for a PBLAS seed routine
- scope the first-wave dependency cone for `PxLATRD`
- prevent the Phase 5 slice from drifting into an unplanned horizontal port

Suggested outputs:

- dependency tree for a named routine
- per-routine file path
- call frequency annotations where known
- grouped report of Fortran auxiliaries vs C kernels

Highest-value use:

- `PDLATRD`, `PSLATRD`, `PCLATRD`, `PZLATRD`

### 2. `i8_boundary_audit.py`

Purpose:

- find width-boundary hazards systematically
- give reviewable evidence of what still narrows to default `INTEGER` or `int`

Suggested checks:

- `INT(...)` conversions in Fortran
- `REAL(...)`/`DBLE(...)` used for workspace or count transport
- `MPI_INT` and related MPI datatype choices
- descriptor fields unpacked into `int`
- BLACS handles narrowed without explicit comments or checks
- PBLAS/LAPACK call boundaries still using default-width descriptors
- process-grid or local-size quantities converted without bounds checks

Suggested outputs:

- Markdown or CSV report grouped by file
- severity buckets:
  - definitely unsafe
  - likely narrowing boundary
  - informational

This script would have caught issues like the single-precision workspace path that passed
through `REAL` and `SGAMN2D`.

### 3. `typed_variant_diff.py`

Purpose:

- compare `s/d/c/z` sibling routines after normalizing type-specific tokens
- identify structural divergence that is not explained by expected precision/complex changes

Value:

- makes review of hand-ported `_I8` routines faster
- helps determine whether shared-core PBLAS refactoring is realistic
- flags cases where one type family received a fix that others missed

Good targets:

- `PxAXPY`, `PxSCAL`, `PxDOT`, `PxDOTC`, `PxGEMV`, `PxSYMV`, `PxHEMV`

### 4. `pblas_signature_inventory.py`

Purpose:

- extract a structured inventory of PBLAS entry points and their interfaces

Suggested fields:

- routine name
- source file
- scalar integer arguments
- descriptor arguments
- BLACS calls
- MPI calls
- BLAS/LAPACK calls
- direct PBLAS dependencies
- type family

Value:

- gives a ranked map of which routines are easiest to port first
- exposes which kernels share the same boundary patterns
- supports the duplication-vs-shared-core decision with concrete data

### 5. `distributed_testgen.py`

Purpose:

- generate small MPI regression drivers from templates

Parameters to vary:

- matrix size
- block size
- grid shape
- `RSRC` / `CSRC`
- upper/lower path
- real vs complex path

Value:

- avoids hand-writing a large family of repetitive tests
- makes it easier to add targeted regression cases as native-I8 PBLAS grows

### 6. `i8_port_scaffold.py`

Purpose:

- generate draft `_i8` sibling files or wrapper skeletons

Use carefully:

- suitable for boilerplate generation
- not suitable as a final refactoring engine

Value:

- reduces repetitive setup work
- makes review smaller by concentrating human effort on real width-boundary logic

## Clang-based tooling

Once the initial kernel pattern is stable, Clang tooling becomes more useful than ad hoc
searches because it understands the C syntax and type system.

### `clang-tidy`

Recommended role:

- continuous safety net during the C refactor

Useful checks:

- implicit narrowing conversions
- suspicious casts
- signed/unsigned mismatches
- dead stores exposed by refactoring
- duplicated branches or unreachable code uncovered by porting

Best use in this project:

- run on the first-wave native-I8 PBLAS kernels continuously
- treat warnings on descriptor and local-size narrowing as review blockers

Potential follow-up:

- add a project-specific check for descriptor unpacking or API-width narrowing without an
  explicit guard

### Clang AST tooling

Recommended role:

- semantic inspection and structured reporting

Useful tasks:

- find functions that take descriptor arrays
- find all calls into BLAS/LAPACK/BLACS/MPI
- locate all assignments and casts into `int`
- identify repeated structural patterns across `s/d/c/z` kernels

Best use in this project:

- produce trustworthy inventories and narrowing reports for the PBLAS C layer
- validate that a chosen shared-core pattern is actually consistent across sibling files

### Clang Transformer

Recommended role:

- controlled mechanical rewriting after the target code shape is chosen

Good targets:

- replacing ad hoc `int` parameters with a chosen widened API type
- wrapping BLAS/LAPACK call boundaries in checked narrowing helpers
- replacing repeated descriptor unpack logic with common helpers
- generating consistent LP64/I8 wrapper layers around a shared implementation

Important constraint:

- do not start here first
- transformer rules are most valuable only after the initial hand-ported kernel pair
  establishes the desired structure

## Recommended workflow

### Stage 1: understand and scope

Build first:

- `pblas_callgraph.py`
- `pblas_signature_inventory.py`
- `i8_boundary_audit.py`

Goal:

- make the dependency cone, narrowing sites, and repeated interface patterns explicit before
  editing large amounts of C

### Stage 2: establish the first native-I8 pattern

Recommended first kernels:

- `PxAXPY_I8`
- `PxSCAL_I8`

Support tooling:

- `clang-tidy` on the edited files
- `typed_variant_diff.py` to compare sibling families

Goal:

- decide whether Phase 5 should continue with duplicated `_i8.c` files or move toward a
  shared-core implementation

### Stage 3: automate mechanical work

After the first pattern is stable:

- use Clang AST tooling to verify consistency
- use Clang Transformer for repetitive rewrites where the pattern is now well-defined
- use `distributed_testgen.py` to expand regression coverage as each new native-I8 kernel lands

## Decision guidance: duplication vs shared core

Tooling can help choose between the two main implementation strategies.

Signals that duplication is still the right short-term path:

- `typed_variant_diff.py` shows nontrivial structural divergence across the four type families
- early kernels require many local special cases
- the first ports are still teaching the right boundary pattern

Signals that a shared core is worth pursuing:

- `typed_variant_diff.py` shows strong structural alignment after type-token normalization
- Clang AST reports show near-identical descriptor and call-boundary structure
- two or more first-wave kernels stabilize on the same widening pattern

## Minimal tooling package for Phase 5

If only a small amount of tooling is built before coding resumes, prioritize:

1. `i8_boundary_audit.py`
2. `pblas_signature_inventory.py`
3. `clang-tidy` integration for the first native-I8 PBLAS files

This is the smallest package that meaningfully improves planning, review, and safety.

## Closing recommendation

Use tooling to make the native-I8 PBLAS work measurable and reviewable, not fully automatic.
The best immediate payoff is:

- Python for inventory and audits
- `clang-tidy` for continuous boundary checking
- ILP64 OpenBLAS and MKL as validation backends

Defer Clang Transformer until the first native-I8 kernel pair proves the target code pattern.
