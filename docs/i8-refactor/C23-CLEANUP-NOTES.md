# Deferred Notes: C23 Cleanup and Warning Strategy

## Purpose

This note records observations from exploratory Clang warning runs during the native-I8
work. The intent is to preserve a future cleanup plan without pulling the current effort
off the I8-critical path.

This is not part of the active I8 phase scope.

## Current observation

Running Clang with a broad warning set on Linux surfaced a large amount of legacy C code
that is not clean under modern C rules, including C23-oriented diagnostics. This is not
surprising for BLACS/PBLAS-era code, but it is substantial enough that it should be
treated as a distinct modernization track rather than folded into the I8 refactor.

Two conclusions follow:

- New `_i8.c` code should be held to a higher warning standard now.
- Broad legacy C23 cleanup should be planned as a separate mechanical effort later.

## What not to do

- Do not pause the native-I8 work to make the full C codebase C23-clean.
- Do not turn `-Weverything` into a project-wide default.
- Do not mix large-scale warning cleanup with correctness-sensitive I8 ports in the same
  patch series unless the changes are tightly scoped and obviously mechanical.

## Recommended near-term policy

For the active I8 work:

- Keep the new PBLAS `_i8.c` wrappers and helpers clean under a curated Clang warning set.
- Use sanitizers and strict warnings as review tools for new code.
- Treat legacy warning cleanup outside touched files as deferred work.

For the broader codebase:

- Capture warning classes and estimate which are mechanical enough for tooling.
- Tackle only the warning families that can be transformed safely and reviewed in bulk.

## High-concern warning families

The first cleanup candidates should be warning classes that are both high-signal and
structurally repetitive:

1. Non-prototype function declarations and definitions
2. Old-style empty parameter lists
3. Missing prototypes for non-static functions
4. Implicit narrowing and sign-conversion issues at integer boundaries
5. Shadowing in small helper scopes
6. Format-string mismatches

These overlap well with the current I8 correctness concerns and are more valuable than
style-only cleanup.

## Recommended warning bundle

The opt-in Clang warning bundle added for the I8 work is a reasonable starting point for
future C cleanup:

- `-Wall`
- `-Wextra`
- `-Wconversion`
- `-Wsign-conversion`
- `-Wshadow`
- `-Wundef`
- `-Wformat=2`
- `-Wcast-qual`
- `-Wwrite-strings`
- `-Wstrict-prototypes`
- `-Wmissing-prototypes`
- `-Wimplicit-fallthrough`
- `-Wnull-dereference`
- `-Wdouble-promotion`

Recommendation: use this curated set to classify and prioritize work. Do not use
`-Weverything` as the working policy for the full tree.

## Tooling strategy

### Use `clang-tidy` for reporting and gating

`clang-tidy` is a good fit for:

- inventorying warning classes in touched C files
- catching new narrowing and cast problems
- preventing regressions once cleanup starts

It is less suitable as the main rewrite engine for old C code.

### Use Clang Transformer or AST tooling only for narrow classes

Clang-based rewriting is appropriate only when the warning class is:

- repetitive
- structurally uniform
- low-risk to change mechanically

Examples:

- empty-parameter declarations
- non-prototype declarations
- straightforward prototype insertion/update work

Recommendation: use one transformation family at a time, not a broad "make it C23"
rewrite.

### Consider `coccinelle` for legacy C modernization

For repetitive C cleanup in old HPC code, `coccinelle` may be a better fit than
`clang-tidy` for some transformations. It is especially attractive where the source is
uniform but not modernized enough to benefit from more semantic tooling.

This should be evaluated when the cleanup becomes an actual phase.

## Proposed deferred plan

When the I8 work reaches a stable checkpoint:

1. Generate a warning histogram from a Linux Clang build log.
2. Group warnings by family and affected directories.
3. Pick one or two high-value mechanical classes first.
4. Decide per class:
   - `clang-tidy` reporting only
   - Clang Transformer rewrite
   - `coccinelle` rewrite
   - manual cleanup
5. Keep the cleanup phase separate from functional I8 changes.

## Suggested first cleanup phase

If a future modernization phase is started, the best initial target is probably:

1. Prototype hygiene
   - non-prototype declarations
   - empty parameter lists
   - missing prototypes

Then:

2. Integer-boundary warnings in C helper/wrapper layers
   - conversions
   - sign changes
   - format mismatches

This ordering keeps the work mechanical first and correctness-relevant second.

## Exit criterion for a future cleanup phase

A first C23-cleanup phase should aim for:

- all new `_i8.c` files clean under the strict Clang warning bundle
- one or two legacy warning families eliminated project-wide
- no growth in warning count for touched C files
- no mixing of warning-only cleanup with unrelated algorithm changes
