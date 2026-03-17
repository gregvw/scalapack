# Resolving the Descriptor ABI Strategy

The key is to stop treating descriptor DTYPE as if it needs to encode both layout semantics and integer width.

  Recommendation
  Use a two-axis model:

  - DTYPE encodes only descriptor semantics/layout.
  - Routine name and argument type encode integer width.

  In this repo, the semantic tags already exist:

  - 1 = 2D block-cyclic
  - 501 = legacy 1D horizontal
  - 502 = legacy 1D vertical

  That is already baked into code like desc_convert.f:14 and many solver paths such as pcptsv.f:407. So the current phase-3 choice of
  BLOCK_CYCLIC_2D_I8 = 501 in SL_i8_params.inc:4 is the wrong abstraction.

  What I would do

  1. Change the I8 2D descriptor tag back to semantic 1.
  2. Treat _I8 entry points as the width distinction.
  3. Add typed helper families for the Fortran boundary:
      - CHK1MAT_I8
      - PCHK1MAT_I8
      - PCHK2MAT_I8
      - DESC_CONVERT_I8 if you need 1D conversions on I8 descriptors
  4. Keep the internal C core width-agnostic after unpacking into widened structs.

  That gives you this rule:

  - LP64 path: PDGEMR2D(..., INTEGER DESCA(9), ...), DESCA(DTYPE_)=1
  - I8 path: PDGEMR2D_I8(..., INTEGER*8 DESCA(9), ...), DESCA(DTYPE_)=1

  Same semantic descriptor, different ABI surface.

  Why this is better

  - It removes the 501 collision immediately.
  - It avoids auditing every existing DTYPE==501/502 branch just to represent width.
  - It matches how the code is already organized: width-specific entry points already exist in pdgemr.c:1050 and pdtrmr.c:972.
  - It scales cleanly when you add I8 1D descriptors later: 501/502 can stay 1D semantic tags for both LP64 and I8 helper families.

  Decision point
  You should explicitly choose between two project goals, because the descriptor strategy differs:

  1. Full ILP64 build
      - One ABI per build.
      - No public _I8 tags needed.
      - DTYPE values stay exactly as today.
  2. Dual surface in one build
      - Legacy LP64 plus selected _I8 entry points.
      - This is what your current branch is doing.
      - In this model, width must live in symbol names and typed helper routines, not in DTYPE.

  Given your current branch, I would stay with option 2.

  Concrete migration

  - Replace BLOCK_CYCLIC_2D_I8=501 with semantic BLOCK_CYCLIC_2D=1 in the new I8 tool layer.
  - Make unpack_desc_i8 accept semantic tags only.
  - Update DESCINIT_I8 and DESCSET_I8 to emit DTYPE=1.
  - Add _I8 validation/conversion helpers before expanding into more Fortran routines.
  - Only after that, start converting the first redistribution-adjacent Fortran callers.



• API Shape

  - Keep descriptor semantics unchanged across widths:
      - 1 = 2D block-cyclic
      - 501 = 1D horizontal
      - 502 = 1D vertical
  - Put width in the symbol, not in DTYPE.
      - DESCINIT / CHK1MAT operate on default INTEGER
      - DESCINIT_I8 / CHK1MAT_I8 operate on INTEGER*8
  - For _I8 entry points, make all public integer arguments INTEGER*8, including INFO. That keeps call sites simple and avoids mixed-kind wrappers.
    Narrow only internally where the value is semantically a BLACS handle, process coordinate, or parameter index.
  - Remove the current BLOCK_CYCLIC_2D_I8 = 501 idea from SL_i8_params.inc:1. The I8 2D descriptor should still carry DTYPE=1, otherwise it conflicts
    with existing 1D logic in desc_convert.f:1 and routines like pcptsv.f:407.

  Helper Set

  - Add I8 twins of the existing boundary helpers:
      - CHK1MAT_I8
      - PCHK1MAT_I8
      - PCHK2MAT_I8
      - DESC_CONVERT_I8
  - Keep their semantics identical to the existing routines in chk1mat.f:1, pchkxmat.f:1, and desc_convert.f:1. The only change should be integer
    width plus explicit narrowing checks at BLACS/MPI boundaries.
  - Add one internal 64-bit comparison/reduction helper for distributed validation.
      - PCHK1MAT_I8 and PCHK2MAT_I8 should not pack values into default INTEGER work arrays the way the current code does.
      - Implement either GLOBCHK_I8 or a small C-backed helper that reduces packed int64_t (value, position) pairs across the BLACS context.

  Rollout

  1. Normalize descriptor tags.
      - Change descinit_i8.f:1 and descset_i8.f:1 so 2D I8 descriptors emit DTYPE=1.
      - Preserve 501/502 only for actual 1D descriptors in both LP64 and I8 helper families.
  2. Build the I8 boundary layer in TOOLS.
      - CHK1MAT_I8 first.
      - Then DESC_CONVERT_I8.
      - Then PCHK1MAT_I8 / PCHK2MAT_I8 once the 64-bit global-check transport is in place.
  3. Test the helper layer before touching major drivers.
      - Wire test_i8_tools.f:1 into CMake.
      - Add focused tests for:
          - valid and invalid I8 2D descriptors
          - I8 2D to I8 1D conversion
          - overflow rejection for CTXT, RSRC, CSRC
          - distributed mismatch detection in PCHK1MAT_I8
  4. Migrate first consumers in a narrow band.
      - Pick small Fortran routines that already depend on PxGEMR2D or PxTRMR2D and have limited workspace logic.
      - Do not start with p?gesvd or eigensolver drivers; they are too entangled.
  5. Only then expand routine-by-routine.
      - Any _I8 routine should use only _I8 helpers.
      - No routine should infer descriptor width from DTYPE.

  Practical rule
  A routine should know descriptor storage width from its own entry point, never from the descriptor contents. DTYPE tells you layout. The symbol
  name tells you ABI.


