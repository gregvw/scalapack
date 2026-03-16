*
*  SL_i8_params.inc -- Shared descriptor parameters for I8 routines
*
*  BLOCK_CYCLIC_2D_I8 = 1 uses the same semantic DTYPE as legacy
*  descriptors (BLOCK_CYCLIC_2D = 1).  Descriptor width is determined
*  by the entry-point name (_I8 suffix), not by the DTYPE value.
*
*  Descriptor index constants (DTYPE_=1 .. LLD_=9) are the same as the
*  legacy layout; only the element type changes (INTEGER*8 vs INTEGER).
*
      INTEGER*8          BLOCK_CYCLIC_2D_I8
      PARAMETER          ( BLOCK_CYCLIC_2D_I8 = 1 )
      INTEGER            CTXT_, DLEN_, DTYPE_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_, CSRC_
      PARAMETER          ( DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5,
     $                     NB_ = 6, RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
