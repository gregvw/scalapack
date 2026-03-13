      SUBROUTINE INFOG2L_I8( GRINDX, GCINDX, DESC, NPROW, NPCOL,
     $                       MYROW, MYCOL, LRINDX, LCINDX, RSRC,
     $                       CSRC )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     .. Scalar Arguments ..
      INTEGER*8          GRINDX, GCINDX, LRINDX, LCINDX
      INTEGER            NPROW, NPCOL, MYROW, MYCOL, RSRC, CSRC
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESC( * )
*     ..
*
*  Purpose
*  =======
*
*  INFOG2L_I8 computes the starting local indexes LRINDX, LCINDX
*  corresponding to the distributed submatrix starting globally at
*  the entry pointed by GRINDX, GCINDX.  This routine returns the
*  coordinates in the grid of the process owning the matrix entry of
*  global indexes GRINDX, GCINDX, namely RSRC and CSRC.
*
*  This is the INTEGER*8 version of INFOG2L, supporting global indices
*  and local indices larger than 2^31-1.
*
*  Arguments
*  =========
*
*  GRINDX    (global input) INTEGER*8
*            The global row starting index of the submatrix.
*
*  GCINDX    (global input) INTEGER*8
*            The global column starting index of the submatrix.
*
*  DESC      (input) INTEGER*8 array of dimension DLEN_.
*            The array descriptor for the underlying distributed matrix.
*
*  NPROW     (global input) INTEGER
*            The total number of process rows over which the distributed
*            matrix is distributed.
*
*  NPCOL     (global input) INTEGER
*            The total number of process columns over which the
*            distributed matrix is distributed.
*
*  MYROW     (local input) INTEGER
*            The row coordinate of the process calling this routine.
*
*  MYCOL     (local input) INTEGER
*            The column coordinate of the process calling this routine.
*
*  LRINDX    (local output) INTEGER*8
*            The local rows starting index of the submatrix.
*
*  LCINDX    (local output) INTEGER*8
*            The local columns starting index of the submatrix.
*
*  RSRC      (global output) INTEGER
*            The row coordinate of the process that possesses the first
*            row and column of the submatrix.
*
*  CSRC      (global output) INTEGER
*            The column coordinate of the process that possesses the
*            first row and column of the submatrix.
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
*     ..
*     .. Local Scalars ..
      INTEGER*8          CBLK, GCCPY, GRCPY, RBLK, NPR8, NPC8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MOD
*     ..
*     .. Executable Statements ..
*
*     Promote grid dimensions to INTEGER*8 for safe mixed-kind arithmetic
*
      NPR8 = NPROW
      NPC8 = NPCOL
*
      GRCPY = GRINDX - 1
      GCCPY = GCINDX - 1
*
      RBLK = GRCPY / DESC(MB_)
      CBLK = GCCPY / DESC(NB_)
*     INT() narrowing is safe: MOD result is in [0, NPROW) or [0, NPCOL),
*     bounded by process grid dimensions, not by matrix dimensions.
      RSRC = INT( MOD( RBLK + DESC(RSRC_), NPR8 ) )
      CSRC = INT( MOD( CBLK + DESC(CSRC_), NPC8 ) )
*
      LRINDX = ( RBLK / NPR8 + 1 ) * DESC(MB_) + 1
      LCINDX = ( CBLK / NPC8 + 1 ) * DESC(NB_) + 1
*
      IF( MOD( MYROW+NPR8-DESC(RSRC_), NPR8 ) .GE.
     $    MOD( RBLK, NPR8 ) ) THEN
         IF( MYROW.EQ.RSRC )
     $      LRINDX = LRINDX + MOD( GRCPY, DESC(MB_) )
         LRINDX = LRINDX - DESC(MB_)
      END IF
*
      IF( MOD( MYCOL+NPC8-DESC(CSRC_), NPC8 ) .GE.
     $    MOD( CBLK, NPC8 ) ) THEN
         IF( MYCOL.EQ.CSRC )
     $      LCINDX = LCINDX + MOD( GCCPY, DESC(NB_) )
         LCINDX = LCINDX - DESC(NB_)
      END IF
*
      RETURN
*
*     End of INFOG2L_I8
*
      END
