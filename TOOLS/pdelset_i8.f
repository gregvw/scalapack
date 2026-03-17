      SUBROUTINE PDELSET_I8( A, IA, JA, DESCA, ALPHA )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine --
*     INTEGER*8 version of PDELSET.
*
*  PDELSET_I8 sets the distributed matrix entry A( IA, JA ) to ALPHA.
*
*     .. Scalar Arguments ..
      INTEGER*8          IA, JA
      DOUBLE PRECISION   ALPHA
*     ..
*     .. Array arguments ..
      INTEGER*8          DESCA( * )
      DOUBLE PRECISION   A( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
*     ..
*     .. Local Scalars ..
      INTEGER*8          IIA, JJA
      INTEGER            IACOL, IAROW, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L_I8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( INT( DESCA( CTXT_ ) ), NPROW, NPCOL,
     $                     MYROW, MYCOL )
*
      CALL INFOG2L_I8( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIA, JJA, IAROW, IACOL )
*
      IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL )
     $   A( IIA+(JJA-1)*DESCA( LLD_ ) ) = ALPHA
*
      RETURN
*
*     End of PDELSET_I8
*
      END
