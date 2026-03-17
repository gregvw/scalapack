      SUBROUTINE PSELSET_I8( A, IA, JA, DESCA, ALPHA )
      IMPLICIT NONE
*
*  INTEGER*8 version of PSELSET.  See PDELSET_I8 for docs.
*
      INTEGER*8          IA, JA
      REAL               ALPHA
      INTEGER*8          DESCA( * )
      REAL               A( * )
*
      INCLUDE 'SL_i8_params.inc'
      INTEGER*8          IIA, JJA
      INTEGER            IACOL, IAROW, MYCOL, MYROW, NPCOL, NPROW
      EXTERNAL           BLACS_GRIDINFO, INFOG2L_I8
      INTRINSIC          INT
*
      CALL BLACS_GRIDINFO( INT( DESCA( CTXT_ ) ), NPROW, NPCOL,
     $                     MYROW, MYCOL )
      CALL INFOG2L_I8( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIA, JJA, IAROW, IACOL )
      IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL )
     $   A( IIA+(JJA-1)*DESCA( LLD_ ) ) = ALPHA
*
      RETURN
      END
