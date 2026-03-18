      SUBROUTINE PCELGET_I8( SCOPE, TOP, ALPHA, A, IA, JA, DESCA )
      IMPLICIT NONE
*
*  INTEGER*8 version of PDELGET.  Retrieves a single distributed
*  matrix element and broadcasts it within the specified scope.
*
*     .. Scalar Arguments ..
      CHARACTER          SCOPE, TOP
      INTEGER*8          IA, JA
      COMPLEX            ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * )
      COMPLEX            A( * )
*     ..
*
*  =====================================================================
*
      INCLUDE 'SL_i8_params.inc'
      COMPLEX            ZERO
      PARAMETER          ( ZERO = (0.0E+0,0.0E+0) )
*
      INTEGER            IACOL, IAROW, ICTXT, MYCOL, MYROW, NPCOL,
     $                   NPROW
      INTEGER*8          IIA, JJA, IOFFA
*
      EXTERNAL           BLACS_GRIDINFO, CGEBR2D, CGEBS2D, INFOG2L_I8
      LOGICAL            LSAME
      EXTERNAL           LSAME
      INTRINSIC          INT
*
      ICTXT = INT( DESCA( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      CALL INFOG2L_I8( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIA, JJA, IAROW, IACOL )
*
      ALPHA = ZERO
*
      IF( LSAME( SCOPE, 'R' ) ) THEN
         IF( MYROW.EQ.IAROW ) THEN
            IF( MYCOL.EQ.IACOL ) THEN
               IOFFA = IIA + (JJA-1)*DESCA( LLD_ )
               CALL CGEBS2D( ICTXT, SCOPE, TOP, 1, 1, A( IOFFA ), 1 )
               ALPHA = A( IOFFA )
            ELSE
               CALL CGEBR2D( ICTXT, SCOPE, TOP, 1, 1, ALPHA, 1,
     $                       IAROW, IACOL )
            END IF
         END IF
      ELSE IF( LSAME( SCOPE, 'C' ) ) THEN
         IF( MYCOL.EQ.IACOL ) THEN
            IF( MYROW.EQ.IAROW ) THEN
               IOFFA = IIA + (JJA-1)*DESCA( LLD_ )
               CALL CGEBS2D( ICTXT, SCOPE, TOP, 1, 1, A( IOFFA ), 1 )
               ALPHA = A( IOFFA )
            ELSE
               CALL CGEBR2D( ICTXT, SCOPE, TOP, 1, 1, ALPHA, 1,
     $                       IAROW, IACOL )
            END IF
         END IF
      ELSE IF( LSAME( SCOPE, 'A' ) ) THEN
         IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL ) THEN
            IOFFA = IIA + (JJA-1)*DESCA( LLD_ )
            CALL CGEBS2D( ICTXT, SCOPE, TOP, 1, 1, A( IOFFA ), 1 )
            ALPHA = A( IOFFA )
         ELSE
            CALL CGEBR2D( ICTXT, SCOPE, TOP, 1, 1, ALPHA, 1,
     $                    IAROW, IACOL )
         END IF
      ELSE
         IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL )
     $      ALPHA = A( IIA + (JJA-1)*DESCA( LLD_ ) )
      END IF
*
      RETURN
      END
