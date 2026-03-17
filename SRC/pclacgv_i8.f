      SUBROUTINE PCLACGV_I8( N, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 version of PCLACGV.
*
*  Conjugates the elements of a distributed complex vector.
*
*     .. Scalar Arguments ..
      INTEGER*8          N, IX, JX, INCX
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCX( * )
      COMPLEX            X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
*     ..
*     .. Local Scalars ..
      INTEGER            IXCOL, IXROW, MYCOL, MYROW, NPCOL, NPROW
      INTEGER*8          ICOFFX, IROFFX, IIX, JJX, IOFFX,
     $                   I8, LDX, NP, NQ
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L_I8
*     ..
*     .. External Functions ..
      INTEGER*8          NUMROC_I8
      EXTERNAL           NUMROC_I8
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MOD, INT
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) RETURN
*
      CALL BLACS_GRIDINFO( INT( DESCX( CTXT_ ) ), NPROW, NPCOL,
     $                     MYROW, MYCOL )
*
      CALL INFOG2L_I8( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIX, JJX, IXROW, IXCOL )
*
      LDX = DESCX( LLD_ )
      IF( INCX.EQ.DESCX( M_ ) ) THEN
*
         IF( MYROW.NE.IXROW ) RETURN
         ICOFFX = MOD( JX-1, DESCX( NB_ ) )
         NQ = NUMROC_I8( N+ICOFFX, DESCX( NB_ ), MYCOL, IXCOL,
     $                   NPCOL )
         IF( MYCOL.EQ.IXCOL ) NQ = NQ - ICOFFX
*
         IF( NQ.GT.0 ) THEN
            IOFFX = IIX + (JJX-1)*LDX
            DO 10 I8 = 1, NQ
               X( IOFFX ) = CONJG( X( IOFFX ) )
               IOFFX = IOFFX + LDX
   10       CONTINUE
         END IF
*
      ELSE IF( INCX.EQ.1 ) THEN
*
         IF( MYCOL.NE.IXCOL ) RETURN
         IROFFX = MOD( IX-1, DESCX( MB_ ) )
         NP = NUMROC_I8( N+IROFFX, DESCX( MB_ ), MYROW, IXROW,
     $                   NPROW )
         IF( MYROW.EQ.IXROW ) NP = NP - IROFFX
*
         IF( NP.GT.0 ) THEN
            IOFFX = IIX + (JJX-1)*LDX
            DO 20 I8 = IOFFX, IOFFX+NP-1
               X( I8 ) = CONJG( X( I8 ) )
   20       CONTINUE
         END IF
*
      END IF
*
      RETURN
*
*     End of PCLACGV_I8
*
      END
