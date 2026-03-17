      SUBROUTINE PSLARFG_I8( N, ALPHA, IAX, JAX, X, IX, JX, DESCX,
     $                       INCX, TAU )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 version of PDLARFG.
*
*  Generates a real elementary reflector H of order n, such that
*  H * sub( X ) = ( alpha, 0 )'.  Uses PSNRM2_I8 and PSSCAL_I8.
*
*     .. Scalar Arguments ..
      INTEGER*8          N, IAX, JAX, IX, JX, INCX
      REAL               ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCX( * )
      REAL               TAU( * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, IXCOL, IXROW, MYCOL, MYROW,
     $                   NPCOL, NPROW, KNT, J4
      INTEGER*8          IIAX, JJAX, INDXTAU, J8
      REAL               BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, SGEBR2D, SGEBS2D,
     $                   PSSCAL_I8, INFOG2L_I8, PSNRM2_I8
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLAPY2
      EXTERNAL           SLAMCH, SLAPY2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, INT
*     ..
*     .. Executable Statements ..
*
      ICTXT = INT( DESCX( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( INCX.EQ.DESCX( M_ ) ) THEN
*
*        sub( X ) is distributed across a process row.
*
         CALL INFOG2L_I8( IX, JAX, DESCX, NPROW, NPCOL, MYROW,
     $                    MYCOL, IIAX, JJAX, IXROW, IXCOL )
*
         IF( MYROW.NE.IXROW )
     $      RETURN
*
         IF( MYCOL.EQ.IXCOL ) THEN
            J8 = IIAX + (JJAX-1)*DESCX( LLD_ )
            CALL SGEBS2D( ICTXT, 'Rowwise', ' ', 1, 1, X( J8 ), 1 )
            ALPHA = X( J8 )
         ELSE
            CALL SGEBR2D( ICTXT, 'Rowwise', ' ', 1, 1, ALPHA, 1,
     $                    MYROW, IXCOL )
         END IF
*
         INDXTAU = IIAX
*
      ELSE
*
*        sub( X ) is distributed across a process column.
*
         CALL INFOG2L_I8( IAX, JX, DESCX, NPROW, NPCOL, MYROW,
     $                    MYCOL, IIAX, JJAX, IXROW, IXCOL )
*
         IF( MYCOL.NE.IXCOL )
     $      RETURN
*
         IF( MYROW.EQ.IXROW ) THEN
            J8 = IIAX + (JJAX-1)*DESCX( LLD_ )
            CALL SGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1, X( J8 ),
     $                    1 )
            ALPHA = X( J8 )
         ELSE
            CALL SGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, ALPHA, 1,
     $                    IXROW, MYCOL )
         END IF
*
         INDXTAU = JJAX
*
      END IF
*
      IF( N.LE.0 ) THEN
         TAU( INDXTAU ) = ZERO
         RETURN
      END IF
*
      CALL PSNRM2_I8( N-1, XNORM, X, IX, JX, DESCX, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
         TAU( INDXTAU ) = ZERO
*
      ELSE
*
         BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = SLAMCH( 'S' )
         RSAFMN = ONE / SAFMIN
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL PSSCAL_I8( N-1, RSAFMN, X, IX, JX, DESCX, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
            CALL PSNRM2_I8( N-1, XNORM, X, IX, JX, DESCX, INCX )
            BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
            TAU( INDXTAU ) = ( BETA-ALPHA ) / BETA
            CALL PSSCAL_I8( N-1, ONE/(ALPHA-BETA), X, IX, JX, DESCX,
     $                      INCX )
*
            ALPHA = BETA
            DO 20 J4 = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU( INDXTAU ) = ( BETA-ALPHA ) / BETA
            CALL PSSCAL_I8( N-1, ONE/(ALPHA-BETA), X, IX, JX, DESCX,
     $                      INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of PSLARFG_I8
*
      END
