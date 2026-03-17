      SUBROUTINE PZLARFG_I8( N, ALPHA, IAX, JAX, X, IX, JX, DESCX,
     $                       INCX, TAU )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 version of PCLARFG.
*
*  Generates a complex elementary reflector H.
*  Uses PDZNRM2_I8, PZDSCAL_I8, PZSCAL_I8.
*
*     .. Scalar Arguments ..
      INTEGER*8          N, IAX, JAX, IX, JX, INCX
      DOUBLE COMPLEX     ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCX( * )
      DOUBLE COMPLEX     TAU( * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, IXCOL, IXROW, MYCOL, MYROW,
     $                   NPCOL, NPROW, KNT, J4
      INTEGER*8          IIAX, JJAX, INDXTAU, J8
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM,
     $                   CR, CI
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, ZGEBR2D, ZGEBS2D,
     $                   PZSCAL_I8, PZDSCAL_I8, INFOG2L_I8,
     $                   PDZNRM2_I8, DLADIV
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3
      EXTERNAL           DLAMCH, DLAPY3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DIMAG, DCMPLX, DBLE, SIGN, INT
*     ..
*     .. Executable Statements ..
*
      ICTXT = INT( DESCX( CTXT_ ) )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( INCX.EQ.DESCX( M_ ) ) THEN
         CALL INFOG2L_I8( IX, JAX, DESCX, NPROW, NPCOL, MYROW,
     $                    MYCOL, IIAX, JJAX, IXROW, IXCOL )
         IF( MYROW.NE.IXROW ) RETURN
         IF( MYCOL.EQ.IXCOL ) THEN
            J8 = IIAX + (JJAX-1)*DESCX( LLD_ )
            CALL ZGEBS2D( ICTXT, 'Rowwise', ' ', 1, 1, X( J8 ), 1 )
            ALPHA = X( J8 )
         ELSE
            CALL ZGEBR2D( ICTXT, 'Rowwise', ' ', 1, 1, ALPHA, 1,
     $                    MYROW, IXCOL )
         END IF
         INDXTAU = IIAX
      ELSE
         CALL INFOG2L_I8( IAX, JX, DESCX, NPROW, NPCOL, MYROW,
     $                    MYCOL, IIAX, JJAX, IXROW, IXCOL )
         IF( MYCOL.NE.IXCOL ) RETURN
         IF( MYROW.EQ.IXROW ) THEN
            J8 = IIAX + (JJAX-1)*DESCX( LLD_ )
            CALL ZGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1, X( J8 ),
     $                    1 )
            ALPHA = X( J8 )
         ELSE
            CALL ZGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, ALPHA, 1,
     $                    IXROW, MYCOL )
         END IF
         INDXTAU = JJAX
      END IF
*
      IF( N.LE.0 ) THEN
         TAU( INDXTAU ) = ZERO
         RETURN
      END IF
*
      CALL PDZNRM2_I8( N-1, XNORM, X, IX, JX, DESCX, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )
*
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
         TAU( INDXTAU ) = ZERO
      ELSE
         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH( 'S' )
         RSAFMN = ONE / SAFMIN
         IF( ABS( BETA ).LT.SAFMIN ) THEN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL PZDSCAL_I8( N-1, RSAFMN, X, IX, JX, DESCX, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
            CALL PDZNRM2_I8( N-1, XNORM, X, IX, JX, DESCX, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
            TAU( INDXTAU ) = DCMPLX( ( BETA-ALPHR ) / BETA,
     $                              -ALPHI / BETA )
            CALL DLADIV( ONE, ZERO, DBLE( ALPHA-BETA ),
     $              DIMAG( ALPHA-BETA ), CR, CI )
            ALPHA = DCMPLX( CR, CI )
            CALL PZSCAL_I8( N-1, ALPHA, X, IX, JX, DESCX, INCX )
            ALPHA = BETA
            DO 20 J4 = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU( INDXTAU ) = DCMPLX( ( BETA-ALPHR ) / BETA,
     $                              -ALPHI / BETA )
            CALL DLADIV( ONE, ZERO, DBLE( ALPHA-BETA ),
     $              DIMAG( ALPHA-BETA ), CR, CI )
            ALPHA = DCMPLX( CR, CI )
            CALL PZSCAL_I8( N-1, ALPHA, X, IX, JX, DESCX, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of PZLARFG_I8
*
      END
