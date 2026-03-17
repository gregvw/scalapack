      SUBROUTINE PSLACPY_I8( UPLO, M, N, A, IA, JA, DESCA,
     $                       B, IB, JB, DESCB )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine --
*     INTEGER*8 version of PSLACPY.
*
*  PSLACPY_I8 copies all or part of a distributed matrix A to another
*  distributed matrix B.  No communication is performed; this is a
*  local copy that decomposes the operation into block-aligned
*  sub-copies handled by PSLACP2_I8.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER*8          IA, IB, JA, JB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER*8          DESCA( * ), DESCB( * )
      REAL               A( * ), B( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INCLUDE 'SL_i8_params.inc'
*     ..
*     .. Local Scalars ..
      INTEGER*8          IAA, IBB, IBLK, IN, ITMP, JAA, JBB,
     $                   JBLK, JN, JTMP, MBA, NBA
*     ..
*     .. External Subroutines ..
      EXTERNAL           PSLACP2_I8
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      MBA = DESCA( MB_ )
      NBA = DESCA( NB_ )
*
      IN = MIN( ((IA+MBA-1)/MBA)*MBA, IA+M-1 )
      JN = MIN( ((JA+NBA-1)/NBA)*NBA, JA+N-1 )
*
      IF( M.LE.( MBA - MOD( IA-1, MBA ) ) .OR.
     $    N.LE.( NBA - MOD( JA-1, NBA ) ) ) THEN
         CALL PSLACP2_I8( UPLO, M, N, A, IA, JA, DESCA,
     $                    B, IB, JB, DESCB )
      ELSE
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL PSLACP2_I8( UPLO, IN-IA+1, N, A, IA, JA, DESCA,
     $                       B, IB, JB, DESCB )
            DO 10 ITMP = IN+1-IA, M-1, MBA
               IBLK = MIN( MBA, M-ITMP )
               IBB = IB + ITMP
               JBB = JB + ITMP
               JAA = JA + ITMP
               CALL PSLACP2_I8( UPLO, IBLK, N-ITMP,
     $                          A, IA+ITMP, JAA, DESCA,
     $                          B, IBB, JBB, DESCB )
   10       CONTINUE
         ELSE IF( LSAME( UPLO, 'L' ) ) THEN
            CALL PSLACP2_I8( UPLO, M, JN-JA+1, A, IA, JA, DESCA,
     $                       B, IB, JB, DESCB )
            DO 20 JTMP = JN+1-JA, N-1, NBA
               JBLK = MIN( NBA, N-JTMP )
               IBB = IB + JTMP
               JBB = JB + JTMP
               IAA = IA + JTMP
               CALL PSLACP2_I8( UPLO, M-JTMP, JBLK,
     $                          A, IAA, JA+JTMP, DESCA,
     $                          B, IBB, JBB, DESCB )
   20       CONTINUE
         ELSE
            IF( M.LE.N ) THEN
               CALL PSLACP2_I8( UPLO, IN-IA+1, N, A, IA, JA, DESCA,
     $                          B, IB, JB, DESCB )
               DO 30 ITMP = IN+1-IA, M-1, MBA
                  IBLK = MIN( MBA, M-ITMP )
                  IBB = IB + ITMP
                  CALL PSLACP2_I8( UPLO, IBLK, N,
     $                             A, IA+ITMP, JA, DESCA,
     $                             B, IBB, JB, DESCB )
   30          CONTINUE
            ELSE
               CALL PSLACP2_I8( UPLO, M, JN-JA+1, A, IA, JA, DESCA,
     $                          B, IB, JB, DESCB )
               DO 40 JTMP = JN+1-JA, N-1, NBA
                  JBLK = MIN( NBA, N-JTMP )
                  JBB = JB + JTMP
                  CALL PSLACP2_I8( UPLO, M, JBLK,
     $                             A, IA, JA+JTMP, DESCA,
     $                             B, IB, JBB, DESCB )
   40          CONTINUE
            END IF
         END IF
*
      END IF
*
      RETURN
*
*     End of PSLACPY_I8
*
      END
