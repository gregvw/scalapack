      SUBROUTINE ZLAMOV_I8( UPLO, M, N, A, LDA, B, LDB )
      IMPLICIT NONE
*
*  INTEGER*8 local matrix copy (DOUBLE COMPLEX).  See DLAMOV_I8 for docs.
*
      CHARACTER          UPLO
      INTEGER*8          M, N, LDA, LDB
      DOUBLE COMPLEX     A( LDA, * ), B( LDB, * )
      INTEGER*8          I, J
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      RETURN
      END
