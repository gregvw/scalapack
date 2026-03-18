      SUBROUTINE NARROW_DESC8( DESC8, DESC4 )
      IMPLICIT NONE
*
*  Copy INTEGER*8 descriptor to default INTEGER descriptor.
*  Used by bridge _I8 drivers when calling legacy PBLAS/LAPACK
*  routines that take INTEGER descriptors.
*
*  Aborts via BLACS_ABORT if any descriptor entry exceeds the
*  default INTEGER range.  Entries that are inherently bounded
*  (DTYPE_, RSRC_, CSRC_, MB_, NB_) are still checked for
*  safety, but M_, N_, and LLD_ are the most likely to overflow
*  on large problems.
*
      INTEGER*8          DESC8( 9 )
      INTEGER            DESC4( 9 )
      INTEGER*8          INTMAX
      PARAMETER          ( INTMAX = 2147483647 )
      INTEGER            K, ICTXT
      INTRINSIC          INT, ABS
      EXTERNAL           BLACS_ABORT
*
      DO 10 K = 1, 9
         IF( ABS( DESC8( K ) ).GT.INTMAX ) THEN
            ICTXT = INT( DESC8( 2 ) )
            CALL BLACS_ABORT( ICTXT, 1 )
         END IF
         DESC4( K ) = INT( DESC8( K ) )
   10 CONTINUE
*
      END
