      SUBROUTINE NARROW_DESC8( DESC8, DESC4 )
      IMPLICIT NONE
*
*  Copy INTEGER*8 descriptor to default INTEGER descriptor.
*  Used by bridge _I8 drivers when calling legacy PBLAS/LAPACK
*  routines that take INTEGER descriptors.
*
      INTEGER*8          DESC8( 9 )
      INTEGER            DESC4( 9 )
      INTEGER            K
      INTRINSIC          INT
*
      DO 10 K = 1, 9
         DESC4( K ) = INT( DESC8( K ) )
   10 CONTINUE
*
      END
