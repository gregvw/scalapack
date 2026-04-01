#include "scalapack-types.h"

#ifdef T3D
#define float double
#endif
#ifdef T3E
#define float double
#endif
#ifdef CRAY
#define float double
#endif
#ifndef Int
#define Int ScaLAPACK_ApiInt
#endif

_Static_assert(sizeof(Int) == sizeof(ScaLAPACK_ApiInt),
               "Int must match the configured ScaLAPACK API integer width.");
