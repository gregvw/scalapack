/*
 *  This file includes the standard C libraries, as well as system dependant
 *  include files.  All BLACS routines include this file.
 */

#ifndef BCONFIG_H
#define BCONFIG_H 1

/*
 * Include files
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "scalapack-types.h"
#include <mpi.h>

/*
 * Integer types used by BLACS
 */
#ifndef Int
#define Int ScaLAPACK_ApiInt
#endif

_Static_assert(sizeof(Int) == sizeof(ScaLAPACK_ApiInt),
               "Int must match the configured ScaLAPACK API integer width.");

/*
 * MPI wrapper definitions for ILP64 support
 * Supports two configurations:
 * 1. MPI 4.0+ with _c variants supporting MPI_Count
 * 2. Standard MPI with int counts (MPI < 4.0)
 *
 * Note: BigMPI is NOT supported because it lacks support for custom
 * reduction operations (MPI_Op_create), which are fundamental to ScaLAPACK.
 */

#if MPI_VERSION >= 4
    /* MPI 4.0+ with _c variants supporting MPI_Count */
    #define MpiInt MPI_Count

    static inline int _MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype,
                               MPI_Count *count) {
      return MPI_Get_count_c(status, datatype, count);
    }

    static inline int _MPI_Isend(const void *buf, MPI_Count count, MPI_Datatype datatype,
                          int dest, int tag, MPI_Comm comm, MPI_Request *request) {
      return MPI_Isend_c(buf, count, datatype, dest, tag, comm, request);
    }

    static inline int _MPI_Irecv(void *buf, MPI_Count count, MPI_Datatype datatype,
                          int source, int tag, MPI_Comm comm, MPI_Request *request) {
      return MPI_Irecv_c(buf, count, datatype, source, tag, comm, request);
    }

    static inline int _MPI_Send(const void *buf, MPI_Count count, MPI_Datatype datatype,
                         int dest, int tag, MPI_Comm comm) {
      return MPI_Send_c(buf, count, datatype, dest, tag, comm);
    }

    static inline int _MPI_Recv(void *buf, MPI_Count count, MPI_Datatype datatype,
                         int source, int tag, MPI_Comm comm, MPI_Status *status) {
      return MPI_Recv_c(buf, count, datatype, source, tag, comm, status);
    }

    static inline int _MPI_Op_create(MPI_User_function_c *user_fn, int commute, MPI_Op *op) {
      return MPI_Op_create_c(user_fn, commute, op);
    }

    static inline int _MPI_Type_create_struct(MPI_Count count,
                                        const MPI_Count array_of_blocklengths[],
                                        const MPI_Aint array_of_displacements[],
                                        const MPI_Datatype array_of_types[],
                                        MPI_Datatype *newtype) {
      int ierr;
      MPI_Count *count_displacements;
      size_t i, nitems, alloc_bytes;

      if (count < 0) return MPI_ERR_COUNT;
      nitems = (size_t) count;
      if (!ScaLAPACK_SizeTMul(nitems, sizeof(MPI_Count), &alloc_bytes)) return MPI_ERR_OTHER;
      count_displacements = (MPI_Count *) malloc(alloc_bytes);
      if ((alloc_bytes > 0) && (count_displacements == NULL)) return MPI_ERR_OTHER;
      for (i = 0; i < (size_t) count; ++i) count_displacements[i] = (MPI_Count) array_of_displacements[i];
      ierr = MPI_Type_create_struct_c(count, array_of_blocklengths,
                                       count_displacements,
                                       array_of_types, newtype);
      free(count_displacements);
      return ierr;
    }

    static inline int _MPI_Type_indexed(MPI_Count count,
                                  const MPI_Count array_of_blocklengths[],
                                  const MPI_Count array_of_displacements[],
                                  MPI_Datatype oldtype, MPI_Datatype *newtype) {
      return MPI_Type_indexed_c(count, array_of_blocklengths,
                                 array_of_displacements, oldtype, newtype);
    }

    static inline int _MPI_Type_vector(MPI_Count count, MPI_Count blocklength,
                                 MPI_Count stride, MPI_Datatype oldtype,
                                 MPI_Datatype *newtype) {
      return MPI_Type_vector_c(count, blocklength, stride, oldtype, newtype);
    }

    static inline int _MPI_Bcast(void *buffer, MPI_Count count, MPI_Datatype datatype,
                          int root, MPI_Comm comm) {
      return MPI_Bcast_c(buffer, count, datatype, root, comm);
    }

    static inline int _MPI_Reduce(const void *sendbuf, void *recvbuf, MPI_Count count,
                           MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
      return MPI_Reduce_c(sendbuf, recvbuf, count, datatype, op, root, comm);
    }

    static inline int _MPI_Allreduce(const void *sendbuf, void *recvbuf, MPI_Count count,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
      return MPI_Allreduce_c(sendbuf, recvbuf, count, datatype, op, comm);
    }

    static inline int _MPI_Pack_size(MPI_Count incount, MPI_Datatype datatype,
                               MPI_Comm comm, ScaLAPACK_BufLen *size) {
      MPI_Count mpi_size;
      int ierr = MPI_Pack_size_c(incount, datatype, comm, &mpi_size);
      if (ierr == MPI_SUCCESS) *size = (ScaLAPACK_BufLen) mpi_size;
      return ierr;
    }

    static inline int _MPI_Pack(const void *inbuf, MPI_Count incount, MPI_Datatype datatype,
                          void *outbuf, ScaLAPACK_BufLen outsize,
                          ScaLAPACK_BufLen *position, MPI_Comm comm) {
      MPI_Count mpi_outsize = (MPI_Count) outsize;
      MPI_Count mpi_position = (MPI_Count) (*position);
      int ierr = MPI_Pack_c(inbuf, incount, datatype, outbuf, mpi_outsize, &mpi_position, comm);
      if (ierr == MPI_SUCCESS) *position = (ScaLAPACK_BufLen) mpi_position;
      return ierr;
    }

    static inline int _MPI_Unpack(const void *inbuf, MPI_Count insize,
                            ScaLAPACK_BufLen *position,
                            void *outbuf, MPI_Count outcount, MPI_Datatype datatype,
                            MPI_Comm comm) {
      MPI_Count mpi_position = (MPI_Count) (*position);
      int ierr = MPI_Unpack_c(inbuf, insize, &mpi_position, outbuf, outcount, datatype, comm);
      if (ierr == MPI_SUCCESS) *position = (ScaLAPACK_BufLen) mpi_position;
      return ierr;
    }

  #else /* MPI version < 4 */
    /* Standard MPI with int counts */
    #define MpiInt int

    static inline int _MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype,
                               int *count) {
      return MPI_Get_count(status, datatype, count);
    }

    static inline int _MPI_Isend(const void *buf, int count, MPI_Datatype datatype,
                          int dest, int tag, MPI_Comm comm, MPI_Request *request) {
      return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
    }

    static inline int _MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
                          int source, int tag, MPI_Comm comm, MPI_Request *request) {
      return MPI_Irecv(buf, count, datatype, source, tag, comm, request);
    }

    static inline int _MPI_Send(const void *buf, int count, MPI_Datatype datatype,
                         int dest, int tag, MPI_Comm comm) {
      return MPI_Send(buf, count, datatype, dest, tag, comm);
    }

    static inline int _MPI_Recv(void *buf, int count, MPI_Datatype datatype,
                         int source, int tag, MPI_Comm comm, MPI_Status *status) {
      return MPI_Recv(buf, count, datatype, source, tag, comm, status);
    }

    static inline int _MPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op) {
      return MPI_Op_create(user_fn, commute, op);
    }

    static inline int _MPI_Type_create_struct(int count, const int array_of_blocklengths[],
                                        const MPI_Aint array_of_displacements[],
                                        const MPI_Datatype array_of_types[],
                                        MPI_Datatype *newtype) {
      return MPI_Type_create_struct(count, array_of_blocklengths,
                                     array_of_displacements,
                                     array_of_types, newtype);
    }

    static inline int _MPI_Type_indexed(int count, const int array_of_blocklengths[],
                                  const int array_of_displacements[],
                                  MPI_Datatype oldtype, MPI_Datatype *newtype) {
      return MPI_Type_indexed(count, array_of_blocklengths,
                               array_of_displacements, oldtype, newtype);
    }

    static inline int _MPI_Type_vector(int count, int blocklength,
                                 int stride, MPI_Datatype oldtype,
                                 MPI_Datatype *newtype) {
      return MPI_Type_vector(count, blocklength, stride, oldtype, newtype);
    }

    static inline int _MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
                          int root, MPI_Comm comm) {
      return MPI_Bcast(buffer, count, datatype, root, comm);
    }

    static inline int _MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
                           MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
      return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    }

    static inline int _MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
      return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    }

    static inline int _MPI_Pack_size(int incount, MPI_Datatype datatype,
                               MPI_Comm comm, ScaLAPACK_BufLen *size) {
      int pack_size;
      int ierr = MPI_Pack_size(incount, datatype, comm, &pack_size);
      if (ierr == MPI_SUCCESS) *size = (ScaLAPACK_BufLen) pack_size;
      return ierr;
    }

    static inline int _MPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype,
                          void *outbuf, ScaLAPACK_BufLen outsize,
                          ScaLAPACK_BufLen *position, MPI_Comm comm) {
      int ierr;
      int pack_outsize, pack_position;

      if (outsize > (ScaLAPACK_BufLen) INT_MAX || *position > (ScaLAPACK_BufLen) INT_MAX)
         return MPI_ERR_COUNT;
      pack_outsize = (int) outsize;
      pack_position = (int) (*position);

      ierr = MPI_Pack(inbuf, incount, datatype, outbuf, pack_outsize, &pack_position, comm);
      if (ierr != MPI_SUCCESS) return ierr;
      *position = (ScaLAPACK_BufLen) pack_position;
      return MPI_SUCCESS;
    }

    static inline int _MPI_Unpack(const void *inbuf, int insize,
                            ScaLAPACK_BufLen *position,
                            void *outbuf, int outcount, MPI_Datatype datatype,
                            MPI_Comm comm) {
      int ierr;
      int unpack_position;

      if (*position > (ScaLAPACK_BufLen) INT_MAX) return MPI_ERR_COUNT;
      unpack_position = (int) (*position);

      ierr = MPI_Unpack(inbuf, insize, &unpack_position, outbuf, outcount, datatype, comm);
      if (ierr != MPI_SUCCESS) return ierr;
      *position = (ScaLAPACK_BufLen) unpack_position;
      return MPI_SUCCESS;
    }

  #endif /* MPI version check */

typedef MpiInt ScaLAPACK_MpiCount;


/*
 * These macros define the naming strategy needed for a fortran
 * routine to call a C routine, and whether to build so they may be
 * called from C or fortran.  For the fortran call C interface, ADD_ assumes that
 * fortran calls expect C routines to have an underscore postfixed to the name
 * (Suns, and the Intel expect this).  NOCHANGE indicates that fortran expects
 * the name called by fortran to be identical to that compiled by C
 * (AIX does this).  UPCASE says it expects C routines called by fortran
 * to be in all upcase (CRAY wants this).  The variable FORTRAN_CALL_C is always
 * set to one of these values.  If the BLACS will be called from C, we define
 * INTFACE to be CALL_C, otherwise, it is set to FORTRAN_CALL_C.
 */
#define ADD_     0
#define NOCHANGE 1
#define UPCASE   2
#define FCISF2C  3
#define C_CALL   4

#ifdef UpCase
#define FORTRAN_CALL_C UPCASE
#endif

#ifdef NoChange
#define FORTRAN_CALL_C NOCHANGE
#endif

#ifdef Add_
#define FORTRAN_CALL_C ADD_
#endif

#ifdef FortranIsF2C
#define FORTRAN_CALL_C FCISF2C
#endif

#ifndef FORTRAN_CALL_C
#define FORTRAN_CALL_C ADD_
#endif

#ifdef CallFromC
#define INTFACE C_CALL
#else
#define INTFACE FORTRAN_CALL_C
#endif

/*
 *  Uncomment these macro definitions, and substitute the topology of your
 *  choice to vary the default topology (TOP = ' ') for broadcast and combines.
#define DefBSTop '1'
#define DefCombTop '1'
 */

/*
 * Uncomment this line if your MPI_Send provides a locally-blocking send
 */
//#define SndIsLocBlk

/*
 * Comment out the following line if your MPI does a data copy on every
 * non-contiguous send
 */
#define MpiBuffGood

/*
 * If your MPI cannot form data types of zero length, uncomment the
 * following definition
 */
/* #define ZeroByteTypeBug */

/*
 *  These macros set the timing and debug levels for the BLACS.  The fastest
 *  code is produced by setting both values to 0.  Higher levels provide
 *  more timing/debug information at the cost of performance.  Present levels
 *  of debug are:
 *  0 : No debug information
 *  1 : Mainly parameter checking.
 *
 *  Present levels of timing are:
 *  0 : No timings taken
 */
#ifndef BlacsDebugLvl
#define BlacsDebugLvl 0
#endif
#ifndef BlacsTimingLvl
#define BlacsTimingLvl 0
#endif

#endif
