#include <mpi.h>

static void scalapack_probe_user_fn(void *invec, void *inoutvec,
                                    MPI_Count *len, MPI_Datatype *datatype)
{
  (void)invec;
  (void)inoutvec;
  (void)len;
  (void)datatype;
}

int main(void)
{
  MPI_Count count = 1;
  MPI_Count blocklengths[1] = { 1 };
  MPI_Count displacements[1] = { 0 };
  MPI_Count struct_displacements[1] = { 0 };
  MPI_Datatype datatype = MPI_INT;
  MPI_Datatype newtype = MPI_DATATYPE_NULL;
  MPI_Status status;
  MPI_Request request = MPI_REQUEST_NULL;
  MPI_Op op = MPI_OP_NULL;
  int root = 0;
  int dest = 0;
  int tag = 0;
  void *buffer = 0;

  (void)MPI_Get_count_c(&status, datatype, &count);
  (void)MPI_Isend_c(buffer, count, datatype, dest, tag, MPI_COMM_WORLD, &request);
  (void)MPI_Irecv_c(buffer, count, datatype, dest, tag, MPI_COMM_WORLD, &request);
  (void)MPI_Send_c(buffer, count, datatype, dest, tag, MPI_COMM_WORLD);
  (void)MPI_Rsend_c(buffer, count, datatype, dest, tag, MPI_COMM_WORLD);
  (void)MPI_Recv_c(buffer, count, datatype, dest, tag, MPI_COMM_WORLD, &status);
  (void)MPI_Sendrecv_c(buffer, count, datatype, dest, tag,
                       buffer, count, datatype, dest, tag,
                       MPI_COMM_WORLD, &status);
  (void)MPI_Op_create_c(scalapack_probe_user_fn, 1, &op);
  (void)MPI_Type_create_struct_c(1, blocklengths, struct_displacements,
                                 &datatype, &newtype);
  (void)MPI_Type_indexed_c(1, blocklengths, displacements, datatype, &newtype);
  (void)MPI_Type_vector_c(count, count, count, datatype, &newtype);
  (void)MPI_Bcast_c(buffer, count, datatype, root, MPI_COMM_WORLD);
  (void)MPI_Reduce_c(buffer, buffer, count, datatype, MPI_SUM, root, MPI_COMM_WORLD);
  (void)MPI_Allreduce_c(buffer, buffer, count, datatype, MPI_SUM, MPI_COMM_WORLD);
  (void)MPI_Pack_size_c(count, datatype, MPI_COMM_WORLD, &count);
  (void)MPI_Pack_c(buffer, count, datatype, buffer, count, &count, MPI_COMM_WORLD);
  (void)MPI_Unpack_c(buffer, count, &count, buffer, count, datatype, MPI_COMM_WORLD);

  return 0;
}
