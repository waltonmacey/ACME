
#include "mpi.h"
#ifndef ZOLTANINTERFACEHPP
#define ZOLTANINTERFACEHPP


#ifdef __cplusplus
extern "C" {
#endif
 void zoltan_partition_problem(
     int *nelem,
     int *xadj,
     int *adjncy,
     double *adjwgt,
     double *vwgt,
     int *nparts,
     MPI_Comm comm,
     double *xcoord,
     double *ycoord,
     double *zcoord,
     int *result_parts,
     int *partmethod);
#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
extern "C" {
#endif
 void zoltan2_print_metrics(
     int *nelem,
     int *xadj,
     int *adjncy,
     double *adjwgt,
     double *vwgt,
     int *nparts,
     MPI_Comm comm,
     int *result_parts);
#ifdef __cplusplus
}
#endif

#endif
