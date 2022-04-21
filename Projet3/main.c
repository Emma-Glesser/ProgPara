#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/*
mpirun  --mca btl  ^openib -np 4 ./main_mpi
*/
int main(int argc, char *argv[])
{
    int rank, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs); // to learn how many proc we are
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // to learn what is my rank

    /* Hello world*/
    printf("rank %d\n",rank);

    MPI_Finalize();
    return 0;
}