#include <stdio.h>
#include <mpi.h>

#define MPI_DEFAULT_TAG 0
#define MPI_DEFAULT_COMM MPI_COMM_WORLD

int
main(int argc, char *argv[])
{
    int n_ranks, current_rank;
    int data_to_send;
    int rank_sum = 0;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    if (current_rank == 0)
        printf("Number of processes: %d.\n", n_ranks);

    // initialise the start of the sum
    data_to_send = current_rank;
    MPI_Allreduce(&data_to_send, &rank_sum, 1, MPI_INT, MPI_SUM,
                  MPI_DEFAULT_COMM);

    printf("Rank %d: sum = %d.\n", current_rank, rank_sum);

    MPI_Finalize();
    return 0;
}

