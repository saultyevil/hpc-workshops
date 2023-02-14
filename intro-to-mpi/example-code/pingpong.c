#include <stdio.h>
#include <mpi.h>

#define MPI_DEFAULT_TAG 0
#define MPI_DEFAULT_COMM MPI_COMM_WORLD

int
main(int argc, char *argv[])
{
    int current_rank, n_ranks;
    int sent_msg = 42;
    int recv_msg;

    MPI_Status recv_status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_DEFAULT_COMM, &current_rank);
    MPI_Comm_size(MPI_DEFAULT_COMM, &n_ranks);

    if (current_rank == 0)
        printf("Number of processes: %d.\n", n_ranks);

    
    if (current_rank == 0)
    {
        MPI_Ssend(&sent_msg, 1, MPI_INT, 1, MPI_DEFAULT_TAG, MPI_DEFAULT_COMM);
        printf("Rank %d: message = %d.\n", current_rank, sent_msg);
        MPI_Recv(&sent_msg, 1, MPI_INT, 1, MPI_DEFAULT_TAG, MPI_DEFAULT_COMM, 
                &recv_status);
        printf("Rank %d: message = %d.\n", current_rank, sent_msg);
    }
    else if (current_rank == 1)
    {
        MPI_Recv(&recv_msg, 1, MPI_INT, 0, MPI_DEFAULT_TAG, MPI_DEFAULT_COMM,
                &recv_status);
        printf("Rank %d: message = %d.\n", current_rank, recv_msg);
        MPI_Ssend(&recv_msg, 1, MPI_INT, 0, MPI_DEFAULT_TAG, MPI_DEFAULT_COMM);
    }

    MPI_Finalize();
    return 0;
}

