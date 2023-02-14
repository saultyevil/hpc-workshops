#include <stdio.h>
#include <mpi.h>

#define MPI_DEFAULT_TAG 0
#define MPI_DEFAULT_COMM MPI_COMM_WORLD

int
main(int argc, char *argv[])
{

/* -------------------------------------------------------------------------- */

    int n_ranks, current_rank;
    int data_to_send, received_data;
    int rank_sum = 0;
    int right_rank, left_rank;

    MPI_Status recv_status;
    MPI_Request rank_request_right, rank_request_left;

/* -------------------------------------------------------------------------- */

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    if (current_rank == 0)
        printf("Number of processes: %d.\n", n_ranks);
    
    // initialise the start of the sum
    data_to_send = (1 + current_rank) * (1 + current_rank); 
    for (int i = 0; i < n_ranks; i++)
    {
        // figure out the neighbouring ranks
        left_rank = current_rank - 1;
        if (left_rank < 0)
            left_rank = n_ranks - 1;
        right_rank = current_rank + 1;
        if (right_rank > (n_ranks - 1))
            right_rank = 0;

        // send (to right) and receive (to left) -- non-blocking 
        MPI_Issend(&data_to_send, 1, MPI_INT, right_rank, MPI_DEFAULT_TAG,
                   MPI_DEFAULT_COMM, &rank_request_right); 
        MPI_Irecv(&received_data, 1, MPI_INT, left_rank, MPI_DEFAULT_TAG,
                  MPI_DEFAULT_COMM, &rank_request_left);

        // stop the non-blocking send and receives
        MPI_Wait(&rank_request_right, &recv_status);
        MPI_Wait(&rank_request_left, &recv_status);

        // sum everything up
        rank_sum += received_data;
        data_to_send = received_data;
    }

    printf("Rank %d: sum = %d.\n", current_rank, rank_sum);

    MPI_Finalize();
    return 0;
}
        
