#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define MPI_DEFAULT_TAG 0
#define N 1000000

double pi_sum(int lower, int upper);

int
main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int current_rank, n_ranks, process_lower, process_upper, chunksize;
    double partial_sum, partial_sums, pi, pi_error;
    MPI_Status recv_status;

    MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    if (current_rank == 0)
    {
        printf("Number of processes used:  %d.\n", n_ranks);
    }

    chunksize = (int) ceil((double) N/n_ranks);
    process_lower = current_rank * chunksize;
    process_upper = (current_rank + 1) * chunksize;
    if (process_upper > N)
        process_upper = N;
  
    partial_sum = pi_sum(process_lower, process_upper);

    if (current_rank == 0)
    {
        pi = partial_sum;
        
        // rank 0 gets things from all other ranks
        for (int rank = 1; rank < n_ranks; rank++)
        {
            MPI_Recv(&partial_sums, 1, MPI_DOUBLE, rank, MPI_DEFAULT_TAG, 
                     MPI_COMM_WORLD, &recv_status);
            pi += partial_sums;
        }

        pi *= 4.0/N;
        pi_error = pi - M_PI;
        
        printf("Pi = %f: error = %f.\n", pi, pi_error); 
    }
    else
    {
        // send partial sum to rank 0
        MPI_Ssend(&partial_sum, 1, MPI_DOUBLE, 0, MPI_DEFAULT_TAG,
                  MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    return 0;
}

double
pi_sum(int lower, int upper)
{

    double sum_value = 0;

    for (int i = lower; i < upper; i++)
    {
        sum_value += 1.0/(1 + ((i - 0.5)/N) * ((i - 0.5)/N));
    }

    return sum_value;
}

