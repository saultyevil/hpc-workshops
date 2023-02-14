#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define TRUE 1
#define FALSE 0

#define MPI_DEFAULT_TAG 0
#define MPI_DEFAULT_COMM MPI_COMM_WORLD

#define NDIMS 2
#define XDIR 1
#define YDIR 0

#define N_NBRS 4
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int print_coords(int *coords, int proc);

int
main(int argc, char *argv[])
{
    int proc, n_procs;
    int data_pass, recv_data, proc_sum = 0;
    int coord_print = FALSE, reorder = 0, disp = 1;
    int nx, ny;

    int coords[NDIMS];
    int dims[NDIMS];
    int dim_period[NDIMS];
    int nbrs[N_NBRS];

    MPI_Status recv_status;
    MPI_Request proc_request;
    MPI_Comm cart_comm;

/* -------------------------------------------------------------------------- */

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    if (argc != 1)
        coord_print = atoi(argv[1]);

    if (proc == 0)
        printf("Number of processes: %d.\n", n_procs);

    /*
     * Initialise the dims array with zeros otherwise MPI_Dims_create will
     * possibly fail when receiving garabge from memory
     */
    for (int i = 0; i < NDIMS; i++)
    {
        dims[i] = 0;
        dim_period[i] = 1;

        if (i == 1)
            dim_period[i] = 0;
    }

    MPI_Dims_create(n_procs, NDIMS, dims);
    MPI_Cart_create(MPI_DEFAULT_COMM, NDIMS, dims, dim_period, reorder,
        &cart_comm);
    MPI_Cart_coords(cart_comm, proc, NDIMS, coords);
    MPI_Comm_rank(cart_comm, &proc);
    if (coord_print == TRUE)
        print_coords(coords, proc);

    nx = dims[0];
    ny = dims[1];

    /*
     * Set the value for each process which will be passed around
     */
    data_pass = proc;
    MPI_Cart_shift(cart_comm, XDIR, disp, &nbrs[LEFT], &nbrs[RIGHT]);
    
    for (int j = 0; j < nx; j++)
    {
        MPI_Issend(&data_pass, 1, MPI_INT, nbrs[RIGHT], MPI_DEFAULT_TAG,
                   cart_comm, &proc_request);
        MPI_Recv(&recv_data, 1, MPI_INT, nbrs[LEFT], MPI_DEFAULT_TAG,
                 cart_comm, &recv_status);
        MPI_Wait(&proc_request, &recv_status);

        proc_sum += recv_data;
        data_pass = recv_data;
    }

    printf("Rank %d: sum = %d.\n", proc, proc_sum);

    MPI_Finalize();
    return 0;
}

int
print_coords(int *coords, int proc)
{
    printf("Rank %d: coords: (", proc);
    for (int k = 0; k < NDIMS; k++)
    {
        printf(" %d ", coords[k]);
    }
    printf(")\n");

    return 0;
}
