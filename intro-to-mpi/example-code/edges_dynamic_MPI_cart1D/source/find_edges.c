#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "find_edges.h"

int
main(int argc, char *argv[])
{
    int i, j, iter;
    int proc, n_procs;
    int nx, ny, nx_proc, ny_proc;
    int *dims, *dim_period, *nbrs;

    double **buff, **master_buff;
    double **old, **new, **edge;

    MPI_Comm cart_comm;
    MPI_Status recv_status;
    MPI_Request proc_request;

    if (argc != 3)
    {
        printf("Incorrect number of CL pars provided. Exepected filename");
        printf(" and n_iterations\n");
        exit(-1);
    }

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(DEFAULT_COMM, &proc);
    MPI_Comm_size(DEFAULT_COMM, &n_procs);

    /*
     * Read in the filename and get the dimensions
     */
    char *filename = argv[1];
    int n_iters = atoi(argv[2]);
    pgmsize(filename, &nx, &ny);
    master_buff = arralloc(sizeof(*master_buff), 2, nx, ny);

    if (proc == MASTER_PROCESS)
    {
        pgmread(filename, &master_buff[0][0], nx, ny);
        printf("\n----------------------\n");
        printf("FILENAME: %s\nRESOLUTION: %d x %d\nN_ITERS: %d\nN_PROCS: %d\n",
            filename, nx, ny, n_iters, n_procs);
        printf("----------------------\n\n");
    }

    /*
     * Calculate the boundaries for each process.
     * dims is initialised as being 0 otherwise MPI_Dims_create will likely
     * return garbage. The directions are non-periodic. Create a cartersian
     * topology and find the neighbouring processes using MPI_Cart_shift.
     */
    int reorder = 0, disp = 1;
    nbrs = malloc(sizeof(*nbrs) * N_NBRS);
    dims = malloc(sizeof(*dims) * NDIMS);
    dim_period = malloc(sizeof(*dim_period) * NDIMS);

    for (i = 0; i < NDIMS; i++)
    {
        dims[i] = 0;
        dim_period[i] = 0;
    }

    MPI_Dims_create(n_procs, NDIMS, dims);
    MPI_Cart_create(DEFAULT_COMM, NDIMS, dims, dim_period, reorder, &cart_comm);
    MPI_Comm_rank(cart_comm, &proc);
    MPI_Cart_shift(cart_comm, XDIR, disp, &nbrs[LEFT], &nbrs[RIGHT]);

    nx_proc = (int) ceil((double) nx/dims[0]);
    ny_proc = ny;

    /*
     * You're a constant source of disappointment, Jack
     * Allocate memory for all of the arrays -- use arralloc because it keeps
     * array elements contiguous
     */
    buff = arralloc(sizeof(*buff), 2, nx_proc, ny_proc);
    old = arralloc(sizeof(*old), 2, nx_proc+2, ny_proc+2);
    new = arralloc(sizeof(*new), 2, nx_proc+2, ny_proc+2);
    edge = arralloc(sizeof(*edge), 2, nx_proc+2, ny_proc+2);

    /*
     * Use MPI_Scatter to scatter the work across all of the processes
     */
    MPI_Scatter(&master_buff[0][0], nx_proc*ny_proc, MPI_DOUBLE, &buff[0][0],
                nx_proc*ny_proc, MPI_DOUBLE, MASTER_PROCESS, cart_comm);

    /*
     * Copy the buffer array (containing the image) to the edge array and
     * initialise the output image as being completely white -- this acts as
     * an initial guess and sets up halo boundary conditions
     */
    for (i = 1; i < nx_proc+1; i++)
    {
        for (j = 1; j < ny_proc+1; j++)
        {
            edge[i][j] = buff[i-1][j-1];
        }
    }

    for (i = 0; i < nx_proc+2; i++)
    {
        for (j = 0; j < ny_proc+2; j++)
        {
            old[i][j] = 255.0;
        }
    }

    if (proc == MASTER_PROCESS)
        printf("\n---- BEGINNING ITERATIONS ----\n\n");

    for (iter = 1; iter <= n_iters; iter++)
    {
        /*
         * Send and recieve halo cells data from left and right
         * neighbouring processes
         */
        MPI_Issend(&old[nx_proc][1], ny_proc, MPI_DOUBLE, nbrs[RIGHT],
                   DEFAULT_TAG, cart_comm, &proc_request);
        MPI_Recv(&old[0][1], ny_proc, MPI_DOUBLE, nbrs[LEFT], DEFAULT_TAG,
                 cart_comm, &recv_status);
        MPI_Wait(&proc_request, &recv_status);

        MPI_Issend(&old[1][1], ny_proc, MPI_DOUBLE, nbrs[LEFT],
                   DEFAULT_TAG, cart_comm, &proc_request);
        MPI_Recv(&old[nx_proc+1][1], ny_proc, MPI_DOUBLE, nbrs[RIGHT],
                 DEFAULT_TAG, cart_comm, &recv_status);
        MPI_Wait(&proc_request, &recv_status);

        for (i = 1; i < nx_proc+1; i++)
        {
            for (j = 1; j < ny_proc+1; j++)
            {
                new[i][j] = 0.25 * (old[i-1][j] + old[i+1][j] + old[i][j-1]
                    + old[i][j+1] - edge[i][j]);
            }
        }

        for (i = 1; i < nx_proc+1; i++)
        {
            for (j = 1; j < ny_proc+1; j++)
            {
                old[i][j] = new[i][j];
            }
        }

        if (iter % FREQ == 0)
            if (proc == MASTER_PROCESS)
                printf("%d iterations complete.\n", iter);
    }

    if (proc == MASTER_PROCESS)
        printf("\n----- END OF ITERATIONS -----\n\n");

    for (i = 1; i < nx_proc+1; i++)
    {
        for (j = 1; j < ny_proc+1; j++)
        {
            buff[i-1][j-1] = old[i][j];
        }
    }

    free(dims);
    free(dim_period);
    free(nbrs);
    free(old);
    free(new);
    free(edge);

    /*
     * Gather the processes buff and send them to the root
     */
    MPI_Gather(&buff[0][0], nx_proc*ny_proc, MPI_DOUBLE, &master_buff[0][0],
               nx_proc*ny_proc, MPI_DOUBLE, MASTER_PROCESS, cart_comm);

    MPI_Finalize();

    free(buff);

    if (proc == MASTER_PROCESS)
    {
        char *out_filename = "output_image.pgm";
        pgmwrite(out_filename, &master_buff[0][0], nx, ny);
    }

    free(master_buff);

    return 0;
}
