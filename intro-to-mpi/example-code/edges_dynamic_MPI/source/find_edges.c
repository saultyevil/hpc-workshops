#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "pgmio.h"
#include "arralloc.h"
#include "find_edges.h"

int
main(int argc, char *argv[])
{
    int i, j, iter;
    int nx, ny, nxP, nyP;
    int current_rank, next_rank, prev_rank, n_ranks;
    MPI_Status recv_status;

    double **buff, **master_buff;
    double **old, **new, **edge;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(DEFAULT_COMM, &current_rank);
    MPI_Comm_size(DEFAULT_COMM, &n_ranks);

    if (argc != 3)
    {
        printf("Incorrect number of CL pars provided. Exepected filename");
        printf(" and n_iterations\n");
        exit(-1);
    }

    char *filename = argv[1];
    int n_iters = atoi(argv[2]);
    pgmsize(filename, &nx, &ny);

    /*
     * Calculate the boundaries for each process
     */
    nxP = nx/n_ranks;
    nyP = ny;

    buff = arralloc(sizeof(*buff), 2, nxP, nyP);
    master_buff = arralloc(sizeof(*master_buff), 2, nx, ny);

    old = arralloc(sizeof(*old), 2, nxP+2, nyP+2);
    new = arralloc(sizeof(*new), 2, nxP+2, nyP+2);
    edge = arralloc(sizeof(*edge), 2, nxP+2, nyP+2);

    if (current_rank == MASTER_PROCESS)
    {
        pgmread(filename, &master_buff[0][0], nx, ny);
        printf("FILENAME: %s\nN_ITERS: %d\nRESOLUTION: %d x %d\nN_PROCS: %d\n",
            filename, n_iters, nx, ny, n_ranks);
    }

    /*
     * Calculate the neighbouring ranks and apply MPI_PROC_NULL if at the edge
     * of the image domain
     */
    if ((next_rank = current_rank + 1) >= n_ranks)
        next_rank = MPI_PROC_NULL;
    if ((prev_rank = current_rank - 1) < 0)
        prev_rank = MPI_PROC_NULL;

    MPI_Scatter(&master_buff[0][0], nxP*nyP, MPI_DOUBLE, &buff[0][0], nxP*nyP,
        MPI_DOUBLE, MASTER_PROCESS, DEFAULT_COMM);

    for (i = 1; i < nxP+1; i++)
    {
        for (j = 1; j < nyP+1; j++)
        {
            edge[i][j] = buff[i-1][j-1];
        }
    }

    for (i = 0; i < nxP+2; i++)
    {
        for (j = 0; j < nyP+2; j++)
        {
            old[i][j] = 255.0;
        }
    }

    if (current_rank == MASTER_PROCESS)
        printf("\n---- BEGINNING ITERATIONS ----\n\n");

    for (iter = 1; iter <= n_iters; iter++)
    {

        MPI_Sendrecv(&old[nxP][1], nyP, MPI_DOUBLE, next_rank, 1,
            &old[0][1], nyP, MPI_DOUBLE, prev_rank, 1, DEFAULT_COMM,
                &recv_status);

        MPI_Sendrecv(&old[1][1], nyP, MPI_DOUBLE, prev_rank, 2,
            &old[nxP+1][1], nyP, MPI_DOUBLE, next_rank, 2, DEFAULT_COMM,
                &recv_status);

        for (i = 1; i < nxP+1; i++)
        {
            for (j = 1; j < nyP+1; j++)
            {
                new[i][j] = 0.25 * (old[i-1][j] + old[i+1][j] + old[i][j-1]
                    + old[i][j+1] - edge[i][j]);
            }
        }

        for (i = 1; i < nxP+1; i++)
        {
            for (j = 1; j < nyP+1; j++)
            {
                old[i][j] = new[i][j];
            }
        }

        if (current_rank == MASTER_PROCESS)
            if (iter % 100 == 0)
                printf("%d iterations complete.\n", iter);
    }

    if (current_rank == MASTER_PROCESS)
        printf("\n----- END OF ITERATIONS -----\n\n");

    for (i = 1; i < nxP+1; i++)
    {
        for (j = 1; j < nyP+1; j++)
        {
            buff[i-1][j-1] = old[i][j];
        }
    }

    MPI_Gather(&buff[0][0], nxP*nyP, MPI_DOUBLE, &master_buff[0][0], nxP*nyP,
        MPI_DOUBLE, MASTER_PROCESS, DEFAULT_COMM);

    if (current_rank == MASTER_PROCESS)
    {
        char *out_filename = "output_image.pgm";
        pgmwrite(out_filename, &master_buff[0][0], nx, ny);
    }

    free(buff);
    free(old);
    free(new);
    free(edge);

    MPI_Finalize();

    return 0;
}
