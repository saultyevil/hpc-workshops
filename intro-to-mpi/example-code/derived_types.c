#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define NX 4
#define NY 4
#define MAX_ITER 4
#define DEFAULT_COMM MPI_COMM_WORLD
#define DEFAULT_TAG 0
#define MASTER_PROCESS 0

int
main(int argc, char *argv[])
{
    int proc, n_procs;
    int ColStart = 2;
    double sendmatrix[NX][NY], recvmatrix[NX][NY];

    MPI_Status recv_status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(DEFAULT_COMM, &proc);
    MPI_Comm_size(DEFAULT_COMM, &n_procs);

    /*
     * Exit the program if only one processes is being used
     */
    if (n_procs < 2)
    {
        if (proc == MASTER_PROCESS)
        {
            printf("At least two processes have to be used.\n");
            MPI_Finalize();
            exit(1);
        }
    }

    int cols_to_send = 1;
    MPI_Datatype mcoltype;
    MPI_Type_vector(NX, cols_to_send, NX, MPI_DOUBLE, &mcoltype);
    MPI_Type_commit(&mcoltype);

    /*
     * Initialize matrix with some values
     */
    for(int i = 0;i < NX;i++)
    {
        for(int j= 0;j < NY;j++)
        {
            sendmatrix[i][j] = (i + 1) * pow(10, j) + proc;
        }
    }

    /*
     * Print the matrix before the ping-pongs
     */
    if (proc == MASTER_PROCESS)
    {
        printf("Matrix before ping-pongs:\n\n");
        for (int i = 0; i < NX; i++)
        {
            for (int j = 0; j < NY; j++)
            {
                printf("[%d, %d] = %f\n", i, j, sendmatrix[i][j]);
            }
        }
    }

    for (int iter = 0; iter < MAX_ITER; iter++)
    {
        if (proc == 0)
        {
            for (int ii = 0; ii < NX; ii++)
            {
                for (int jj = ColStart; jj < (ColStart + NY); jj++)
                {
                    MPI_Ssend(&(sendmatrix[0][ColStart-1]), 1, mcoltype, 1,
                              DEFAULT_TAG, DEFAULT_COMM);
                    MPI_Recv(&(recvmatrix[0][ColStart-1]), 1, mcoltype, 1,
                             DEFAULT_TAG, DEFAULT_COMM, &recv_status);
                }
            }
        }
        else if (proc == 1)
        {
            for (int ii = 0; ii < NX; ii++)
            {
                for (int jj = ColStart; jj < (ColStart + NY); jj++)
                {
                    MPI_Recv(&(recvmatrix[0][ColStart-1]), 1, mcoltype, 0,
                             DEFAULT_TAG, DEFAULT_COMM, &recv_status);
                    MPI_Ssend(&(sendmatrix[0][ColStart-1]), 1, mcoltype, 0,
                              DEFAULT_TAG, DEFAULT_COMM);
                }
            }
        }

        for (int iii = 0; iii < NX; iii++)
        {
            for (int jjj = 0; jjj < NY; jjj++)
            {
                sendmatrix[iii][jjj-1] = recvmatrix[iii][jjj-1];
            }
        }

    }

    /*
     * Print the matrix after the ping-pongs
     */
    if (proc == MASTER_PROCESS)
    {
        printf("\nMatrix after ping-pongs:\n\n");
        for (int i = 0; i < NX; i++)
        {
            for (int j = 0; j < NY; j++)
            {
                printf("[%d, %d] = %f\n", i, j, sendmatrix[i][j]);
            }
        }
    }

    MPI_Finalize();
    return 0;
}
