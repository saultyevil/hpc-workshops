#include <stdio.h>
#include <mpi.h>


#define NDIMS 2
#define XDIR 1
#define YDIR 0

#define N_NBRS 4
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int
main(int argc, char *argv[])
{

/* -------------------------------------------------------------------------- */

    int proc, n_procs;
    int reorder = 0, disp = 1;
    int dims[NDIMS], dim_periods[NDIMS], coords[NDIMS];
    int nbrs[N_NBRS];
    
    MPI_Comm cart_comm;

/* -------------------------------------------------------------------------- */

    MPI_Init(NULL, NULL);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    for (int i = 0; i < NDIMS; i++)
    {
        dims[i] = 0;
        dim_periods[i] = 0;
    }

    /*
     * Create the cartesian topology
     */
    MPI_Dims_create(n_procs, NDIMS, dims);
    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, dim_periods, reorder, 
                    &cart_comm);
    MPI_Cart_coords(cart_comm, proc, NDIMS, coords);

    /*
     * Find the up, down, left and right neighbours for each proc
     */
    MPI_Cart_shift(cart_comm, XDIR, disp, &nbrs[LEFT], &nbrs[RIGHT]);
    MPI_Cart_shift(cart_comm, YDIR, disp, &nbrs[UP], &nbrs[DOWN]);

    printf("Rank %d (", proc);
    for (int j = 0; j < NDIMS; j++)
    {
        printf(" %d ", coords[j]);
    }
    printf("): up %d down %d left %d right %d\n",
        nbrs[UP], nbrs[DOWN], nbrs[LEFT], nbrs[RIGHT]);
        
    MPI_Finalize();

    return 0;
}
