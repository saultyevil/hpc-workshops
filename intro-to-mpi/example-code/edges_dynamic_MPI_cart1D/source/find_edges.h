void pgmsize (char *filename, int *nx, int *ny);
void pgmread (char *filename, void *vx, int nx, int ny);
void pgmwrite(char *filename, void *vx, int nx, int ny);

void *arralloc(size_t size, int ndim, ...);

#define DEFAULT_COMM MPI_COMM_WORLD
#define DEFAULT_TAG 0
#define MASTER_PROCESS 0

#define FREQ 100

#define NDIMS 1
#define XDIR 0

#define N_NBRS 4
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3
