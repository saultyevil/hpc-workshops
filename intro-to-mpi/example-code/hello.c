#include <stdio.h>
#include <mpi.h>

int
main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int rank, size, namelen;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Get_processor_name(hostname, &namelen);

	printf("Hello world from rank %d with hostname %s!", rank, hostname);
 	printf("There are %d processes.\n", size);

	if (rank == 0)
	{
		printf("This is a special message from the master process :-).\n");
	}

	MPI_Finalize();

	return 0;
}
