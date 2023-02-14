#include <stdlib.h>
#include <stdio.h>

#include "pgmio.h"

int
main(int argc, char *argv[])
{
    int i, j, iter;
    int nx, ny;

    if (argc != 3)
    {
        printf("Incorrect number of CL pars provided. Exepected filename and ");
        printf("n_iterations\n");
        exit(-1);
    }

    char *filename = argv[1];
    int n_iters = atoi(argv[2]);

    pgmsize(filename, &nx, &ny);

    double buff[nx][ny];
    double old[nx+2][ny+2];
    double new[nx+2][ny+2];
    double edge[nx+2][ny+2];

    pgmread(filename, buff, nx, ny);
    printf("FILENAME: %s\nN_ITERS: %d\n", filename, n_iters);

    printf("\nCopying image to buffer...\n");
    for (i = 1; i < nx+1; i++)
    {
        for (j = 1; j < ny+1; j++)
        {
            edge[i][j] = buff[i-1][j-1];
        }
    }

    for (i = 0; i < nx+2; i++)
    {
        for (j = 0; j < ny+2; j++)
        {
            old[i][j] = 255.0;
        }
    }

    printf("\n---- BEGINNING ITERATIONS ----\n\n");

    for (iter = 1; iter <= n_iters; iter++)
    {
        for (i = 1; i < nx+1; i++)
        {
            for (j = 1; j < ny+1; j++)
            {
                new[i][j] = 0.25 * (old[i-1][j] + old[i+1][j] + old[i][j-1]
                    + old[i][j+1] - edge[i][j]);
            }
        }

        for (i = 1; i < nx+1; i++)
        {
            for (j = 1; j < ny+1; j++)
            {
                old[i][j] = new[i][j];
            }
        }

        if (iter % 100 == 0)
            printf("%d iterations complete.\n", iter);
    }

    printf("\n----- END OF ITERATIONS -----\n\n");

    // for (i = 1; i < nx+1; i++)
    //     for (j = 1; j < ny+1; j++)
    //         printf("old[%d, %d] = %f\n", i, j, old[i][j]);

    printf("\nCopying new image to buffer...\n");
    for (i = 1; i < nx+1; i++)
    {
        for (j = 1; j < ny+1; j++)
        {
            buff[i-1][j-1] = old[i][j];
        }
    }

    char *out_filename = "output_image.pgm";
    pgmwrite(out_filename, buff, nx, ny);

    return 0;
}
