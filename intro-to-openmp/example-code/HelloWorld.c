#include <omp.h>
#include <stdio.h>

int main(void)
{
    #pragma omp parallel
    {
        printf("Hello from thread: %i\n", omp_get_thread_num());
    }

}
