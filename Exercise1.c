#include <stdio.h>
#include <omp.h>

int main() {

    printf("Amount of thread at THIS time: %d\n", omp_get_num_threads());

#pragma omp parallel
    {
        printf("Amount of thread NOW: %d\n", omp_get_num_threads());
        int ID = omp_get_thread_num();
        printf("hello(%d) ", ID);
        printf("world(%d)\n", ID);
    };

    printf("Amount of thread at the end: %d\n", omp_get_num_threads());

    return 0;
}
