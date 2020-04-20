#include <stdio.h>
#include <omp.h>

int main() {

    printf("Amount of thread at THIS time: %d\n", omp_get_num_threads());

    int id;

    // using "shared" and "default(none)" seems good practice since my compiler was showing some warning
    #pragma omp parallel shared(id) default(none)// private(id)
    {
        printf("Amount of thread NOW: %d\n", omp_get_num_threads());

        id = omp_get_thread_num();
        // usleep(id); // We can try to mess up with the threads by sleeping a different time for each thread
        printf("hello(%d) ", id);

        printf("world(%d)\n", id);
    };

    printf("Amount of thread at the end: %d\n", omp_get_num_threads());

    return 0;
}
