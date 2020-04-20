#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

#define N 1000000
#define T 30
#define OUTPUT 6

// We have to declare x outside main otherwise this makes the pile of execution explode
double x[N]; // Put in the tas


// Copied from stream.c
double mysecond()
{
    struct timeval tp;
    struct timezone tzp;

    gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

int main() {

    int thN = omp_get_max_threads();
    printf("#Thread = %d\n\n", thN);

    double avgTime[OUTPUT] = {0};
    double minTime[OUTPUT] = {[0 ... 4] = 20};
    double maxTime[OUTPUT] = {0, 0, 0, 0, 0};

    double maxval[OUTPUT] = {0.0};
    int maxloc[OUTPUT] = {0};

    double t[OUTPUT+1] = {0};
    double time[OUTPUT] = {0};

    char* label[OUTPUT] = {"Serial", "//", "Crit", "Avoid", "Padding", "loc Str"};

    int thLoc[thN];
    double thVal[thN];

    // Avoid False sharing with padding
    int pad = 128;
    int thLocPadded[thN*pad];
    double thValPadded[thN*pad];

    // Avoid False sharing with local storage
    double maxVlocal[thN];
    int maxLlocal[thN];

    // Fancy progress bar
    for (int i = 0; i<T; i++) printf("|");
    printf("\n");

    for (int k = 0; k < T; k++) {

        usleep(20000); // Let some temporal space between tries (20ms) in case the system is busy or something

        srand((unsigned int) mysecond()*1000); // seed
        for (int i = 0; i < N; i++) {
            // Generate random number between 0 and 1
            x[i] = ((double) (rand()) / RAND_MAX) * ((double) (rand()) / RAND_MAX) * ((double) (rand()) / RAND_MAX) * 1000;
        }

        t[0] = mysecond();

        // === Serial ===

        for (int i = 0; i < N; i++) {
            if (x[i] > maxval[0]) {
                maxval[0] = x[i];
                maxloc[0] = i;
            }
        }

        t[1] = mysecond();

        // === Parallel ===

        // using "shared" and "default(none)" seems good practice since my compiler was showing some warning
        // This should lead to some error due to writting at the same time a variable... Can't see any on my 4 threads
        #pragma omp parallel for shared(x, maxval, maxloc) default(none)
        for (int i = 0; i < N; i++) {

            if (x[i] > maxval[1]) {
                maxval[1] = x[i];
                maxloc[1] = i;
            }
        }

        t[2] = mysecond();

        // === Crit ===

        #pragma omp parallel for shared(x, maxval, maxloc) default(none)
            // This should be way more slow than before
            for (int i = 0; i < N; i++) {
                if (x[i] > maxval[2]) // This if statement avoid staying blocked waiting for the critical part to be free
                // Indeed if x[i] is lower than the maximum, it doesn't need to modify the maximum
                // Therefore we don't wait if we already know that we don't need to rewrite the maximum and this improve by quite a lot
                {
                    #pragma omp critical
                    {
                        if (x[i] > maxval[2]) {
                            maxval[2] = x[i];
                            maxloc[2] = i;
                        }
                    }
                }
            }

        t[3] = mysecond();

        // === Avoid ===

        #pragma omp parallel shared(x, thVal, thLoc) default(none)
        {
            int id = omp_get_thread_num();
            // Without initialization it works better than with ...
            // thVal[id] = -1e30;
            #pragma omp for
            for (int i = 0; i < N; i++) {
                if (x[i] > thVal[id]) {
                    thVal[id] = x[i];
                    thLoc[id] = i;
                }
            }
        }

        // Single thread part
        maxval[3] = thVal[0];
        maxloc[3] = thLoc[0];
        for(int j=1; j<thN; j++){
            if(thVal[j] > maxval[3]){
                maxval[3] = thVal[j];
                maxloc[3] = thLoc[j];
            }
        }

        t[4] = mysecond();


        // === W/o FS : padding ===

        #pragma omp parallel shared(x, thValPadded, thLocPadded, pad) default(none)
        {
            int id = omp_get_thread_num()*pad;
            // Without initialization it works better than with ...
            // thValPadded[id] = -1e30;

            #pragma omp for
            for (int i = 0; i < N; i++) {
                if (x[i] > thValPadded[id]) {
                    thValPadded[id] = x[i];
                    thLocPadded[id] = i;
                }
            }

        }

        // Single thread part
        maxval[4] = thValPadded[0];
        maxloc[4] = thLocPadded[0];
        for(int j=1; j<thN; j++){
            if(thValPadded[j*pad] > maxval[4]){
                maxval[4] = thValPadded[j*pad];
                maxloc[4] = thLocPadded[j*pad];
            }
        }


        t[5] = mysecond();

        // === W/o FS : Local store ===
        // Sometimes this one gives the wrong answer for no reason found

        #pragma omp parallel shared(x, maxVlocal, maxLlocal) default(none)
        {
            int id = omp_get_thread_num();
            double thValMaxLocal = -1e30;
            int thLocMaxLocal = 0;

            #pragma omp for
            for (int i = 0; i < N; i++) {
                if (x[i] > thValMaxLocal) {
                    thValMaxLocal = x[i];
                    thLocMaxLocal = i;
                }
            }

            maxVlocal[id] = thValMaxLocal;
            maxLlocal[id] = thLocMaxLocal;
        }

        // Single thread part
        maxval[5] = maxVlocal[0];
        maxloc[5] = maxLlocal[0];
        for(int j=1; j<thN; j++){
            if(thVal[j] > maxval[5]){
                maxval[5] = thVal[j];
                maxloc[5] = thLoc[j];
            }
        }

        t[6] = mysecond();



        // Just compute the time stats
        for(int m = 0; m < OUTPUT; m++){
            time[m] = t[m+1] - t[m];
            avgTime[m] += time[m]/T;

            if(time[m] > maxTime[m]) maxTime[m] = time[m];
            else if (time[m] < minTime[m]) minTime[m] = time[m];

            // Print if the parallel version get a wrong answer
            if(maxloc[0] != maxloc[m]) printf("maxVal:%f!=%f [%d]\n", maxval[0], maxval[m], m);
        }

        // Fancy progress bar
        printf(".");

    }

    printf("\n\nLabel\t|Avg(s)\t\t|Min(s)\t\t|Max(s)\n");

    for(int m = 0; m < OUTPUT; m++){
        printf("%s\t|%f\t|%f\t|%f\n", label[m], avgTime[m], minTime[m], maxTime[m]);
    }
    return 0;
}

