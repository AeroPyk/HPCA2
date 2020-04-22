#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define N 500 // number of particle
#define DELTA_T 0.05
#define CYC 100
#define DIM 2
#define G 6.6742e-4 // instead of e-11
// Have to put a higher.
// It's either that or we get masses way higher because this force is mostly applicable to big object


void init();
void sA();
void rA();
void mover(double delta_t);
void printBodies();
void doubleCheck();
double round(double x);

void reducedAlgo(int q);

void simpleAlgo(int q);

typedef double vect_t[DIM];

vect_t pos[N];
vect_t old_pos[N];
vect_t vel[N];
vect_t old_vel[N];
vect_t forces[N];

vect_t old_res[N]; // to compare final result between both version

double mass[N];

int main(){
    printf("Parallel version\n");

    init();
    double start_time = omp_get_wtime();

    for (int c = 0; c < CYC; ++c) {
        sA(); // update forces
        mover(DELTA_T); // update velocity and position
    }

    double run_time = omp_get_wtime() - start_time;

    printBodies();
    printf("Run time simple: %f\n", run_time);

    // We use old velocity and position to compare and make sure we have the same results
#pragma omp parallel for // Let's optimize even though it's not in the recorded time
    for (int q = 0; q < N; ++q) {
        for (int d = 0; d < DIM; ++d) {
            old_res[q][d] = pos[q][d]; // back up position to compare both version (correctness)
            pos[q][d] = old_pos[q][d];
            vel[q][d] = old_vel[q][d];
        }
    }

    start_time = omp_get_wtime();

    for (int c = 0; c < CYC; ++c) {
        rA(); // update forces
        mover(DELTA_T); // update velocity and position
    }

    run_time = omp_get_wtime() - start_time;

    // printBodies();
    printf("Run time reduced: %f\n", run_time);

    // This is sadly not really reliable since both methods calculate slight different values that make them diverge from one another
    // as much as there're cycles.
    doubleCheck();

    return 0;
}

void init(){
#pragma omp parallel for
    for(int q=0; q < N; q++){
        for (int d = 0; d < DIM; ++d) {
            forces[q][d] = 0;
        }
    }

    // We'll avoid to mess up with rand which is not thread safe
    for (int q=0; q < N; q++){

        for(int d=0; d < DIM; d++){
            pos[q][d] = (rand() / (double)(RAND_MAX)) * 2 - 1;
            old_pos[q][d] = pos[q][d];
            vel[q][d] = (rand() / (double)(RAND_MAX)) * 2 - 1;
            old_vel[q][d] = vel[q][d];
        }

        mass[q] = fabs((rand() / (double)(RAND_MAX)) * 2 - 1);

    }

}

void rA(){

#pragma omp parallel for
    for(int q=0; q < N; q++){
        for (int d = 0; d < DIM; ++d) {
            forces[q][d] = 0;
        }
    }

#pragma omp parallel for
    for(int q=0; q < N; q++) reducedAlgo(q);
}
void sA(){

#pragma omp parallel for
    for(int q=0; q < N; q++){
        for (int d = 0; d < DIM; ++d) {
            forces[q][d] = 0;
        }
    }
#pragma omp parallel for // This is a safe optimization because the sensible part in this loop is forces[q][d]
// but since each thread will have a specific bunch of distinct q values there won't be any problem such as writing
// same memory space at same time
    for(int q=0; q < N; q++) simpleAlgo(q);
}

void simpleAlgo(int q){

// #pragma omp parallel for // this wouldn't be safe optimization because 2 threads could write forces[q][d] at the same time
    for (int k=0; k < N; k++){

        if(k != q){
            double diff[DIM] = {0};
            double dist = 0;

            for(int d = 0; d < DIM; d++){
                diff[d] = pos[q][d] - pos[k][d]; // reading pos: thread safe
                dist += pow(diff[d], 2);
            }
            dist = sqrt(dist);
            double dist3 = pow(dist, 3);

            for(int d = 0; d < DIM; d++){
                double F = (G*mass[q]*mass[k]/dist3*diff[d]);

                forces[q][d] -= F; // writting forces: potential danger on threading
            }

        }

    }

}

void reducedAlgo(int q){

    // printf("======Th=%d\n", omp_get_num_threads());
    vect_t localForcesK[N] = {0};
    vect_t localForcesQ[N] = {0};

    for (int k=0; k < N; k++){

        if(k > q){
            double diff[DIM] = {0};
            double dist = 0;

            for(int d = 0; d < DIM; d++){
                diff[d] = pos[q][d] - pos[k][d];
                dist += pow(diff[d], 2);
            }
            dist = sqrt(dist);
            double dist3 = pow(dist, 3);


            for(int d = 0; d < DIM; d++){
                double F = (G*mass[q]*mass[k]/dist3*diff[d]);

                localForcesQ[q][d] -= F; // writting forces: potential danger on threading
                // In particular this one, indeed another thread could write this same value
                // if we decide to parallelize the whole function
                localForcesK[k][d] += F; // writting forces: potential danger on threading

            }

        }

    }

    // Other possibility: store localForcesK in a global array (means creating #thread array of the same size of localForcesK
    // And then use threads to sum up everything in a loop over k (in another parallel section).
    // That would be memory consuming but might be faster


    #pragma omp critical // Here each thread is supposed to queue to add its localForcesK values
    {
        // printf("Th=%d\n", omp_get_num_threads());
        for (int k = 0; k < N; k++) {
            if(k > q) { // We don't have time to sum up 0s
                for (int d = 0; d < DIM; d++) {
                    forces[k][d] += localForcesK[k][d];
                }
            }
        }
        // We have to be careful with forces[q][d] as well!
        for (int d = 0; d < DIM; d++) {
            forces[q][d] += localForcesQ[q][d];
        }
        /*
         * At the beginning I let it in the loop over k because I though it would be safe for a thread to write it
         * And then I put the forces[k][d] in the critical region because I was convinced it was the only threat
         * However if a thread is faster than another it can start writing on the forces[k][d] and then interfer with
         * a slower thread that would still be computing on the loop on forces[q][d] so I used a second local array
         *
         */
    }


    // I wonder if the trade would be worth it between memory and time

}

void mover(double delta_t){

#pragma omp parallel for // Same reason than before: vel[q][d] and pos[q][d] are accessed only by 1 thread
    // for a particular value of q
    for(int q=0; q < N; q++){

        for(int d=0; d < DIM; d++){
            vel[q][d] += delta_t/mass[q]*forces[q][d];
            pos[q][d] += delta_t*vel[q][d];
        }

    }
}

void printBodies(){
    char dim[3] = {'X', 'Y', 'Z'};
    printf("\nID,");
    for (int d = 0; d < DIM; ++d) {
        printf("%c,", dim[d]);
    }
    for (int d = 0; d < DIM; ++d) {
        printf("f%c,", dim[d]);
    }
    for (int d = 0; d < DIM; ++d) {
        printf("v%c,", dim[d]);
    }
    printf("mass\n");

    for (int q = 0; q < __min(N,5); ++q) {
        printf("%d,", q);
        for (int d = 0; d < DIM; ++d) {
            printf("%f,", pos[q][d]);
        }
        for (int d = 0; d < DIM; ++d) {
            printf("%f,", forces[q][d]);
        }
        for (int d = 0; d < DIM; ++d) {
            printf("%f,", vel[q][d]);
        }
        printf("%f\n", mass[q]);
    }

}

void doubleCheck(){

    double eps = 0.0001;

#pragma omp parallel for
    for (int q = 0; q < 10; ++q) {
        for (int d = 0; d < DIM; ++d) {
            if (fabs(old_res[q][d] - pos[q][d]) > eps) printf("Error:%f!=%f [%d][%d]\n", old_res[q][d], pos[q][d], q, d);
        }
    }

}

double round(double x){
    return (double) ((int) (x * 10000))/10000;
}