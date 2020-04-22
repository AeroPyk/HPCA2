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
    printf("Serial version\n");

    init();
    double start_time = omp_get_wtime();

    for (int c = 0; c < CYC; ++c) {
        sA(); // update forces
        mover(DELTA_T); // update velocity and position
    }

    double run_time = omp_get_wtime() - start_time;

    // printBodies();
    printf("Run time simple: %f\n", run_time);

    // We use old velocity and position to compare and make sure we have the same results

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

    printBodies();
    printf("Run time reduced: %f\n", run_time);

    doubleCheck();

    return 0;
}

void init(){

    for(int q=0; q < N; q++){
        for (int d = 0; d < DIM; ++d) {
            forces[q][d] = 0;
        }
    }

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

    for(int q=0; q < N; q++){
        for (int d = 0; d < DIM; ++d) {
            forces[q][d] = 0;
        }
    }

    for(int q=0; q < N; q++) reducedAlgo(q);
}
void sA(){

    for(int q=0; q < N; q++){
        for (int d = 0; d < DIM; ++d) {
            forces[q][d] = 0;
        }
    }

    for(int q=0; q < N; q++) simpleAlgo(q);
}

void simpleAlgo(int q){

    for (int k=0; k < N; k++){

        if(k != q){
            double diff[DIM] = {0};
            double dist = 0;

            for(int d = 0; d < DIM; d++){
                diff[d] = pos[q][d] - pos[k][d];
                dist += pow(diff[d], 2);
            }
            dist = sqrt(dist);
            double dist3 = pow(dist, 3);

            for(int d = 0; d < DIM; d++){
                double F = G*mass[q]*mass[k]/dist3*diff[d];

                forces[q][d] -= F;
            }

        }

    }

}

void reducedAlgo(int q){

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
                double F = G*mass[q]*mass[k]/dist3*diff[d];

                forces[q][d] -= F;
                forces[k][d] += F;
            }

        }

    }

}

void mover(double delta_t){

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