#include <iostream>
//#include <mpi.h>
#include <cmath>
#include <sys/time.h>
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#define A 200000
#define DX 1 // размеры параллелепипеда
#define DY 1
#define DZ 1
#define NX 128 //количество узлов
#define NY 10
#define NZ 10
const double EPSILON = 0.000001;
const double HX = (double)DX/(NX-1);
const double HY = (double)DY/(NY-1);
const double HZ = (double)DZ/(NZ-1);

double count_xi(int x0, int i){
    return (double(x0) + i*HX);
}

double count_yi(int y0, int i){
    return (double(y0) + i*HY);
}

double count_zi(int z0, int i){
    return (double(z0) + i*HZ);
}


double count_fi(double x, double y, double z){
    return (x*x + y*y + z*z);
}

void formula(double** fi, double** prev_fi, int i, int j, int k){
    fi[i][j*NZ + k] = (prev_fi[i+1][j*NZ + k] + prev_fi[i-1][j*NZ + k])/(HX*HX) +
            (prev_fi[i][(j+1)*NZ + k] + prev_fi[i][(j-1)*NZ + k])/(HY*HY) +
            (prev_fi[i][j*NZ + (k+1)] + prev_fi[i][j*NZ + (k-1)])/(HZ*HZ);
    fi[i][j*NZ + k] -= (6 - A*prev_fi[i][j*NZ+k]);
    fi[i][j*NZ + k] /= (2/(HX*HX) + 2/(HY*HY) + 2/(HZ*HZ) + A);
}

double count_res_for_epsilon(double** fi, double** prev_fi, int x_for_proc){
    double max = 0;
    double tmp = 0;
    for(int i = 1; i < x_for_proc-1; i++){
        for(int j = 1; j < NY-1; j++){
            for(int k = 1; k < NZ-1; k++){
                tmp = fabs(fi[i][j*NZ + k] - prev_fi[i][j*NZ + k]);
                if(tmp > max){
                    max = tmp;
                }
            }
        }
    }
    return max;
}

double count_accuracy(double** fi,int rank, int size, int x_for_proc, int x0, int y0, int z0){
    double max = 0;
    double tmp = 0;
    for(int i = 0; i < x_for_proc;i++){
        for(int j = 0; j < NY; j++){
            for(int k = 0; k < NZ; k++){
                if((rank == 0 && i == 0) || (rank == size-1 && i == x_for_proc-1)){
                    continue;
                }
                tmp = fabs(fi[i][j*NZ+k] - count_fi(count_xi(x0, i + rank*(x_for_proc-2)), count_yi(y0, j), count_zi(z0, k)));
                if(tmp > max){
                    max = tmp;
                }
            }
        }
    }
    return max;
}

void fill_block(double** fi, int x_for_proc, int size, int rank,int x0, int y0, int z0){
    for(int i = 1; i < x_for_proc-2; i++){
        for(int j = 0; j < NY; j++){
            for(int k = 0; k < NZ; k++){
                if((rank == 0 && i == 1) || j == 0 || k == 0 ||
                (rank == size-1 && ((i == x_for_proc-2) == NX/size)) ||
                j == NY-1 || k == NZ-1){
                    fi[i][j*NZ + k] = count_fi(count_xi(x0, i - 1 + rank*(x_for_proc-2)), count_yi(y0, j), count_zi(z0, k));
                }else{
                    fi[i][j*NZ + k] = 0;
                }
            }
        }
    }
}

void delete_func(int x_for_proc, double** fi, double** prev_fi, double* sending_first_border, double* sending_last_border, double* prev_proc_last_border, double* next_proc_first_border){
    for(int i = 0; i < x_for_proc;i++){
        delete[] fi[i];
        delete[] prev_fi[i];
    }
    delete[] fi;
    delete[] prev_fi;
    delete[] sending_last_border;
    delete[] sending_first_border;
    delete[] next_proc_first_border;
    delete[] prev_proc_last_border;
}

void algorythm(int argc, char* argv[]){

    MPI_Init(&argc,&argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int dims[2] = { 1, size };
    int periods[2] = { 0, 0 };
    int coords[2] = { 0, 0 };
    int prev_proc, next_proc;
    int reorder = 1;
    MPI_Comm dec_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &dec_comm);
    MPI_Cart_get(dec_comm, 2, dims, periods, coords);
    MPI_Cart_shift(dec_comm, 1, 1, &prev_proc, &next_proc);

    int x0 = -1;
    int y0 = -1;
    int z0 = -1;
    int x_for_proc = (NX/size) + 2;
    double** fi = new double*[x_for_proc];
    for(int i = 0; i < x_for_proc; i++){
        fi[i] = new double[NY*NZ];
    }
    fill_block(fi, x_for_proc, size,rank,x0,y0,z0);
    double** prev_fi = new double*[x_for_proc];
    for(int i = 0; i < x_for_proc; i++){
        prev_fi[i] = new double[NY*NZ];
    }
    if(rank != 0){
        MPI_Sendrecv(fi[1],NY*NZ, MPI_DOUBLE, prev_proc, 1, fi[x_for_proc-1], NY*NZ, MPI_DOUBLE, prev_proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if(rank != size-1){
        MPI_Sendrecv(fi[x_for_proc-2], NY*NZ, MPI_DOUBLE, next_proc, 1, fi[0], NY*NZ, MPI_DOUBLE, next_proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    double max_eps = 0;
    double local_max_eps = 0;
    int iterations = 0;
    double max_accuracy = 0;
    double local_max_accuracy = 0;
    double* sending_last_border = new double[NY*NZ];
    double* sending_first_border = new double[NY*NZ];
    double* prev_proc_last_border = new double[NY*NZ];
    double* next_proc_first_border = new double[NY*NZ];

    do{
        max_eps = 0;
        local_max_eps = 0;

        for(int i =0; i < x_for_proc; i++){
            for(int j=0;j < NY; j++){
                for(int k=0;k < NZ; k++){
                    prev_fi[i][j*NZ + k] = fi[i][j*NZ + k];
                }
            }
        }


        if(rank != 0){
            for(int i = 0; i < NY; i++){
                for(int j = 0; j < NZ; j++){
                    prev_fi[0][i*NZ + j] = fi[0][i*NZ + j];
                }
            }
        }


        if(rank != size-1){
            for(int i = 0;i < NY; i++){
                for (int j = 0; j < NZ; j++) {
                    prev_fi[x_for_proc-1][i*NZ + j] = fi[x_for_proc-1][i*NZ + j];
                }
            }
        }


        if(rank != 0){
            for(int i = 0;i < NY; i++){
                for (int j = 0; j < NZ; j++) {
                    formula(fi, prev_fi, 1, i, j);
                }
            }
        }

        if(rank != size-1){
            for(int i = 0;i < NY; i++){
                for (int j = 0; j < NZ; j++) {
                    formula(fi, prev_fi, x_for_proc-2, i, j);
                }
            }
        }

        for(int i = 1; i < NY-1; i++){
            for(int j = 1; j < NZ-1; j++){
                if(rank != 0){
                    sending_first_border[i*NZ + j] = fi[1][i*NZ + j];
                }
                if(rank != size-1){
                    sending_last_border[i*NZ + j] = fi[x_for_proc-2][i*NZ + j];
                }
            }
        }

        MPI_Request req[4];
        MPI_Status status[4];

        if (rank != 0){
            MPI_Isend(sending_first_border, NY*NZ, MPI_DOUBLE, prev_proc, 1, MPI_COMM_WORLD, &req[0]);
        }

        if(rank != size-1){
            MPI_Isend(sending_last_border, NY*NZ, MPI_DOUBLE, next_proc, 1, MPI_COMM_WORLD, &req[1]);
        }

        for(int i = 2; i < x_for_proc - 2; i++){
            for(int j = 1; j < NY-1; j++){
                for(int k = 1; k < NZ-1; k++){
                    formula(fi, prev_fi, i, j, k);
                }
            }
        }


        if(rank != 0){
            MPI_Irecv(prev_proc_last_border, NY*NZ, MPI_DOUBLE, prev_proc, 1, MPI_COMM_WORLD, &req[2]);
        }

        if(rank != size-1){
            MPI_Irecv(next_proc_first_border, NY*NZ, MPI_DOUBLE, next_proc, 1, MPI_COMM_WORLD, &req[3]);
        }

        if(rank!=0){
            MPI_Wait(&req[0], &status[0]);
            MPI_Wait(&req[2], &status[2]);
        }
        if(rank!=size-1){
            MPI_Wait(&req[1], &status[1]);
            MPI_Wait(&req[3], &status[3]);
        }

        local_max_eps = count_res_for_epsilon(fi, prev_fi, x_for_proc);
        MPI_Allreduce(&local_max_eps, &max_eps, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        iterations++;

        for(int i = 1; i < NY-1; i++){
            for(int j = 1; j < NZ-1; j++){
                if(rank != 0){
                    fi[x_for_proc-1][i*NZ+j] = prev_proc_last_border[i*NZ + j];
                }
                if(rank != size-1){
                    fi[0][i*NZ + j] = next_proc_first_border[i*NZ + j];
                }
            }
        }

    }while(max_eps >= EPSILON);

    local_max_eps = count_res_for_epsilon(fi, prev_fi, x_for_proc);
    local_max_accuracy = count_accuracy(fi, rank, size, x_for_proc, x0, y0, z0);
    MPI_Allreduce(&local_max_eps, &max_eps, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max_accuracy, &max_accuracy, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if(rank == 0){
        std::cout << "ACCURACY: " << max_accuracy << '\n';
        std::cout << "ITERATIONS: " << iterations << '\n';
    }

    delete_func(x_for_proc, fi,prev_fi,sending_first_border, sending_last_border, prev_proc_last_border, next_proc_first_border);
    MPI_Finalize();
}


int main(int argc, char* argv[]) {

    struct timeval start, stop;
    gettimeofday(&start, NULL);

    algorythm(argc, argv);

    gettimeofday(&stop, NULL);
    std::cout << "PROG TIME = " << stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) * 1e-6 << '\n';
    return 0;
}