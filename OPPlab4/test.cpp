#include "mpi.h"
//#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>

const double EPS = 0.00001;
const int a = 100000;

const int DX = 2;
const int DY = 2;
const int DZ = 2;

const int NX = 10;
const int NY = 10;
const int NZ = 128;

const double HX = (double)DX / (NX - 1);
const double HY = (double)DY / (NY - 1);
const double HZ = (double)DZ / (NZ - 1);

//calculating coords;
double xi(int x0, int i) {
    return ((double)x0 + i * HX);
}

double yj(int y0, int j) {
    return ((double)y0 + j * HY);
}

double zk(int z0, int k) {
    return ((double)z0 + k * HZ);
}
//end of calculating coords

//formula of the iteration proccess;
void iteration(double** phi, double** prevPhi, int i, int j, int k) {
    phi[k][j * NX + i] = (prevPhi[k + 1][j * NX + i] + prevPhi[k - 1][j * NX + i]) / (HZ * HZ) +
                         (prevPhi[k][(j + 1) * NX + i] + prevPhi[k][(j - 1) * NX + i]) / (HY * HY) +
                         (prevPhi[k][j * NX + i + 1] + prevPhi[k][j * NX + i - 1]) / (HX * HX);
    phi[k][j * NX + i] -= (6 - a * prevPhi[k][j * NX + i]);		//sum(phi) - rho;
    phi[k][j * NX + i] /= (2 / (HX * HX) + 2 / (HY * HY) + 2 / (HZ * HZ) + a);

}


//return max_eps
double calcEps(double** phi, double** prevPhi, int linesForProc) {
    double max_eps = 0.0;
    double cur_eps;
    for (int k = 1; k < linesForProc - 1; k++) {
        for (int j = 1; j < NY - 1; j++) {
            for (int i = 1; i < NX - 1; i++) {
                cur_eps = fabs(phi[k][j * NX + i] - prevPhi[k][j * NX + i]);
                if (cur_eps > max_eps) max_eps = cur_eps;
            }
        }
    }

    return max_eps;
}

double calcPhi(double x, double y, double z) {
    return (x * x + y * y + z * z);
}

double calcAcc(int x0, int y0, int z0, int zForProc, double** phi, int rank, int size) {
    double max_acc = 0.0;
    double cur_acc = 0.0;

    for (int k = 0; k < zForProc; k++) {
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                if (rank == 0 && k == 0) continue;
                if (rank == size-1 && k == zForProc-1) continue;
                cur_acc = fabs(phi[k][j * NX + i] - calcPhi(xi(x0, i), yj(y0, j), zk(z0, k + rank * (zForProc - 2))));

                if (cur_acc > max_acc) max_acc = cur_acc;
            }
        }
    }
    return max_acc;
}

int main(int argc, char** argv) {


    MPI_Init(&argc, &argv);
    int size, rank, rank_dec;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf( "AAAAAAAAAA\n");
    int dims[2] = { 1, size };
    int periods[2] = { 0, 0 };
    int coords[2] = { 0, 0 };
    int prevx, nextx;
    int reorder = 1;
    MPI_Comm comm_dec;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm_dec);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_dec);
    //получили координаты процесса
    MPI_Cart_get(comm_dec, 2, dims, periods, coords);
    MPI_Cart_shift(comm_dec, 1, 1, &prevx, &nextx);
    struct timeval tv1, tv2;
    int x0 = -1, y0 = -1, z0 = -1;
    int zForProc = (NZ / size) + 2; //oof
    double** phi = (double**)malloc(zForProc * sizeof(double*));
    for (int i = 0; i < zForProc; i++) {
        phi[i] = (double*)malloc(NX * NY * sizeof(double));
    }
    //printf("\n%d\n", size);


    for (int k = 1; k < zForProc - 2; k++) {
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) { //k=1 в первом условии потому что с этого момента идёт собственная строка блока
                if ((rank_dec == 0 && k == 1) || j == 0 || i == 0 || ((rank_dec == size - 1) && ((k == zForProc - 2) == (NZ / size)))
                    || j == NY - 1 || i == NX - 1) {

                    phi[k][j * NX + i] = calcPhi(xi(x0, i), yj(y0, j), zk(z0, k - 1 + rank_dec * (zForProc - 2)));
                    printf("%lf k=%d, j=%d, i=%d\n", phi[k][j*NX + i], k, j, i);
                }
                else phi[k][j * NX + i] = 0;
            }
        }
    }


    double** prevPhi = (double**)malloc(zForProc * sizeof(double*));
    for (int i = 0; i < zForProc; i++) {
        prevPhi[i] = (double*)malloc(NX * NY * sizeof(double));
    }

    gettimeofday(&tv1, NULL);

    //первый собственный "слой" процесса передаем в последний(технический) "слой" предыдущего
    if (rank_dec != 0) {
        MPI_Sendrecv(phi[1], NX * NY, MPI_DOUBLE, prevx, 12345,
                     phi[zForProc - 1], NX * NY, MPI_DOUBLE, prevx, 12345, comm_dec, MPI_STATUS_IGNORE);
    }

    //последний собственный "слой" процесса передаем в первый(технический) "слой" последующего
    if (rank_dec != size - 1) {
        MPI_Sendrecv(phi[zForProc - 2], NX * NY, MPI_DOUBLE, nextx, 12345,
                     phi[0], NX * NY, MPI_DOUBLE, nextx, 12345, comm_dec, MPI_STATUS_IGNORE);
    }
    //printf("Its ok\n");
    double allMax = 0.0;
    int iter = 0;
    double locMax = 0.0;
    double loc_acc = 0.0;
    double glob_acc = 0.0;

    double* lastLine = (double*)malloc(NX * NY * sizeof(double));
    double* firstLine = (double*)malloc(NX * NY * sizeof(double));

    double* getFirst = (double*)malloc(NX * NY * sizeof(double));
    double* getLast = (double*)malloc(NX * NY * sizeof(double));
    do {
        allMax = 0.0;
        locMax = 0.0;

        //скопировали предыдущие значения вместе с границами
        for (int k = 0; k < zForProc; k++) {
            for (int j = 0; j < NY; j++) {
                for (int i = 0; i < NX; i++) {
                    prevPhi[k][j * NX + i] = phi[k][j * NX + i];
                }
            }
        }

        if (rank_dec != 0)
            for (int j = 0; j < NY; j++) {
                for (int i = 0; i < NX; i++) {
                    prevPhi[0][j * NX + i] = phi[0][j * NX + i];
                }
            }


        if (rank_dec != size - 1)
            for (int j = 0; j < NY; j++) {
                for (int i = 0; i < NX; i++) {
                    prevPhi[zForProc - 1][j * NX + i] = phi[zForProc - 1][j * NX + i];
                }
            }


        //вычисляем нижнюю границу
        if (rank_dec != size - 1)
            for (int k = zForProc - 2; k < zForProc - 1; k++) {
                for (int j = 1; j < NY - 1; j++) {
                    for (int i = 1; i < NX - 1; i++) {
                        iteration(phi, prevPhi, i, j, k);
                    }
                }
            }
        //теперь верхнюю
        if (rank_dec != 0)
            for (int k = 1; k < 2; k++) {
                for (int j = 1; j < NY - 1; j++) {
                    for (int i = 1; i < NX - 1; i++) {
                        iteration(phi, prevPhi, i, j, k);
                    }
                }
            }

        for (int j = 1; j < NY - 1; j++)
            for (int i = 1; i < NX - 1; i++) {
                if (rank_dec != 0) lastLine[j * NX + i] = phi[1][j * NX + i];
                if (rank_dec != size - 1) firstLine[j * NX + i] = phi[zForProc - 2][j * NX + i];
            }
        MPI_Request req[4];
        MPI_Status status[4];

        //посылаем нижнюю
        if (rank_dec != 0) {
            MPI_Isend(lastLine, NX * NY, MPI_DOUBLE, prevx, 12345, comm_dec, &req[0]);
        }
        //посылаем верхнюю:
        if (rank_dec != size - 1) {
            MPI_Isend(firstLine, NX * NY, MPI_DOUBLE, nextx, 12345, comm_dec, &req[1]);
        }

        //вычисляем остальные строки блока
        for (int k = 2; k < zForProc - 2; k++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int i = 1; i < NX - 1; i++) {
                    iteration(phi, prevPhi, i, j, k);
                }
            }
        }

        //получаем границы:
        //нижняя:
        if (rank_dec != 0) {
            MPI_Irecv(getLast, NX * NY, MPI_DOUBLE, prevx, 12345, comm_dec, &req[2]);
        }
        //верхняя:
        if (rank_dec != size - 1) {
            MPI_Irecv(getFirst, NX * NY, MPI_DOUBLE, nextx, 12345, comm_dec, &req[3]);
        }


        if (rank_dec != 0) {
            MPI_Wait(&req[0], &status[0]);
            MPI_Wait(&req[2], &status[2]);
        }
        if (rank_dec != size - 1) {
            MPI_Wait(&req[1], &status[1]);
            MPI_Wait(&req[3], &status[3]);
        }

        locMax = calcEps(phi, prevPhi, zForProc);
        printf("EPS: %lf\n", locMax);
        MPI_Allreduce(&locMax, &allMax, 1, MPI_DOUBLE, MPI_MAX, comm_dec);
        printf("MAX_EPS: %lf\n", allMax);
        iter++;

        for (int j = 1; j < NY - 1; j++)
            for (int i = 1; i < NX - 1; i++) {
                if (rank_dec != 0) phi[zForProc - 1][j * NX + i] = getLast[j * NX + i];
                if (rank_dec != size - 1) phi[0][j * NX + i] = getFirst[j * NX + i];
            }
        //printf("Cur eps: %lf\n", allMax);
    } while (allMax >= EPS);

    gettimeofday(&tv2, NULL);

    double time = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec) * 1e-6;

    locMax = calcEps(phi, prevPhi, zForProc);
    loc_acc = calcAcc(x0, y0, z0, zForProc, phi, rank_dec, size);
    MPI_Allreduce(&loc_acc, &glob_acc, 1, MPI_DOUBLE, MPI_MAX, comm_dec);
    MPI_Allreduce(&locMax, &allMax, 1, MPI_DOUBLE, MPI_MAX, comm_dec);
    if (rank_dec == 0) {
        printf("Accuracity: %lf\n", glob_acc);
        printf("Time: %lf\n", time);
        printf("Iterations: %d\nMax eps: %1.20f\n", iter, allMax);
    }

    MPI_Finalize();
    return 0;
}
