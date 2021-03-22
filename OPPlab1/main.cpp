#include <iostream>
#include <cmath>
#include <malloc.h>
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
//#include <cblas>
#define CLEAR 0
#define MPI1 1
#define MPI2 2

void fill_matrix(double* Matrix, int N){



       for (int i = 0; i < N; i++) {
           for (int j = 0; j < N; j++) {
               if (i == j) {
                   Matrix[i * N + j] = 2.0;
               } else if (i < j) {
                   Matrix[i * N + j] = 1.0;
               } else {
                   Matrix[i * N + j] = 1.0;
               }
           }
       }

}

void print_matrix(double *Matrix, int N){
        int checker = 0;
        for (int i = 0; i < N * N; i++) {
            if (checker == N) {
                std::cout << '\n';
                checker = 0;
            }
            std::cout << Matrix[i];
            checker++;
        }
        std::cout << '\n';

}



void mul_f(double *A, double *B, int N, int M, double *res){


        for (int i = 0; i < N; i++) {
            double tmp_string = 0;
            for (int j = 0; j < N; j++) {
                tmp_string += A[i * N + j] * B[j];
            }
            res[i] = tmp_string;
        }

}

void MPI_mul_f(double *A, double *B, int N, int M, double *res, int proc_number){
    for(int i = 0; i < M;i++){
        double tmp_string = 0;
        for(int j = 0; j<N; j++){
            tmp_string += A[i*N + j] * B[j];
        }
        res[M*proc_number+i] = tmp_string;
    }
    std::cout << "MPI res" << '\n';
    for (int i = 0; i < N; i++){
        std::cout << res[i] << '\n';
    }
    std::cout << '\n';
}

double scalar(double *a, double *b, int N){
        double res = 0;
        for (int i = 0; i < N; i++) {
            res += a[i] * b[i];
        }
        return res;

}

void const_mul(double *a, double value, int N, double* res){
    for(int i = 0; i < N; i++){
        res[i] = a[i] * value;
    }
}

double vector_length(const double *vector, int N){

        double res[N];

        for (int i = 0; i < N; i++) {
            res[i] = vector[i] * vector[i];
        }
        double sum = 0;
        for (int i = 0; i < N; i++) {
            sum += res[i];
        }
        return sqrt(sum);

}

void one_iter_f(double *Matrix, double *b, double *x, double* r, int N){
    double matrix_res[N];
    double mul_res[N];
    //cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, 1, N, 1.0, Matrix,N, b, N, 0.0, res, N);
    double z[N];
    double epsilon = 0.00001;

        std::cout << vector_length(r, N) << " " << vector_length(b, N) << "\n";

        while (vector_length(r, N) / vector_length(b, N) >= epsilon) {

            mul_f(Matrix, x, N, 1, matrix_res);
            std::cout << '\n' << "matrix_res" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << matrix_res[i] << '\n';
            }
            for (int i = 0; i < N; i++) {
                r[i] = b[i] - matrix_res[i];
            }
            std::cout << '\n' << "r" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << r[i] << '\n';
            }
            std::cout<<'\n';
            for (int i = 0; i < N; i++) {
                z[i] = r[i];
            }

            std::cout << '\n' << "z" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << z[i] << '\n';
            }
            std::cout<<'\n';

            mul_f(Matrix, z, N, 1, matrix_res);
            std::cout << '\n' << "matrix_res" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << matrix_res[i] << '\n';
            }
            std::cout<<'\n';
            double alpha = scalar(r, r, N) / scalar(matrix_res, z, N);
            std::cout << '\n' << "alpha" << '\n';
            std::cout << alpha << '\n';
            std::cout<<'\n';

            double* next_x = new double[N];
            const_mul(z, alpha, N, mul_res);
            std::cout << '\n' << "mul_res" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << mul_res[i] << '\n';
            }
            std::cout<<'\n';
            for (int i = 0; i < N; i++) {
                next_x[i] = x[i] + mul_res[i];
            }


            std::cout << '\n' << "next_x" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << next_x[i] << '\n';
            }
            std::cout<<'\n';


            mul_f(Matrix, z, N, 1, matrix_res);
            std::cout << '\n' << "matrix_res" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << matrix_res[i] << '\n';
            }
            std::cout<<'\n';
            const_mul(matrix_res, alpha, N, mul_res);

            std::cout << '\n' << "mul_res" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << mul_res[i] << '\n';
            }
            std::cout<<'\n';
            double* next_r = new double [N];
            for (int i = 0; i < N; i++) {
                next_r[i] = r[i] - mul_res[i];
            }
            std::cout << '\n' << "next_r" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << next_r[i] << '\n';
            }
            std::cout<<'\n';

            double beta = scalar(next_r, next_r, N) / scalar(r, r, N);
            std::cout << '\n' << "beta" << '\n';
            std::cout << beta << '\n';
            std::cout<<'\n';
            double* next_z = new double [N];
            const_mul(z, beta, N, mul_res);
            std::cout << '\n' << "mul_res" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << mul_res[i] << '\n';
            }
            std::cout<<'\n';
            for (int i = 0; i < N; i++) {
                next_z[i] = next_r[i] + mul_res[i];
            }
            std::cout << '\n' << "next_z" << '\n';
            for(int i = 0; i < N; i++){
                std::cout << next_z[i] << '\n';
            }
            std::cout<<'\n';

            for (int i = 0; i < N; i++) {
                z[i] = next_z[i];
                x[i] = next_x[i];
                r[i] = next_r[i];
            }
            std::cout <<"epsilon checker: " << vector_length(r, N) / vector_length(b, N) << "\n";
        }

}


void standart_prog(int argc, char** argv){
    int N;
    std::cin >> N;

    double * Matrix= new double [N*N];

    fill_matrix(Matrix, N);
    //print_matrix(Matrix, N);

    double* b = (double*)malloc(sizeof (double)*N);
    double* x = new double [N];
    double* r = new double [N];
    for(int i = 0; i < N; i++){
        b[i] = N+1;
        x[i] = 0;
        r[i] = 1;
    }

    std::cout << '\n' << "b  x  r" << '\n';
    for(int i = 0; i < N; i++){
        std:: cout << b[i] << "  " << x[i] << "  " << r[i] << '\n';
    }

    one_iter_f(Matrix, b, x, r, N);

    std::cout << "result x" << '\n';
    for(int i = 0; i < N; i++){
        std:: cout << x[i] << '\n';
    }

    std::cout << '\n';
}
//Mpi_comm_size
//Mpi_comm_rank
//Mpi_send
//Mpi_recieve
/*
void mpi_1(int argc, char* argv[]){
   MPI_Init(&argc, &argv);
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if(rank == 0){
       int N;
       std::cin >> N;

       double Matrix[N*N];
       double* b = new double(N);
       double* x = new double(N);
       double* r = new double(N);

       fill_matrix(Matrix, N);

       for(int i = 0; i < N; i++){
           b[i] = i+1;
           x[i] = 0;
           r[i] = 1;
       }

       int ex_N = N;
       int ex_size = size;
       int skip = 0;
       for(int i = 0; i < size; i++){
           int M = ex_N/ex_size;
           MPI_Send(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
           MPI_Send(&M, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
           MPI_Send(Matrix+N*skip, N*M,MPI_DOUBLE, i,0,MPI_COMM_WORLD);
           MPI_Send(x, N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
           MPI_Send(b, N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
           MPI_Send(r,N,MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
           skip += M;
           ex_N -= M;
           ex_size--;
       }
   }
    int N, M;
    MPI_Recv(&N, 1, MPI_INT, 0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&M, 1, MPI_INT, 0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    double Matrix[N*M];
    double* x = new double(N);
    double* b = new double(N);
    double* r = new double(N);
    double epsilon = 0.00001;
    MPI_Recv(&Matrix, N*M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&b, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r,N,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    double* matrix_res = new double(N);
    double* mul_res = new double(N);
    for(int i = 0; i < N; i++){
        matrix_res[i] = 0;
    }



    while(vector_length(r, N) / vector_length(b, N) >= epsilon) {
        MPI_mul_f(Matrix, x, M, 1, matrix_res, rank);

        double final_matrix_res[N];
        MPI_Allreduce(&matrix_res, &final_matrix_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for (int i = 0; i < N; i++) {
            r[i] = b[i] - final_matrix_res[i];
        }


        std::cout << '\n' << "r" << '\n';
        for(int i = 0;i < N; i++){
            std::cout << r[i] << '\n';
        }
        std::cout << '\n';
        //
        double z[N];
        for (int i = 0; i < N; i++) {
            z[i] = r[i];
        }


        mul_f(Matrix, z, N, 1, final_matrix_res);
        double alpha = scalar(r, r, N) / scalar(final_matrix_res, z, N);

        double* next_x = new double(N);
        const_mul(z, alpha, N, mul_res);
        for (int i = 0; i < N; i++) {
            next_x[i] = x[i] + mul_res[i];
        }


        mul_f(Matrix, z, N, 1, matrix_res);
        const_mul(final_matrix_res, alpha, N, mul_res);
        double* next_r = new double(N);
        for (int i = 0; i < N; i++) {
            next_r[i] = r[i] - mul_res[i];
        }

        double beta = scalar(next_r, next_r, N) / scalar(r, r, N);

        double* next_z = new double(N);
        const_mul(z, beta, N, mul_res);
        for (int i = 0; i < N; i++) {
            next_z[i] = next_r[i] + mul_res[i];
        }

        for (int i = 0; i < N; i++) {
            z[i] = next_z[i];
            x[i] = next_x[i];
            r[i] = next_r[i];
        }
        std::cout << vector_length(r, N) / vector_length(b, N) << "\n";
    }

    for(int i = 0; i < N; i++){
        std::cout << x[i] << '\n';
    }
    MPI_Finalize();
}

void mpi_2(int argc, char* argv[]){
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){
        int N;
        std::cin >> N;

        double Matrix[N*N];
        double b[N];
        double x[N];
        double r[N];

        fill_matrix(Matrix, N);

        for(int i = 0; i < N; i++){
            b[i] = i+1;
            x[i] = 0;
            r[i] = 1;
        }

        int ex_N = N;
        int ex_size = size;
        int skip = 0;
        for(int i = 0; i < size; i++){
            int M = ex_N/ex_size;
            MPI_Send(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&M, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(Matrix+N*skip, N*M,MPI_DOUBLE, i,0,MPI_COMM_WORLD);
            MPI_Send(x+M*skip, M, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send(b+skip, M, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            skip += M;
            ex_N -= M;
            ex_size--;
        }
    }
    int N, M;
    MPI_Recv(&N, 1, MPI_INT, 0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&M, 1, MPI_INT, 0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    double Matrix[N*M];
    double x[N];
    double b[M];
    double r[N];
    double epsilon = 0.00001;
    MPI_Recv(&Matrix, N*M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&b, M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    double matrix_res[N];
    double mul_res[N];
    for(int i = 0; i < N; i++){
        matrix_res[i] = 0;
    }



    while(vector_length(r, N) / vector_length(b, N) >= epsilon) {
            MPI_mul_f(Matrix, x, M, 1, matrix_res, rank);

        MPI_Allreduce(&x, &matrix_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for (int i = 0; i < N; i++) {
            r[i] = b[i] - matrix_res[i];
        }

        //
        double z[N];
        for (int i = 0; i < N; i++) {
            z[i] = r[i];
        }

        mul_f(Matrix, z, N, 1, matrix_res);
        double alpha = scalar(r, r, N) / scalar(matrix_res, z, N);

        double next_x[N];
        const_mul(z, alpha, N, mul_res);
        for (int i = 0; i < N; i++) {
            next_x[i] = x[i] + mul_res[i];
        }


        mul_f(Matrix, z, N, 1, matrix_res);
        const_mul(matrix_res, alpha, N, mul_res);
        double next_r[N];
        for (int i = 0; i < N; i++) {
            next_r[i] = r[i] - mul_res[i];
        }

        double beta = scalar(next_r, next_r, N) / scalar(r, r, N);

        double next_z[N];
        const_mul(z, beta, N, mul_res);
        for (int i = 0; i < N; i++) {
            next_z[i] = next_r[i] + mul_res[i];
        }

        for (int i = 0; i < N; i++) {
            z[i] = next_z[i];
            x[i] = next_x[i];
            r[i] = next_r[i];
        }
        std::cout << vector_length(r, N) / vector_length(b, N) << "\n";
    }

    for(int i = 0; i < N; i++){
        std::cout << x[i] << '\n';
    }
    MPI_Finalize();
}
*/

int main(int argc, char** argv) {

    standart_prog(argc,argv);
    //mpi_1(argc,argv);
    return 0;
}
