#include <iostream>
#include <cmath>
#include <malloc.h>
//#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"

void fill_matrix(double* Matrix, int N){



       for (int i = 0; i < N; i++) {
           for (int j = 0; j < N; j++) {
               if (i == j) {
                   Matrix[i * N + j] = double(i*i);
               } else if (i < j) {
                   Matrix[i * N + j] = double(i+j);
               } else {
                   Matrix[i * N + j] = double(i+j);
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
            std::cout << Matrix[i] << " ";
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

void MPI_mul_f(double *A, double *B, int N, int M, double *res, int proc_number, int collected_M){
    for(int i = 0; i < M;i++){
        double tmp_string = 0;
        for(int j = 0; j<N; j++){
            tmp_string += A[i*N + j] * B[j];
        }
        res[collected_M+i] = tmp_string;
    }
}

void MPI2_mul_f(double *A, double *B, int N, int M, double *res){
    for(int i = 0; i < N;i++){
        double tmp_string = 0;
        for(int j = 0; j<M; j++){
            tmp_string += A[j*N + i] * B[j];
        }
        res[i] = tmp_string;
    }
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

        double* res = new double[N];

        for (int i = 0; i < N; i++) {
            res[i] = vector[i] * vector[i];
        }
        double sum = 0;
        for (int i = 0; i < N; i++) {
            sum += res[i];
        }
        delete[] res;
        return sqrt(sum);
}


double MPI_vector_length(const double *vector, int N){

    double* res = new double[N];

    for (int i = 0; i < N; i++) {
        res[i] = vector[i] * vector[i];
    }
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += res[i];
    }
    delete[] res;
    return sum;
}

void one_iter_f(double *Matrix, double *b, double *x, double* r, int N){
    double* matrix_res = new double [N];
    double* mul_res = new double [N];
    double* z = new double [N];
    double epsilon = 0.00001;
    int iter_counter = 0;

    mul_f(Matrix, x, N, 1, matrix_res);

    for (int i = 0; i < N; i++) {
        r[i] = b[i] - matrix_res[i];
    }

    for (int i = 0; i < N; i++) {
        z[i] = r[i];
    }

    while (double (vector_length(r, N) / vector_length(b, N)) >= epsilon) {



            mul_f(Matrix, z, N, 1, matrix_res);

            double alpha = scalar(r, r, N) / scalar(matrix_res, z, N);

            double* next_x = new double[N];
            const_mul(z, alpha, N, mul_res);

            for (int i = 0; i < N; i++) {
                next_x[i] = x[i] + mul_res[i];
            }


            mul_f(Matrix, z, N, 1, matrix_res);

            const_mul(matrix_res, alpha, N, mul_res);


            double* next_r = new double [N];
            for (int i = 0; i < N; i++) {
                next_r[i] = r[i] - mul_res[i];
            }

            double beta = scalar(next_r, next_r, N) / scalar(r, r, N);
            double* next_z = new double [N];
            const_mul(z, beta, N, mul_res);


            for (int i = 0; i < N; i++) {
                next_z[i] = next_r[i] + mul_res[i];
            }



            for (int i = 0; i < N; i++) {
                z[i] = next_z[i];
                x[i] = next_x[i];
                r[i] = next_r[i];
            }
            iter_counter++;
            delete[] next_z;
            delete[] next_r;
            delete[] next_x;
        }


        std:: cout << "iter counter: " << iter_counter << '\n';

        delete[] x;
        delete[] z;
        delete[] r;
        delete[] matrix_res;
        delete[] mul_res;

}


void standart_prog(int argc, char** argv){
    int N = 4000;

    double * Matrix= new double [N*N];

    fill_matrix(Matrix, N);

    double* b = new double [N];
    double* x = new double [N];
    double* r = new double [N];
    double* u = new double [N];
    for(int i = 0; i < N; i++){
        b[i] = 0;
        x[i] = rand()%100;
        r[i] = 1;
        u[i] = rand()%100;

    }

    mul_f(Matrix, u, N, 1, b);

    one_iter_f(Matrix, b, x, r, N);

    delete[] Matrix;

}

void mpi_1(int argc, char* argv[]){
   MPI_Init(&argc, &argv);
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   int M;
   int N = 4000;

    M = N/size;
    int collected_M = M;
    for(int i = 1; i < rank+1; i ++){
        M = (N-collected_M)/(size-i);
        collected_M += M;
    }
    collected_M -= M;



    double* Matrix = new double [N*M];
    int straighter = 0;
    int adder;
    if(rank == 0){
        adder = 0;
    }else{
        adder = -1;
    }
    for(int i = 0; i < N*M; i++){
        if(i+N*collected_M == N*(collected_M + straighter) + (collected_M + straighter)){
            Matrix[i] = double((collected_M+straighter)*(collected_M+straighter));
            straighter++;
        }else{
            if(i%N == 0){
                adder++;
            }
            Matrix[i] = double(i%N+collected_M+adder);
        }
    }

    double* x = new double[N];
    double* b = new double[N];
    double* r = new double[N];
    double* u = new double[N];
    for(int i = 0; i < N; i++){     
        b[i] = 0;
        x[i] = rand()%100;
        r[i] = 1;
        u[i] = rand()%100 ;
    }
    double epsilon = 0.00001;


    double* matrix_res = new double[N];
    double* mul_res = new double[N];
    for(int i = 0; i < N; i++){
        matrix_res[i] = 0;
    }

    MPI_mul_f(Matrix, u, N, M, matrix_res, rank, collected_M);
    double* final_matrix_res = new double[N];
    MPI_Allreduce(matrix_res, final_matrix_res, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < N; i++) {
        b[i] = final_matrix_res[i];
    }

    delete[] u;



    int iter_counter = 0;

    MPI_mul_f(Matrix, x, N, M, matrix_res, rank, collected_M);


    MPI_Allreduce(matrix_res, final_matrix_res, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < N; i++) {
        r[i] = b[i] - final_matrix_res[i];
    }

    double* z = new double[N];
    for (int i = 0; i < N; i++) {
        z[i] = r[i];
    }


    while(double(vector_length(r, N) / vector_length(b, N)) >= epsilon) {


        MPI_mul_f(Matrix, z, N, M, matrix_res,rank, collected_M);
        MPI_Allreduce(matrix_res, final_matrix_res, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double alpha = scalar(r, r, N) / scalar(final_matrix_res, z, N);

        double* next_x = new double[N];
        const_mul(z, alpha, N, mul_res);
        for (int i = 0; i < N; i++) {
            next_x[i] = x[i] + mul_res[i];
        }


        MPI_mul_f(Matrix, z, N, M, matrix_res, rank, collected_M);
        MPI_Allreduce(matrix_res, final_matrix_res, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        const_mul(final_matrix_res, alpha, N, mul_res);
        double* next_r = new double[N];
        for (int i = 0; i < N; i++) {
            next_r[i] = r[i] - mul_res[i];
        }

        double beta = scalar(next_r, next_r, N)/scalar(r, r, N);

        double* next_z = new double[N];
        const_mul(z, beta, N, mul_res);
        for (int i = 0; i < N; i++) {
            next_z[i] = next_r[i] + mul_res[i];
        }

        for (int i = 0; i < N; i++) {
            z[i] = next_z[i];
            x[i] = next_x[i];
            r[i] = next_r[i];
        }

        delete[] next_x;
        delete[] next_r;
        delete[] next_z;

        iter_counter++;

    }

    if(rank == 0) {
        std::cout << "iter_counter: " << iter_counter << '\n';
    }
    delete[] Matrix;
    delete[] r;
    delete[] x;
    delete[] matrix_res;
    delete[] mul_res;
    delete[] final_matrix_res;
    delete[] b;
    delete[] z;
    MPI_Finalize();
}


void mpi_2(int argc, char* argv[]){
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int M;
    int N = 4000;

    M = N/size;
    int collected_M = M;
    for(int i = 1; i < rank+1; i ++){
        M = (N-collected_M)/(size-i);
        collected_M += M;
    }
    collected_M -= M;



    double* Matrix = new double [N*M];
    int straighter = 0;
    int adder;
    if(rank == 0){
        adder = 0;
    }else{
        adder = -1;
    }
    for(int i = 0; i < N*M; i++){
        if(i+N*collected_M == N*(collected_M + straighter) + (collected_M + straighter)){
            Matrix[i] = double((collected_M+straighter)*(collected_M+straighter));
            straighter++;
        }else{
            if(i%N == 0){
                adder++;
            }
            Matrix[i] = double(i%N+collected_M+adder);
        }
    }


    double* x = new double[M];
    double* b = new double[M];
    double* r = new double[M];
    double* u = new double[M];
    for(int i = 0; i < M; i++){
        b[i] = 0;
        x[i] = rand()%100;
        r[i] = 1;
        u[i] = rand()%100;
    }
    double epsilon = 0.00001;

    double* matrix_res = new double[N];
    double* mul_res = new double[M];
    double* final_mul_res = new double[N];
    for(int i = 0; i < N; i++){
        matrix_res[i] = 0;
    }

    MPI2_mul_f(Matrix, u, N, M, matrix_res);
    double* matrix_ans = new double [N];
    MPI_Allreduce(matrix_res, matrix_ans, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(int i = 0; i < M; i++){
        b[i] = matrix_ans[i + collected_M];
    }



    MPI2_mul_f(Matrix, x, N, M, matrix_res);

    double* final_matrix_res = new double[N];
    MPI_Allreduce(matrix_res, final_matrix_res, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double* p_final_matrix_res = new double[M];
    for(int i = 0; i < M; i++){
        p_final_matrix_res[i] = final_matrix_res[i + collected_M];
    }

    for (int i = 0; i < M; i++) {
        r[i] = b[i] - p_final_matrix_res[i];
    }

    double* z = new double[M];
    for (int i = 0; i < M; i++) {
        z[i] = r[i];
    }

    delete[] matrix_ans;
    delete[] u;
    double vl_top = MPI_vector_length(r,M);
    double vl_bot = MPI_vector_length(b, M);
    double final_vl_top = 0;
    double final_vl_bot = 0;
    MPI_Allreduce(&vl_top, &final_vl_top, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&vl_bot, &final_vl_bot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    int iter_counter = 0;

    while(double(sqrt(final_vl_top) / sqrt(final_vl_bot)) >= epsilon) {

        MPI2_mul_f(Matrix, z, N, M, matrix_res);
        MPI_Allreduce(matrix_res, final_matrix_res, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for(int i = 0; i < M; i++){
            p_final_matrix_res[i] = final_matrix_res[i + collected_M];
        }

        double alpha_top = scalar(r,r,M);
        double alpha_bot = scalar(p_final_matrix_res,z,M);
        double final_alpha_top = 0;
        double final_alpha_bot = 0;
        MPI_Allreduce(&alpha_top, &final_alpha_top, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&alpha_bot, &final_alpha_bot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double alpha = final_alpha_top / final_alpha_bot;

        double* next_x = new double[M];
        const_mul(z, alpha, M, mul_res);
        for (int i = 0; i < M; i++) {
            next_x[i] = x[i] + mul_res[i];
        }


        MPI2_mul_f(Matrix, z, N, M, matrix_res);
        MPI_Allreduce(matrix_res, final_matrix_res, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for(int i = 0; i < M; i++){
            p_final_matrix_res[i] = final_matrix_res[i + collected_M];
        }

        const_mul(p_final_matrix_res, alpha, M, mul_res);
        double* next_r = new double[M];
        for (int i = 0; i < M; i++) {
            next_r[i] = r[i] - mul_res[i];
        }


        double beta_top = scalar(next_r,next_r,M);
        double beta_bot = scalar(r,r,M);
        double final_beta_top = 0;
        double final_beta_bot = 0;
        MPI_Allreduce(&beta_top, &final_beta_top, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&beta_bot, &final_beta_bot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double beta = final_beta_top / final_beta_bot;

        double* next_z = new double[M];
        const_mul(z, beta, M, mul_res);
        for (int i = 0; i < M; i++) {
            next_z[i] = next_r[i] + mul_res[i];
        }

        for (int i = 0; i < M; i++) {
            z[i] = next_z[i];
            x[i] = next_x[i];
            r[i] = next_r[i];
        }
        vl_top = MPI_vector_length(r,M);
        vl_bot = MPI_vector_length(b, M);
        final_vl_top = 0;
        final_vl_bot = 0;
        MPI_Allreduce(&vl_top, &final_vl_top, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&vl_bot, &final_vl_bot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        iter_counter++;
        delete[] next_r;
        delete[] next_z;
        delete[] next_x;

    }


    double* buf_x = new double [N];
    for(int i = 0; i < N; i++){
        buf_x[i] = 0;
    }
    for(int i = 0; i < M;i++){
        buf_x[collected_M + i] = x[i];
    }
    double * final_buf_x = new double [N];

    MPI_Allreduce(buf_x,final_buf_x, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(rank == 0) {
        std::cout << "iter_counter: " << iter_counter << '\n';
    }
    delete[] mul_res;
    delete[] matrix_res;
    delete[] z;
    delete[] final_buf_x;
    delete[] buf_x;
    delete[] final_mul_res;
    delete[] x;
    delete[] final_matrix_res;
    delete[] p_final_matrix_res;
    delete[] r;
    delete[] b;
    delete[] Matrix;
    MPI_Finalize();
}

int main(int argc, char** argv) {
    srand(100);
    clock_t start, stop;
    start = clock();

    standart_prog(argc,argv);
    //mpi_1(argc,argv);
    //mpi_2(argc, argv);
    stop = clock();
    std::cout << "PROG TIME = " << double((stop - start))/CLOCKS_PER_SEC<< '\n';
    return 0;
}
