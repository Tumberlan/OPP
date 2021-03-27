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


double MPI_vector_length(const double *vector, int N){

    double res[N];

    for (int i = 0; i < N; i++) {
        res[i] = vector[i] * vector[i];
    }
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += res[i];
    }
    return sum;
}

void one_iter_f(double *Matrix, double *b, double *x, double* r, int N){
    double matrix_res[N];
    double mul_res[N];
    double z[N];
    double epsilon = 0.00001;

        std::cout << "vl "<< (vector_length(r, N)/vector_length(b,N)) << "\n";

        while (vector_length(r, N) / vector_length(b, N) >= epsilon) {

            mul_f(Matrix, x, N, 1, matrix_res);

            std::cout<<"matrix_res"<<'\n';
            for(int i = 0; i< N; i++){
                std::cout << matrix_res[i] << '\n';
            }
            for (int i = 0; i < N; i++) {
                r[i] = b[i] - matrix_res[i];
            }

            for (int i = 0; i < N; i++) {
                z[i] = r[i];
            }

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
            std::cout << "top beta" << scalar(next_r,next_r,N) << '\n';
            std::cout << "beta" << beta << '\n';
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
/*
void mpi_1(int argc, char* argv[]){
   MPI_Init(&argc, &argv);
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   int M;
   int N =4;

    M = N/size;
    int collected_M = M;
    for(int i = 1; i < rank+1; i ++){
        M = (N-collected_M)/(size-i);
        collected_M += M;
    }
    collected_M -= M;


    //MPI_Recv(&M, 1, MPI_INT, 0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    double* Matrix = new double [N*M];
    int straighter = 0;
    for(int i = 0; i < N*M; i++){
        if(i+N*collected_M == N*(collected_M + straighter) + (collected_M + straighter)){
            Matrix[i] = 2.0;
            straighter++;
        }else{
            Matrix[i] = 1.0;
        }
    }

    double* x = new double(N);
    double* b = new double(N);
    double* r = new double(N);
    for(int i = 0; i < N; i++){
        b[i] = N+1;
        x[i] = 0;
        r[i] = 1;
    }
    double epsilon = 0.00001;
   // MPI_Recv(&Matrix, N*M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    double* matrix_res = new double[N];
    double* mul_res = new double[N];
    for(int i = 0; i < N; i++){
        matrix_res[i] = 0;
    }



    while(vector_length(r, N) / vector_length(b, N) >= epsilon) {
        MPI_mul_f(Matrix, x, M, 1, matrix_res, rank, collected_M);

        double* final_matrix_res = new double[N];
        MPI_Allreduce(matrix_res, final_matrix_res, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for (int i = 0; i < N; i++) {
            r[i] = b[i] - final_matrix_res[i];
        }

        double* z = new double[N];
        for (int i = 0; i < N; i++) {
            z[i] = r[i];
        }


        MPI_mul_f(Matrix, z, N, M, matrix_res,rank, collected_M);
        MPI_Allreduce(matrix_res, final_matrix_res, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double alpha = scalar(r, r, N) / scalar(final_matrix_res, z, N);

        double* next_x = new double(N);
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

        double beta = scalar(next_r, next_r, N) / scalar(r, r, N);

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
        std::cout << vector_length(r, N) / vector_length(b, N) << "\n";
    }

    for(int i = 0; i < N; i++){
        std::cout << x[i] << '\n';
    }
    MPI_Finalize();
}
*/

void mpi_2(int argc, char* argv[]){
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int M;
    int N = 4;

    M = N/size;
    int collected_M = M;
    for(int i = 1; i < rank+1; i ++){
        M = (N-collected_M)/(size-i);
        collected_M += M;
    }
    collected_M -= M;


    //MPI_Recv(&M, 1, MPI_INT, 0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    double* Matrix = new double [N*M];
    int straighter = 0;
    for(int i = 0; i < N*M; i++){
        if(i+N*collected_M == N*(collected_M + straighter) + (collected_M + straighter)){
            Matrix[i] = 2.0;
            straighter++;
        }else{
            Matrix[i] = 1.0;
        }
    }

    double* x = new double[M];
    double* b = new double[M];
    double* r = new double[M];
    for(int i = 0; i < M; i++){
        b[i] = N+1;
        x[i] = 0;
        r[i] = 1;
    }
    double epsilon = 0.00001;
    // MPI_Recv(&Matrix, N*M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    double* matrix_res = new double[N];
    double* mul_res = new double[N];
    double* final_mul_res = new double[N];
    for(int i = 0; i < N; i++){
        matrix_res[i] = 0;
    }



    double vl_top = MPI_vector_length(r,M);
    double vl_bot = MPI_vector_length(b, M);
    double* vl_top_buf = new double[size];
    double* vl_bot_buf = new double[size];
    double final_vl_top = 0;
    double final_vl_bot = 0;
    MPI_Allreduce(&vl_top, vl_top_buf, M, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&vl_bot, vl_bot_buf, M, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(int i = 0; i < size; i++){
        final_vl_top += vl_top_buf[i];
        final_vl_bot += vl_bot_buf[i];
    }
    int iter_counter = 0;

    while(sqrt(vl_top) / sqrt(vl_bot) >= epsilon) {
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


        MPI2_mul_f(Matrix, z, N, M, matrix_res);
        MPI_Allreduce(matrix_res, final_matrix_res, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for(int i = 0; i < M; i++){
            p_final_matrix_res[i] = final_matrix_res[i + collected_M];
        }

        double alpha_top = scalar(r,r,M);
        double alpha_bot = scalar(p_final_matrix_res,z,M);
        double* alpha_top_buf = new double[1];
        double* alpha_bot_buf = new double[1];
        double final_alpha_top = 0;
        double final_alpha_bot = 0;
        MPI_Allreduce(&alpha_top, alpha_top_buf, M, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&alpha_bot, alpha_bot_buf, M, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for(int i = 0; i < 1; i++){
            final_alpha_top += alpha_top_buf[i];
            final_alpha_bot += alpha_bot_buf[i];
        }
        double alpha = final_alpha_top / final_alpha_bot;

        double* next_x = new double[M];
        const_mul(z, alpha, M, mul_res);
        for (int i = 0; i < N; i++) {
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
        double* beta_top_buf = new double[1];
        double* beta_bot_buf = new double[1];
        double final_beta_top = 0;
        double final_beta_bot = 0;
        MPI_Allreduce(&beta_top, beta_top_buf, M, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&beta_bot, beta_bot_buf, M, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for(int i = 0; i < size; i++){
            final_beta_top += beta_top_buf[i];
            final_beta_bot += beta_bot_buf[i];
        }
        double beta = beta_top / beta_bot;

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
        MPI_Allreduce(&vl_top, vl_top_buf, M, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&vl_bot, vl_bot_buf, M, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for(int i = 0; i < size; i++){
            final_vl_top += vl_top_buf[i];
            final_vl_bot += vl_bot_buf[i];
        }
        iter_counter++;
        delete[] alpha_bot_buf;
        delete[] alpha_top_buf;
        delete[] beta_bot_buf;
        delete[] beta_top_buf;
        delete[] p_final_matrix_res;
        delete[] z;
        delete[] next_r;
        delete[] next_z;
        delete[] next_x;
    }


    double* buf_x = new double [N];
    for(int i = 0; i < M;i++){
        buf_x[collected_M + i] = x[i];
    }
    double * final_buf_x = new double [N];

    MPI_Allreduce(buf_x,final_buf_x, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    std::cout<< "answer x" << '\n';
    for(int i = 0; i < N; i++){
        std::cout << final_buf_x[i] << '\n';
    }

    std::cout<< "iter_counter: " << iter_counter << '\n';
    delete[] mul_res;
    delete[] matrix_res;
    delete[] final_buf_x;
    delete[] buf_x;
    delete[] final_mul_res;
    delete[] x;
    delete[] r;
    delete[] b;
    delete[] Matrix;
    delete[] vl_top_buf;
    delete[] vl_bot_buf;
    MPI_Finalize();
}


int main(int argc, char** argv) {

    standart_prog(argc,argv);
    //mpi_1(argc,argv);
    return 0;
}
