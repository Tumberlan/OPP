#include <iostream>
#include <cmath>
#include <malloc.h>
//#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

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


#pragma omp parallel for
        for (int i = 0; i < N; i++) {
            double tmp_string = 0;
            for (int j = 0; j < N; j++) {
                tmp_string += A[i * N + j] * B[j];
            }
            res[i] = tmp_string;
        }


}

double scalar(double *a, double *b, int N){
        double res = 0;
#pragma omp parallel reduction (+: res)
    {
        double local_res = 0;
#pragma omp for
        for (int i = 0; i < N; i++) {
            local_res += a[i] * b[i];
        }
    #pragma omp atomic
        res += local_res;
    }
        return res;

}


void const_mul(double *a, double value, int N, double* res){
#pragma omp parallel for
    for(int i = 0; i < N; i++){
        res[i] = a[i] * value;
    }
}

double vector_length(const double *vector, int N){


        double sum = 0;
#pragma omp parallel reduction (+: sum)
    {
        double local_sum = 0;
#pragma omp for
        for (int i = 0; i < N; i++) {
            local_sum += vector[i] * vector[i];
        }
#pragma omp atomic
        sum += local_sum;
    }
        return sqrt(sum);
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

    double* matrix_res = new double [N];
    double* mul_res = new double [N];
    double* z = new double [N];
    double epsilon = 0.00001;
    int iter_counter = 0;

    mul_f(Matrix, x, N, 1, matrix_res);

#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        r[i] = b[i] - matrix_res[i];
    }
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        z[i] = r[i];
    }

    while (double (vector_length(r, N) / vector_length(b, N)) >= epsilon) {

        mul_f(Matrix, z, N, 1, matrix_res);

        double alpha = scalar(r, r, N) / scalar(matrix_res, z, N);

        double* next_x = new double[N];
        const_mul(z, alpha, N, mul_res);

#pragma omp parallel for
        for (int i = 0; i < N; i++) {
            next_x[i] = x[i] + mul_res[i];
        }


        mul_f(Matrix, z, N, 1, matrix_res);

        const_mul(matrix_res, alpha, N, mul_res);


        double* next_r = new double [N];
#pragma omp parallel for
        for (int i = 0; i < N; i++) {
            next_r[i] = r[i] - mul_res[i];
        }

        double beta = scalar(next_r, next_r, N) / scalar(r, r, N);
        double* next_z = new double [N];
        const_mul(z, beta, N, mul_res);

#pragma omp parallel for
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

    delete[] Matrix;

}

int main(int argc, char** argv) {
    srand(100);
    clock_t start, stop;
    start = clock();

    standart_prog(argc,argv);

    stop = clock();
    std::cout << "PROG TIME = " << double((stop - start))/CLOCKS_PER_SEC<< '\n';
    return 0;
}
