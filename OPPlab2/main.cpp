#include <iostream>
#include <cmath>
#include <malloc.h>
//#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#define FROM_NANOS 1000000000

void fill_matrix(double* Matrix, int N){


#pragma omp parallel for
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



void prog_1(int argc, char** argv){
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

void prog_2(int argc, char** argv){



        int N = 4000;

        double *Matrix = new double[N * N];

        fill_matrix(Matrix, N);



        double *b = new double[N];
        double *x = new double[N];
        double *r = new double[N];
        double *u = new double[N];
        double *next_z = new double[N];
        double *next_x = new double[N];
        double *next_r = new double[N];
        double *matrix_res = new double[N];
        double *mul_res = new double[N];
        double *z = new double[N];
        double alpha = 0;
        double beta = 0;
        double sum_top = 0;
        double sum_bot = 0;
        double res_top = 0;
        double res_bot = 0;
        double epsilon = 0.00001;
        int iter_counter = 0;


        for (int i = 0; i < N; i++) {
            b[i] = 0;
            x[i] = rand() % 100;
            r[i] = 1;
            u[i] = rand() % 100;
        }

#pragma omp parallel
    {
//mul_f(Matrix, u, N, 1, b);
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
        for (int i = 0; i < N; i++) {
            double tmp_string = 0;
            for (int j = 0; j < N; j++) {
                tmp_string += Matrix[i * N + j] * u[j];
            }
            b[i] = tmp_string;
        }



        //mul_f(Matrix, x, N, 1, matrix_res);
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
        for (int i = 0; i < N; i++) {
            double tmp_string = 0;
            for (int j = 0; j < N; j++) {
                tmp_string += Matrix[i * N + j] * x[j];
            }
            matrix_res[i] = tmp_string;
        }


//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
        for (int i = 0; i < N; i++) {
            r[i] = b[i] - matrix_res[i];
        }

//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
        for (int i = 0; i < N; i++) {
            z[i] = r[i];
        }


        sum_top = 0;
        sum_bot = 0;
#pragma omp reduction (+: sum_top)
        {
            double local_sum_top = 0;
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
            for (int i = 0; i < N; i++) {
                local_sum_top += r[i] * r[i];
            }
#pragma omp atomic
            sum_top += local_sum_top;
        }
#pragma omp reduction (+: sum_bot)
        {
            double local_sum_bot = 0;
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
            for (int i = 0; i < N; i++) {
                local_sum_bot += b[i] * b[i];
            }
#pragma omp atomic
            sum_bot += local_sum_bot;
        }
#pragma omp barrier

        while (double(sqrt(sum_top) / sqrt(sum_bot)) >= epsilon) {

            //mul_f(Matrix, z, N, 1, matrix_res);
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
            for (int i = 0; i < N; i++) {
                double tmp_string = 0;
                for (int j = 0; j < N; j++) {
                    tmp_string += Matrix[i * N + j] * z[j];
                }
                matrix_res[i] = tmp_string;
            }


            res_top = 0;
            res_bot = 0;
#pragma omp reduction (+: res_top)
            {
                double local_res_top = 0;
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
                for (int i = 0; i < N; i++) {
                    local_res_top += r[i] * r[i];
                }
#pragma omp atomic
                res_top += local_res_top;
            }

#pragma omp reduction (+: res_bot)
            {
                double local_res_bot = 0;
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
                for (int i = 0; i < N; i++) {
                    local_res_bot += matrix_res[i] * z[i];
                }
#pragma omp atomic
                res_bot += local_res_bot;
            }
#pragma omp barrier

            alpha = res_top / res_bot;


            //const_mul(z, alpha, N, mul_res);
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
            for (int i = 0; i < N; i++) {
                mul_res[i] = z[i] * alpha;
            }


//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
            for (int i = 0; i < N; i++) {
                next_x[i] = x[i] + mul_res[i];
            }

            //mul_f(Matrix, z, N, 1, matrix_res);
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
            for (int i = 0; i < N; i++) {
                double tmp_string = 0;
                for (int j = 0; j < N; j++) {
                    tmp_string += Matrix[i * N + j] * z[j];
                }
                matrix_res[i] = tmp_string;
            }


            //const_mul(matrix_res, alpha, N, mul_res);
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
            for (int i = 0; i < N; i++) {
                mul_res[i] = matrix_res[i] * alpha;
            }



//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
            for (int i = 0; i < N; i++) {
                next_r[i] = r[i] - mul_res[i];
            }


            res_top = 0;
            res_bot = 0;
#pragma omp reduction (+: res_top)
            {
                double local_res_top = 0;
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
                for (int i = 0; i < N; i++) {
                    local_res_top += next_r[i] * next_r[i];
                }
#pragma omp atomic
                res_top += local_res_top;
            }

#pragma omp reduction (+: res_bot)
            {
                double local_res_bot = 0;
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
                for (int i = 0; i < N; i++) {
                    local_res_bot += r[i] * r[i];
                }
#pragma omp atomic
                res_bot += local_res_bot;
            }
#pragma omp barrier

            beta = res_top / res_bot;
            //const_mul(z, beta, N, mul_res);
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
            for (int i = 0; i < N; i++) {
                mul_res[i] = z[i] * beta;
            }


//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
            for (int i = 0; i < N; i++) {
                next_z[i] = next_r[i] + mul_res[i];
            }


            for (int i = 0; i < N; i++) {
                z[i] = next_z[i];
                x[i] = next_x[i];
                r[i] = next_r[i];
            }
            iter_counter++;
            sum_top = 0;
            sum_bot = 0;
#pragma omp reduction (+: sum_top)
            {
                double local_sum_top = 0;
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
                for (int i = 0; i < N; i++) {
                    local_sum_top += r[i] * r[i];
                }
#pragma omp atomic
                sum_top += local_sum_top;
            }

#pragma omp reduction (+: sum_bot)
            {
                double local_sum_bot = 0;
//#pragma omp for schedule(static)
//#pragma omp for schedule(static, 2)
//#pragma omp for schedule(dynamic)
//#pragma omp for schedule(dynamic, 2)
//#pragma omp for schedule(guided)
#pragma omp for schedule(guided, 2)     
                for (int i = 0; i < N; i++) {
                    local_sum_bot += b[i] * b[i];
                }
#pragma omp atomic
                sum_bot += local_sum_bot;
            }
#pragma omp barrier

        }
    }

        std::cout << "iter counter: " << iter_counter << '\n';

        delete[] x;
        delete[] z;
        delete[] r;
        delete[] matrix_res;
        delete[] mul_res;

        delete[] Matrix;

}

int main(int argc, char** argv) {
    srand(100);


    struct timespec mt1, mt2;

    long int tt;


    clock_gettime (CLOCK_REALTIME, &mt1);

    prog_1(argc,argv);
    prog_2(argc, argv);

    clock_gettime (CLOCK_REALTIME, &mt2);
    tt=FROM_NANOS*(mt2.tv_sec - mt1.tv_sec)+(mt2.tv_nsec - mt1.tv_nsec);
    std::cout << "PROG TIME = " << (double)tt/FROM_NANOS<< '\n';
    return 0;
}
