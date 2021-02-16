#include <iostream>
#include <cmath>
//#include <cblas>


void fill_matrix(double* Matrix, int N){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(i == j){
                Matrix[i*N + j] = 2.0;
            }else if(i < j){
                Matrix[i*N + j] = 1.0;
            }else{
                Matrix[i*N + j] = 0.0;
            }
        }
    }
}

void print_matrix(double *Matrix, int N){
    int checker = 0;
    for(int i = 0; i < N*N; i++){
        if(checker == N){
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
        for(int j = 0; j < N; j++){
            tmp_string += A[i*N + j]*B[j];
        }
        res[i] = tmp_string;
    }
}

double scalar(double *a, double *b, int N){
    double res = 0;
    for(int i = 0; i < N; i++){
        res += a[i]*b[i];
    }
    return res;
}

void const_mul(double *a, double value, int N, double* res){
    for(int i = 0; i < N; i++){
        res[i] = a[i] * value;
    }
}

double vector_length(double *vector, int N){
    auto res = new double [N];
    for(int i = 0; i < N; i++){
        res[i] = vector[i]*vector[i];
    }
    double sum = 0;
    for(int i = 0; i < N; i++) {
        sum += res[i];
    }
    return sqrt(sum);
}

void one_iter_f(double *Matrix, double *b, double *x, double* r, int N){
    auto* matrix_res = new double [N];
    auto* mul_res = new double [N];
    //cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, 1, N, 1.0, Matrix,N, b, N, 0.0, res, N);

    mul_f(Matrix, x, N, 1, matrix_res);
    std::cout << "mul res" << '\n';
    for(int i = 0; i < N; i++){
        r[i] = b[i] - matrix_res[i];
    }
    std::cout << "r" << '\n';
    for(int i = 0; i < N; i++){
        std::cout << r[i] << '\n';
    }
    auto z = new double [N];
    for(int i = 0; i < N; i++){
        z[i] = r[i];
    }
    std::cout << "z" << '\n';
    for(int i = 0; i < N; i++){
        std::cout << z[i] << '\n';
    }

    mul_f(Matrix, z, N, 1, matrix_res);
    double alpha = scalar(r, r, N)/scalar(matrix_res, z, N);
    std::cout << "alpha " << alpha << '\n';

    auto next_x = new double [N];
    const_mul(z, alpha, N, mul_res);
    for(int i = 0; i < N; i++){
        next_x[i] = x[i] + mul_res[i];
    }
    std::cout << "next_x" <<'\n';
    for(int i = 0; i<N; i++){
        std::cout << next_x[i] << '\n';
    }


    mul_f(Matrix, z, N, 1, matrix_res);
    const_mul(matrix_res, alpha, N, mul_res);
    auto next_r = new double[N];
    for(int i = 0; i < N; i++){
        next_r[i] = r[i] - mul_res[i];
    }
    std::cout << "next_r" << '\n';
    for(int i = 0; i < N; i++){
        std::cout << next_r[i] << '\n';
    }

    double beta = scalar(next_r, next_r, N)/scalar(r,r,N);
    std::cout << "beta " << beta << '\n';

    auto next_z = new double[N];
    const_mul(z,beta,N,mul_res);
    for(int i = 0; i < N; i++){
        next_z[i] = next_r[i] + mul_res[i];
    }
    std::cout << "next_z" << '\n';
    for(int i = 0; i < N; i++){
        std::cout << next_z[i] << '\n';
    }
}


int main() {

    int N;
    std::cin >> N;

    auto* Matrix = new double[N*N];

    fill_matrix(Matrix, N);
    print_matrix(Matrix, N);

    auto* b = new double[N];
    auto* x = new double[N];
    auto* r = new double [N];
    for(int i = 0; i < N; i++){
        b[i] = N+1;
        x[i] = 1;
        r[i] = 0;
    }

    std::cout << '\n' << "b  x  r" << '\n';
    for(int i = 0; i < N; i++){
        std:: cout << b[i] << "  " << x[i] << "  " << r[i] << '\n';
    }

    one_iter_f(Matrix, b, x, r, N);

    return 0;
}
