#include <iostream>
//#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#define N 4096
#define M 4096


void init_Matrix(double* Matrix){
    for (int i = 0;i < N; i++){
        for(int j = 0; j<M; j++){
            if(i == j){
                Matrix[i*M + j] = 2;
            }else{
                Matrix[i*M + j] = 1;
            }
        }
    }
}

void print_Matrix(double* Matrix){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            std::cout << Matrix[i*N + j] << " ";
        }
        std::cout << "\n";
    }
}

void print_local_Matrix(double* Matrix, int h, int w){
    for(int i = 0; i < h; i++){
        for(int j = 0; j < w; j++){
            std::cout << Matrix[i*w + j] << " ";
        }
        std::cout << "\n";
    }
}

double row_x_column_res(double *A, double *B, int first_matrix_start, int second_matrix_start){
    double res = 0;
    for (int i = 0; i < N; i++){
        res += A[first_matrix_start+i]*B[second_matrix_start+i];
    }
    return res;
}
void Matrix_mul(double* A, double* B, double* result, int rows, int columns){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            result[i*columns+j] = row_x_column_res(A,B, i*N, j*M);
        }
    }
}


void programm(int argc, char* argv[]){
    MPI_Comm comm_cart, column_comm, row_comm;
    MPI_Init(&argc,&argv);

    int rank, size;


    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    int proc_matrix_x = 4;
    int proc_matrix_y = 4;

    double* A = NULL;
    double* B = NULL;
    double* C = NULL;

    if(rank == 0){
        A = new double [N*M];
        init_Matrix(A);
        B = new double [N*M];
        init_Matrix(B);
        C = new double [N*M];
    }


    MPI_Comm_split(MPI_COMM_WORLD, rank/proc_matrix_x, rank, &row_comm);
    MPI_Comm_split(MPI_COMM_WORLD, rank%proc_matrix_x,rank,&column_comm);
    MPI_Comm_split(MPI_COMM_WORLD, rank<proc_matrix_x*proc_matrix_y,rank, &comm_cart);

    int proc_number_of_rows = N/proc_matrix_x;
    int proc_number_of_columns = M/proc_matrix_y;

    MPI_Datatype column;
    MPI_Type_vector(N, proc_number_of_columns, M, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);

    double* local_A = new double [proc_number_of_rows*N];
    double* local_B =new double [proc_number_of_columns*M];
    double* ex_local_B =new double [proc_number_of_columns*M];

    if(rank % proc_matrix_y == 0){
        MPI_Scatter(A, proc_number_of_rows*N, MPI_DOUBLE, local_A, proc_number_of_rows*N, MPI_DOUBLE, 0, column_comm);
    }
    MPI_Barrier(comm_cart);
    if(rank == 0){
        for(int i = 0; i < proc_number_of_columns; i++) {
            for (int j = 0; j < M; j++) {
                local_B[i * M + j] = B[j * M + i];
            }
        }
        for(int i = 1; i < proc_matrix_y; i++){
            MPI_Send(B+i*proc_number_of_columns,1, column, i, 10, row_comm);
        }
    }
    if(rank > 0 && rank < proc_matrix_y ){
        MPI_Recv(ex_local_B, proc_number_of_columns*N, MPI_DOUBLE, 0, 10, row_comm, MPI_STATUS_IGNORE);
        for(int i = 0; i < proc_number_of_columns; i++){
            for (int j = 0; j < M; j++){
                local_B[i*M + j] = ex_local_B[j*proc_number_of_columns+i];
            }
        }
    }
    MPI_Bcast(local_A, proc_number_of_rows*N, MPI_DOUBLE, 0, row_comm);

    MPI_Bcast(local_B, proc_number_of_columns*M, MPI_DOUBLE, 0, column_comm);

    double* local_C = new double [proc_number_of_rows*proc_number_of_columns];

    Matrix_mul(local_A, local_B, local_C, proc_number_of_rows, proc_number_of_columns);

    MPI_Barrier(comm_cart);
    //print_local_Matrix(local_C, proc_number_of_rows, proc_number_of_columns);
    MPI_Gather(local_C, proc_number_of_rows*proc_number_of_columns, MPI_DOUBLE, C, proc_number_of_rows*proc_number_of_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* result_C = new double [N*M];


    if(rank == 0) {
        int result_idx = 0;
        for (int i = 0; i < size; i++) {
            int x = i % proc_matrix_y;
            int y = i / proc_matrix_y;
            for (int j = y * proc_number_of_rows; j < y * proc_number_of_rows + proc_number_of_rows; j++) {
                for (int k = x * proc_number_of_columns; k < (x + 1) * proc_number_of_columns; k++) {
                    result_C[j * N + k] = C[result_idx];
                    result_idx++;
                }
            }
        }
    }


    MPI_Type_free(&column);
    MPI_Finalize();

    if (rank == 0){
        //print_Matrix(result_C);
        delete[] A;
        delete[] B;
        delete[] C;
    }
    delete[] local_A;
    delete[] local_B;
    delete[] ex_local_B;
    delete[] local_C;
    delete[] result_C;
}




int main(int argc, char* argv[]) {
    srand(100);
    clock_t start, stop;
    start = clock();

    programm(argc, argv);
    stop = clock();
    std::cout << "PROG TIME = " << double((stop - start))/CLOCKS_PER_SEC<< '\n';
    return 0;
}
