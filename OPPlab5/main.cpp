#include <iostream>
#include <cmath>
//#include <mpi.h>
#include <pthread.h>
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#include <sys/time.h>

#define SUCCESS 0
#define REQUEST_TAG 10
#define ANSWER_TAG 20
#define TASK_SENT 0
#define NEED_TASKS 1
#define NO_TASKS 2
#define TURN_OFF 3
#define ITERATIONS 5
#define TASKS_NUMBER 250
#define L 3000

pthread_mutex_t task_mutex;
pthread_mutex_t start_idx_mutex;
pthread_mutex_t given_task_idx_mutex;

int rank, size, iteration = 0;
int* task_list;
int given_task_idx, start_idx;
long double result = 0.0;

int get_task_amount(){
    return size / (rank+1);
}
void calculate(int taskIdx) {
    for (int i = 0; i < task_list[taskIdx]; ++i) {
        result += sqrt(i)/(L*L);
    }
}

void reset_task_list() {

    pthread_mutex_lock(&task_mutex);
    for (int i = 0; i < size*TASKS_NUMBER; i++){
        task_list[i] = std::abs(50 - i%TASKS_NUMBER)* std::abs(rank-ITERATIONS%size)*L;
    }
    pthread_mutex_unlock(&task_mutex);
}

bool get_task(int from, int* amount_of_get, int* from_idx){
    pthread_mutex_lock(&task_mutex);
    int message_send = NEED_TASKS;
    int message_recieve;
    pthread_mutex_unlock(&task_mutex);

    MPI_Send(&message_send, 1, MPI_INT, from, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&message_recieve, 1, MPI_INT, from, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (message_recieve == NO_TASKS){
        return false;
    }
    pthread_mutex_lock(&task_mutex);
    MPI_Recv(amount_of_get,1,MPI_INT, from, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&task_list[from*TASKS_NUMBER], *amount_of_get, MPI_INT, from, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    *from_idx = from*TASKS_NUMBER;
    pthread_mutex_unlock(&task_mutex);
    return true;
}

void* do_task(void*) {
    int task_amount, global_task_amount;
    double max_time, min_time;
    double global_disbalance = 0;
    double global_res;
    double disbalance;
    while (iteration < ITERATIONS) {
        struct timeval start, stop;
        gettimeofday(&start, NULL);
        reset_task_list();
        task_amount = 0;

        pthread_mutex_lock(&start_idx_mutex);
        start_idx = TASKS_NUMBER * rank;
        pthread_mutex_unlock(&start_idx_mutex);
        pthread_mutex_lock(&given_task_idx_mutex);
        given_task_idx = 0;
        pthread_mutex_unlock(&given_task_idx_mutex);

        while (start_idx < TASKS_NUMBER * (rank + 1) - given_task_idx) {
            calculate(start_idx);
            pthread_mutex_lock(&start_idx_mutex);
            start_idx++;
            pthread_mutex_unlock(&start_idx_mutex);
            task_amount++;
        }

        bool areNewTasks;
        do {
            areNewTasks = false;
            for (int currRank = 0; currRank < size; ++currRank) {
                int receivedTasks, otherTaskIdx;
                if (currRank != rank && get_task(currRank, &receivedTasks, &otherTaskIdx)) {
                    for (int i = 0; i < receivedTasks; ++i) {
                        calculate(otherTaskIdx);
                        otherTaskIdx++;
                        task_amount++;
                    }
                    areNewTasks = true;
                    break;
                }
            }
        } while (areNewTasks);

        gettimeofday(&stop, NULL);
        double task_time;
        task_time = stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) * 1e-6;
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Reduce(&task_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&task_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&task_amount, &global_task_amount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        disbalance = max_time-min_time;
        double disb_procentage;
        disb_procentage = disbalance/max_time*100;
        if(rank == 0) {
            std::cout << "----------RESULT " << iteration << " ----------" << '\n';
            std::cout << rank << " done " << global_task_amount << " tasks " <<   '\n';
            std::cout << rank << " has result " << result << '\n';
            std::cout << rank << " ITERATIONS TIME = " << task_time << '\n';

            std::cout << "DISBALANCE " << disbalance << '\n';
            std::cout << "PROCENTAGE OF DISBALANCE " << disb_procentage <<"%"<<'\n';
            global_disbalance += (max_time - min_time) / max_time * 100;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        iteration++;
    }

    int sending_message = TURN_OFF;
    MPI_Send(&sending_message, 1, MPI_INT, rank, REQUEST_TAG, MPI_COMM_WORLD);

    MPI_Reduce(&result, &global_res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << ("----------RESULTS-----------") << '\n';
        std::cout << "GLOBAL RES: " << result << '\n';
        std::cout << "PROCENTAGE OF GLOBAL DISBALANCE: " << global_disbalance / ITERATIONS << "%" << '\n';
    }
    return NULL;
}

void* send_task(void*) {

    MPI_Status recv_status;
    int sending_message, receive_message;

    while (iteration < ITERATIONS) {

        MPI_Recv(&receive_message, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &recv_status);
        int sender = recv_status.MPI_SOURCE;
        bool task_is_needed = receive_message == NEED_TASKS && sender != rank;
        pthread_mutex_lock(&start_idx_mutex);
        bool nothing_to_send = start_idx > (rank + 1) * TASKS_NUMBER - given_task_idx;
        pthread_mutex_unlock(&start_idx_mutex);
        if (task_is_needed && nothing_to_send) {
            sending_message = NO_TASKS;
            MPI_Send(&sending_message, 1, MPI_INT, sender, ANSWER_TAG, MPI_COMM_WORLD);
            continue;
        }
        if (task_is_needed) {
            int given_tasks_amount = get_task_amount();
            pthread_mutex_lock(&given_task_idx_mutex);
            given_task_idx += given_tasks_amount;
            pthread_mutex_unlock(&given_task_idx_mutex);
            pthread_mutex_lock(&task_mutex);
            int *send_tasks = &task_list[(rank + 1) * TASKS_NUMBER - given_task_idx];
            pthread_mutex_unlock(&task_mutex);
            sending_message = TASK_SENT;
            MPI_Send(&sending_message, 1, MPI_INT, sender, ANSWER_TAG, MPI_COMM_WORLD);
            MPI_Send(&given_tasks_amount, 1, MPI_INT, sender, ANSWER_TAG, MPI_COMM_WORLD);
            MPI_Send(send_tasks, given_tasks_amount, MPI_INT, sender, ANSWER_TAG, MPI_COMM_WORLD);
            continue;
        }
        if (receive_message == TURN_OFF) {
            break;
        }
    }
    return NULL;
}

int main(int argc, char *argv[]) {

    struct timeval strt, stp;
    gettimeofday(&strt, NULL);

    int providedLevel;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &providedLevel);

    if (providedLevel != MPI_THREAD_MULTIPLE) {
        fprintf(stderr, "Error on pthread init");
        MPI_Finalize();
        return 1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    task_list = new int [TASKS_NUMBER * size]();

    pthread_attr_t attrs;
    int attrsInitRes = pthread_attr_init(&attrs);
    if (attrsInitRes != SUCCESS) {
        perror("Error on initializing attrs");
        MPI_Finalize();
        return 1;
    }
    int mutexInitRes = pthread_mutex_init(&task_mutex, NULL);
    if (mutexInitRes != SUCCESS) {
        perror("Error on initializing mutex");
        MPI_Finalize();
        return 1;
    }
    mutexInitRes = pthread_mutex_init(&given_task_idx_mutex, NULL);
    if (mutexInitRes != SUCCESS) {
        perror("Error on initializing mutex");
        MPI_Finalize();
        return 1;
    }
    mutexInitRes = pthread_mutex_init(&start_idx_mutex, NULL);
    if (mutexInitRes != SUCCESS) {
        perror("Error on initializing mutex");
        MPI_Finalize();
        return 1;
    }
    int setStateRes = pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);
    if (setStateRes != SUCCESS) {
        perror("Error on setting detach state");
        MPI_Finalize();
        return 1;
    }

    pthread_t threads[2];
    pthread_create(&threads[0], &attrs, send_task, NULL);
    pthread_create(&threads[1], &attrs, do_task, NULL);
    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);

    pthread_mutex_destroy(&task_mutex);
    pthread_mutex_destroy(&start_idx_mutex);
    pthread_mutex_destroy(&given_task_idx_mutex);
    pthread_attr_destroy(&attrs);
    delete[] task_list;

    MPI_Finalize();
    gettimeofday(&stp, NULL);
    std::cout << "PROG_TIME: " << stp.tv_sec - strt.tv_sec + (stp.tv_usec - strt.tv_usec) * 1e-6 << '\n';
    return 0;
}
