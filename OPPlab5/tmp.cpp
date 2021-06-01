#include <iostream>
//#include "mpi.h"
#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>
#include <stddef.h>
#include <cmath>
#include <sys/time.h>

#define ITERATIONS 5
#define REQUEST_TAG 10
#define ANSWER_TAG 20
#define NO_TASKS 0
#define SUCCESS 1
#define NEED_TASKS 2
#define TURN_OFF 3
#define L 1000
#define TASK_NUMBER 50 
int* task_list;
double local_res = 0.0;
int size, rank;
int iter_counter = 0;
int start_task_numb = 0;
int given_tasks_amount = 0;
pthread_mutex_t task_mutex;


int get_task_amount(){
    return (int)(5.0 * (double)size / (double)(rank+1));
}

void calculate(int taskIdx) {
    for (int i = 0; i < task_list[taskIdx]; ++i) {
         local_res += exp(sin(i));
    }
}

void reset_task_list() {

    //std::cout << rank << " is resetting task list" << '\n';
    pthread_mutex_lock(&task_mutex);

    for (int i = 0; i < size*TASK_NUMBER; i++){
        task_list[i] = abs(50 - i%TASK_NUMBER)* abs(rank-ITERATIONS%size)*L;
    }
    pthread_mutex_unlock(&task_mutex);
}

int get_task(int from, int* amount_of_get, int* from_idx){
    pthread_mutex_lock(&task_mutex);
    int message_send = NEED_TASKS;
    int message_recieve;
    pthread_mutex_unlock(&task_mutex);

    MPI_Send(&message_send, 1, MPI_INT, from, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&message_recieve, 1, MPI_INT, from, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (message_recieve == NO_TASKS){
        return NO_TASKS;
    }
    pthread_mutex_lock(&task_mutex);
    MPI_Recv(amount_of_get,1,MPI_INT, from, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&task_list[from*TASK_NUMBER], *amount_of_get, MPI_INT, from, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    *from_idx = from*TASK_NUMBER;
    pthread_mutex_unlock(&task_mutex);
    return SUCCESS;
}

void *do_task(void *args){
    long double operation_res = 0;
    int done_task_number = 0;
    while (iter_counter < ITERATIONS){
        reset_task_list();
        struct timeval start, stop;
        gettimeofday(&start, NULL);
        pthread_mutex_lock(&task_mutex);
        given_tasks_amount = 0;
        pthread_mutex_unlock(&task_mutex);
        done_task_number = 0;

        pthread_mutex_lock(&task_mutex);
        start_task_numb = rank*TASK_NUMBER;
        while (start_task_numb< (rank+1)*TASK_NUMBER - given_tasks_amount - 1){
            done_task_number += task_list[start_task_numb];
            pthread_mutex_unlock(&task_mutex);
            for (long int i = 0; i < task_list[start_task_numb]; i++){
                pthread_mutex_lock(&task_mutex);
                calculate(i);
                pthread_mutex_unlock(&task_mutex);
            }
            pthread_mutex_lock(&task_mutex);
            start_task_numb++;

        }
        pthread_mutex_unlock(&task_mutex);

        pthread_mutex_lock(&task_mutex);
        bool check = false;
        pthread_mutex_unlock(&task_mutex);
        int amount_of_get;
        pthread_mutex_lock(&task_mutex);
        int from_idx = 0;
        pthread_mutex_unlock(&task_mutex);
        do {
            check = false;
            for (int i = 0; i < size; i++) {
                if (i!=rank) {
                    if (get_task(i, &amount_of_get, &from_idx) == SUCCESS) {
                        check = true;
                        for(int i = from_idx; i < from_idx + amount_of_get; i++){
                            done_task_number += task_list[i];
                            for (long int j = 0; j < task_list[i]; j++){
                                pthread_mutex_lock(&task_mutex);
                                calculate(j);
                                pthread_mutex_unlock(&task_mutex);
                            }
                        }
                        break;
                    }
                }
            }
        }while(check);
        std::cout << "-K.O.-" << '\n';
        pthread_mutex_lock(&task_mutex);
        iter_counter++;
        pthread_mutex_unlock(&task_mutex);
        gettimeofday(&stop, NULL);
        double task_time;
        task_time = stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) * 1e-6;
        std::cout << "PIPIPUPU" << '\n';
        double max_time = 0;
        double min_time = 0;
        MPI_Allreduce(&task_time,&max_time, 1, MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(&task_time,&min_time, 1, MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        std::cout << "CHECK))))" << '\n';
        double disbalance;
        pthread_mutex_lock(&task_mutex);
        disbalance = max_time-min_time;
        pthread_mutex_unlock(&task_mutex);
        double disb_procentage;
        pthread_mutex_lock(&task_mutex);
        disb_procentage = disbalance/max_time*100;
        pthread_mutex_unlock(&task_mutex);
        std::cout << "----------RESULTS "<< rank << " ----------";
        std::cout << rank << " done " << done_task_number<< " tasks " <<   '\n';
        std::cout << rank << " has result " << operation_res << '\n';
        std::cout << rank << " iterations TIME = " << stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) * 1e-6 << '\n';
        if(rank == 0){
            std::cout << "DISBALANCE " << disbalance << '\n';
            std::cout << "PROCENTAGE OF DISBALANCE " << disb_procentage <<'\n';
        }
    }
    pthread_mutex_lock(&task_mutex);
    int message_send = TURN_OFF;
    pthread_mutex_unlock(&task_mutex);
    MPI_Send(&message_send, 1, MPI_INT, rank, REQUEST_TAG, MPI_COMM_WORLD);
    return NULL;
}

void *send_task(void *args){
    int message_send;
    int message_recieve;
    while (iter_counter < ITERATIONS) {
        MPI_Status status;
        MPI_Recv (&message_recieve, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);

        if (message_recieve == TURN_OFF){
            //std::cout << "turn off "<< rank << '\n';
            break;
        }
        //std::cout << rank << " is sending tasks to " << status.MPI_SOURCE << '\n';
        std::cout << rank << " START NUMBER: " << start_task_numb << " GIVE NUMBER " << given_tasks_amount << '\n';
        std::cout << rank << " LEFT: " << (rank+1)*TASK_NUMBER - given_tasks_amount << '\n';
        std::cout << "BIG LAST CHECK: RANK " << rank << " SOURCE " << status.MPI_SOURCE << '\n';
        if(rank != status.MPI_SOURCE){
            std::cout << "SOURCE: " << status.MPI_SOURCE << '\n';
            pthread_mutex_lock(&task_mutex);
            if (start_task_numb > (rank+1)*TASK_NUMBER - given_tasks_amount - 1) {
                //std::cout << rank << " has no tasks to send" <<'\n';
                std::cout<< "YEAH"<<'\n';
                pthread_mutex_unlock(&task_mutex);
                message_send = NO_TASKS;
                MPI_Send(&message_send, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
                continue;
            }

            pthread_mutex_unlock(&task_mutex);
            int number_want_to_give = get_task_amount();
            int *tasks_to_send = &task_list[(rank+1) * TASK_NUMBER - given_tasks_amount - 1];
            pthread_mutex_lock(&task_mutex);
            given_tasks_amount += number_want_to_give;
            pthread_mutex_unlock(&task_mutex);
            message_send = SUCCESS;
            MPI_Send(&message_send, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            MPI_Send(&number_want_to_give, 1, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
            MPI_Send(&tasks_to_send, number_want_to_give, MPI_INT, status.MPI_SOURCE, ANSWER_TAG, MPI_COMM_WORLD);
        }
    }
    return NULL;
}


int main(int argc, char** argv){
    struct timeval start, stop;
    gettimeofday(&start, NULL);
    int prov;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov);
    if (prov != MPI_THREAD_MULTIPLE) {
        //std::cout << prov << '\n';
        MPI_Finalize();
        return 1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    //std::cout << "Creating new tasks for " << rank << '\n';

    task_list = new int[size * TASK_NUMBER];
    pthread_mutex_init(&task_mutex, NULL);

    pthread_attr_t attrs;
    if (pthread_attr_init(&attrs) != 0){
        perror("Error in setting attributes");
        abort();
    }
    if (pthread_attr_setdetachstate (&attrs, PTHREAD_CREATE_JOINABLE) != 0){
        perror("Error in setting attributes");
        abort();
    }

    pthread_t threads[2];
    pthread_create(&threads[0], &attrs, send_task, NULL);
    pthread_create(&threads[1], &attrs, do_task, NULL);
    pthread_attr_destroy(&attrs);
    for (int i = 0; i < 2; i++){
        if (0 != pthread_join(threads[i], NULL)){
            perror("Can't join a thread");
            abort();
        }
    }

    pthread_mutex_destroy(&task_mutex);

    delete[] task_list;

    MPI_Finalize();
    gettimeofday(&stop, NULL);

    std::cout << "PROG TIME = " << stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) * 1e-6 << '\n';

    return 0;
}

