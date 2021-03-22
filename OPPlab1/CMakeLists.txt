cmake_minimum_required(VERSION 3.17)
project(OPPlab1)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_COMPILER "C:/Program Files/Microsoft MPI/Bin/mpiexec")
find_package(MPI)
include_directories(SYSTEM ${MSMPI_INC})

add_executable(OPPlab1 main.cpp)

target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
