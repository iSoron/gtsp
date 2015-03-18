cmake_minimum_required(VERSION 3.1)
project(project)

include_directories(/opt/cplex-12.4/distribution/cplex/include/ilcplex)

set(CMAKE_BINARY_DIR bin)
set(CMAKE_SOURCE_DIR src)

set(CMAKE_C_FLAGS "-O3 -g -Wall -pedantic -g --std=c11 -Winline")

set(SOURCE_FILES
    src/lp.c
    src/lp.h
    src/main.c
    src/main.h
    src/util.c
    src/util.h
    src/graph.c
    src/graph.h
    src/tsp.h
    src/tsp.c
    src/branch_and_cut.h
    src/branch_and_cut.c)

add_executable(project ${SOURCE_FILES} )

target_link_libraries(project cplex m pthread)