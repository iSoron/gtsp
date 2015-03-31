//
// Created by isoron on 18/03/15.
//

#ifndef _PROJECT_GTSP_H_
#define _PROJECT_GTSP_H_

#include "lp.h"
#include "graph.h"
#include "branch_and_cut.h"

struct CLUSTER 
{
	int size;
	int* set;
};

struct GTSP
{
    struct Graph *graph;

    int *clusters;
    int cluster_count;

    double *x_coordinates;
    double *y_coordinates;
    int** dist_matrix;
    struct CLUSTER *vertex_set;
};

int GTSP_create_random_problem(
        int node_count, int cluster_count, int grid_size, struct GTSP *data);

int inital_tour_value(struct GTSP *data); 

void GTSP_free(struct GTSP *data);

int GTSP_init_data(struct GTSP *data);

int GTSP_init_lp(struct LP *lp, struct GTSP *data);

int GTSP_add_cutting_planes(struct LP *lp, struct GTSP *data);

int GTSP_write_problem(struct GTSP *data, char *filename);

int GTSP_write_solution(struct GTSP *data, char *filename, double *x);

int GTSP_main(int argc, char **argv);

int optimize_vertex_in_cluster(struct TOUR * tour, struct GTSP *data);

int two_opt(struct TOUR* tour, struct GTSP *data);

int K_opt(int* tour, struct GTSP *data);

int tour_length(int* tour, struct GTSP* data);

void print_tour(int* tour, struct GTSP* data);

int list_length(struct TOUR* tour, struct GTSP* data);

void print_list(struct TOUR * tour, struct GTSP* data);

extern double *OPTIMAL_X;
extern double FLOW_CPU_TIME;

#endif //_PROJECT_GTSP_H_
