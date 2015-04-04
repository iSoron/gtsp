#ifndef _PROJECT_GTSP_H_
#define _PROJECT_GTSP_H_

#include "lp.h"
#include "graph.h"
#include "branch_and_cut.h"

struct Tour
{
    int vertex;
    int next;
    int prev;
};

struct Cluster
{
	int size;
	int* nodes;
};

struct GTSP
{
    struct Graph *graph;

    int** dist_matrix;

    int cluster_count;
    int *node_to_cluster;
    struct Cluster *clusters;
};

int GTSP_create_random_problem(
        int node_count, int cluster_count, int grid_size, struct GTSP *data);

int inital_tour_value(struct GTSP *data, int *value, double *x);

void GTSP_free(struct GTSP *data);

int GTSP_init_data(struct GTSP *data);

int GTSP_init_lp(struct LP *lp, struct GTSP *data);

int GTSP_add_cutting_planes(struct LP *lp, struct GTSP *data);

int GTSP_write_problem(struct GTSP *data, char *filename);

int GTSP_write_solution(struct GTSP *data, char *filename, double *x);

int optimize_vertex_in_cluster(struct Tour * tour, struct GTSP *data);

int two_opt(struct Tour * tour, struct GTSP *data);

int K_opt(struct Tour* tour, struct GTSP *data);

int tour_length(int* tour, struct GTSP* data);

void print_tour(int* tour, struct GTSP* data);

int list_length(struct Tour * tour, struct GTSP* data);

void print_list(struct Tour * tour, struct GTSP* data);

int build_tour_from_x(struct GTSP *data, struct Tour *tour, double *x);

int GTSP_solution_found(struct BNC *bnc, struct GTSP *data, double *x);

int GTSP_check_solution(struct GTSP *data, double *x);

int GTSP_read_solution(struct GTSP *gtsp, char *filename, double **p_x);

#endif //_PROJECT_GTSP_H_
