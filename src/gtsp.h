#ifndef _PROJECT_GTSP_H_
#define _PROJECT_GTSP_H_

#include "lp.h"
#include "graph.h"

struct GTSP
{
    struct Graph *graph;

    int *clusters;
    int cluster_count;
};

int GTSP_create_random_problem(
        int node_count, int cluster_count, int grid_size, struct GTSP *data);

void GTSP_free(struct GTSP *data);

int GTSP_init_data(struct GTSP *data);

int GTSP_init_lp(struct LP *lp, struct GTSP *data);

int GTSP_add_cutting_planes(struct LP *lp, struct GTSP *data);

int GTSP_write_problem(struct GTSP *data, char *filename);

int GTSP_write_solution(struct GTSP *data, char *filename, double *x);

int GTSP_main(int argc, char **argv);

extern double *OPTIMAL_X;
extern double FLOW_CPU_TIME;

#endif //_PROJECT_GTSP_H_
