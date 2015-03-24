//
// Created by isoron on 18/03/15.
//

#ifndef _PROJECT_GTSP_H_
#define _PROJECT_GTSP_H_

#include "lp.h"
#include "graph.h"

struct GTSP
{
    struct Graph *graph;

    int *clusters;
    int cluster_count;

    double *x_coordinates;
    double *y_coordinates;
};

int GTSP_create_random_problem(
        int node_count, int cluster_count, int grid_size, struct GTSP *data);

void GTSP_free(struct GTSP *data);

int GTSP_init_data(struct GTSP *data);

int GTSP_init_lp(struct LP *lp, struct GTSP *data);

int GTSP_write_input_data(struct GTSP *data, char *filename);

int GTSP_write_solution(struct GTSP *data, char *filename, double *x);

int GTSP_add_cutting_planes(struct LP *lp, struct GTSP *data);

int GTSP_main(int argc, char **argv);

#endif //_PROJECT_GTSP_H_
