//
// Created by isoron on 18/03/15.
//

#ifndef _PROJECT_GTSP_H_
#define _PROJECT_GTSP_H_

#include "lp.h"

struct Edge
{
    int from;
    int to;
    int weight;
};

struct GTSP
{
    int node_count;
    int edge_count;
    struct Edge *edges;

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

int GTSP_write_data(struct GTSP *data, char *filename);

int GTSP_write_solution(struct GTSP *data, char *filename, double *x);

#endif //_PROJECT_GTSP_H_
