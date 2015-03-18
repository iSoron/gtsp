#ifndef _PROJECT_TSP_H_
#define _PROJECT_TSP_H_

#include "lp.h"
#include "graph.h"

struct TSPData
{
    int node_count;
    int edge_count;
    int *edge_list;
    int *edge_weights;
};

int TSP_init_data(struct TSPData *data);

void TSP_free_data(struct TSPData *data);

int TSP_find_violated_subtour_elimination_cut
        (struct LP *lp, struct TSPData *data);

int TSP_is_graph_connected(
        struct Graph *G,
        double *x,
        int *island_count,
        int *island_sizes,
        int *island_start,
        int *island_nodes);

int TSP_find_closest_neighbor_tour(
        int start,
        int node_count,
        int edge_count,
        int *edges,
        int *elen,
        int *path_length);

int TSP_add_subtour_elimination_cut(struct LP *lp, int delta_length, int *delta);

int TSP_read_problem(char *filename, struct TSPData *data);

int TSP_add_cutting_planes(struct LP *lp, struct TSPData *data);

int TSP_init_lp(struct LP *lp, struct TSPData *data);

double TSP_find_initial_solution(struct TSPData *data);

#endif
