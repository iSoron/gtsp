//
// Created by isoron on 3/17/15.
//

#ifndef ___TSP_H_
#define ___TSP_H_

#include "lp.h"
#include "graph.h"

int add_all_subtours(int ncount, int ecount, int *elist, struct LP *lp);

int TSP_find_violated_subtour_elimination_cut
        (int ncount, int ecount, int *elist, struct LP *lp);

int TSP_is_graph_connected(
        struct Graph *G, double *x, int *island_count, int *island_sizes,
        int *island_start, int *island_nodes);

int TSP_find_closest_neighbor_tour(
        int start, int node_count, int edge_count, int *edges, int *elen,
        int *path_length);

int TSP_add_subtour_elimination_cut(struct LP *lp, int deltacount, int *delta);

int TSP_read_problem(
        char *filename, int *p_ncount, int *p_ecount, int **p_elist,
        int **p_elen);

int TSP_add_cutting_planes(int ncount, int ecount, int *elist, struct LP *lp);

int TSP_init_lp(
        int node_count, struct LP *lp, int edge_count, int *edge_weights,
        int *edge_list);

double TSP_find_initial_solution
        (int *edge_weights, int *edge_list, int node_count, int edge_count);

#endif //_PROJECT_TSP_H_
