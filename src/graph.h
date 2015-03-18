#ifndef __GRAPH_H_
#define __GRAPH_H_

#include "main.h"

void graph_dfs(
        int n, struct Graph *G, double *x, int *icount, int *island);

void graph_init(struct Graph *G);

void graph_free(struct Graph *G);

int graph_build(int node_count, int edge_count, int *edge_list, struct Graph *G);

int euclid_edgelen(int i, int j, double *x, double *y);

void get_delta(
        int nsize, int *nlist, int ecount, int *elist, int *deltacount,
        int *delta, int *marks);

#endif //___GRAPH_H_
