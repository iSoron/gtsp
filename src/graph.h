#ifndef _PROJECT_GRAPH_H_
#define _PROJECT_GRAPH_H_

#include "main.h"

struct AdjObj
{
    /* Index of neighbor node */
    int n;

    /* Index of adj joining neighbor */
    int e;
};

struct Node
{
    int deg;
    struct AdjObj *adj;
    int mark;
};

struct Graph
{
    int node_count;
    int edge_count;
    struct Node *node_list;
    struct AdjObj *adj_space;
};

void graph_dfs(
        int n, struct Graph *G, double *x, int *icount, int *island);

void graph_init(struct Graph *G);

void graph_free(struct Graph *G);

int graph_build
        (int node_count, int edge_count, int *edge_list, struct Graph *G);

void get_delta(
        int nsize,
        int *nlist,
        int ecount,
        int *elist,
        int *deltacount,
        int *delta,
        int *marks);

#endif
