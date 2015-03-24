#ifndef _PROJECT_FLOW_H_
#define _PROJECT_FLOW_H_

#include "graph.h"

int flow_find_augmenting_path(
        const struct Graph *graph,
        const double *residual_caps,
        struct Node *from,
        struct Node *to,
        int *path_length,
        struct Edge **path_edges,
        double *path_capacity);

int flow_find_max_flow(
        const struct Graph *digraph,
        const double *capacities,
        struct Node *from,
        struct Node *to,
        double *flow,
        double *value);

int flow_mark_reachable_nodes(
        const struct Graph *graph, double *residual_caps, struct Node *from);

int flow_main(int argc, char **argv);

#include "graph.h"

#endif //_PROJECT_FLOW_H_
