//
// Created by isoron on 19/03/15.
//

#ifndef _PROJECT_FLOW_H_
#define _PROJECT_FLOW_H_

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

#include "graph.h"

#endif //_PROJECT_FLOW_H_
