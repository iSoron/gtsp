#include <malloc.h>
#include <float.h>
#include "gtsp.h"
#include "flow.h"
#include "util.h"

int flow_find_max_flow(
        const struct Graph *digraph,
        const double *capacities,
        struct Node *from,
        struct Node *to,
        double *flow,
        double *value)
{
    int rval = 0;

    int path_length;
    struct Edge **path_edges = 0;
    double path_capacity;

    double *residual_caps = 0;

    residual_caps = (double *) malloc(digraph->edge_count * sizeof(double));
    abort_if(!residual_caps, "could not allocate residual_caps");

    path_edges = (struct Edge **) malloc(
            digraph->edge_count * sizeof(struct Edge *));
    abort_if(!path_edges, "could not allocate path_edges");

    for (int i = 0; i < digraph->edge_count; i++)
    {
        flow[i] = 0;
        residual_caps[i] = capacities[i];
        abort_if(!digraph->edges[i].reverse,
                "digraph must have reverse edge information");
    }

    *value = 0;

    while (1)
    {
        flow_find_augmenting_path(digraph, residual_caps, from, to,
                &path_length, path_edges, &path_capacity);

        if (path_length == 0) break;

        (*value) += path_capacity;

        for (int i = 0; i < path_length; i++)
        {
            struct Edge *e = &digraph->edges[path_edges[i]->index];

            residual_caps[e->index] -= path_capacity;
            residual_caps[e->reverse->index] += path_capacity;

            flow[e->index] += path_capacity;
            flow[e->reverse->index] -= path_capacity;
        }
    }

    CLEANUP:
    if (path_edges) free(path_edges);
    if (residual_caps) free(residual_caps);
    return rval;
}

int flow_find_augmenting_path(
        const struct Graph *graph,
        const double *residual_caps,
        struct Node *from,
        struct Node *to,
        int *path_length,
        struct Edge **path_edges,
        double *path_capacity)
{
    int rval = 0;

    struct Node **queue = 0;
    int queue_start = 0;
    int queue_end = 0;

    struct Node **parents = 0;
    struct Edge **parent_edges = 0;

    int node_count = graph->node_count;

    queue = (struct Node **) malloc(node_count * sizeof(struct Node *));
    parents = (struct Node **) malloc(node_count * sizeof(struct Node *));
    parent_edges = (struct Edge **) malloc(node_count * sizeof(struct Edge *));

    abort_if(!queue, "could not allocate queue");
    abort_if(!parents, "could not allocate parents");
    abort_if(!parent_edges, "could not allocate parent_edges");

    for (int i = 0; i < node_count; i++)
        graph->nodes[i].mark = 0;

    int found = 0;
    queue[queue_end++] = from;

    while (queue_end > queue_start)
    {
        struct Node *n = queue[queue_start++];

        n->mark = 2;

        for (int i = 0; i < n->degree; i++)
        {
            struct Node *neighbor = n->adj[i].neighbor;
            struct Edge *edge = n->adj[i].edge;

            if (neighbor->mark > 0) continue;
            if (residual_caps[edge->index] < LP_EPSILON) continue;

            parents[neighbor->index] = n;
            parent_edges[neighbor->index] = edge;

            queue[queue_end++] = neighbor;
            neighbor->mark = 1;

            if (neighbor == to)
            {
                found = 1;
                break;
            }
        }

        if (found) break;
    }

    *path_length = 0;
    *path_capacity = DBL_MAX;

    if (queue_end == queue_start) goto CLEANUP;

    struct Node *n = to;
    while (n != from)
    {
        struct Edge *edge = parent_edges[n->index];
        path_edges[*path_length] = edge;

        double c = residual_caps[edge->index];
        if (c < *path_capacity) *path_capacity = c;

        n = parents[n->index];
        (*path_length)++;
    }

    CLEANUP:
    if (parents) free(parents);
    if (parent_edges) free(parent_edges);
    if (queue) free(queue);
    return rval;
}