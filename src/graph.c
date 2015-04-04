#include <malloc.h>
#include "graph.h"
#include "util.h"

void graph_init(struct Graph *graph)
{
    if (!graph) return;

    graph->nodes = 0;
    graph->adj = 0;
    graph->node_count = 0;
    graph->edge_count = 0;
    graph->x_coordinates = 0;
    graph->y_coordinates = 0;
}

void graph_free(struct Graph *graph)
{
    if (!graph) return;

    if (graph->edges) free(graph->edges);
    if (graph->nodes) free(graph->nodes);
    if (graph->adj) free(graph->adj);
    if (graph->x_coordinates) free(graph->x_coordinates);
    if (graph->y_coordinates) free(graph->y_coordinates);
}

int graph_build(
        int node_count,
        int edge_count,
        int *edges,
        int is_directed,
        struct Graph *graph)
{
    int rval = 0;
    struct Node *n;
    struct Adjacency *p;

    graph->edges = (struct Edge *) malloc(edge_count * sizeof(struct Edge));
    graph->nodes = (struct Node *) malloc(node_count * sizeof(struct Node));
    graph->adj = (struct Adjacency *) malloc(
            2 * edge_count * sizeof(struct Adjacency));

    abort_if(!graph->edges, "could not allocate G->edges\n");
    abort_if(!graph->nodes, "could not allocate G->nodes");
    abort_if(!graph->adj, "could not allocate G->adj");

    for (int i = 0; i < node_count; i++)
    {
        graph->nodes[i].index = i;
        graph->nodes[i].degree = 0;
    }

    for (int i = 0; i < edge_count; i++)
    {
        int a = edges[2 * i];
        int b = edges[2 * i + 1];
        graph->nodes[a].degree++;
        if (!is_directed) graph->nodes[b].degree++;

        graph->edges[i].reverse = 0;
        graph->edges[i].index = i;
        graph->edges[i].from = &graph->nodes[a];
        graph->edges[i].to = &graph->nodes[b];
    }

    p = graph->adj;
    for (int i = 0; i < node_count; i++)
    {
        graph->nodes[i].adj = p;
        p += graph->nodes[i].degree;
        graph->nodes[i].degree = 0;
    }

    for (int i = 0; i < edge_count; i++)
    {
        int a = edges[2 * i];
        int b = edges[2 * i + 1];

        n = &graph->nodes[a];
        n->adj[n->degree].neighbor = &graph->nodes[b];
        n->adj[n->degree].edge = &graph->edges[i];
        n->degree++;

        if (!is_directed)
        {
            n = &graph->nodes[b];
            n->adj[n->degree].neighbor = &graph->nodes[a];
            n->adj[n->degree].edge = &graph->edges[i];
            n->degree++;
        }
    }

    graph->node_count = node_count;
    graph->edge_count = edge_count;

    CLEANUP:
    if (rval)
    {
        if (graph->edges) free(graph->edges);
        if (graph->nodes) free(graph->nodes);
        if (graph->adj) free(graph->adj);
    }
    return rval;
}

int get_cut_edges_from_marks(
        struct Graph *graph, int *cut_edges_count, struct Edge **cut_edges)
{
    *cut_edges_count = 0;

    for (int i = 0; i < graph->edge_count; ++i)
    {
        struct Edge *e = &graph->edges[i];
        struct Node *from = e->from;
        struct Node *to = e->to;
        if (from->mark != to->mark)
            cut_edges[(*cut_edges_count)++] = e;
    }

    return 0;
}

int graph_dump(const struct Graph *graph)
{
    (void) graph;
#if LOG_LEVEL >= LOG_LEVEL_DEBUG
    log_debug("node_count: %d edge_count: %d\n", graph->node_count,
            graph->edge_count);

    for (int i = 0; i < graph->node_count; i++)
    {
        struct Node *n = &graph->nodes[i];
        log_debug("%3d degree: %d mark: %d\n", n->index, n->degree, n->mark);
    }

    for (int i = 0; i < graph->edge_count; i++)
    {
        struct Edge *e = &graph->edges[i];
        log_debug("%3d (%d, %d) weight: %d ", e->index, e->from->index,
                e->to->index, e->weight);
            if (e->reverse) printf("reverse: %d ", e->reverse->index);
            printf("\n");

    }
    #endif
    return 0;
}
