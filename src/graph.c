#include <malloc.h>
#include "graph.h"
#include "util.h"
#include "lp.h"

void graph_init(struct Graph *graph)
{
    if (!graph) return;

    graph->nodes = 0;
    graph->adj = 0;
    graph->node_count = 0;
    graph->edge_count = 0;
}

void graph_free(struct Graph *graph)
{
    if (!graph) return;

    if (graph->nodes) free(graph->nodes);
    if (graph->adj) free(graph->adj);
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
        n->adj[n->degree].neighbor_index = b;
        n->adj[n->degree].edge_index = i;
        n->adj[n->degree].neighbor = &graph->nodes[b];
        n->adj[n->degree].edge = &graph->edges[i];
        n->degree++;

        if (!is_directed)
        {
            n = &graph->nodes[b];
            n->adj[n->degree].neighbor_index = a;
            n->adj[n->degree].edge_index = i;
            n->adj[n->degree].neighbor = &graph->nodes[b];
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

void graph_dfs(
        int n, struct Graph *G, double *x, int *island_size, int *island_nodes)
{
    *(island_nodes + (*island_size)) = n;
    (*island_size)++;

    struct Node *pn = &G->nodes[n];
    pn->mark = 1;

    for (int i = 0; i < pn->degree; i++)
    {
        if (x[pn->adj[i].edge_index] > LP_EPSILON)
        {
            int neighbor = pn->adj[i].neighbor_index;

            if (G->nodes[neighbor].mark == 0)
                graph_dfs(neighbor, G, x, island_size, island_nodes);
        }
    }
}

int graph_build_directed_from_undirected(
        const struct Graph *graph, struct Graph *digraph)
{
    int rval = 0;

    int *edges = 0;

    edges = (int *) malloc(4 * graph->edge_count * sizeof(int));
    abort_if(!edges, "could not allocate edges");

    for (int i = 0; i < graph->edge_count; i++)
    {
        struct Edge *e = &graph->edges[i];
        edges[4 * i] = edges[4 * i + 3] = e->from->index;
        edges[4 * i + 1] = edges[4 * i + 2] = e->to->index;
    }

    rval = graph_build(graph->node_count, 2 * graph->edge_count, edges, 1,
            digraph);
    abort_if(rval, "graph_build failed");

    for (int i = 0; i < graph->edge_count; i++)
    {
        digraph->edges[2 * i].reverse = &digraph->edges[i * 2 + 1];
        digraph->edges[2 * i + 1].reverse = &digraph->edges[i * 2];
    }

    CLEANUP:
    if (!edges) free(edges);
    return rval;
}

void get_delta(
        int island_node_count,
        int *island_nodes,
        int edge_count,
        int *edges,
        int *delta_count,
        int *delta,
        int *marks)
{
    for (int i = 0; i < island_node_count; i++)
        marks[island_nodes[i]] = 1;

    int k = 0;
    for (int i = 0; i < edge_count; i++)
        if (marks[edges[2 * i]] + marks[edges[2 * i + 1]] == 1)
            delta[k++] = i;

    *delta_count = k;

    for (int i = 0; i < island_node_count; i++)
        marks[island_nodes[i]] = 0;
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
        if (from->mark && !to->mark)
            cut_edges[(*cut_edges_count)++] = e;
    }

    return 0;
}
