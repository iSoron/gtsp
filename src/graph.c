#include <malloc.h>
#include "main.h"
#include "graph.h"
#include "util.h"
#include "lp.h"

void graph_dfs(
        int n, struct Graph *G, double *x, int *island_size, int *island_nodes)
{
    *(island_nodes + (*island_size)) = n;
    (*island_size)++;

    struct Node *pn = &G->node_list[n];
    pn->mark = 1;

    for (int i = 0; i < pn->deg; i++)
    {
        if (x[pn->adj[i].e] > LP_EPSILON)
        {
            int neighbor = pn->adj[i].n;

            if (G->node_list[neighbor].mark == 0)
                graph_dfs(neighbor, G, x, island_size, island_nodes);
        }
    }
}

void graph_init(struct Graph *G)
{
    if (!G) return;

    G->node_list = 0;
    G->adj_space = 0;
    G->node_count = 0;
    G->edge_count = 0;
}

void graph_free(struct Graph *G)
{
    if (!G) return;

    if (G->node_list) free(G->node_list);
    if (G->adj_space) free(G->adj_space);
}

int graph_build(int node_count, int edge_count, int *edge_list, struct Graph *G)
{
    int rval = 0;
    struct Node *n;
    struct AdjObj *p;

    G->node_list = (struct Node *) malloc(node_count * sizeof(struct Node));
    G->adj_space =
            (struct AdjObj *) malloc(2 * edge_count * sizeof(struct AdjObj));

    abort_if(!G->node_list, "could not allocate G->node_list");
    abort_if(!G->adj_space, "could not allocate G->adj_space");

    for (int i = 0; i < node_count; i++)
        G->node_list[i].deg = 0;

    for (int i = 0; i < edge_count; i++)
    {
        int a = edge_list[2 * i];
        int b = edge_list[2 * i + 1];
        G->node_list[a].deg++;
        G->node_list[b].deg++;
    }

    p = G->adj_space;
    for (int i = 0; i < node_count; i++)
    {
        G->node_list[i].adj = p;
        p += G->node_list[i].deg;
        G->node_list[i].deg = 0;
    }

    for (int i = 0; i < edge_count; i++)
    {
        int a = edge_list[2 * i];
        int b = edge_list[2 * i + 1];
        n = &G->node_list[a];
        n->adj[n->deg].n = b;
        n->adj[n->deg].e = i;
        n->deg++;
        n = &G->node_list[b];
        n->adj[n->deg].n = a;
        n->adj[n->deg].e = i;
        n->deg++;
    }

    G->node_count = node_count;
    G->edge_count = edge_count;

    CLEANUP:
    if (rval)
    {
        if (G->node_list) free(G->node_list);
        if (G->adj_space) free(G->adj_space);
    }
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