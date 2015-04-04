#include "gtsp-subtour.h"
#include <assert.h>
#include "util.h"
#include "flow.h"

static void deactivate_cluster_node(double *capacities, struct Node *cluster_node)
{
    for (int i = 0; i < cluster_node->degree; i++)
    {
        struct Adjacency *adj = &cluster_node->adj[i];
        struct Edge *e = adj->edge;

        capacities[e->index] = 0;
    }
}

static void activate_cluster_node(double *capacities, struct Node *cluster_node)
{
    for (int i = 0; i < cluster_node->degree; i++)
    {
        struct Adjacency *adj = &cluster_node->adj[i];
        struct Edge *e = adj->edge;

        capacities[e->index] = 1e10;
        capacities[e->reverse->index] = 1e10;
    }
}

static int build_flow_digraph(
        struct GTSP *data, double *x, struct Graph *digraph, double *capacities)
{
    int rval = 0;

    int *digraph_edges = 0;

    int node_count = data->graph->node_count;
    struct Graph *graph = data->graph;

    int digraph_node_count = node_count + data->cluster_count;
    int digraph_edge_count = 4 * graph->edge_count + 2 * graph->node_count;

    digraph_edges = (int *) malloc(2 * digraph_edge_count * sizeof(int));
    abort_if(!digraph_edges, "could not allocate digraph_edges");

    // Create four directed edges for each edge of the original graph
    int ke = 0;
    int kc = 0;
    for (int i = 0; i < graph->edge_count; i++)
    {
        if (x[i] < EPSILON) continue;

        struct Edge *e = &graph->edges[i];
        int from = e->from->index;
        int to = e->to->index;

        digraph_edges[ke++] = from;
        digraph_edges[ke++] = to;
        capacities[kc++] = x[i];

        digraph_edges[ke++] = to;
        digraph_edges[ke++] = from;
        capacities[kc++] = 0;

        digraph_edges[ke++] = to;
        digraph_edges[ke++] = from;
        capacities[kc++] = x[i];

        digraph_edges[ke++] = from;
        digraph_edges[ke++] = to;
        capacities[kc++] = 0;
    }

    // Create an extra node for each cluster and connect it to the vertices
    // of the cluster through some edge with very high capacity
    for (int i = 0; i < node_count; i++)
    {
        struct Node *n = &graph->nodes[i];
        int cl = data->node_to_cluster[n->index];

        digraph_edges[ke++] = n->index;
        digraph_edges[ke++] = node_count + cl;
        capacities[kc++] = 0;

        digraph_edges[ke++] = node_count + cl;
        digraph_edges[ke++] = n->index;
        capacities[kc++] = 0;
    }

    assert(ke <= 2 * digraph_edge_count);
    assert(kc <= digraph_edge_count);

    digraph_edge_count = kc;

    rval = graph_build(digraph_node_count, kc, digraph_edges, 1, digraph);
    abort_if(rval, "graph_build failed");

    for (int i = 0; i < digraph_edge_count; i += 2)
    {
        digraph->edges[i].reverse = &digraph->edges[i + 1];
        digraph->edges[i + 1].reverse = &digraph->edges[i];
    }

    CLEANUP:
    if (digraph_edges) free(digraph_edges);
    return rval;
}

static int add_subtour_cut(
        struct LP *lp, struct Edge **cut_edges, int cut_edges_count)
{
    int rval = 0;

    char sense = 'G';
    double rhs = 2.0;
    int newnz = cut_edges_count;

    int *rmatind = 0;
    double *rmatval = 0;

    rmatind = (int *) malloc(newnz * sizeof(int));
    abort_if(!rmatind, "could not allocate rmatind");

    rmatval = (double *) malloc(newnz * sizeof(double));
    abort_if(!rmatval, "could not allocate rmatval");

    for (int i = 0; i < cut_edges_count; i++)
    {
        rmatind[i] = cut_edges[i]->index;
        rmatval[i] = 1.0;
    }

    log_verbose("Generated cut:\n");
    for (int i = 0; i < newnz; i++)
            log_verbose("    %8.2f x%d\n", rmatval[i], rmatind[i]);
    log_verbose("    %c %.2lf\n", sense, rhs);

    if (OPTIMAL_X)
    {
        double sum = 0;
        for (int i = 0; i < newnz; i++)
            sum += rmatval[i] * OPTIMAL_X[rmatind[i]];
        abort_if(sum <= rhs - EPSILON, "cannot add invalid cut");
    }

    struct Row *cut = 0;
    cut = (struct Row *) malloc(sizeof(struct Row));
    abort_if(!cut, "could not allocate cut");

    cut->nz = newnz;
    cut->sense = sense;
    cut->rhs = rhs;
    cut->rmatval = rmatval;
    cut->rmatind = rmatind;

    rval = LP_add_cut(lp, cut);
    abort_if(rval, "LP_add_cut failed");

    CLEANUP:
    return rval;
}

int GTSP_find_exact_subtour_cuts(struct LP *lp, struct GTSP *data)
{
    int rval = 0;

    double *x = 0;
    double *flow = 0;
    double *capacities = 0;
    struct Edge **cut_edges = 0;

    double initial_time = get_user_time();

    int added_cuts_count = 0;
    struct Graph *graph = data->graph;
    const int edge_count = graph->edge_count;
    const int node_count = graph->node_count;
    const int cluster_count = data->cluster_count;

    int num_cols = LP_get_num_cols(lp);
    int digraph_edge_count = 4 * edge_count + 2 * node_count;
    int original_cut_pool_size = lp->cut_pool_size;

    x = (double *) malloc(num_cols * sizeof(double));
    abort_if(!x, "could not allocate x");

    rval = LP_get_x(lp, x);
    abort_if(rval, "LP_get_x failed");

#if LOG_LEVEL >= LOG_LEVEL_DEBUG
    rval = GTSP_write_solution(data, "gtsp-frac.out", x);
    abort_if(rval, "GTSP_write_solution failed");
#endif

    struct Graph digraph;
    graph_init(&digraph);

    flow = (double *) malloc(digraph_edge_count * sizeof(double));
    capacities = (double *) malloc(digraph_edge_count * sizeof(double));
    cut_edges = (struct Edge **) malloc(edge_count * sizeof(struct Edge *));

    abort_if(!flow, "could not allocate flow");
    abort_if(!capacities, "could not allocate capacities");
    abort_if(!cut_edges, "could not allocate cut_edges");

    rval = build_flow_digraph(data, x, &digraph, capacities);
    abort_if(rval, "build_flow_digraph failed");

    int cuts_count = 0;

    for (int i = 0; i < cluster_count; i++)
    {
        for (int j = i + 1; j < cluster_count; j++)
        {
            struct Node *from = &digraph.nodes[node_count + i];
            struct Node *to = &digraph.nodes[node_count + j];

            double flow_value;
            int cut_edges_count;

            activate_cluster_node(capacities, from);
            activate_cluster_node(capacities, to);

            log_verbose("Sending flow from cluster %d to cluster %d\n", i, j);

            rval = flow_find_max_flow(&digraph, capacities, from, to, flow,
                    &flow_value);

            abort_if(rval, "flow_find_max_flow failed");

            log_verbose("    flow value = %.4lf\n", flow_value);

            deactivate_cluster_node(capacities, from);
            deactivate_cluster_node(capacities, to);

            if (flow_value >= 2 - EPSILON) continue;

            log_verbose("Marked nodes:\n");
            for (int k = 0; k < node_count; k++)
            {
                graph->nodes[k].mark = digraph.nodes[k].mark;
                if (digraph.nodes[k].mark) log_verbose("    %d\n", k);
            }

            rval = get_cut_edges_from_marks(graph, &cut_edges_count, cut_edges);
            abort_if(rval, "get_cut_edges_from_marks failed");

            log_verbose("Cut edges:\n");
            for (int k = 0; k < cut_edges_count; k++)
                    log_verbose("  %d %d (%d)\n", cut_edges[k]->from->index,
                            cut_edges[k]->to->index, cut_edges[k]->index);

            rval = add_subtour_cut(lp, cut_edges, cut_edges_count);
            abort_if(rval, "add_subtour_cut failed");

            cuts_count++;
        }
    }

    added_cuts_count = lp->cut_pool_size - original_cut_pool_size;
    log_debug("    %d cluster-to-cluster\n", added_cuts_count);

    SUBTOUR_COUNT += added_cuts_count;
    SUBTOUR_TIME += get_user_time() - initial_time;

    CLEANUP:
    graph_free(&digraph);
    if (capacities) free(capacities);
    if (cut_edges) free(cut_edges);
    if (flow) free(flow);
    if (x) free(x);
    return rval;
}
