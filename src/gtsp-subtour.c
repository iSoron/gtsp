#include "gtsp-subtour.h"
#include <assert.h>
#include <float.h>
#include "util.h"
#include "flow.h"

extern double FLOW_CPU_TIME;

int static build_flow_digraph(
        struct GTSP *data, double *x, struct Graph *digraph, double *capacities)
{
    int rval = 0;

    int *digraph_edges = 0;

    int node_count = data->graph->node_count;
    struct Graph *graph = data->graph;

    int digraph_node_count = node_count + data->cluster_count + 1;
    int digraph_edge_count = 4 * graph->edge_count + 2 * graph->node_count +
            2 * data->cluster_count;

    digraph_edges = (int *) malloc(2 * digraph_edge_count * sizeof(int));
    abort_if(!digraph_edges, "could not allocate digraph_edges");

    // Create four directed edges for each edge of the original graph
    int ke = 0;
    int kc = 0;
    for (int i = 0; i < graph->edge_count; i++)
    {
        if (x[node_count + i] < LP_EPSILON) continue;

        struct Edge *e = &graph->edges[i];
        int from = e->from->index;
        int to = e->to->index;

        digraph_edges[ke++] = from;
        digraph_edges[ke++] = to;
        capacities[kc++] = x[node_count + i];

        digraph_edges[ke++] = to;
        digraph_edges[ke++] = from;
        capacities[kc++] = 0;

        digraph_edges[ke++] = to;
        digraph_edges[ke++] = from;
        capacities[kc++] = x[node_count + i];

        digraph_edges[ke++] = from;
        digraph_edges[ke++] = to;
        capacities[kc++] = 0;
    }

    // Create an extra node for each cluster and connect it to the vertices
    // of the cluster through some edge with very high capacity
    for (int i = 0; i < node_count; i++)
    {
        struct Node *n = &graph->nodes[i];
        int cl = data->clusters[n->index];

        digraph_edges[ke++] = n->index;
        digraph_edges[ke++] = node_count + cl;
        capacities[kc++] = 0;

        digraph_edges[ke++] = node_count + cl;
        digraph_edges[ke++] = n->index;
        capacities[kc++] = 0;
    }

    // Create an extra root node and connect it to each cluster node through
    // some edge with zero capacity
    for (int i = 0; i < data->cluster_count; i++)
    {
        digraph_edges[ke++] = node_count + i;
        digraph_edges[ke++] = node_count + data->cluster_count;
        capacities[kc++] = 0;

        digraph_edges[ke++] = node_count + data->cluster_count;
        digraph_edges[ke++] = node_count + i;
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

int static add_subtour_cut(
        struct LP *lp,
        struct Graph *graph,
        struct Node *from,
        struct Node *to,
        struct Edge **cut_edges,
        int cut_edges_count,
        int type)
{
    int rval = 0;

    char sense = 'G';
    double rhs = 2.0 - 2.0 * type;
    int newnz = cut_edges_count + type;

    int *rmatind = 0;
    double *rmatval = 0;

    rmatind = (int *) malloc(newnz * sizeof(int));
    abort_if(!rmatind, "could not allocate rmatind");

    rmatval = (double *) malloc(newnz * sizeof(double));
    abort_if(!rmatval, "could not allocate rmatval");

    for (int i = 0; i < cut_edges_count; i++)
    {
        rmatind[i] = cut_edges[i]->index + graph->node_count;
        rmatval[i] = 1.0;
    }

    if (type >= 1)
    {
        rmatind[cut_edges_count] = from->index;
        rmatval[cut_edges_count] = -2.0;
    }

    if (type >= 2)
    {
        rmatind[cut_edges_count + 1] = to->index;
        rmatval[cut_edges_count + 1] = -2.0;
    }

    log_verbose("Generated cut:\n");
    for (int i = 0; i < newnz; i++)
            log_verbose("%8.2f x%d\n", rmatval[i], rmatind[i]);
    log_verbose("    %c %.2lf\n", sense, rhs);

    if (OPTIMAL_X)
    {
        double sum = 0;
        for (int i = 0; i < newnz; i++)
            sum += rmatval[i] * OPTIMAL_X[rmatind[i]];
        abort_if(sum <= rhs - LP_EPSILON, "cannot add invalid cut");
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

int find_exact_subtour_cuts(
        struct LP *lp, struct GTSP *data, double min_cut_violation)
{
    int rval = 0;

    double *x = 0;
    double *capacities = 0;

    int added_cuts_count = 0;
    struct Graph *graph = data->graph;

    int num_cols = LP_get_num_cols(lp);
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
    int digraph_edge_count = 4 * graph->edge_count + 2 * graph->node_count +
            2 * data->cluster_count;

    int original_cut_pool_size = lp->cut_pool_size;

    capacities = (double *) malloc(digraph_edge_count * sizeof(double));
    abort_if(!capacities, "could not allocate capacities");

    rval = build_flow_digraph(data, x, &digraph, capacities);
    abort_if(rval, "build_flow_digraph failed");

    // Constraints (2.1)
    rval = find_exact_subtour_cuts_cluster_to_cluster(lp, data, &digraph,
            capacities, min_cut_violation);
    abort_if(rval, "find_exact_subtour_cuts_cluster_to_cluster failed");

    added_cuts_count = lp->cut_pool_size - original_cut_pool_size;
    if (added_cuts_count > 0)
    {
        log_debug("Added %d cluster-to-cluster subtour cuts\n",
                added_cuts_count);
        goto CLEANUP;
    }

    // Constraints (2.2)
    original_cut_pool_size = lp->cut_pool_size;
    rval = find_exact_subtour_cuts_node_to_cluster(lp, data, x, &digraph,
            capacities, min_cut_violation);
    abort_if(rval, "find_exact_subtour_cuts_node_to_cluster failed");

    added_cuts_count = lp->cut_pool_size - original_cut_pool_size;
    if (added_cuts_count > 0)
    {
        log_debug("Added %d node-to-cluster subtour cuts\n", added_cuts_count);
        goto CLEANUP;
    }

    // Constraints (2.3)
    original_cut_pool_size = lp->cut_pool_size;
    rval = find_exact_subtour_cuts_node_to_node(lp, data, x, &digraph,
            capacities, min_cut_violation);
    abort_if(rval, "find_exact_subtour_cuts_node_to_node failed");

    added_cuts_count = lp->cut_pool_size - original_cut_pool_size;
    if (added_cuts_count > 0)
    {
        log_debug("Added %d node-to-node subtour cuts\n", added_cuts_count);
        goto CLEANUP;
    }

    CLEANUP:
    graph_free(&digraph);
    if (capacities) free(capacities);
    if (x) free(x);
    return rval;
}

int find_exact_subtour_cuts_node_to_node(
        struct LP *lp,
        struct GTSP *data,
        double *x,
        struct Graph *digraph,
        double *capacities,
        double min_cut_violation)
{
    int rval = 0;

    struct Edge **cut_edges = 0;
    double *flow = 0;

    struct Graph *graph = data->graph;
    int *clusters = data->clusters;

    cut_edges = (struct Edge **) malloc(
            graph->edge_count * sizeof(struct Edge *));
    flow = (double *) malloc(digraph->edge_count * sizeof(double));

    abort_if(!cut_edges, "could not allocate cut_edges");
    abort_if(!flow, "could not allocate flow");

    int max_x_index = 0;
    double max_x = DBL_MIN;

    for (int i = 0; i < graph->node_count; i++)
    {
        struct Node *n = &graph->nodes[i];
        if (x[n->index] > max_x)
        {
            max_x = x[n->index];
            max_x_index = i;
        }
    }

    int i = max_x_index;

    for (int j = 0; j < graph->node_count; j++)
    {
        if (i == j) continue;
        if (clusters[i] == clusters[j]) continue;
        if (x[i] + x[j] - 1 <= LP_EPSILON) continue;

        struct Node *from = &digraph->nodes[i];
        struct Node *to = &digraph->nodes[j];

        int cut_edges_count;
        double flow_value;

        rval = flow_find_max_flow(digraph, capacities, from, to, flow,
                &flow_value);
        abort_if(rval, "flow_find_max_flow failed");

        if (flow_value >= 2 * (x[i] + x[j] - 1) - min_cut_violation)
            continue;

        log_verbose("Marked nodes:\n");
        for (int k = 0; k < graph->node_count; k++)
        {
            graph->nodes[k].mark = digraph->nodes[k].mark;
            if (digraph->nodes[k].mark) log_verbose("    %d\n", k);
        }

        rval = get_cut_edges_from_marks(graph, &cut_edges_count, cut_edges);
        abort_if(rval, "get_cut_edges_from_marks failed");

        log_verbose("Cut edges:\n");
        for (int k = 0; k < cut_edges_count; k++)
                log_verbose("  %d %d (%d)\n", cut_edges[k]->from->index,
                        cut_edges[k]->to->index, cut_edges[k]->index);

        rval = add_subtour_cut(lp, graph, from, to, cut_edges, cut_edges_count,
                2);
        abort_if(rval, "add_subtour_cut failed");
    }

    CLEANUP:
    if (flow) free(flow);
    if (cut_edges) free(cut_edges);
    return rval;
}

int find_exact_subtour_cuts_node_to_cluster(
        struct LP *lp,
        struct GTSP *data,
        double *x,
        struct Graph *digraph,
        double *capacities,
        double min_cut_violation)
{
    int rval = 0;

    int cuts_count = 0;
    struct Edge **cut_edges = 0;
    double *flow = 0;

    struct Graph *graph = data->graph;
    int *clusters = data->clusters;

    cut_edges = (struct Edge **) malloc(
            graph->edge_count * sizeof(struct Edge *));
    flow = (double *) malloc(digraph->edge_count * sizeof(double));
    abort_if(!cut_edges, "could not allocate cut_edges");
    abort_if(!flow, "could not allocate flow");

    for (int i = 0; i < graph->node_count; i++)
    {
        for (int j = 0; j < data->cluster_count; j++)
        {
            if (clusters[i] == j) continue;
            if (x[i] < LP_EPSILON) continue;

            struct Node *from = &digraph->nodes[i];
            struct Node *to = &digraph->nodes[graph->node_count + j];

            log_verbose(
                    "Sending flow from node %d to cluster %d (must be >= %.4lf)\n",
                    i, j, 2 * x[i]);

            activate_cluster_node(capacities, to);

            double flow_value;
            int cut_edges_count;

            rval = flow_find_max_flow(digraph, capacities, from, to, flow,
                    &flow_value);
            abort_if(rval, "flow_find_max_flow failed");

            log_verbose("    flow value = %.4lf\n", flow_value);

            deactivate_cluster_node(capacities, to);

            if (flow_value + min_cut_violation >= 2 * x[i])
                continue;

            log_verbose("Marked nodes:\n");
            for (int k = 0; k < graph->node_count; k++)
            {
                graph->nodes[k].mark = digraph->nodes[k].mark;
                if (graph->nodes[k].mark) log_verbose("    %d\n", k);
            }

            rval = get_cut_edges_from_marks(graph, &cut_edges_count, cut_edges);
            abort_if(rval, "get_cut_edges_from_marks failed");

            log_verbose("Cut edges:\n");
            for (int k = 0; k < cut_edges_count; k++)
            {
                struct Edge *e = cut_edges[k];
                assert(e->from->mark != e->to->mark);
                log_verbose("  %d (%d) %d (%d) [%d]\n", e->from->index,
                        e->from->mark, e->to->index, e->to->mark, e->index);
            }

            rval = add_subtour_cut(lp, graph, from, 0, cut_edges,
                    cut_edges_count, 1);
            abort_if(rval, "add_subtour_cut failed");

            cuts_count++;
        }
    }

    CLEANUP:
    if (cut_edges) free(cut_edges);
    if (flow) free(flow);
    return rval;
}

int find_exact_subtour_cuts_cluster_to_cluster(
        struct LP *lp,
        struct GTSP *data,
        struct Graph *digraph,
        double *capacities,
        double min_cut_violation)
{
    int rval = 0;

    double *flow = 0;
    struct Edge **cut_edges = 0;

    int cuts_count = 0;

    struct Graph *graph = data->graph;

    cut_edges = (struct Edge **) malloc(
            graph->edge_count * sizeof(struct Edge *));
    flow = (double *) malloc(digraph->edge_count * sizeof(double));
    abort_if(!cut_edges, "could not allocate cut_edges");
    abort_if(!flow, "could not allocate flow");

    struct Node *root_node = &digraph->nodes[graph->node_count +
            data->cluster_count];

    for (int i = 0; i < data->cluster_count; i++)
    {
        for (int j = i + 1; j < data->cluster_count; j++)
        {
            struct Node *from = &digraph->nodes[graph->node_count + i];
            struct Node *to = &digraph->nodes[graph->node_count + j];

            double flow_value;
            int cut_edges_count;

            activate_cluster_node(capacities, from);
            activate_cluster_node(capacities, to);
            deactivate_cluster_node(capacities, root_node);

            log_verbose("Sending flow from cluster %d to cluster %d\n", i, j);

            rval = flow_find_max_flow(digraph, capacities, from, to, flow,
                    &flow_value);

            abort_if(rval, "flow_find_max_flow failed");

            log_verbose("    flow value = %.4lf\n", flow_value);

            deactivate_cluster_node(capacities, from);
            deactivate_cluster_node(capacities, to);

            if (flow_value >= 2 - min_cut_violation) continue;

            log_verbose("Marked nodes:\n");
            for (int k = 0; k < graph->node_count; k++)
            {
                graph->nodes[k].mark = digraph->nodes[k].mark;
                if (digraph->nodes[k].mark) log_verbose("    %d\n", k);
            }

            rval = get_cut_edges_from_marks(graph, &cut_edges_count, cut_edges);
            abort_if(rval, "get_cut_edges_from_marks failed");

            log_verbose("Cut edges:\n");
            for (int k = 0; k < cut_edges_count; k++)
                    log_verbose("  %d %d (%d)\n", cut_edges[k]->from->index,
                            cut_edges[k]->to->index, cut_edges[k]->index);

            rval = add_subtour_cut(lp, graph, 0, 0, cut_edges, cut_edges_count,
                    0);
            abort_if(rval, "add_subtour_cut failed");

            cuts_count++;
        }
    }

    CLEANUP:
    if (cut_edges) free(cut_edges);
    if (flow) free(flow);
    return rval;
}

void deactivate_cluster_node(double *capacities, struct Node *cluster_node)
{
    for (int i = 0; i < cluster_node->degree; i++)
    {
        struct Adjacency *adj = &cluster_node->adj[i];
        struct Edge *e = adj->edge;

        capacities[e->index] = 0;
    }
}

void activate_cluster_node(double *capacities, struct Node *cluster_node)
{
    for (int i = 0; i < cluster_node->degree; i++)
    {
        struct Adjacency *adj = &cluster_node->adj[i];
        struct Edge *e = adj->edge;

        capacities[e->index] = 1e10;
        capacities[e->reverse->index] = 1e10;
    }
}