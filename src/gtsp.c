#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <getopt.h>
#include "gtsp.h"
#include "geometry.h"
#include "util.h"
#include "flow.h"
#include "branch_and_cut.h"

static double *OPTIMAL_X = 0;

int GTSP_init_data(struct GTSP *data)
{
    int rval = 0;

    data->clusters = 0;
    data->cluster_count = 0;
    data->x_coordinates = 0;
    data->y_coordinates = 0;

    data->graph = (struct Graph *) malloc(sizeof(struct Graph));
    abort_if(!data->graph, "could not allocate data->graph");

    graph_init(data->graph);

    CLEANUP:
    return rval;
}

void GTSP_free(struct GTSP *data)
{
    if (!data) return;

    graph_free(data->graph);
    free(data->graph);

    if (data->clusters) free(data->clusters);
    if (data->x_coordinates) free(data->x_coordinates);
    if (data->y_coordinates) free(data->y_coordinates);
}

int GTSP_create_random_problem(
        int node_count, int cluster_count, int grid_size, struct GTSP *data)
{
    int rval = 0;

    int *edges = 0;
    int *weights = 0;
    int *clusters = 0;

    double *x_coords = 0;
    double *y_coords = 0;

    struct Graph *graph = 0;

    int edge_count = (node_count * (node_count - 1)) / 2;

    graph = (struct Graph *) malloc(sizeof(struct Graph));
    abort_if(!graph, "could not allocate graph\n");

    graph_init(graph);

    edges = (int *) malloc(2 * edge_count * sizeof(int));
    weights = (int *) malloc(edge_count * sizeof(int));
    clusters = (int *) malloc(node_count * sizeof(int));

    abort_if(!edges, "could not allocate data->edges\n");
    abort_if(!weights, "could not allocate weights\n");
    abort_if(!clusters, "could not allocate clusters\n");

    x_coords = (double *) malloc(node_count * sizeof(double));
    y_coords = (double *) malloc(node_count * sizeof(double));

    abort_if(!x_coords, "could not allocate x_coords\n");
    abort_if(!y_coords, "could not allocate y_coords\n");

    rval = generate_random_clusters_2d(node_count, cluster_count, grid_size,
            x_coords, y_coords, clusters);
    abort_if(rval, "generate_random_clusters_2d failed");

    int curr_edge = 0;
    for (int i = 0; i < edge_count; i++)
        for (int j = i + 1; j < node_count; j++)
        {
            edges[curr_edge * 2] = i;
            edges[curr_edge * 2 + 1] = j;
            weights[curr_edge] =
                    get_euclidean_distance(x_coords, y_coords, i, j);

            curr_edge++;
        }

    rval = graph_build(node_count, edge_count, edges, 0, graph);
    abort_if(rval, "graph_build failed");

    for (int i = 0; i < edge_count; i++)
        graph->edges[i].weight = weights[i];

    data->graph = graph;
    data->clusters = clusters;
    data->cluster_count = cluster_count;
    data->x_coordinates = x_coords;
    data->y_coordinates = y_coords;

    CLEANUP:
    if (weights) free(weights);
    if (edges) free(edges);
    if (rval)
    {
        if (clusters) free(clusters);
    }
    return rval;
}

int GTSP_init_lp(struct LP *lp, struct GTSP *data)
{
    int rval = 0;

    int node_count = data->graph->node_count;
    int edge_count = data->graph->edge_count;
    int cluster_count = data->cluster_count;
    int *clusters = data->clusters;
    struct Edge *edges = data->graph->edges;

    for (int i = 0; i < node_count; i++)
    {
        rval = LP_new_row(lp, 'E', 0.0);
        abort_if(rval, "LP_new_row failed");
    }

    for (int i = 0; i < cluster_count; i++)
    {
        rval = LP_new_row(lp, 'E', 1.0);
        abort_if(rval, "LP_new_row failed");
    }

    double lb = 0.0;
    double ub = 1.0;
    int cmatbeg = 0;

    for (int i = 0; i < node_count; i++)
    {
        double obj = 0.0;
        double cmatval[] = {-2.0, 1.0};
        int cmatind[] = {i, node_count + clusters[i]};

        rval = LP_add_cols(lp, 1, 2, &obj, &cmatbeg, cmatind, cmatval, &lb,
                &ub);
        abort_if(rval, "LP_add_cols failed");
    }

    for (int i = 0; i < edge_count; i++)
    {
        double obj = (double) edges[i].weight;
        double cmatval[] = {1.0, 1.0};
        int cmatind[] = {edges[i].from->index, edges[i].to->index};

        rval = LP_add_cols(lp, 1, 2, &obj, &cmatbeg, cmatind, cmatval, &lb,
                &ub);
        abort_if(rval, "LP_add_cols failed");
    }

    CLEANUP:
    return rval;
}

int GTSP_add_subcluster_cut(
        struct LP *lp,
        struct Graph *graph,
        struct Edge *e,
        struct Edge **cut_edges,
        int cut_edges_count)
{
    int rval = 0;

    char sense = 'G';
    double rhs = 0.0;
    int newnz = cut_edges_count + 1;

    int rmatbeg = 0;
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

    rmatind[cut_edges_count] = graph->node_count + e->index;
    rmatval[cut_edges_count] = -1.0;

    log_debug("Generated cut:\n");
    for (int i = 0; i < newnz; i++)
            log_debug("%8.2f x%d\n", rmatval[i], rmatind[i]);
    log_debug("    %c %.2lf\n", sense, rhs);

    if (OPTIMAL_X)
    {
        double sum = 0;
        for (int i = 0; i < newnz; i++)
            sum += rmatval[i] * OPTIMAL_X[rmatind[i]];
        abort_if(sum <= rhs - LP_EPSILON, "cannot add invalid cut");
    }

    rval = LP_add_rows(lp, 1, newnz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
    abort_if(rval, "LP_add_rows failed");

    CLEANUP:
    if (rmatval) free(rmatval);
    if (rmatind) free(rmatind);
    return rval;
}

int GTSP_add_subtour_elimination_cut(
        struct LP *lp,
        struct Graph *graph,
        struct Node *from,
        struct Node *to,
        struct Edge **cut_edges,
        int cut_edges_count)
{
    int rval = 0;

    char sense = 'G';
    double rhs = -2.0;
    int newnz = cut_edges_count + 2;

    int rmatbeg = 0;
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

    rmatind[cut_edges_count] = from->index;
    rmatval[cut_edges_count] = -2.0;

    rmatind[cut_edges_count + 1] = to->index;
    rmatval[cut_edges_count + 1] = -2.0;

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

    rval = LP_add_rows(lp, 1, newnz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
    abort_if(rval, "LP_add_rows failed");

    CLEANUP:
    if (rmatval) free(rmatval);
    if (rmatind) free(rmatind);
    return rval;
}

int GTSP_add_subtour_elimination_cut_2(
        struct LP *lp,
        struct Graph *graph,
        struct Node *from,
        struct Node *to,
        struct Edge **cut_edges,
        int cut_edges_count)
{
    int rval = 0;

    char sense = 'G';
    double rhs = 0.0;
    int newnz = cut_edges_count + 1;

    int rmatbeg = 0;
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

    rmatind[cut_edges_count] = from->index;
    rmatval[cut_edges_count] = -2.0;

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

    rval = LP_add_rows(lp, 1, newnz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
    abort_if(rval, "LP_add_rows failed");

    CLEANUP:
    if (rmatval) free(rmatval);
    if (rmatind) free(rmatind);
    return rval;
}

int GTSP_add_subtour_elimination_cut_3(
        struct LP *lp,
        struct Graph *graph,
        struct Node *from,
        struct Node *to,
        struct Edge **cut_edges,
        int cut_edges_count)
{
    int rval = 0;

    char sense = 'G';
    double rhs = 2.0;
    int newnz = cut_edges_count;

    int rmatbeg = 0;
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

    rval = LP_add_rows(lp, 1, newnz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
    abort_if(rval, "LP_add_rows failed");

    CLEANUP:
    if (rmatval) free(rmatval);
    if (rmatind) free(rmatind);
    return rval;
}

int GTSP_build_flow_digraph(struct GTSP *data, double *x, struct Graph *digraph, double *capacities)
{
    int rval = 0;

    int *digraph_edges = 0;
    struct Graph *graph = data->graph;

    int node_count = graph->node_count;

    int digraph_edge_count = 4 * graph->edge_count + 2 * graph->node_count +
            2 * data->cluster_count;
    int digraph_node_count = node_count + data->cluster_count + 1;

    digraph_edges = (int *) malloc(2 * digraph_edge_count * sizeof(int));

    // Create four directed edges for each edge of the original graph
    int ke = 0;
    int kc = 0;
    for (int i = 0; i < graph->edge_count; i++)
    {
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
        capacities[kc++] = 1e10;

        digraph_edges[ke++] = node_count + cl;
        digraph_edges[ke++] = n->index;
        capacities[kc++] = 1e10;
    }

    // Create an extra node and connect it to each cluster node through
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

    assert(ke == 2 * digraph_edge_count);
    assert(kc == digraph_edge_count);

    rval = graph_build(digraph_node_count, digraph_edge_count, digraph_edges, 1,
            digraph);
    abort_if(rval, "graph_build failed");

    for (int i = 0; i < digraph_edge_count; i += 2)
    {
        digraph->edges[i].reverse = &digraph->edges[i + 1];
        digraph->edges[i + 1].reverse = &digraph->edges[i];
    }


    CLEANUP:
    return rval;
}

int GTSP_find_exact_subtour_elimination_cuts(
        struct LP *lp, struct GTSP *data, int *added_cuts_count)
{
    int rval = 0;

    double *x = 0;
    double *flow = 0;
    double *capacities = 0;

    struct Edge **cut_edges = 0;

    int *clusters = data->clusters;
    struct Graph *graph = data->graph;
    int node_count = graph->node_count;


    int num_cols = LP_get_num_cols(lp);
    x = (double *) malloc(num_cols * sizeof(double));
    abort_if(!x, "could not allocate x");

    rval = LP_get_x(lp, x);
    abort_if(rval, "LP_get_x failed");

    struct Graph digraph;
    graph_init(&digraph);

    int digraph_edge_count = 4 * graph->edge_count + 2 * graph->node_count +
            2 * data->cluster_count;

    flow = (double *) malloc(digraph_edge_count * sizeof(double));
    capacities = (double *) malloc(digraph_edge_count * sizeof(double));
    cut_edges =
            (struct Edge **) malloc(digraph_edge_count * sizeof(struct Edge *));

    abort_if(!flow, "could not allocate flow");
    abort_if(!capacities, "could not allocate capacities");
    abort_if(!cut_edges, "could not allocate cut_edges");


    rval = GTSP_build_flow_digraph(data, x, &digraph, capacities);
    abort_if(rval, "GTSP_build_flow_digraph failed");


    // Constraints (2.3)
    {
        int max_x_index = 0;
        double max_x = DBL_MIN;

        for (int i = 0; i < node_count; i++)
        {
            struct Node *n = &graph->nodes[i];
            if (x[n->index] > max_x)
            {
                max_x = x[n->index];
                max_x_index = i;
            }
        }

        int i = max_x_index;

        for (int j = 0; j < node_count; j++)
        {
            if (i == j) continue;

            if (clusters[i] == clusters[j]) continue;
            if (x[i] + x[j] - 1 <= LP_EPSILON) continue;

            struct Node *from = &digraph.nodes[i];
            struct Node *to = &digraph.nodes[j];

            log_verbose("Calculating max flow from node %d to node %to\n",
                    from->index, to->index);
            double flow_value;
            rval = flow_find_max_flow(&digraph, capacities, from, to, flow,
                    &flow_value);
            abort_if(rval, "flow_find_max_flow failed");

            log_verbose("    %.2lf\n", flow_value);

            if (flow_value >= 2 * (x[i] + x[j] - 1) - LP_EPSILON) continue;


            log_verbose("violation: %.2lf >= %.2lf\n", flow_value,
                    2 * (x[i] + x[j] - 1));

            int cut_edges_count;
            rval = get_cut_edges_from_marks(&digraph, &cut_edges_count,
                    cut_edges);
            abort_if(rval, "get_cut_edges_from_marks failed");

            log_verbose("Adding cut for i=%d j=%d, cut edges:\n", i, j);
            int c = 0;
            for (int k = 0; k < cut_edges_count / 2; k++)
            {
                int idx = cut_edges[k * 2]->index / 4;
                if (idx > graph->edge_count) continue;

                cut_edges[c++] = &graph->edges[idx];
                log_verbose("    %d %d\n", cut_edges[c - 1]->from->index,
                        cut_edges[c - 1]->to->index);
            }

            rval = GTSP_add_subtour_elimination_cut(lp, graph, from, to,
                    cut_edges, c);
            abort_if(rval, "GTSP_add_subtour_elimination_cut failed");

            (*added_cuts_count)++;
            if (*added_cuts_count > 10) goto CLEANUP;
        }
    }

    // Constraints (2.2)
    for (int i = 0; i < node_count; i++)
    {
        for (int j = 0; j < data->cluster_count; j++)
        {
            if (clusters[i] == j) continue;
            if (x[i] < LP_EPSILON) continue;

            struct Node *from = &digraph.nodes[i];
            struct Node *to = &digraph.nodes[node_count + j];

//            for (int k = 0; k < graph->node_count; k++)
//            {
//                struct Node *n = &graph->nodes[k];
//                if (clusters[n->index] != clusters[i]) continue;
//
//                int offset = 4 * graph->edge_count + 2 * k;
//                capacities[offset] = 0;
//                capacities[offset + 1] = 0;
//            }

            log_verbose("Calculating max flow from node %d to cluster %to\n", i,
                    j);
            double flow_value;
            rval = flow_find_max_flow(&digraph, capacities, from, to, flow,
                    &flow_value);
            abort_if(rval, "flow_find_max_flow failed");

            log_verbose("    %.2lf\n", flow_value);

            if (flow_value >= 2 * x[i] - LP_EPSILON) continue;

            log_verbose("violation: %.2lf >= %.2lf\n", flow_value, 2 * x[i]);

            int cut_edges_count;
            rval = get_cut_edges_from_marks(&digraph, &cut_edges_count,
                    cut_edges);
            abort_if(rval, "get_cut_edges_from_marks failed");

            log_verbose("Adding cut for i=%d j=%d, cut edges:\n", i, j);
            int c = 0;
            for (int k = 0; k < cut_edges_count / 2; k++)
            {
                int idx = cut_edges[k * 2]->index / 4;
                if (idx > graph->edge_count) continue;

                cut_edges[c++] = &graph->edges[idx];
                log_verbose("    %d %d\n", cut_edges[c - 1]->from->index,
                        cut_edges[c - 1]->to->index);
            }

            rval = GTSP_add_subtour_elimination_cut_2(lp, graph, from, to,
                    cut_edges, c);
            abort_if(rval, "GTSP_add_subtour_elimination_cut failed");

            for (int k = 0; k < graph->node_count; k++)
            {
                int offset = 4 * graph->edge_count + 2 * k;
                capacities[offset] = 1e10;
                capacities[offset + 1] = 1e10;
            }

            (*added_cuts_count)++;
            if (*added_cuts_count > 10) goto CLEANUP;
        }
    }

    // Constraints (2.1)
    for (int i = 0; i < data->cluster_count; i++)
    {
        for (int j = i + 1; j < data->cluster_count; j++)
        {
            struct Node *from = &digraph.nodes[node_count + i];
            struct Node *to = &digraph.nodes[node_count + j];

            log_verbose("Calculating max flow from cluster %d to cluster %to\n",
                    i, j);
            double flow_value;
            rval = flow_find_max_flow(&digraph, capacities, from, to, flow,
                    &flow_value);
            abort_if(rval, "flow_find_max_flow failed");

            log_verbose("    %.2lf\n", flow_value);

            if (flow_value >= 2 - LP_EPSILON) continue;

            log_verbose("violation: %.2lf >= 2\n", flow_value);

            int cut_edges_count;
            rval = get_cut_edges_from_marks(&digraph, &cut_edges_count,
                    cut_edges);
            abort_if(rval, "get_cut_edges_from_marks failed");

            log_verbose("Adding cut for i=%d j=%d, cut edges:\n", i, j);
            int c = 0;
            for (int k = 0; k < cut_edges_count / 2; k++)
            {
                int idx = cut_edges[k * 2]->index / 4;
                if (idx > graph->edge_count) continue;

                cut_edges[c++] = &graph->edges[idx];
                log_verbose("    %d %d\n", cut_edges[c - 1]->from->index,
                        cut_edges[c - 1]->to->index);
            }

            rval = GTSP_add_subtour_elimination_cut_3(lp, graph, from, to,
                    cut_edges, c);
            abort_if(rval, "GTSP_add_subtour_elimination_cut3 failed");

            (*added_cuts_count)++;
            if (*added_cuts_count > 10) goto CLEANUP;
        }
    }

    // subcluster
    struct Node *root = &digraph.nodes[digraph.node_count - 1];
    for (int e_index = 0; e_index < graph->edge_count; e_index++)
    {
        struct Edge *e = &graph->edges[e_index];
        double x_e = x[node_count + e_index];
        if (x_e < LP_EPSILON) continue;

        struct Node *from = &digraph.nodes[e->from->index];
        struct Node *to = &digraph.nodes[e->to->index];

        if (x[from->index] > 1 - LP_EPSILON && x[to->index] > 1 - LP_EPSILON)
            continue;

        capacities[4 * e_index] = capacities[4 * e_index + 2] = 0;

        int cluster_from_index = data->clusters[from->index];
        int cluster_to_index = data->clusters[to->index];

        for (int k = 0; k < data->cluster_count; k++)
        {
            if (cluster_from_index == k) continue;
            if (cluster_to_index == k) continue;

            int offset = 4 * graph->edge_count + 2 * node_count + 2 * k;
            capacities[offset] = 1e10;
            capacities[offset + 1] = 1e10;
        }

        for (int k = 0; k < graph->node_count; k++)
        {
            struct Node *n = &graph->nodes[k];
            if (clusters[n->index] != cluster_from_index &&
                    clusters[n->index] != cluster_to_index)
                continue;

            int offset = 4 * graph->edge_count + 2 * k;
            capacities[offset] = 0;
            capacities[offset + 1] = 0;
        }

        // First direction
        log_debug("Calculating max flow from (%d,%d) to root\n", from->index,
                to->index);
        double flow_value;
        rval = flow_find_max_flow(&digraph, capacities, from, root, flow,
                &flow_value);
        abort_if(rval, "flow_find_max_flow failed");

        log_debug("    %.2lf\n", flow_value);

        if (flow_value + LP_EPSILON < x_e)
        {
            log_debug("violation: %.2lf > %.2lf\n", flow_value, x_e);

            int cut_edges_count;
            rval = get_cut_edges_from_marks(&digraph, &cut_edges_count,
                    cut_edges);
            abort_if(rval, "get_cut_edges_from_marks failed");

            log_debug("Adding cut for i=%d j=root, cut edges:\n", from->index);
            int c = 0;
            for (int k = 0; k < cut_edges_count / 2; k++)
            {
                int idx = cut_edges[k * 2]->index / 4;
                if (idx == e_index) continue;
                if (idx >= graph->edge_count) continue;

                cut_edges[c++] = &graph->edges[idx];
                log_debug("    %d %d\n", cut_edges[c - 1]->from->index,
                        cut_edges[c - 1]->to->index);
            }

            rval = GTSP_add_subcluster_cut(lp, graph, e, cut_edges, c);
            abort_if(rval, "GTSP_add_subcluster_cut failed");

            (*added_cuts_count)++;
            if (*added_cuts_count > 10) goto CLEANUP;

        } else
        {
            // Reverse direction
            log_debug("Trying reverse edge:\n", to->index, from->index);

            rval = flow_find_max_flow(&digraph, capacities, to, root, flow,
                    &flow_value);
            abort_if(rval, "flow_find_max_flow failed");

            log_debug("    %.2lf\n", flow_value);

            if (flow_value + LP_EPSILON < x_e)
            {
                log_debug("violation: %.2lf > %.2lf\n", flow_value, x_e);

                int cut_edges_count;
                rval = get_cut_edges_from_marks(&digraph, &cut_edges_count,
                        cut_edges);
                abort_if(rval, "get_cut_edges_from_marks failed");

                log_debug("Adding cut for i=%d j=root, cut edges:\n",
                        from->index);
                int c = 0;
                for (int k = 0; k < cut_edges_count / 2; k++)
                {
                    int idx = cut_edges[k * 2]->index / 4;
                    if (idx == e_index) continue;
                    if (idx >= graph->edge_count) continue;

                    cut_edges[c++] = &graph->edges[idx];
                    log_debug("    %d %d\n", cut_edges[c - 1]->from->index,
                            cut_edges[c - 1]->to->index);
                }

                rval = GTSP_add_subcluster_cut(lp, graph, e, cut_edges, c);
                abort_if(rval, "GTSP_add_subcluster_cut failed");

                (*added_cuts_count)++;
                if (*added_cuts_count > 10) goto CLEANUP;
            }
        }

        capacities[4 * e_index] = x_e;
        capacities[4 * e_index + 2] = x_e;

        for (int k = 0; k < graph->node_count; k++)
        {
            int offset = 4 * graph->edge_count + 2 * k;
            capacities[offset] = 1e10;
            capacities[offset + 1] = 1e10;
        }

        for (int k = 0; k < data->cluster_count; k++)
        {
            int offset = 4 * graph->edge_count + 2 * graph->node_count;
            capacities[offset + 2 * k] = 0;
            capacities[offset + 2 * k + 1] = 0;
        }
    }

    CLEANUP:
    graph_free(&digraph);
    if (flow) free(flow);
    if (cut_edges) free(cut_edges);
    if (capacities) free(capacities);
    if (x) free(x);
    return rval;
}

int GTSP_add_cutting_planes(struct LP *lp, struct GTSP *data)
{
    int rval = 0;
    double *x = 0;

    int num_cols = LP_get_num_cols(lp);
    x = (double *) malloc(num_cols * sizeof(double));
    abort_if(!x, "could not allocate x");

    while (1)
    {
        int added_cuts_count = 0;

        log_debug("Finding subtour cuts...\n");

        rval = GTSP_find_exact_subtour_elimination_cuts(lp, data,
                &added_cuts_count);
        abort_if(rval, "GTSP_find_exact_subtour_elimination_cuts failed");

        if (added_cuts_count > 0)
        {
            log_debug("Found %d subtour elimination cuts using exact "
                    "separation\n", added_cuts_count);
        } else break;

        log_debug("Reoptimizing...\n");
        int is_infeasible;
        rval = LP_optimize(lp, &is_infeasible);
        abort_if(rval, "LP_optimize failed");

        if (is_infeasible) break;

        double objval;
        rval = LP_get_obj_val(lp, &objval);
        abort_if(rval, "LP_get_obj_val failed");

        rval = LP_get_x(lp, x);
        abort_if(rval, "LP_get_x failed");

        log_debug("Current solution:\n");
        for (int i = 0; i < data->graph->node_count; i++)
            if (x[i] > LP_EPSILON) log_debug("    node %d = %.2f\n", i, x[i]);

        for (int i = 0; i < data->graph->edge_count; i++)
        {
            struct Edge *e = &data->graph->edges[i];
            int idx = e->index + data->graph->node_count;
            if (x[idx] > LP_EPSILON)
            {
                log_debug("    edge (%d, %d) = %.2f\n", e->from->index,
                        e->to->index, x[idx]);
            }
        }

        log_debug("    obj val = %f\n", objval);
    }

    CLEANUP:
    if (x) free(x);
    return rval;
}

int GTSP_write_input_data(struct GTSP *data, char *filename)
{
    int rval = 0;

    FILE *file;

    file = fopen(filename, "w");
    abort_if(!file, "could not open file");

    fprintf(file, "%d %d\n", data->graph->node_count, data->cluster_count);

    for (int i = 0; i < data->graph->node_count; i++)
    {
        fprintf(file, "%.2lf %.2lf %d\n", data->x_coordinates[i],
                data->y_coordinates[i], data->clusters[i]);
    }

    CLEANUP:
    if (file) fclose(file);
    return rval;
}

int GTSP_write_solution(struct GTSP *data, char *filename, double *x)
{
    int rval = 0;

    struct Edge *edges = data->graph->edges;
    int node_count = data->graph->node_count;
    int edge_count = data->graph->edge_count;

    FILE *file;
    file = fopen(filename, "w");
    abort_if(!file, "could not open file");

    int positive_edge_count = 0;
    for (int i = 0; i < edge_count; i++)
        if (x[i + node_count] > LP_EPSILON)
            positive_edge_count++;

    fprintf(file, "%d %d\n", node_count, edge_count);

    fprintf(file, "%d\n", positive_edge_count);

    for (int i = 0; i < edge_count; i++)
        if (x[i + node_count] > LP_EPSILON)
            fprintf(file, "%d %d %.4lf\n", edges[i].from->index,
                    edges[i].to->index, x[i + node_count]);

    CLEANUP:
    if (file) fclose(file);
    return rval;
}

int get_edge_num(int node_count, int from, int to)
{
    int idx = node_count;

    for (int k = 0; k < from; k++)
        idx += node_count - k - 1;

    idx += to - from - 1;

    return idx;
}

int GTSP_read_x(char *filename, double **p_x)
{
    int rval = 0;

    int node_count;
    int edge_count;

    double *x;

    FILE *file;

    log_info("Reading optimal solution from file %s\n", filename);

    file = fopen(filename, "r");
    abort_if(!file, "could not open file");

    rval = fscanf(file, "%d %d", &node_count, &edge_count);
    abort_if(rval != 2, "invalid input format (node and edge count)");

    int num_cols = node_count + edge_count;

    x = (double *) malloc(num_cols * sizeof(double));
    abort_if(!x, "could not allocate x");

    for (int i = 0; i < node_count + edge_count; i++) x[i] = 0.0;

    rval = fscanf(file, "%d", &edge_count);
    abort_if(rval != 1, "invalid input format (positive edge count)");

    for (int i = 0; i < edge_count; i++)
    {
        int from, to, edge;
        rval = fscanf(file, "%d %d", &from, &to);
        abort_if(rval != 2, "invalid input format (edge endpoints)");

        if (from > to) swap(from, to);

        edge = get_edge_num(node_count, from, to);
        abort_if(edge > num_cols, "invalid edge");

        x[from] += 0.5;
        x[to] += 0.5;
        x[edge] = 1;
    }

    for (int i = 0; i < num_cols; i++)
    {
        if (x[i] <= LP_EPSILON) continue;
        log_debug(" x%-3d = %.2f\n", i, x[i]);
    }

    *p_x = x;
    rval = 0;

    CLEANUP:
    return rval;
}

static const struct option options_tab[] =
        {{"help", no_argument, 0, 'h'}, {"nodes", required_argument, 0, 'n'},
                {"clusters", required_argument, 0, 'm'},
                {"grid-size", required_argument, 0, 'g'},
                {"optimal", required_argument, 0, 'x'},
                {"seed", required_argument, 0, 's'},
                {(char *) 0, (int) 0, (int *) 0, (int) 0}};

static int input_node_count = 20;
static int input_cluster_count = 5;
static int grid_size = 100;

static void GTSP_print_usage(char **argv)
{
    printf("wrong usage\n");
}

static int GTSP_parse_args(int argc, char **argv)
{
    int rval = 0;

    opterr = 0;

    while (1)
    {
        int c = 0;
        int option_index = 0;
        c = getopt_long(argc, argv, "n:m:g:x:s:", options_tab, &option_index);

        if (c < 0) break;

        switch (c)
        {
            case 'n':
                input_node_count = atoi(optarg);
                break;

            case 'm':
                input_cluster_count = atoi(optarg);
                break;

            case 'g':
                grid_size = atoi(optarg);
                break;

            case 'x':
                rval = GTSP_read_x(optarg, &OPTIMAL_X);
                abort_if(rval, "GTSP_read_x failed");
                break;

            case 's':
                SEED = (unsigned) atoi(optarg);
                break;

            case ':':
                fprintf(stderr, "option '-%c' requires an argument\n", optopt);
                return 1;

            case '?':
            default:
                fprintf(stderr, "option '-%c' is invalid\n", optopt);
                return 1;

        }
    }

    CLEANUP:
    return rval;
}

int GTSP_solution_found(struct GTSP *data, double *x)
{
    int rval = 0;

    log_info("Writting solution to file gtsp.out\n");
    rval = GTSP_write_solution(data, "gtsp.out", x);
    abort_if(rval, "GTSP_write_solution failed");

    CLEANUP:
    return rval;
}

int GTSP_main(int argc, char **argv)
{
    int rval = 0;

    struct BNC bnc;
    struct GTSP data;

    SEED = (unsigned int) get_real_time() % 1000000;

    rval = GTSP_init_data(&data);
    abort_if(rval, "GTSP_init_data failed");

    rval = BNC_init(&bnc);
    abort_if(rval, "BNC_init failed");

    rval = GTSP_parse_args(argc, argv);
    if (rval) return 1;

    srand(SEED);

    log_info("Generating random GTSP instance...\n");
    log_info("    seed = %d\n", SEED);
    log_info("    input_node_count = %d\n", input_node_count);
    log_info("    input_cluster_count = %d\n", input_cluster_count);
    log_info("    grid_size = %d\n", grid_size);

    rval = GTSP_create_random_problem(input_node_count, input_cluster_count,
            grid_size, &data);
    abort_if(rval, "GTSP_create_random_problem failed");

    log_info("Writing random instance to file gtsp.in\n");
    rval = GTSP_write_input_data(&data, "gtsp.in");
    abort_if(rval, "GTSP_write_problem failed");

    bnc.best_obj_val = DBL_MAX;
    bnc.problem_data = (void *) &data;
    bnc.problem_init_lp = (int (*)(struct LP *, void *)) GTSP_init_lp;
    bnc.problem_add_cutting_planes =
            (int (*)(struct LP *, void *)) GTSP_add_cutting_planes;
    bnc.problem_solution_found =
            (int (*)(void *, double *)) GTSP_solution_found;

    if (OPTIMAL_X)
    {
        log_info("Optimal solution is available. Cuts will be checked.\n");

        double opt_val = 0.0;
        for (int i = 0; i < data.graph->edge_count; i++)
        {
            struct Edge *e = &data.graph->edges[i];
            opt_val += OPTIMAL_X[i + input_node_count] * e->weight;
        }

        log_info("    opt = %.2lf\n", opt_val);
    }

    log_info("Initializing LP...\n");
    rval = BNC_init_lp(&bnc);
    abort_if(rval, "BNC_init_lp failed");

    log_info("Writing LP to file gtsp.lp...\n");
    rval = LP_write(bnc.lp, "gtsp.lp");
    abort_if(rval, "LP_write failed");

    log_info("Starting branch-and-cut solver...\n");
    rval = BNC_solve(&bnc);
    abort_if(rval, "BNC_solve_node failed");

    abort_if(!bnc.best_x, "problem has no feasible solution");

    log_info("Optimal integral solution:\n");
    log_info("    obj value = %.2lf **\n", bnc.best_obj_val);

    log_info("Branch-and-bound nodes: %d\n", BNC_NODE_COUNT);
    log_info("Max-flow computations: %d\n", FLOW_MAX_FLOW_COUNT);

    CLEANUP:
    GTSP_free(&data);
    BNC_free(&bnc);
    return rval;
}