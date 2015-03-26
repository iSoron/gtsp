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
#include "math.h"

static double *OPTIMAL_X = 0;

static int get_edge_num(int node_count, int from, int to)
{
    int idx = node_count;

    for (int k = 0; k < from; k++)
        idx += node_count - k - 1;

    idx += to - from - 1;

    return idx;
}

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
            if (clusters[i] == clusters[j]) continue;

            edges[curr_edge * 2] = i;
            edges[curr_edge * 2 + 1] = j;
            weights[curr_edge] = get_euclidean_distance(x_coords, y_coords, i,
                    j);

            curr_edge++;
        }

    edge_count = curr_edge;

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

static int add_subtour_cut(
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

    rval = LP_add_rows(lp, 1, newnz, &rhs, &sense, &rmatbeg, rmatind, rmatval);
    abort_if(rval, "LP_add_rows failed");

    CLEANUP:
    if (rmatval) free(rmatval);
    if (rmatind) free(rmatind);
    return rval;
}

static int build_flow_digraph(
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
//        if (x[node_count + i] < LP_EPSILON) continue;

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

//    digraph_edge_count = kc;

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

int find_exact_subtour_cuts_cluster_to_cluster(
        struct LP *lp,
        struct GTSP *data,
        struct Graph *digraph,
        double *capacities,
        int *added_cuts_count)
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

            if (flow_value >= 2 - LP_EPSILON) continue;

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

            (*added_cuts_count)++;
            cuts_count++;
        }
    }

    CLEANUP:
    if (cut_edges) free(cut_edges);
    if (flow) free(flow);
    return rval;
}

int find_exact_subtour_cuts_node_to_cluster(
        struct LP *lp,
        struct GTSP *data,
        double *x,
        struct Graph *digraph,
        double *capacities,
        int *added_cuts_count)
{
    int rval = 0;

    int cuts_count = 0;
    struct Edge **cut_edges = 0;
    double *flow = 0;

    struct Graph *graph = data->graph;
    int *clusters = data->clusters;

    cut_edges = (struct Edge **) malloc(
            digraph->edge_count * sizeof(struct Edge *));
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

            log_verbose("Sending flow from node %d to cluster %d (must be >= %.4lf)\n", i, j, 2*x[i]);

            activate_cluster_node(capacities, to);

            double flow_value;
            int cut_edges_count;

            rval = flow_find_max_flow(digraph, capacities, from, to, flow,
                    &flow_value);
            abort_if(rval, "flow_find_max_flow failed");

            log_verbose("    flow value = %.4lf\n", flow_value);

            deactivate_cluster_node(capacities, to);

            if (flow_value + LP_EPSILON >= 2 * x[i])
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
                log_verbose("  %d (%d) %d (%d) [%d]\n", e->from->index, e->from->mark,
                        e->to->index, e->to->mark, e->index);
            }

            rval = add_subtour_cut(lp, graph, from, 0, cut_edges,
                    cut_edges_count, 1);
            abort_if(rval, "add_subtour_cut failed");

            (*added_cuts_count)++;
            cuts_count++;
        }
    }

    CLEANUP:
    if (cut_edges) free(cut_edges);
    if (flow) free(flow);
    return rval;
}

int find_exact_subtour_cuts_node_to_node(
        struct LP *lp,
        struct GTSP *data,
        double *x,
        struct Graph *digraph,
        double *capacities,
        int *added_cuts_count)
{
    int rval = 0;

    struct Edge **cut_edges = 0;
    double *flow = 0;

    struct Graph *graph = data->graph;
    int *clusters = data->clusters;

    cut_edges = (struct Edge **) malloc(
            digraph->edge_count * sizeof(struct Edge *));
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

        if (flow_value >= 2 * (x[i] + x[j] - 1) - LP_EPSILON)
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

        (*added_cuts_count)++;
    }

    CLEANUP:
    if (flow) free(flow);
    if (cut_edges) free(cut_edges);
    return rval;
}

int find_exact_subtour_cuts(
        struct LP *lp, struct GTSP *data, int *total_added_cuts)
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

    capacities = (double *) malloc(digraph_edge_count * sizeof(double));
    abort_if(!capacities, "could not allocate capacities");

    rval = build_flow_digraph(data, x, &digraph, capacities);
    abort_if(rval, "build_flow_digraph failed");

    // Constraints (2.1)
    rval = find_exact_subtour_cuts_cluster_to_cluster(lp, data, &digraph,
            capacities, &added_cuts_count);
    abort_if(rval, "find_exact_subtour_cuts_cluster_to_cluster failed");

    log_debug("Added %d cluster-to-cluster subtour cuts\n", added_cuts_count);
    (*total_added_cuts) += added_cuts_count;

    if (added_cuts_count > 0) goto CLEANUP;

    // Constraints (2.2)
    rval = find_exact_subtour_cuts_node_to_cluster(lp, data, x, &digraph,
            capacities, &added_cuts_count);
    abort_if(rval, "find_exact_subtour_cuts_node_to_cluster failed");

    log_debug("Added %d node-to-cluster subtour cuts\n", added_cuts_count);
    (*total_added_cuts) += added_cuts_count;

    // Constraints (2.3)
    rval = find_exact_subtour_cuts_node_to_node(lp, data, x, &digraph,
            capacities, &added_cuts_count);
    abort_if(rval, "find_exact_subtour_cuts_node_to_node failed");

    log_debug("Added %d node-to-node subtour cuts\n", added_cuts_count);
    (*total_added_cuts) += added_cuts_count;

    CLEANUP:
    graph_free(&digraph);
    if (capacities) free(capacities);
    if (x) free(x);
    return rval;
}

int GTSP_add_cutting_planes(struct LP *lp, struct GTSP *data)
{
    int rval = 0;

    int round = 0;
    int added_cuts_count = 0;

    while (1)
    {
        round++;
        int added_cuts_this_round = 0;

        log_debug("Finding subtour cuts, round %d...\n", round);

        rval = find_exact_subtour_cuts(lp, data, &added_cuts_this_round);
        abort_if(rval, "find_exact_subtour_cuts failed");

        if (added_cuts_this_round == 0)
        {
            log_debug("No more subtour cuts found.\n");
            break;
        }

        log_debug("Reoptimizing...\n");

        int is_infeasible;
        rval = LP_optimize(lp, &is_infeasible);
        abort_if(rval, "LP_optimize failed");

        double obj_val;
        rval = LP_get_obj_val(lp, &obj_val);
        abort_if(rval, "LP_get_obj_val failed");

        log_debug("    obj val = %.4lf\n", obj_val);

        if (is_infeasible) break;

        added_cuts_count += added_cuts_this_round;
    }

    CLEANUP:
    return rval;
}

int GTSP_write_problem(struct GTSP *data, char *filename)
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

int GTSP_read_solution(char *filename, double **p_x)
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

static const struct option options_tab[] = {{"help", no_argument, 0, 'h'},
        {"nodes", required_argument, 0, 'n'},
        {"clusters", required_argument, 0, 'm'},
        {"grid-size", required_argument, 0, 'g'},
        {"optimal", required_argument, 0, 'x'},
        {"seed", required_argument, 0, 's'},
        {(char *) 0, (int) 0, (int *) 0, (int) 0}};

static int input_node_count = -1;
static int input_cluster_count = -1;
static int grid_size = 100;

static void GTSP_print_usage()
{
    printf("Parameters:\n");
    printf("%4s %-13s %s\n", "-n", "--nodes", "number of nodes");
    printf("%4s %-13s %s\n", "-m", "--clusters", "number of clusters");
    printf("%4s %-13s %s\n", "-s", "--seed", "random seed");
    printf("%4s %-13s %s\n", "-g", "--grid-size",
            "size of the box used for generating random points");
    printf("%4s %-13s %s\n", "-x", "--optimal",
            "file containg valid solution (used to assert validity of cuts)");
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
                rval = GTSP_read_solution(optarg, &OPTIMAL_X);
                abort_if(rval, "GTSP_read_solution failed");
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

    if (input_cluster_count < 0)
    {
        input_cluster_count = (int) ceil(input_node_count / 5.0);
        if (input_cluster_count < 3) input_cluster_count = 3;
    }

    if (input_node_count < 0)
    {
        printf("You must specify the number of nodes.\n");
        rval = 1;
    }

    if (input_cluster_count > input_node_count)
    {
        printf("Number of clusters must be at most number of nodes.\n");
        rval = 1;
    }

    if (rval)
    {
        GTSP_print_usage();
        rval = 1;
    }

    CLEANUP:
    return rval;
}

int GTSP_solution_found(struct GTSP *data, double *x)
{
    int rval = 0;

    log_info("Writting integral solution to file gtsp.out\n");
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
    rval = GTSP_write_problem(&data, "gtsp.in");
    abort_if(rval, "GTSP_write_problem failed");

    bnc.best_obj_val = DBL_MAX;
    bnc.problem_data = (void *) &data;
    bnc.problem_init_lp = (int (*)(struct LP *, void *)) GTSP_init_lp;
    bnc.problem_add_cutting_planes = (int (*)(
            struct LP *, void *)) GTSP_add_cutting_planes;
    bnc.problem_solution_found = (int (*)(
            void *, double *)) GTSP_solution_found;

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
