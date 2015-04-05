/* Copyright (c) 2015 Armin Sadeghi
 * Copyright (c) 2015 Alinson Xavier
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include "gtsp.h"
#include "geometry.h"
#include "util.h"
#include "gtsp-subtour.h"
#include "gtsp-comb.h"
#include "gtsp-cols.h"

int GTSP_init_data(struct GTSP *data)
{
    int rval = 0;

    data->node_to_cluster = 0;
    data->cluster_count = 0;

    data->graph = (struct Graph *) malloc(sizeof(struct Graph));
    abort_if(!data->graph, "could not allocate data->graph");

    data->clusters = (struct Cluster *) malloc(sizeof(struct Cluster));
    abort_if(!data->clusters, "could not allocate data->clusters");

    graph_init(data->graph);

    CLEANUP:
    return rval;
}

void GTSP_free(struct GTSP *data)
{
    if (!data) return;

    for (int i = 0; i < data->graph->node_count; i++)
        free(data->dist_matrix[i]);

    for (int i = 0; i < data->cluster_count; i++)
        free(data->clusters[i].nodes);

    if (data->clusters) free(data->clusters);
    if (data->dist_matrix) free(data->dist_matrix);
    if (data->node_to_cluster) free(data->node_to_cluster);

    graph_free(data->graph);
    free(data->graph);
}

int GTSP_create_random_problem(
        int node_count, int cluster_count, int grid_size, struct GTSP *data)
{
    int rval = 0;
    int *edges = 0;
    int *weights = 0;
    int *node_to_cluster = 0;
    struct Cluster *clusters = 0;
    int *edge_map = 0;

    int **dist_matrix = 0;

    double *x_coords = 0;
    double *y_coords = 0;

    struct Graph *graph = 0;

    int edge_count = (node_count * (node_count - 1)) / 2;
    graph = (struct Graph *) malloc(sizeof(struct Graph));
    abort_if(!graph, "could not allocate graph\n");

    graph_init(graph);

    edges = (int *) malloc(2 * edge_count * sizeof(int));
    weights = (int *) malloc(edge_count * sizeof(int));
    node_to_cluster = (int *) malloc(node_count * sizeof(int));
    abort_if(!data->graph, "could not allocate data->graph");
    abort_if(!edges, "could not allocate data->edges\n");
    abort_if(!weights, "could not allocate weights\n");
    abort_if(!node_to_cluster, "could not allocate node_to_cluster\n");

    x_coords = (double *) malloc(node_count * sizeof(double));
    y_coords = (double *) malloc(node_count * sizeof(double));

    abort_if(!x_coords, "could not allocate x_coords\n");
    abort_if(!y_coords, "could not allocate y_coords\n");

    dist_matrix = (int **) malloc(node_count * sizeof(int *));
    abort_if(!dist_matrix, "could not allocate dist_matrix\n");

    for (int i = 0; i < node_count; i++)
    {
        dist_matrix[i] = (int *) malloc(node_count * sizeof(int));
        abort_iff(!dist_matrix[i], "could not allocate dist_matrix[%d]\n", i);
    }

    rval = generate_random_clusters_2d(node_count, cluster_count, grid_size,
            x_coords, y_coords, node_to_cluster);
    abort_if(rval, "generate_random_clusters_2d failed");

    rval = generate_dist_matrix(node_count, x_coords, y_coords, dist_matrix);
    abort_if(rval, "generate_distance_matrix_2d failed");

    clusters = (struct Cluster *) malloc(
            cluster_count * sizeof(struct Cluster));
    abort_if(!clusters, "could not allocate clusters");

    for (int i = 0; i < cluster_count; i++)
        clusters[i].size = 0;

    for (int i = 0; i < node_count; i++)
        clusters[node_to_cluster[i]].size += 1;

    for (int i = 0; i < cluster_count; i++)
    {
        clusters[i].nodes = (int *) malloc(clusters[i].size * sizeof(int));
        abort_iff(!clusters[i].nodes, "could not allocate clusters[%d].nodes",
                i);
    }

    int current_vertex = 0;
    for (int j = 0; j < cluster_count; j++)
    {
        current_vertex = 0;
        for (int i = 0; i < node_count; i++)
            if (node_to_cluster[i] == j)
            {
                clusters[j].nodes[current_vertex] = i;
                current_vertex += 1;
            }
    }

    int curr_edge = 0;
    for (int i = 0; i < edge_count; i++)
        for (int j = i + 1; j < node_count; j++)
        {
            if (node_to_cluster[i] == node_to_cluster[j]) continue;

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
    data->node_to_cluster = node_to_cluster;
    data->cluster_count = cluster_count;
    graph->x_coordinates = x_coords;
    graph->y_coordinates = y_coords;
    data->dist_matrix = dist_matrix;
    data->clusters = clusters;

    CLEANUP:
    if (edge_map) free(edge_map);
    if (weights) free(weights);
    if (edges) free(edges);
    if (rval)
    {
        if (clusters) free(clusters);
        if (node_to_cluster) free(node_to_cluster);
        if (dist_matrix) free(dist_matrix);
    }
    return rval;
}

int GTSP_init_lp(struct LP *lp, struct GTSP *data)
{
    int rval = 0;

    int edge_count = data->graph->edge_count;
    struct Edge *edges = data->graph->edges;

    double lb = 0.0;
    double ub = 1.0;
    int empty = 0;
    double emptyd = 0.0;

    int k = 0;
    for (int i = 0; i < edge_count; i++)
    {
        if (edges[i].column < 0) continue;
        edges[i].column = k++;

        double obj = (double) edges[i].weight;

        rval = LP_add_cols(lp, 1, 0, &obj, &empty, &empty, &emptyd, &lb, &ub);
        abort_if(rval, "LP_add_cols failed");
    }

    CLEANUP:
    return rval;
}

int GTSP_add_cutting_planes(struct LP *lp, struct GTSP *data)
{
    int rval = 0;
    int current_round = 0;

    while (1)
    {
        if (current_round > 0)
        {
            int is_infeasible;
            rval = LP_optimize(lp, &is_infeasible);
            abort_if(rval, "LP_optimize failed");

            if (is_infeasible) break;
        }

        current_round++;
        abort_if(get_user_time() - INITIAL_TIME >= MAX_TOTAL_TIME,
                "time limit exceeded");

        int original_cut_pool_size;
        int added_cuts_count;

        // Generalized subtour cuts
        original_cut_pool_size = lp->cut_pool_size;
        log_debug("Finding subtour cuts, round %d...\n", current_round);

        rval = GTSP_find_exact_subtour_cuts(lp, data);
        abort_if(rval, "GTSP_find_exact_subtour_cuts failed");

        added_cuts_count = lp->cut_pool_size - original_cut_pool_size;
        if (added_cuts_count > 0)
            continue;

        // Generalized comb cuts
        original_cut_pool_size = lp->cut_pool_size;
        log_debug("Finding comb cuts, round %d...\n", current_round);

        rval = GTSP_find_comb_cuts(lp, data);
        abort_if(rval, "GTSP_find_comb_cuts failed");

        added_cuts_count = lp->cut_pool_size - original_cut_pool_size;
        if (added_cuts_count > 0)
            continue;

        // Column generation
        int original_cols_count = LP_get_num_cols(lp);

        rval = GTSP_find_columns(lp, data);
        abort_if(rval, "GTSP_find_columns failed");

        int added_cols_count = LP_get_num_cols(lp) - original_cols_count;

        if (added_cols_count > 0)
            continue;

        break;
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

    const struct Graph *graph = data->graph;

    fprintf(file, "%d %d %d\n", graph->node_count, data->cluster_count,
            graph->edge_count);

    for (int i = 0; i < graph->node_count; i++)
    {
        fprintf(file, "%.2lf %.2lf %d\n", graph->x_coordinates[i],
                graph->y_coordinates[i], data->node_to_cluster[i]);
    }

    CLEANUP:
    if (file) fclose(file);
    return rval;
}

int GTSP_print_solution(struct GTSP *data, double *x)
{
    struct Edge *edges = data->graph->edges;
    int edge_count = data->graph->edge_count;

    for (int i = 0; i < edge_count; i++)
    {
        int col = edges[i].column;
        if (col < 0) continue;
        if (x[col] > EPSILON)
            log_info("    %-3d %-3d %8.4lf\n", edges[i].from->index,
                    edges[i].to->index, x[col]);
    }

    return 0;
}

int GTSP_write_solution(struct GTSP *data, char *filename, double *x)
{
    int rval = 0;

    struct Edge *edges = data->graph->edges;
    int edge_count = data->graph->edge_count;

    FILE *file;

    file = fopen(filename, "w");
    abort_if(!file, "could not open file");

    int positive_edge_count = 0;
    for (int i = 0; i < edge_count; i++)
    {
        struct Edge *e = &edges[i];

        if (e->column < 0) continue;
        if (x[e->column] < EPSILON) continue;

        positive_edge_count++;
    }

    fprintf(file, "%d\n", positive_edge_count);

    for (int i = 0; i < edge_count; i++)
    {
        struct Edge *e = &edges[i];

        if (e->column < 0) continue;
        if (x[e->column] < EPSILON) continue;

        fprintf(file, "%d %d %.4lf\n", edges[i].from->index, edges[i].to->index,
                x[e->column]);
    }

    CLEANUP:
    if (file) fclose(file);
    return rval;
}

int GTSP_read_solution(struct GTSP *gtsp, char *filename, double **p_x)
{
    int rval = 0;

    int node_count;
    int edge_count;
    int *edge_map = 0;

    double *x;

    FILE *file;

    log_info("Reading optimal solution from file %s\n", filename);

    file = fopen(filename, "r");
    abort_if(!file, "could not open file");

    rval = fscanf(file, "%d %d", &node_count, &edge_count);
    abort_if(rval != 2, "invalid input format (node and edge count)");

    x = (double *) malloc(edge_count * sizeof(double));
    abort_if(!x, "could not allocate x");

    for (int i = 0; i < edge_count; i++) x[i] = 0.0;

    rval = fscanf(file, "%d", &edge_count);
    abort_if(rval != 1, "invalid input format (positive edge count)");

    edge_map = (int *) malloc(node_count * node_count * sizeof(int));
    abort_if(!edge_map, "could not allocate edge_map");

    rval = GTSP_build_edge_map(gtsp, edge_map);
    abort_if(rval, "GTSP_build_edge_map failed");

    for (int i = 0; i < edge_count; i++)
    {
        int from, to, edge;
        double val;
        rval = fscanf(file, "%d %d %lf", &from, &to, &val);
        abort_if(rval != 3, "invalid input format (edge endpoints)");

        edge = edge_map[from * node_count + to];

        x[edge] = val;
    }

    for (int i = 0; i < edge_count; i++)
    {
        if (x[i] <= EPSILON) continue;
        log_debug("    x%-5d = %.6f\n", i, x[i]);
    }

    *p_x = x;
    rval = 0;

    CLEANUP:
    if (file) fclose(file);
    if (edge_map) free(edge_map);
    return rval;
}

int GTSP_check_solution(struct GTSP *data, double *x)
{
    int rval = 0;
    int *cluster_mark = 0;

    struct Node **stack = 0;
    int stack_top = 0;

    struct Graph *graph = data->graph;
    const int node_count = graph->node_count;
    const int edge_count = graph->edge_count;

    cluster_mark = (int *) malloc(data->cluster_count * sizeof(int));
    abort_if(!cluster_mark, "could not allocate cluster_mark");

    stack = (struct Node **) malloc(graph->node_count * sizeof(struct Node *));
    abort_if(!stack, "could not allocate stack");

    for (int i = 0; i < edge_count; i++)
    {
        int col = graph->edges[i].column;
        if (col < 0) continue;

        abort_iff(x[col] < 1.0 - EPSILON && x[col] > EPSILON,
                "solution is not integral: x%d = %.4lf", i, x[col]);

        abort_iff(x[col] > 1.0 + EPSILON || x[col] < 0.0 - EPSILON,
                "value out of bounds: x%d = %.4lf", i, x[col]);
    }

    for (int i = 0; i < node_count; i++)
        graph->nodes[i].mark = 0;

    for (int i = 0; i < data->cluster_count; i++)
        cluster_mark[i] = 0;

    int initial;
    for (initial = 0; initial < edge_count; initial++)
        if (x[initial] > 1.0 - EPSILON) break;

    abort_if(initial == edge_count, "no initial node");

    for (int i = 0; i < edge_count; i++)
    {
        if (graph->edges[i].column == initial)
        {
            stack[stack_top++] = graph->edges[i].from;
            graph->edges[i].from->mark = 1;
            break;
        }
    }

    while (stack_top > 0)
    {
        struct Node *n = stack[--stack_top];
        cluster_mark[data->node_to_cluster[n->index]]++;

        for (int i = 0; i < n->degree; i++)
        {
            struct Adjacency *adj = &n->adj[i];
            struct Node *neighbor = adj->neighbor;

            if (neighbor->mark) continue;

            if (adj->edge->column < 0) continue;
            if (x[adj->edge->column] < EPSILON) continue;

            stack[stack_top++] = neighbor;
            neighbor->mark = 1;
        }
    }

    for (int i = 0; i < data->cluster_count; i++)
        abort_iff(cluster_mark[i] < 1, "cluster %d not visited", i);

    log_info("    solution is valid\n");

    CLEANUP:
    if (stack) free(stack);
    if (cluster_mark) free(cluster_mark);
    return rval;
}

int GTSP_solution_found(struct BNC *bnc, struct GTSP *data, double *x)
{
    int rval = 0;

    if (strlen(SOLUTION_FILENAME) > 0)
    {
        log_info("Writing solution to file %s\n", SOLUTION_FILENAME);
        rval = GTSP_write_solution(data, SOLUTION_FILENAME, x);
        abort_if(rval, "GTSP_write_solution failed");
    }

    log_info("Checking solution...\n");
    rval = GTSP_check_solution(data, x);
    abort_if(rval, "GTSP_check_solution failed");

    CLEANUP:
    return rval;
}

int GTSP_build_edge_map(struct GTSP *gtsp, int *edge_map)
{
    int node_count = gtsp->graph->node_count;

    int k = 0;
    for (int i = 0; i < node_count; i++)
    {
        for (int j = i + 1; j < node_count; j++)
        {
            if (gtsp->node_to_cluster[i] == gtsp->node_to_cluster[j]) continue;
            edge_map[i * node_count + j] = k;
            edge_map[j * node_count + i] = k;
            k++;
        }
    }

    return 0;
}
