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
#include <assert.h>
#include "gtsp.h"
#include "geometry.h"
#include "util.h"
#include "gtsp-subtour.h"
#include "gtsp-comb.h"

int large_neighborhood_search(
        struct Tour *tour, struct GTSP *data, int *tour_cost);

int build_edge_map(struct GTSP *gtsp, int *edge_map);

double *OPTIMAL_X = 0;

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
    int cluster_count = data->cluster_count;
    int *clusters = data->node_to_cluster;
    struct Edge *edges = data->graph->edges;

    for (int i = 0; i < cluster_count; i++)
    {
        rval = LP_new_row(lp, 'E', 2.0);
        abort_if(rval, "LP_new_row failed");
    }

    double lb = 0.0;
    double ub = 1.0;
    int cmatbeg = 0;

    for (int i = 0; i < edge_count; i++)
    {
        struct Node *from = edges[i].from;
        struct Node *to = edges[i].to;

        double obj = (double) edges[i].weight;
        double cmatval[] = {1.0, 1.0};
        int cmatind[] = {clusters[from->index], clusters[to->index]};

        rval = LP_add_cols(lp, 1, 2, &obj, &cmatbeg, cmatind, cmatval, &lb,
                &ub);
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

        original_cut_pool_size = lp->cut_pool_size;
        log_debug("Finding subtour cuts, round %d...\n", current_round);

        rval = GTSP_find_exact_subtour_cuts(lp, data);
        abort_if(rval, "GTSP_find_exact_subtour_cuts failed");

        added_cuts_count = lp->cut_pool_size - original_cut_pool_size;
        if (added_cuts_count > 0)
            continue;

        original_cut_pool_size = lp->cut_pool_size;
        log_debug("Finding comb cuts, round %d...\n", current_round);

        rval = GTSP_find_comb_cuts(lp, data);
        abort_if(rval, "GTSP_find_comb_cuts failed");

        added_cuts_count = lp->cut_pool_size - original_cut_pool_size;
        if (added_cuts_count > 0)
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
        if (x[i] > EPSILON)
            log_info("    %-3d %-3d %8.4lf\n", edges[i].from->index,
                    edges[i].to->index, x[i]);

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
        if (x[i] > EPSILON)
            positive_edge_count++;

    fprintf(file, "%d\n", positive_edge_count);

    for (int i = 0; i < edge_count; i++)
        if (x[i] > EPSILON)
            fprintf(file, "%d %d %.4lf\n", edges[i].from->index,
                    edges[i].to->index, x[i]);

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

    rval = build_edge_map(gtsp, edge_map);
    abort_if(rval, "build_edge_map failed");

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

int build_edge_map(struct GTSP *gtsp, int *edge_map)
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
        abort_iff(x[i] < 1.0 - EPSILON && x[i] > EPSILON,
                "solution is not integral: x%d = %.4lf", i, x[i]);

        abort_iff(x[i] > 1.0 + EPSILON || x[i] < 0.0 - EPSILON,
                "value out of bounds: x%d = %.4lf", i, x[i]);
    }

    for (int i = 0; i < node_count; i++)
        graph->nodes[i].mark = 0;

    for (int i = 0; i < data->cluster_count; i++)
        cluster_mark[i] = 0;

    int initial;
    for (initial = 0; initial < edge_count; initial++)
        if (x[initial] > 1.0 - EPSILON) break;

    abort_if(initial == edge_count, "no initial node");

    stack[stack_top++] = graph->edges[initial].from;
    graph->edges[initial].from->mark = 1;

    while (stack_top > 0)
    {
        struct Node *n = stack[--stack_top];
        cluster_mark[data->node_to_cluster[n->index]]++;

        for (int i = 0; i < n->degree; i++)
        {
            struct Adjacency *adj = &n->adj[i];
            struct Node *neighbor = adj->neighbor;

            if (neighbor->mark) continue;
            if (x[adj->edge->index] < EPSILON) continue;

            stack[stack_top++] = neighbor;
            neighbor->mark = 1;
        }
    }

    for (int i = 0; i < data->cluster_count; i++)
    {
        abort_iff(cluster_mark[i] > 1, "cluster %d visited multiple times", i);
        abort_iff(cluster_mark[i] < 1, "cluster %d not visited", i);
    }

    log_info("    solution is valid\n");

    CLEANUP:
    if (stack) free(stack);
    if (cluster_mark) free(cluster_mark);
    return rval;
}

int GTSP_solution_found(struct BNC *bnc, struct GTSP *data, double *x)
{
    UNUSED(bnc);

    int rval = 0;
    int tour_cost;
    double *best_val = &bnc->best_obj_val;

    struct Tour *tour;
    tour = (struct Tour *) malloc(data->cluster_count * sizeof(struct Tour));

    if (strlen(SOLUTION_FILENAME) > 0)
    {
        log_info("Writing solution to file %s\n", SOLUTION_FILENAME);
        rval = GTSP_write_solution(data, SOLUTION_FILENAME, x);
        abort_if(rval, "GTSP_write_solution failed");
    }

    rval = build_tour_from_x(data, tour, x);
    abort_if(rval, "build_tour_from_x failed");

    rval = large_neighborhood_search(tour, data, &tour_cost);
    abort_if(rval, "large_neighborhood_search failed");

    if (tour_cost + EPSILON < *best_val)
    {
        log_info("Local search improve the integral solution\n");
        log_info("         obj val = %f\n", *best_val);
        log_info("         after LNS = %d\n", tour_cost);
        *best_val = tour_cost;
    }

    log_info("Checking solution...\n");
    rval = GTSP_check_solution(data, x);
    abort_if(rval, "GTSP_check_solution failed");

    CLEANUP:
    return rval;
}

int build_x_from_tour(struct GTSP *data, struct Tour *tour, double *x)
{
    int rval = 0;
    int *edge_map = 0;

    int node_count = data->graph->node_count;
    int edge_count = data->graph->edge_count;

    edge_map = (int *) malloc(node_count * node_count * sizeof(int));
    abort_if(!edge_map, "could not allocate edge_map");

    rval = build_edge_map(data, edge_map);
    abort_if(rval, "build_edge_map failed");

    for (int i = 0; i < edge_count; i++)
        x[i] = 0.0;

    int next_vertex = tour[0].next;
    int current_vertex = tour[0].vertex;
    for (int i = 0; i < data->cluster_count; i++)
    {
        int from = current_vertex;
        int to = tour[next_vertex].vertex;
        current_vertex = tour[next_vertex].vertex;
        next_vertex = tour[next_vertex].next;

        x[edge_map[from * node_count + to]] = 1.0;
    }

    CLEANUP:
    if (edge_map) free(edge_map);
    return rval;
}

int inital_tour_value(struct GTSP *data, int *tour_cost, double *x)
{
    int rval = 0;

    int cluster_count = data->cluster_count;

    struct Tour *tour = 0;
    int *uncovered_sets = 0;
    int *cluster_in_tour = 0;

    tour = (struct Tour *) malloc((cluster_count + 1) * sizeof(struct Tour));
    uncovered_sets = (int *) malloc((cluster_count - 1) * sizeof(int));
    cluster_in_tour = (int *) malloc(cluster_count * sizeof(int));
    abort_if(!tour, "could not allocate tour");
    abort_if(!uncovered_sets, "could not allocate uncovered_sets");
    abort_if(!cluster_in_tour, "could not allocate cluster_in_tour");

    int cluster_num = 0;
    for (int i = 0; i < cluster_count; i++)
    {
        cluster_in_tour[i] = 0;
        if (data->node_to_cluster[0] != i)
        {
            uncovered_sets[cluster_num] = i;
            cluster_num++;
        }
    }

    int new_vertex = 1, cost;

    tour[0].vertex = 0;
    cluster_in_tour[0] = 1;

    while (new_vertex < data->cluster_count)
    {
        int min_dist_vertex = -1;
        int min_cost = INT_MAX;

        for (int i = 1; i < data->graph->node_count; i++)
        {
            if (cluster_in_tour[data->node_to_cluster[i]]) continue;

            for (int k = 0; k < new_vertex; k++)
            {
                cost = data->dist_matrix[i][tour[k].vertex];
                if (cost < min_cost)
                {
                    min_cost = cost;
                    min_dist_vertex = i;
                }
            }
        }

        assert(min_dist_vertex >= 0);
        int cluster_to_insert = data->node_to_cluster[min_dist_vertex];
        cluster_in_tour[cluster_to_insert] = 1;

        min_cost = INT_MAX;

        int insertion_cost = -1, best_pose = -1, best_vertex = -1;

        for (int k = 0; k < data->clusters[cluster_to_insert].size; k++)
        {
            for (int j = 0; j < new_vertex; j++)
            {
                int vertex_to_insert = data->clusters[cluster_to_insert]
                        .nodes[k];
                if (new_vertex == 1)
                {
                    int vertex1 = tour[0].vertex;
                    insertion_cost =
                            data->dist_matrix[vertex_to_insert][vertex1] +
                                    data->dist_matrix[vertex1][vertex_to_insert];
                }
                else
                {
                    int vertex1 = tour[j].vertex;
                    int vertex2 = tour[tour[j].next].vertex;
                    insertion_cost =
                            data->dist_matrix[vertex1][vertex_to_insert] +
                                    data->dist_matrix[vertex_to_insert][vertex2] -
                                    data->dist_matrix[vertex1][vertex2];
                }
                if (insertion_cost < min_cost)
                {
                    min_cost = insertion_cost;
                    best_pose = j;
                    best_vertex = vertex_to_insert;
                }
            }
        }

        tour[new_vertex].vertex = best_vertex;
        tour[new_vertex].prev = best_pose;

        if (new_vertex == 1)
        {
            tour[new_vertex].next = best_pose;
            tour[best_pose].prev = new_vertex;
        }
        else
        {
            int temp_vertex = tour[best_pose].next;
            tour[new_vertex].next = temp_vertex;
            tour[temp_vertex].prev = new_vertex;
        }
        tour[best_pose].next = new_vertex;

        new_vertex += 1;
    }

    tour[data->cluster_count].vertex = 0;

    rval = large_neighborhood_search(tour, data, tour_cost);
    abort_if(rval, "large_neighborhood_search failed");

    log_info("Initial upper-bound: %d \n", *tour_cost);

    rval = build_x_from_tour(data, tour, x);
    abort_if(rval, "build_x_from_tour failed");

    CLEANUP:
    if (tour) free(tour);
    if (cluster_in_tour) free(cluster_in_tour);
    if (uncovered_sets) free(uncovered_sets);
    return rval;
}

int optimize_vertex_in_cluster(struct Tour *tour, struct GTSP *data)
{
    int current_cluster;
    int insertion_cost;

    int **dist_matrix = data->dist_matrix;
    int cluster_count = data->cluster_count;
    struct Cluster *vertex_set = data->clusters;

    for (int i = 0; i < cluster_count; i++)
    {
        int vertex = tour[i].vertex;
        int prev_vertex = tour[tour[i].prev].vertex;
        int next_vertex = tour[tour[i].next].vertex;

        current_cluster = data->node_to_cluster[vertex];

        insertion_cost = dist_matrix[prev_vertex][vertex] +
                dist_matrix[vertex][next_vertex];

        for (int j = 0; j < vertex_set[current_cluster].size; j++)
        {
            int vertex_in_cluster = vertex_set[current_cluster].nodes[j];
            int cost = dist_matrix[vertex_in_cluster][prev_vertex] +
                    dist_matrix[vertex_in_cluster][next_vertex];
            if (insertion_cost > cost)
            {
                insertion_cost = cost;
                tour[i].vertex = vertex_in_cluster;
            }
        }
    }

    return 0;
}

int two_opt(struct Tour *tour, struct GTSP *data)
{
    int **dist_matrix = data->dist_matrix;

    for (int i = 0; i < data->cluster_count; i++)
    {
        int v1 = tour[i].vertex;
        int v2 = tour[tour[i].prev].vertex;
        int v3 = tour[tour[i].next].vertex;
        int v4 = tour[tour[tour[i].next].next].vertex;

        int current_cost = dist_matrix[v2][v1] + dist_matrix[v3][v4];
        int temp_cost = dist_matrix[v2][v3] + dist_matrix[v1][v4];

        if (current_cost > temp_cost)
        {
            int temp_next = tour[i].next;
            int temp_prev = tour[i].prev;

            tour[i].next = tour[temp_next].next;
            tour[i].prev = temp_next;

            tour[tour[temp_next].next].prev = i;

            tour[temp_next].next = i;
            tour[temp_next].prev = temp_prev;

            tour[temp_prev].next = temp_next;

        }
    }

    return 0;
}

int large_neighborhood_search(
        struct Tour *tour, struct GTSP *data, int *tour_cost)
{
    int rval = 0;

    int cluster_count = data->cluster_count;
    int *clusters = data->node_to_cluster;
    int **dist_matrix = data->dist_matrix;
    struct Cluster *vertex_set = data->clusters;

    //LNS starts
    for (int iter = 0; iter < 2000; iter++)
    {
        //Delete a random vertex

        int delete_vertex = rand() % (cluster_count - 1) + 1;

        int prev_vertex = tour[delete_vertex].prev;
        int next_vertex = tour[delete_vertex].next;

        tour[prev_vertex].next = next_vertex;
        tour[next_vertex].prev = prev_vertex;

        int cluster_to_insert = clusters[tour[delete_vertex].vertex];

        int best_pose = -1;
        int best_vertex = -1;
        int min_cost = INT_MAX;

        for (int i = 0; i < vertex_set[cluster_to_insert].size; i++)
        {
            int vertex_to_insert = vertex_set[cluster_to_insert].nodes[i];

            int next_edge = tour[0].next;
            for (int j = 1; j < cluster_count; j++)
            {
                int vertex1 = tour[next_edge].vertex;
                int vertex2 = tour[tour[next_edge].next].vertex;

                int insert_cost = dist_matrix[vertex1][vertex_to_insert] +
                        dist_matrix[vertex_to_insert][vertex2] -
                        dist_matrix[vertex1][vertex2];

                if (insert_cost < min_cost)
                {
                    min_cost = insert_cost;
                    best_pose = next_edge;
                    best_vertex = vertex_to_insert;
                }

                next_edge = tour[next_edge].next;
            }
        }

        assert(best_pose >= 0);
        assert(best_vertex >= 0);

        next_vertex = tour[best_pose].next;
        tour[delete_vertex].prev = best_pose;
        tour[delete_vertex].vertex = best_vertex;
        tour[delete_vertex].next = next_vertex;
        tour[best_pose].next = delete_vertex;
        tour[next_vertex].prev = delete_vertex;

        rval = optimize_vertex_in_cluster(tour, data);
        abort_if(rval, "optimize_vertex_in_cluster failed");
    }

    rval = K_opt(tour, data);
    abort_if(rval, "two_opt failed");

    rval = two_opt(tour, data);
    abort_if(rval, "two_opt failed");

    *tour_cost = list_length(tour, data);

    CLEANUP:
    //if (vertex_seq) free(vertex_seq);
    return rval;
}

int tour_length(int *tour, struct GTSP *data)
{
    int tour_cost = 0;
    for (int i = 0; i < data->cluster_count; i++)
    {
        if (i == data->cluster_count - 1)
            tour_cost += data->dist_matrix[tour[i]][tour[0]];
        else
            tour_cost += data->dist_matrix[tour[i]][tour[i + 1]];

    }
    return tour_cost;
}

int list_length(struct Tour *tour, struct GTSP *data)
{
    int tour_cost = 0;
    for (int i = 0; i < data->cluster_count; i++)
    {
        int vertex1 = tour[i].vertex;
        int vertex2 = tour[tour[i].next].vertex;
        tour_cost += data->dist_matrix[vertex1][vertex2];
    }
    return tour_cost;
}

void print_tour(int *tour, struct GTSP *data)
{
    for (int i = 0; i < data->cluster_count; i++)
    {
        printf("%d\t", tour[i]);
    }

    printf("\n");
}

void print_list(struct Tour *tour, struct GTSP *data)
{
    printf("%d\t", tour[0].vertex);
    int vertex_next = tour[0].next;

    for (int i = 1; i < data->cluster_count; i++)
    {
        printf("%d\t", tour[vertex_next].vertex);
        vertex_next = tour[vertex_next].next;
    }

    printf("\n");
}

int K_opt(struct Tour *tour, struct GTSP *data)
{

    int cluster_count = data->cluster_count;
    int l;

    for (int i = 0; i < cluster_count - 3; i++)
    {
        int location_in_path = 0;

        for (int k = 0; k < i; k++)
            location_in_path = tour[location_in_path].next;

        int current_vertex = tour[location_in_path].vertex;

        for (int k = 3; k < cluster_count - i - 2; k++)
        {
            int first_next_location = tour[location_in_path].next;
            int first_next_vertex = tour[tour[location_in_path].next].vertex;
            int next_location = tour[location_in_path].next;
            for (l = 0; l < k; l++)
            {
                next_location = tour[next_location].next;
            }
            int next_vertex = tour[next_location].vertex;
            int next_next_location = tour[next_location].next;

            int next_next_vertex = tour[next_next_location].vertex;

            if (next_next_vertex == current_vertex) break;

            int cost1 = data->dist_matrix[current_vertex][first_next_vertex] +
                    data->dist_matrix[next_vertex][next_next_vertex];
            int cost2 = data->dist_matrix[current_vertex][next_vertex] +
                    data->dist_matrix[first_next_vertex][next_next_vertex];

            if (cost2 < cost1)
            {
                int tmp_location = next_location;
                for (int m = 0; m < l + 1; m++)
                {
                    int tmp_vertex = tour[tmp_location].next;
                    tour[tmp_location].next = tour[tmp_location].prev;
                    tour[tmp_location].prev = tmp_vertex;
                    tmp_location = tour[tmp_location].next;
                }

                tour[location_in_path].next = next_location;
                tour[next_location].prev = location_in_path;

                tour[first_next_location].next = next_next_location;

                tour[next_next_location].prev = first_next_location;

            }
        }
    }
    return 0;
}

int build_tour_from_x(struct GTSP *data, struct Tour *tour, double *x)
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

    for (int i = 0; i < node_count; i++)
        graph->nodes[i].mark = 0;

    for (int i = 0; i < data->cluster_count; i++)
        cluster_mark[i] = 0;

    int initial;
    for (initial = 0; initial < edge_count; initial++)
        if (x[initial] > 1.0 - EPSILON) break;

    initial = graph->edges[initial].from->index;

    abort_if(initial == edge_count, "no initial node");

    stack[stack_top++] = &graph->nodes[initial];
    graph->nodes[initial].mark = 1;

    tour[0].vertex = graph->nodes[initial].index;
    int next_vertex = 1;

    while (stack_top > 0)
    {
        struct Node *n = stack[--stack_top];
        cluster_mark[data->node_to_cluster[n->index]]++;

        for (int i = 0; i < n->degree; i++)
        {
            struct Adjacency *adj = &n->adj[i];
            struct Node *neighbor = adj->neighbor;

            if (neighbor->mark) continue;
            if (x[adj->edge->index] < EPSILON) continue;

            stack[stack_top++] = neighbor;
            tour[next_vertex].vertex = neighbor->index;
            next_vertex++;
            neighbor->mark = 1;

        }
    }

    for (int i = 0; i < data->cluster_count; i++)
    {
        if (i == 0)
        {
            tour[i].prev = data->cluster_count - 1;
        }
        else
        {
            tour[i].prev = i - 1;
        }
        if (i == data->cluster_count - 1)
        {
            tour[i].next = 0;
        }
        else
        {
            tour[i].next = i + 1;
        }
    }

    CLEANUP:
    if (stack) free(stack);
    if (cluster_mark) free(cluster_mark);
    return rval;
}

