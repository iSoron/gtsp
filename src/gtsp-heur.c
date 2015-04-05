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

#include <assert.h>
#include "util.h"
#include "gtsp.h"
#include "gtsp-heur.h"

static int large_neighborhood_search(
        struct Tour *tour, struct GTSP *data, int *tour_cost);

static int build_tour_from_x(struct GTSP *data, struct Tour *tour, double *x);

static int build_x_from_tour(struct GTSP *data, struct Tour *tour, double *x);

static int optimize_vertex_in_cluster(struct Tour *tour, struct GTSP *data);

static int two_opt(struct Tour *tour, struct GTSP *data);

//static int tour_length(int *tour, struct GTSP *data);

static int get_tour_length(struct Tour *tour, struct GTSP *data);

static int K_opt(struct Tour *tour, struct GTSP *data);

int GTSP_find_initial_tour(struct GTSP *data, int *tour_cost, double *x)
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
        assert(min_dist_vertex < data->graph->node_count);

        int cluster_to_insert = data->node_to_cluster[min_dist_vertex];
        cluster_in_tour[cluster_to_insert] = 1;

        min_cost = INT_MAX;

        int insertion_cost = -1, best_pose = -1, best_vertex = -1;

        for (int k = 0; k < data->clusters[cluster_to_insert].size; k++)
        {
            for (int j = 0; j < new_vertex; j++)
            {
                int vertex_to_insert = data->clusters[cluster_to_insert].nodes[k];
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

int GTSP_improve_solution(struct BNC *bnc, struct GTSP *data, double *x)
{
    int rval = 0;

    int tour_cost;
    double *best_val = &bnc->best_obj_val;
    double **best_x = &bnc->best_x;

    struct Tour *tour;
    tour = (struct Tour *) malloc(
            (data->cluster_count + 1) * sizeof(struct Tour));

    rval = build_tour_from_x(data, tour, x);
    abort_if(rval, "build_tour_from_x failed");

    for (int i = 0; i < data->cluster_count; i++)
    {
        log_debug("    %d\n", tour[i]);
    }

    rval = large_neighborhood_search(tour, data, &tour_cost);
    abort_if(rval, "large_neighborhood_search failed");

    if (tour_cost + EPSILON < *best_val)
    {
        log_info("Local search improved the integral solution\n");
        log_info("         before = %f\n", *best_val);
        log_info("         after  = %f\n", tour_cost);

        build_x_from_tour(data, tour, x);

        *best_val = tour_cost;
        *best_x = x;
    }

    CLEANUP:
    return rval;
}

int static optimize_vertex_in_cluster(struct Tour *tour, struct GTSP *data)
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

int static two_opt(struct Tour *tour, struct GTSP *data)
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

static int large_neighborhood_search(
        struct Tour *tour, struct GTSP *data, int *tour_cost)
{
    int rval = 0;

    int cluster_count = data->cluster_count;
    int *node_to_cluster = data->node_to_cluster;
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

        int cluster_to_insert = node_to_cluster[tour[delete_vertex].vertex];

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

    *tour_cost = get_tour_length(tour, data);

    CLEANUP:
    //if (vertex_seq) free(vertex_seq);
    return rval;
}

int static get_tour_length(struct Tour *tour, struct GTSP *data)
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

int static K_opt(struct Tour *tour, struct GTSP *data)
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

int static build_tour_from_x(struct GTSP *data, struct Tour *tour, double *x)
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

    int initial = -1;

    for (int i = 0; i < edge_count; i++)
    {
        int col = graph->edges[i].column;
        if (col < 0) continue;

        if (x[col] > 1.0 - EPSILON)
        {
            initial = i;
            break;
        }
    }

    initial = graph->edges[initial].from->index;

    abort_if(initial < 0, "no initial node");

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

            if (adj->edge->column < 0) continue;
            if (x[adj->edge->column] < EPSILON) continue;

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

static int build_x_from_tour(struct GTSP *data, struct Tour *tour, double *x)
{
    int rval = 0;
    int *edge_map = 0;

    int node_count = data->graph->node_count;
    int edge_count = data->graph->edge_count;

    edge_map = (int *) malloc(node_count * node_count * sizeof(int));
    abort_if(!edge_map, "could not allocate edge_map");

    rval = GTSP_build_edge_map(data, edge_map);
    abort_if(rval, "GTSP_build_edge_map failed");

    for (int i = 0; i < edge_count; i++)
        x[i] = 0.0;

    int k = 0;
    int next_vertex = tour[0].next;
    int current_vertex = tour[0].vertex;
    for (int i = 0; i < data->cluster_count; i++)
    {
        int from = current_vertex;
        int to = tour[next_vertex].vertex;
        current_vertex = tour[next_vertex].vertex;
        next_vertex = tour[next_vertex].next;
        struct Edge *e = &data->graph->edges[edge_map[from * node_count + to]];

        if (e->column < 0) e->column = k++;
        x[e->column] = 1.0;
    }

    CLEANUP:
    if (edge_map) free(edge_map);
    return rval;
}

