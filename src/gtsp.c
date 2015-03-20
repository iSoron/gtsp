#include <stdio.h>
#include <stdlib.h>
#include "gtsp.h"
#include "geometry.h"
#include "util.h"

int GTSP_init_data(struct GTSP *data)
{
    data->clusters = 0;
    data->cluster_count = 0;
    data->x_coordinates = 0;
    data->y_coordinates = 0;
    graph_init(data->graph);

    return 0;
}

void GTSP_free(struct GTSP *data)
{
    if (!data) return;
    if (data->graph) graph_free(data->graph);
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

    graph = (struct Graph*) malloc(sizeof(struct Graph));
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

    data->graph = graph;
    data->clusters = clusters;
    data->cluster_count = cluster_count;
    data->x_coordinates = x_coords;
    data->y_coordinates = y_coords;

    CLEANUP:
    if (rval)
    {
        if (edges) free(edges);
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

int GTSP_write_data(struct GTSP *data, char *filename)
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
        if (x[i + node_count] > 0.5)
            positive_edge_count++;

    fprintf(file, "%d\n", positive_edge_count);

    for (int i = 0; i < edge_count; i++)
        if (x[i + node_count] > 0.5)
            fprintf(file, "%d %d\n", edges[i].from->index, edges[i].to->index);

    CLEANUP:
    if (file) fclose(file);
    return rval;
}

