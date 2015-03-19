#include <stdio.h>
#include <stdlib.h>
#include "gtsp.h"
#include "geometry.h"
#include "util.h"

int GTSP_init_data(struct GTSP *data)
{
    data->node_count = 0;
    data->edge_count = 0;
    data->edges = 0;
    data->clusters = 0;
    data->cluster_count = 0;
    data->x_coordinates = 0;
    data->y_coordinates = 0;
    return 0;
}

void GTSP_free(struct GTSP *data)
{
    if (!data) return;
    if (data->edges) free(data->edges);
    if (data->clusters) free(data->clusters);
    if (data->x_coordinates) free(data->x_coordinates);
    if (data->y_coordinates) free(data->y_coordinates);
}

int GTSP_create_random_problem(
        int node_count, int cluster_count, int grid_size, struct GTSP *data)
{
    int rval = 0;

    struct Edge *edges = 0;
    int *clusters = 0;

    double *x_coords = 0;
    double *y_coords = 0;

    int edge_count = (node_count * (node_count - 1)) / 2;

    edges = (struct Edge *) malloc(edge_count * sizeof(struct Edge));
    clusters = (int *) malloc(node_count * sizeof(int));
    ABORT_IF (!edges, "could not allocate data->edges\n");
    ABORT_IF (!clusters, "could not allocate clusters\n");

    x_coords = (double *) malloc(node_count * sizeof(double));
    y_coords = (double *) malloc(node_count * sizeof(double));
    ABORT_IF (!x_coords, "could not allocate x_coords\n");
    ABORT_IF (!y_coords, "could not allocate y_coords\n");

    rval = generate_random_clusters_2d(node_count, cluster_count, grid_size,
            x_coords, y_coords, clusters);
    ABORT_IF(rval, "generate_random_clusters_2d failed");

    int current_edge = 0;
    for (int i = 0; i < edge_count; i++)
        for (int j = i + 1; j < node_count; j++)
        {
            edges[current_edge].from = i;
            edges[current_edge].to = j;
            edges[current_edge].weight =
                    get_euclidean_distance(x_coords, y_coords, i, j);

            current_edge++;
        }

    data->node_count = node_count;
    data->edge_count = edge_count;
    data->edges = edges;
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

    int node_count = data->node_count;
    int edge_count = data->edge_count;
    int cluster_count = data->cluster_count;
    int *clusters = data->clusters;
    struct Edge *edges = data->edges;

    for (int i = 0; i < node_count; i++)
    {
        rval = LP_new_row(lp, 'E', 0.0);
        ABORT_IF(rval, "LP_new_row failed");
    }

    for (int i = 0; i < cluster_count; i++)
    {
        rval = LP_new_row(lp, 'E', 1.0);
        ABORT_IF(rval, "LP_new_row failed");
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
        ABORT_IF(rval, "LP_add_cols failed");
    }

    for (int i = 0; i < edge_count; i++)
    {
        double obj = (double) edges[i].weight;
        double cmatval[] = {1.0, 1.0};
        int cmatind[] = {edges[i].from, edges[i].to};

        rval = LP_add_cols(lp, 1, 2, &obj, &cmatbeg, cmatind, cmatval, &lb,
                &ub);
        ABORT_IF(rval, "LP_add_cols failed");
    }

    CLEANUP:
    return rval;
}

int GTSP_write_data(struct GTSP *data, char *filename)
{
    int rval = 0;

    FILE *file;

    file = fopen(filename, "w");
    ABORT_IF(!file, "could not open file");

    fprintf(file, "%d %d\n", data->node_count, data->cluster_count);

    for (int i = 0; i < data->node_count; i++)
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

    struct Edge *edges = data->edges;
    int node_count = data->node_count;

    FILE *file;
    file = fopen(filename, "w");
    ABORT_IF(!file, "could not open file");

    int positive_edge_count = 0;
    for (int i = 0; i < data->edge_count; i++)
        if (x[i + node_count] > 0.5)
            positive_edge_count++;

    fprintf(file, "%d\n", positive_edge_count);

    for (int i = 0; i < data->edge_count; i++)
        if (x[i + node_count] > 0.5)
            fprintf(file, "%d %d\n", edges[i].from, edges[i].to);

    CLEANUP:
    if (file) fclose(file);
    return rval;
}