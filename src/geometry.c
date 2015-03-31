#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "geometry.h"
#include "util.h"

/* function for creating a random set of points in unit square */

int generate_random_points_2d(
        int node_count,
        int grid_size,
        double *x_coordinates,
        double *y_coordinates)
{
    int rval = 0, i, j, winner, x, y;
    int **hit = 0, *hit_count = 0;

    hit = (int **) malloc(grid_size * sizeof(int *));
    abort_if(!hit, "could not allocate hit");

    for (i = 0; i < grid_size; i++)
        hit[i] = 0;

    hit_count = (int *) malloc(grid_size * sizeof(int));
    abort_if(!hit_count, "could not allocate hit_count");

    for (i = 0; i < grid_size; i++)
        hit_count[i] = 0;

    for (i = 0; i < node_count; i++)
    {
        winner = 0;
        do
        {
            x = rand() % grid_size;
            y = rand() % grid_size;

            /* check to see if (x,y) is a duplicate point */
            for (j = 0; j < hit_count[x]; j++)
                if (hit[x][j] == y) break;

            if (j == hit_count[x])
            {
                void *tmp_ptr = (void *) hit[x];
                tmp_ptr = realloc(tmp_ptr, (hit_count[x] + 1) * sizeof(int));
                abort_if(!tmp_ptr, "could not reallocate hit_count");

                hit[x] = (int *) tmp_ptr;
                hit[x][hit_count[x]] = y;
                hit_count[x]++;
                winner = 1;
            }
        } while (!winner);

        x_coordinates[i] = (double) x;
        y_coordinates[i] = (double) y;
    }

    CLEANUP:

    if (hit)
    {
        for (i = 0; i < grid_size; i++)
            if (hit[i]) free(hit[i]);

        free(hit);
    }
    if (hit_count) free(hit_count);
    return rval;
}

int generate_random_clusters_2d(
        int node_count,
        int cluster_count,
        int grid_size,
        double *x_coordinates,
        double *y_coordinates,
        int *node_to_cluster)
{
    int rval = 0;

    rval = generate_random_points_2d(node_count, grid_size, x_coordinates,
            y_coordinates);
    abort_if(rval, "generate_random_points_2d failed");

    for (int i = 0; i < cluster_count; i++)
        node_to_cluster[i] = i;

    for (int i = cluster_count; i < node_count; i++)
    {
        int closest_point = 0;
        int closest_distance = INT_MAX;

        for (int j = 0; j < cluster_count; j++)
        {
            int distance =
                    get_euclidean_distance(x_coordinates, y_coordinates, i, j);

            if (distance < closest_distance)
            {
                closest_distance = distance;
                closest_point = j;
            }
        }

        node_to_cluster[i] = closest_point;
    }

    CLEANUP:
    return rval;
}

int generate_dist_matrix(
		int node_count,
        double *x_coordinates,
        double *y_coordinates, int** dist_matrix)
{
	int i,j;
    for (i = 0; i < node_count; i++){
		for (j = 0; j < node_count; j++){
			dist_matrix[i][j] = 
				get_euclidean_distance(x_coordinates, y_coordinates, i, j);
		}
	}
	return 0;
}

int get_euclidean_distance(
        double *x_coordinates,
        double *y_coordinates,
        int p1_index,
        int p2_index)
{
    double t1 = x_coordinates[p1_index] - x_coordinates[p2_index];
    double t2 = y_coordinates[p1_index] - y_coordinates[p2_index];
    return (int) (sqrt(t1 * t1 + t2 * t2) + 0.5);
}
