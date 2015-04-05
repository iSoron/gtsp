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

#ifndef _PROJECT_GTSP_H_
#define _PROJECT_GTSP_H_

#include "lp.h"
#include "graph.h"
#include "branch-and-cut.h"

struct Tour
{
    int vertex;
    int next;
    int prev;
};

struct Cluster
{
    /*
     * Number of nodes inside the cluster.
     */
    int size;

    /*
     * List of nodes inside the cluster.
     */
    int *nodes;
};

struct GTSP
{
    struct Graph *graph;

    int **dist_matrix;

    /*
     * Mapping between a node and the cluster which contains it. If a node
     * has index i and is contained in cluster k, then node_to_cluster[i] = k.
     */
    int *node_to_cluster;

    int cluster_count;
    struct Cluster *clusters;
};

int GTSP_create_random_problem(
        int node_count, int cluster_count, int grid_size, struct GTSP *data);

void GTSP_free(struct GTSP *data);

int GTSP_init_data(struct GTSP *data);

int GTSP_init_lp(struct LP *lp, struct GTSP *data);

int GTSP_add_cutting_planes(struct LP *lp, struct GTSP *data);

int GTSP_write_problem(struct GTSP *data, char *filename);

int GTSP_write_solution(struct GTSP *data, char *filename, double *x);

int GTSP_print_solution(struct GTSP *data, double *x);

int GTSP_solution_found(struct BNC *bnc, struct GTSP *data, double *x);

int GTSP_check_solution(struct GTSP *data, double *x);

int GTSP_read_solution(struct GTSP *gtsp, char *filename, double **p_x);

int GTSP_build_edge_map(struct GTSP *gtsp, int *edge_map);


#endif //_PROJECT_GTSP_H_
