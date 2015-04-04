/* Copyright (c) 2015 Alinson Xavier
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

#ifndef _PROJECT_GEOMETRY_H_
#define _PROJECT_GEOMETRY_H_

int generate_random_points_2d(
        int node_count,
        int grid_size,
        double *x_coordinates,
        double *y_coordinates);

int generate_random_clusters_2d(
        int node_count,
        int cluster_count,
        int grid_size,
        double *x_coordinates,
        double *y_coordinates,
        int *node_to_cluster);

int get_euclidean_distance(
        double *x_coordinates,
        double *y_coordinates,
        int p1_index,
        int p2_index);

int generate_dist_matrix(
        int node_count,
        double *x_coordinates,
        double *y_coordinates,
        int **dist_matrix);

#endif //_PROJECT_GEOMETRY_H_
