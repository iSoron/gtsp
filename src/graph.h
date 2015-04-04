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

#ifndef _PROJECT_GRAPH_H_
#define _PROJECT_GRAPH_H_

#include "main.h"

struct Adjacency
{
    struct Edge *edge;
    struct Node *neighbor;
};

struct Node
{
    int mark;

    int index;
    int degree;

    struct Adjacency *adj;
};

struct Edge
{
    int index;
    int weight;

    struct Node *from;
    struct Node *to;

    struct Edge *reverse;
};

struct Graph
{
    int node_count;
    int edge_count;

    struct Edge *edges;
    struct Node *nodes;

    double *x_coordinates;
    double *y_coordinates;

    struct Adjacency *adj;
};

void graph_init(struct Graph *graph);

void graph_free(struct Graph *graph);

int graph_build(
        int node_count,
        int edge_count,
        int *edges,
        int is_directed,
        struct Graph *graph);

int get_cut_edges_from_marks(
        struct Graph *graph, int *cut_edges_count, struct Edge **cut_edges);

int graph_dump(const struct Graph *graph);

#endif
