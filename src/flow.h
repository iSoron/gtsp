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

#ifndef _PROJECT_FLOW_H_
#define _PROJECT_FLOW_H_

#include "graph.h"

int flow_find_augmenting_path(
        const struct Graph *graph,
        const double *residual_caps,
        struct Node *from,
        struct Node *to,
        int *path_length,
        struct Edge **path_edges,
        double *path_capacity);

int flow_find_max_flow(
        const struct Graph *digraph,
        const double *capacities,
        struct Node *from,
        struct Node *to,
        double *flow,
        double *value);

int flow_mark_reachable_nodes(
        const struct Graph *graph, double *residual_caps, struct Node *from);

#endif //_PROJECT_FLOW_H_
