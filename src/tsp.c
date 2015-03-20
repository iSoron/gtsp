#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "tsp.h"
#include "util.h"
#include "geometry.h"

int TSP_init_data(struct TSPData *data)
{
    data->node_count = 0;
    data->edge_count = 0;
    data->edge_weights = 0;
    data->edge_list = 0;

    return 0;
}

void TSP_free_data(struct TSPData *data)
{
    if(!data) return;
    if (data->edge_list) free(data->edge_list);
    if (data->edge_weights) free(data->edge_weights);
}

int TSP_init_lp(struct LP *lp, struct TSPData *data)
{
    int node_count = data->node_count;
    int edge_count = data->edge_count;
    int *edge_list = data->edge_list;
    int *edge_weights = data->edge_weights;

    int rval = 0;

    /* Build a row for each degree equation */
    for (int i = 0; i < node_count; i++)
    {
        rval = LP_new_row(lp, 'E', 2.0);
        abort_if(rval, "LP_new_row failed");
    }

    /* Build a column for each edge_index of the graph */
    double lb = 0.0;
    double ub = 1.0;
    int cmatbeg = 0;
    double cmatval[] = {1.0, 1.0};
    for (int j = 0; j < edge_count; j++)
    {
        double obj = (double) edge_weights[j];
        int cmatind[] = {edge_list[2 * j], edge_list[2 * j + 1]};

        rval = LP_add_cols(lp, 1, 2, &obj, &cmatbeg, cmatind, cmatval, &lb,
                &ub);

        abort_if(rval, "LP_add_cols failed");
    }

    CLEANUP:
    return rval;
}

int TSP_find_violated_subtour_elimination_cut(
        struct LP *lp,
        struct TSPData *data)
{
    int ncount = data->node_count;
    int edge_count = data->edge_count;
    int *edges = data->edge_list;

    int rval = 0;
    int is_infeasible = 0;

    double *x = 0;
    int *delta = 0;
    int *marks = 0;
    int *island_nodes = 0;
    int *island_start = 0;
    int *island_sizes = 0;

    struct Graph G;
    graph_init(&G);

    rval = LP_optimize(lp, &is_infeasible);
    abort_if(rval, "LP_optimize failed");
    abort_if(is_infeasible, "LP is infeasible");

    rval = graph_build(ncount, edge_count, edges, 0, &G);
    abort_if(rval, "graph_build failed");

    x = (double *) malloc(edge_count * sizeof(double));
    delta = (int *) malloc(edge_count * sizeof(int));
    marks = (int *) malloc(ncount * sizeof(int));
    abort_if(!x, "Could not allocate memory for x");
    abort_if(!delta, "Could not allocate memory for delta");
    abort_if(!marks, "Could not allocate memory for marks");

    island_nodes = (int *) malloc(ncount * sizeof(int));
    island_start = (int *) malloc(ncount * sizeof(int));
    island_sizes = (int *) malloc(ncount * sizeof(int));
    abort_if(!island_nodes, "Could not allocate memory for island_nodes");
    abort_if(!island_start, "Could not allocate memory for island_start");
    abort_if(!island_sizes, "Could not allocate memory for island_sizes");

    for (int i = 0; i < ncount; i++)
        marks[i] = 0;

    rval = LP_get_x(lp, x);
    abort_if(rval, "LP_get_x failed");

    int round = 0;
    int delta_count = 0;
    int island_count = 0;

    while (!TSP_is_graph_connected(&G, x, &island_count, island_sizes,
            island_start, island_nodes))
    {
        log_verbose("Adding %d BNC_solve_node inequalities...\n", island_count);
        for (int i = 0; i < island_count; i++)
        {
            get_delta(island_sizes[i], island_nodes + island_start[i],
                    edge_count, edges, &delta_count, delta, marks);

            rval = TSP_add_subtour_elimination_cut(lp, delta_count, delta);
        }

        log_verbose("Reoptimizing (round %d)...\n", ++round);
        abort_if(rval, "TSP_add_subtour_elimination_cut failed");

        rval = LP_optimize(lp, &is_infeasible);
        abort_if(rval, "LP_optimize failed");

        if(is_infeasible) goto CLEANUP;

        double objval = 0;
        rval = LP_get_obj_val(lp, &objval);
        abort_if(rval, "LP_get_obj_val failed");

        rval = LP_get_x(lp, x);
        abort_if(rval, "LP_get_x failed");
    }

    log_verbose("    graph is connected\n");

    CLEANUP:
    graph_free(&G);
    if (x) free(x);
    if (island_nodes) free(island_nodes);
    if (delta) free(delta);
    if (marks) free(marks);
    return rval;
}

int TSP_add_subtour_elimination_cut(struct LP *lp, int delta_length, int *delta)
{
    int rval = 0;
    char sense = 'G';
    double rhs = 2.0;
    int rmatbeg = 0;

    double *rmatval;
    int *rmatind = delta;

    rmatval = (double *) malloc(delta_length * sizeof(double));
    abort_if(!rmatval, "out of memory for rmatval");

    for (int i = 0; i < delta_length; i++)
        rmatval[i] = 1.0;

    rval = LP_add_rows(lp, 1, delta_length, &rhs, &sense, &rmatbeg, rmatind,
            rmatval);

    abort_if(rval, "LP_add_rows failed");

    CLEANUP:
    if (rmatval) free(rmatval);
    return rval;
}

int TSP_is_graph_connected(
        struct Graph *G,
        double *x,
        int *island_count,
        int *island_sizes,
        int *island_start,
        int *island_nodes)
{
    for (int i = 0; i < G->node_count; i++)
    {
        G->nodes[i].mark = 0;
        island_nodes[i] = -1;
    }

    int k = 0, current_island = 0;

    for (int i = 0; i < G->node_count; i++)
    {
        if (G->nodes[i].mark != 0) continue;

        island_sizes[current_island] = 0;

        graph_dfs(i, G, x, island_sizes + current_island, island_nodes + k);

        island_start[current_island] = k;
        k += island_sizes[current_island];

        current_island++;
    }

    (*island_count) = current_island;

    return (*island_count == 1);
}

int TSP_find_closest_neighbor_tour(
        int start,
        int node_count,
        int edge_count,
        int *edges,
        int *elen,
        int *path_length)
{
    int rval;
    int current_node = start;

    struct Graph G;
    graph_init(&G);

    rval = graph_build(node_count, edge_count, edges, 0, &G);
    abort_if(rval, "graph_build failed");

    for (int j = 0; j < node_count; j++)
        G.nodes[j].mark = 0;

    for (int j = 0; j < node_count; j++)
    {
        if (j == node_count - 1)
            G.nodes[start].mark = 0;

        struct Node *pn = &G.nodes[current_node];
        pn->mark = 1;

        int closest_neighbor = -1;
        int closest_edge_length = 10000000;

        for (int i = 0; i < pn->degree; i++)
        {
            int edge = pn->adj[i].edge_index;
            int neighbor = pn->adj[i].neighbor_index;

            if (G.nodes[neighbor].mark == 1) continue;
            if (elen[edge] > closest_edge_length) continue;

            closest_neighbor = neighbor;
            closest_edge_length = elen[edge];
        }

        *path_length += closest_edge_length;
        current_node = closest_neighbor;
    }

    CLEANUP:
    graph_free(&G);
    return rval;
}

int TSP_read_problem(char *filename, struct TSPData *data)
{
    int *p_node_count = &data->node_count;
    int *p_edge_count = &data->edge_count;
    int **p_edge_list = &data->edge_list;
    int **p_edge_weights = &data->edge_weights;

    struct _IO_FILE *f = (struct _IO_FILE *) NULL;
    int i, j, end1, end2, w, rval = 0, node_count, edge_count;
    int *edge_list = (int *) NULL, *edge_weights = (int *) NULL;
    double *x = (double *) NULL, *y = (double *) NULL;

    if (filename)
    {
        if ((f = fopen(filename, "r")) == NULL)
        {
            fprintf(stderr, "Unable to open %s for input\n", filename);
            rval = 1;
            goto CLEANUP;
        }
    }

    if (filename && GEOMETRIC_DATA == 0)
    {
        if (fscanf(f, "%d %d", &node_count, &edge_count) != 2)
        {
            fprintf(stderr, "Input file %s has invalid format\n", filename);
            rval = 1;
            goto CLEANUP;
        }

        printf("Nodes: %d  Edges: %d\n", node_count, edge_count);
        fflush(stdout);

        edge_list = (int *) malloc(2 * edge_count * sizeof(int));
        if (!edge_list)
        {
            fprintf(stderr, "out of memory for edge_list\n");
            rval = 1;
            goto CLEANUP;
        }

        edge_weights = (int *) malloc(edge_count * sizeof(int));
        if (!edge_weights)
        {
            fprintf(stderr, "out of memory for edge_weights\n");
            rval = 1;
            goto CLEANUP;
        }

        for (i = 0; i < edge_count; i++)
        {
            if (fscanf(f, "%d %d %d", &end1, &end2, &w) != 3)
            {
                fprintf(stderr, "%s has invalid input format\n", filename);
                rval = 1;
                goto CLEANUP;
            }
            edge_list[2 * i] = end1;
            edge_list[2 * i + 1] = end2;
            edge_weights[i] = w;
        }
    }
    else
    {
        if (filename)
        {
            if (fscanf(f, "%d", &node_count) != 1)
            {
                fprintf(stderr, "Input file %s has invalid format\n", filename);
                rval = 1;
                goto CLEANUP;
            }
        }
        else
        {
            node_count = NODE_COUNT_RAND;
        }

        x = (double *) malloc(node_count * sizeof(double));
        y = (double *) malloc(node_count * sizeof(double));
        if (!x || !y)
        {
            fprintf(stdout, "out of memory for x or y\n");
            rval = 1;
            goto CLEANUP;
        }

        if (filename)
        {
            for (i = 0; i < node_count; i++)
            {
                if (fscanf(f, "%lf %lf", &x[i], &y[i]) != 2)
                {
                    fprintf(stderr, "%s has invalid input format\n", filename);
                    rval = 1;
                    goto CLEANUP;
                }
            }
        }
        else
        {
            rval = generate_random_points_2d(node_count, GRID_SIZE_RAND, x, y);
            if (rval)
            {
                fprintf(stderr, "generate_random_points_2d failed\n");
                goto CLEANUP;
            }

            printf("%d\n", node_count);
            for (i = 0; i < node_count; i++)
            {
                printf("%.0f %.0f\n", x[i], y[i]);
            }
            printf("\n");
        }

        edge_count = (node_count * (node_count - 1)) / 2;
        log_verbose("Complete graph: %d nodes, %d edges\n", node_count,
                edge_count);

        edge_list = (int *) malloc(2 * edge_count * sizeof(int));
        if (!edge_list)
        {
            fprintf(stderr, "out of memory for edge_list\n");
            rval = 1;
            goto CLEANUP;
        }

        edge_weights = (int *) malloc(edge_count * sizeof(int));
        if (!edge_weights)
        {
            fprintf(stderr, "out of memory for edge_weights\n");
            rval = 1;
            goto CLEANUP;
        }

        edge_count = 0;
        for (i = 0; i < node_count; i++)
        {
            for (j = i + 1; j < node_count; j++)
            {
                edge_list[2 * edge_count] = i;
                edge_list[2 * edge_count + 1] = j;
                edge_weights[edge_count] = get_euclidean_distance(x, y, i, j);
                edge_count++;
            }
        }
    }

    *p_node_count = node_count;
    *p_edge_count = edge_count;
    *p_edge_list = edge_list;
    *p_edge_weights = edge_weights;

    CLEANUP:
    if (f) fclose(f);
    if (x) free(x);
    if (y) free(y);
    return rval;
}

int TSP_add_cutting_planes(struct LP *lp, struct TSPData *data)
{
    int rval = 0;

    rval = TSP_find_violated_subtour_elimination_cut(lp, data);
    abort_if (rval, "TSP_find_violated_subtour_elimination_cut failed\n");

    CLEANUP:
    return rval;
}

double TSP_find_initial_solution(struct TSPData *data)
{
    double best_val = 1e99;

    log_verbose("Finding closest neighbor_index tour\n");
    for (int i = 0; i < data->node_count; i++)
    {
        int path_length = 0;

        TSP_find_closest_neighbor_tour(i, data->node_count, data->edge_count,
                data->edge_list, data->edge_weights, &path_length);

        if (best_val > path_length) best_val = path_length;
    }

    log_verbose("    length = %lf\n", best_val);

    return best_val;
}