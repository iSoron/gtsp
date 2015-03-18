#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "tsp.h"
#include "util.h"

int TSP_init_lp(
        int node_count, struct LP *lp, int edge_count, int *edge_weights,
        int *edge_list)
{
    int rval = 0;

    /* Build a row for each degree equation */
    for (int i = 0; i < node_count; i++)
    {
        rval = lp_new_row(lp, 'E', 2.0);
        ABORT_IF(rval, "lp_new_row failed\n");
    }

    /* Build a column for each edge of the graph */
    double lb = 0.0;
    double ub = 1.0;
    int cmatbeg = 0;
    double cmatval[] = {1.0, 1.0};
    for (int j = 0; j < edge_count; j++)
    {
        double obj = (double) edge_weights[j];
        int cmatind[] = {edge_list[2 * j], edge_list[2 * j + 1]};

        rval = lp_add_cols(lp, 1, 2, &obj, &cmatbeg, cmatind, cmatval, &lb,
                &ub);

        ABORT_IF(rval, "lp_add_cols failed\n");
    }

    CLEANUP:
    return rval;
}

int TSP_find_violated_subtour_elimination_cut(
        int ncount, int edge_count, int *edges, struct LP *lp)
{
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

    rval = lp_optimize(lp, &is_infeasible);
    ABORT_IF(rval, "lp_optimize failed\n");
    ABORT_IF(is_infeasible, "LP is infeasible\n");

    rval = graph_build(ncount, edge_count, edges, &G);
    ABORT_IF(rval, "graph_build failed\n");

    x = (double *) malloc(edge_count * sizeof(double));
    delta = (int *) malloc(edge_count * sizeof(int));
    marks = (int *) malloc(ncount * sizeof(int));
    ABORT_IF(!x, "Could not allocate memory for x");
    ABORT_IF(!delta, "Could not allocate memory for delta");
    ABORT_IF(!marks, "Could not allocate memory for marks");

    island_nodes = (int *) malloc(ncount * sizeof(int));
    island_start = (int *) malloc(ncount * sizeof(int));
    island_sizes = (int *) malloc(ncount * sizeof(int));
    ABORT_IF(!island_nodes, "Could not allocate memory for island_nodes");
    ABORT_IF(!island_start, "Could not allocate memory for island_start");
    ABORT_IF(!island_sizes, "Could not allocate memory for island_sizes");

    for (int i = 0; i < ncount; i++)
        marks[i] = 0;

    rval = lp_get_x(lp, x);
    ABORT_IF(rval, "lp_get_x failed\n");

    int round = 0;
    int delta_count = 0;
    int island_count = 0;

    while (!TSP_is_graph_connected(&G, x, &island_count, island_sizes,
            island_start, island_nodes))
    {
        time_printf("Adding %d bnc_solve_node inequalities...\n", island_count);
        for (int i = 0; i < island_count; i++)
        {
            get_delta(island_sizes[i], island_nodes + island_start[i],
                    edge_count, edges, &delta_count, delta, marks);

            rval = TSP_add_subtour_elimination_cut(lp, delta_count, delta);
        }

        time_printf("Reoptimizing (round %d)...\n", ++round);
        ABORT_IF(rval, "TSP_add_subtour_elimination_cut failed");

        rval = lp_optimize(lp, &is_infeasible);
        ABORT_IF(rval, "lp_optimize failed\n");
        ABORT_IF(is_infeasible, "LP is infeasible\n");

        double objval = 0;
        rval = lp_get_obj_val(lp, &objval);
        ABORT_IF(rval, "lp_get_obj_val failed\n");

        rval = lp_get_x(lp, x);
        ABORT_IF(rval, "lp_get_x failed\n");
    }

    time_printf("    graph is TSP_is_graph_connected\n");

    CLEANUP:
    graph_free(&G);
    if (x) free(x);
    if (island_nodes) free(island_nodes);
    if (delta) free(delta);
    if (marks) free(marks);
    return rval;
}

int TSP_add_subtour_elimination_cut(struct LP *lp, int deltacount, int *delta)
{
    int rval = 0;
    char sense = 'G';
    double rhs = 2.0;
    int rmatbeg = 0;

    double *rmatval;
    int *rmatind = delta;

    rmatval = (double *) malloc(deltacount * sizeof(double));
    ABORT_IF(!rmatval, "out of memory for rmatval\n");

    for (int i = 0; i < deltacount; i++)
        rmatval[i] = 1.0;

    rval = lp_add_rows(lp, 1, deltacount, &rhs, &sense, &rmatbeg, rmatind,
            rmatval);

    ABORT_IF(rval, "lp_add_rows failed");

    CLEANUP:
    if (rmatval) free(rmatval);
    return rval;
}

int TSP_is_graph_connected(
        struct Graph *G, double *x, int *island_count, int *island_sizes,
        int *island_start, int *island_nodes)
{
    for (int i = 0; i < G->node_count; i++)
    {
        G->node_list[i].mark = 0;
        island_nodes[i] = -1;
    }

    int k = 0, current_island = 0;

    for (int i = 0; i < G->node_count; i++)
    {
        if (G->node_list[i].mark != 0) continue;

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
        int start, int node_count, int edge_count, int *edges, int *elen,
        int *path_length)
{
    int rval;
    int current_node = start;

    struct Graph G;
    graph_init(&G);

    rval = graph_build(node_count, edge_count, edges, &G);
    ABORT_IF(rval, "graph_build failed\n");

    for (int j = 0; j < node_count; j++)
        G.node_list[j].mark = 0;

    for (int j = 0; j < node_count; j++)
    {
        if (j == node_count - 1)
            G.node_list[start].mark = 0;

        struct Node *pn = &G.node_list[current_node];
        pn->mark = 1;

        int closest_neighbor = -1;
        int closest_edge_length = 10000000;

        for (int i = 0; i < pn->deg; i++)
        {
            int edge = pn->adj[i].e;
            int neighbor = pn->adj[i].n;

            if (G.node_list[neighbor].mark == 1) continue;
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

int TSP_read_problem(
        char *filename, int *p_ncount, int *p_ecount, int **p_elist,
        int **p_elen)
{
    struct _IO_FILE *f = (struct _IO_FILE *) NULL;
    int i, j, end1, end2, w, rval = 0, ncount, ecount;
    int *elist = (int *) NULL, *elen = (int *) NULL;
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

    if (filename && geometric_data == 0)
    {
        if (fscanf(f, "%d %d", &ncount, &ecount) != 2)
        {
            fprintf(stderr, "Input file %s has invalid format\n", filename);
            rval = 1;
            goto CLEANUP;
        }

        printf("Nodes: %d  Edges: %d\n", ncount, ecount);
        fflush(stdout);

        elist = (int *) malloc(2 * ecount * sizeof(int));
        if (!elist)
        {
            fprintf(stderr, "out of memory for elist\n");
            rval = 1;
            goto CLEANUP;
        }

        elen = (int *) malloc(ecount * sizeof(int));
        if (!elen)
        {
            fprintf(stderr, "out of memory for elen\n");
            rval = 1;
            goto CLEANUP;
        }

        for (i = 0; i < ecount; i++)
        {
            if (fscanf(f, "%d %d %d", &end1, &end2, &w) != 3)
            {
                fprintf(stderr, "%s has invalid input format\n", filename);
                rval = 1;
                goto CLEANUP;
            }
            elist[2 * i] = end1;
            elist[2 * i + 1] = end2;
            elen[i] = w;
        }
    }
    else
    {
        if (filename)
        {
            if (fscanf(f, "%d", &ncount) != 1)
            {
                fprintf(stderr, "Input file %s has invalid format\n", filename);
                rval = 1;
                goto CLEANUP;
            }
        }
        else
        {
            ncount = ncount_rand;
        }

        x = (double *) malloc(ncount * sizeof(double));
        y = (double *) malloc(ncount * sizeof(double));
        if (!x || !y)
        {
            fprintf(stdout, "out of memory for x or y\n");
            rval = 1;
            goto CLEANUP;
        }

        if (filename)
        {
            for (i = 0; i < ncount; i++)
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
            rval = CO759_build_xy(ncount, x, y, gridsize_rand);
            if (rval)
            {
                fprintf(stderr, "CO759_build_xy failed\n");
                goto CLEANUP;
            }

            printf("%d\n", ncount);
            for (i = 0; i < ncount; i++)
            {
                printf("%.0f %.0f\n", x[i], y[i]);
            }
            printf("\n");
        }

        ecount = (ncount * (ncount - 1)) / 2;
        time_printf("Complete graph: %d nodes, %d edges\n", ncount, ecount);

        elist = (int *) malloc(2 * ecount * sizeof(int));
        if (!elist)
        {
            fprintf(stderr, "out of memory for elist\n");
            rval = 1;
            goto CLEANUP;
        }

        elen = (int *) malloc(ecount * sizeof(int));
        if (!elen)
        {
            fprintf(stderr, "out of memory for elen\n");
            rval = 1;
            goto CLEANUP;
        }

        ecount = 0;
        for (i = 0; i < ncount; i++)
        {
            for (j = i + 1; j < ncount; j++)
            {
                elist[2 * ecount] = i;
                elist[2 * ecount + 1] = j;
                elen[ecount] = euclid_edgelen(i, j, x, y);
                ecount++;
            }
        }
    }

    *p_ncount = ncount;
    *p_ecount = ecount;
    *p_elist = elist;
    *p_elen = elen;

    CLEANUP:
    if (f) fclose(f);
    if (x) free(x);
    if (y) free(y);
    return rval;
}

int TSP_add_cutting_planes(int ncount, int ecount, int *elist, struct LP *lp)
{
    int rval = 0;

    rval = TSP_find_violated_subtour_elimination_cut(ncount, ecount, elist, lp);
    ABORT_IF (rval, "TSP_find_violated_subtour_elimination_cut failed\n");

    CLEANUP:
    return rval;
}

double TSP_find_initial_solution(
        int *edge_weights, int *edge_list, int node_count, int edge_count)
{
    double best_val = 1e99;

    time_printf("Finding closest neighbor tour\n");
    for (int i = 0; i < node_count; i++)
    {
        int path_length = 0;

        TSP_find_closest_neighbor_tour(i, node_count, edge_count, edge_list,
                edge_weights, &path_length);

        if (best_val > path_length) best_val = path_length;
    }

    time_printf("    length = %lf\n", best_val);

    return best_val;
}