#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lp.h"
#include "branch_and_cut.h"
#include "tsp.h"
#include "util.h"

int bnc_is_integral(double *x, int length);

int bnc_find_best_branch_var(double *x, int length);

int bnc_solve_node(
        struct LP *lp, double *best_val, int ncount, int ecount, int *elist,
        int depth)
{
    int rval = 0;
    double *x = (double *) NULL;

    time_printf("Optimizing...\n");

    int is_infeasible;
    rval = lp_optimize(lp, &is_infeasible);
    ABORT_IF (rval, "lp_optimize failed\n");

    if (is_infeasible)
    {
        time_printf("Branch pruned by infeasibility.\n");
        goto CLEANUP;
    }

    double objval;
    rval = lp_get_obj_val(lp, &objval);
    ABORT_IF (rval, "lp_get_obj_val failed\n");

    time_printf("    objective value = %.2f\n", objval);

    if (objval > *best_val)
    {
        time_printf("Branch pruned by bound (%.2lf > %.2lf).\n", objval,
                *best_val);
        rval = 0;
        goto CLEANUP;
    }

    x = (double *) malloc(ecount * sizeof(double));
    ABORT_IF(!x, "could not allocate x\n");

    rval = lp_get_x(lp, x);
    ABORT_IF(rval, "lp_get_x failed\n");

    rval = TSP_add_cutting_planes(ncount, ecount, elist, lp);
    ABORT_IF(rval, "TSP_add_cutting_planes failed\n");

    rval = lp_optimize(lp, &is_infeasible);
    ABORT_IF (rval, "lp_optimize failed\n");

    rval = lp_get_obj_val(lp, &objval);
    ABORT_IF(rval, "lp_get_obj_val failed\n");

    rval = lp_get_x(lp, x);
    ABORT_IF(rval, "lp_get_x failed\n");

    if (bnc_is_integral(x, ecount))
    {
        time_printf("    solution is integral\n");

        if (objval < *best_val)
        {
            *best_val = objval;
            time_printf("Found a better integral solution:\n");
            time_printf("    objval = %.2lf **\n", objval);
        }
    } else
    {
        time_printf("    solution is fractional\n");
        rval = bnc_branch_node(lp, x, ncount, ecount, depth, best_val, elist);
        ABORT_IF(rval, "bnc_branch_node failed\n");
    }

    CLEANUP:
    if (x) free(x);
    return rval;
}

int bnc_branch_node(
        struct LP *lp, double *x, int ncount, int ecount, int depth,
        double *best_val, int *elist)
{
    int rval = 0;

    int best_index = bnc_find_best_branch_var(x, ecount);

    time_printf("Branching on variable x%d = %.6lf (depth %d)...\n", best_index,
            x[best_index], depth);

    time_printf("Fixing variable x%d to one...\n", best_index);
    rval = lp_change_bound(lp, best_index, 'L', 1.0);
    ABORT_IF(rval, "lp_change_bound failed\n");

    rval = bnc_solve_node(lp, best_val, ncount, ecount, elist, depth + 1);
    ABORT_IF(rval, "bnc_solve_node failed\n");

    rval = lp_change_bound(lp, best_index, 'L', 0.0);
    ABORT_IF(rval, "lp_change_bound failed\n");

    time_printf("Fixing variable x%d to zero...\n", best_index);
    rval = lp_change_bound(lp, best_index, 'U', 0.0);
    ABORT_IF(rval, "lp_change_bound failed\n");

    rval = bnc_solve_node(lp, best_val, ncount, ecount, elist, depth + 1);
    ABORT_IF(rval, "bnc_solve_node failed\n");

    rval = lp_change_bound(lp, best_index, 'U', 1.0);
    ABORT_IF(rval, "lp_change_bound failed\n");

    time_printf("Finished branching on variable %d\n", best_index);

    CLEANUP:
    return rval;
}

int bnc_find_best_branch_var(double *x, int length)
{
    int best_index = 0;
    double best_index_frac = 1.0;

    for (int j = 0; j < length; j++)
    {
        if (fabs(x[j] - 0.5) < best_index_frac)
        {
            best_index = j;
            best_index_frac = fabs(x[j] - 0.5);
        }
    }

    return best_index;
}

int bnc_is_integral(double *x, int length)
{
    int all_integral = 1;

    for (int j = 0; j < length; j++)
    {
        if (x[j] > LP_EPSILON && x[j] < 1.0 - LP_EPSILON)
        {
            all_integral = 0;
            break;
        }
    }

    return all_integral;
}

int bnc_init_lp(
        struct LP *lp, int node_count, int edge_count, int *edge_list,
        int *edge_weights)
{
    int rval = 0;
    time_printf("Initializing LP...\n");

    rval = lp_open(lp);
    ABORT_IF(rval, "lp_open failed\n");

    rval = lp_create(lp, "subtour");
    ABORT_IF(rval, "lp_create failed\n");

    rval = TSP_init_lp(node_count, lp, edge_count, edge_weights, edge_list);
    ABORT_IF(rval, "TSP_init_lp failed\n");

    rval = lp_write(lp, "subtour.lp");
    ABORT_IF(rval, "lp_write failed\n");

    CLEANUP:
    return rval;
}