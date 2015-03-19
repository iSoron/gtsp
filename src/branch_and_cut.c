#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lp.h"
#include "branch_and_cut.h"
#include "util.h"

static int BNC_solve_node(struct BNC *bnc, int depth);

static int BNC_branch_node(struct BNC *bnc, double *x, int depth);

static int BNC_is_integral(double *x, int num_cols);

static int BNC_find_best_branching_var(double *x, int num_cols);

int BNC_init(struct BNC *bnc)
{
    int rval = 0;

    bnc->lp = 0;
    bnc->problem_data = 0;
    bnc->problem_init_lp = 0;
    bnc->problem_add_cutting_planes = 0;

    bnc->best_x = 0;
    bnc->best_obj_val = 0;

    bnc->lp = (struct LP *) malloc(sizeof(struct LP));
    ABORT_IF(!bnc->lp, "could not allocate bnc->lp\n");

    CLEANUP:
    return rval;
}

void BNC_free(struct BNC *bnc)
{
    if (!bnc) return;
    if (bnc->lp)
    {
        LP_free(bnc->lp);
        free(bnc->lp);
    }
}

int BNC_init_lp(struct BNC *bnc)
{
    int rval = 0;
    time_printf("Initializing LP...\n");

    rval = LP_open(bnc->lp);
    ABORT_IF(rval, "LP_open failed\n");

    rval = LP_create(bnc->lp, "subtour");
    ABORT_IF(rval, "LP_create failed\n");

    rval = bnc->problem_init_lp(bnc->lp, bnc->problem_data);
    ABORT_IF(rval, "problem_init_lp failed\n");

    rval = LP_write(bnc->lp, "subtour.lp");
    ABORT_IF(rval, "LP_write failed\n");

    CLEANUP:
    return rval;
}

int BNC_solve(struct BNC *bnc)
{
    return BNC_solve_node(bnc, 1);
}

static int BNC_solve_node(struct BNC *bnc, int depth)
{
    struct LP *lp = bnc->lp;
    double *best_val = &bnc->best_obj_val;

    int rval = 0;
    double *x = (double *) NULL;

    time_printf("Optimizing...\n");

    int is_infeasible;
    rval = LP_optimize(lp, &is_infeasible);
    ABORT_IF (rval, "LP_optimize failed\n");

    if (is_infeasible)
    {
        time_printf("Branch pruned by infeasibility.\n");
        goto CLEANUP;
    }

    double objval;
    rval = LP_get_obj_val(lp, &objval);
    ABORT_IF (rval, "LP_get_obj_val failed\n");

    time_printf("    objective value = %.2f\n", objval);

    if (objval > *best_val)
    {
        time_printf("Branch pruned by bound (%.2lf > %.2lf).\n", objval,
                *best_val);
        rval = 0;
        goto CLEANUP;
    }

    int num_cols = LP_get_num_cols(lp);

    x = (double *) malloc(num_cols * sizeof(double));
    ABORT_IF(!x, "could not allocate x\n");

    rval = LP_get_x(lp, x);
    ABORT_IF(rval, "LP_get_x failed\n");

    if (bnc->problem_add_cutting_planes)
    {
        rval = bnc->problem_add_cutting_planes(lp, bnc->problem_data);
        ABORT_IF(rval, "problem_add_cutting_planes failed\n");
    }

    rval = LP_optimize(lp, &is_infeasible);
    ABORT_IF (rval, "LP_optimize failed\n");

    rval = LP_get_obj_val(lp, &objval);
    ABORT_IF(rval, "LP_get_obj_val failed\n");

    rval = LP_get_x(lp, x);
    ABORT_IF(rval, "LP_get_x failed\n");

    if (BNC_is_integral(x, num_cols))
    {
        time_printf("    solution is integral\n");

        if (objval < *best_val)
        {
            *best_val = objval;
            bnc->best_x = x;
            x = 0;

            time_printf("Found a better integral solution:\n");
            time_printf("    objval = %.2lf **\n", objval);
        }
    }
    else
    {
        time_printf("    solution is fractional\n");
        rval = BNC_branch_node(bnc, x, depth);
        ABORT_IF(rval, "BNC_branch_node failed\n");
    }

    CLEANUP:
    if (x) free(x);
    return rval;
}

static int BNC_branch_node(struct BNC *bnc, double *x, int depth)
{
    int rval = 0;

    struct LP *lp = bnc->lp;

    int num_cols = LP_get_num_cols(lp);
    int best_branch_var = BNC_find_best_branching_var(x, num_cols);

    time_printf("Branching on variable x%d = %.6lf (depth %d)...\n",
            best_branch_var, x[best_branch_var], depth);

    time_printf("Fixing variable x%d to one...\n", best_branch_var);
    rval = LP_change_bound(lp, best_branch_var, 'L', 1.0);
    ABORT_IF(rval, "LP_change_bound failed\n");

    rval = BNC_solve_node(bnc, depth + 1);
    ABORT_IF(rval, "BNC_solve_node failed\n");

    rval = LP_change_bound(lp, best_branch_var, 'L', 0.0);
    ABORT_IF(rval, "LP_change_bound failed\n");

    time_printf("Fixing variable x%d to zero...\n", best_branch_var);
    rval = LP_change_bound(lp, best_branch_var, 'U', 0.0);
    ABORT_IF(rval, "LP_change_bound failed\n");

    rval = BNC_solve_node(bnc, depth + 1);
    ABORT_IF(rval, "BNC_solve_node failed\n");

    rval = LP_change_bound(lp, best_branch_var, 'U', 1.0);
    ABORT_IF(rval, "LP_change_bound failed\n");

    time_printf("Finished branching on variable %d\n", best_branch_var);

    CLEANUP:
    return rval;
}

static int BNC_is_integral(double *x, int num_cols)
{
    for (int i = 0; i < num_cols; i++)
        if (x[i] > LP_EPSILON && x[i] < 1.0 - LP_EPSILON)
            return 0;

    return 1;
}

static int BNC_find_best_branching_var(double *x, int num_cols)
{
    int best_index = 0;
    double best_index_frac = 1.0;

    for (int i = 0; i < num_cols; i++)
    {
        if (fabs(x[i] - 0.5) < best_index_frac)
        {
            best_index = i;
            best_index_frac = fabs(x[i] - 0.5);
        }
    }

    return best_index;
}
