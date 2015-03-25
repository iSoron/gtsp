#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lp.h"
#include "branch_and_cut.h"
#include "util.h"

int BNC_NODE_COUNT = 0;

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
    bnc->problem_solution_found = 0;

    bnc->best_x = 0;
    bnc->best_obj_val = 0;

    bnc->lp = (struct LP *) malloc(sizeof(struct LP));
    abort_if(!bnc->lp, "could not allocate bnc->lp");

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
    if (bnc->best_x) free(bnc->best_x);
}

int BNC_init_lp(struct BNC *bnc)
{
    int rval = 0;

    rval = LP_open(bnc->lp);
    abort_if(rval, "LP_open failed");

    rval = LP_create(bnc->lp, "subtour");
    abort_if(rval, "LP_create failed");

    rval = bnc->problem_init_lp(bnc->lp, bnc->problem_data);
    abort_if(rval, "problem_init_lp failed");

    rval = LP_write(bnc->lp, "subtour.lp");
    abort_if(rval, "LP_write failed");

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

    BNC_NODE_COUNT++;

    int rval = 0;
    double *x = (double *) NULL;

    log_debug("Optimizing...\n");

    int is_infeasible;
    rval = LP_optimize(lp, &is_infeasible);
    abort_if(rval, "LP_optimize failed\n");

    if (is_infeasible)
    {
        log_debug("Branch pruned by infeasibility.\n");
        goto CLEANUP;
    }

    double objval = 0;
    rval = LP_get_obj_val(lp, &objval);
    abort_if(rval, "LP_get_obj_val failed\n");

    log_debug("    obj value = %.2f\n", objval);

    if (objval > *best_val + LP_EPSILON)
    {
        log_debug("Branch pruned by bound (%.2lf > %.2lf).\n", objval,
                *best_val);
        rval = 0;
        goto CLEANUP;
    }

    int num_cols = LP_get_num_cols(lp);

    x = (double *) malloc(num_cols * sizeof(double));
    abort_if(!x, "could not allocate x");

    rval = LP_get_x(lp, x);
    abort_if(rval, "LP_get_x failed");

    if (bnc->problem_add_cutting_planes)
    {
//        log_info("Adding problem cutting planes...\n");
        rval = bnc->problem_add_cutting_planes(lp, bnc->problem_data);
        abort_if(rval, "problem_add_cutting_planes failed");
    }

    rval = LP_optimize(lp, &is_infeasible);
    abort_if(rval, "LP_optimize failed\n");

    rval = LP_get_obj_val(lp, &objval);
    abort_if(rval, "LP_get_obj_val failed");

    rval = LP_get_x(lp, x);
    abort_if(rval, "LP_get_x failed");

    if (BNC_is_integral(x, num_cols))
    {
        log_debug("    solution is integral\n");

        if (objval + LP_EPSILON < *best_val)
        {
            if (bnc->best_x) free(bnc->best_x);
            *best_val = objval;
            bnc->best_x = x;
            x = 0;

            log_info("Found a better integral solution:\n");
            log_info("    obj val = %.2lf **\n", objval);

            if (bnc->problem_solution_found)
                bnc->problem_solution_found(bnc->problem_data, bnc->best_x);
        }
    } else
    {
        log_debug("    solution is fractional\n");
        rval = BNC_branch_node(bnc, x, depth);
        abort_if(rval, "BNC_branch_node failed");
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

    log_debug("Branching on variable x%d = %.6lf (depth %d)...\n",
            best_branch_var, x[best_branch_var], depth);

    log_debug("Fixing variable x%d to one...\n", best_branch_var);
    rval = LP_change_bound(lp, best_branch_var, 'L', 1.0);
    abort_if(rval, "LP_change_bound failed");

    rval = BNC_solve_node(bnc, depth + 1);
    abort_if(rval, "BNC_solve_node failed");

    rval = LP_change_bound(lp, best_branch_var, 'L', 0.0);
    abort_if(rval, "LP_change_bound failed");

    log_debug("Fixing variable x%d to zero...\n", best_branch_var);
    rval = LP_change_bound(lp, best_branch_var, 'U', 0.0);
    abort_if(rval, "LP_change_bound failed");

    rval = BNC_solve_node(bnc, depth + 1);
    abort_if(rval, "BNC_solve_node failed");

    rval = LP_change_bound(lp, best_branch_var, 'U', 1.0);
    abort_if(rval, "LP_change_bound failed");

    log_debug("Finished branching on variable %d\n", best_branch_var);

    CLEANUP:
    return rval;
}

static int BNC_is_integral(double *x, int num_cols)
{
//    return 1;

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
