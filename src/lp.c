#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>
#include "lp.h"
#include "util.h"

#define LP_EPSILON 0.000001
#define MAX_CUT_POOL_SIZE 100000

int LP_open(struct LP *lp)
{
    int rval = 0;

    lp->cplex_lp = (CPXLPptr) NULL;
    lp->cplex_env = CPXopenCPLEX(&rval);
    abort_if(rval, "CPXopenCPLEX failed");

    lp->cut_pool = (struct Row **) malloc(
            MAX_CUT_POOL_SIZE * sizeof(struct Row *));
    abort_if(!lp->cut_pool, "could not allocate cut_pool");

    CLEANUP:
    return rval;
}

void LP_free(struct LP *lp)
{
    if (!lp) return;
    if (!lp->cplex_env) return;

    if (lp->cplex_lp)
        CPXfreeprob(lp->cplex_env, &(lp->cplex_lp));

    CPXcloseCPLEX(&lp->cplex_env);
    lp->cplex_env = 0;
    if (lp->cut_pool) free(lp->cut_pool);
}

int LP_create(struct LP *lp, const char *name)
{
    int rval = 0;
    char nambuf[MAX_NAME_LENGTH];

    abort_if(!lp->cplex_env, "cplex_env is null");

    strncpy(nambuf, name, MAX_NAME_LENGTH);
    nambuf[MAX_NAME_LENGTH - 1] = '\0';

    lp->cplex_lp = CPXcreateprob(lp->cplex_env, &rval, nambuf);
    abort_if(rval, "CPXcreateprob failed");

    CLEANUP:
    return rval;
}

int LP_new_row(struct LP *lp, char sense, double rhs)
{
    int rval = 0;

    rval = CPXnewrows(lp->cplex_env, lp->cplex_lp, 1, &rhs, &sense, 0, 0);

    abort_if(rval, "CPXnewrows failed");

    CLEANUP:
    return rval;
}

int LP_add_rows(
        struct LP *lp,
        int newrows,
        int newnz,
        double *rhs,
        char *sense,
        int *rmatbeg,
        int *rmatind,
        double *rmatval)
{
    int rval = 0;

    rval = CPXaddrows(lp->cplex_env, lp->cplex_lp, 0, newrows, newnz, rhs,
            sense, rmatbeg, rmatind, rmatval, 0, 0);

    abort_if(rval, "CPXaddrows failed");

    CLEANUP:
    return rval;
}

int LP_add_cols(
        struct LP *lp,
        int newcols,
        int newnz,
        double *obj,
        int *cmatbeg,
        int *cmatind,
        double *cmatval,
        double *lb,
        double *ub)
{
    int rval = 0;

    rval = CPXaddcols(lp->cplex_env, lp->cplex_lp, newcols, newnz, obj, cmatbeg,
            cmatind, cmatval, lb, ub, (char **) NULL);
    abort_if(rval, "CPXaddcols failed");

    CLEANUP:
    return rval;
}

int LP_change_bound(struct LP *lp, int col, char lower_or_upper, double bnd)
{
    int rval = 0;

    rval = CPXchgbds(lp->cplex_env, lp->cplex_lp, 1, &col, &lower_or_upper,
            &bnd);
    abort_if(rval, "CPXchgbds failed");

    CLEANUP:
    return rval;
}

extern int LP_OPTIMIZE_COUNT;
extern double LP_CPU_TIME;

int LP_optimize(struct LP *lp, int *infeasible)
{
    LP_OPTIMIZE_COUNT++;

    int rval = 0, solstat;

    *infeasible = 0;

    double current = get_current_time();
    rval = CPXdualopt(lp->cplex_env, lp->cplex_lp);
    abort_if(rval, "CPXdualopt failed");
    LP_CPU_TIME += get_current_time() - current;

    solstat = CPXgetstat(lp->cplex_env, lp->cplex_lp);
    if (solstat == CPX_STAT_INFEASIBLE)
    {
        *infeasible = 1;
        goto CLEANUP;
    }
    else
    {
        abort_if(solstat != CPX_STAT_OPTIMAL && solstat != CPXMIP_OPTIMAL &&
                solstat != CPX_STAT_OPTIMAL_INFEAS, "Invalid solution status");
    }

    CLEANUP:
    return rval;
}

int LP_remove_slacks(struct LP *lp, int first_row, double max_slack)
{
    int rval = 0;

    double *slacks = 0;
    int *should_remove = 0;

    int numrows = CPXgetnumrows(lp->cplex_env, lp->cplex_lp);
    if (numrows < 5000) return 0;

    should_remove = (int *) malloc((numrows + 1) * sizeof(int));
    abort_if(!should_remove, "could not allocate should_remove");

    slacks = (double *) malloc(numrows * sizeof(double));
    abort_if(!slacks, "could not allocate slacks");

    rval = CPXgetslack(lp->cplex_env, lp->cplex_lp, slacks, 0, numrows - 1);
    abort_if(rval, "CPXgetslack failed");

    for (int i = 0; i < numrows; i++)
        should_remove[i] = (slacks[i] < -max_slack);
    should_remove[numrows] = 0;

    log_debug("Deleting constraints...\n");
    int start = 0;
    int end = -1;
    int count = 0;
    for (int i = first_row; i < numrows; i++)
    {
        if (should_remove[i])
        {
            end++;
        }
        else
        {
            if (end >= start)
            {
                rval = CPXdelrows(lp->cplex_env, lp->cplex_lp, start - count,
                        end - count);
                abort_if(rval, "CPXdelrows failed");
                log_verbose("    %d %d (%d)\n", start, end, end - start + 1);

                count += end - start + 1;
            }

            start = i + 1;
            end = i;
        }
    }

    log_info("Removed %d of %d constraints\n", count, numrows);

    CLEANUP:
    if (should_remove) free(should_remove);
    if (slacks) free(slacks);
    return rval;
}

int LP_get_obj_val(struct LP *lp, double *obj)
{
    int rval = 0;

    rval = CPXgetobjval(lp->cplex_env, lp->cplex_lp, obj);
    abort_if(rval, "CPXgetobjval failed");

    CLEANUP:
    return rval;
}

int LP_get_x(struct LP *lp, double *x)
{
    int rval = 0;

    int ncols = CPXgetnumcols(lp->cplex_env, lp->cplex_lp);
    abort_if(!ncols, "No columns in LP");

    rval = CPXgetx(lp->cplex_env, lp->cplex_lp, x, 0, ncols - 1);
    abort_if(rval, "CPXgetx failed");

    CLEANUP:
    return rval;
}

int LP_get_num_cols(struct LP *lp)
{
    return CPXgetnumcols(lp->cplex_env, lp->cplex_lp);
}

int LP_write(struct LP *lp, const char *fname)
{
    int rval = 0;

    char nambuf[MAX_NAME_LENGTH];
    strncpy(nambuf, fname, MAX_NAME_LENGTH);
    nambuf[MAX_NAME_LENGTH - 1] = '\0';

    rval = CPXwriteprob(lp->cplex_env, lp->cplex_lp, nambuf, "RLP");
    abort_if(rval, "CPXwriteprob failed");

    CLEANUP:
    return rval;
}

#define return_if_neq(a, b) \
    if((a)<(b)) return -1; \
    if((a)>(b)) return 1;

#define return_if_neq_epsilon(a, b) \
    if((a+LP_EPSILON)<(b)) return -1; \
    if((a-LP_EPSILON)>(b)) return 1;

int compare_cuts(struct Row *cut1, struct Row *cut2)
{
    return_if_neq(cut1->nz, cut2->nz);

    for (int i = 0; i < cut1->nz; i++)
    {
        return_if_neq(cut1->rmatind[i], cut2->rmatind[i]);
        return_if_neq_epsilon(cut1->rmatval[i], cut2->rmatval[i]);
    }

    return 0;
}

int LP_add_cut(struct LP *lp, struct Row *cut)
{
    int rval = 0;

    rval = LP_update_hash(cut);
    abort_if(rval, "LP_update_hash failed");

    for (int i = 0; i < lp->cut_pool_size; i++)
    {
        if (lp->cut_pool[i]->hash != cut->hash) continue;
        if (!compare_cuts(lp->cut_pool[i], cut))
        {
            free(cut->rmatval);
            free(cut->rmatind);
            free(cut);
            return 0;
        }
    }

    lp->cut_pool[lp->cut_pool_size++] = cut;

    int rmatbeg = 0;
    rval = LP_add_rows(lp, 1, cut->nz, &cut->rhs, &cut->sense, &rmatbeg,
            cut->rmatind, cut->rmatval);
    abort_if(rval, "LP_add_rows failed");

    CLEANUP:
    return rval;
}

int LP_update_hash(struct Row *cut)
{
    unsigned long hash = 0;

    for (int i = 0; i < cut->nz; i++)
    {
        hash += cut->rmatind[i] * 65521;
        hash %= 4294967291;
    }

    cut->hash = hash;

    return 0;
}