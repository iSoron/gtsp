#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>
#include "lp.h"
#include "util.h"

#define LP_EPSILON 0.000001

int LP_open(struct LP *lp)
{
    int rval = 0;

    lp->cplex_lp = (CPXLPptr) NULL;
    lp->cplex_env = CPXopenCPLEX(&rval);
    if (rval)
    {
        fprintf(stderr, "CPXopenCPLEX failed, return code %d\n", rval);
        goto CLEANUP;
    }

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

int LP_optimize(struct LP *lp, int *infeasible)
{
    int rval = 0, solstat;

    *infeasible = 0;

    rval = CPXdualopt(lp->cplex_env, lp->cplex_lp);
    abort_if(rval, "CPXdualopt failed");

    solstat = CPXgetstat(lp->cplex_env, lp->cplex_lp);
    if (solstat == CPX_STAT_INFEASIBLE)
    {
        *infeasible = 1;
    }
    else
    {
        abort_if(solstat != CPX_STAT_OPTIMAL
                && solstat != CPX_STAT_OPTIMAL_INFEAS,
                "Invalid solution status");
    }

    CLEANUP:
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
