#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>
#include "lp.h"
#include "util.h"

#define LP_EPSILON 0.000001

int lp_open(struct LP *lp)
{
    int rval = 0;

    lp->cplex_lp = (CPXLPptr)NULL;
    lp->cplex_env = CPXopenCPLEX(&rval);
    if (rval)
    {
        fprintf(stderr, "CPXopenCPLEX failed, return code %d\n", rval);
        goto CLEANUP;
    }

    CLEANUP:
    return rval;
}

void lp_free(struct LP *lp)
{
    if (!lp) return;
    if (!lp->cplex_env) return;

    if (lp->cplex_lp)
        CPXfreeprob(lp->cplex_env, &(lp->cplex_lp));

    CPXcloseCPLEX(&lp->cplex_env);
    lp->cplex_env = 0;
}

int lp_create(struct LP *lp, const char *name)
{
    int rval = 0;
    char nambuf[MAX_NAME_LENGTH];

    ABORT_IF(!lp->cplex_env, "cplex_env is null\n");

    strncpy(nambuf, name, MAX_NAME_LENGTH);
    nambuf[MAX_NAME_LENGTH - 1] = '\0';

    lp->cplex_lp = CPXcreateprob(lp->cplex_env, &rval, nambuf);
    ABORT_IF(rval, "CPXcreateprob failed\n");

    CLEANUP:
    return rval;
}

int lp_new_row(struct LP *lp, char sense, double rhs)
{
    int rval = 0;

    rval = CPXnewrows(lp->cplex_env, lp->cplex_lp, 1, &rhs, &sense,
            (double *) NULL, (char **) NULL);

    ABORT_IF(rval, "CPXnewrows failed\n");

    CLEANUP:
    return rval;
}

int lp_add_rows(
        struct LP *lp, int newrows, int newnz, double *rhs, char *sense,
        int *rmatbeg, int *rmatind, double *rmatval)
{
    int rval = 0;

    rval = CPXaddrows(lp->cplex_env, lp->cplex_lp, 0, newrows, newnz, rhs,
            sense, rmatbeg, rmatind, rmatval, (char **) NULL, (char **) NULL);

    ABORT_IF(rval, "CPXaddrows failed\n");

    CLEANUP:
    return rval;
}

int lp_add_cols(
        struct LP *lp, int newcols, int newnz, double *obj, int *cmatbeg,
        int *cmatind, double *cmatval, double *lb, double *ub)
{
    int rval = 0;

    rval = CPXaddcols(lp->cplex_env, lp->cplex_lp, newcols, newnz, obj, cmatbeg,
            cmatind, cmatval, lb, ub, (char **) NULL);

    ABORT_IF(rval, "CPXaddcols failed\n");

    CLEANUP:
    return rval;
}

int lp_change_bound(struct LP *lp, int col, char lower_or_upper, double bnd)
{
    int rval = 0;

    rval = CPXchgbds(lp->cplex_env, lp->cplex_lp, 1, &col, &lower_or_upper,
            &bnd);
    ABORT_IF(rval, "CPXchgbds failed\n");

    CLEANUP:
    return rval;
}

int lp_optimize(struct LP *lp, int *infeasible)
{
    int rval = 0, solstat;

    *infeasible = 0;

    rval = CPXdualopt(lp->cplex_env, lp->cplex_lp);
    ABORT_IF(rval, "CPXdualopt failed\n");

    solstat = CPXgetstat(lp->cplex_env, lp->cplex_lp);
    if (solstat == CPX_STAT_INFEASIBLE)
    {
        if (infeasible)
            *infeasible = 1;

    }
    else
    {
        ABORT_IF(solstat != CPX_STAT_OPTIMAL
                && solstat != CPX_STAT_OPTIMAL_INFEAS,
                "Invalid solution status");
    }

    CLEANUP:
    return rval;
}

int lp_get_obj_val(struct LP *lp, double *obj)
{
    int rval = 0;

    rval = CPXgetobjval(lp->cplex_env, lp->cplex_lp, obj);
    ABORT_IF(rval, "CPXgetobjval failed\n");

    CLEANUP:
    return rval;
}

int lp_get_x(struct LP *lp, double *x)
{
    int rval = 0;

    int ncols = CPXgetnumcols(lp->cplex_env, lp->cplex_lp);
    ABORT_IF(!ncols, "No columns in LP\n");

    rval = CPXgetx(lp->cplex_env, lp->cplex_lp, x, 0, ncols - 1);
    ABORT_IF(rval, "CPXgetx failed\n");

    CLEANUP:
    return rval;
}

int lp_write(struct LP *lp, const char *fname)
{
    int rval = 0;

    char nambuf[MAX_NAME_LENGTH];
    strncpy(nambuf, fname, MAX_NAME_LENGTH);
    nambuf[MAX_NAME_LENGTH - 1] = '\0';

    rval = CPXwriteprob(lp->cplex_env, lp->cplex_lp, nambuf, "RLP");
    ABORT_IF(rval, "CPXwriteprob failed\n");

    CLEANUP:
    return rval;
}
