#ifndef __MAIN_
#define __MAIN_

#include <cplex.h>

#define LP_EPSILON 0.000001

struct LP
{
    CPXENVptr cplex_env;
    CPXLPptr cplex_lp;
};

static const int MAX_NAME_LENGTH = 100;

int lp_open(struct LP *lp);

int lp_create(struct LP *lp, const char *name);

int lp_write(struct LP *lp, const char *fname);

void lp_free(struct LP *lp);

int lp_new_row(struct LP *lp, char sense, double rhs);

int lp_add_rows(struct LP *lp, int newrows, int newnz, double *rhs, char *sense,
        int *rmatbeg, int *rmatind, double *rmatval);

int lp_add_cols(struct LP *lp, int newcols, int newnz, double *obj,
        int *cmatbeg, int *cmatind, double *cmatval, double *lb, double *ub);

int lp_change_bound(struct LP *lp, int col, char lower_or_upper, double bnd);

int lp_optimize(struct LP *lp, int *infeasible);

int lp_get_obj_val(struct LP *lp, double *obj);

int lp_get_x(struct LP *lp, double *x);

#endif