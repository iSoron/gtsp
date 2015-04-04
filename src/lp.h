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

#ifndef _PROJECT_LP_H_
#define _PROJECT_LP_H_

#include <cplex.h>
#include "params.h"

struct LP
{
    CPXENVptr cplex_env;
    CPXLPptr cplex_lp;

    int cut_pool_size;
    struct Row **cut_pool;
};

struct Row
{
    unsigned long hash;
    int cplex_row_index;
    int age;

    int nz;
    char sense;
    double rhs;
    double *rmatval;
    int *rmatind;
};

static const int MAX_NAME_LENGTH = 100;

int LP_open(struct LP *lp);

int LP_create(struct LP *lp, const char *name);

int LP_write(struct LP *lp, const char *fname);

void LP_free(struct LP *lp);

int LP_new_row(struct LP *lp, char sense, double rhs);

int LP_add_rows(
        struct LP *lp,
        int newrows,
        int newnz,
        double *rhs,
        char *sense,
        int *rmatbeg,
        int *rmatind,
        double *rmatval);

int LP_add_cols(
        struct LP *lp,
        int newcols,
        int newnz,
        double *obj,
        int *cmatbeg,
        int *cmatind,
        double *cmatval,
        double *lb,
        double *ub);

int LP_change_bound(struct LP *lp, int col, char lower_or_upper, double bnd);

int LP_optimize(struct LP *lp, int *infeasible);

int LP_get_obj_val(struct LP *lp, double *obj);

int LP_get_x(struct LP *lp, double *x);

int LP_get_num_cols(struct LP *lp);

int LP_add_cut(struct LP *lp, struct Row *cut);

#endif