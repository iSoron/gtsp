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

#include "gtsp-cols.h"
#include "util.h"

int GTSP_find_columns(struct LP *lp, struct GTSP *data)
{
    int rval = 0;

    double *y = 0;
    double initial_time = get_user_time();
    int num_rows = LP_get_num_rows(lp);

    y = (double *) malloc(num_rows * sizeof(double));
    abort_if(!y, "could not allocate y");

    rval = LP_get_y(lp, y);
    abort_if(rval, "LP_get_y failed");

    int violated_count = 0;
    int edge_count = data->graph->edge_count;

    log_debug("Finding new columns...\n");

    for (int i = 0; i < edge_count; i++)
    {
        struct Edge *e = &data->graph->edges[i];
        if (e->column >= 0) continue;

        double sum = 0;

        for (int j = 0; j < lp->cut_pool_size; j++)
        {
            struct Row *cut = lp->cut_pool[j];
            if (!cut->edges[e->index]) continue;

            sum += y[cut->cplex_row_index];
        }

        if (sum - EPSILON < e->weight) continue;

        rval = GTSP_add_column(lp, e);
        abort_if(rval, "GTSP_add_column failed");

        violated_count++;
    }

    log_debug("    %d of %d edges are violated\n", violated_count, edge_count);

    COLUMNS_TIME += get_user_time() - initial_time;

    CLEANUP:
    if (y) free(y);
    return rval;
}

int GTSP_add_column(struct LP *lp, struct Edge *e)
{
    int rval = 0;

    int *cmatind = 0;
    double *cmatval = 0;

    int num_rows = LP_get_num_rows(lp);
    int num_cols = LP_get_num_cols(lp);

    cmatind = (int *) malloc(num_rows * sizeof(int));
    cmatval = (double *) malloc(num_rows * sizeof(double));

    abort_if(!cmatind, "could not allocate cmatbeg");
    abort_if(!cmatval, "could not allocate cmatval");

    int nz = 0;
    int cmatbeg = 0;
    double lb = 0.0;
    double ub = 1.0;
    double obj = e->weight;

    for (int j = 0; j < lp->cut_pool_size; j++)
    {
        struct Row *cut = lp->cut_pool[j];
        if (!cut->edges[e->index]) continue;

        cmatind[nz] = cut->cplex_row_index;
        cmatval[nz] = cut->edge_val;
        nz++;
    }

    e->column = num_cols;
    LP_add_cols(lp, 1, nz, &obj, &cmatbeg, cmatind, cmatval, &lb, &ub);

    CLEANUP:
    if (cmatind) free(cmatind);
    if (cmatval) free(cmatval);
    return rval;
}

