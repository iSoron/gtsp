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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>
#include <assert.h>
#include "lp.h"
#include "util.h"
#include "gtsp-subtour.h"

static int compress_cut_pool(struct LP *lp)
{
    int delete_count = 0;
    for (int i = 0; i < lp->cut_pool_size; i++)
    {
        struct Row *cut = lp->cut_pool[i];

        if (cut->cplex_row_index < 0)
        {
            free(cut->edges);
            free(cut);
            delete_count++;
        }
        else
        {
            lp->cut_pool[i - delete_count] = lp->cut_pool[i];
        }
    }

    lp->cut_pool_size -= delete_count;

    return 0;
}

static int remove_old_cuts(struct LP *lp)
{
    int rval = 0;

    int *should_remove = 0;

    int numrows = CPXgetnumrows(lp->cplex_env, lp->cplex_lp);
    log_verbose("    numrows=%d\n", numrows);

    should_remove = (int *) malloc((numrows + 1) * sizeof(int));
    abort_if(!should_remove, "could not allocate should_remove");

    for (int i = 0; i < numrows; i++)
        should_remove[i] = 0;
    should_remove[numrows] = 0;

    log_verbose("Old cplex row index:\n");
    for (int i = 0; i < lp->cut_pool_size; i++)
            log_verbose("    %d\n", lp->cut_pool[i]->cplex_row_index);

    log_verbose("Should remove:\n");
    for (int i = 0; i < lp->cut_pool_size; i++)
    {
        struct Row *cut = lp->cut_pool[i];
        if (cut->age <= MAX_CUT_AGE) continue;
        if (cut->cplex_row_index < 0) continue;

        should_remove[cut->cplex_row_index] = 1;
        log_verbose("    %d\n", cut->cplex_row_index);
    }

    // Update cut->cplex_row_index
    int count = 0;
    for (int i = 0; i < numrows; i++)
    {
        if (!should_remove[i]) continue;

        for (int j = 0; j < lp->cut_pool_size; j++)
        {
            struct Row *cut = lp->cut_pool[j];

            if (cut->cplex_row_index == i - count)
                cut->cplex_row_index = -1;
            else if (cut->cplex_row_index > i - count)
                cut->cplex_row_index--;
        }

        count++;
    }

    log_verbose("New cplex row index:\n");
    for (int i = 0; i < lp->cut_pool_size; i++)
            log_verbose("    %d\n", lp->cut_pool[i]->cplex_row_index);

    // Remove from CPLEX, finding the largest intervals of contiguous rows
    // that should be removed.
    int start = 0;
    int end = -1;
    count = 0;
    for (int i = 0; i < numrows + 1; i++)
    {
        if (should_remove[i])
            end++;
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

    log_debug("    removed %d old cuts\n", count);

    rval = CPXdualopt(lp->cplex_env, lp->cplex_lp);
    abort_if(rval, "CPXoptimize failed");

    log_debug("Compressing cut pool...\n");
    compress_cut_pool(lp);

    long nz = 0;
    long size = 0;
    for (int i = 0; i < lp->cut_pool_size; i++)
    {
        size += sizeof(struct Row);
        struct Row *cut = lp->cut_pool[i];

        nz += cut->nz;
        size += cut->nz * sizeof(double);
        size += cut->nz * sizeof(int);
        size += cut->edge_count * sizeof(char);
    }

    log_debug("    %ld cuts (%ld nz, %ld MiB)\n", lp->cut_pool_size, nz,
            size / 1024 / 1024);

    if (size > CUT_POOL_MAX_MEMORY)
        CUT_POOL_MAX_MEMORY = size;

    CLEANUP:
    if (should_remove) free(should_remove);
    return rval;
}

static int update_cut_ages(struct LP *lp)
{
    int rval = 0;

    double *slacks = 0;

    int numrows = CPXgetnumrows(lp->cplex_env, lp->cplex_lp);

    slacks = (double *) malloc(numrows * sizeof(double));
    abort_if(!slacks, "could not allocate slacks");

    rval = CPXgetslack(lp->cplex_env, lp->cplex_lp, slacks, 0, numrows - 1);
    abort_if(rval, "CPXgetslack failed");

    for (int i = 0; i < lp->cut_pool_size; i++)
    {
        struct Row *cut = lp->cut_pool[i];
        if (cut->cplex_row_index < 0) continue;
        assert(cut->cplex_row_index < numrows);

        if (slacks[cut->cplex_row_index] < -EPSILON)
            cut->age++;
        else
            cut->age = 0;
    }

    CLEANUP:
    if (slacks) free(slacks);
    return rval;
}

int LP_open(struct LP *lp)
{
    int rval = 0;

    lp->cplex_lp = (CPXLPptr) NULL;
    lp->cplex_env = CPXopenCPLEX(&rval);
    abort_if(rval, "CPXopenCPLEX failed");

    lp->cut_pool_size = 0;
    lp->cut_pool = (struct Row **) malloc(
            MAX_CUT_POOL_SIZE * sizeof(struct Row *));
    abort_if(!lp->cut_pool, "could not allocate cut_pool");

    CLEANUP:
    return rval;
}

void LP_free_cut_pool(struct LP *lp)
{
    for (int i = 0; i < lp->cut_pool_size; i++)
    {
        struct Row *cut = lp->cut_pool[i];
        free(cut->edges);
        free(cut);
    }

    if (lp->cut_pool) free(lp->cut_pool);
}

void LP_free(struct LP *lp)
{
    if (!lp) return;
    if (!lp->cplex_env) return;

    if (lp->cplex_lp)
        CPXfreeprob(lp->cplex_env, &(lp->cplex_lp));

    CPXcloseCPLEX(&lp->cplex_env);
    lp->cplex_env = 0;

    LP_free_cut_pool(lp);
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
    LP_SOLVE_COUNT++;

    int rval = 0, solstat;

    *infeasible = 0;

    int numrows = CPXgetnumrows(lp->cplex_env, lp->cplex_lp);
    int numcols = CPXgetnumcols(lp->cplex_env, lp->cplex_lp);

    if (numrows > LP_MAX_ROWS) LP_MAX_ROWS = numrows;
    if (numrows > LP_MAX_COLS) LP_MAX_COLS = numcols;

    log_debug("Optimizing LP (%d rows %d cols)...\n", numrows, numcols);

    double initial_time = get_user_time();

    rval = CPXdualopt(lp->cplex_env, lp->cplex_lp);
    abort_if(rval, "CPXdualopt failed");

    LP_SOLVE_TIME += get_user_time() - initial_time;

    solstat = CPXgetstat(lp->cplex_env, lp->cplex_lp);
    if (solstat == CPX_STAT_INFEASIBLE)
    {
        log_debug("    infeasible\n");

        *infeasible = 1;
        goto CLEANUP;
    }
    else
    {
        abort_if(solstat != CPX_STAT_OPTIMAL && solstat != CPXMIP_OPTIMAL &&
                solstat != CPX_STAT_OPTIMAL_INFEAS, "Invalid solution status");
    }

    double objval;
    rval = LP_get_obj_val(lp, &objval);
    abort_if(rval, "LP_get_obj_val failed");

    log_debug("    obj val = %.4lf\n", objval);
    log_debug("    time = %.4lf\n", get_user_time() - initial_time);

    initial_time = get_user_time();
    rval = update_cut_ages(lp);
    abort_if(rval, "update_cut_ages failed");

    if (LP_SOLVE_COUNT % MAX_CUT_AGE == 0)
    {
        log_debug("Removing old cuts...\n");
        rval = remove_old_cuts(lp);
        abort_if(rval, "LP_remove_old_cuts failed");
    }

    CUT_POOL_TIME += get_user_time() - initial_time;

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

int LP_get_y(struct LP *lp, double *y)
{
    int rval = 0;

    int nrows = CPXgetnumrows(lp->cplex_env, lp->cplex_lp);
    abort_if(!nrows, "No rows in LP");

    rval = CPXgetpi(lp->cplex_env, lp->cplex_lp, y, 0, nrows - 1);
    abort_iff(rval, "CPXgetpi failed (errno = %d)", rval);

    CLEANUP:
    return rval;
}

int LP_get_num_cols(struct LP *lp)
{
    return CPXgetnumcols(lp->cplex_env, lp->cplex_lp);
}

int LP_get_num_rows(struct LP *lp)
{
    return CPXgetnumrows(lp->cplex_env, lp->cplex_lp);
}

int LP_write(struct LP *lp, const char *fname)
{
    int rval = 0;
    char nambuf[MAX_NAME_LENGTH];

    FILE *f = fopen(fname, "w");
    abort_iff(!f, "could not open file %s", fname);
    fclose(f);

    strncpy(nambuf, fname, MAX_NAME_LENGTH);
    nambuf[MAX_NAME_LENGTH - 1] = '\0';

    log_info("Writing LP to file %s...\n", fname);
    rval = CPXwriteprob(lp->cplex_env, lp->cplex_lp, nambuf, "RLP");
    abort_if(rval, "CPXwriteprob failed");

    CLEANUP:
    return rval;
}

#define return_if_neq(a, b) \
    if((a)<(b)) return -1; \
    if((a)>(b)) return 1;

#define return_if_neq_epsilon(a, b) \
    if((a+EPSILON)<(b)) return -1; \
    if((a-EPSILON)>(b)) return 1;

int compare_cuts(struct Row *cut1, struct Row *cut2)
{

    for (int i = 0; i < cut1->edge_count; i++)
    {
        return_if_neq(cut1->edges[i], cut2->edges[i]);
    }

    return 0;
}

static int update_hash(struct Row *cut)
{
    unsigned long hash = 0;

    for (int i = 0; i < cut->edge_count; i++)
    {
        hash += cut->edges[i] * 65521;
        hash %= 4294967291;
    }

    cut->hash = hash;

    return 0;
}

int LP_add_cut(struct LP *lp, struct Row *cut)
{
    int rval = 0;
    double initial_time = get_user_time();

    rval = update_hash(cut);
    abort_if(rval, "LP_update_hash failed");

    for (int i = 0; i < lp->cut_pool_size; i++)
    {
        if (lp->cut_pool[i]->cplex_row_index < 0) continue;
        if (lp->cut_pool[i]->hash != cut->hash) continue;
        if (!compare_cuts(lp->cut_pool[i], cut))
        {
            log_verbose("Discarding duplicate cut (same as cplex row %d)\n",
                    lp->cut_pool[i]->cplex_row_index);

            free(cut->rmatval);
            free(cut->rmatind);
            free(cut->edges);
            free(cut);
            return 0;
        }
    }

    abort_if(lp->cut_pool_size > MAX_CUT_POOL_SIZE, "Cut pool is too large");
    lp->cut_pool[lp->cut_pool_size++] = cut;

    int rmatbeg = 0;
    rval = LP_add_rows(lp, 1, cut->nz, &cut->rhs, &cut->sense, &rmatbeg,
            cut->rmatind, cut->rmatval);
    abort_if(rval, "LP_add_rows failed");

    cut->cplex_row_index = CPXgetnumrows(lp->cplex_env, lp->cplex_lp) - 1;
    cut->age = 0;

    free(cut->rmatval);
    free(cut->rmatind);
    cut->rmatind = 0;
    cut->rmatval = 0;

    CUT_POOL_TIME += get_user_time() - initial_time;
    CLEANUP:
    return rval;
}
