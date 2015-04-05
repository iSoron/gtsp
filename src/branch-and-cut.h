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

#ifndef _PROJECT_BRANCH_AND_CUT_H_
#define _PROJECT_BRANCH_AND_CUT_H_

#include "lp.h"

struct BNC
{

    struct LP *lp;

    /*
     * Best integral solution found so far.
     */
    double *best_x;

    /*
     * Value of the best integral solution found so far.
     */
    double best_obj_val;

    /*
     * Pointer to problem-specific data. This pointer is passed to other
     * problem-specific functions, such as problem_init_lp and
     * problem_add_cutting_planes.
     */
    void *problem_data;

    /*
     * This callback is called at the beginning of the branch-and-cut search.
     * It should initialize the LP, create the initial rows and columns.
     */
    int (*problem_init_lp)(struct LP *lp, void *data);

    /*
     * This callback is called after an LP relaxation is solved and before
     * checking whether the solution is integral. It should find and add
     * any problem-specific cutting planes to the LP.
     */
    int (*problem_add_cutting_planes)(struct LP *lp, void *data);

    /*
     * This callback is called after a better integral solution has been found.
     * This solution satisfies all the problem-specific cuts.
     */
    int (*problem_solution_found)(struct BNC *, void *data, double *x);
};

int BNC_init(struct BNC *bnc);

int BNC_solve(struct BNC *bnc);

int BNC_init_lp(struct BNC *bnc);

void BNC_free(struct BNC *bnc);

#endif //_PROJECT_BRANCH_AND_CUT_H_
