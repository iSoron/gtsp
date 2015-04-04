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

#ifndef _PROJECT_MAIN_H_
#define _PROJECT_MAIN_H_

extern unsigned int SEED;

extern double *OPTIMAL_X;
extern double SUBTOUR_TIME;
extern double COMBS_TIME;

extern double LP_SOLVE_TIME;
extern int LP_SOLVE_COUNT;
extern int LP_MAX_ROWS;
extern int LP_MAX_COLS;

extern double CUT_POOL_TIME;
extern long CUT_POOL_MAX_MEMORY;

extern int FLOW_MAX_FLOW_COUNT;

extern double TOTAL_TIME;
extern double INITIAL_TIME;
extern double ROOT_VALUE;

extern int SUBTOUR_COUNT;
extern int COMBS_COUNT;

extern char LP_FILENAME[100];
extern char SOLUTION_FILENAME[100];
extern char FRAC_SOLUTION_FILENAME[100];

#endif