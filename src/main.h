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