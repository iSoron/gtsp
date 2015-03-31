#ifndef _PROJECT_MAIN_H_
#define _PROJECT_MAIN_H_

extern char *INPUT_FILENAME;
extern unsigned int SEED;
extern int GEOMETRIC_DATA;
extern int NODE_COUNT_RAND;
extern int GRID_SIZE_RAND;

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
extern double ROOT_VALUE;

extern int SUBTOUR_CLUSTER_CLUSTER_COUNT;
extern int SUBTOUR_NODE_CLUSTER_COUNT;
extern int SUBTOUR_NODE_NODE_COUNT;

extern int COMBS_COUNT;

#endif