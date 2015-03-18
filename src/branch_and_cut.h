#ifndef _PROJECT_BRANCH_AND_CUT_H_
#define _PROJECT_BRANCH_AND_CUT_H_

#include "lp.h"

int bnc_solve_node(
        struct LP *lp, double *best_val, int ncount, int ecount, int *elist,
        int depth);

int bnc_branch_node(
        struct LP *lp, double *x, int ncount, int ecount, int depth,
        double *current_val, int *elist);

int bnc_init_lp(
        struct LP *lp, int node_count, int edge_count, int *edge_list, int *edge_weights);

#endif //_PROJECT_BRANCH_AND_CUT_H_
