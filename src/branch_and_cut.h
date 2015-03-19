#ifndef _PROJECT_BRANCH_AND_CUT_H_
#define _PROJECT_BRANCH_AND_CUT_H_

#include "lp.h"

struct BNC
{
    struct LP *lp;

    double *best_x;
    double best_obj_val;

    int *problem_data;

    int (*problem_init_lp)(struct LP *, void *);

    int (*problem_add_cutting_planes)(struct LP *, void *);
};

int BNC_init(struct BNC *bnc);

int BNC_solve(struct BNC *bnc);

int BNC_init_lp(struct BNC *bnc);

void BNC_free(struct BNC *bnc);

#endif //_PROJECT_BRANCH_AND_CUT_H_
