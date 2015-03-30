#ifndef _PROJECT_BRANCH_AND_CUT_H_
#define _PROJECT_BRANCH_AND_CUT_H_

#include "lp.h"

struct TOUR {
	int vertex;
	int next;
	int prev;
	};

struct BNC
{
    struct LP *lp;

    double *best_x;
    double best_obj_val;

    double *optimal_x;

    int *problem_data;

    int (*problem_init_lp)(struct LP *, void *);

    int (*problem_add_cutting_planes)(struct LP *, void *);

    int (*problem_solution_found)(void *data, double *x);
};

int BNC_init(struct BNC *bnc);

int BNC_solve(struct BNC *bnc);

int BNC_init_lp(struct BNC *bnc);

void BNC_free(struct BNC *bnc);

int re_optimize_integral(struct BNC *bnc);

//int optimize_vertex_in_cluster(struct BNC *bnc, double best_val);

extern int BNC_NODE_COUNT;


#endif //_PROJECT_BRANCH_AND_CUT_H_
