#ifndef _PROJECT_UTIL_H_
#define _PROJECT_UTIL_H_

#define ABORT_IF(cond, msg) if(cond) { \
    fprintf(stderr, msg); rval = 1; goto CLEANUP; }

double get_current_time(void);

int build_random_2d_points
        (int node_count, double *x_list, double *y_list, int grid_size);

double get_current_time(void);

void time_printf(const char *fmt, ...);

void next_set(int sz, int *set);

#endif
