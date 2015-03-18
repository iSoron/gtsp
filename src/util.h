#ifndef __CO759_UTIL_H
#define __CO759_UTIL_H

#define ABORT_IF(cond, msg) if(cond) { \
    fprintf(stderr, msg); rval = 1; goto CLEANUP; }

double util_get_current_time(void);

int CO759_build_xy(int ncount, double *xlist, double *ylist, int gridsize);

double util_get_current_time(void);

void time_printf(const char *fmt, ...);

void next_set(int sz, int *Set);

#endif  /* __CO759_UTIL_H */
