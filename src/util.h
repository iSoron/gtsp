#ifndef _PROJECT_UTIL_H_
#define _PROJECT_UTIL_H_

#define ABORT_IF(cond, msg) if(cond) { \
    fprintf(stderr, msg); rval = 1; goto CLEANUP; }

double get_current_time(void);

double get_real_time();

void time_printf(const char *fmt, ...);

void next_set(int sz, int *set);

#endif
