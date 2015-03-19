#include <stdio.h>
#include <sys/resource.h>
#include <stdarg.h>
#include <time.h>
#include "util.h"

double get_current_time()
{
    struct rusage ru;

    getrusage(RUSAGE_SELF, &ru);

    return ((double) ru.ru_utime.tv_sec)
            + ((double) ru.ru_utime.tv_usec) / 1000000.0;
}

double get_real_time()
{
    return (double) time (0);
}

static double initial_time = 0;

void time_printf(const char *fmt, ...)
{
    if (initial_time == 0)
        initial_time = get_current_time();

    printf("[%10.2lf] ", get_current_time() - initial_time);

    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);

    fflush(stdout);
}

void next_set(int sz, int *set)
{
    int i;
    for (i = 0; i < sz - 1 && set[i] + 1 == set[i + 1]; i++) set[i] = i;
    set[i] = set[i] + 1;
}