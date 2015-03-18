#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <stdarg.h>
#include "util.h"

double get_current_time(void)
{
    struct rusage ru;

    getrusage(RUSAGE_SELF, &ru);

    return ((double) ru.ru_utime.tv_sec)
            + ((double) ru.ru_utime.tv_usec) / 1000000.0;
}

/* function for creating a random set of points in unit square */

int build_random_2d_points(
        int node_count, double *x_list, double *y_list, int grid_size)
{
    int rval = 0, i, j, winner, x, y;
    int **hit = (int **) NULL, *hitcount = (int *) NULL;

    printf("Random %d point set, grid_size = %d\n", node_count, grid_size);
    fflush(stdout);

    hit = (int **) malloc(grid_size * sizeof(int *));
    if (!hit)
    {
        fprintf(stderr, "out of memory for hit\n");
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < grid_size; i++) hit[i] = (int *) NULL;

    hitcount = (int *) malloc(grid_size * sizeof(int));
    if (!hitcount)
    {
        fprintf(stderr, "out of memory for hitcount\n");
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < grid_size; i++) hitcount[i] = 0;

    for (i = 0; i < node_count; i++)
    {
        winner = 0;
        do
        {
            x = (int) (rand() % grid_size);
            y = (int) (rand() % grid_size);

            /* check to see if (x,y) is a duplicate point */

            for (j = 0; j < hitcount[x]; j++)
            {
                if (hit[x][j] == y) break;
            }
            if (j == hitcount[x])
            {
                void *tmp_ptr = (void *) hit[x];
                tmp_ptr = realloc(tmp_ptr, (hitcount[x] + 1) * sizeof(int));
                if (!tmp_ptr)
                {
                    fprintf(stderr, "out of member in realloc of hit\n");
                    rval = 1;
                    goto CLEANUP;
                }
                hit[x] = (int *) tmp_ptr;
                hit[x][hitcount[x]] = y;
                hitcount[x]++;
                winner = 1;
            }
            if (!winner)
            {
                printf("X");
                fflush(stdout);
            }
        } while (!winner);
        x_list[i] = (double) x;
        y_list[i] = (double) y;
    }

    CLEANUP:

    printf("\n");

    if (hit)
    {
        for (i = 0; i < grid_size; i++)
        {
            if (hit[i]) free(hit[i]);
        }
        free(hit);
    }
    if (hitcount) free(hitcount);
    return rval;
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