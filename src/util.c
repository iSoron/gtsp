/****************************************************************************/
/*                                                                          */
/*                      Utility Functions for CO759                         */
/*                                                                          */
/****************************************************************************/

#include "main.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <stdarg.h>
#include "util.h"

double util_get_current_time(void)
{
    struct rusage ru;

    getrusage (RUSAGE_SELF, &ru);

    return ((double) ru.ru_utime.tv_sec) +
            ((double) ru.ru_utime.tv_usec)/1000000.0;
}

/* function for creating a random set of points in unit square */

int CO759_build_xy (int ncount, double *xlist, double *ylist, int gridsize)
{
    int rval = 0, i, j, winner, x, y;
    int **hit = (int **) NULL, *hitcount = (int *) NULL;

    printf ("Random %d point set, gridsize = %d\n", ncount, gridsize);
    fflush (stdout);

    hit =  (int **) malloc (gridsize * sizeof (int *));
    if (!hit) {
        fprintf (stderr, "out of memory for hit\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < gridsize; i++) hit[i] = (int *) NULL;

    hitcount = (int *) malloc (gridsize * sizeof (int));
    if (!hitcount) {
        fprintf (stderr, "out of memory for hitcount\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < gridsize; i++) hitcount[i] = 0;

    for (i = 0; i < ncount; i++) {
        winner = 0;
        do {
            x = (int) (rand () % gridsize);
            y = (int) (rand () % gridsize);

            /* check to see if (x,y) is a duplicate point */

            for (j = 0; j < hitcount[x]; j++) {
                if (hit[x][j] == y) break;
            }
            if (j == hitcount[x]) {
                void *tmp_ptr = (void *) hit[x];
                tmp_ptr = realloc (tmp_ptr, (hitcount[x]+1)*sizeof (int));
                if (!tmp_ptr) {
                    fprintf (stderr, "out of member in realloc of hit\n");
                    rval = 1; goto CLEANUP;
                }
                hit[x] = (int *) tmp_ptr;
                hit[x][hitcount[x]] = y;
                hitcount[x]++;
                winner = 1;
            }
            if (!winner) {
                printf ("X"); fflush (stdout);
            }
        } while (!winner);
        xlist[i] = (double) x;
        ylist[i] = (double) y;
    }

CLEANUP:

    printf ("\n");

    if (hit) {
        for (i = 0; i < gridsize; i++) {
            if (hit[i]) free (hit[i]);
        }
        free (hit);
    }
    if (hitcount) free (hitcount);
    return rval;
}

double initial_time = 0;

void time_printf(const char *fmt, ...)
{
    if(initial_time == 0)
        initial_time = util_get_current_time();

    printf("[%10.2lf] ", util_get_current_time() - initial_time);

    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);

    fflush(stdout);
}

void next_set (int sz, int *Set)
{
   int i;
   for (i=0; i < sz-1 && Set[i]+1 == Set[i+1]; i++) Set[i] = i;
   Set[i] = Set[i]+1;
}