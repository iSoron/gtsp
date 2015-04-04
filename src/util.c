/* Copyright (c) 2015 Alinson Xavier
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <sys/resource.h>
#include <stdarg.h>
#include <time.h>
#include "util.h"
#include "main.h"

double get_user_time()
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

void time_printf(const char *fmt, ...)
{
    if (INITIAL_TIME == 0)
        INITIAL_TIME = get_user_time();

    printf("[%10.2lf] ", get_user_time() - INITIAL_TIME);

    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);

    fflush(stdout);
}

void pause()
{
    printf("Prese [Enter] to continue...");
    fflush(stdout);
    getchar();
}