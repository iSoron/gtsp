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

#ifndef PROJECT_PARAMS_H
#define PROJECT_PARAMS_H

/*
 * Error margin for floating point comparisons.
 */
#define EPSILON 0.000001

/*
 * Cut pool settings.
 */
#define MAX_CUT_AGE 5
#define MAX_CUT_POOL_SIZE 1000000

/*
 * Available log levels, in decreasing level of verboseness, are:
 *
 *      LOG_LEVEL_VERBOSE
 *      LOG_LEVEL_DEBUG
 *      LOG_LEVEL_INFO
 *      LOG_LEVEL_WARNING
 *      LOG_LEVEL_ERROR
 */
#define LOG_LEVEL LOG_LEVEL_INFO

/*
 * Time limit for the computation (user time, in seconds).
 */
#define MAX_TOTAL_TIME 3600

#endif //PROJECT_PARAMS_H
