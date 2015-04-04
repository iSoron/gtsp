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
