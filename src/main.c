#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "lp.h"
#include "util.h"
#include "main.h"
#include "tsp.h"
#include "branch_and_cut.h"

char *fname = (char *) NULL;
int seed = 0;
int geometric_data = 0;
int ncount_rand = 0;
int gridsize_rand = 100;
int use_all_subtours = 0;

int main(int ac, char **av)
{
    int rval = 0;
    int *edge_weights = 0;
    int *edge_list = 0;
    int *tlist = 0;

    seed = (int) util_get_current_time();

    rval = parseargs(ac, av);
    ABORT_IF(rval, "Failed to parse arguments.\n");
    ABORT_IF(!fname && !ncount_rand, "Must specify a problem file or use -k for"
            " random prob\n");

    printf("Seed = %d\n", seed);
    srand(seed);

    if (fname)
    {
        printf("Problem name: %s\n", fname);
        if (geometric_data) printf("Geometric data\n");
    }

    int node_count = 0, edge_count = 0;
    rval = TSP_read_problem(fname, &node_count, &edge_count, &edge_list,
            &edge_weights);
    ABORT_IF(rval, "TSP_read_problem failed\n");
    ABORT_IF(use_all_subtours && node_count > 20, "Too many nodes to add all"
            " subtours\n");

    tlist = (int *) malloc((node_count) * sizeof(int));

    ABORT_IF(!tlist, "out of memory for tlist\n");

    double best_val = TSP_find_initial_solution(edge_weights, edge_list,
            node_count, edge_count);

    initial_time = util_get_current_time();

    struct LP *lp;
    lp = (struct LP *) malloc(sizeof(struct LP *));

    rval = bnc_init_lp(lp, node_count, edge_count, edge_list, edge_weights);
    ABORT_IF(rval, "bnc_init_lp failed\n");

    rval = bnc_solve_node(lp, &best_val, node_count, edge_count, edge_list, 1);
    ABORT_IF(rval, "bnc_solve_node failed\n");

    time_printf("Optimal integral solution:\n");
    time_printf("    objective value = %.2lf **\n", best_val);

    printf("\nRunning Time: %.2f seconds\n",
            util_get_current_time() - initial_time);
    fflush(stdout);

    CLEANUP:
    if (tlist) free(tlist);
    if (edge_list) free(edge_list);
    if (edge_weights) free(edge_weights);
    return rval;
}

int parseargs(int ac, char **av)
{
    int c;

    if (ac == 1)
    {
        usage(av[0]);
        return 1;
    }

    while ((c = getopt(ac, av, "ab:gk:s:")) != EOF)
    {
        switch (c)
        {
            case 'a':
                use_all_subtours = 1;
                break;
            case 'b':
                gridsize_rand = atoi(optarg);
                break;
            case 'g':
                geometric_data = 1;
                break;
            case 'k':
                ncount_rand = atoi(optarg);
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case '?':
            default:
                usage(av[0]);
                return 1;
        }
    }

    if (optind < ac) fname = av[optind++];

    if (optind != ac)
    {
        usage(av[0]);
        return 1;
    }

    return 0;
}

void usage(char *f)
{
    fprintf(stderr, "Usage: %s [-see below-] [prob_file]\n"
            "   -a    add all subtours cuts at once\n"
            "   -b d  gridsize d for random problems\n"
            "   -g    prob_file has x-y coordinates\n"
            "   -k d  generate problem with d cities\n"
            "   -s d  random seed\n", f);
}

