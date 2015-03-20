#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <float.h>
#include "lp.h"
#include "util.h"
#include "main.h"
#include "tsp.h"
#include "branch_and_cut.h"
#include "gtsp.h"
#include "flow.h"

char *INPUT_FILENAME = 0;
unsigned int SEED = 0;
int GEOMETRIC_DATA = 0;
int NODE_COUNT_RAND = 0;
int GRID_SIZE_RAND = 100;

static int parse_arguments_tsp(int ac, char **av);

static void print_usage_tsp(char *f);

int test_max_flow()
{
    int rval = 0;

    int *edges = 0;
    double *capacities = 0;
    double *flow = 0;
    double flow_value;

    FILE *f = fopen("tmp/flow.in", "r");
    abort_if(!f, "could not open input file");

    struct Graph graph;
    graph_init(&graph);

    int node_count, edge_count;

    rval = fscanf(f, "%d %d ", &node_count, &edge_count);
    abort_if(rval != 2, "invalid input format");

    edges = (int *) malloc(4 * edge_count * sizeof(int));
    abort_if(!edges, "could not allocate edges\n");

    capacities = (double *) malloc(2 * edge_count * sizeof(double));
    abort_if(!capacities, "could not allocate capacities");

    for (int i = 0; i < edge_count; i++)
    {
        int from, to, cap;

        rval = fscanf(f, "%d %d %d ", &from, &to, &cap);
        abort_if(rval != 3, "invalid input format");

        edges[i*4] = edges[i*4+3] = from;
        edges[i*4+1] = edges[i*4 + 2] = to;
        capacities[2 * i] = cap;
        capacities[2 * i + 1] = 0;
    }

    rval = graph_build(node_count, 2 * edge_count, edges, 1, &graph);
    abort_if(rval, "graph_build failed");

    for (int i = 0; i < edge_count; i++)
    {
        graph.edges[2*i].reverse = &graph.edges[2*i+1];
        graph.edges[2*i+1].reverse = &graph.edges[2*i];
    }

    flow = (double *) malloc(graph.edge_count * sizeof(double));
    abort_if(!flow, "could not allocate flow");

    struct Node *from = &graph.nodes[0];
    struct Node *to = &graph.nodes[graph.node_count - 1];

    rval = flow_find_max_flow(&graph, capacities, from, to, flow, &flow_value);
    abort_if(rval, "flow_find_max_flow failed");

    log_info("Optimal flow has value %f\n", flow_value);
    for (int i = 0; i < graph.edge_count; i++)
    {
        struct Edge *e = &graph.edges[i];
        if(flow[e->index] <= 0) continue;

        log_info("  %d %d %6.2f / %6.2f\n", e->from->index, e->to->index, flow[e->index], capacities[e->index]);
    }

    CLEANUP:
    if (capacities) free(capacities);
    if (edges) free(edges);
    if (flow) free(flow);
    return rval;
}

int main_tsp(int ac, char **av)
{
    int rval = 0;

    SEED = (unsigned int) get_real_time();

    struct BNC bnc;
    struct TSPData data;

    rval = TSP_init_data(&data);
    abort_if(rval, "TSP_init_data failed");

    rval = BNC_init(&bnc);
    abort_if(rval, "BNC_init failed");

    rval = parse_arguments_tsp(ac, av);
    abort_if(rval, "Failed to parse arguments.");

    printf("Seed = %d\n", SEED);
    srand(SEED);

    rval = TSP_read_problem(INPUT_FILENAME, &data);
    abort_if(rval, "TSP_read_problem failed");

    bnc.best_obj_val = TSP_find_initial_solution(&data);
    bnc.problem_data = (void *) &data;
    bnc.problem_init_lp = (int (*)(struct LP *, void *)) TSP_init_lp;
    bnc.problem_add_cutting_planes =
            (int (*)(struct LP *, void *)) TSP_add_cutting_planes;

    rval = BNC_init_lp(&bnc);
    abort_if(rval, "BNC_init_lp failed");

    rval = BNC_solve(&bnc);
    abort_if(rval, "BNC_solve_node failed");

    log_info("Optimal integral solution:\n");
    log_info("    obj value = %.2lf **\n", bnc.best_obj_val);

    CLEANUP:
    BNC_free(&bnc);
    TSP_free_data(&data);
    return rval;
}

int main_gtsp(int ac, char **av)
{
    int rval = 0;

    SEED = (unsigned int) get_real_time();
    srand(SEED);

    int node_count = 50;
    int cluster_count = node_count / 5;
    int grid_size = 100;

    struct BNC bnc;
    struct GTSP data;

    rval = GTSP_init_data(&data);
    abort_if(rval, "GTSP_init_data failed");

    rval = GTSP_create_random_problem(node_count, cluster_count, grid_size,
            &data);
    abort_if(rval, "GTSP_create_random_problem failed");

    rval = GTSP_write_data(&data, "gtsp.in");
    abort_if(rval, "GTSP_write_problem failed");

    rval = BNC_init(&bnc);
    abort_if(rval, "BNC_init failed");

    log_info("Setting seed = %d\n", SEED);
    srand(SEED);

    bnc.best_obj_val = DBL_MAX;
    bnc.problem_data = (void *) &data;
    bnc.problem_init_lp = (int (*)(struct LP *, void *)) GTSP_init_lp;

    rval = BNC_init_lp(&bnc);
    abort_if(rval, "BNC_init_lp failed");

    log_info("Starting branch-and-cut solver...\n");
    rval = BNC_solve(&bnc);
    abort_if(rval, "BNC_solve_node failed");

    log_info("Optimal integral solution:\n");
    log_info("    obj value = %.2lf **\n", bnc.best_obj_val);

    rval = GTSP_write_solution(&data, "gtsp.out", bnc.best_x);
    abort_if(rval, "GTSP_write_solution failed");

    CLEANUP:
    GTSP_free(&data);
    BNC_free(&bnc);
    return rval;

}

static int parse_arguments_tsp(int ac, char **av)
{
    int rval = 0;

    int c;

    if (ac == 1)
    {
        print_usage_tsp(av[0]);
        return 1;
    }

    while ((c = getopt(ac, av, "ab:gk:s:")) != EOF)
    {
        switch (c)
        {
            case 'a':;
                break;
            case 'b':
                GRID_SIZE_RAND = atoi(optarg);
                break;
            case 'g':
                GEOMETRIC_DATA = 1;
                break;
            case 'k':
                NODE_COUNT_RAND = atoi(optarg);
                break;
            case 's':
                SEED = (unsigned) atoi(optarg);
                break;
            case '?':
            default:
                print_usage_tsp(av[0]);
                return 1;
        }
    }

    if (optind < ac) INPUT_FILENAME = av[optind++];

    if (optind != ac)
    {
        print_usage_tsp(av[0]);
        return 1;
    }

    abort_if(!INPUT_FILENAME && !NODE_COUNT_RAND,
            "Must specify an input file or use -k for random problem\n");

    CLEANUP:
    return rval;
}

int main(int ac, char **av)
{
    return test_max_flow();
//    return main_gtsp(ac, av);
//  return main_tsp(ac, av);
}

static void print_usage_tsp(char *f)
{
    fprintf(stderr, "Usage: %s [-see below-] [prob_file]\n"
            "   -a    add all subtours cuts at once\n"
            "   -b d  gridsize d for random problems\n"
            "   -g    prob_file has x-y coordinates\n"
            "   -k d  generate problem with d cities\n"
            "   -s d  random SEED\n", f);
}

