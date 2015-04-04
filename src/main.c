#include <getopt.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "main.h"
#include "gtsp.h"
#include "util.h"

unsigned int SEED = 0;

double SUBTOUR_TIME = 0;
double COMBS_TIME = 0;

double CUT_POOL_TIME = 0;
long CUT_POOL_MAX_MEMORY = 0;

double LP_SOLVE_TIME = 0;

int LP_MAX_COLS = 0;
int LP_MAX_ROWS = 0;
int LP_SOLVE_COUNT = 0;

double TOTAL_TIME = 0;
double ROOT_VALUE = 0;

int SUBTOUR_COUNT = 0;
int COMBS_COUNT = 0;

double INITIAL_TIME = 0;

int BNC_NODE_COUNT = 0;

static int input_node_count = -1;
static int input_cluster_count = -1;
static int grid_size = 100;

static const struct option options_tab[] = {{"help", no_argument, 0, 'h'},
        {"nodes", required_argument, 0, 'n'},
        {"clusters", required_argument, 0, 'm'},
        {"grid-size", required_argument, 0, 'g'},
        {"optimal", required_argument, 0, 'x'},
        {"seed", required_argument, 0, 's'},
        {(char *) 0, (int) 0, (int *) 0, (int) 0}};

static char input_x_filename[1000] = {0};

static void print_usage(int argc, char **argv)
{
    printf("Usage: %s [OPTION]...\n", argv[0]);
    printf("Solves the Generalized Traveling Salesman problem for the input graph.\n\n");

    printf("Parameters:\n");
    printf("%4s %-13s %s\n", "-n", "--nodes", "number of nodes");
    printf("%4s %-13s %s\n", "-m", "--clusters", "number of clusters");
    printf("%4s %-13s %s\n", "-s", "--seed", "random seed");
    printf("%4s %-13s %s\n", "-g", "--grid-size",
            "size of the box used for generating random points");
    printf("%4s %-13s %s\n", "-x", "--optimal",
            "file containg valid solution (used to assert validity of cuts)");
}

static int parse_args(int argc, char **argv)
{
    int rval = 0;

    opterr = 0;

    while (1)
    {
        int c = 0;
        int option_index = 0;
        c = getopt_long(argc, argv, "n:m:g:x:s:h:", options_tab, &option_index);

        if (c < 0) break;

        switch (c)
        {
            case 'n':
                input_node_count = atoi(optarg);
                break;

            case 'm':
                input_cluster_count = atoi(optarg);
                break;

            case 'g':
                grid_size = atoi(optarg);
                break;

            case 'x':
                strcpy(input_x_filename, optarg);
                break;

            case 's':
                SEED = (unsigned) atoi(optarg);
                break;

            case 'h':
                print_usage(argc, argv);
                exit(0);

            case ':':
                fprintf(stderr, "%s: option '-%c' requires an argument\n",
                        argv[0], optopt);
                rval = 1;
                goto CLEANUP;

            case '?':
            default:
                fprintf(stderr, "%s: option '-%c' is invalid\n", argv[0],
                        optopt);
                rval = 1;
                goto CLEANUP;

        }
    }

    if (input_cluster_count < 0)
    {
        input_cluster_count = (int) ceil(input_node_count / 5.0);
        if (input_cluster_count < 3) input_cluster_count = 3;
    }

    if (input_node_count < 0)
    {
        fprintf(stderr, "You must specify the number of nodes.\n");
        rval = 1;
    }
    else if (input_cluster_count > input_node_count)
    {
        printf("Number of clusters must be at most number of nodes.\n");
        rval = 1;
    }

    CLEANUP:
    if (rval)
        fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
    return rval;
}

int main(int argc, char **argv)
{
    int rval = 0;

    struct BNC bnc;
    struct GTSP data;

    double *initial_x = 0;
    double initial_time = get_user_time();

    SEED = (unsigned int) get_real_time() % 1000;

    rval = GTSP_init_data(&data);
    abort_if(rval, "GTSP_init_data failed");

    rval = BNC_init(&bnc);
    abort_if(rval, "BNC_init failed");

    rval = parse_args(argc, argv);
    if (rval) return 1;

    srand(SEED);

    log_info("Generating random GTSP instance...\n");
    log_info("    seed = %d\n", SEED);
    log_info("    input_node_count = %d\n", input_node_count);
    log_info("    input_cluster_count = %d\n", input_cluster_count);
    log_info("    grid_size = %d\n", grid_size);

    char instance_name[1000];
    sprintf(instance_name, "gtsp-m%d-n%d-s%d", input_cluster_count,
            input_node_count, SEED);

    rval = GTSP_create_random_problem(input_node_count, input_cluster_count,
            grid_size, &data);
    abort_if(rval, "GTSP_create_random_problem failed");

    char filename[100];
    sprintf(filename, "input/%s.in", instance_name);
    log_info("Writing random instance to file %s\n", filename);
    rval = GTSP_write_problem(&data, filename);
    abort_if(rval, "GTSP_write_problem failed");

#if LOG_LEVEL >= LOG_LEVEL_DEBUG
    log_info("Writing random instance to file gtsp.in\n");
    rval = GTSP_write_problem(&data, "gtsp.in");
#endif

    int init_val = 0;

    initial_x = (double *) malloc(
            (data.graph->node_count + data.graph->edge_count) * sizeof(double));
    abort_if(!initial_x, "could not allocate initial_x");

    rval = inital_tour_value(&data, &init_val, initial_x);
    abort_if(rval, "initial_tour_value failed");

    rval = GTSP_solution_found(&bnc, &data, initial_x);
    abort_if(rval, "check_sol failed");

    bnc.best_x = initial_x;
    bnc.best_obj_val = init_val;
    bnc.problem_data = (void *) &data;
    bnc.problem_init_lp = (int (*)(struct LP *, void *)) GTSP_init_lp;
    bnc.problem_add_cutting_planes = (int (*)(
            struct LP *, void *)) GTSP_add_cutting_planes;
    bnc.problem_solution_found = (int (*)(
            struct BNC *, void *, double *)) GTSP_solution_found;

    double opt_val = 0.0;

    if (strlen(input_x_filename) == 0)
    {
        sprintf(input_x_filename, "optimal/%s.out", instance_name);

        FILE *file = fopen(input_x_filename, "r");

        if (!file)
            input_x_filename[0] = 0;
        else
            fclose(file);
    }

    if (strlen(input_x_filename) > 0)
    {
        rval = GTSP_read_solution(&data, input_x_filename, &OPTIMAL_X);
        abort_if(rval, "GTSP_read_solution failed");

        log_info("Optimal solution is available. Cuts will be checked.\n");

        for (int i = 0; i < data.graph->edge_count; i++)
        {
            struct Edge *e = &data.graph->edges[i];
            opt_val += OPTIMAL_X[i] * e->weight;
        }

        log_info("    opt = %.2lf\n", opt_val);
    }

    log_info("Initializing LP...\n");
    rval = BNC_init_lp(&bnc);
    abort_if(rval, "BNC_init_lp failed");

//    log_info("Writing LP to file gtsp.lp...\n");
//    rval = LP_write(bnc.lp, "gtsp.lp");
//    abort_if(rval, "LP_write failed");

    log_info("Starting branch-and-cut solver...\n");
    rval = BNC_solve(&bnc);
    abort_if(rval, "BNC_solve_node failed");

    abort_if(!bnc.best_x, "problem has no feasible solution");

    log_info("Optimal integral solution:\n");
    log_info("    obj value = %.2lf **\n", bnc.best_obj_val);

    if (OPTIMAL_X)
    {
        abort_iff(bnc.best_obj_val - EPSILON > opt_val,
                "Solution is not optimal: %.4lf > %.4lf", bnc.best_obj_val,
                opt_val);
    }

    TOTAL_TIME = get_user_time() - initial_time;

    log_info("Branch-and-bound nodes: %d\n", BNC_NODE_COUNT);
    log_info("LP optimize calls: %d\n", LP_SOLVE_COUNT);
    log_info("LP solving time: %.2lf\n", LP_SOLVE_TIME);
    log_info("LP cut pool management time: %.2lf\n", CUT_POOL_TIME);

    // Write statistics to a file
    FILE *file = fopen("stats.tab", "a");
    abort_if(!file, "could not open stats.tab");

    struct stat st;
    stat("stats.tab", &st);

    if (st.st_size == 0)
    {
        fprintf(file, "%-20s  ", "instance");
        fprintf(file, "%-8s  ", "time");
        fprintf(file, "%-8s  ", "subt-t");
        fprintf(file, "%-8s  ", "combs-t");
        fprintf(file, "%-8s  ", "pool-t");
        fprintf(file, "%-8s  ", "pool-m");
        fprintf(file, "%-8s  ", "lp-count");
        fprintf(file, "%-8s  ", "lp-time");
        fprintf(file, "%-8s  ", "lp-rows");
        fprintf(file, "%-8s  ", "lp-cols");
        fprintf(file, "%-8s  ", "init-v");
        fprintf(file, "%-8s  ", "opt-v");
        fprintf(file, "%-8s  ", "root-v");
        fprintf(file, "%-8s  ", "nodes");
        fprintf(file, "%-8s  ", "subt-cc");
        fprintf(file, "%-8s  ", "subt-nc");
        fprintf(file, "%-8s  ", "subt-nn");
        fprintf(file, "%-8s  ", "combs");

        fprintf(file, "\n");
    }

    fprintf(file, "%-20s  ", instance_name);
    fprintf(file, "%-8.2lf  ", TOTAL_TIME);
    fprintf(file, "%-8.2lf  ", SUBTOUR_TIME);
    fprintf(file, "%-8.2lf  ", COMBS_TIME);
    fprintf(file, "%-8.2lf  ", CUT_POOL_TIME);
    fprintf(file, "%-8ld  ", CUT_POOL_MAX_MEMORY / 1024 / 1024);
    fprintf(file, "%-8d  ", LP_SOLVE_COUNT);
    fprintf(file, "%-8.2lf  ", LP_SOLVE_TIME);
    fprintf(file, "%-8d  ", LP_MAX_ROWS);
    fprintf(file, "%-8d  ", LP_MAX_COLS);
    fprintf(file, "%-8d  ", init_val);
    fprintf(file, "%-8.0lf  ", bnc.best_obj_val);
    fprintf(file, "%-8.0lf  ", ROOT_VALUE);
    fprintf(file, "%-8d  ", BNC_NODE_COUNT);
    fprintf(file, "%-8d  ", SUBTOUR_COUNT);
    fprintf(file, "%-8d  ", COMBS_COUNT);
    fprintf(file, "\n");
    fclose(file);

    CLEANUP:
    if (OPTIMAL_X) free(OPTIMAL_X);
    GTSP_free(&data);
    BNC_free(&bnc);
    return rval;
}


