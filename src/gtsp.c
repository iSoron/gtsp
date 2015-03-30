#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <getopt.h>
#include <math.h>
#include "gtsp.h"
#include "geometry.h"
#include "util.h"
#include "flow.h"
#include "branch_and_cut.h"
#include "gtsp-subtour.h"

double *OPTIMAL_X = 0;

int GTSP_init_data(struct GTSP *data)
{
    int rval = 0;

    data->clusters = 0;
    data->cluster_count = 0;

    data->graph = (struct Graph *) malloc(sizeof(struct Graph));
    abort_if(!data->graph, "could not allocate data->graph");

    graph_init(data->graph);

    CLEANUP:
    return rval;
}

void GTSP_free(struct GTSP *data)
{
    if (!data) return;

    graph_free(data->graph);
    free(data->graph);

    if (data->clusters) free(data->clusters);
}

int GTSP_create_random_problem(
        int node_count, int cluster_count, int grid_size, struct GTSP *data)
{
    int rval = 0;

    int *edges = 0;
    int *weights = 0;
    int *clusters = 0;

    double *x_coords = 0;
    double *y_coords = 0;

    struct Graph *graph = 0;

    int edge_count = (node_count * (node_count - 1)) / 2;

    graph = (struct Graph *) malloc(sizeof(struct Graph));
    abort_if(!graph, "could not allocate graph\n");

    graph_init(graph);

    edges = (int *) malloc(2 * edge_count * sizeof(int));
    weights = (int *) malloc(edge_count * sizeof(int));
    clusters = (int *) malloc(node_count * sizeof(int));

    abort_if(!edges, "could not allocate data->edges\n");
    abort_if(!weights, "could not allocate weights\n");
    abort_if(!clusters, "could not allocate clusters\n");

    x_coords = (double *) malloc(node_count * sizeof(double));
    y_coords = (double *) malloc(node_count * sizeof(double));

    abort_if(!x_coords, "could not allocate x_coords\n");
    abort_if(!y_coords, "could not allocate y_coords\n");

    rval = generate_random_clusters_2d(node_count, cluster_count, grid_size,
            x_coords, y_coords, clusters);
    abort_if(rval, "generate_random_clusters_2d failed");

    int curr_edge = 0;
    for (int i = 0; i < edge_count; i++)
        for (int j = i + 1; j < node_count; j++)
        {
            if (clusters[i] == clusters[j]) continue;

            edges[curr_edge * 2] = i;
            edges[curr_edge * 2 + 1] = j;
            weights[curr_edge] = get_euclidean_distance(x_coords, y_coords, i,
                    j);

            curr_edge++;
        }

    edge_count = curr_edge;

    rval = graph_build(node_count, edge_count, edges, 0, graph);
    abort_if(rval, "graph_build failed");

    for (int i = 0; i < edge_count; i++)
        graph->edges[i].weight = weights[i];

    data->graph = graph;
    data->clusters = clusters;
    data->cluster_count = cluster_count;
    graph->x_coordinates = x_coords;
    graph->y_coordinates = y_coords;

    CLEANUP:
    if (weights) free(weights);
    if (edges) free(edges);
    if (rval)
    {
        if (clusters) free(clusters);
    }
    return rval;
}

int GTSP_init_lp(struct LP *lp, struct GTSP *data)
{
    int rval = 0;

    int node_count = data->graph->node_count;
    int edge_count = data->graph->edge_count;
    int cluster_count = data->cluster_count;
    int *clusters = data->clusters;
    struct Edge *edges = data->graph->edges;

    for (int i = 0; i < node_count; i++)
    {
        rval = LP_new_row(lp, 'E', 0.0);
        abort_if(rval, "LP_new_row failed");
    }

    for (int i = 0; i < cluster_count; i++)
    {
        rval = LP_new_row(lp, 'E', 1.0);
        abort_if(rval, "LP_new_row failed");
    }

    double lb = 0.0;
    double ub = 1.0;
    int cmatbeg = 0;

    for (int i = 0; i < node_count; i++)
    {
        double obj = 0.0;
        double cmatval[] = {-2.0, 1.0};
        int cmatind[] = {i, node_count + clusters[i]};

        rval = LP_add_cols(lp, 1, 2, &obj, &cmatbeg, cmatind, cmatval, &lb,
                &ub);
        abort_if(rval, "LP_add_cols failed");
    }

    for (int i = 0; i < edge_count; i++)
    {
        double obj = (double) edges[i].weight;
        double cmatval[] = {1.0, 1.0};
        int cmatind[] = {edges[i].from->index, edges[i].to->index};

        rval = LP_add_cols(lp, 1, 2, &obj, &cmatbeg, cmatind, cmatval, &lb,
                &ub);
        abort_if(rval, "LP_add_cols failed");
    }

    CLEANUP:
    return rval;
}

int GTSP_add_cutting_planes(struct LP *lp, struct GTSP *data)
{
    int rval = 0;
    int current_round = 0;

    while (1)
    {
        if (current_round > 0)
        {
            int is_infeasible;
            rval = LP_optimize(lp, &is_infeasible);
            abort_if(rval, "LP_optimize failed");

            if (is_infeasible) break;
        }

        current_round++;

        int original_cut_pool_size;
        int added_cuts_count;

        original_cut_pool_size = lp->cut_pool_size;
        log_debug("Finding subtour cuts, round %d...\n", current_round);

        rval = find_exact_subtour_cuts(lp, data, LP_EPSILON);
        abort_if(rval, "find_exact_subtour_cuts failed");

        added_cuts_count = lp->cut_pool_size - original_cut_pool_size;
        if (added_cuts_count > 0)
            continue;

#ifdef ENABLE_COMB_INEQUALITIES
        original_cut_pool_size = lp->cut_pool_size;
        log_debug("Finding comb cuts, round %d...\n", current_round);

        rval = find_comb_cuts(lp, data);
        abort_if(rval, "find_comb_cuts failed");

        added_cuts_count = lp->cut_pool_size - original_cut_pool_size;
        if (added_cuts_count > 0)
            continue;
#endif

        break;
    }

    CLEANUP:
    return rval;
}

int GTSP_write_problem(struct GTSP *data, char *filename)
{
    int rval = 0;

    FILE *file;

    file = fopen(filename, "w");
    abort_if(!file, "could not open file");

    const struct Graph *graph = data->graph;

    fprintf(file, "%d %d\n", graph->node_count, data->cluster_count);

    for (int i = 0; i < graph->node_count; i++)
    {
        fprintf(file, "%.2lf %.2lf %d\n", graph->x_coordinates[i],
                graph->y_coordinates[i], data->clusters[i]);
    }

    CLEANUP:
    if (file) fclose(file);
    return rval;
}

int GTSP_write_solution(struct GTSP *data, char *filename, double *x)
{
    int rval = 0;

    struct Edge *edges = data->graph->edges;
    int node_count = data->graph->node_count;
    int edge_count = data->graph->edge_count;

    FILE *file;
    file = fopen(filename, "w");
    abort_if(!file, "could not open file");

    int positive_edge_count = 0;
    for (int i = 0; i < edge_count; i++)
        if (x[i + node_count] > LP_EPSILON)
            positive_edge_count++;

    fprintf(file, "%d %d\n", node_count, edge_count);

    fprintf(file, "%d\n", positive_edge_count);

    for (int i = 0; i < edge_count; i++)
        if (x[i + node_count] > LP_EPSILON)
            fprintf(file, "%d %d %.4lf\n", edges[i].from->index,
                    edges[i].to->index, x[i + node_count]);

    CLEANUP:
    if (file) fclose(file);
    return rval;
}

int GTSP_read_solution(struct GTSP *gtsp, char *filename, double **p_x)
{
    int rval = 0;

    int node_count;
    int edge_count;
    int *edge_map = 0;

    double *x;

    FILE *file;

    log_info("Reading optimal solution from file %s\n", filename);

    file = fopen(filename, "r");
    abort_if(!file, "could not open file");

    rval = fscanf(file, "%d %d", &node_count, &edge_count);
    abort_if(rval != 2, "invalid input format (node and edge count)");

    int num_cols = node_count + edge_count;

    x = (double *) malloc(num_cols * sizeof(double));
    abort_if(!x, "could not allocate x");

    for (int i = 0; i < node_count + edge_count; i++) x[i] = 0.0;

    rval = fscanf(file, "%d", &edge_count);
    abort_if(rval != 1, "invalid input format (positive edge count)");

    edge_map = (int *) malloc(node_count * node_count * sizeof(int));
    abort_if(!edge_map, "could not allocate edge_map");

    int k = node_count;
    for (int i = 0; i < node_count; i++)
    {
        for (int j = i + 1; j < node_count; j++)
        {
            if (gtsp->clusters[i] == gtsp->clusters[j]) continue;
            edge_map[i * node_count + j] = k;
            edge_map[j * node_count + i] = k;
            k++;
        }
    }

    for (int i = 0; i < edge_count; i++)
    {
        int from, to, edge;
        double val;
        rval = fscanf(file, "%d %d %lf", &from, &to, &val);
        abort_if(rval != 3, "invalid input format (edge endpoints)");

        edge = edge_map[from * node_count + to];
        abort_if(edge > num_cols, "invalid edge");

        x[from] += val / 2;
        x[to] += val / 2;
        x[edge] = val;
    }

    for (int i = 0; i < num_cols; i++)
    {
        if (x[i] <= LP_EPSILON) continue;
        log_debug("    x%-5d = %.6f\n", i, x[i]);
    }

    *p_x = x;
    rval = 0;

    CLEANUP:
    if (edge_map) free(edge_map);
    return rval;
}

int GTSP_check_solution(struct GTSP *data, double *x)
{
    int rval = 0;
    int *cluster_mark = 0;

    struct Node **stack = 0;
    int stack_top = 0;

    struct Graph *graph = data->graph;
    const int node_count = graph->node_count;
    const int edge_count = graph->edge_count;

    cluster_mark = (int *) malloc(data->cluster_count * sizeof(int));
    abort_if(!cluster_mark, "could not allocate cluster_mark");

    stack = (struct Node **) malloc(graph->node_count * sizeof(struct Node *));
    abort_if(!stack, "could not allocate stack");

    for (int i = 0; i < node_count + edge_count; i++)
    {
        abort_iff(x[i] < 1.0 - LP_EPSILON && x[i] > LP_EPSILON,
                "solution is not integral: x%d = %.4lf", i, x[i]);

        abort_iff(x[i] > 1.0 + LP_EPSILON || x[i] < 0.0 - LP_EPSILON,
                "value out of bounds: x%d = %.4lf", i, x[i]);
    }

    for (int i = 0; i < node_count; i++)
        graph->nodes[i].mark = 0;

    for (int i = 0; i < data->cluster_count; i++)
        cluster_mark[i] = 0;

    int initial;
    for (initial = 0; initial < node_count; initial++)
        if (x[initial] > 1.0 - LP_EPSILON) break;

    abort_if(initial == node_count, "no initial node");

    stack[stack_top++] = &graph->nodes[initial];
    graph->nodes[initial].mark = 1;

    while (stack_top > 0)
    {
        struct Node *n = stack[--stack_top];
        cluster_mark[data->clusters[n->index]]++;

        for (int i = 0; i < n->degree; i++)
        {
            struct Adjacency *adj = &n->adj[i];
            struct Node *neighbor = adj->neighbor;

            if (neighbor->mark) continue;
            if (x[node_count + adj->edge->index] < LP_EPSILON) continue;

            stack[stack_top++] = neighbor;
            neighbor->mark = 1;
        }
    }

    for (int i = 0; i < data->cluster_count; i++)
        abort_if(cluster_mark[i] != 1, "cluster not visited exactly one time");

    log_info("    solution is valid\n");

    CLEANUP:
    if (stack) free(stack);
    if (cluster_mark) free(cluster_mark);
    return rval;
}

int GTSP_solution_found(struct GTSP *data, double *x)
{
    int rval = 0;

    char filename[100];

    sprintf(filename, "tmp/gtsp-m%d-n%d-s%d.out", data->cluster_count,
            data->graph->node_count, SEED);

    log_info("Writting solution to file %s\n", filename);
    rval = GTSP_write_solution(data, filename, x);
    abort_if(rval, "GTSP_write_solution failed");

    log_info("Checking solution...\n");
    rval = GTSP_check_solution(data, x);
    abort_if(rval, "GTSP_check_solution failed");

    CLEANUP:
    return rval;
}

static const struct option options_tab[] = {{"help", no_argument, 0, 'h'},
        {"nodes", required_argument, 0, 'n'},
        {"clusters", required_argument, 0, 'm'},
        {"grid-size", required_argument, 0, 'g'},
        {"optimal", required_argument, 0, 'x'},
        {"seed", required_argument, 0, 's'},
        {(char *) 0, (int) 0, (int *) 0, (int) 0}};
static int input_node_count = -1;
static int input_cluster_count = -1;
static int grid_size = 100;

static char input_x_filename[1000] = {0};

static void GTSP_print_usage()
{
    printf("Parameters:\n");
    printf("%4s %-13s %s\n", "-n", "--nodes", "number of nodes");
    printf("%4s %-13s %s\n", "-m", "--clusters", "number of clusters");
    printf("%4s %-13s %s\n", "-s", "--seed", "random seed");
    printf("%4s %-13s %s\n", "-g", "--grid-size",
            "size of the box used for generating random points");
    printf("%4s %-13s %s\n", "-x", "--optimal",
            "file containg valid solution (used to assert validity of cuts)");
}

static int GTSP_parse_args(int argc, char **argv)
{
    int rval = 0;

    opterr = 0;

    while (1)
    {
        int c = 0;
        int option_index = 0;
        c = getopt_long(argc, argv, "n:m:g:x:s:", options_tab, &option_index);

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

            case ':':
                fprintf(stderr, "option '-%c' requires an argument\n", optopt);
                rval = 1;
                goto CLEANUP;

            case '?':
            default:
                fprintf(stderr, "option '-%c' is invalid\n", optopt);
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
        printf("You must specify the number of nodes.\n");
        rval = 1;
    }

    if (input_cluster_count > input_node_count)
    {
        printf("Number of clusters must be at most number of nodes.\n");
        rval = 1;
    }

    if (rval)
    {
        GTSP_print_usage();
        rval = 1;
    }

    CLEANUP:
    return rval;
}

double FLOW_CPU_TIME = 0;
double LP_SOLVE_TIME = 0;
double LP_CUT_POOL_TIME = 0;
int LP_OPTIMIZE_COUNT = 0;

int GTSP_main(int argc, char **argv)
{
    int rval = 0;

    struct BNC bnc;
    struct GTSP data;

    SEED = (unsigned int) get_real_time() % 1000000;

    rval = GTSP_init_data(&data);
    abort_if(rval, "GTSP_init_data failed");

    rval = BNC_init(&bnc);
    abort_if(rval, "BNC_init failed");

    rval = GTSP_parse_args(argc, argv);
    if (rval) return 1;

    srand(SEED);

    log_info("Generating random GTSP instance...\n");
    log_info("    seed = %d\n", SEED);
    log_info("    input_node_count = %d\n", input_node_count);
    log_info("    input_cluster_count = %d\n", input_cluster_count);
    log_info("    grid_size = %d\n", grid_size);

    rval = GTSP_create_random_problem(input_node_count, input_cluster_count,
            grid_size, &data);
    abort_if(rval, "GTSP_create_random_problem failed");

    char filename[100];
    sprintf(filename, "input/gtsp-m%d-n%d-s%d.in", input_cluster_count,
            input_node_count, SEED);
    log_info("Writing random instance to file %s\n", filename);
    rval = GTSP_write_problem(&data, filename);
    abort_if(rval, "GTSP_write_problem failed");

    bnc.best_obj_val = DBL_MAX;
    bnc.problem_data = (void *) &data;
    bnc.problem_init_lp = (int (*)(struct LP *, void *)) GTSP_init_lp;
    bnc.problem_add_cutting_planes = (int (*)(
            struct LP *, void *)) GTSP_add_cutting_planes;
    bnc.problem_solution_found = (int (*)(
            void *, double *)) GTSP_solution_found;

    double opt_val = 0.0;

    if (strlen(input_x_filename) == 0)
    {
        sprintf(input_x_filename, "optimal/gtsp-m%d-n%d-s%d.out",
                input_cluster_count, input_node_count, SEED);

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
            opt_val += OPTIMAL_X[i + input_node_count] * e->weight;
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
        abort_iff(bnc.best_obj_val > opt_val,
                "Solution is not optimal: %.4lf > %.4lf", bnc.best_obj_val,
                opt_val);
    }

    log_info("Branch-and-bound nodes: %d\n", BNC_NODE_COUNT);
    log_info("Max-flow calls: %d\n", FLOW_MAX_FLOW_COUNT);
    log_info("Max-flow computation time: %.2lf\n", FLOW_CPU_TIME);
    log_info("LP optimize calls: %d\n", LP_OPTIMIZE_COUNT);
    log_info("LP solving time: %.2lf\n", LP_SOLVE_TIME);
    log_info("LP cut pool management time: %.2lf\n", LP_CUT_POOL_TIME);

    CLEANUP:
    GTSP_free(&data);
    BNC_free(&bnc);
    return rval;
}
