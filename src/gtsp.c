#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <getopt.h>
#include "gtsp.h"
#include "geometry.h"
#include "util.h"
#include "flow.h"
#include "branch_and_cut.h"
#include "math.h"
#include "gtsp-subtour.h"

double *OPTIMAL_X = 0;

static int get_edge_num(int node_count, int from, int to)
{
    int idx = node_count;

    for (int k = 0; k < from; k++)
        idx += node_count - k - 1;

    idx += to - from - 1;

    return idx;
}

int GTSP_init_data(struct GTSP *data)
{
    int rval = 0;

    data->clusters = 0;
    data->cluster_count = 0;
    data->x_coordinates = 0;
    data->y_coordinates = 0;

    data->graph = (struct Graph *) malloc(sizeof(struct Graph));
    abort_if(!data->graph, "could not allocate data->graph");

	data->vertex_set = (struct CLUSTER *) malloc(sizeof(struct CLUSTER));
	
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
    if (data->x_coordinates) free(data->x_coordinates);
    if (data->y_coordinates) free(data->y_coordinates);
}

int GTSP_create_random_problem(
        int node_count, int cluster_count, int grid_size, struct GTSP *data)
{
    int rval = 0;
	int i = 0;
    int *edges = 0;
    int *weights = 0;
    int *clusters = 0;
	
	int ** dist_matrix = 0;

    double *x_coords = 0;
    double *y_coords = 0;

    struct Graph *graph = 0;

    int edge_count = (node_count * (node_count - 1)) / 2;
http://ns-webapp-b.private.uwaterloo.ca/wifisetup/index.html?cmd=login&switchip=fd74:6b6a:8eca:402::6&mac=44:6d:57:16:be:2a&ip=2620:101:f000:701::5d3c:2c75&essid=uw-wifi-setup-no-encrypt&apname=ENG-AP-E5-5--E&apgroup=E5&url=http%3A%2F%2Fwww%2Egstatic%2Ecom%2Fgenerate_204
    graph = (struct Graph *) malloc(sizeof(struct Graph));
    abort_if(!graph, "could not allocate graph\n");

    graph_init(graph);

    edges = (int *) malloc(2 * edge_count * sizeof(int));
    weights = (int *) malloc(edge_count * sizeof(int));
    clusters = (int *) malloc(node_count * sizeof(int));
abort_if(!data->graph, "could not allocate data->graph");
    abort_if(!edges, "could not allocate data->edges\n");
    abort_if(!weights, "could not allocate weights\n");
    abort_if(!clusters, "could not allocate clusters\n");

    x_coords = (double *) malloc(node_count * sizeof(double));
    y_coords = (double *) malloc(node_count * sizeof(double));

    abort_if(!x_coords, "could not allocate x_coords\n");
    abort_if(!y_coords, "could not allocate y_coords\n");

	dist_matrix = (int **) malloc(node_count * sizeof(int*));
	for(i=0; i<node_count; i++)
		dist_matrix[i] = (int *) malloc(node_count * sizeof(int));
	abort_if(!dist_matrix, "could not allocate dist_matrix\n");
    
    
    rval = generate_random_clusters_2d(node_count, cluster_count, grid_size,
            x_coords, y_coords, clusters);
    abort_if(rval, "generate_random_clusters_2d failed");

	rval = generate_dist_matrix(node_count,
            x_coords, y_coords, dist_matrix);
    abort_if(rval, "generate_distance_matrix_2d failed");
    
	struct CLUSTER *cluster_member;
	cluster_member = (struct CLUSTER *) malloc(cluster_count * sizeof(struct CLUSTER));
	for (int j=0; j<cluster_count; j++){
		cluster_member[j].size = 0;
		for(int i=0; i<node_count; i++)
			if(clusters[i] == j)
				cluster_member[j].size+= 1;
	}
	for (int j=0; j<cluster_count; j++)
		cluster_member[j].set = (int *) malloc(cluster_member[j].size * sizeof(int));
		
	int current_vertex = 0;	
	for (int j=0; j<cluster_count; j++){
		current_vertex = 0;
		for(int i=0; i<node_count; i++)
			if(clusters[i] == j){
				cluster_member[j].set[current_vertex] = i;

				current_vertex += 1;
			}
	}
	
			
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
    data->x_coordinates = x_coords;
    data->y_coordinates = y_coords;
	data->dist_matrix = dist_matrix;
	data->vertex_set = cluster_member;
		
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

    int round = 0;

    int violation_total = 3;
    int violation_current = 0;
    double violations[] = {1.0, 0.1, LP_EPSILON};

    while (1)
    {
        round++;

        log_debug("Finding subtour cuts, round %d, violation %.4lf...\n", round,
                violations[violation_current]);

        int original_cut_pool_size = lp->cut_pool_size;
        rval = find_exact_subtour_cuts(lp, data, violations[violation_current]);
        abort_if(rval, "find_exact_subtour_cuts failed");

        if (lp->cut_pool_size - original_cut_pool_size == 0)
        {
            if (++violation_current < violation_total)
            {
                log_debug("No cuts found. Decreasing minimum cut violation.\n");
                continue;
            }
            else
            {
                log_debug("No additional cuts found.\n");
                break;
            }
        }

        int is_infeasible;
        rval = LP_optimize(lp, &is_infeasible);
        abort_if(rval, "LP_optimize failed");

        if (is_infeasible) break;
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

    fprintf(file, "%d %d\n", data->graph->node_count, data->cluster_count);

    for (int i = 0; i < data->graph->node_count; i++)
    {
        fprintf(file, "%.2lf %.2lf %d\n", data->x_coordinates[i],
                data->y_coordinates[i], data->clusters[i]);
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

int GTSP_read_solution(char *filename, double **p_x)
{
    int rval = 0;

    int node_count;
    int edge_count;

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

    for (int i = 0; i < edge_count; i++)
    {
        int from, to, edge;
        rval = fscanf(file, "%d %d", &from, &to);
        abort_if(rval != 2, "invalid input format (edge endpoints)");

        if (from > to) swap(from, to);

        edge = get_edge_num(node_count, from, to);
        abort_if(edge > num_cols, "invalid edge");

        x[from] += 0.5;
        x[to] += 0.5;
        x[edge] = 1;
    }

    for (int i = 0; i < num_cols; i++)
    {
        if (x[i] <= LP_EPSILON) continue;
        log_debug(" x%-3d = %.2f\n", i, x[i]);
    }

    *p_x = x;
    rval = 0;

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
                rval = GTSP_read_solution(optarg, &OPTIMAL_X);
                abort_if(rval, "GTSP_read_solution failed");
                break;

            case 's':
                SEED = (unsigned) atoi(optarg);
                break;

            case ':':
                fprintf(stderr, "option '-%c' requires an argument\n", optopt);
                return 1;

            case '?':
            default:
                fprintf(stderr, "option '-%c' is invalid\n", optopt);
                return 1;

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

int GTSP_solution_found(struct GTSP *data, double *x)
{
    int rval = 0;

    log_info("Writting integral solution to file gtsp.out\n");
    rval = GTSP_write_solution(data, "gtsp.out", x);
    abort_if(rval, "GTSP_write_solution failed");

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
	int init_val ;
	
	init_val = inital_tour_value(&data);
    log_info("Writing random instance to file gtsp.in\n");
    rval = GTSP_write_problem(&data, "gtsp.in");
    abort_if(rval, "GTSP_write_problem failed");

    bnc.best_obj_val = init_val;
    bnc.problem_data = (void *) &data;
    bnc.problem_init_lp = (int (*)(struct LP *, void *)) GTSP_init_lp;
    bnc.problem_add_cutting_planes = (int (*)(
            struct LP *, void *)) GTSP_add_cutting_planes;
    bnc.problem_solution_found = (int (*)(
            void *, double *)) GTSP_solution_found;

    if (OPTIMAL_X)
    {
        log_info("Optimal solution is available. Cuts will be checked.\n");

        double opt_val = 0.0;
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

    log_info("Writing LP to file gtsp.lp...\n");
    rval = LP_write(bnc.lp, "gtsp.lp");
    abort_if(rval, "LP_write failed");

    log_info("Starting branch-and-cut solver...\n");
    rval = BNC_solve(&bnc);
    abort_if(rval, "BNC_solve_node failed");

    abort_if(!bnc.best_x, "problem has no feasible solution");

    log_info("Optimal integral solution:\n");
    log_info("    obj value = %.2lf **\n", bnc.best_obj_val);

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

int inital_tour_value(struct GTSP *data)
{
	int cluster_count = data->cluster_count;
	int * tour;
	int * uncovered_sets;
	int min_vertex = -1;
	int min_cost = 100000000;
	int tour_cost = 0;
	
	int* cluster_in_tour;
	cluster_in_tour = (int *) malloc(cluster_count*sizeof(int));
	
	tour = (int *) malloc(cluster_count*sizeof(int));
	uncovered_sets = (int *) malloc((cluster_count-1)*sizeof(int));
	
	int cluster_num = 0;
	for(int i =0; i< cluster_count; i++){
		cluster_in_tour[i] = 0;
		if(data->clusters[0] != i){
			uncovered_sets[cluster_num] = i;
			cluster_num += 1;
		}
	}
	int new_vertex = 1;
	tour[0] = 0;
	cluster_in_tour[0] = 1;
	while(new_vertex <= data->cluster_count){
		min_vertex = -1;
		min_cost = 100000000;
		for (int i = 1; i < data->graph->node_count; i++) {
			if (cluster_in_tour[data->clusters[i]] == 0){	
					for (int k = 0; k < new_vertex; k++) {
						int cost      = data->dist_matrix[i][tour[k]];
						if (cost < min_cost) {
							min_cost    = cost;
							min_vertex  = i;
						}
				}
			}
		}
		tour[new_vertex] = min_vertex;
		cluster_in_tour[data->clusters[min_vertex]] = 1;
		new_vertex += 1;
	}
	
	

	
	
	tour_cost = Larg_neighborhood_search(tour, data);	
	//tour_cost = optimize_vertex_in_cluster(tour, data);	
	log_info("Initial upper-bound: %d \n", tour_cost);
	return tour_cost;		
}



int optimize_vertex_in_cluster(struct TOUR * tour, struct GTSP *data)
{
	int i = 0 , j, current_cluster, tour_cost;
	int insertion_cost = 1000000;
	int rval = 0;	

	if(rval)
		printf("Larg_neighborhood_search stopped unexpectedly");
		
	for(i = 0; i < data->cluster_count ; i++){
		int vertex = tour[i].vertex;
		int prev_vertex = tour[tour[i].prev].vertex;
		int next_vertex = tour[tour[i].next].vertex;
		
		current_cluster = data->clusters[vertex];
		insertion_cost = data->dist_matrix[prev_vertex][vertex] +
			data->dist_matrix[vertex][next_vertex];
		
		for(j = 0; j < data->vertex_set[current_cluster].size; j++){
				int vertex_in_cluster = data->vertex_set[current_cluster].set[j];
				int cost = data->dist_matrix[vertex_in_cluster][prev_vertex] +
					data->dist_matrix[vertex_in_cluster][next_vertex];
				if (insertion_cost > cost){
					insertion_cost = cost;
					tour[i].vertex = vertex_in_cluster;
				}
		}
	}
	
	return 0;
}

int two_opt(struct TOUR* tour, struct GTSP *data){
	
	int rval = 0, i;
	for (i = 0; i < data->cluster_count; i++){
		int vertex1 = tour[i].vertex;
		int vertex2 = tour[tour[i].prev].vertex;
		int vertex3 = tour[tour[i].next].vertex;
		int vertex4 = tour[tour[tour[i].next].next].vertex;

		int current_cost = data->dist_matrix[vertex2][vertex1] + 
		data->dist_matrix[vertex3][vertex4];
		
		int temp_cost = data->dist_matrix[vertex2][vertex3] + 
		data->dist_matrix[vertex1][vertex4];
		if(current_cost > temp_cost){
			int temp_next = tour[i].next;
			int temp_prev = tour[i].prev;
			
			tour[i].next = tour[temp_next].next;
			tour[i].prev = temp_next;
			
			tour[tour[temp_next].next].prev = i;
			
			tour[temp_next].next = i;
			tour[temp_next].prev = temp_prev;
			
			tour[temp_prev].next = temp_next;
			
		}
	}
	return rval;
}
/*
int K_opt(int* tour, struct GTSP *data){
	int rval = 0, i, k, I, j;
	int current_cost, temp_cost, J, temp_vertex;
	int tour_length = data->cluster_count;
	
	for (i = 1; i < tour_length - 2; i++){
			I = (i+k)%(tour_length);
			if (I == tour_length - 1){
				current_cost = data->dist_matrix[tour[i-1]][tour[i]] + 
				data->dist_matrix[tour[I]][tour[0]];
				temp_cost = data->dist_matrix[tour[i-1]][tour[I]] + 
			data->dist_matrix[tour[i]][tour[0]];
			}else{				
				current_cost = data->dist_matrix[tour[i-1]][tour[i]] + 
				data->dist_matrix[tour[I]][tour[I+1]];
				temp_cost = data->dist_matrix[tour[i-1]][tour[I]] + 
			data->dist_matrix[tour[i]][tour[I+1]];
			}	
			
			if(current_cost > temp_cost){
				log_info("K_opt improved the bound\n");
				for(j = k; j > 0 ; j--){
					if(i + j > tour_length - 1)
						J = i + j - tour_length;
						temp_vertex = tour[i + k - j];
						tour[i + k - j] = tour[J];
						tour[J] = temp_vertex;
					}
			}
	}
	return rval;
}*/

int Larg_neighborhood_search(int* tour, struct GTSP *data){
	int i, best_vertex, best_pose, tour_cost, rval; 
	struct TOUR *vertex_seq;
	vertex_seq = (struct TOUR*) malloc(data->cluster_count*sizeof(struct TOUR));
	
	//Construct the list	
	for(i = 0; i < data->cluster_count; i++){
		vertex_seq[i].vertex = tour[i];
		if ( i == 0){
			vertex_seq[i].prev = data->cluster_count-1;
		}else{
			vertex_seq[i].prev = i - 1;
		}
		if ( i == data->cluster_count-1){
			vertex_seq[i].next = 0;
		}else{
			vertex_seq[i].next = i + 1;
		}
	}

	
	//LNS starts
	for(int iter = 0; iter < 500; iter++){
		//Delete a vertex
		int delete_vertex = rand()%(data->cluster_count - 1) + 1;

		int prev_vertex = vertex_seq[delete_vertex].prev;
		int next_vertex = vertex_seq[delete_vertex].next;
		
		vertex_seq[prev_vertex].next = next_vertex;
		vertex_seq[next_vertex].prev = prev_vertex;
		
		int cluster_to_insert = data->clusters[vertex_seq[delete_vertex].vertex];
		int min_cost = 10000000;
		for(i =0 ; i < data->vertex_set[cluster_to_insert].size ; i++){
			int vertex_to_insert = data->vertex_set[cluster_to_insert].set[i];
			
			int next_edge = vertex_seq[0].next;
			for(int j = 1; j < data->cluster_count ; j++){
				int vertex1 = vertex_seq[next_edge].vertex;
				int vertex2 = vertex_seq[vertex_seq[next_edge].next].vertex;
	
				int insert_cost =  data->dist_matrix[vertex1][vertex_to_insert] +
					data->dist_matrix[vertex_to_insert][vertex2] -
					data->dist_matrix[vertex1][vertex2];
				if(insert_cost < min_cost){
					min_cost = insert_cost;
					best_pose = next_edge;
					best_vertex = vertex_to_insert;
				}
				next_edge = vertex_seq[next_edge].next;
			}
		}

		next_vertex = vertex_seq[best_pose].next;
		vertex_seq[delete_vertex].prev = best_pose;
		vertex_seq[delete_vertex].vertex = best_vertex;
		vertex_seq[delete_vertex].next = next_vertex;
		vertex_seq[best_pose].next =  delete_vertex;
		vertex_seq[next_vertex].prev = delete_vertex;
		
		rval = optimize_vertex_in_cluster(vertex_seq, data);
		
	}
	
	//TWO OPT MOVE
	rval = two_opt(vertex_seq, data);
	
	tour_cost = list_length(vertex_seq, data);	
	
	return tour_cost;
}

int tour_length(int* tour, struct GTSP* data){
	int tour_cost = 0;
		for(int i = 0; i< data->cluster_count ; i++){
			if(i==data->cluster_count-1){
				tour_cost += data->dist_matrix[tour[i]][tour[0]];
			}else{
				tour_cost += data->dist_matrix[tour[i]][tour[i+1]];
				}
			
		}
	return tour_cost;	
}

int list_length(struct TOUR *tour, struct GTSP* data){
	int tour_cost = 0;
		for(int i = 0; i< data->cluster_count ; i++){
				int vertex1 = tour[i].vertex;
				int vertex2 = tour[tour[i].next].vertex;
				tour_cost += data->dist_matrix[vertex1][vertex2];
		}
	return tour_cost;	
}

void print_tour(int* tour, struct GTSP* data){

		for(int i = 0; i< data->cluster_count ; i++){
			printf("%d\t", tour[i]);
			}
			printf("\n");
}

void print_list(struct TOUR * tour, struct GTSP* data){
		printf("%d\t", tour[0].vertex);
		int vertex_next = tour[0].next;
		for(int i = 1; i< data->cluster_count ; i++){
				printf("%d\t", tour[vertex_next].vertex);
				vertex_next = tour[vertex_next].next;
		}
		printf("\n");	
}


