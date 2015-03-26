#include <malloc.h>
#include <float.h>
#include "flow.h"
#include "gtsp.h"
#include "util.h"

int FLOW_MAX_FLOW_COUNT = 0;

int flow_mark_reachable_nodes(
        const struct Graph *graph, double *residual_caps, struct Node *from)
{
    int rval = 0;

    struct Node **stack;
    int stack_top = 0;
    int *parents = 0;

    stack = (struct Node **) malloc(graph->node_count * sizeof(struct Node *));
    abort_if(!stack, "could not allocate stack");

    parents = (int *) malloc(graph->node_count * sizeof(int ));
    abort_if(!parents, "could not allocate parents");

    stack[stack_top++] = from;
    from->mark = 1;

    while (stack_top > 0)
    {
        struct Node *n = stack[--stack_top];

        for (int j = 0; j < n->degree; j++)
        {
            struct Edge *e = n->adj[j].edge;
            struct Node *neighbor = n->adj[j].neighbor;

            if (neighbor->mark) continue;
            if (residual_caps[e->index] <= 0) continue;

            stack[stack_top++] = neighbor;
            neighbor->mark = 1;
            parents[neighbor->index] = n->index;
        }

    }

    log_verbose("Reachable nodes:\n");
    for (int i = 0; i < graph->node_count; i++)
        if (graph->nodes[i].mark)
            log_verbose("    %d from %d\n", graph->nodes[i].index, parents[i]);

    CLEANUP:
    if(parents) free(parents);
    if (stack) free(stack);
    return rval;
}

extern double FLOW_CPU_TIME;

int flow_find_max_flow(
        const struct Graph *digraph,
        const double *capacities,
        struct Node *from,
        struct Node *to,
        double *flow,
        double *value)
{
    int rval = 0;

    FLOW_MAX_FLOW_COUNT++;
    double initial_time = get_current_time();

    for (int i = 0; i < digraph->node_count; i++)
        digraph->nodes[i].mark = 0;

    log_verbose("Input graph:\n");

    #if LOG_LEVEL >= LOG_LEVEL_VERBOSE
    graph_dump(digraph);
    #endif

    log_verbose("Solving flow problem:\n");

    log_verbose("%d %d\n", digraph->node_count, digraph->edge_count);
    log_verbose("%d %d\n", from->index, to->index);
    for (int i = 0; i < digraph->edge_count; i++)
    {
        log_verbose("%d %d %.4lf\n", digraph->edges[i].from->index,
                digraph->edges[i].to->index, capacities[i]);
    }

    int path_length;
    struct Edge **path_edges = 0;
    double path_capacity;

    double *residual_caps = 0;

    residual_caps = (double *) malloc(digraph->edge_count * sizeof(double));
    abort_if(!residual_caps, "could not allocate residual_caps");

    path_edges = (struct Edge **) malloc(
            digraph->edge_count * sizeof(struct Edge *));
    abort_if(!path_edges, "could not allocate path_edges");

    for (int i = 0; i < digraph->edge_count; i++)
    {
        flow[i] = 0;
        residual_caps[i] = capacities[i];
        abort_if(!digraph->edges[i].reverse,
                "digraph must have reverse edge information");
        abort_if(digraph->edges[i].reverse->reverse != &digraph->edges[i],
                "invalid reverse edge");
    }

    *value = 0;

    while (1)
    {
        flow_find_augmenting_path(digraph, residual_caps, from, to,
                &path_length, path_edges, &path_capacity);

        if (path_length == 0) break;

        log_verbose("Found augmenting path of capacity %.4lf:\n",
                path_capacity);

        (*value) += path_capacity;

        for (int i = 0; i < path_length; i++)
        {
            struct Edge *e = &digraph->edges[path_edges[i]->index];

            log_verbose("  %d %d (%d)\n", e->from->index, e->to->index, e->index);

            residual_caps[e->index] -= path_capacity;
            residual_caps[e->reverse->index] += path_capacity;

            flow[e->index] += path_capacity;
            flow[e->reverse->index] -= path_capacity;
        }

        log_verbose("New residual capacities:\n");
        for (int i = 0; i < digraph->edge_count; i++)
        {
            #if LOG_LEVEL >= LOG_LEVEL_VERBOSE
            struct Edge *e = &digraph->edges[i];
            #endif

            if (residual_caps[i] < LP_EPSILON) continue;

            log_verbose("%d %d %.4lf (%d)\n", e->from->index, e->to->index, e->index,
                    residual_caps[e->index]);
        }
    }

    log_verbose("No more paths found.\n");

    rval = flow_mark_reachable_nodes(digraph, residual_caps, from);
    abort_if(rval, "flow_mark_reachable_nodes failed");

    FLOW_CPU_TIME += get_current_time() - initial_time;

    CLEANUP:
    if (path_edges) free(path_edges);
    if (residual_caps) free(residual_caps);
    return rval;
}

int flow_find_augmenting_path(
        const struct Graph *graph,
        const double *residual_caps,
        struct Node *from,
        struct Node *to,
        int *path_length,
        struct Edge **path_edges,
        double *path_capacity)
{
    int rval = 0;

    struct Node **queue = 0;
    int queue_start = 0;
    int queue_end = 0;

    struct Node **parents = 0;
    struct Edge **parent_edges = 0;

    int node_count = graph->node_count;

    queue = (struct Node **) malloc(node_count * sizeof(struct Node *));
    parents = (struct Node **) malloc(node_count * sizeof(struct Node *));
    parent_edges = (struct Edge **) malloc(node_count * sizeof(struct Edge *));

    abort_if(!queue, "could not allocate queue");
    abort_if(!parents, "could not allocate parents");
    abort_if(!parent_edges, "could not allocate parent_edges");

    for (int i = 0; i < node_count; i++)
        graph->nodes[i].mark = 0;

    int found = 0;
    queue[queue_end++] = from;

    while (queue_end > queue_start)
    {
        struct Node *n = queue[queue_start++];

        n->mark = 1;

        for (int i = 0; i < n->degree; i++)
        {
            struct Node *neighbor = n->adj[i].neighbor;
            struct Edge *edge = n->adj[i].edge;

            if (neighbor->mark > 0) continue;
            if (residual_caps[edge->index] < LP_EPSILON) continue;

            parents[neighbor->index] = n;
            parent_edges[neighbor->index] = edge;

            queue[queue_end++] = neighbor;
            neighbor->mark = 1;

            if (neighbor == to)
            {
                found = 1;
                break;
            }
        }

        if (found) break;
    }

    *path_length = 0;
    *path_capacity = DBL_MAX;

    if (queue_end == queue_start) goto CLEANUP;

    struct Node *n = to;
    while (n != from)
    {
        struct Edge *edge = parent_edges[n->index];
        path_edges[*path_length] = edge;

        double c = residual_caps[edge->index];
        if (c < *path_capacity) *path_capacity = c;

        n = parents[n->index];
        (*path_length)++;
    }

    CLEANUP:
    if (parents) free(parents);
    if (parent_edges) free(parent_edges);
    if (queue) free(queue);
    return rval;
}

int flow_main(int argc, char **argv)
{
    UNUSED(argc);
    UNUSED(argv);

    int rval = 0;

    int *edges = 0;
    double *capacities = 0;
    double *flow = 0;
    double flow_value;

    struct Edge **cut_edges = 0;

    FILE *f = fopen("tmp/flow.in", "r");
    abort_if(!f, "could not open input file");

    struct Graph graph;
    graph_init(&graph);

    int node_count, edge_count;

    rval = fscanf(f, "%d %d ", &node_count, &edge_count);
    abort_if(rval != 2, "invalid input format (node count, edge count)");

    int sink, source;
    rval = fscanf(f, "%d %d", &source, &sink);
    abort_if(rval != 2, "invalid input format (source, sink)");

    edges = (int *) malloc(4 * edge_count * sizeof(int));
    abort_if(!edges, "could not allocate edges\n");

    capacities = (double *) malloc(2 * edge_count * sizeof(double));
    abort_if(!capacities, "could not allocate capacities");

    for (int i = 0; i < edge_count; i++)
    {
        int from, to;
        double cap;

        rval = fscanf(f, "%d %d %lf ", &from, &to, &cap);
        abort_if(rval != 3, "invalid input format (edge specification)");

        edges[i * 4] = edges[i * 4 + 3] = from;
        edges[i * 4 + 1] = edges[i * 4 + 2] = to;
        capacities[2 * i] = cap;
        capacities[2 * i + 1] = 0;
    }

    rval = graph_build(node_count, 2 * edge_count, edges, 1, &graph);
    abort_if(rval, "graph_build failed");

    for (int i = 0; i < edge_count; i++)
    {
        graph.edges[2 * i].reverse = &graph.edges[2 * i + 1];
        graph.edges[2 * i + 1].reverse = &graph.edges[2 * i];
    }

    flow = (double *) malloc(graph.edge_count * sizeof(double));
    abort_if(!flow, "could not allocate flow");

    struct Node *from = &graph.nodes[source];
    struct Node *to = &graph.nodes[sink];

    rval = flow_find_max_flow(&graph, capacities, from, to, flow, &flow_value);
    abort_if(rval, "flow_find_max_flow failed");

    log_info("Optimal flow has value %f\n", flow_value);
    for (int i = 0; i < graph.edge_count; i++)
    {
        struct Edge *e = &graph.edges[i];
        if (flow[e->index] <= 0) continue;

        log_info("  %d %d %6.2f / %6.2f\n", e->from->index, e->to->index,
                flow[e->index], capacities[e->index]);
    }

    log_info("Nodes reachable from origin on residual graph:\n");
    for (int i = 0; i < graph.node_count; i++)
    {
        struct Node *n = &graph.nodes[i];
        if (n->mark)
            log_info("  %d\n", n->index);
    }

    int cut_edges_count = 0;
    cut_edges =
            (struct Edge **) malloc(graph.edge_count * sizeof(struct Edge *));
    abort_if(!cut_edges, "could not allocate cut_edges");

    rval = get_cut_edges_from_marks(&graph, &cut_edges_count, cut_edges);
    abort_if(rval, "get_cut_edges_from_marks failed");

    log_info("Min cut edges:\n");
    for (int i = 0; i < cut_edges_count; i++)
    {
        struct Edge *e = cut_edges[i];
        if (capacities[e->index] <= 0) continue;
        log_info("  %d %d\n", e->from->index, e->to->index);
    }

    CLEANUP:
    if (cut_edges) free(cut_edges);
    if (capacities) free(capacities);
    if (edges) free(edges);
    if (flow) free(flow);
    return rval;
}