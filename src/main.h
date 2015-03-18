#ifndef __MAIN_H_
#define __MAIN_H_

struct AdjObj
{
    int n;
    /* index of neighbor node */
    int e;   /* index of adj joining neighbor */
};

struct Node
{
    int deg;
    struct AdjObj *adj;
    int mark;
};

struct Graph
{
    int node_count;
    int edge_count;
    struct Node *node_list;
    struct AdjObj *adj_space;
};

void usage(char *f);

int parseargs(int ac, char **av);

extern char *fname;
extern int seed;
extern int geometric_data;
extern int ncount_rand;
extern int gridsize_rand;
extern int use_all_subtours;

extern double initial_time;

#endif