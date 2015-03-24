#include <getopt.h>
#include "main.h"
#include "gtsp.h"
#include "tsp.h"
#include "flow.h"

char *INPUT_FILENAME = 0;
unsigned int SEED = 0;
int GEOMETRIC_DATA = 0;
int NODE_COUNT_RAND = 0;
int GRID_SIZE_RAND = 100;

static const struct option options_tab[] = {
        {"help", no_argument, 0, 'h'}, {"tsp", no_argument, 0, 't'},
        {"gtsp", no_argument, 0, 'g'}, {"flow", no_argument, 0, 'f'},
        {(char *) 0, (int) 0, (int *) 0, (int) 0}
};

void GTSP_print_usage(char **argv)
{
    printf("wrong usage\n");
}

int main(int argc, char **argv)
{
    int c = 0;
    int option_index = 0;

    c = getopt_long(argc, argv, "htgf", options_tab, &option_index);

    if (c < 0)
    {
        GTSP_print_usage(argv);
        return 1;
    }

    switch (c)
    {
        case 'h':
            GTSP_print_usage(argv);
            return 1;

        case 'f':
            return flow_main(argc, argv);

        case 't':
            return TSP_main(argc, argv);

        case 'g':
            return GTSP_main(argc, argv);
    }

}

