#include "include/matmul.h"
#include "include/timer.h"
#include "include/tweets_spatial_analysis.h"
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[]) {
    string ellipse_file("../data/tweets_median_working_filtered.csv");
    string adj_file("./data/tweets_median_working_adjacency_matrix.csv");
    string adj_ordered_file("./data/tweets_median_working_adjacency_matrix_ordered.csv");
    string dis_file("./data/tweets_median_working_distance_matrix.csv");

    int c, rows = 0, cols = 0, n_procs = 100;
    while ((c = getopt(argc, argv, "e:a:o:d:r:c:n:t:")) != -1) {
        switch (c) {
        case 'e':
            if (optarg)
                ellipse_file = optarg;
            break;
        case 'a':
            if (optarg)
                adj_file = optarg;
            break;
        case 'o':
            if (optarg)
                adj_ordered_file = optarg;
            break;
        case 'd':
            if (optarg)
                dis_file = optarg;
            break;
        case 'r':
            if (optarg)
                rows = atoi(optarg);
            break;
        case 'c':
            if (optarg)
                cols = atoi(optarg);
            break;
        case 'n':
            if (optarg)
                n_procs = atoi(optarg);
            break;
        }
    }

    if (n_procs == 1) {
        // test_matmul(rows);
    } else {
        // test_parallel_matmul(rows, &argc, &argv);
    }
    // build_overlap_matrix(ellipse_file, adj_file, rows);
    // find_components(adj_file);
    // build_distance_matrix(adj_ordered_file, dis_file);
    build_distance_matrix_parallel(adj_ordered_file, dis_file, &argc, &argv);

    return 0;
}
