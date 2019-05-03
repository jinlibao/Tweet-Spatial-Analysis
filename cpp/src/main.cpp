#include "include/matmul.h"
#include "include/timer.h"
#include "include/tweets_spatial_analysis.h"
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[]) {
    string mat_A("./data/data_100_A.csv");
    string mat_B("./data/data_100_B.csv");
    string mat_C("./data/data_100_C.csv");
    string ellipse_file("../data/tweets_median_working_filtered.csv");
    string adj_file("./data/tweets_median_working_adjacency_matrix.csv");
    string adj_ordered_file("./data/tweets_median_working_adjacency_matrix_ordered.csv");
    string dis_file("./data/tweets_median_working_distance_matrix.csv");

    int c, rows, cols, n_procs, job;
    rows = cols = 0;
    n_procs = 1;
    job = 0;
    while ((c = getopt(argc, argv, "A:B:C:e:a:o:d:r:c:n:t:j:")) != -1) {
        switch (c) {
        case 'A':
            if (optarg)
                mat_A = optarg;
            break;
        case 'B':
            if (optarg)
                mat_B = optarg;
            break;
        case 'C':
            if (optarg)
                mat_C = optarg;
            break;
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
        case 'j':
            if (optarg)
                job = atoi(optarg);
            break;
        }
    }

    if (job == 0) {
        if (n_procs == 1) {
            test_APD(mat_A);
            // test_matmul(rows);
            // test_matmul(rows, 1, mat_A, mat_B, mat_C);
            // test_matsq(rows);
            // test_matsq(rows, 1, mat_A, mat_C);
        } else {
            test_APD_parallel(mat_A, 1, &argc, &argv);
            // test_parallel_matmul(rows, &argc, &argv);
            // test_parallel_matmul(rows, &argc, &argv, 1, mat_A, mat_B, mat_C);
            // test_parallel_matsq(rows, &argc, &argv);
            // test_parallel_matsq(rows, &argc, &argv, 1, mat_A, mat_C);
        }
    }
    
    if (job == 1) {
        build_overlap_matrix(ellipse_file, adj_file, rows);
        find_components(adj_file);
    }
    if (job == 2) {
        build_distance_matrix(adj_ordered_file, dis_file);
    }
    if (job == 3) {
        build_distance_matrix_parallel(adj_ordered_file, dis_file, &argc, &argv);
    }

    return 0;
}
