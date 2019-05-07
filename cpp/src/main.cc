#include "include/tweets_spatial_analysis.h"
#include <cstdlib>

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
        build_overlap_matrix<float>(ellipse_file, adj_file, rows);
        find_components<float>(adj_file);
    }
    if (job == 1) {
        if (n_procs == 1) {
            build_distance_matrix<float>(adj_ordered_file, dis_file);
        } else {
            build_distance_matrix_parallel<float>(adj_ordered_file, dis_file, &argc, &argv);
        }
    }
    if (job == 2) {
        if (n_procs == 1) {
            test_APD<float>(mat_A);
        } else {
            // test_APD_parallel<float>(mat_A, 1, &argc, &argv);
            test_APD_parallel_non_recursive<float>(mat_A, 1, &argc, &argv);
        }
    }
    if (job == 3) {
        if (n_procs == 1) {
            test_matmul<float>(rows);
        } else {
            test_parallel_matmul<float>(rows, &argc, &argv);
        }
    }
    if (job == 4) {
        if (n_procs == 1) {
            test_matmul<float>(rows, 1, mat_A, mat_B, mat_C);
        } else {
            test_parallel_matmul<float>(rows, &argc, &argv, 1, mat_A, mat_B, mat_C);
        }
    }
    if (job == 5) {
        if (n_procs == 1) {
            test_matsq<float>(rows);
        } else {
            test_parallel_matsq<float>(rows, &argc, &argv);
        }
    }
    if (job == 6) {
        if (n_procs == 1) {
            test_matsq<float>(rows, 1, mat_A, mat_C);
        } else {
            test_parallel_matsq<float>(rows, &argc, &argv, 1, mat_A, mat_C);
        }
    }
    if (job == 7) {
        if (n_procs == 1) {
            build_overlap_matrix<float>(ellipse_file, adj_file, rows);
            find_components<float>(adj_file);
            build_distance_matrix<float>(adj_ordered_file, dis_file);
        } else {
            build_distance_matrix_parallel<float>(adj_ordered_file, dis_file, &argc, &argv);
        }
    }
    if (job == 8) {
        if (n_procs == 1) {
            build_overlap_matrix<short>(ellipse_file, adj_file, rows);
            find_components<short>(adj_file);
            build_distance_matrix<short>(adj_ordered_file, dis_file);
        } else {
            build_distance_matrix_parallel<short>(adj_ordered_file, dis_file, &argc, &argv);
        }
    }

    return 0;
}
