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
    string outlier_file("");
    string suc_file("./data/tweets_median_working_successor_matrix.csv");
    string id_file("./data/tweets_median_working_adjacency_matrix_ordered_id.csv");

    long unsigned id_from = 379005817LU, id_to = 134577044LU;

    int c, rows, cols, job;
    rows = cols = 0;
    job = 0;
    while ((c = getopt(argc, argv, "A:B:C:O:e:a:o:d:r:c:j:s:i:f:t:")) != -1) {
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
        case 'O':
            if (optarg)
                outlier_file = optarg;
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
        case 'j':
            if (optarg)
                job = atoi(optarg);
            break;
        case 's':
            if (optarg)
                suc_file = optarg;
            break;
        case 'i':
            if (optarg)
                id_file = optarg;
            break;
        case 'f':
            if (optarg)
                id_from = (long unsigned)atol(optarg);
            break;
        case 't':
            if (optarg)
                id_to = (long unsigned)atol(optarg);
            break;
        }
    }

    int node, n_procs;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &node);
    MPI_Comm_size(comm, &n_procs);

    if (job == 0) {
        APSP<float>(rows, ellipse_file, adj_file, adj_ordered_file, dis_file, outlier_file, node, n_procs);
    }
    if (job == 1) {
        build_adjacency_matrix<float>(ellipse_file, adj_file, rows);
    }
    if (job == 2) {
        find_components<float>(adj_file, outlier_file);
    }
    if (job == 3) {
        build_distance_matrix<float>(adj_ordered_file, dis_file, node, n_procs);
    }
    if (job == 4) {
        build_successor_matrix<float>(adj_ordered_file, dis_file, node, n_procs);
    }
    if (job == 5) {
        test_APD_recursive<float>(mat_A, 1, node, n_procs);
    }
    if (job == 6) {
        test_APD<float>(mat_A, 1, node, n_procs);
    }
    if (job == 7) {
        test_matrix_multiplication<float>(rows, node, n_procs);
    }
    if (job == 8) {
        test_matrix_multiplication<float>(rows, node, n_procs, 1, mat_A, mat_B, mat_C);
    }
    if (job == 9) {
        test_matrix_square<float>(rows, node, n_procs);
    }
    if (job == 10) {
        test_matrix_square<float>(rows, node, n_procs, 1, mat_A, mat_C);
    }
    if (job == 11) {
        get_shortest_path_by_id(suc_file, id_file, id_from, id_to);
    }

    MPI_Finalize();
    return 0;
}
