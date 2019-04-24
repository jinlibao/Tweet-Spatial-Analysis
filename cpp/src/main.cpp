#include "include/timer.h"
#include "include/tweets_spatial_analysis.h"
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[])
{
    string ellipse_file("../data/tweets_median_working_filtered.csv");
    string adj_file("./data/tweets_median_working_adjacency_matrix.csv");
    string adj_ordered_file("./data/tweets_median_working_adjacency_matrix_ordered.csv");
    string dis_file("./data/tweets_median_working_distance_matrix.csv");

    int c, rows = 0;
    while ((c = getopt(argc, argv, "e:a:o:d:r:")) != -1) {
        switch (c) {
        case 'e':
            if (optarg) ellipse_file = optarg;
            break;
        case 'a':
            if (optarg) adj_file = optarg;
            break;
        case 'o':
            if (optarg) adj_ordered_file = optarg;
            break;
        case 'd':
            if (optarg) dis_file = optarg;
            break;
        case 'r':
            if (optarg) rows = atoi(optarg);
            break;
        }
    }

    shino::precise_stopwatch stopwatch;

    //test_matmul();
    build_overlap_matrix(ellipse_file, adj_file, rows);
    find_components(adj_file);
    build_distance_matrix(adj_ordered_file, dis_file);

    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    cout << "Wall clock time elapsed: " << elapsed_time << " milliseconds" << endl;

    return 0;
}


