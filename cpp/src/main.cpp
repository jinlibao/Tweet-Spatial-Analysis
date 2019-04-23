#include "include/timer.h"
#include "include/tweets_spatial_analysis.h"
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[])
{
    string ellipse_file("../data/tweets_median_working_filtered.csv");
    string input_file("./data/tweets_median_working_adjacency_matrix.csv");
    string output_file("./data/tweets_median_working_distance_matrix.csv");

    int c, rows = 0;
    while ((c = getopt(argc, argv, "e:i:o:r:")) != -1) {
        switch (c) {
        case 'e':
            if (optarg) ellipse_file = optarg;
            break;
        case 'i':
            if (optarg) input_file = optarg;
            break;
        case 'o':
            if (optarg) output_file = optarg;
            break;
        case 'r':
            if (optarg) rows = atoi(optarg);
            break;
        }
    }

    shino::precise_stopwatch stopwatch;

    build_overlap_matrix(ellipse_file, input_file, rows);

    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    cout << "Wall clock time elapsed: " << elapsed_time << " milliseconds" << endl;

    return 0;
}


