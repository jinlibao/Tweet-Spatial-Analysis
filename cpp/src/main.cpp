#include "include/timer.h"
#include "include/tweets_spatial_analysis.h"
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[])
{
    string input_file("./data/adjacency_matrix.csv"), output_file("./data/distance_matrix.csv");

    int c;
    while ((c = getopt(argc, argv, "i::o::")) != -1) {
        switch (c) {
        case 'i':
            if (optarg) input_file = optarg;
            break;
        case 'o':
            if (optarg) output_file = optarg;
            break;
        }
    }

    shino::precise_stopwatch stopwatch;

    imat E;
    E.load(input_file, csv_ascii);
    imat F = APD(E);
    F.save(output_file, csv_ascii);

    auto elapsed_time =
        stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    cout << "Wall clock time elapsed: " << elapsed_time << " milliseconds"
         << endl;

    return 0;
}
