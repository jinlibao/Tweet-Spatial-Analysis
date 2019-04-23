#include "include/ellipse.h"
#include "include/timer.h"
#include "include/tweets_spatial_analysis.h"
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[])
{
    string input_file("./data/adjacency_matrix.csv"), output_file("./data/distance_matrix.csv");

    int c;
    while ((c = getopt(argc, argv, "i:o:")) != -1) {
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

    double A[4][5] = {
        {0, 0, 2, 1, 0},
        {2, 1, 2, 1, 0},
        {2, 1, 2, 1, 0.5},
        {2, 1, 2, 1, 0.5}
    };

    double B[4][5] = {
        {0, 0, 3, 2, 0},
        {-1, 2, 3, 2, 0.5},
        {-1, 2, 3, 1.5, 0.5},
        {-1, 2, 3, 1.5, 1.6}
    };

    for (int i = 0; i < 4; ++i) {
        Ellipse e1 = Ellipse(A[i][0], A[i][1], A[i][2], A[i][3], A[i][4]);
        Ellipse e2 = Ellipse(B[i][0], B[i][1], B[i][2], B[i][3], B[i][4]);
        if (e1.overlap(e2)) {
            cout << "true" << endl;
        } else {
            cout << "false" << endl;
        }
    }
    
    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    cout << "Wall clock time elapsed: " << elapsed_time << " milliseconds" << endl;

    return 0;
}
