#include <armadillo>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace arma;

imat APD(imat A);
long long recursion_depth = 0;

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

    imat E;
    E.load(input_file, csv_ascii);
    imat F = APD(E);
    // cout << F << endl;
    F.save(output_file, csv_ascii);

    return 0;
}

imat APD(imat A)
{
    int n = A.n_rows;
    imat Z = A * A;
    imat B(n, n, fill::zeros);
    int cnt = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && (A(i, j) == 1 || Z(i, j) > 0)) {
                B(i, j) = 1;
                ++cnt;
            }
            else {
                B(i, j) = 0;
            }
        }
    }

    ++recursion_depth;
    cout << "Recursion depth: " << recursion_depth << endl;

    imat D(n, n, fill::ones);
    D = -1 * D;
    if (cnt == (n - 1) * n) {
        D = 2 * B - A;
        return D;
    }
    imat T = APD(B);

    cout << "Finished recursion " << --recursion_depth << endl;
    imat X = T * A;
    vec deg(n, fill::zeros);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            deg(i) += A(i, j);
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (X(i, j) >= T(i, j) * deg(j)) {
                D(i, j) = 2 * T(i, j);
            }
            else {
                D(i, j) = 2 * T(i, j) - 1;
            }
        }
    }
    return D;
}
