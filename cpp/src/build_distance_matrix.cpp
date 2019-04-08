#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

mat APD(mat A);

int main()
{
    mat E;
    E.load("./data/adjacency_matrix.csv", csv_ascii);
    //cout << E << endl;
    mat F = APD(E);
    cout << F << endl;
    F.save("./data/distance_matrix.csv", csv_ascii);

    return 0;
}

mat APD(mat A)
{
    int n = A.n_rows;
    mat Z = A * A;
    mat B(n, n, fill::zeros);
    int cnt = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && (A(i, j) == 1 || Z(i, j) > 0)) {
                B(i, j) = 1;
                ++cnt;
            } else {
                B(i, j) = 0;
            }
        }
    }

    mat D(n, n, fill::ones);
    D = -1 * D;
    if (cnt == (n - 1) * n) {
        D = 2 * B - A;
        return D;
    }
    mat T = APD(B);
    mat X = T * A;
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
