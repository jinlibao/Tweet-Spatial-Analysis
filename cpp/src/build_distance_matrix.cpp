#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

mat APD(mat A);

int main()
{
    mat B(3, 3, fill::zeros);
    B(0, 1) = 1;
    B(1, 0) = 1;
    B(1, 2) = 1;
    B(2, 1) = 1;

    cout << 2 * B << endl;
    cout << B.n_rows << ' ' << B.n_cols << endl;
    mat D = APD(B);
    cout << D << endl;

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
            cnt += (i != j);
            if (i != j && (A(i, j) == 1 || Z(i, j) > 0)) {
                B(i, j) = 1;
            }
            else {
                B(i, j) = 0;
            }
        }
    }

    mat D(n, n, fill::zeros);
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
            } else {
                D(i, j) = 2 * T(i, j) - 1;
            }
        }
    }
    return D;
}
