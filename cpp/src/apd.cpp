#include "include/tweets_spatial_analysis.h"

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

    imat D(n, n, fill::ones);
    D = -1 * D;
    if (cnt == (n - 1) * n) {
        D = 2 * B - A;
        return D;
    }
    imat T = APD(B);

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
