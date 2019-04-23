#include "include/tweets_spatial_analysis.h"

void build_overlap_matrix(string input_file, string output_file, int rows)
{
    mat A;
    imat Adj(rows, rows, fill::zeros);
    vector<int> remove_idx;
    vector<pair<unsigned int, unsigned long>> id;

    A.load(input_file, csv_ascii);
    if (rows == 0) rows = A.n_rows;
    for (int i = 0; i < rows; ++i) {
        if (A(i, 2) == 0 || A(i, 3) == 0) {
            remove_idx.push_back(i);
            //Adj.row(i) = -1 * ones<ivec>(rows).t();
            continue;
        }
        for (int j = 0; j < i; ++j) {
            Ellipse e1 = Ellipse(A(i, 0), A(i, 1), A(i, 2), A(i, 3), A(i, 4));
            Ellipse e2 = Ellipse(A(j, 0), A(j, 1), A(j, 2), A(j, 3), A(j, 4));
            Adj(i, j) = e1.overlap(e2);
        }
        id.push_back({(unsigned int)i, (unsigned long)A(i, 5)});
    }
    Adj = Adj + Adj.t();

    // remove rows and cols that contain -1
    reverse(remove_idx.begin(), remove_idx.end());
    for (auto& i : remove_idx) {
        Adj.shed_row(i);
        Adj.shed_col(i);
    }
    rows = Adj.n_rows;
    Mat<unsigned long> id_mat(rows, 2, fill::zeros);
    for (int i = 0; i < rows; ++i) {
        id_mat(i, 0) = i;
        id_mat(i, 1) = id[i].second;
    }

    // output
    string output_file_id = output_file;
    output_file_id.replace(output_file_id.end() - 4, output_file_id.end(), "_id.csv");
    id_mat.save(output_file_id, csv_ascii);
    Adj.save(output_file, csv_ascii);
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
