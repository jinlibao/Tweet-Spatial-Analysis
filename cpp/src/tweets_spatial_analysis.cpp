#include "include/tweets_spatial_analysis.h"
#include <cstdio>

void build_overlap_matrix(string ellipse_file, string adj_file, long rows)
{
    mat A;
    vector<int> remove_idx;
    vector<pair<unsigned int, unsigned long>> id;

    A.load(ellipse_file, csv_ascii);
    if (rows == 0) rows = A.n_rows;
    Mat<short> Adj(rows, rows, fill::zeros);
    cout << ellipse_file << ": " << rows << endl;
    long m = rows * rows / 2, mm = 0;
    for (int i = 0; i < rows; ++i) {
        if (A(i, 2) == 0 || A(i, 3) == 0) {
            remove_idx.push_back(i);
            //Adj.row(i) = -1 * ones<ivec>(rows).t();
            continue;
        }
        for (int j = 0; j < i; ++j) {
            Ellipse e1 = Ellipse(A(i, 0), A(i, 1), A(i, 2), A(i, 3), A(i, 4));
            Ellipse e2 = Ellipse(A(j, 0), A(j, 1), A(j, 2), A(j, 3), A(j, 4));
            if (e1.overlap(e2))
                Adj(i, j) = 1;
            else
                Adj(i, j) = 0;
            ++mm;
            if (mm % 1000000 == 0) {
                printf("%7.4f%% completed\n", (double)mm / m * 100);
            }
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
    char suffix[100] = "";
    sprintf(suffix, "_%ld.csv", rows);
    adj_file.replace(adj_file.end() - 4, adj_file.end(), suffix);
    string adj_id_file = adj_file;
    adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
    id_mat.save(adj_id_file, csv_ascii);
    Adj.save(adj_file, csv_ascii);
}

void find_components(string adj_file)
{
    Mat<short> A;
    Mat<long unsigned> id_mat;
    A.load(adj_file, csv_ascii);
    string adj_id_file = adj_file;
    adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
    id_mat.load(adj_id_file);
    vector<pair<int, vector<int>>> components = bfs(A);
    sort(components.begin(), components.end());
    reverse(components.begin(), components.end());

    int rows = A.n_rows;
    int cur = 0;
    vector<pair<long unsigned, long unsigned>> id;
    for (int i = 0; i < rows; i++) {
        id.push_back({(long unsigned)id_mat(i, 0), (long unsigned)id_mat(i, 1)});
    }
    Mat<short> B(rows, rows, fill::zeros);
    for (auto& c : components) {
        int m = c.first;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                B(cur + i, cur + j) = A(c.second[i], c.second[j]);
            }
            id[c.second[i]].first = cur + i;
        }
        cur += m;
    }
    sort(id.begin(), id.end());
    for (int i = 0; i < rows; ++i) {
        id_mat(i, 0) = i;
        id_mat(i, 1) = id[i].second;
    }

    // output
    adj_file.replace(adj_file.end() - 4, adj_file.end(), "_ordered.csv");
    adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_ordered.csv");
    id_mat.save(adj_id_file, csv_ascii);
    B.save(adj_file, csv_ascii);
}

vector<pair<int, vector<int>>> bfs(Mat<short>& A)
{
    int rows = A.n_rows;
    Col<short> visited(rows, fill::zeros);
    vector<pair<int, vector<int>>> components;
    for (int i = 0; i < rows; ++i) {
        if (!visited(i)) {
            visited(i) = 1;
            vector<int> component({i});
            queue<int> q;
            for (int j = 0; j < rows; ++j) {
                if (!visited(j) && A(i, j) > 0) {
                    visited(j) = 1;
                    q.push(j);
                    component.push_back(j);
                }
            }
            while (!q.empty()) {
                int j = q.front();
                q.pop();
                for (int k = 0; k < rows; ++k) {
                    if (!visited(k) && A(j, k) > 0) {
                        visited(k) = 1;
                        component.push_back(k);
                        q.push(k);
                    }
                }
            }
            components.push_back({(int)component.size(), component});
        }
    }
    return components;
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