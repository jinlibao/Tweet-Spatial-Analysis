#include "include/tweets_spatial_analysis.h"
#include "include/matmul.h"
#include "include/timer.h"
#include <cstdio>

void build_overlap_matrix(string ellipse_file, string adj_file, long rows) {
    shino::precise_stopwatch stopwatch;
    mat A;
    vector<int> remove_idx;
    vector<pair<unsigned int, long unsigned>> id;

    A.load(ellipse_file, csv_ascii);
    if (rows == 0)
        rows = A.n_rows;
    Mat<short> Adj(rows, rows, fill::zeros);
    cout << ellipse_file << ": " << rows << endl;
    long m = rows * rows / 2, mm = 0;
    for (int i = 0; i < rows; ++i) {
        if (A(i, 2) == 0 || A(i, 3) == 0) {
            remove_idx.push_back(i);
            // Adj.row(i) = -1 * ones<ivec>(rows).t();
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
        id.push_back({(unsigned int)i, (long unsigned)A(i, 5)});
    }
    Adj = Adj + Adj.t();

    // remove rows and cols that contain -1
    reverse(remove_idx.begin(), remove_idx.end());
    for (auto &i : remove_idx) {
        Adj.shed_row(i);
        Adj.shed_col(i);
    }
    rows = Adj.n_rows;
    Mat<long unsigned> id_mat(rows, 2, fill::zeros);
    for (int i = 0; i < rows; ++i) {
        id_mat(i, 0) = i;
        id_mat(i, 1) = id[i].second;
    }

    // output
    // char suffix[100] = "";
    // sprintf(suffix, "_%ld.csv", rows);
    // adj_file.replace(adj_file.end() - 4, adj_file.end(), suffix);
    string adj_id_file = adj_file;
    adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
    id_mat.save(adj_id_file, csv_ascii);
    Adj.save(adj_file, csv_ascii);

    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    cout << "Wall clock time elapsed: " << elapsed_time << " milliseconds" << endl;
}

void build_distance_matrix(string adj_file, string dis_file) {
    shino::precise_stopwatch stopwatch;
    Mat<short> A;
    Mat<long unsigned> id_mat;
    cout << "Read from " << adj_file << endl;
    A.load(adj_file, csv_ascii);
    int rows = A.n_rows;
    string adj_id_file = adj_file;
    adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
    cout << "Read from " << adj_id_file << endl;
    id_mat.load(adj_id_file);

    vector<pair<unsigned int, long unsigned>> id;
    for (int i = 0; i < rows; i++) {
        id.push_back({(long unsigned)id_mat(i, 0), (long unsigned)id_mat(i, 1)});
    }

    vector<pair<int, int>> idx;
    for (int i = 0; i < rows;) {
        int start = i;
        while (i < rows && id[start].first == id[i].first) {
            ++i;
        }
        int end = i - 1;
        idx.push_back({start, end});
        // cout << start << ' ' << end << endl;
    }
    Mat<short> AA = -1 * ones<Mat<short>>(rows, rows);
    for (int i = 0; i < (int)idx.size(); ++i) {
        int r1 = idx[i].first;
        int c1 = idx[i].first;
        int r2 = idx[i].second;
        int c2 = idx[i].second;

        cout << "Applying APD to Block matrix " << i << endl;
        Mat<short> D;
        D = APD(A.submat(r1, c1, r2, c2));
        for (int j = r1; j < r2 + 1; ++j) {
            for (int k = c1; k < c2 + 1; ++k) {
                AA(j, k) = D(j - r1, k - c1);
            }
        }
    }

    cout << "Write to " << dis_file << endl;
    AA.save(dis_file, csv_ascii);

    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    cout << "Wall clock time elapsed: " << elapsed_time << " milliseconds" << endl;
}

void find_components(string adj_file) {
    shino::precise_stopwatch stopwatch;
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
    vector<vector<long unsigned>> id;
    for (int i = 0; i < rows; i++) {
        id.push_back({(long unsigned)id_mat(i, 0), (long unsigned)id_mat(i, 1)});
    }
    Mat<short> B(rows, rows, fill::zeros);
    int k = 0;
    for (auto &c : components) {
        int m = c.first;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                B(cur + i, cur + j) = A(c.second[i], c.second[j]);
            }
            id[c.second[i]][0] = cur + i;
            id[c.second[i]].insert(id[c.second[i]].begin() + 1, k);
        }
        cur += m;
        ++k;
    }
    sort(id.begin(), id.end());
    for (int i = 0; i < rows; ++i) {
        id_mat(i, 0) = id[i][1];
        id_mat(i, 1) = id[i][2];
    }

    // output
    adj_file.replace(adj_file.end() - 4, adj_file.end(), "_ordered.csv");
    adj_id_file = adj_file;
    adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
    id_mat.save(adj_id_file, csv_ascii);
    B.save(adj_file, csv_ascii);

    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    cout << "Wall clock time elapsed: " << elapsed_time << " milliseconds" << endl;
}

vector<pair<int, vector<int>>> bfs(Mat<short> &A) {
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

Mat<short> APD(const Mat<short> &A) {
    int n = A.n_rows;
    Mat<short> Z = A * A;
    Mat<short> B(n, n, fill::zeros);
    long cnt = 0;
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

    Mat<short> D(n, n, fill::ones);
    D = -1 * D;
    if (cnt == (long)(n - 1) * n) {
        D = 2 * B - A;
        return D;
    }
    Mat<short> T = APD(B);

    Mat<short> X = T * A;
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

Mat<short> APD_parallel(const Mat<short> &A, int rank, int n_procs) {
    int n = A.n_rows;
    Mat<short> Z;

    if (n > 100) {
        Z = parallel_matmul(A, A, rank, n_procs);
    } else {
        Z = A * A;
    }
    Mat<short> B(n, n, fill::zeros);
    long cnt = 0;
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

    Mat<short> D(n, n, fill::ones);
    D = -1 * D;
    if (cnt == (long)(n - 1) * n) {
        D = 2 * B - A;
        return D;
    }
    Mat<short> T = APD_parallel(B, rank, n_procs);

    Mat<short> X;
    if (n > 100) {
        X = parallel_matmul(T, A, rank, n_procs);
    } else {
        X = T * A;
    }
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

void build_distance_matrix_parallel(string adj_file, string dis_file, int *argc, char ***argv) {
    int n_procs, rank;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(argc, argv);
    MPI_Comm_size(comm, &n_procs);
    MPI_Comm_rank(comm, &rank);

    shino::precise_stopwatch stopwatch;

    Mat<short> A;
    Mat<long unsigned> id_mat;
    if (rank == 0) {
        cout << "Read from " << adj_file << endl;
    }
    A.load(adj_file, csv_ascii);
    int rows = A.n_rows;
    string adj_id_file = adj_file;
    adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
    if (rank == 0) {
        cout << "Read from " << adj_id_file << endl;
    }
    id_mat.load(adj_id_file);

    vector<pair<unsigned int, long unsigned>> id;
    for (int i = 0; i < rows; i++) {
        id.push_back({(long unsigned)id_mat(i, 0), (long unsigned)id_mat(i, 1)});
    }

    vector<pair<int, int>> idx;
    for (int i = 0; i < rows;) {
        int start = i;
        while (i < rows && id[start].first == id[i].first) {
            ++i;
        }
        int end = i - 1;
        idx.push_back({start, end});
        // cout << start << ' ' << end << endl;
    }

    Mat<short> AA = -1 * ones<Mat<short>>(rows, rows);
    for (int i = 0; i < (int)idx.size(); ++i) {
        int r1 = idx[i].first;
        int c1 = idx[i].first;
        int r2 = idx[i].second;
        int c2 = idx[i].second;

        if (rank == 0) {
            cout << "Applying APD to Block matrix " << i << endl;
        }

        Mat<short> D;
        D = APD_parallel(A.submat(r1, c1, r2, c2), rank, n_procs);
        if (rank == 0) {
            for (int j = r1; j < r2 + 1; ++j) {
                for (int k = c1; k < c2 + 1; ++k) {
                    AA(j, k) = D(j - r1, k - c1);
                }
            }
        }
    }

    if (rank == 0) {
        cout << "Write to " << dis_file << endl;
        AA.save(dis_file, csv_ascii);

        auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        cout << "Wall clock time elapsed: " << elapsed_time << " milliseconds" << endl;
    }

    MPI_Finalize();
}
