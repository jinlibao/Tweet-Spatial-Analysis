#include "include/tweets_spatial_analysis.h"
#include "include/stopwatch.h"
#include <cstdio>
#include <stack>

template void build_overlap_matrix<short>(string input_file, string output_file, long rows = 0);
template void find_components<short>(string adj_file);
template vector<pair<int, vector<int>>> bfs<short>(Mat<short>& A);
template void build_distance_matrix<short>(string adj_file, string dis_file);
template void build_distance_matrix_parallel<short>(string adj_file, string dis_file, int *argc, char ***argv);
template Mat<short> APD<short>(const Mat<short>& A);
template Mat<short> APD_parallel<short>(const Mat<short> &A, int rank, int n_procs);
template Mat<short> APD_parallel_non_recursive<short>(const Mat<short> &A, int rank, int n_procs);
template void test_APD<short>(string mat_file);
template void test_APD_parallel<short>(string mat_file, int mode, int *argc, char ***argv);
template void test_APD_parallel_non_recursive<short>(string mat_file, int mode, int *argc, char ***argv);

template void build_overlap_matrix<int>(string input_file, string output_file, long rows = 0);
template void find_components<int>(string adj_file);
template vector<pair<int, vector<int>>> bfs<int>(Mat<int>& A);
template void build_distance_matrix<int>(string adj_file, string dis_file);
template void build_distance_matrix_parallel<int>(string adj_file, string dis_file, int *argc, char ***argv);
template Mat<int> APD<int>(const Mat<int>& A);
template Mat<int> APD_parallel<int>(const Mat<int> &A, int rank, int n_procs);
template Mat<int> APD_parallel_non_recursive<int>(const Mat<int> &A, int rank, int n_procs);
template void test_APD<int>(string mat_file);
template void test_APD_parallel<int>(string mat_file, int mode, int *argc, char ***argv);
template void test_APD_parallel_non_recursive<int>(string mat_file, int mode, int *argc, char ***argv);

template void build_overlap_matrix<long>(string input_file, string output_file, long rows = 0);
template void find_components<long>(string adj_file);
template vector<pair<int, vector<int>>> bfs<long>(Mat<long>& A);
template void build_distance_matrix<long>(string adj_file, string dis_file);
template void build_distance_matrix_parallel<long>(string adj_file, string dis_file, int *argc, char ***argv);
template Mat<long> APD<long>(const Mat<long>& A);
template Mat<long> APD_parallel<long>(const Mat<long> &A, int rank, int n_procs);
template Mat<long> APD_parallel_non_recursive<long>(const Mat<long> &A, int rank, int n_procs);
template void test_APD<long>(string mat_file);
template void test_APD_parallel<long>(string mat_file, int mode, int *argc, char ***argv);
template void test_APD_parallel_non_recursive<long>(string mat_file, int mode, int *argc, char ***argv);

template void build_overlap_matrix<float>(string input_file, string output_file, long rows = 0);
template void find_components<float>(string adj_file);
template vector<pair<int, vector<int>>> bfs<float>(Mat<float>& A);
template void build_distance_matrix<float>(string adj_file, string dis_file);
template void build_distance_matrix_parallel<float>(string adj_file, string dis_file, int *argc, char ***argv);
template Mat<float> APD<float>(const Mat<float>& A);
template Mat<float> APD_parallel<float>(const Mat<float> &A, int rank, int n_procs);
template Mat<float> APD_parallel_non_recursive<float>(const Mat<float> &A, int rank, int n_procs);
template void test_APD<float>(string mat_file);
template void test_APD_parallel<float>(string mat_file, int mode, int *argc, char ***argv);
template void test_APD_parallel_non_recursive<float>(string mat_file, int mode, int *argc, char ***argv);

template void build_overlap_matrix<double>(string input_file, string output_file, long rows = 0);
template void find_components<double>(string adj_file);
template vector<pair<int, vector<int>>> bfs<double>(Mat<double>& A);
template void build_distance_matrix<double>(string adj_file, string dis_file);
template void build_distance_matrix_parallel<double>(string adj_file, string dis_file, int *argc, char ***argv);
template Mat<double> APD<double>(const Mat<double>& A);
template Mat<double> APD_parallel<double>(const Mat<double> &A, int rank, int n_procs);
template Mat<double> APD_parallel_non_recursive<double>(const Mat<double> &A, int rank, int n_procs);
template void test_APD<double>(string mat_file);
template void test_APD_parallel<double>(string mat_file, int mode, int *argc, char ***argv);
template void test_APD_parallel_non_recursive<double>(string mat_file, int mode, int *argc, char ***argv);

template <class T>
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
    cout << "Wall clock time elapsed: " << elapsed_time << " ms" << endl;
}

template <class T>
void build_distance_matrix(string adj_file, string dis_file) {
    shino::precise_stopwatch stopwatch;
    Mat<T> A;
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
    Mat<int> AA = -1 * ones<Mat<int>>(rows, rows);
    for (int i = 0; i < (int)idx.size(); ++i) {
        int r1 = idx[i].first;
        int c1 = idx[i].first;
        int r2 = idx[i].second;
        int c2 = idx[i].second;

        cout << "Applying APD to Block matrix " << i << endl;
        Mat<T> D = APD<T>(A.submat(r1, c1, r2, c2));
        for (int j = r1; j < r2 + 1; ++j) {
            for (int k = c1; k < c2 + 1; ++k) {
                AA(j, k) = (int)D(j - r1, k - c1);
            }
        }
    }

    cout << "Write to " << dis_file << endl;
    AA.save(dis_file, csv_ascii);

    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    cout << "Wall clock time elapsed: " << elapsed_time << " ms" << endl;
}

template <class T>
void find_components(string adj_file) {
    shino::precise_stopwatch stopwatch;
    Mat<T> A;
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
                B(cur + i, cur + j) = (short)A(c.second[i], c.second[j]);
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
    cout << "Wall clock time elapsed: " << elapsed_time << " ms" << endl;
}

template <class T>
vector<pair<int, vector<int>>> bfs(Mat<T> &A) {
    int rows = A.n_rows;
    Col<T> visited(rows, fill::zeros);
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

template <class T>
Mat<T> APD(const Mat<T> &A) {
    long n = A.n_rows;
    Mat<T> Z = A * A;
    Mat<T> B(n, n, fill::zeros);
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

    Mat<T> D(n, n, fill::ones);
    D = -1 * D;
    if (cnt == (long)(n - 1) * n) {
        D = 2 * B - A;
        return D;
    }
    Mat<T> S = APD<T>(B);

    Mat<T> X = S * A;
    if (n > 1000) {
        //S.save("/gscratch/ljin1/data/twitter/csv/matmul_test/mat_S.csv", csv_ascii);
        //A.save("/gscratch/ljin1/data/twitter/csv/matmul_test/mat_A.csv", csv_ascii);
        //X.save("/gscratch/ljin1/data/twitter/csv/matmul_test/mat_X.csv", csv_ascii);
    }

    vec deg(n, fill::zeros);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            deg(i) += A(i, j);
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (X(i, j) >= S(i, j) * deg(j)) {
                D(i, j) = 2 * S(i, j);
            } else {
                D(i, j) = 2 * S(i, j) - 1;
            }
        }
    }
    return D;
}

template <class T>
Mat<T> APD_parallel(const Mat<T> &A, int node, int n_procs) {
    long n = A.n_rows;
    Mat<T> Z;

    if (n >= 100) {
        Z = parallel_matsq(A, A, node, n_procs);
    } else {
        Z = A * A;
    }

    Mat<T> B(n, n, fill::zeros);
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

    Mat<T> D(n, n, fill::ones);
    D = -1 * D;
    if (cnt == (n - 1) * n) {
        D = 2 * B - A;
        return D;
    }

    Mat<T> S = APD_parallel<T>(B, node, n_procs);

    Mat<T> X;
    if (n >= 100) {
        X = parallel_matmul(S, A, node, n_procs);
        //S.save("/gscratch/ljin1/data/twitter/csv/matmul_test/mat_S_" + to_string(node) + ".csv", csv_ascii);
        //A.save("/gscratch/ljin1/data/twitter/csv/matmul_test/mat_A_" + to_string(node) + ".csv", csv_ascii);
        //X.save("/gscratch/ljin1/data/twitter/csv/matmul_test/mat_X_" + to_string(node) + ".csv", csv_ascii);
    } else {
        X = S * A;
    }

    vec deg(n, fill::zeros);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            deg(i) += A(i, j);
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (X(i, j) >= S(i, j) * deg(j)) {
                D(i, j) = 2 * S(i, j);
            } else {
                D(i, j) = 2 * S(i, j) - 1;
            }
        }
    }
    return D;
}

template <class T>
Mat<T> APD_parallel_non_recursive(const Mat<T> &AA, int node, int n_procs) {
    long n = AA.n_rows;

    stack<Mat<T>> A_stack;
    Mat<T> A = AA;
    Mat<T> D;
    while (true) {
        Mat<T> Z;

        if (n >= 100) {
            Z = parallel_matsq(A, A, node, n_procs);
        } else {
            Z = A * A;
        }

        Mat<T> B(n, n, fill::zeros);
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

        if (cnt == (n - 1) * n) {
            D = 2 * B - A;
            break;
        }

        A_stack.push(A);
        A = B;
    }
    while (!A_stack.empty()) {
        A = A_stack.top();
        A_stack.pop();

        Mat<T> S = D;
        Mat<T> X;
        if (n >= 100) {
            X = parallel_matmul(S, A, node, n_procs);
            //S.save("/gscratch/ljin1/data/twitter/csv/matmul_test/mat_S_" + to_string(node) + ".csv", csv_ascii);
            //A.save("/gscratch/ljin1/data/twitter/csv/matmul_test/mat_A_" + to_string(node) + ".csv", csv_ascii);
            //X.save("/gscratch/ljin1/data/twitter/csv/matmul_test/mat_X_" + to_string(node) + ".csv", csv_ascii);
        } else {
            X = S * A;
        }

        vec deg(n, fill::zeros);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                deg(i) += A(i, j);
            }
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (X(i, j) >= S(i, j) * deg(j)) {
                    D(i, j) = 2 * S(i, j);
                } else {
                    D(i, j) = 2 * S(i, j) - 1;
                }
            }
        }
    }
    return D;
}

template <class T>
void build_distance_matrix_parallel(string adj_file, string dis_file, int *argc, char ***argv) {
    shino::precise_stopwatch stopwatch;
    int n_procs, node;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(argc, argv);
    MPI_Comm_size(comm, &n_procs);
    MPI_Comm_rank(comm, &node);

    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    if (node == 0) {
        printf("CPU %d: MPI initialized (%u ms)\n", node, elapsed_time);
    }

    Mat<T> A;
    Mat<long unsigned> id_mat;
    if (node == 0) {
        elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        printf("CPU %d: Reading from %s... (%u ms)\n", node, adj_file.c_str(), elapsed_time);
    }
    A.load(adj_file, csv_ascii);

    int rows = A.n_rows;
    string adj_id_file = adj_file;
    adj_id_file.replace(adj_id_file.end() - 4, adj_id_file.end(), "_id.csv");
    if (node == 0) {
        elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        printf("CPU %d: Reading from %s... (%u ms)\n", node, adj_id_file.c_str(), elapsed_time);
    }
    id_mat.load(adj_id_file);

    if (node == 0) {
        elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        printf("CPU %d: Unpacking id... (%u ms)\n", node, elapsed_time);
    }
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
    }

    Mat<int> AA = -1 * ones<Mat<int>>(rows, rows);
    for (int i = 0; i < (int)idx.size(); ++i) {
        int r1 = idx[i].first;
        int c1 = idx[i].first;
        int r2 = idx[i].second;
        int c2 = idx[i].second;

        if (node == 0) {
            elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
            printf("CPU %d: Applying APD to block matrix %d (%d-by-%d) (%u ms)\n", node, i, r2 - r1 + 1, c2 - c1 + 1, elapsed_time);
        }

        //Mat<T> D = APD_parallel<T>(A.submat(r1, c1, r2, c2), node, n_procs);
        Mat<T> D = APD_parallel_non_recursive<T>(A.submat(r1, c1, r2, c2), node, n_procs);
        if (node == 0) {
            //elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
            //printf("CPU %d: Assemble the distance matrix... (%u ms)\n", node, elapsed_time);
            for (int j = r1; j < r2 + 1; ++j) {
                for (int k = c1; k < c2 + 1; ++k) {
                    AA(j, k) = (int)D(j - r1, k - c1);
                }
            }
        }
    }

    if (node == 0) {
        //dis_file.replace(dis_file.end() - 4, dis_file.end(), "_parallel_recursive.csv");
        dis_file.replace(dis_file.end() - 4, dis_file.end(), "_parallel_non_recursive.csv");
        elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        printf("CPU %d: Writing to %s... (%u ms)\n", node, dis_file.c_str(), elapsed_time);
        AA.save(dis_file, csv_ascii);

        cout << "Wall clock time elapsed: " << elapsed_time << " ms" << endl;
    }
    MPI_Finalize();
}

template <class T>
void test_APD(string mat_file) {
    Mat<T> A;
    A.load(mat_file, csv_ascii);
    Mat<T> D = APD<T>(A);
    D.save(mat_file.replace(mat_file.end() - 4, mat_file.end(), "_D.csv"), csv_ascii);
}

template <class T>
void test_APD_parallel(string mat_file, int mode, int *argc, char ***argv) {
    int node, n_procs;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(argc, argv);
    MPI_Comm_rank(comm, &node);
    MPI_Comm_size(comm, &n_procs);

    Mat<T> A;
    A.load(mat_file, csv_ascii);
    Mat<T> D = APD_parallel<T>(A, node, n_procs);
    if (mode == 2) {
        D.save(mat_file.replace(mat_file.end() - 4, mat_file.end(), "_D_parallel_" + to_string(node) +  ".csv"), csv_ascii);
    } else {
        if (node == 0) {
            D.save(mat_file.replace(mat_file.end() - 4, mat_file.end(), "_D_parallel.csv"), csv_ascii);
        }
    }
    MPI_Finalize();
}

template <class T>
void test_APD_parallel_non_recursive(string mat_file, int mode, int *argc, char ***argv) {
    int node, n_procs;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(argc, argv);
    MPI_Comm_rank(comm, &node);
    MPI_Comm_size(comm, &n_procs);

    Mat<T> A;
    A.load(mat_file, csv_ascii);
    Mat<T> D = APD_parallel_non_recursive<T>(A, node, n_procs);
    if (mode == 2) {
        D.save(mat_file.replace(mat_file.end() - 4, mat_file.end(), "_D_parallel_" + to_string(node) +  ".csv"), csv_ascii);
    } else {
        if (node == 0) {
            D.save(mat_file.replace(mat_file.end() - 4, mat_file.end(), "_D_parallel.csv"), csv_ascii);
        }
    }
    MPI_Finalize();
}

