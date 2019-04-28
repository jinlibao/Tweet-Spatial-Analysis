#include "include/matmul.h"
#include "include/timer.h"

Mat<short> parallel_matmul(const Mat<short> &A, const Mat<short> &B, int rank, int n_procs) {
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;
    long rows = A.n_rows, cols = B.n_cols;
    pair<vector<pair<int, int>>, vector<pair<int, int>>> ranges = configure_cpu(rows, cols, n_procs);
    vector<pair<int, int>> row_idx = ranges.first;
    vector<pair<int, int>> col_idx = ranges.second;
    int n_procs_col = (int)col_idx.size();
    int r1, r2, c1, c2;
    Mat<short> C(rows, cols, fill::zeros);
    shino::precise_stopwatch stopwatch;

    if (rank == 0) {
        r1 = row_idx[0].first;
        r2 = row_idx[0].second;
        c1 = col_idx[0].first;
        c2 = col_idx[0].second;
        C.submat(r1, c1, r2, c2) = A.rows(r1, r2) * B.cols(c1, c2);
        int size = (r2 - r1 + 1) * (c2 - c1 + 1);

        for (int k = 0; k < n_procs - 1; ++k) {
            short *c;
            c = new short[size];
            // MPI_Alloc_mem(sizeof(short)*size, MPI_INFO_NULL, &c);
            MPI_Recv(c, size, MPI_SHORT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
            int tag = status.MPI_TAG;
            int i = tag / n_procs_col;
            int j = tag % n_procs_col;
            r1 = row_idx[i].first;
            r2 = row_idx[i].second;
            c1 = col_idx[j].first;
            c2 = col_idx[j].second;
            Mat<short> C_sub(r2 - r1 + 1, c2 - c1 + 1, fill::zeros);
            array2mat(&c, C_sub);
            // printf("CPU %d: Receive C(%d, %d) from CPU %d (%d, %d): %d-by-%d, %d\n", rank, i, j,  i * n_procs_col + j, i, j,
            // (int)C_sub.n_rows, (int)C_sub.n_cols, (int)(C_sub.n_rows * C_sub.n_cols));
            C.submat(r1, c1, r2, c2) = C_sub;
        }

        short *c = new short[rows * cols];
        mat2array(&c, C);
        MPI_Barrier(comm);
        MPI_Bcast(c, rows * cols, MPI_SHORT, 0, comm);

        auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        cout << "Number of CPUs: " << n_procs << ", rows: " << rows;
        cout << ", wall clock time elapsed to perform parallel_matmul: " << elapsed_time << " milliseconds" << endl;
    }

    if (rank > 0) {
        int i = rank / n_procs_col;
        int j = rank % n_procs_col;
        r1 = row_idx[i].first;
        r2 = row_idx[i].second;
        c1 = col_idx[j].first;
        c2 = col_idx[j].second;
        long size = (r2 - r1 + 1) * (c2 - c1 + 1);
        Mat<short> TC = A.rows(r1, r2) * B.cols(c1, c2);
        short *c;
        mat2array(&c, TC);
        // printf("CPU %d: Send C(%d, %d) to CPU 0 (0, 0)\n", rank, i, j);
        MPI_Send(c, size, MPI_SHORT, 0, rank, comm);
        delete[] c;

        c = new short[rows * cols];
        MPI_Barrier(comm);
        MPI_Bcast(c, rows * cols, MPI_SHORT, 0, comm);
        array2mat(&c, C);
    }

    return C;
}

void mat2array(short **A, const Mat<short> &B) {
    int n_rows = B.n_rows, n_cols = B.n_cols;
    (*A) = new short[n_rows * n_cols];
    for (int i = 0; i < n_rows; ++i)
        for (int j = 0; j < n_cols; ++j)
            (*A)[i * n_cols + j] = B(i, j);
}

void array2mat(short **A, Mat<short> &B) {
    int n_rows = B.n_rows, n_cols = B.n_cols;
    for (int i = 0; i < n_rows; ++i)
        for (int j = 0; j < n_cols; ++j)
            B(i, j) = (*A)[i * n_cols + j];
    delete[](*A);
}

pair<vector<pair<int, int>>, vector<pair<int, int>>> configure_cpu(int rows, int cols, int n_procs, int n_procs_row, int n_procs_col) {
    // n_procs_row x n_procs_col grid processors
    if (n_procs_row == 0 || n_procs_col == 0) {
        n_procs_row = floor(sqrt(n_procs));
        n_procs_col = floor(sqrt(n_procs));
    }
    int rows_per_proc = floor((double)rows / n_procs_row);
    int cols_per_proc = floor((double)cols / n_procs_col);
    int t_row = rows - rows_per_proc * n_procs_row;
    int t_col = cols - cols_per_proc * n_procs_col;
    vector<pair<int, int>> row_idx;
    vector<pair<int, int>> col_idx;
    for (int i = 0; i < n_procs_row; ++i) {
        if (i < t_row) {
            row_idx.push_back({i * (rows_per_proc + 1), (i + 1) * (rows_per_proc + 1) - 1});
        } else {
            row_idx.push_back({t_row * (rows_per_proc + 1) + (i - t_row) * rows_per_proc,
                               t_row * (rows_per_proc + 1) + (i - t_row + 1) * rows_per_proc - 1});
        }
    }
    for (int i = 0; i < n_procs_col; ++i) {
        if (i < t_col) {
            col_idx.push_back({i * (cols_per_proc + 1), (i + 1) * (cols_per_proc + 1) - 1});
        } else {
            col_idx.push_back({t_col * (cols_per_proc + 1) + (i - t_col) * cols_per_proc,
                               t_col * (cols_per_proc + 1) + (i - t_col + 1) * cols_per_proc - 1});
        }
    }
    // printf("rows: %d, cols: %d\n", rows, cols);
    // printf("n_procs: %d, n_procs_row: %d, n_procs_col: %d\n", n_procs, n_procs_row, n_procs_col);
    // printf("rows_per_proc: %d, cols_per_proc: %d\n", rows_per_proc, cols_per_proc);
    // printf("t_row: %d, t_col: %d\n", t_row, t_col);
    // for (int i = 0; i < (int)row_idx.size(); ++i)
    //     printf("%d: %d - %d\n", i, row_idx[i].first, row_idx[i].second);
    // for (int i = 0; i < (int)col_idx.size(); ++i)
    //     printf("%d: %d - %d\n", i, col_idx[i].first, col_idx[i].second);
    return {row_idx, col_idx};
}

void test_parallel_matmul(int rows, int *argc, char ***argv) {
    int n_procs, rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    shino::precise_stopwatch stopwatch;

    MPI_Init(argc, argv);
    MPI_Comm_size(comm, &n_procs);
    MPI_Comm_rank(comm, &rank);
    Mat<short> A(rows, rows, fill::randn);
    Mat<short> B(rows, rows, fill::randn);
    Mat<short> C = parallel_matmul(A, B, rank, n_procs);
    if (rank == 0) {
        auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        cout << "Number of CPUs: " << n_procs << ", rows: " << rows << endl;
        cout << "Wall clock time elapsed: " << elapsed_time << " milliseconds" << endl;
    }
    MPI_Finalize();
}

void test_matmul(int rows) {
    Mat<short> A(rows, rows, fill::randn);
    Mat<short> B(rows, rows, fill::randn);
    shino::precise_stopwatch stopwatch;
    Mat<short> C = A * B;
    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    cout << "Number of CPUs: 1, rows: " << rows << endl;
    cout << "Wall clock time elapsed: " << elapsed_time << " milliseconds" << endl;
}
