#include "include/matmul.h"

Mat<short> parallel_matmul(Mat<short> A, Mat<short> B, int *argc, char ***argv)
{
    int n_procs, rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;

    MPI_Init(argc, argv);
    MPI_Comm_size(comm, &n_procs);
    MPI_Comm_rank(comm, &rank);

    int a_rows = A.n_rows, a_cols = A.n_cols;
    int b_rows = B.n_rows, b_cols = B.n_cols;
    int rows = a_rows, cols = b_cols;
    pair<vector<pair<int, int>>, vector<pair<int, int>>> ranges = configure_cpu(rows, cols, n_procs);
    vector<pair<int, int>> row_idx = ranges.first;
    vector<pair<int, int>> col_idx = ranges.second;
    int n_procs_row = (int)row_idx.size();
    int n_procs_col = (int)col_idx.size();
    int r1, r2, c1, c2;
    Mat<short> C(rows, cols, fill::zeros);
    double start_time, end_time;

    if (rank == 0) {
        start_time = MPI_Wtime();
        vector<pair<short*, long>> arrays_row;
        vector<pair<short*, long>> arrays_col;

        for (int i = 0; i < n_procs_row; ++i) {
            r1 = row_idx[i].first;
            r2 = row_idx[i].second;
            short *a;
            mat2array(&a, A.rows(r1, r2));
            arrays_row.push_back({a, (r2 - r1 + 1) * a_cols});
        }

        for (int j = 0; j < n_procs_col; ++j) {
            c1 = col_idx[j].first;
            c2 = col_idx[j].second;
            short *b;
            mat2array(&b, B.cols(c1, c2));
            arrays_col.push_back({b, (c2 - c1 + 1) * b_rows});
        }

        //cout << "Starting calculation...\n";
        for (int i = 0; i < n_procs_row; ++i) {
            r1 = row_idx[i].first;
            r2 = row_idx[i].second;
            for (int j = 0; j < n_procs_col; ++j) {
                c1 = col_idx[j].first;
                c2 = col_idx[j].second;
                if (i + j == 0) continue;
                //printf("CPU %d: Send row block %d of A to CPU %d (%d, %d)\n", rank, i, i * n_procs_col + j, i, j);
                MPI_Send(arrays_row[i].first, arrays_row[i].second, MPI_SHORT, i * n_procs_col + j, 1, comm);
                //printf("CPU %d: Send col block %d of B to CPU %d (%d, %d)\n", rank, j, i * n_procs_col + j, i, j);
                MPI_Send(arrays_col[j].first, arrays_col[j].second, MPI_SHORT, i * n_procs_col + j, 2, comm);
            }
        }

        for (auto r : arrays_row) {
            delete[] r.first;
        }

        for (auto c : arrays_col) {
            delete[] c.first;
        }

        r1 = row_idx[0].first;
        r2 = row_idx[0].second;
        c1 = col_idx[0].first;
        c2 = col_idx[0].second;
        C.submat(r1, c1, r2, c2) = A.rows(r1, r2) * B.cols(c1, c2);

        // int n_recv = 1;
        // while (n_recv++ < n_procs) {
        //     int size = (row_idx[0].second - row_idx[0].first + 1) * (col_idx[0].second - col_idx[0].first + 1);
        //     short *c = new short[size];
        //     MPI_Recv(c, size, MPI_SHORT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
        //     int tag = status.MPI_TAG;
        //     int i = tag / n_procs_col;
        //     int j = tag % n_procs_col;
        //     r1 = row_idx[i].first;
        //     r2 = row_idx[i].second;
        //     c1 = col_idx[j].first;
        //     c2 = col_idx[j].second;
        //     Mat<short> C_sub(r2 - r1 + 1, c2 - c1 + 1, fill::zeros);
        //     array2mat(&c, C_sub);
        //     printf("CPU %d: Receive C(%d, %d) from CPU %d (%d, %d): %d-by-%d, %d\n", rank, i, j,  i * n_procs_col + j, i, j, (int)C_sub.n_rows, (int)C_sub.n_cols, (int)(C_sub.n_rows * C_sub.n_cols));
        //     C.submat(r1, c1, r2, c2) = C_sub;
        // }

        for (int i = 0; i < n_procs_row; ++i) {
            r1 = row_idx[i].first;
            r2 = row_idx[i].second;
            for (int j = 0; j < n_procs_col; ++j) {
                if (i + j == 0) continue;
                c1 = col_idx[j].first;
                c2 = col_idx[j].second;
                long size = (c2 - c1 + 1) * (r2 - r1 + 1);
                short *c = new short[size];
                Mat<short> C_sub(r2 - r1 + 1, c2 - c1 + 1, fill::zeros);
                //printf("CPU %d: Receive C(%d, %d) from CPU %d (%d, %d): %d-by-%d, %d\n", rank, i, j,  i * n_procs_col + j, i, j, (int)C_sub.n_rows, (int)C_sub.n_cols, (int)(C_sub.n_rows * C_sub.n_cols));
                MPI_Recv(c, size, MPI_SHORT, i * n_procs_col + j, i * n_procs_col + j, comm, &status);
                array2mat(&c, C_sub);
                C.submat(r1, c1, r2, c2) = C_sub;
            }
        }
        end_time = MPI_Wtime();
        cout << "Elapsed time for calculating matrix D: " << end_time - start_time << endl;
    }

    if (rank > 0) {
        int i = rank / n_procs_col;
        int j = rank % n_procs_col;
        r1 = row_idx[i].first;
        r2 = row_idx[i].second;
        c1 = col_idx[j].first;
        c2 = col_idx[j].second;
        long a_size = (r2 - r1 + 1) * cols;
        long b_size = (c2 - c1 + 1) * rows;
        long size = (r2 - r1 + 1) * (c2 - c1 + 1);
        short *a = new short[a_size];
        short *b = new short[b_size];
        //short *a = (short *)malloc(sizeof(short) * a_size);
        //short *b = (short *)malloc(sizeof(short) * b_size);
        //printf("CPU %d: Receive row block %d of A from CPU %d (%d, %d)\n", rank, i, 0, 0, 0);
        MPI_Recv(a, a_size, MPI_SHORT, 0, 1, comm, &status);
        //printf("CPU %d: Receive col block %d of B from CPU %d (%d, %d)\n", rank, j, 0, 0, 0);
        MPI_Recv(b, b_size, MPI_SHORT, 0, 2, comm, &status);
        Mat<short> TA(r2 - r1 + 1, cols, fill::zeros);
        Mat<short> TB(rows, c2 - c1 + 1, fill::zeros);
        array2mat(&a, TA);
        array2mat(&b, TB);
        Mat<short> TC = TA * TB;
        //cout << "CPU " << rank << ": \n" << TC << endl;
        //printf("CPU %d: C(%d, %d): %d-by-%d\n", rank, i, j, (int)TC.n_rows, (int)TC.n_cols);
        short *c;
        mat2array(&c, TC);

        //printf("CPU %d: Send C(%d, %d) to CPU 0 (0, 0)\n", rank, i, j);
        MPI_Send(c, size, MPI_SHORT, 0, rank, comm);
        delete[] c;
    }
    MPI_Finalize();
    return C;
}

void mat2array(short **A, const Mat<short>& B)
{
    int n_rows = B.n_rows, n_cols = B.n_cols;
    (*A) = new short[n_rows * n_cols];
    for (int i = 0; i < n_rows; ++i)
        for (int j = 0; j < n_cols; ++j)
            (*A)[i * n_cols + j] = B(i, j);
}

void array2mat(short **A, Mat<short>& B)
{
    int n_rows = B.n_rows, n_cols = B.n_cols;
    for (int i = 0; i < n_rows; ++i)
        for (int j = 0; j < n_cols; ++j)
            B(i, j) = (*A)[i * n_cols + j];
    delete[] (*A);
}

pair<vector<pair<int, int>>, vector<pair<int, int>>> configure_cpu(int rows, int cols, int n_procs, int n_procs_row, int n_procs_col)
{
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
            row_idx.push_back({t_row * (rows_per_proc + 1) + (i - t_row) * rows_per_proc, t_row * (rows_per_proc + 1) + (i - t_row + 1) * rows_per_proc - 1});
        }
    }
    for (int i = 0; i < n_procs_col; ++i) {
        if (i < t_col) {
            col_idx.push_back({i * (cols_per_proc + 1), (i + 1) * (cols_per_proc + 1) - 1});
        } else {
            col_idx.push_back({t_col * (cols_per_proc + 1) + (i - t_col) * cols_per_proc, t_col * (cols_per_proc + 1) + (i - t_col + 1) * cols_per_proc - 1});
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
