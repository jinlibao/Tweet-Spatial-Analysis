#include "include/matmul.h"
#include "include/stopwatch.h"

template Mat<short> matrix_multiplication<short>(const Mat<short> &A, const Mat<short> &B, int node, int n_procs);
template Mat<short> matrix_square<short>(const Mat<short> &A, int node, int n_procs);
template void mat2array<short>(short **A, const Mat<short> &B);
template void array2mat<short>(short **A, Mat<short> &B);
template void mat2array_square<short>(short **A, const Mat<short> &B);
template void array2mat_square<short>(short **A, Mat<short> &B);
template void get_mpi_type<short>(MPI_Datatype *mpi_type);
template void test_matrix_multiplication<short>(int rows, int node, int n_procs, int mode = 0, string mat_A = "", string mat_B = "", string mat_C = "");
template void test_matrix_square<short>(int rows, int node, int n_procs, int mode = 0, string mat_A = "", string mat_C = "");

template Mat<int> matrix_multiplication<int>(const Mat<int> &A, const Mat<int> &B, int node, int n_procs);
template Mat<int> matrix_square<int>(const Mat<int> &A, int node, int n_procs);
template void mat2array<int>(int **A, const Mat<int> &B);
template void array2mat<int>(int **A, Mat<int> &B);
template void mat2array_square<int>(int **A, const Mat<int> &B);
template void array2mat_square<int>(int **A, Mat<int> &B);
template void get_mpi_type<int>(MPI_Datatype *mpi_type);
template void test_matrix_multiplication<int>(int rows, int node, int n_procs, int mode = 0, string mat_A = "", string mat_B = "", string mat_C = "");
template void test_matrix_square<int>(int rows, int node, int n_procs, int mode = 0, string mat_A = "", string mat_C = "");

template Mat<long> matrix_multiplication<long>(const Mat<long> &A, const Mat<long> &B, int node, int n_procs);
template Mat<long> matrix_square<long>(const Mat<long> &A, int node, int n_procs);
template void mat2array<long>(long **A, const Mat<long> &B);
template void array2mat<long>(long **A, Mat<long> &B);
template void mat2array_square<long>(long **A, const Mat<long> &B);
template void array2mat_square<long>(long **A, Mat<long> &B);
template void get_mpi_type<long>(MPI_Datatype *mpi_type);
template void test_matrix_multiplication<long>(int rows, int node, int n_procs, int mode = 0, string mat_A = "", string mat_B = "", string mat_C = "");
template void test_matrix_square<long>(int rows, int node, int n_procs, int mode = 0, string mat_A = "", string mat_C = "");

template Mat<float> matrix_multiplication<float>(const Mat<float> &A, const Mat<float> &B, int node, int n_procs);
template Mat<float> matrix_square<float>(const Mat<float> &A, int node, int n_procs);
template void mat2array<float>(float **A, const Mat<float> &B);
template void array2mat<float>(float **A, Mat<float> &B);
template void mat2array_square<float>(float **A, const Mat<float> &B);
template void array2mat_square<float>(float **A, Mat<float> &B);
template void get_mpi_type<float>(MPI_Datatype *mpi_type);
template void test_matrix_multiplication<float>(int rows, int node, int n_procs, int mode = 0, string mat_A = "", string mat_B = "", string mat_C = "");
template void test_matrix_square<float>(int rows, int node, int n_procs, int mode = 0, string mat_A = "", string mat_C = "");

template Mat<double> matrix_multiplication<double>(const Mat<double> &A, const Mat<double> &B, int node, int n_procs);
template Mat<double> matrix_square<double>(const Mat<double> &A, int node, int n_procs);
template void mat2array<double>(double **A, const Mat<double> &B);
template void array2mat<double>(double **A, Mat<double> &B);
template void mat2array_square<double>(double **A, const Mat<double> &B);
template void array2mat_square<double>(double **A, Mat<double> &B);
template void get_mpi_type<double>(MPI_Datatype *mpi_type);
template void test_matrix_multiplication<double>(int rows, int node, int n_procs, int mode = 0, string mat_A = "", string mat_B = "", string mat_C = "");
template void test_matrix_square<double>(int rows, int node, int n_procs, int mode = 0, string mat_A = "", string mat_C = "");

template <class T>
void mat2array(T **A, const Mat<T> &B) {
    long n_rows = B.n_rows, n_cols = B.n_cols;
    for (long i = 0; i < n_rows; ++i)
        for (long j = 0; j < n_cols; ++j)
            (*A)[i * n_cols + j] = B(i, j);
}

template <class T>
void array2mat(T **A, Mat<T> &B) {
    long n_rows = B.n_rows, n_cols = B.n_cols;
    for (long i = 0; i < n_rows; ++i)
        for (long j = 0; j < n_cols; ++j)
            B(i, j) = (*A)[i * n_cols + j];
}

template <class T>
void mat2array_square(T **A, const Mat<T> &B) {
    long n_rows = B.n_rows;
    for (long i = 0; i < n_rows; ++i)
        for (long j = 0; j <= i; ++j)
            (*A)[i * (i + 1) / 2 + j] = B(i, j);
}

template <class T>
void array2mat_square(T **A, Mat<T> &B) {
    long n_rows = B.n_rows;
    for (long i = 0; i < n_rows; ++i) {
        for (long j = 0; j <= i; ++j) {
            B(i, j) = (*A)[i * (i + 1) / 2 + j];
            B(j, i) = (*A)[i * (i + 1) / 2 + j];
        }
    }
}

template<class T>
void get_mpi_type(MPI_Datatype *mpi_type) {
    auto name = typeid(T).name();
    if (strcmp(name, "c") == 0) {
        *mpi_type = MPI_CHAR;
    } else if (strcmp(name, "h") == 0) {
        *mpi_type = MPI_UNSIGNED_CHAR;
    } else if (strcmp(name, "s") == 0) {
        *mpi_type = MPI_SHORT;
    } else if (strcmp(name, "t") == 0) {
        *mpi_type = MPI_UNSIGNED_SHORT;
    } else if (strcmp(name, "i") == 0) {
        *mpi_type = MPI_INT;
    } else if (strcmp(name, "j") == 0) {
        *mpi_type = MPI_UNSIGNED;
    } else if (strcmp(name, "l") == 0) {
        *mpi_type = MPI_LONG;
    } else if (strcmp(name, "m") == 0) {
        *mpi_type = MPI_UNSIGNED_LONG;
    } else if (strcmp(name, "x") == 0) {
        *mpi_type = MPI_LONG_LONG;
    } else if (strcmp(name, "y") == 0) {
        *mpi_type = MPI_UNSIGNED_LONG_LONG;
    } else if (strcmp(name, "f") == 0) {
        *mpi_type = MPI_FLOAT;
    } else if (strcmp(name, "d") == 0) {
        *mpi_type = MPI_DOUBLE;
    } else if (strcmp(name, "e") == 0) {
        *mpi_type = MPI_LONG_DOUBLE;
    } else {
        *mpi_type = MPI_SHORT;
    }
}

pair<vector<pair<int, int>>, vector<pair<int, int>>> configure_cpu_triangular(int rows, int cols, int n_procs) {
    // lower trigular (including diagonal) part of m x m grid processors (total m x (m + 1) / 2 processors)
    int n_procs_row = floor(sqrt(n_procs * 2 + 1.0 / 4) - 1.0 / 2);
    n_procs_row = n_procs_row > rows ? rows : n_procs_row;
    int n_procs_col = n_procs_row;
    int n_procs_required = n_procs_row * (n_procs_row + 1) / 2;
    int rows_per_proc = floor((double)rows / n_procs_row);
    int cols_per_proc = floor((double)cols / n_procs_col);
    int t_row = rows - rows_per_proc * n_procs_row;
    int t_col = cols - cols_per_proc * n_procs_col;

    vector<pair<int, int>> row_idx, col_idx;
    for (int k = 0; k < n_procs_required; ++k) {
        int i = floor(sqrt(2 * k + 1.0 / 4) - 1.0 / 2);
        int j = k - i * (i + 1) / 2;

        if (i < t_row) {
            row_idx.push_back({i * (rows_per_proc + 1), (i + 1) * (rows_per_proc + 1) - 1});
        } else {
            row_idx.push_back({t_row * (rows_per_proc + 1) + (i - t_row) * rows_per_proc, t_row * (rows_per_proc + 1) + (i - t_row + 1) * rows_per_proc - 1});
        }

        if (j < t_col) {
            col_idx.push_back({j * (cols_per_proc + 1), (j + 1) * (cols_per_proc + 1) - 1});
        } else {
            col_idx.push_back({t_col * (cols_per_proc + 1) + (j - t_col) * cols_per_proc, t_col * (cols_per_proc + 1) + (j - t_col + 1) * cols_per_proc - 1});
        }
    }
    // printf("rows: %d, cols: %d\n", rows, cols);
    // printf("n_procs: %d, n_procs_required: %d, n_procs_row: %d, n_procs_col: %d\n", n_procs, n_procs_required, n_procs_row, n_procs_col);
    // printf("rows_per_proc: %d, cols_per_proc: %d\n", rows_per_proc, cols_per_proc);
    // printf("t_row: %d, t_col: %d\n", t_row, t_col);
    // for (int i = 0; i < (int)col_idx.size(); ++i) {
    //     printf("row_idx[%d]: %d - %d\n", i, row_idx[i].first, row_idx[i].second);
    //     printf("col_idx[%d]: %d - %d\n", i, col_idx[i].first, col_idx[i].second);
    // }
    return {row_idx, col_idx};
}

pair<vector<pair<int, int>>, vector<pair<int, int>>> configure_cpu_rectangular(int rows, int cols, int n_procs, int n_procs_row, int n_procs_col) {
    // n_procs_row x n_procs_col grid processors
    if (n_procs_row == 0 || n_procs_col == 0) {
        n_procs_row = floor(sqrt(n_procs));
        n_procs_col = floor(sqrt(n_procs));
    }
    n_procs_row = n_procs_row > rows ? rows : n_procs_row;
    n_procs_col = n_procs_col > cols ? cols : n_procs_col;
    int rows_per_proc = floor((double)rows / n_procs_row);
    int cols_per_proc = floor((double)cols / n_procs_col);
    int t_row = rows - rows_per_proc * n_procs_row;
    int t_col = cols - cols_per_proc * n_procs_col;

    vector<pair<int, int>> row_idx;
    for (int i = 0; i < n_procs_row; ++i) {
        if (i < t_row) {
            row_idx.push_back({i * (rows_per_proc + 1), (i + 1) * (rows_per_proc + 1) - 1});
        } else {
            row_idx.push_back({t_row * (rows_per_proc + 1) + (i - t_row) * rows_per_proc, t_row * (rows_per_proc + 1) + (i - t_row + 1) * rows_per_proc - 1});
        }
    }

    vector<pair<int, int>> col_idx;
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

template <class T>
Mat<T> matrix_multiplication(const Mat<T> &A, const Mat<T> &B, int node, int n_procs) {
    shino::precise_stopwatch stopwatch;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;
    MPI_Datatype mpi_type;
    get_mpi_type<T>(&mpi_type);
    long rows = A.n_rows, cols = B.n_cols;

    Mat<T> C(rows, cols, fill::zeros);

    if (n_procs == 1) {
        C = A * B;
        return C;
    }

    pair<vector<pair<int, int>>, vector<pair<int, int>>> ranges = configure_cpu_rectangular(rows, cols, n_procs);
    vector<pair<int, int>> row_idx = ranges.first;
    vector<pair<int, int>> col_idx = ranges.second;
    int n_procs_row = (int)row_idx.size();
    int n_procs_col = (int)col_idx.size();
    int n_procs_required = n_procs_row * n_procs_col;

    // int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype * newtype)
    long limit = INT_MAX;
    long mpi_count = ceil((double)rows * cols / limit);
    long mpi_unit = ceil((double)rows * cols / mpi_count);
    MPI_Datatype MPI_LARGE;
    if (rows * cols > limit) {
        MPI_Type_contiguous(mpi_unit, mpi_type, &MPI_LARGE);
        MPI_Type_commit(&MPI_LARGE);
        if (node == 0) {
            printf("CPU %d: rows: %ld, cols: %ld, limit: %ld, mpi_count: %ld, mpi_unit: %ld\n", node, rows, cols, limit, mpi_count, mpi_unit);
        }
    }

    if (node == 0) {
        long r1 = row_idx[0].first;
        long r2 = row_idx[0].second;
        long c1 = col_idx[0].first;
        long c2 = col_idx[0].second;
        long size = (r2 - r1 + 1) * (c2 - c1 + 1);
        try {
            C.submat(r1, c1, r2, c2) = A.rows(r1, r2) * B.cols(c1, c2);
        } catch (const std::exception &e) {
            cout << e.what() << endl;
            printf("CPU %d: row_idx[%d]: r1 - r2: %ld - %ld, col_idx[%d] c1 - c2: %ld - %ld, A.n_rows: %ld, B.n_cols: %ld\n", node, 0, r1, r2, 0, c1, c2, rows, cols);
        }
        auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        printf("CPU %d: CPU used: %d, Receive C(%d, %d) from CPU %d (%d, %d): %ld-by-%ld, %ld", node, n_procs_required, 0, 0, 0, 0, 0, r2 - r1 + 1, c2 - c1 + 1, size);
        cout << " (" << elapsed_time << " ms)\n";

        for (int k = 1; k < n_procs_required; ++k) {
            T *c = new T[size];
            MPI_Recv(c, size, mpi_type, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
            int tag = status.MPI_TAG;
            int i = tag / n_procs_col;
            int j = tag % n_procs_col;
            long r1 = row_idx[i].first;
            long r2 = row_idx[i].second;
            long c1 = col_idx[j].first;
            long c2 = col_idx[j].second;
            Mat<T> C_sub(r2 - r1 + 1, c2 - c1 + 1, fill::zeros);
            array2mat<T>(&c, C_sub);
            delete[] c;
            C.submat(r1, c1, r2, c2) = C_sub;
            elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
            printf("CPU %d: CPU used: %d, Receive C(%d, %d) from CPU %d (%d, %d): %ld-by-%ld, %ld", node, n_procs_required, i, j, i * n_procs_col + j, i, j, (long)C_sub.n_rows, (long)C_sub.n_cols, (long)(C_sub.n_rows * C_sub.n_cols));
            cout << " (" << elapsed_time << " ms)\n";
        }

        T *c;
        size = rows * cols;
        MPI_Barrier(comm);
        if (size > limit) {
            c = new T[mpi_count * mpi_unit];
            mat2array(&c, C);
            MPI_Bcast(c, mpi_count, MPI_LARGE, 0, comm);
        } else {
            c = new T[size];
            mat2array(&c, C);
            MPI_Bcast(c, size, mpi_type, 0, comm);
        }
        delete[] c;
    }

    if (node > 0 && node < n_procs_required) {
        int i = node / n_procs_col;
        int j = node % n_procs_col;
        long r1 = row_idx[i].first;
        long r2 = row_idx[i].second;
        long c1 = col_idx[j].first;
        long c2 = col_idx[j].second;
        long size = (r2 - r1 + 1) * (c2 - c1 + 1);
        Mat<T> TC;
        try {
            TC = A.rows(r1, r2) * B.cols(c1, c2);
        } catch (const std::exception &e) {
            cout << e.what() << endl;
            printf("CPU %d: row_idx[%d]: r1 - r2: %ld - %ld, col_idx[%d] c1 - c2: %ld - %ld, A.n_rows: %ld, B.n_cols: %ld\n", node, i, r1, r2, j, c1, c2, rows, cols);
        }
        T *c = new T[size];
        mat2array(&c, TC);
        // printf("CPU %d: Send C(%d, %d) to CPU 0 (0, 0)\n", node, i, j);
        MPI_Send(c, size, mpi_type, 0, node, comm);
        delete[] c;
    }

    if (node > 0) {
        T *c;
        long size = rows * cols;
        MPI_Barrier(comm);
        if (size > limit) {
            c = new T[mpi_count * mpi_unit];
            MPI_Bcast(c, mpi_count, MPI_LARGE, 0, comm);
        } else {
            c = new T[size];
            MPI_Bcast(c, size, mpi_type, 0, comm);
        }
        array2mat<T>(&c, C);
        delete[] c;
    }
    return C;
}

template <class T>
Mat<T> matrix_square(const Mat<T> &A, int node, int n_procs) {
    shino::precise_stopwatch stopwatch;

    MPI_Datatype mpi_type;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Status status;
    get_mpi_type<T>(&mpi_type);
    long rows = A.n_rows;
    long cols = A.n_cols;
    Mat<T> C(rows, cols, fill::zeros);

    if (n_procs == 1) {
        C = A * A;
        return C;
    }

    pair<vector<pair<int, int>>, vector<pair<int, int>>> ranges = configure_cpu_triangular(rows, cols, n_procs);
    vector<pair<int, int>> row_idx = ranges.first;
    vector<pair<int, int>> col_idx = ranges.second;
    int n_procs_required = (int)col_idx.size();

    if (node == 0) {
        long r1 = row_idx[node].first;
        long r2 = row_idx[node].second;
        long c1 = col_idx[node].first;
        long c2 = col_idx[node].second;
        long size = (r2 - r1 + 1) * (c2 - c1 + 1);
        C.submat(r1, c1, r2, c2) = A.rows(r1, r2) * A.cols(c1, c2);
        auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        printf("CPU_required: %d: CPU %d: Receive C(%d, %d) from CPU %d (%d, %d): %ld-by-%ld, %ld", n_procs_required, node, 0, 0, 0, 0, 0, r2 - r1 + 1, c2 - c1 + 1, size);
        cout << " (" << elapsed_time << " ms)\n";

        for (int k = 1; k < n_procs_required; ++k) {
            T *c = new T[size];
            MPI_Recv(c, size, mpi_type, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
            int tag = status.MPI_TAG;
            long r1 = row_idx[tag].first;
            long r2 = row_idx[tag].second;
            long c1 = col_idx[tag].first;
            long c2 = col_idx[tag].second;
            Mat<T> C_sub(r2 - r1 + 1, c2 - c1 + 1, fill::zeros);
            array2mat<T>(&c, C_sub);
            C.submat(r1, c1, r2, c2) = C_sub;
            int i = floor(sqrt(2 * tag + 1.0 / 4) - 1.0 / 2);
            int j = tag - i * (i + 1) / 2;
            if (i != j)
                C.submat(c1, r1, c2, r2) = C_sub.t();
            elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
            printf("CPU_required: %d: CPU %d: Receive C(%d, %d) from CPU %d (%d, %d): %ld-by-%ld, %ld", n_procs_required, node, i, j, tag, i, j, (long)C_sub.n_rows, (long)C_sub.n_cols, (long)(C_sub.n_rows * C_sub.n_cols));
            cout << " (" << elapsed_time << " ms)\n";
            delete[] c;
        }

        size = rows * (cols + 1) / 2;
        T *c = new T[size];
        mat2array_square<T>(&c, C);
        MPI_Barrier(comm);
        MPI_Bcast(c, size, mpi_type, 0, comm);
        delete[] c;

        // cout << "Number of CPUs: " << n_procs << ", rows: " << rows;
        // cout << ", wall clock time elapsed to perform matrix_multiplication: " << elapsed_time << " ms" << endl;
    }

    if (node > 0 && node < n_procs_required) {
        long r1 = row_idx[node].first;
        long r2 = row_idx[node].second;
        long c1 = col_idx[node].first;
        long c2 = col_idx[node].second;
        // printf("%02d: r1: %ld, r2: %ld, c1: %ld, c2: %ld\n", node, r1, r2, c1, c2);
        Mat<T> TC;
        TC = A.rows(r1, r2) * A.cols(c1, c2);
        long size = (r2 - r1 + 1) * (c2 - c1 + 1);
        T *c = new T[size];
        mat2array(&c, TC);
        // printf("CPU %d: Send C(%d, %d) to CPU 0 (0, 0)\n", node, i, j);
        MPI_Send(c, size, mpi_type, 0, node, comm);
        delete[] c;
    }

    if (node > 0) {
        long size = rows * (cols + 1) / 2;
        T *c = new T[size];
        MPI_Barrier(comm);
        MPI_Bcast(c, size, mpi_type, 0, comm);
        array2mat_square<T>(&c, C);
        delete[] c;
    }
    return C;
}

template <class T>
void test_matrix_multiplication(int rows, int node, int n_procs, int mode, string mat_A, string mat_B, string mat_C) {
    shino::precise_stopwatch stopwatch;

    if (node == 0) {
        cout << "-----------------------------------" << endl;
        cout << "Number of CPUs: " << n_procs << endl;
        cout << "Loading/generating A and B..." << endl;
        if (mat_A.size() > 0) {
            cout << "Reading from " << mat_A << endl;
            cout << "Reading from " << mat_B << endl;
        }
    }

    Mat<T> A;
    Mat<T> B;
    if (mode == 1) {
        A.load(mat_A, csv_ascii);
        B.load(mat_B, csv_ascii);
    } else {
        A.randn(rows, rows);
        B.randn(rows, rows);
    }
    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();

    if (node == 0) {
        cout << "Performing matrix multiplication A * B..."
             << " (" << elapsed_time << " ms)" << endl;
    }
    Mat<T> C = matrix_multiplication<T>(A, B, node, n_procs);
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    if (mode == 2) {
        mat_C.replace(mat_C.end() - 4, mat_C.end(), "_matrix_multiplication_" + to_string(node) + ".csv");
        C.save(mat_C, csv_ascii);
    } else if (mode == 1) {
        if (node == 0) {
            mat_C.replace(mat_C.end() - 4, mat_C.end(), "_matrix_multiplication.csv");
            cout << "Writing to " << mat_C << " (" << elapsed_time << " ms)" << endl;
            C.save(mat_C, csv_ascii);
        }
    }

    if (node == 0) {
        elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        cout << "Number of CPUs: " << n_procs << ", rows: " << C.n_rows;
        cout << ", total wall clock time elapsed: " << elapsed_time << " ms" << endl;
        cout << "-----------------------------------" << endl;
    }
}

template <class T>
void test_matrix_square(int rows, int node, int n_procs, int mode, string mat_A, string mat_C) {
    shino::precise_stopwatch stopwatch;

    if (node == 0) {
        cout << "-----------------------------------" << endl;
        cout << "Number of CPUs: " << n_procs << endl;
        cout << "Loading/generating A..." << endl;
        if (mat_A.size() > 0) {
            cout << "Reading from " << mat_A << endl;
        }
    }

    Mat<T> A;
    if (mode == 1) {
        A.load(mat_A, csv_ascii);
    } else {
        A.randn(rows, rows);
    }
    auto elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();

    if (node == 0) {
        cout << "Performing matrix squaring A * A..."
             << " (" << elapsed_time << " ms)" << endl;
    }
    Mat<T> C = matrix_square<T>(A, node, n_procs);
    elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
    if (node == 0) {
        if (mode == 1) {
            mat_C.replace(mat_C.end() - 4, mat_C.end(), "_matrix_square.csv");
            cout << "Writing to " << mat_C << " (" << elapsed_time << " ms)" << endl;
            C.save(mat_C, csv_ascii);
        }
        elapsed_time = stopwatch.elapsed_time<unsigned int, std::chrono::milliseconds>();
        cout << "Number of CPUs: " << n_procs << ", rows: " << C.n_rows;
        cout << ", total wall clock time elapsed: " << elapsed_time << " ms" << endl;
        cout << "-----------------------------------" << endl;
    }
}

