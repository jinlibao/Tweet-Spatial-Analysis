#ifndef MATMUL_H
#define MATMUL_H

#include <armadillo>
#include <climits>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <mpi.h>
#include <typeinfo>
#include <utility>
#include <vector>

using namespace arma;
using namespace std;

// template function declaration
template <class T> Mat<T> matrix_multiplication(const Mat<T> &A, const Mat<T> &B, int node, int n_procs);
template <class T> Mat<T> matrix_square(const Mat<T> &A, int node, int n_procs);
template <class T> void mat2array(T **A, const Mat<T> &B);
template <class T> void array2mat(T **A, Mat<T> &B);
template <class T> void mat2array_square(T **A, const Mat<T> &B);
template <class T> void array2mat_square(T **A, Mat<T> &B);
template <class T> void get_mpi_type(MPI_Datatype *mpi_type);
template <class T> void test_matrix_multiplication(int rows, int *argc, char ***argv, int mode = 0, string mat_A = "", string mat_B = "", string mat_C = "");
template <class T> void test_matrix_square(int rows, int *argc, char ***argv, int mode = 0, string mat_A = "", string mat_C = "");
pair<vector<pair<int, int>>, vector<pair<int, int>>> configure_cpu_triangular(int rows, int cols, int n_procs);
pair<vector<pair<int, int>>, vector<pair<int, int>>> configure_cpu_rectangular(int rows, int cols, int n_procs, int n_procs_row = 0, int n_procs_col = 0);

#endif
