#ifndef MATMUL_H
#define MATMUL_H

#include <armadillo>
#include <climits>
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <utility>
#include <vector>

using namespace arma;
using namespace std;

Mat<short> parallel_matmul(const Mat<short> &A, const Mat<short> &B, int node, int n_procs);
Mat<short> parallel_matsq(const Mat<short> &A, const Mat<short> &B, int node, int n_procs);
void mat2array(short **A, const Mat<short> &B);
void array2mat(short **A, Mat<short> &B);
void mat2array_square(short **A, const Mat<short> &B);
void array2mat_square(short **A, Mat<short> &B);
void test_matmul(int rows, int mode = 0, string mat_A = "", string mat_B = "", string mat_C = "");
void test_matsq(int rows, int mode = 0, string mat_A = "", string mat_C = "");
void test_parallel_matmul(int rows, int *argc, char ***argv, int mode = 0, string mat_A = "", string mat_B = "", string mat_C = "");
void test_parallel_matsq(int rows, int *argc, char ***argv, int mode = 0, string mat_A = "", string mat_C = "");
pair<vector<pair<int, int>>, vector<pair<int, int>>> configure_cpu_triangular(int rows, int cols, int n_procs);
pair<vector<pair<int, int>>, vector<pair<int, int>>> configure_cpu_rectangular(int rows, int cols, int n_procs, int n_procs_row = 0,
                                                                               int n_procs_col = 0);

#endif
