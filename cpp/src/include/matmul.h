#ifndef MATMUL_H
#define MATMUL_H

#include <armadillo>
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <utility>
#include <vector>

using namespace arma;
using namespace std;

void mat2array(short **A, const Mat<short> &B);
void array2mat(short **A, Mat<short>& B);
Mat<short> parallel_matmul(const Mat<short> &A, const Mat<short> &B, int rank, int n_procs);
pair<vector<pair<int, int>>, vector<pair<int, int>>> configure_cpu(int rows, int cols, int n_procs, int n_procs_row = 0, int n_procs_col = 0);
void test_matmul(int rows);
void test_parallel_matmul(int rows, int *argc, char*** argv);

#endif
