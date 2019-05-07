#ifndef TWEETS_SPATIAL_ANALYSIS_H
#define TWEETS_SPATIAL_ANALYSIS_H

#include "ellipse.h"
#include "matmul.h"
#include <armadillo>
#include <queue>
#include <vector>

using namespace arma;
using namespace std;

template <class T> void build_overlap_matrix(string input_file, string output_file, long rows = 0);
template <class T> void find_components(string adj_file);
template <class T> vector<pair<int, vector<int>>> bfs(Mat<T> &A);
template <class T> void build_distance_matrix(string adj_file, string dis_file);
template <class T> void build_distance_matrix_parallel(string adj_file, string dis_file, int *argc, char ***argv);
template <class T> Mat<T> APD(const Mat<T> &A);
template <class T> Mat<T> APD_parallel(const Mat<T> &A, int rank, int n_procs);
template <class T> Mat<T> APD_parallel_non_recursive(const Mat<T> &A, int rank, int n_procs);
template <class T> void test_APD(string mat_file);
template <class T> void test_APD_parallel(string mat_file, int mode, int *argc, char ***argv);
template <class T> void test_APD_parallel_non_recursive(string mat_file, int mode, int *argc, char ***argv);

#endif
