#ifndef TWEETS_SPATIAL_ANALYSIS_H
#define TWEETS_SPATIAL_ANALYSIS_H

#include "ellipse.h"
#include "matmul.h"
#include <armadillo>
#include <queue>
#include <vector>

using namespace arma;
using namespace std;

template <class T> void build_adjacency_matrix(string input_file, string output_file, long rows = 0);
template <class T> vector<pair<int, vector<int>>> bfs(Mat<T> &A);
template <class T> void find_components(string adj_file, string outlier_file = "");
template <class T> Mat<T> APD_recursive(const Mat<T> &A, int rank, int n_procs);
template <class T> Mat<T> APD(const Mat<T> &A, int rank, int n_procs);
template <class T> void build_distance_matrix(string adj_file, string dis_file, int node, int n_procs, bool non_recursive = true);
template <class T> Mat<T> BPWM(const Mat<T> &A, const Mat<T> &B, int node, int n_procs);
template <class T> Mat<T> compute_successor_matrix(const Mat<T> &A, const Mat<T> &D, int node, int n_procs);
template <class T> void build_successor_matrix(string adj_file, string dis_file, int node, int n_procs);
template <class T> void APSP(long rows, string ellipse_file, string adj_file, string adj_ordered_file, string dis_file, string outlier_file, int node, int n_procs);
template <class T> void test_APD_recursive(string mat_file, int mode, int node, int n_procs);
template <class T> void test_APD(string mat_file, int mode, int node, int n_procs);

#endif
