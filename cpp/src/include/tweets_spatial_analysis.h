#ifndef TWEETS_SPATIAL_ANALYSIS_H
#define TWEETS_SPATIAL_ANALYSIS_H

#include <armadillo>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <vector>
#include "ellipse.h"
#include "matmul.h"

using namespace arma;
using namespace std;

template <class T>
void build_adjacency_matrix(string input_file, string output_file, long rows = 0);
template <class T>
vector<pair<int, vector<int>>> bfs(Mat<T> &A);
template <class T>
void find_components(string adj_file, string outlier_file = "");
template <class T>
Mat<T> APD_recursive(const Mat<T> &A, int rank, int n_procs);
template <class T>
Mat<T> APD(const Mat<T> &A, int rank, int n_procs);
template <class T>
void build_distance_matrix(string adj_file, string dis_file, int node, int n_procs, bool non_recursive = true);
template <class T>
Mat<T> BPWM(const Mat<T> &A, const Mat<T> &B, int node, int n_procs);
template <class T>
Mat<T> get_witness_matrix(const Mat<T> &A, const Mat<T> &D, int node, int n_procs, int i, int r, string dis_file);
template <class T>
void build_witness_matrix(string adj_file, string dis_file, int r, int node, int n_procs);
template <class T>
Mat<T> compute_successor_matrix(const Mat<T> &A, const Mat<T> &D, int node, int n_procs, int i, string name);
template <class T>
void build_successor_matrix(string adj_file, string dis_file, int node, int n_procs);
template <class T>
void APSP(long rows, string ellipse_file, string adj_file, string adj_ordered_file, string dis_file, string outlier_file, int node, int n_procs);
template <class T>
void test_APD_recursive(string mat_file, int mode, int node, int n_procs);
template <class T>
void test_APD(string mat_file, int mode, int node, int n_procs);
vector<long unsigned> convert_index_path_to_id_path(Mat<long unsigned> &id_mat, vector<int> index_path);
void convert_paths_length_to_csv(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths);
void convert_paths_length_to_csv(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths, string filename);
void convert_paths_to_csv(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths);
void convert_paths_to_csv(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths, string filename);
void convert_paths_to_gml(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths, string filename);
void convert_paths_to_gml(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths);
void convert_paths_to_json(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths);
void convert_paths_to_json(Mat<long unsigned> &id_mat, vector<vector<int>> &index_paths, string filename);
void convert_paths_to_json(vector<vector<int>> &index_paths);
void convert_paths_to_json(vector<vector<int>> &index_paths, string filename);
vector<vector<int>> find_all_shortest_index_paths(Mat<int> &S, Mat<long unsigned> &id_mat);
pair<vector<int>, vector<long unsigned>> get_shortest_path_by_id(Mat<int> &S, Mat<long unsigned> &id_mat, long unsigned id_from, long unsigned id_to);
vector<int> get_shortest_path_by_index(Mat<int> &S, Mat<long unsigned> &id_mat, int from, int to);
void print_shortest_path(vector<int> &index_path, Mat<long unsigned> &id_mat);
void print_shortest_path(vector<int> &index_path, vector<long unsigned> &id_path);
void test_find_all_shortest_index_paths(string suc_file, string id_file);
void test_get_shortest_path_by_id(string suc_file, string id_file, long unsigned id_from, long unsigned id_to);
void test_get_shortest_path_by_index(string suc_file, string id_file, int from, int to);
void find_node_degree(string adj_file);

#endif
