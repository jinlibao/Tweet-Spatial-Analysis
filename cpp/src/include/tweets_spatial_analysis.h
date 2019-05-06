#ifndef TWEETS_SPATIAL_ANALYSIS_H
#define TWEETS_SPATIAL_ANALYSIS_H

#include "ellipse.h"
#include <armadillo>
#include <vector>
#include <queue>

using namespace arma;
using namespace std;

void build_overlap_matrix(string input_file, string output_file, long rows = 0);
void find_components(string adj_file);
vector<pair<int, vector<int>>> bfs(Mat<float>& A);
void build_distance_matrix(string adj_file, string dis_file);
void build_distance_matrix_parallel(string adj_file, string dis_file, int *argc, char ***argv);
Mat<float> APD(const Mat<float>& A);
Mat<float> APD_parallel(const Mat<float> &A, int rank, int n_procs);
void test_APD(string mat_file);
void test_APD_parallel(string mat_file, int mode, int *argc, char ***argv);

#endif
