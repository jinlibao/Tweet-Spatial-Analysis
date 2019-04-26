#ifndef TWEETS_SPATIAL_ANALYSIS_H
#define TWEETS_SPATIAL_ANALYSIS_H

#include "ellipse.h"
#include "matmul.h"
#include <armadillo>
#include <vector>
#include <queue>

using namespace arma;
using namespace std;

Mat<short> APD(Mat<short> A);
void build_overlap_matrix(string input_file, string output_file, long rows = 0);
void find_components(string adj_file);
vector<pair<int, vector<int>>> bfs(Mat<short>& A);
void build_distance_matrix(string adj_file, string dis_file);
void test_matmul(int rows);
void test_parallel_matmul(int rows, int *argc, char*** argv);

#endif
