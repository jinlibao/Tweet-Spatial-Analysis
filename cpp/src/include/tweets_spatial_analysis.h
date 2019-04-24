#ifndef TWEETS_SPATIAL_ANALYSIS_H
#define TWEETS_SPATIAL_ANALYSIS_H

#include "ellipse.h"
#include <armadillo>
#include <vector>
#include <queue>

using namespace arma;
using namespace std;

imat APD(imat A);
void build_overlap_matrix(string input_file, string output_file, long rows = 0);
void find_components(string adj_file);
vector<pair<int, vector<int>>> bfs(Mat<short>& A);

#endif
