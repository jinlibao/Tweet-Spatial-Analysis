#ifndef TWEETS_SPATIAL_ANALYSIS_H
#define TWEETS_SPATIAL_ANALYSIS_H

#include "ellipse.h"
#include <armadillo>
#include <vector>

using namespace arma;
using namespace std;

imat APD(imat A);
void build_overlap_matrix(string input_file, string output_file, int rows = 0);

#endif
