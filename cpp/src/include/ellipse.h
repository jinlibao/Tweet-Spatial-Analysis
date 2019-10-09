#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <utility>

using namespace std;

const double PI = atan(1) * 4;

class Ellipse {
 public:
  Ellipse(double x, double y, double a, double b, double p);
  Ellipse(const Ellipse& e);
  bool overlap(const Ellipse& e);
  double f(double c1, double c2, double c3, double d1, double d2, double d3, double t);
  double g(double c1, double c2, double c3, double d1, double d2, double d3, double t);
  double dist(double c1, double c2, double c3, double d1, double d2, double d3, double t);
  double newton(double c1, double c2, double c3, double d1, double d2, double d3, double t, double tol, int maxIter);

 private:
  double x, y, a, b, p;
};

#endif
