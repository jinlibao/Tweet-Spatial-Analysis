#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <cmath>
#include <utility>
#include <algorithm>

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
    pair<double, double> newton(double c1, double c2, double c3, double d1, double d2, double d3, double t, double tol, int maxIter);

private:
    double x, y, a, b, p;
};

Ellipse::Ellipse(double x0, double y0, double a0, double b0, double p0)
{
    x = x0;
    y = y0;
    a = a0;
    b = b0;
    p = p0;
}

Ellipse::Ellipse(const Ellipse& e)
{
    x = e.x;
    y = e.y;
    a = e.a;
    b = e.b;
    p = e.p;
}

bool Ellipse::overlap(const Ellipse& e)
{
    double x1, y1, a1, b1, p1, x2, y2, a2, b2, p2;
    if (max(a, b) > max(e.a, e.b)) {
        x1 = e.x;
        y1 = e.y;
        a1 = e.a;
        b1 = e.b;
        p1 = e.p;
        x2 = x;
        y2 = y;
        a2 = a;
        b2 = b;
        p2 = p;
    } else {
        x1 = x;
        y1 = y;
        a1 = a;
        b1 = b;
        p1 = p;
        x2 = e.x;
        y2 = e.y;
        a2 = e.a;
        b2 = e.b;
        p2 = e.p;
    }

    double c1 = (a1 / a2) * cos(p1 - p2);
    double c2 =-(b1 / a2) * sin(p1 - p2);
    double c3 = ((x1 - x2) * cos(p2) + (y1 - y2) * sin(p2)) / a2;
    double d1 = (a1 / b2) * sin(p1 - p2);
    double d2 = (b1 / b2) * cos(p1 - p2);
    double d3 = (-(x1 - x2) * sin(p2) + (y1 - y2) * cos(p2)) / b2;

    pair<double, double> t = Ellipse::newton(c1, c2, c3, d1, d2, d3, 0.0, 1e-5, 100000);
    double dist1 = Ellipse::dist(c1, c2, c3, d1, d2, d3, t.first);
    double dist2 = Ellipse::dist(c1, c2, c3, d1, d2, d3, t.second);
    if (min(dist1, dist2) <= 1.0) {
        return true;
    } else {
        return false;
    }
}

double Ellipse::f(double c1, double c2, double c3, double d1, double d2, double d3, double t)
{
    return 2 * ((c1 * cos(t) + c2 * sin(t) + c3) * (-c1 * sin(t) + c2 * cos(t)) + 
                (d1 * cos(t) + d2 * sin(t) + d3) * (-d1 * sin(t) + d2 * cos(t)));
}

double Ellipse::g(double c1, double c2, double c3, double d1, double d2, double d3, double t)
{
    return 2 * ((c1 * cos(t) + c2 * sin(t) + c3) * (-c1 * cos(t) - c2 * sin(t)) + pow(-c1 * sin(t) + c2 * cos(t), 2) + 
                (d1 * cos(t) + d2 * sin(t) + d3) * (-d1 * cos(t) - d2 * sin(t)) + pow(-d1 * sin(t) + d2 * cos(t), 2));
}

pair<double, double> Ellipse::newton(double c1, double c2, double c3, double d1, double d2, double d3, double t, double tol, int maxIter)
{
    double t1, t2;
    int iter = 0;
    while (fabs(Ellipse::f(c1, c2, c3, d1, d2, d3, t)) >= tol && iter < maxIter) {
        if (Ellipse::g(c1, c2, c3, d1, d2, d3, t) == 0) {
            t = t - Ellipse::f(c1, c2, c3, d1, d2, d3, t) * 0.01;
        } else {
            t = t - Ellipse::f(c1, c2, c3, d1, d2, d3, t) / Ellipse::g(c1, c2, c3, d1, d2, d3, t);
        }
        ++iter;
    }
    t1 = (t / (2 * PI) - floor(t / (2 * PI))) * 2 * PI;
    if (t1 > PI) {
        t2 = t1 - PI;
    } else {
        t2 = t1 + PI;
    }
    return {t1, t2};
}

double Ellipse::dist(double c1, double c2, double c3, double d1, double d2, double d3, double t)
{
    return pow(c1 * cos(t) + c2 * sin(t) + c3, 2) +  pow(d1 * cos(t) + d2 * sin(t) + d3, 2);
}

#endif
