from pyproj import Proj, transform, Geod
from math import *
import matplotlib.pyplot as plt
import numpy as np

def lon_lat_to_east_north(lon, lat):
    try:
        east, north = transform(Proj(init='epsg:4326'), Proj(init='epsg:3857'), lon, lat)
        return east, north
    except:
        return None, None

def are_two_ellipses_overlapping(x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2):
    c1 =  (a1 / a2) * cos(phi1 - phi2)
    c2 = -(b1 / a2) * sin(phi1 - phi2)
    c3 =  ((x1 - x2) * cos(phi2) + (y1 - y2) * sin(phi2)) / a2
    d1 =  (a1 / b2) * sin(phi1 - phi2)
    d2 = -(b1 / b2) * cos(phi1 - phi2)
    d3 = (-(x1 - x2) * sin(phi2) + (y1 - y2) * cos(phi2)) / b2
    t1, t2 = newton(f, g, c1, c2, c3, d1, d2, d3, 0)
    dist1 = f(c1, c2, c3, d1, d2, d3, t1)
    dist2 = f(c1, c2, c3, d1, d2, d3, t2)
    if min(dist1, dist2) <= 1:
        return True
    else:
        return False

def plot_ellipses(x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2):
    c1 =  (a1 / a2) * cos(phi1 - phi2)
    c2 = -(b1 / a2) * sin(phi1 - phi2)
    c3 =  ((x1 - x2) * cos(phi2) + (y1 - y2) * sin(phi2)) / a2
    d1 =  (a1 / b2) * sin(phi1 - phi2)
    d2 = -(b1 / b2) * cos(phi1 - phi2)
    d3 = (-(x1 - x2) * sin(phi2) + (y1 - y2) * cos(phi2)) / b2

    pi = atan(1) * 4
    t = np.linspace(0, 2 * pi, 100)
    ex1 = a1 * np.cos(t) * cos(phi1) - b1 * np.sin(t) * sin(phi1) + x1
    ey1 = a1 * np.cos(t) * sin(phi1) + b1 * np.sin(t) * cos(phi1) + y1
    ex2 = a2 * np.cos(t) * cos(phi2) - b2 * np.sin(t) * sin(phi2) + x2
    ey2 = a2 * np.cos(t) * sin(phi2) + b2 * np.sin(t) * cos(phi2) + y2
    ex3 = c1 * np.cos(t) + c2 * np.sin(t) + c3
    ey3 = d1 * np.cos(t) + d2 * np.sin(t) + d3
    ex4 = np.cos(t)
    ey4 = np.sin(t)
    plt.figure()
    plt.plot(ex1, ey1, ex2, ey2, ex3, ey3, ex4, ey4)
    plt.axis('equal')
    plt.show()

def f(c1, c2, c3, d1, d2, d3, t):
    return (c1 * cos(t) + c2 * sin(t) + c3) ** 2 + \
           (d1 * cos(t) + d2 * sin(t) + d3) ** 2

def g(c1, c2, c3, d1, d2, d3, t):
    return 2 * ((c1 * cos(t) + c2 * sin(t) + c3) * (-c1 * sin(t) + c2 * cos(t)) + \
                (d1 * cos(t) + d2 * sin(t) + d3) * (-d1 * sin(t) + d2 * cos(t)))

def h(c1, c2, c3, d1, d2, d3, t):
    return 2 * ((c1 * cos(t) + c2 * sin(t) + c3) * (-c1 * cos(t) - c2 * sin(t)) + \
                (-c1 * sin(t) + c2 * cos(t)) ** 2 + \
                (d1 * cos(t) + d2 * sin(t) + d3) * (-d1 * cos(t) - d2 * sin(t)) + \
                (-d1 * sin(t) + d2 * cos(t)) ** 2)

def newton(f, g, c1, c2, c3, d1, d2, d3, t0, tol=1e-6, maxIter=100):
    pi = atan(1) * 4
    iter = 0
    t = t0
    while iter < maxIter and fabs(f(c1, c2, c3, d1, d2, d3, t)) >= tol:
        if g(c1, c2, c3, d1, d2, d3, t) == 0:
            t = t - f(c1, c2, c3, d1, d2, d3, t) / 1
        else:
            t = t - f(c1, c2, c3, d1, d2, d3, t) / g(c1, c2, c3, d1, d2, d3, t)
        iter += 1
    t1 = (t / (2 * pi) - floor(t / (2 * pi))) * 2 * pi
    if t1 > pi:
        t2 = t1 - pi
    else:
        t2 = t1 + pi
    return (t1, t2)

# Point and ellipse (rotated) position test: algorithm
# https://stackoverflow.com/questions/7946187/point-and-ellipse-rotated-position-test-algorithm
def is_point_inside_ellipse(xo, yo, a, b, theta, xp, yp):
    #theta = pi / 2.0 - theta
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    xp_minus_xo = xp - xo
    yp_minus_yo = yp - yo
    part1 = ((cos_theta * xp_minus_xo + sin_theta * yp_minus_yo)**2) / a**2
    part2 = ((sin_theta * xp_minus_xo - cos_theta * yp_minus_yo)**2) / b**2

    if part1 + part2 <= 1.0:
        return True
    else:
        return False

def calculate_ellipse_a_and_b_from_area_and_ratio(area, xy_ratio):
    # Code to calculate a and b ellispe parameters from the area and the x/y ratio:
    # The original area data is defined in km squared
    # Convert a and b into meters.
    a = sqrt((area * xy_ratio) / pi) * 1000.0
    b = sqrt(area / (xy_ratio * pi)) * 1000.0
    return a, b, xy_ratio

def calculate_distance(x0, y0, x1, y1):
    return sqrt((x0 - x1)**2 + (y0 - y1)**2)

# Haversine Formula in Python (Bearing and Distance between two GPS points)
# https://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
def calculate_distance_non_projected_havershine(lat1, lon1, lat2, lon2):
    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6372.8 # 6371
    return c * r


