from pyproj import Proj, transform, Geod

from math import *

def lon_lat_to_east_north(lon, lat):
    try:
        east, north = transform(Proj(init='epsg:4326'), Proj(init='epsg:3857'), lon, lat)
        return east, north
    except:
        return None, None

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
    return a, b

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


