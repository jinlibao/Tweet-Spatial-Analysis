from pyproj import Proj, transform, Geod
from math import *
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def lon_lat_to_east_north(lon, lat):
    try:
        east, north = transform(Proj(init='epsg:4326'), Proj(init='epsg:3857'), lon, lat)
        return east, north
    except:
        return None, None

def are_two_ellipses_overlapping(x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2):
    # Make sure the second ellipse has large scaling constant
    if max(a1, b1) > max(a2, b2):
        x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2 = x2, y2, a2, b2, phi2, x1, y1, a1, b1, phi1
    c1 =  (a1 / a2) * cos(phi1 - phi2)
    c2 = -(b1 / a2) * sin(phi1 - phi2)
    c3 =  ((x1 - x2) * cos(phi2) + (y1 - y2) * sin(phi2)) / a2
    d1 =  (a1 / b2) * sin(phi1 - phi2)
    d2 =  (b1 / b2) * cos(phi1 - phi2)
    d3 = (-(x1 - x2) * sin(phi2) + (y1 - y2) * cos(phi2)) / b2

    iter_n = 1
    pi = atan(1) * 4
    t = newton(f, g, c1, c2, c3, d1, d2, d3, 0)
    t_old = t
    while dist(c1, c2, c3, d1, d2, d3, t) > 1 and iter_n < 4:
        t = t_old - pi / 2 * iter_n
        t = newton(f, g, c1, c2, c3, d1, d2, d3, t)
        iter_n += 1

    if dist(c1, c2, c3, d1, d2, d3, t) <= 1:
        return True
    else:
        return False

def plot_shortest_path(data, filename='shortest_path.pdf'):
    plt.style.use('ggplot')
    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)

    color = [
        (209, 85, 62), (77, 136, 185), (150, 142, 208), (119, 119, 119),
        (242, 195, 110), (151, 185, 85), (244, 184, 185), (157, 209, 199),
        (255, 255, 188), (189, 186, 215), (235, 135, 119), (138, 176, 208),
        (242, 183, 112), (188, 221, 120), (245, 207, 228),
        (31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40),
        (148, 103, 189), (140, 86, 75), (227, 119, 194), (127, 127, 127),
        (188, 189, 34), (23, 190, 207), (174, 199, 232), (255, 187, 120),
        (152, 223, 138), (255, 152, 150), (197, 176, 213),
        (196, 156, 148)
    ]
    # Scale the RGB values to the [0, 1] range
    for i in range(len(color)):
        r, g, b = color[i]
        color[i] = (r / 255., g / 255., b / 255.)

    pi = atan(1) * 4
    t = np.linspace(0, 2 * pi, 100)
    with PdfPages(filename) as pdf:
        fig = plt.figure()
        handles = []
        for i, datum in enumerate(data):
            x, y, a, b, phi, id_name = datum
            ex = a * np.cos(t) * np.cos(phi) - b * np.sin(t) * np.sin(phi) + x
            ey = a * np.cos(t) * np.sin(phi) + b * np.sin(t) * np.cos(phi) + y
            # plt.plot(ex, ey)
            # plt.annotate(str(int(id_name)), (x, y))
            handle, = plt.plot(ex, ey, linewidth=0.5, label=str(int(id_name)), color=color[i])
            handles.append(handle)
        plt.axis('equal')
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.legend(handles=handles, loc='best')
        plt.title('Shortest path from {} to {}'.format(str(int(data[0][5])), str(int(data[-1][5]))))
        plt.show(block=False)
        pdf.savefig(fig)

def plot_ellipses(x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2, filename="plot.pdf"):

    if max(a1, b1) > max(a2, b2):
        x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2 = x2, y2, a2, b2, phi2, x1, y1, a1, b1, phi1

    c1 =  (a1 / a2) * cos(phi1 - phi2)
    c2 = -(b1 / a2) * sin(phi1 - phi2)
    c3 =  ((x1 - x2) * cos(phi2) + (y1 - y2) * sin(phi2)) / a2
    d1 =  (a1 / b2) * sin(phi1 - phi2)
    d2 =  (b1 / b2) * cos(phi1 - phi2)
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
    ex = [[ex1, ex2], [ex3, ex4]]
    ey = [[ey1, ey2], [ey3, ey4]]
    label = ['Before Inverse Transformation', 'After Inverse Transformation']
    style = [['r:', 'b:'], ['r-', 'b-']]

    dv = dist(c1, c2, c3, d1, d2, d3, t)
    fv = f(c1, c2, c3, d1, d2, d3, t)
    gv = g(c1, c2, c3, d1, d2, d3, t)
    d_list = [dv, fv, gv]
    func_name = ['$\mathrm{dist}(t)$', '$\mathrm{dist}\'(t)$', '$\mathrm{dist}\'\'(t)$']

    plt.style.use('ggplot')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    with PdfPages(filename) as pdf:
        for i in range(len(ex)):
            fig = plt.figure()
            plt.plot(ex[i][0], ey[i][0], 'r', ex[i][1], ey[i][1], 'b')
            plt.axis('equal')
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.title(label[i])
            plt.show(block=False)
            pdf.savefig(fig)
            plt.close()

        fig = plt.figure()
        handles = []
        for i in range(len(ex)):
            # fig = plt.figure()
            cur_h_0, = plt.plot(ex[i][0], ey[i][0], style[i][0], label=label[i])
            handles.append(cur_h_0)
            cur_h_1, = plt.plot(ex[i][1], ey[i][1], style[i][1], label=label[i])
            handles.append(cur_h_1)
            plt.axis('equal')
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.title('Before vs. After Inverse Transformation')
        plt.legend(handles=handles, loc='best')
        plt.show(block=False)
        pdf.savefig(fig)
        plt.close()

        fig = plt.figure()
        handles = []
        for i in range(len(d_list)):
            handle, = plt.plot(t, d_list[i], label=func_name[i])
            handles.append(handle)
        plt.xlabel('$t$')
        plt.ylabel('$y$')
        plt.legend(handles=handles, loc='best')
        plt.title('{}, {}, {} vs. $t$'.format(func_name[0], func_name[1], func_name[2]))
        plt.show(block=False)
        pdf.savefig(fig)
        plt.close()

def dist(c1, c2, c3, d1, d2, d3, t):
    return (c1 * np.cos(t) + c2 * np.sin(t) + c3) ** 2 + \
           (d1 * np.cos(t) + d2 * np.sin(t) + d3) ** 2

def f(c1, c2, c3, d1, d2, d3, t):
    return 2 * ((c1 * np.cos(t) + c2 * np.sin(t) + c3) * (-c1 * np.sin(t) + c2 * np.cos(t)) + \
                (d1 * np.cos(t) + d2 * np.sin(t) + d3) * (-d1 * np.sin(t) + d2 * np.cos(t)))

def g(c1, c2, c3, d1, d2, d3, t):
    return 2 * ((c1 * np.cos(t) + c2 * np.sin(t) + c3) * (-c1 * np.cos(t) - c2 * np.sin(t)) + \
                (-c1 * np.sin(t) + c2 * np.cos(t)) ** 2 + \
                (d1 * np.cos(t) + d2 * np.sin(t) + d3) * (-d1 * np.cos(t) - d2 * np.sin(t)) + \
                (-d1 * np.sin(t) + d2 * np.cos(t)) ** 2)

def newton(f, g, c1, c2, c3, d1, d2, d3, t0, tol=1e-10, maxIter=100):
    cur_iter = 0
    t = t0
    while fabs(f(c1, c2, c3, d1, d2, d3, t)) >= tol and cur_iter < maxIter:
        if g(c1, c2, c3, d1, d2, d3, t) == 0:
            t = t - f(c1, c2, c3, d1, d2, d3, t) * 0.01
        else:
            t = t - f(c1, c2, c3, d1, d2, d3, t) / g(c1, c2, c3, d1, d2, d3, t)
        cur_iter += 1

    t = (t / (2 * pi) - math.floor(t / (2 * pi))) * 2 * pi
    return t

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

