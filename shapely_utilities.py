
from shapely.geometry import Point, Polygon
from shapely.ops import cascaded_union
import shapely.affinity

import math

# Based upon:
# Drawing ellipse with shapely?
# https://gis.stackexchange.com/questions/243459/drawing-ellipse-with-shapely
# Draw an ellipse using Shapely:
# https://stackoverflow.com/questions/13105915/draw-an-ellipse-using-shapely
# The angle passed in is measured counter-clockwise from the x-axis, not clockwise
def define_points_for_ellipse(x, y, a, b, angle):
    #print(str(x) + " " + str(y) + " " + str(a) + " " + str(b) + " " + str(angle))
    circ = shapely.geometry.Point(x, y).buffer(1)
    ell = shapely.affinity.scale(circ, a, b)
    ellr = shapely.affinity.rotate(ell, angle, use_radians=True)
    #return (list(ellr.exterior.coords.xy))
    return list(ellr.exterior.coords)

def define_polygon_for_ellipse(x, y, a, b, angle):
    # 1st elem = center point (x,y) coordinates
    # Let create a circle of radius 1 around center point:
    circ = shapely.geometry.Point(x, y).buffer(1)

    # 2nd elem = the two semi-axis values (along x, along y)
    # Let create the ellipse along x and y:
    ell = shapely.affinity.scale(circ, a, b)

    # 3rd elem = angle in degrees between x-axis of the Cartesian base and the corresponding semi-axis
    # Let rotate the ellipse (clockwise, x axis pointing right):
    # The angle of rotation can be specified in either degrees (default) or radians by setting use_radians=True.
    # Positive angles are counter-clockwise and negative are clockwise rotations.
    ellr = shapely.affinity.rotate(ell, angle, use_radians=True)
    print("Ellipse Area: " + str(ellr.area))
    coords =  list(ellr.exterior.coords)
    polygon = Polygon(coords)
    print("Polygon Ellipse Area: " + str(polygon.area))

def dissolve_ellipses(ellipse_list):
    polygons = []
    for row_idx in range(0, len(ellipse_list['x'])):
        coords = define_points_for_ellipse(ellipse_list['x'][row_idx], ellipse_list['y'][row_idx],
                                                     ellipse_list['a'][row_idx], ellipse_list['b'][row_idx],
                                                     ellipse_list['angle'][row_idx])
        polygon = Polygon(coords)
        polygons.append(polygon)

    dissolve = cascaded_union(polygons)

    new_data = dict()
    new_data['x'] = dissolve.exterior.coords.xy[0].tolist()
    new_data['y'] = dissolve.exterior.coords.xy[1].tolist()
    new_data['area'] = dissolve.area
    return new_data
    #new_data['x'], new_data['y'] = list(dissolve.exterior.coords.xy)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def main():
    define_polygon_for_ellipse(0, 0, math.sqrt(1.0 / math.pi), math.sqrt(1.0 / math.pi), 0)

if __name__ == '__main__':
    main()