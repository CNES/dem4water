from math import atan, cos, sin

import geopandas as gpd
import numpy as np
import rasterio as rio
import shapely
from shapely.geometry import LineString, Point, box
from shapely.ops import split


def compute_distance(point1, point2, thresholding=None):
    """."""
    if thresholding == "max":
        dist = np.ceil(
            np.sqrt(
                (point1[0] - point2[0]) * (point1[0] - point2[0])
                + (point1[1] - point2[1]) * (point1[1] - point2[1])
            )
        )
    elif thresholding == "min":
        dist = np.floor(
            np.sqrt(
                (point1[0] - point2[0]) * (point1[0] - point2[0])
                + (point1[1] - point2[1]) * (point1[1] - point2[1])
            )
        )
    else:
        dist = np.sqrt(
            (point1[0] - point2[0]) * (point1[0] - point2[0])
            + (point1[1] - point2[1]) * (point1[1] - point2[1])
        )
    return dist


raster = "/home/btardy/Documents/activites/WATER/GDP/Berthier_chain1/masked_mnt.tif"
raster = "/home/btardy/Documents/activites/WATER/GDP/Etang-de-Sault/waterbody_bin.tif"
with rio.open(raster) as image:
    bounds = image.bounds
    geom = box(*bounds)
    gdf_mnt = gpd.GeoDataFrame({"id": 1, "geometry": [geom]}, crs=image.crs)

    print(bounds)
    distance = (
        compute_distance([bounds.top, bounds.left], [bounds.bottom, bounds.right]) + 10
    ) / 2
cutline_base = (
    "/home/btardy/Documents/activites/WATER/GDP/Berthier_chain1/cutline_base.geojson"
)
cutline_base = (
    "/home/btardy/Documents/activites/WATER/GDP/Etang-de-Sault/cutline_base.geojson"
)

gdf = gpd.read_file(cutline_base)
# distance = 300
coords = list(gdf.geometry.values[0].coords)
print(coords)


point_a = Point(coords[0])
point_b = Point(coords[-1])

middle_point = Point((point_a.x + point_b.x) / 2, (point_a.y + point_b.y) / 2)

theta = atan((point_b.y - point_a.y) / (point_b.x - point_a.x))

point_c = Point(
    (middle_point.x - distance * sin(theta)), (middle_point.y + distance * cos(theta))
)
point_d = Point(
    (middle_point.x + distance * sin(theta)), (middle_point.y - distance * cos(theta))
)

line = LineString([point_c, point_d])
gdf.geometry = [line]
gdf.to_file("/home/btardy/Documents/activites/WATER/GDP/test_mediatrice.geojson")

inter = split(geom, line)  #  gdf.overlay(gdf_mnt, how="intersection")
inter = [g for g in inter.geoms]
gdf2 = gpd.GeoDataFrame({"sub": range(len(inter))}, geometry=inter, crs=gdf.crs)
gdf2.to_file("/home/btardy/Documents/activites/WATER/GDP/test_mediatrice_inter.geojson")


def oversampling_polygon_boundary(gdf_wb, max_dist):
    """Ensure that each point of boundary is lower than the max_dist parameter.

    Parameters
    ----------
    gdf_wb:
         geodataframe containing only one polygon
    max_dist:
        the maximum distance in meter allowed between two points
    """
    # ensure projection is in meters
    gdf_w = gdf_wb.to_crs(2154)
    list_of_points = list(gdf_w.geometry.values[0].exterior.coords)
    new_list_of_points = []
    while len(list_of_points) >= 2:
        ori_point = list_of_points.pop(0)
        new_list_of_points.append(ori_point)

        target_point = list_of_points.pop(0)
        line = LineString([ori_point, target_point])
        if max_dist < line.length:
            # Too spaced generate point
            new_point = line.interpolate(max_dist)
            list_of_points = [new_point, target_point] + list_of_points
        else:
            # Close enough
            list_of_points = [target_point] + list_of_points
    gdf_w.geometry = [LineString(new_list_of_points)]
    gdf_w = gdf_w.to_crs(gdf_wb.crs)
    return gdf_w


def find_perpendicular_bisector(cutline_base, bisector_lenght, bisector_file):
    """Draw the perpendicular bisector of a unique cutline.

    Parameters
    ----------
    cutline_base:
        geojson or equivalent file containing the cutline
    bisector_lenght:
        length in meter from the cutline to the point ending the bisector
        on a side, i.e
    """
    gdf = gpd.read_file(cutline_base)
    coords = list(gdf.geometry.values[0].coords)
    point_a = Point(coords[0])
    point_b = Point(coords[-1])
    middle_point = Point((point_a.x + point_b.x) / 2, (point_a.y + point_b.y) / 2)
    theta = atan((point_b.y - point_a.y) / (point_b.x - point_a.x))
    point_c = Point(
        (middle_point.x - bisector_lenght * sin(theta)),
        (middle_point.y + bisector_lenght * cos(theta)),
    )
    point_d = Point(
        (middle_point.x + bisector_lenght * sin(theta)),
        (middle_point.y - bisector_lenght * cos(theta)),
    )

    line = LineString([point_c, point_d])
    gdf.geometry = [line]
    gdf.to_file(bisector_file)
