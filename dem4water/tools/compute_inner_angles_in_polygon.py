import math

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Point


def getAngle(a, b, c):
    ang = math.degrees(
        math.atan2(c[1] - b[1], c[0] - b[0]) - math.atan2(a[1] - b[1], a[0] - b[0])
    )
    return ang


def get_angle(a, b, c):
    ab = math.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]))
    ac = math.sqrt((a[0] - c[0]) * (a[0] - c[0]) + (a[1] - c[1]) * (a[1] - c[1]))
    cb = math.sqrt((c[0] - b[0]) * (c[0] - b[0]) + (c[1] - b[1]) * (c[1] - b[1]))
    ang = math.acos((ac * ac + ab * ab - cb * cb) / (2 * ac * ab))
    return np.degrees(ang)


gdf = gpd.read_file(
    "/home/btardy/Documents/activites/WATER/GDP/extract/Berthier/berthier_bd.geojson"
)
#    "/home/btardy/Documents/activites/WATER/GDP/bertier_new_cut/bd_clean_simpl.geojson"
gdf = gdf.to_crs(2154)
gdf_simpl = gpd.read_file(
    "/home/btardy/Documents/activites/WATER/GDP/extract/Berthier/berthier_bd.geojson"
)
gdf_simpl = gdf_simpl.to_crs(2154)
gdf_simpl.geometry = gdf_simpl.geometry.simplify(1)
coords_points = gdf_simpl.geometry.values[0].exterior.coords
simple_coords = list(coords_points)[:-1]

list_a = simple_coords[:]
list_b = simple_coords[1:] + [simple_coords[0]]
list_c = [simple_coords[-1]] + simple_coords[:-1]

print(len(list_c), len(list_a), len(list_b))
list_angle = []
list_rupt = []
for a, b, c in zip(list_a, list_b, list_c):
    angle = get_angle(a, c, b)
    print(a, b, c, angle)
    list_angle.append(angle)
    if angle < 150:
        list_rupt.append(a)
    # ac = np.array(list_c) - np.array(list_a)
    # ab = np.array(list_b) - np.array(list_a)
    # print(ac)
    # print(ab)
    # cosine_angle = np.dot(ac, ab) / (np.linalg.norm(ac) * np.linalg.norm(ab))
    # angle = np.arccos(cosine_angle)

    # print(np.degrees(angle))
# print(min(list_angle))
# print(max(list_angle))
# points_x = [x[0] for x in list_a]
# points_y = [x[1] for x in list_a]
# points = [Point(x, y) for x, y in zip(points_x, points_y)]
# gdf = gpd.GeoDataFrame(
#     {"point_x": points_x, "point_y": points_y, "angle": list_angle},
#     geometry=points,
#     crs=2154,
# )
# gdf.to_file(
#     "/home/btardy/Documents/activites/WATER/GDP/bertier_new_cut/Points_angle.geojson"
# )
coords_points = list(gdf.geometry.values[0].exterior.coords)
segments = []
seg = []
for point in coords_points:
    seg.append(point)
    if point in list_rupt:
        segments.append(seg)
        seg = []

lines = []
for seg in segments:
    if len(seg) > 2:
        lines.append(LineString(seg))

gdf = gpd.GeoDataFrame(
    {"segments": range(len(lines))},
    geometry=lines,
    crs=2154,
)
gdf.to_file(
    "/home/btardy/Documents/activites/WATER/GDP/bertier_new_cut/segments_test.geojson"
)


gdf_gdp = gpd.read_file(
    "/home/btardy/Documents/activites/WATER/GDP/bertier_new_cut/gdp_vector.geojson"
)
gdf_gdp = gdf_gdp.to_crs(2154)
inter = gpd.sjoin(gdf, gdf_gdp, how="inner")
print(inter)
inter.to_file(
    "/home/btardy/Documents/activites/WATER/GDP/bertier_new_cut/base_cutline.geojson"
)
