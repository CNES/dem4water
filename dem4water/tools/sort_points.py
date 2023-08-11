import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point


def sort_coordinates(list_of_xy_coords):
    cx, cy = list_of_xy_coords.mean(0)
    x, y = list_of_xy_coords.T
    angles = np.arctan2(x - cx, y - cy)
    indices = np.argsort(angles)
    return list_of_xy_coords[indices]


in_points = "/home/btardy/Documents/activites/WATER/GDP/tests_bds/export_points_to_be_sorted.geojson"
out_file = "/home/btardy/Documents/activites/WATER/GDP/tests_bds/test_sort_points_lines.geojson"

algo = "inpe"
in_points = f"/home/btardy/Documents/activites/WATER/GDP/test_python/{algo}/test_intersect_lines.geojson"
out_file = f"/home/btardy/Documents/activites/WATER/GDP/test_python/{algo}/test_intersect_lines_sort.geojson"

# in_points = (
#     "/home/btardy/Documents/activites/WATER/GDP/test_python/test_poly_to_points.geojson"
# )
# out_file = "/home/btardy/Documents/activites/WATER/GDP/test_python/test_poly_to_points_sort.geojson"


def sort_and_filter_points(in_points, out_file):
    gdf = gpd.GeoDataFrame().from_file(in_points)
    print(gdf)

    unique_object = gdf.ident.unique()
    # {"ident": list(gdf.ident)}
    df = pd.DataFrame()

    print(unique_object)
    ident_list = []
    coords_x_list = []
    coords_y_list = []
    for uni in unique_object:
        vals = gdf.loc[gdf.ident == uni]
        # ident = uni * np.ones(len(vals.index))
        # print(ident)
        xxx = vals.geometry.x
        yyy = vals.geometry.y
        coords = np.array([(x, y) for x, y in zip(xxx, yyy)])
        sorted_c = sort_coordinates(coords)
        coords_x = [sorted_c[0][0], sorted_c[-1][0]]
        coords_y = [sorted_c[0][1], sorted_c[-1][1]]
        # coords_x = [x[0] for x in sorted_c]
        # coords_y = [x[1] for x in sorted_c]
        coords_x_list += coords_x
        coords_y_list += coords_y
        ident_list += [uni, uni]

    df["x"] = coords_x_list
    df["y"] = coords_y_list
    df["ident"] = ident_list
    df["idx"] = range(len(df.index))
    print(df)
    gdf2 = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.x, df.y), crs=gdf.crs)
    gdf2.to_file(out_file.replace(".geojson", "points.geojson"))
    lines = gdf2.groupby(["ident"])["geometry"].apply(lambda x: LineString(x.tolist()))

    # store as a GeodataFrame and add 'ID' as a column (currently stored as the 'index')
    lines = gpd.GeoDataFrame(lines, geometry="geometry", crs=gdf.crs)
    lines.reset_index(inplace=True)
    lines["length"] = lines.geometry.length
    # lines = lines.loc[lines["length"] > 20]
    print(lines)
    lines.to_file(out_file)


sort_and_filter_points(in_points, out_file)
