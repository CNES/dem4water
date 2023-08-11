import os

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString


def poly_to_points(feature, poly):
    """."""
    return {feature: poly.exterior.coords}


def convert_geodataframe_poly_to_points(in_gdf, id_col, ident_gdp):
    """."""
    in_gdf["points"] = in_gdf.apply(
        lambda p: poly_to_points(p[id_col], p["geometry"]), axis=1
    )
    list_points = list(in_gdf.points)

    ident = []
    coordx = []
    coordy = []
    fidx = []
    for i, poly in enumerate(list_points):
        for _, points in poly.items():
            for idx, point in enumerate(points[:-1]):
                fidx.append(idx)
                ident.append(ident_gdp)
                coordx.append(point[0])
                coordy.append(point[1])
    df_coords = pd.DataFrame()
    df_coords["x"] = coordx
    df_coords["y"] = coordy
    df_coords[id_col] = ident
    df_coords["gdp_point"] = fidx  # range(len(df_coords.index))

    gdf1 = gpd.GeoDataFrame(
        df_coords, geometry=gpd.points_from_xy(df_coords.x, df_coords.y), crs=in_gdf.crs
    )
    return gdf1


def sort_coordinates(list_of_xy_coords):
    """."""
    cx, cy, _ = list_of_xy_coords.mean(0)
    x, y, _ = list_of_xy_coords.T
    angles = np.arctan2(x - cx, y - cy)
    indices = np.argsort(angles)
    return list_of_xy_coords[indices]


# xxx = [0, 4, 1, 2, 3]
# yyy = [0, 4, 1, 2, 3]
# zzz = [0, 1, 2, 3, 4]
# list_of_xy_coords = np.array(list(zip(xxx, yyy, zzz)))

# list_of_xy_coords = sort_coordinates(list_of_xy_coords)
# input(list_of_xy_coords)
gdf_gdp = gpd.GeoDataFrame().from_file(
    "/home/btardy/Documents/activites/WATER/GDP/test_chain_2/filtered_poly.geojson"
)

gdf_wb = gpd.GeoDataFrame().from_file(
    "/home/btardy/Documents/activites/WATER/GDP/test_chain_2/contour_points.geojson"
)

gdf = gpd.sjoin_nearest(gdf_gdp, gdf_wb)
# gdf.to_file("/home/btardy/Documents/activites/WATER/GDP/test_chain_2/inter.geojson")
print(gdf_wb)
print(gdf)
# input("wait")

list_df = []
for ident in gdf.gdp_unique_id.unique():
    gdf_work = gdf[gdf.gdp_unique_id == ident]
    gdf_work.iloc[[0]].to_file(
        "/home/btardy/Documents/activites/WATER/GDP/test_chain_2/gdf_work.geojson"
    )
    if len(gdf_work) == 1:
        print(f"No enought points to start cutline search for gdp n°{ident}")
        continue
    print(gdf_work)
    wb_ident = gdf_work.wb_unique_id.unique()
    if len(wb_ident) > 1:
        print("Impossible to match the gdp with a unique water body.")
    else:
        wb_ident = wb_ident[0]
    max_indices_wb = gdf_wb[gdf_wb.wb_unique_id == wb_ident].id_point.max()
    min_indices_gdp = gdf_work.id_point.min()
    max_indices_gdp = gdf_work.id_point.max()
    if min_indices_gdp < 100 and max_indices_gdp > max_indices_wb - 100:
        print("passe par l'origine")
        gdf_points = convert_geodataframe_poly_to_points(
            gdf_work.iloc[[0]], "gdp_unique_id", ident
        )
        print("*" * 10)
        print(gdf_points)
        gdf_points_contour = gpd.GeoDataFrame(
            gdf_work[["wb_unique_id", "id_point"]],
            geometry=gpd.points_from_xy(gdf_work.x, gdf_work.y),
            crs=gdf_work.crs,
        )
        inter = gpd.sjoin_nearest(gdf_points_contour, gdf_points)
        gdf_temp = inter.sort_values("gdp_point")
        gdf_work = gdf_temp.reset_index(drop=True)
        gdf_temp = gpd.GeoDataFrame(
            gdf_temp[["gdp_unique_id", "wb_unique_id", "id_point"]],
            geometry=gpd.points_from_xy(gdf_temp.x, gdf_temp.y),
            crs=gdf_work.crs,
        )
        # print("*" * 10)
        # print(inter)
        # print("*" * 10)
        # print(list(inter.id_point))
        # inter.to_file(
        #     "/home/btardy/Documents/activites/WATER/GDP/test_chain_2/gdf_work_1.geojson"
        # )
        # input("wait")
        list_df.append(gdf_temp)
    else:
        print("no zero")
        gdf_work = gdf_work.sort_values("id_point")
        gdf_work = gdf_work.reset_index(drop=True)
        gdf_temp = gpd.GeoDataFrame(
            gdf_work[["gdp_unique_id", "wb_unique_id", "id_point"]],
            geometry=gpd.points_from_xy(gdf_work.x, gdf_work.y),
            crs=gdf_work.crs,
        )
        list_df.append(gdf_temp)
gdf_n = gpd.GeoDataFrame(pd.concat(list_df, ignore_index=True), crs=list_df[0].crs)
lines = gdf_n.groupby(["gdp_unique_id"])["geometry"].apply(
    lambda x: LineString(x.tolist())
)
lines = gpd.GeoDataFrame(lines, geometry="geometry", crs=gdf.crs)
lines.reset_index(inplace=True)
work_dir = "/home/btardy/Documents/activites/WATER/GDP/test_chain_2/"
lines.to_file(os.path.join(work_dir, "cutline.geojson"))
