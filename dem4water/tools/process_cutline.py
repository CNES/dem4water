import os

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from shapely.geometry import LineString

from dem4water.tools.compute_grandient_dot_product import compute_gradient_product
from dem4water.tools.polygonize_raster import PolygonizeParams, polygonize
from dem4water.tools.rasterize_vectors import RasterizarionParams, rasterize
from dem4water.tools.remove_holes_in_shapes import close_holes


# ######################################
# Preprocess database
# #######################################
def remove_mutlipolygon(in_vector, epsg, buffer_size=10):
    """Apply a buffer in both sign to close small gap between polygons."""
    gdf = gpd.GeoDataFrame().from_file(in_vector)
    gdf = gdf.to_crs(epsg)
    gdf["geometry"] = gdf.geometry.buffer(buffer_size)
    gdf["geometry"] = gdf.geometry.buffer(-buffer_size)
    gdf = gdf.explode(ignore_index=True)
    return gdf


def preprocess_water_body(in_vector, epsg):
    """."""
    # Fuse multipolygon
    gdf_wb = remove_mutlipolygon(in_vector, epsg, 10)
    # Remove holes
    gdf_wb.geometry = gdf_wb.geometry.apply(lambda p: close_holes(p))

    return gdf_wb


def clear_polygonize(in_vector, epsg):
    """."""
    gdf = gpd.GeoDataFrame().from_file(in_vector)
    gdf = gdf.to_crs(epsg)
    gdf = gdf.loc[gdf["DN"] == 1]
    return gdf


def filter_by_convex_hull(gdf, work_dir, id_col):
    """."""
    gdf_convex_hull = gdf.copy()
    gdf_convex_hull.geometry = gdf_convex_hull.geometry.convex_hull
    gdf_convex_hull["area"] = gdf_convex_hull.geometry.area
    gdf_convex_hull = gdf_convex_hull.sort_values("area", ascending=False)
    gdf_convex_hull = gdf_convex_hull.reset_index(drop=True)

    id_line = list(gdf_convex_hull[id_col])
    view_id = []
    not_to_see = []
    for line in id_line:
        if line not in not_to_see:
            df_work = gdf_convex_hull.loc[gdf_convex_hull[id_col] == line]
            intersect = gpd.sjoin(gdf, df_work, how="inner")
            view_id.append(line)
            not_to_see += list(intersect[f"{id_col}_left"])

            # print(intersect)
            # intersect.to_file(os.path.join(work_dir, "intersect.geojson"))
    gdf_filtered = gdf[gdf[id_col].isin(view_id)]
    # gdf_filtered = gdf_filtered.drop(["area"], axis=1)
    gdf_filtered.to_file(os.path.join(work_dir, "filtered_poly.geojson"))
    return gdf_filtered


# ##################################################
# Handle points
# ##################################################


def poly_to_points(feature, poly):
    """."""
    return {feature: poly.exterior.coords}


def convert_geodataframe_poly_to_points(in_gdf, id_col, out_id_name, ident_gdp=None):
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
                if ident_gdp is None:
                    ident.append(i)
                else:
                    ident.append(ident_gdp)
                coordx.append(point[0])
                coordy.append(point[1])
    df_coords = pd.DataFrame()
    df_coords["x"] = coordx
    df_coords["y"] = coordy
    df_coords[id_col] = ident
    df_coords[out_id_name] = fidx  # range(len(df_coords.index))
    gdf = gpd.GeoDataFrame(
        df_coords, geometry=gpd.points_from_xy(df_coords.x, df_coords.y), crs=in_gdf.crs
    )
    return gdf


def draw_lines(gdf_wb_points, gdf_gdp_poly):
    gdf = gpd.sjoin_nearest(gdf_gdp_poly, gdf_wb_points)
    list_df = []
    for ident in gdf.gdp_unique_id.unique():
        gdf_work = gdf[gdf.gdp_unique_id == ident]
        if len(gdf_work) == 1:
            print(f"No enought points to start cutline search for gdp n°{ident}")
            continue
        wb_ident = gdf_work.wb_unique_id.unique()
        if len(wb_ident) > 1:
            print("Impossible to match the gdp with a unique water body.")
        else:
            wb_ident = wb_ident[0]
            max_indices_wb = gdf_wb_points[
                gdf_wb_points.wb_unique_id == wb_ident
            ].id_point.max()
        min_indices_gdp = gdf_work.id_point.min()
        max_indices_gdp = gdf_work.id_point.max()
        if min_indices_gdp < 100 and max_indices_gdp > max_indices_wb - 100:
            print("passe par l'origine")
            gdf_points = convert_geodataframe_poly_to_points(
                gdf_work.iloc[[0]], "gdp_unique_id", "gdp_point", ident
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
    # Remove useless points by simplying each line
    lines.geometry = lines.simplify(1)
    lines.reset_index(inplace=True)
    return lines


# ##############################################
# Find higher altitude to complete cutline
# ##############################################


def search_line(init_point, gdf_lines, mnt_raster, buffer_size):
    """."""
    gdf = gpd.GeoDataFrame(
        {"init_point": 1},
        geometry=gpd.points_from_xy(init_point[0], init_point[1]),
        crs=gdf_lines.crs,
    )
    gdf.geometry = gdf.geometry.buffer(buffer_size)
    geoms = gdf.geometry.values
    with rasterio.open(mnt_raster) as mnt:
        out_image, out_transform = mask(mnt, geoms, crop=True)
        no_data = mnt.nodata
        data = out_image[0, :, :]
        row, col = np.where(data != no_data)
        vals = np.extract(data != no_data, data)
        x_coords, y_coords = rasterio.transform.xy(
            out_transform, row, cols, offset="center"
        )
        gdf = gpd.GeoDataFrame(
            {"row": row, "col": col, "atlitude": vals},
            geometry=gpd.points_from_xy(x_coords, y_coords),
            crs=gdf_lines.crs,
        )


def fill_cutline(gdf_line):
    """."""

    # Extraire chaque ligne du dataframe
    for i in gdf_line.index:
        line = list(gdf_line.iloc[i].geometry.coords)
        left_point = line[0]
        right_point = line[-1]
    # Prendre les 2 points extrème de chaque cotés

    # Chercher dans un radius de N le point le plus haut

    # Assurer qu'on ne revient pas en arrière


# ##############################################
# Process
# ##############################################
def prepare_inputs(in_vector, mnt_raster, work_dir, gdp_buffer_size):
    """."""
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)
    with rasterio.open(mnt_raster) as src:
        epsg = src.crs.to_epsg()
    # 1 Handle mutlipolygons
    # 2 remove holes
    gdf_bd = preprocess_water_body(in_vector, epsg)
    gdf_bd.to_file(os.path.join(work_dir, "bd_clean.geojson"))
    # 4 Rasterize the waterbody filtered
    params_raster = RasterizarionParams(
        mode="binary",
        binary_foreground_value=1,
        background_value=0,
        column_field=None,
        dtype="uint8",
    )
    waterbody_bin = os.path.join(work_dir, "waterbody_bin.tif")
    rasterize(gdf_bd, mnt_raster, waterbody_bin, params_raster)
    # 5 run GDP
    gdp_raster = os.path.join(work_dir, "gdp_raster.tif")
    compute_gradient_product(waterbody_bin, mnt_raster, gdp_raster)
    # 6 Vectorize
    params_poly = PolygonizeParams(
        layer_name="output", field_name="DN", driver="GeoJson", overwrite=True
    )
    gdp_vector = os.path.join(work_dir, "gdp_vector.geojson")
    polygonize(gdp_raster, gdp_vector, params_poly)
    gdf_gdp = clear_polygonize(gdp_vector, epsg)
    # 7 Relier les polygones GDP proches
    # TODO compute area w.r.t the waterbody whole area
    gdf_gdp = gdf_gdp.loc[gdf_gdp.geometry.area > 2000]
    # find intersection according the max distance
    gdf_gdp.geometry = gdf_gdp.geometry.buffer(gdp_buffer_size)
    # then fuse
    gdf_gdp = gdf_gdp.dissolve(by="DN")
    # remove the buffer
    gdf_gdp.geometry = gdf_gdp.geometry.buffer(-gdp_buffer_size)
    # then separe into unique object
    gdf_gdp = gdf_gdp.explode(ignore_index=True)
    # dissolve remove the DN column
    gdf_gdp["gdp_unique_id"] = range(len(gdf_gdp.index))
    gdf_gdp.to_file(os.path.join(work_dir, "gdp_fusion_nearest.geojson"))

    gdf_gdp = filter_by_convex_hull(gdf_gdp, work_dir, "gdp_unique_id")
    # 8 Convertir la BD en points
    gdf_bd["wb_unique_id"] = range(len(gdf_bd.index))
    gdf_wb_points = convert_geodataframe_poly_to_points(
        gdf_bd, "wb_unique_id", "id_point"
    )
    gdf_wb_points.to_file(os.path.join(work_dir, "contour_points.geojson"))

    # 9 Draw line
    lines = draw_lines(gdf_wb_points, gdf_gdp)
    lines.to_file(os.path.join(work_dir, "cutline.geojson"))


prepare_inputs(
    "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/extract_marne_inpe_noholes.geojson",
    "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/dem_extract_Marne-Giffaumont.tif",
    "/home/btardy/Documents/activites/WATER/GDP/Marne_chain",
    100,
)

prepare_inputs(
    "/home/btardy/Documents/activites/WATER/GDP/extract/Laparan/laparan_bd.geojson",
    "/home/btardy/Documents/activites/WATER/GDP/extract/Laparan/dem_extract_Laparan.tif",
    "/home/btardy/Documents/activites/WATER/GDP/Laparan_chain",
    100,
)
prepare_inputs(
    "/home/btardy/Documents/activites/WATER/GDP/extract/Vaufrey/vaufrey_bd.geojson",
    "/home/btardy/Documents/activites/WATER/GDP/extract/Vaufrey/dem_extract_Vaufrey.tif",
    "/home/btardy/Documents/activites/WATER/GDP/Vaufrey_chain",
    100,
)
prepare_inputs(
    "/home/btardy/Documents/activites/WATER/GDP/extract/Naussac/naussac_bd.geojson",
    "/home/btardy/Documents/activites/WATER/GDP/extract/Naussac/dem_extract_Naussac.tif",
    "/home/btardy/Documents/activites/WATER/GDP/Naussac_chain",
    100,
)
