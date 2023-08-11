import os

# import alphashape
# import cartopy.crs as ccrs
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from shapely.geometry import LineString

from dem4water.tools.compute_grandient_dot_product import \
    compute_gradient_product
# from dem4water.tools.concave_hull import concave_hull
from dem4water.tools.polygonize_raster import PolygonizeParams, polygonize
from dem4water.tools.rasterize_vectors import RasterizarionParams, rasterize
from dem4water.tools.remove_holes_in_shapes import close_holes


def sort_coordinates(list_of_xy_coords):
    """."""
    cx, cy = list_of_xy_coords.mean(0)
    x, y = list_of_xy_coords.T
    angles = np.arctan2(x - cx, y - cy)
    indices = np.argsort(angles)
    return list_of_xy_coords[indices]


def sort_and_filter_points(gdf):
    """."""
    # gdf = gpd.GeoDataFrame().from_file(in_points)
    # print(gdf)

    unique_object = gdf.ident.unique()
    # {"ident": list(gdf.ident)}
    df_coords = pd.DataFrame()

    # print(unique_object)
    ident_list = []
    coords_x_list = []
    coords_y_list = []
    for uni in unique_object:
        vals = gdf.loc[gdf.ident == uni]
        # ident = uni * np.ones(len(vals.index))
        # print(ident)
        xxx = vals.geometry.x
        yyy = vals.geometry.y
        coords = np.array(list(zip(xxx, yyy)))
        sorted_c = sort_coordinates(coords)
        index = int(len(coords) / 2)
        coords_x = [sorted_c[0][0], sorted_c[index][0]]
        coords_y = [sorted_c[0][1], sorted_c[index][1]]
        # coords_x = [x[0] for x in sorted_c]
        # coords_y = [x[1] for x in sorted_c]
        coords_x_list += coords_x
        coords_y_list += coords_y
        ident_list += [uni, uni]

    df_coords["x"] = coords_x_list
    df_coords["y"] = coords_y_list
    df_coords["ident"] = ident_list
    df_coords["idx"] = range(len(df_coords.index))
    # print(df_coords)
    gdf2 = gpd.GeoDataFrame(
        df_coords, geometry=gpd.points_from_xy(df_coords.x, df_coords.y), crs=gdf.crs
    )
    return gdf2


def sort_points(gdf):
    unique_object = gdf.ident.unique()
    # {"ident": list(gdf.ident)}
    df_coords = pd.DataFrame()

    # print(unique_object)
    ident_list = []
    coords_x_list = []
    coords_y_list = []
    fidx = []
    for uni in unique_object:
        vals = gdf.loc[gdf.ident == uni]
        # ident = uni * np.ones(len(vals.index))
        # print(ident)
        xxx = vals.geometry.x
        yyy = vals.geometry.y
        coords = np.array(list(zip(xxx, yyy)))
        sorted_c = sort_coordinates(coords)
        coords_x = list(sorted_c[:, 0])
        coords_y = list(sorted_c[:, 1])
        ident_list += [uni] * len(coords_x)
        coords_x_list += coords_x
        coords_y_list += coords_y
        fidx += range(len(coords_x))
    df_coords["x"] = coords_x_list
    df_coords["y"] = coords_y_list
    df_coords["ident"] = ident_list
    df_coords["idx"] = range(len(df_coords.index))
    df_coords["ident_uni"] = fidx
    # print(df_coords)
    gdf2 = gpd.GeoDataFrame(
        df_coords, geometry=gpd.points_from_xy(df_coords.x, df_coords.y), crs=gdf.crs
    )
    return gdf2


def extract_points_of_interest_from_poly(gdf):
    """."""
    unique_object = gdf.ident.unique()
    df_coords = pd.DataFrame()
    ident_list = []
    coords_x_list = []
    coords_y_list = []
    for uni in unique_object:
        vals = gdf.loc[gdf.ident == uni]
        xxx = vals.geometry.x
        yyy = vals.geometry.y
        coords = np.array(list(zip(xxx, yyy)))
        index = int(len(coords) / 2)
        coords_x = [coords[0][0], coords[index][0]]
        coords_y = [coords[0][1], coords[index][1]]
        coords_x_list += coords_x
        coords_y_list += coords_y
        ident_list += [uni, uni]
    df_coords["x"] = coords_x_list
    df_coords["y"] = coords_y_list
    df_coords["ident"] = ident_list
    df_coords["idx"] = range(len(df_coords.index))
    gdf2 = gpd.GeoDataFrame(
        df_coords, geometry=gpd.points_from_xy(df_coords.x, df_coords.y), crs=gdf.crs
    )
    return gdf2


def poly_to_points(feature, poly):
    """."""
    return {feature: poly.exterior.coords}


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
    gdf = gpd.GeoDataFrame().from_file(in_vector)
    gdf = gdf.to_crs(epsg)
    gdf = gdf.loc[gdf["DN"] == 1]
    return gdf


def convert_geodataframe_poly_to_points(in_gdf, col):
    """."""
    in_gdf["points"] = in_gdf.apply(
        lambda p: poly_to_points(p[col], p["geometry"]), axis=1
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
                ident.append(i)
                coordx.append(point[0])
                coordy.append(point[1])
    df_coords = pd.DataFrame()
    df_coords["x"] = coordx
    df_coords["y"] = coordy
    df_coords["ident"] = ident
    df_coords["fidx"] = fidx  # range(len(df_coords.index))

    gdf = gpd.GeoDataFrame(
        df_coords, geometry=gpd.points_from_xy(df_coords.x, df_coords.y), crs=in_gdf.crs
    )
    return gdf


def find_nearest_points(gdf_wb, gdf_gdp):
    gdf_n = []
    for i in range(len(gdf_gdp.index)):
        df_n = gpd.sjoin_nearest(gdf_wb, gdf_gdp.iloc[[i]], distance_col="distance")
        # print(df_n)
        df_n = df_n.loc[df_n["distance"] == min(list(df_n["distance"]))]
        gdf_n.append(df_n)
        # print("\n")
        # print(df_n)
        # for p in list(df_n.geometry):
        #     print(p)
        # dist = list(df_n["distance"])
        # print(min(dist))
    gdf_n = gpd.GeoDataFrame(pd.concat(gdf_n, ignore_index=True), crs=gdf_n[0].crs)
    gdf_n = gdf_n.sort_values("fidx", ascending=True)
    gdf_n = gdf_n.reset_index(drop=True)
    return gdf_n


def get_lines_per_gdp_points(gdf_points_selected, gdf_wb_contour):
    """."""
    uni_lines = gdf_points_selected["ident_right"].unique()
    list_of_dataframe = []
    for line in uni_lines:
        gdf_work = gdf_points_selected.loc[gdf_points_selected["ident_right"] == line]
        # print(gdf_work)
        fidx_id = list(gdf_work["fidx"])
        print(fidx_id[0], fidx_id[-1] + 1)
        id_selected = range(fidx_id[0], fidx_id[-1] + 1)
        gdf_filter = gdf_wb_contour[gdf_wb_contour["fidx"].isin(id_selected)]
        gdf_filter["line"] = line
        if gdf_filter.shape[0] > 1:
            list_of_dataframe.append(gdf_filter)
        # print(gdf_filter)
    gdf_final = gpd.GeoDataFrame(
        pd.concat(list_of_dataframe, ignore_index=True), crs=list_of_dataframe[0].crs
    )
    return gdf_final


def filter_by_convex_hull(gdf, work_dir):
    """."""
    gdf_convex_hull = gdf.copy()
    gdf_convex_hull.geometry = gdf_convex_hull.geometry.convex_hull
    gdf_convex_hull["area"] = gdf_convex_hull.geometry.area
    gdf_convex_hull = gdf_convex_hull.sort_values("area", ascending=False)
    gdf_convex_hull = gdf_convex_hull.reset_index(drop=True)

    id_line = list(gdf_convex_hull["uni_id"])
    view_id = []
    not_to_see = []
    for line in id_line:
        if line not in not_to_see:
            df_work = gdf_convex_hull.loc[gdf_convex_hull["uni_id"] == line]
            intersect = gpd.sjoin(gdf, df_work, how="inner")
            view_id.append(line)
            not_to_see += list(intersect["uni_id_left"])

            # print(intersect)
            intersect.to_file(os.path.join(work_dir, "intersect.geojson"))
    gdf_filtered = gdf[gdf["uni_id"].isin(view_id)]
    # gdf_filtered = gdf_filtered.drop(["area"], axis=1)
    gdf_filtered.to_file(os.path.join(work_dir, "filtered_poly.geojson"))
    return gdf_filtered


def compute_cutline_from_gdp_poly(gdf_wb_contour, gdf_gdp):
    inter = gpd.sjoin_nearest(gdf_gdp, gdf_wb_contour)
    uni_poly = gdf_gdp["uni_id"].unique()
    ident = []
    coordx = []
    coordy = []
    idx = []
    fidx = []
    for id_poly in uni_poly:
        df_work = inter.loc[inter["uni_id"] == id_poly]
        if len(df_work.index) > 1:
            ident += [id_poly] * len(df_work.index)
            idx += range(len(df_work.index))
            coordx += list(df_work.x)
            coordy += list(df_work.y)
            fidx += list(df_work.fidx)
    df_coords = pd.DataFrame()
    df_coords["x"] = coordx
    df_coords["y"] = coordy
    df_coords["ident"] = ident
    df_coords["idx"] = idx
    df_coords["fidx"] = fidx
    gdf = gpd.GeoDataFrame(
        df_coords,
        geometry=gpd.points_from_xy(df_coords.x, df_coords.y),
        crs=gdf_gdp.crs,
    )
    gdf.to_file(
        "/home/btardy/Documents/activites/WATER/GDP/test_chain/filterd_points.geojson"
    )
    gdf = gdf.drop_duplicates("geometry")
    gdf = sort_points(gdf)  # fous le bordel
    gdf.to_file(
        "/home/btardy/Documents/activites/WATER/GDP/test_chain/filterd_points_no_dup.geojson"
    )
    # gdf = concave_hull(gdf)
    # gdf.geometry = gdf.geometry.concave_hull()
    # gdf_proj = gdf.to_crs(ccrs.AlbersEqualArea().proj4_init)
    # alpha_shape = alphashape.alphashape(gdf_proj)
    # print(type(alpha_shape))
    gdf.to_file(
        "/home/btardy/Documents/activites/WATER/GDP/test_chain/concave_hull.geojson"
    )
    lines = gdf.groupby(["ident"])["geometry"].apply(lambda x: LineString(x.tolist()))
    lines = gpd.GeoDataFrame(lines, geometry="geometry", crs=gdf.crs)
    lines.reset_index(inplace=True)
    lines["length"] = lines.geometry.length
    return lines


def prepare_inputs(in_vector, mnt_raster, work_dir, gdp_buffer_size):
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
    gdf_gdp["uni_id"] = range(len(gdf_gdp.index))
    gdf_gdp.to_file(os.path.join(work_dir, "gdp_fusion_nearest.geojson"))

    gdf_gdp = filter_by_convex_hull(gdf_gdp, work_dir)

    # 8 Extract points from GDP
    gdf_gdp1 = convert_geodataframe_poly_to_points(gdf_gdp, "uni_id")
    # gdf_gdp_points = extract_points_of_interest_from_poly(gdf_gdp)
    # gdf_gdp_points = sort_and_filter_points(gdf_gdp_points)
    gdf_gdp1.to_file(os.path.join(work_dir, "sorted_points.geojson"))

    # 9 Follow contour between points A&B to draw cutline
    gdf_bd["ident"] = range(len(gdf_bd.index))
    gdf_wb_points = convert_geodataframe_poly_to_points(gdf_bd, "id")
    gdf_wb_points.to_file(os.path.join(work_dir, "contour_points.geojson"))
    # print(gdf_gdp_points.columns)
    # gdf_nearest_points = find_nearest_points(gdf_wb_points, gdf_gdp_points)
    # print(gdf_nearest_points)
    # gdf_nearest_points.to_file(os.path.join(work_dir, "nearest_points.geojson"))
    # gdf_final = get_lines_per_gdp_points(gdf_nearest_points, gdf_wb_points)
    # gdf_final.to_file(os.path.join(work_dir, "cutline_points.geojson"))
    # lines = gdf_final.groupby(["line"])["geometry"].apply(
    #     lambda x: LineString(x.tolist())
    # )
    # lines = gpd.GeoDataFrame(lines, geometry="geometry", crs=gdf_final.crs)
    # lines.reset_index(inplace=True)
    # lines["length"] = lines.geometry.length
    lines = compute_cutline_from_gdp_poly(gdf_wb_points, gdf_gdp)
    lines.to_file(os.path.join(work_dir, "cutline.geojson"))


prepare_inputs(
    "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/extract_marne_inpe_noholes.geojson",
    "/home/btardy/Documents/activites/WATER/GDP/extract/Marne-Giffaumont/dem_extract_Marne-Giffaumont.tif",
    "/home/btardy/Documents/activites/WATER/GDP/test_chain",
    100,
)
