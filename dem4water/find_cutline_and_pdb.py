#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""This module search the PDB and the cutline."""


import glob

# import json
import logging
import math
import os
import sys
from itertools import groupby
from math import atan, cos, sin

import geopandas as gpd
import numpy as np
import pandas as pd

# import pyproj
import rasterio
import shapely
from rasterio import Affine

# from rasterio import transform as TR
from rasterio.mask import mask
from shapely.geometry import LineString, Point, box
from shapely.ops import split

from dem4water.tools.compute_grandient_dot_product import compute_gradient_product
from dem4water.tools.polygonize_raster import polygonize
from dem4water.tools.rasterize_vectors import RasterizarionParams, rasterize
from dem4water.tools.remove_holes_in_shapes import close_holes

# from typing import Optional


# ########################################
# Database utils
# ########################################
def remove_mutlipolygon(in_vector, epsg, buffer_size=10):
    """Apply a buffer in both sign to close small gap between polygons."""
    gdf = gpd.GeoDataFrame().from_file(in_vector)
    # print(epsg)
    gdf = gdf.to_crs(epsg)
    gdf["geometry"] = gdf.geometry.buffer(buffer_size)
    gdf["geometry"] = gdf.geometry.buffer(-buffer_size)
    gdf = gdf.explode(ignore_index=True)
    return gdf


def preprocess_water_body(in_vector, epsg, simplify=False, buff_simplify=200):
    """Simplify the water body to GDP."""
    # Fuse multipolygon
    gdf_wb = remove_mutlipolygon(in_vector, epsg, 10)
    # Remove holes
    gdf_wb.geometry = gdf_wb.geometry.apply(lambda p: close_holes(p))
    if simplify:
        gdf_wb.geometry = gdf_wb.geometry.buffer(buff_simplify)
        gdf_wb.geometry = gdf_wb.geometry.buffer(-buff_simplify)
    if len(gdf_wb.index) > 1:
        logging.error("WARNING: multiple water body detected, process only the larger")
        gdf_wb = gdf_wb[gdf_wb.geometry.area == max(gdf_wb.geometry.area)]
    return gdf_wb


# ##############################################
# GDP utils
# ##############################################
def clear_polygonize(in_vector, epsg):
    """Remove all 0 polygons from a vector file."""
    gdf = gpd.GeoDataFrame().from_file(in_vector)
    gdf = gdf.to_crs(epsg)
    if gdf.empty:
        return None
    gdf = gdf.loc[gdf["DN"] == 1]
    # input(gdf.crs)
    return gdf


def merge_close_gdp(gdf_gdp, max_dist):
    """Merge gdp entity if they are close enough."""
    # find intersection according the max distance
    gdf_gdp.geometry = gdf_gdp.geometry.buffer(max_dist)
    # then fuse
    gdf_gdp = gdf_gdp.dissolve(by="DN")
    # remove the buffer
    gdf_gdp.geometry = gdf_gdp.geometry.buffer(-max_dist)
    # then separe into unique object
    gdf_gdp = gdf_gdp.explode(ignore_index=True)
    # dissolve remove the DN column
    gdf_gdp["gdp_unique_id"] = range(len(gdf_gdp.index))
    return gdf_gdp


def filter_by_convex_hull(gdf, id_col):
    """Filter insider polygon from decreasing area."""
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
    gdf_filtered = gdf[gdf[id_col].isin(view_id)]
    return gdf_filtered


# ####################################################
# Format conversion
# ####################################################
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
    list_of_points = [
        Point(x[0], x[1]) for x in gdf_w.geometry.values[0].exterior.coords
    ]

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
    # gdf_w.geometry = [LineString(new_list_of_points)]
    # gdf_w = gdf_w.to_crs(gdf_wb.crs)
    gdf_f = gpd.GeoDataFrame(
        {"id_point": range(len(new_list_of_points))},
        geometry=new_list_of_points,
        crs=2154,
    )
    gdf_f = gdf_f.to_crs(gdf_wb.crs)
    return gdf_f


# ########################################################
# Points extraction
# ########################################################
def convert_coord_to_pixel(col, row, transform):
    """Convert coordinates according a transform."""
    return (col, row) * transform


def extract_points_from_raster(raster, poly_dataframe):
    """Extract all points inside a polygon."""
    geoms = poly_dataframe.geometry.values
    with rasterio.open(raster) as src:
        out_image, out_transform = mask(src, geoms, crop=True)
        # reference the pixel centre
        transf = out_transform * Affine.translation(0.5, 0.5)
        transformer = rasterio.transform.AffineTransformer(transf)

        no_data = src.nodata
        data = out_image[0, :, :]
        row, col = np.where(data != no_data)
        values = np.extract(data != no_data, data)
        df_p = pd.DataFrame({"col": col, "row": row, "val": values})
        df_p["x"] = df_p.apply(
            lambda row: transformer.xy(row.row, row.col)[0]
            # convert_coord_to_pixel(row.row, row.col, transf)[0]
            ,
            axis=1,
        )
        df_p["y"] = df_p.apply(
            lambda row: transformer.xy(row.row, row.col)[1]
            # convert_coord_to_pixel(row.row, row.col, transf)[1]
            ,
            axis=1,
        )
        gdf = gpd.GeoDataFrame(
            df_p,
            geometry=df_p.apply(lambda row: Point(row["x"], row["y"]), axis=1),
            crs=poly_dataframe.crs,
        )
        return gdf


# #######################################################
# PDB
# ######################################################


def find_pdb(poly_dataframe, raster):
    """Find pdb point."""
    gdf_pdb = extract_points_from_raster(raster, poly_dataframe)
    min_dem = gdf_pdb["val"].idxmin()
    value_min = gdf_pdb.loc[min_dem]
    pdb = Point(value_min.x, value_min.y)
    return pdb, gdf_pdb["val"].min()


def find_insider(wb_poly, poly_gdp):
    """Find insider point."""
    poly = poly_gdp.copy()
    poly.geometry = poly.geometry.buffer(30)
    inter = gpd.sjoin(poly, wb_poly)
    insider = inter.geometry.representative_point()
    return insider


# #######################################################
# Find points from polygon to draw base lines
# #######################################################
def manage_gdp_and_origin(gdf, max_indices_wb, ident, index_line):
    """."""
    id_points = list(gdf.id_point)
    group = (n - i for i, n in enumerate(id_points))
    dist = list(group)
    max_dist = dist[-1]
    dist = (0 if x > max_dist * 0.1 else 1 for x in dist)
    ident_point = [g for _, (*g,) in groupby(id_points, lambda _: next(dist))]
    if len(ident_point) == 2:
        # print("Two points", "ident", ident)
        id_point_reoder = ident_point[-1] + ident_point[0]
        x_coords = []
        y_coords = []
        for val_point in id_point_reoder:
            # print(val_point)
            x_coords.append(gdf[gdf.id_point == val_point].geometry.x)
            y_coords.append(gdf[gdf.id_point == val_point].geometry.y)
        df_coords = pd.DataFrame()

        df_coords["x_cut"] = x_coords
        df_coords["y_cut"] = y_coords
        df_coords["gdp_unique_id"] = index_line
        df_coords["gdf_sub_line"] = ident
        ident += 1
        gdf_temp = gpd.GeoDataFrame(
            df_coords,
            geometry=gpd.points_from_xy(df_coords["x_cut"], df_coords["y_cut"]),
            crs=gdf.crs,
        )
    else:
        print("segments founds ", len(ident_point))
        if ident_point[-1][-1] == max_indices_wb:
            print("Join to 0")
            init_segment = ident_point[-1] + ident_point[0]
            ident_point = [init_segment] + ident_point[1:-1]
        gdf_temp, ident = draw_multisegment(gdf, ident_point, ident, index_line)
    return gdf_temp, ident


def draw_multisegment(gdf, ident_point, ident, index_line):
    """."""
    list_df = []
    for seg in ident_point:
        x_coords = []
        y_coords = []
        for val_point in seg:
            x_coords.append(gdf[gdf.id_point == val_point].geometry.x)
            y_coords.append(gdf[gdf.id_point == val_point].geometry.y)
        df_coords = pd.DataFrame()
        df_coords["x_cut"] = x_coords
        df_coords["y_cut"] = y_coords
        df_coords["gdp_unique_id"] = index_line
        df_coords["gdf_sub_line"] = ident
        ident += 1
        gdf_temp = gpd.GeoDataFrame(
            df_coords,
            geometry=gpd.points_from_xy(df_coords["x_cut"], df_coords["y_cut"]),
            crs=gdf.crs,
        )
        list_df.append(gdf_temp)
    gdf_temp = gpd.GeoDataFrame(
        pd.concat(list_df, ignore_index=True), crs=list_df[0].crs
    )
    return gdf_temp, ident


def manage_no_origin(gdf, ident, index_line):
    """."""
    id_points = list(gdf.id_point)
    group = (n - i for i, n in enumerate(id_points))
    dist = list(group)

    max_dist = dist[-1]
    dist = (0 if x > max_dist * 0.1 else 1 for x in dist)
    ident_point = [g for _, (*g,) in groupby(id_points, lambda _: next(dist))]
    return draw_multisegment(gdf, ident_point, ident, index_line)


def find_points_for_base_lines(
    gdf_wb_points, gdf_gdp_poly, range_buff, ident, index_line
):
    """."""
    gdf_work = gdf_gdp_poly.copy()
    # print(gdf_wb_points.crs, gdf_work.crs)
    for buffer_size in range(0, range_buff, 10):
        gdf_work.geometry = gdf_gdp_poly.geometry.buffer(buffer_size)
        gdf = gpd.sjoin(gdf_wb_points, gdf_work, how="inner", predicate="within")
        if gdf.empty:
            logging.info(
                f"No intersection found between gdp {ident} and the water body"
                f" for buffer size {buffer_size}. Try to increase"
            )
            continue
        if gdf.shape[0] > 2:
            logging.info(
                f"Buffer size {buffer_size} provide a line with"
                f" {gdf.shape[0]} points."
            )
            break
    if gdf.empty:
        logging.info(
            f"ERROR: No intersection found between gdp {ident} and the water body"
            f" for buffer size {buffer_size}. Try to increase"
        )
    # print(gdf)
    min_indices_gdp = gdf.id_point.min()
    max_indices_gdp = gdf.id_point.max()
    max_indices_wb = gdf_wb_points.id_point.max()
    if min_indices_gdp < 100 and max_indices_gdp > max_indices_wb - 100:
        print("passe par l'origine")
        gdf_line = manage_gdp_and_origin(gdf, max_indices_wb, ident, index_line)
    else:
        print("Pas origine")
        gdf_line = manage_no_origin(gdf, ident, index_line)

    return gdf_line


# ##################################################
# Approach using angles
# ##################################################
def get_angle(a, b, c):
    """Provide the angle A between the points A,B,C."""
    ab = math.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]))
    ac = math.sqrt((a[0] - c[0]) * (a[0] - c[0]) + (a[1] - c[1]) * (a[1] - c[1]))
    cb = math.sqrt((c[0] - b[0]) * (c[0] - b[0]) + (c[1] - b[1]) * (c[1] - b[1]))
    ang = math.acos((ac * ac + ab * ab - cb * cb) / (2 * ac * ab))
    return np.degrees(ang)


def find_base_line_using_segments(
    gdf_wb, gdf_gdp, ident, index_line, buffer_size_max=50, step_buff=5, angle_thres=130
):
    """."""
    gdf_simpl = gdf_wb.copy()
    # gdf_simpl = gdf_simpl.to_crs(2154)
    # 1. Generate segments
    # Reduce the number of angles to try to find direction changes
    gdf_simpl.geometry = gdf_simpl.geometry.simplify(1)
    coords_points = gdf_simpl.geometry.values[0].exterior.coords
    simple_coords = list(coords_points)[:-1]
    list_a = simple_coords[:]
    list_b = simple_coords[1:] + [simple_coords[0]]
    list_c = [simple_coords[-1]] + simple_coords[:-1]
    list_rupt = []
    list_angle = []
    for a, b, c in zip(list_a, list_b, list_c):
        list_angle.append(get_angle(a, b, c))
        if get_angle(a, b, c) < angle_thres:
            list_rupt.append(a)
    # merge_last_segment = False
    #     merge_last_segment = True

    points_x = [x[0] for x in list_a]
    points_y = [x[1] for x in list_a]
    points = [Point(x, y) for x, y in zip(points_x, points_y)]
    gdf = gpd.GeoDataFrame(
        {
            "point_x": points_x,
            "point_y": points_y,
            "i": range(len(list_angle)),
            "angle": list_angle,
        },
        geometry=points,
        crs=gdf_wb.crs,
    )
    # gdf.to_file(
    #     "/home/btardy/Documents/activites/WATER/GDP/bertier_new_cut/Points_angle.geojson"
    # )
    # gdf_wb1 = gdf_wb.to_crs(2154)
    coords_points_ori = list(gdf_wb.geometry.values[0].exterior.coords)
    segments = []
    seg = []
    for point in coords_points_ori:
        seg.append(point)
        if point in list_rupt:
            segments.append(seg)
            seg = []
    ori_lines = []
    for seg in segments:
        if len(seg) > 2:
            ori_lines.append(LineString(seg))

    gdf = gpd.GeoDataFrame(
        {"segments": range(len(ori_lines))},
        geometry=ori_lines,
        crs=gdf_wb.crs,
    )
    # gdf.to_file(
    #     "/home/btardy/Documents/activites/WATER/GDP/bertier_new_cut/segments.geojson"
    # )
    for buff in range(0, buffer_size_max + step_buff, step_buff):
        gdf_gdp_buff = gdf_gdp.copy()
        gdf_gdp_buff.geometry = gdf_gdp_buff.geometry.buffer(buff)
        inter = gpd.sjoin(gdf, gdf_gdp_buff, how="inner")
        # handle the case of mutliple segment intersecting the same GDP
        if inter.empty:
            logging.info("Search base for cutline failed. No intersection found")
            logging.info(
                f"No intersection for buffer size {buff}. "
                f"Try higher (max {buffer_size_max})"
            )
        else:
            logging.info(f"Intersection found for buffer {buff}. Process")
            break
        if buffer_size_max == buff:
            return None, ident
    if not len(inter.index) == 1:
        logging.info("More than one segment intersect the same GDP. Try to fuse them")

    lines = list(inter.geometry.values)
    id_seg = list(inter.segments.values)
    # print(id_seg, len(ori_lines))
    coords_new_seg = []
    if id_seg[0] == 0 and id_seg[-1] == len(ori_lines) - 1:
        print("fuse origine")
        coords_new_seg = list(lines[-1].coords) + list(lines[0].coords)
        for seg in lines[1:-1]:
            coords_new_seg += seg.coords
    else:
        for seg in lines:
            coords_new_seg += seg.coords

    # line = LineString(coords_new_seg)
    # inter = gpd.GeoDataFrame(
    #     {"segments": range(len([line]))},
    #     geometry=[line],
    #     crs=gdf_wb.crs,
    # )
    df_coords = pd.DataFrame()
    df_coords["x_cut"] = [x[0] for x in coords_new_seg]
    df_coords["y_cut"] = [x[1] for x in coords_new_seg]
    df_coords["gdf_unique_id"] = index_line
    df_coords["gdf_sub_line"] = ident
    ident += 1
    inter = gpd.GeoDataFrame(
        df_coords,
        geometry=gpd.points_from_xy(df_coords["x_cut"], df_coords["y_cut"]),
        crs=gdf_wb.crs,
    )
    return inter, ident


# def fuse_segment(segments, water_body):
#     """."""


# #######################################################
# Find extent to base cutline
# #######################################################
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


def find_perpendicular_bisector(point_a, point_b, bisector_lenght):
    """Draw the perpendicular bisector of a unique cutline.

    Parameters
    ----------
    cutline_base:
        geojson or equivalent file containing the cutline
    bisector_lenght:
        length in meter from the cutline to the point ending the bisector
        on a side, i.e
    """
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
    # gdf.geometry = [line]
    # gdf.to_file(bisector_file)
    return line


def cut_area_according_perpendicular_bisector(ref_image, coords_cutline, crs):
    """."""
    with rasterio.open(ref_image) as image:
        bounds = image.bounds
        poly_image = box(*bounds)
        distance_max = (
            compute_distance([bounds.top, bounds.left], [bounds.bottom, bounds.right])
            + 10
        )

        point_a = coords_cutline[0]
        point_b = coords_cutline[-1]
        per_bisect = find_perpendicular_bisector(point_a, point_b, distance_max)
        search_area = split(poly_image, per_bisect)
        # search_area must have only two parts
        search_area = list(search_area.geoms)
        # print(search_area)
        first_poly = search_area[0]
        second_poly = search_area[1]
        if first_poly.contains(point_a):
            geom = [first_poly, second_poly]
        else:
            geom = [second_poly, first_poly]
        gdf = gpd.GeoDataFrame(
            {"search_loc": ["left", "right"]}, geometry=geom, crs=crs
        )
        return gdf


def search_next_point(mnt_raster, shapes, work_dir, direction, index):
    """."""
    with rasterio.open(mnt_raster) as dem:
        image, transform = rasterio.mask.mask(dem, shapes, crop=True)
        transformer = rasterio.transform.AffineTransformer(transform)
        alt_max = np.amax(image)
        coord_max = np.where(image == alt_max)
        new_coord = transformer.xy(coord_max[1], coord_max[2])
        new_point = Point(new_coord[0][0], new_coord[1][0])
        out_profile = dem.profile.copy()

        out_profile.update(
            {"width": image.shape[2], "height": image.shape[1], "transform": transform}
        )
        with rasterio.open(
            work_dir + f"/mask_{direction}_{index}.tif", "w", **out_profile
        ) as dst:
            dst.write(image)
        return new_point, alt_max


def add_point_to_list(value, list_values, direction):
    """Add point on start or end of list."""
    if direction == "left":
        # print("value", value)
        list_values = [value] + list_values
    else:
        list_values.append(value)
    return list_values


def follow_direction(
    cutline_coords,
    direction,
    mask_polygon,
    search_area_size,
    mnt_raster,
    work_dir,
    index,
):
    """."""
    starting_point = cutline_coords[0] if direction == "left" else cutline_coords[-1]
    print("starting_point", index, starting_point, mask_polygon.crs)
    gdf_point = gpd.GeoDataFrame(
        {"point": ["a"]}, geometry=[starting_point], crs=mask_polygon.crs
    )
    gdf_point.geometry = gdf_point.geometry.buffer(search_area_size)

    # Define the area where finding a point
    mask_polygon.geometry = mask_polygon.geometry.buffer(1)
    search_area = mask_polygon.overlay(gdf_point, how="intersection")
    search_area = search_area.explode(ignore_index=True)
    # print(search_area)
    if search_area.empty:
        logging.info(f"Searching point {direction} is stopped. Empty search area")
        return cutline_coords, None, mask_polygon, True
    search_area.to_file(work_dir + f"/pot_search_{direction}_{index}.geojson")

    # crop the mnt
    shapes = [p for p in search_area.geometry.values]  # if p.contains(starting_point)]
    search_area = gpd.GeoDataFrame(
        {"i": range(len(shapes))}, geometry=shapes, crs=mask_polygon.crs
    )
    if search_area.empty:
        logging.info(
            f"Searching point {direction} is stopped. Empty search area after filtering"
        )
        return cutline_coords, None, mask_polygon, True
    search_area.to_file(work_dir + f"/search_{direction}_{index}.geojson")

    # find the max in the area
    new_point, alt_max = search_next_point(
        mnt_raster, shapes, work_dir, direction, index
    )
    # print(alt_max)
    # add the point to cutline
    if alt_max < -100:
        logging.info("Invalid altitude found. Stop search")
        return cutline_coords, None, mask_polygon, True
    if new_point in cutline_coords:
        logging.info("Cutline extent: point already found stop search")
        return cutline_coords, None, mask_polygon, True
    cutline_coords = add_point_to_list(new_point, cutline_coords, direction)
    # remove the previous area to the search area ?
    # search_area.geometry = search_area.geometry.buffer(search_area_size)
    mask_polygon = mask_polygon.overlay(
        search_area,
        how="difference",
    )
    return cutline_coords, alt_max, mask_polygon, False


def find_extent_to_line_with_vector(
    mnt_raster, gdf_cutline, water_body, radius_search, maximum_alt, work_dir
):
    """.

    Parameters
    ----------
    mnt_raster: the mnt area around the reservoir
    cutline: a geodataframe containing one line
    """
    # 1. Find the perpendicular bisector according the line
    # gdf_cutline = gpd.read_file(cutline)
    # print(gdf_cutline)
    coords_cutline = list(gdf_cutline.geometry.values)
    base = gpd.GeoDataFrame(
        {"i": [1]}, geometry=[LineString(coords_cutline)], crs=gdf_cutline.crs
    )
    base.to_file(work_dir + "/cutline_base.geojson")
    gdf_split = cut_area_according_perpendicular_bisector(
        mnt_raster, coords_cutline, gdf_cutline.crs
    )
    gdf_split.to_file(work_dir + "/split_area.geojson")
    # 2. Generate a left and right masked mnt
    # Remove the water body to the areas of searching points
    gdf_split_w_wb = gdf_split.overlay(water_body, how="difference")
    # TODO: ensure that there is always two polygons after removing waterplan

    # 3. For each direction search point until a stop criterion is reached
    mask_polygon = gdf_split_w_wb.loc[[0]]
    stop = False
    index = 0
    while not stop:
        coords_cutline, alt, mask_polygon, stop = follow_direction(
            coords_cutline,
            "left",
            mask_polygon,
            radius_search,
            mnt_raster,
            work_dir,
            index,
        )
        # print(coords_cutline)
        mask_polygon.to_file(work_dir + f"/left_mask_{index}.geojson")
        index += 1
        if alt is not None and alt > maximum_alt:
            logging.info(f"Cutline extent: max found on left side at {alt}")
            break
    mask_polygon = gdf_split_w_wb.loc[[1]]
    stop = False
    index = 0
    while not stop:
        coords_cutline, alt, mask_polygon, stop = follow_direction(
            coords_cutline,
            "right",
            mask_polygon,
            radius_search,
            mnt_raster,
            work_dir,
            index,
        )

        mask_polygon.to_file(work_dir + f"/right_mask_{index}.geojson")
        index += 1
        if alt is not None and alt > maximum_alt:
            logging.info(f"Cutline extent: max found on right side at {alt}")
            break
    return coords_cutline


# #######################################################
# Extend but with raster approach
# #######################################################
def search_point(
    init_point,
    prev_point,
    cutline,
    search_radius,
    mnt_raster,
    shapes,
    alt,
    alt_max,
    direction,
    points_save,
    number_of_added_points,
):
    # new_point, alt_max = search_next_point(mnt_raster, shapes, work_dir, direction, index)
    coords_cutline = list(cutline.coords)
    with rasterio.open(mnt_raster) as dem:
        mnt_array, mnt_transform = rasterio.mask.mask(dem, shapes, crop=True)
        mnt_array = mnt_array[0, :, :]
        x_grid, y_grid = np.meshgrid(
            np.arange(mnt_array.shape[0]), np.arange(mnt_array.shape[1]), indexing="ij"
        )
        center_x, center_y = rasterio.transform.rowcol(
            mnt_transform, init_point[0], init_point[1]
        )
        prev_point_x, prev_point_y = rasterio.transform.rowcol(
            mnt_transform, prev_point[0], prev_point[1]
        )
        mask_radius = np.ceil(
            np.sqrt(
                (prev_point_x - center_x) * (prev_point_x - center_x)
                + (prev_point_y - center_y) * (prev_point_y - center_y)
            )
        )
        disc_search = (
            (x_grid - center_x) ** 2 + (y_grid - center_y) ** 2
        ) <= search_radius**2
        disc_mask = (
            (x_grid - prev_point_x) ** 2 + (y_grid - prev_point_y) ** 2
        ) >= mask_radius**2
        circle = np.logical_and(disc_search, disc_mask)
        mnt_array[~circle] = 0
        alt.append(np.amax(mnt_array))
        l_indices = np.where(mnt_array == [np.amax(mnt_array)])
        new_x, new_y = rasterio.transform.xy(mnt_transform, l_indices[0], l_indices[1])
        points_save.append((new_x, new_y))
        stop = False
        if alt[-1] > alt_max and number_of_added_points > 0:
            print(alt[-1], alt_max)
            print("Warn: Altitude max reached")
            logging.info(
                f"Altitude maximum reached on {direction} side. Stop searching points"
            )
            stop = True
        if (new_x[0], new_y[0]) in coords_cutline:
            print("point already found")
            logging.info(
                "The current point was previously added to the line."
                f" Stop searching points on {direction} side."
            )
            stop = True
        else:
            if len(alt) > 2:
                if alt[-1] <= alt[-2]:
                    # stop = True
                    print("Warn: altitude decrease")
                    logging.info(f"Altitude decrease on {direction} side.")
        if not stop:
            if direction == "left":
                left_line = LineString([(new_x[0], new_y[0]), coords_cutline[0]])
                if left_line.length > 1000:
                    print("Left Point too far")
                    stop = True
                    return cutline, stop, alt, points_save, number_of_added_points
                # inter = shapely.intersection(left_line, water_body.geometry.values[0])
                # print("inter", inter)
                self_inter = shapely.intersection(
                    left_line,
                    cutline,
                )
                # If another point than the origin is found as intersection
                if isinstance(self_inter, shapely.geometry.MultiPoint):
                    print("Left Self intersection detected")
                    stop = True
                    # input(s)
                else:
                    cutline = LineString([(new_x[0], new_y[0])] + coords_cutline)
                    number_of_added_points += 1
            else:
                right_line = LineString([(new_x[0], new_y[0]), coords_cutline[-1]])

                if right_line.length > 1000:
                    print("Right Point too far")
                    stop = True
                    return cutline, stop, alt, points_save, number_of_added_points
                # inter = shapely.intersection(right_line, water_body.geometry.values[0])
                # print("inter", inter)
                self_inter = shapely.intersection(
                    right_line,
                    cutline,
                )
                if isinstance(self_inter, shapely.geometry.MultiPoint):
                    print("Right: Self intersection detected")
                    stop = True
                    # input(s)
                else:
                    cutline = LineString(coords_cutline + [(new_x[0], new_y[0])])
                    number_of_added_points += 1
        return cutline, stop, alt, points_save, number_of_added_points


def find_extent_to_line(
    gdf_cutline, mnt_raster, water_body, search_radius, alt_max, work_dir
):
    """."""
    # 1. Find the perpendicular bisector according the line
    # gdf_cutline = gpd.read_file(cutline)
    # print(gdf_cutline)
    coords_cutline = list(gdf_cutline.geometry.values)
    cutline = LineString(coords_cutline)
    base = gpd.GeoDataFrame({"i": [1]}, geometry=[cutline], crs=gdf_cutline.crs)
    base.to_file(work_dir + "/cutline_base.geojson")
    gdf_split = cut_area_according_perpendicular_bisector(
        mnt_raster, coords_cutline, gdf_cutline.crs
    )
    gdf_split.to_file(work_dir + "/split_area.geojson")
    # 2. Generate a left and right masked mnt
    # Remove the water body to the areas of searching points
    gdf_split_w_wb = gdf_split.overlay(water_body, how="difference")
    mask_polygon_start = list(
        gdf_split_w_wb[gdf_split_w_wb["search_loc"] == "left"].geometry
    )
    # Search points to complete line by the start of the base
    stop = False
    number_of_added_points = 0
    alt_seen = []
    points_save = []
    while not stop:
        coords_cut = list(cutline.coords)
        left_point = coords_cut[0]
        prev_left_point = coords_cut[1]
        cutline, stop, alt_seen, points_save, number_of_added_points = search_point(
            left_point,
            prev_left_point,
            cutline,
            search_radius,
            mnt_raster,
            mask_polygon_start,
            alt_seen,
            alt_max,
            "left",
            points_save,
            number_of_added_points,
        )
    logging.info(f"Number of points added on left side : {number_of_added_points}")
    mask_polygon_end = list(
        gdf_split_w_wb[gdf_split_w_wb["search_loc"] == "right"].geometry
    )
    # Search points to complete line from the end of the base
    stop = False
    number_of_added_points = 0
    alt_seen = []
    points_save = []
    while not stop:
        coords_cut = list(cutline.coords)
        right_point = coords_cut[-1]
        prev_right_point = coords_cut[-2]
        cutline, stop, alt_seen, points_save, number_of_added_points = search_point(
            right_point,
            prev_right_point,
            cutline,
            search_radius,
            mnt_raster,
            mask_polygon_end,
            alt_seen,
            alt_max,
            "right",
            points_save,
            number_of_added_points,
        )
    logging.info(f"Number of points added on right side : {number_of_added_points}")
    # gdf_cutline_final = gpd.GeoDataFrame(
    #     {"id_cutline": [1]}, geometry=[cutline], crs=gdf_cutline.crs
    # )
    # return gdf_cutline_final
    return list(cutline.coords)


# #######################################################
# Main function
# #######################################################
def find_pdb_and_cutline(
    database_file,
    dem_raster,
    work_dir,
    gdp_buffer_size,
    radius_search_size,
    maximum_alt,
    debug=False,
):
    """Search the PDB and cutline from the database and the DEM.

    Parameters
    ----------
    database_file:
    dem_raster:
    """
    logging_format = (
        "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
    )
    if debug:
        logging.basicConfig(
            stream=sys.stdout, level=logging.DEBUG, format=logging_format
        )
    else:
        logging.basicConfig(
            stream=sys.stdout, level=logging.INFO, format=logging_format
        )
    logging.info("Starting search cutline and PDB using GDP")

    # 1. Manage the particularity of database:
    # - Fuse multipolygon
    # - Remove holes inside water body
    # - Smooth waterbody (optional)
    # - Rasterize the waterbody
    with rasterio.open(dem_raster) as dem:
        epsg = dem.crs.to_epsg()
        gdf_wb = preprocess_water_body(database_file, epsg)
        gdf_wb.to_file(os.path.join(work_dir, "bd_clean.geojson"))
        params_raster = RasterizarionParams(
            mode="binary",
            binary_foreground_value=1,
            background_value=0,
            column_field=None,
            dtype="uint8",
            nodata=0,
        )
        waterbody_bin = os.path.join(work_dir, "waterbody_bin.tif")
        rasterize(gdf_wb, dem, params_raster, waterbody_bin)
        # 2. Compute the gradient dot product using dem and cleaned water body
        # - The output is a raster then vectorize it over the whole area
        # - Remove the 0 polygons which means no GDP found
        gdp_raster = os.path.join(work_dir, "gdp_raster.tif")
        compute_gradient_product(waterbody_bin, dem_raster, gdp_raster)
        gdp_vector = os.path.join(work_dir, "gdp_vector.geojson")
        polygonize(gdp_raster, gdp_vector)
        gdf_gdp = clear_polygonize(gdp_vector, epsg)
        if gdf_gdp is None:
            logging.info("ERROR: when computing GDP, no slope found")
            logging.info("ERROR: Stopping the chain")
            return None
        # 3. Fuse close GDP
        gdf_gdp = merge_close_gdp(gdf_gdp, gdp_buffer_size)
        if gdf_gdp.empty:
            logging.info("ERROR: Empty dataframe after merging GDP.")
            return None
        gdf_gdp.to_file(os.path.join(work_dir, "gdp_fusion_nearest.geojson"))
        # 4. Filter by convex hull to remove all insider GDP
        gdf_gdp = filter_by_convex_hull(gdf_gdp, "gdp_unique_id")
        # 5. TODO: filter by area to remove small GDP
        # print(gdf_gdp)
        gdf_gdp = gdf_gdp[gdf_gdp.geometry.area > 200]
        gdf_gdp = gdf_gdp.reset_index(drop=True)
        # print(gdf_gdp)
        # 6 Convert the water body into a sett of point with regular sampling
        wb_contour_points = oversampling_polygon_boundary(gdf_wb, dem.res[0])
        list_ident = []
        list_line = []
        ident = 0
        list_pdb = []
        list_alt_pdb = []
        for index in gdf_gdp.index:
            row = gdf_gdp.loc[[index]]
            pdb, alt = find_pdb(row, dem_raster)
            list_pdb.append(pdb)
            if alt >= -100:
                list_alt_pdb.append(alt)
            else:
                list_alt_pdb.append(np.nan)
        index_min = np.nanargmin(list_alt_pdb)
        print(list_alt_pdb)
        print("index min ", index_min)
        for index in [index_min]:  # gdf_gdp.index:
            row = gdf_gdp.loc[[index]]
            # 7. For each GDP find a PDB
            # pdb = find_pdb(row, dem_raster)
            # 8. For each GDP found a insider point
            # insider = find_insider(gdf_wb, row)
            # 9. For each GDP draw baselines
            print(row)
            wb_contour_points.to_file(work_dir + "/point_contour.geojson")
            row.to_file(work_dir + f"/gdp_{index}.geojson")
            # gdf_line, ident = find_points_for_base_lines(
            #     wb_contour_points, row, gdp_buffer_size, ident, index
            # )
            print("bef", ident)
            gdf_line, ident = find_base_line_using_segments(
                gdf_wb, row, ident, index, angle_thres=150
            )
            print("aft", ident)
            if gdf_line is None:
                continue
            # gdf_line.to_file(work_dir + f"/base_cutline_{index}.geojson")
            # 10. For each baseline find extents
            # coords_new_line = find_extent_to_line(
            #     dem_raster, gdf_line, gdf_wb, radius_search_size, maximum_alt, work_dir
            # )
            coords_new_line = find_extent_to_line(
                gdf_line, dem_raster, gdf_wb, radius_search_size, maximum_alt, work_dir
            )
            list_line.append(LineString(coords_new_line))
            list_ident.append(ident)

        # 11. Merge all cutlines into one file
        if not list_line:
            return None
        gdf_final = gpd.GeoDataFrame(
            {"ident_line": list_ident}, geometry=list_line, crs=gdf_wb.crs
        )
        gdf_final.to_file(os.path.join(work_dir, "cutline.geojson"))
        return os.path.join(work_dir, "cutline.geojson")


# df = gpd.read_file(
#     "/home/btardy/Documents/activites/WATER/GDP/extract/Berthier/berthier_bd.geojson"
# )
# alt = df.DAM_LVL_M.values[0]
# find_pdb_and_cutline(
#     database_file="/home/btardy/Documents/activites/WATER/GDP/extract/Berthier/berthier_bd.geojson",
#     dem_raster="/home/btardy/Documents/activites/WATER/GDP/extract/Berthier/dem_extract_Berthier.tif",
#     work_dir="/home/btardy/Documents/activites/WATER/GDP/bertier_new_cut",
#     gdp_buffer_size=50,
#     radius_search_size=20,
#     maximum_alt=alt + 20,
#     debug=False,
# )
working_dir = "/work/CAMPUS/etudes/hydro_aval/dem4water/work_benjamin/cutlines_v11"
out_file = "/work/CAMPUS/etudes/hydro_aval/dem4water/work_benjamin/cutlines_v11"
db_full = (
    "/work/CAMPUS/etudes/hydro_aval/MTE_2022_Reservoirs/livraisons"
    "/dams_database/db_CS/db_340_dams/340-retenues-pourLoiZSV_V6_sans_tampon_corrections.geojson"
)
gdf_db = gpd.GeoDataFrame().from_file(db_full)
extract_folder = (
    "/work/CAMPUS/etudes/hydro_aval/dem4water/work_benjamin/340_MAE_first_03/extracts"
)

# print(gdf_db.columns)
# print(list(gdf_db.DAM_LVL_M))
# print(list(gdf_db.DEPTH_M))
for DAM, ALTI, HAUTEUR in zip(
    list(gdf_db.DAM_NAME), list(gdf_db.DAM_LVL_M), list(gdf_db.DEPTH_M)
):
    # if dam not in ["Laparan", "Marne Giffaumont", "Alesani", "Borfloc h"]:
    #     continue
    print(DAM, ALTI)
    GDF_T = gdf_db.loc[gdf_db.DAM_NAME == DAM]
    DAM = DAM.replace(" ", "-")
    WDIR = os.path.join(working_dir, DAM)
    if not os.path.exists(WDIR):
        os.mkdir(WDIR)

    if not os.path.exists(os.path.join(WDIR, "cutline.geojson")):
        DAM_DB = os.path.join(WDIR, f"bd_{DAM}.geojson")
        GDF_T.to_file(DAM_DB)
        EXTRACT = glob.glob(extract_folder + f"/{DAM}/dem*.tif")[0]
        # print(extract)
        # out_extract = os.path.join(wdir, "dem_reproj.tif")
        # cmd = f"gdalwarp {extract} {out_extract} -t_srs 'EPSG:2154'"
        # os.system(cmd)
        if ALTI is None:
            ALTI = 10000
        # res = prepare_inputs(dam_db, extract, wdir, 100, alti + 20, 30)  # buffer size
        RES = find_pdb_and_cutline(
            DAM_DB,
            EXTRACT,
            WDIR,
            gdp_buffer_size=100,
            radius_search_size=30,
            maximum_alt=ALTI + 20,
            debug=False,
        )
        if RES is not None:
            CUTLINE = os.path.join(WDIR, "cutline.geojson")
            OUT_CUTLINE = os.path.join(working_dir, f"{DAM}_cutline.geojson")
            os.system(f"cp {CUTLINE} {OUT_CUTLINE}")
