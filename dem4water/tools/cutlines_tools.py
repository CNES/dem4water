#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""Provide all tools for find cut lines."""
import logging
from math import atan, cos, sin

import geopandas as gpd
import numpy as np
import rasterio
import shapely
from shapely import LineString, Point, box
from shapely.ops import split


def compute_distance(point1, point2, thresholding=None):
    """Compute distance between two points.

    Thresholding
    """
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


def find_perpendicular_bisector(point_a, point_b, bisector_length):
    """Draw the perpendicular bisector of a unique cutline.

    Parameters
    ----------
    point_a:
        side point the line
    point_b:
        opposite side point of point_a
    bisector_length:
        length in meter from the cutline to the point ending the bisector on a side
    """
    middle_point = Point((point_a.x + point_b.x) / 2, (point_a.y + point_b.y) / 2)
    theta = atan((point_b.y - point_a.y) / (point_b.x - point_a.x))
    point_c = Point(
        (middle_point.x - bisector_length * sin(theta)),
        (middle_point.y + bisector_length * cos(theta)),
    )
    point_d = Point(
        (middle_point.x + bisector_length * sin(theta)),
        (middle_point.y - bisector_length * cos(theta)),
    )

    line = LineString([point_c, point_d])
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


def is_point_valid(masked_dem, dem_transform, reservoir_shape, prev_point, alt):
    search = True
    while search:
        # Select a point
        alt_max = np.amax(masked_dem)
        if alt_max == -10000:
            logging.warning(
                "Searching point reach a no data value. "
                "The area is not valid. Trying another area."
            )
            return None, None, alt
        # The altitude seems correct convert to point
        l_indices = np.where(masked_dem == [alt_max])
        new_x, new_y = rasterio.transform.xy(
            dem_transform, list(l_indices[0]), list(l_indices[1])
        )
        # Ensure the selected point not cross water
        line = LineString([prev_point, (new_x, new_y)])
        mid_point = line.interpolate(line.length / 2)
        if not reservoir_shape.contains(mid_point):
            alt.append(alt_max)
            return new_x, new_y, alt
        logging.warning("Candidate cross the reservoir. Find other candidate.")
        # if invalid iterate over the masked dem until
        masked_dem[l_indices] = -10000


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
    small_erode_water,
):
    """"""
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
        # 0 is not a valid no data value as some dam can be in lower altitude
        mnt_array[~circle] = -10000
        # TODO: search valid point
        is_point_valid(mnt_array, mnt_transform, small_erode_water, prev_point, alt)
        alt.append(np.amax(mnt_array))
        l_indices = np.where(mnt_array == [np.amax(mnt_array)])
        new_x, new_y = rasterio.transform.xy(
            mnt_transform, list(l_indices[0]), list(l_indices[1])
        )
        points_save.append((new_x, new_y))
        stop = False
        if (new_x[0], new_y[0]) in coords_cutline:
            print("point already found")
            logging.info(
                "The current point was previously added to the line."
                f" Try {search_radius} on {direction} side"
                # f" Stop searching points on {direction} side."
            )
            # stop = True
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
                    # stop = True
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
                    # stop = True
                    # input(s)
                else:
                    cutline = LineString(coords_cutline + [(new_x[0], new_y[0])])
                    number_of_added_points += 1

        if alt[-1] > alt_max and number_of_added_points > 0:
            print(alt[-1], alt_max)
            print("Warn: Altitude max reached")
            logging.info(
                f"Altitude maximum reached on {direction} side. Stop searching points"
            )
            stop = True
        return cutline, stop, alt, points_save, number_of_added_points


def find_extent_to_line(
    gdf_cutline, mnt_raster, water_body, search_radius_max, alt_max, work_dir
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
    small_erode_water = water_body.geometry.buffer(-0.5).values[0]
    mask_polygon_start = list(
        gdf_split_w_wb[gdf_split_w_wb["search_loc"] == "left"].geometry
    )
    # Search points to complete line by the start of the base
    stop = False
    number_of_added_points = 0
    alt_seen = []
    points_save = []
    search_radius = 5
    coords_cut = list(cutline.coords)

    left_point = coords_cut[0]
    prev_left_point = coords_cut[1]
    while not stop:
        # coords_cut = list(cutline.coords)
        # left_point = coords_cut[0]
        # prev_left_point = coords_cut[1]
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
            small_erode_water,
        )
        search_radius += 5
        if search_radius > search_radius_max and not stop:
            logging.info("Max radius reached before alt_max on left side")
            stop = True
    logging.info(f"Number of points added on left side : {number_of_added_points}")
    mask_polygon_end = list(
        gdf_split_w_wb[gdf_split_w_wb["search_loc"] == "right"].geometry
    )
    # Search points to complete line from the end of the base
    stop = False
    number_of_added_points = 0
    alt_seen = []
    points_save = []
    search_radius = 5
    coords_cut = list(cutline.coords)

    right_point = coords_cut[-1]
    prev_right_point = coords_cut[-2]
    while not stop:
        # coords_cut = list(cutline.coords)
        # right_point = coords_cut[-1]
        # prev_right_point = coords_cut[-2]
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
            small_erode_water,
        )
        search_radius += 5

        if search_radius > search_radius_max:
            logging.info("Max radius reached before alt_max on right side")
            stop = True
    logging.info(f"Number of points added on right side : {number_of_added_points}")
    # gdf_cutline_final = gpd.GeoDataFrame(
    #     {"id_cutline": [1]}, geometry=[cutline], crs=gdf_cutline.crs
    # )
    # return gdf_cutline_final
    return list(cutline.coords)
