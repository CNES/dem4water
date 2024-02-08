#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""This module provides tools to manipulate the watermap."""

import logging
import os

import geopandas as gpd
import numpy as np
import rasterio as rio

from dem4water.tools.save_raster import save_image



def create_water_mask(watermap, water_thres=0.15):
    """Find the water body."""
    binary_water_map = watermap.replace(".tif", "_binary.tif")
    raster_watermap = rio.open(watermap)
    raster_water_map = raster_watermap.read()
    raster_w_map_profile = raster_watermap.profile
    binary_watermap = np.where(raster_water_map > water_thres, 1, 0)

    save_image(binary_watermap, raster_w_map_profile, binary_water_map)
    shp_w_map = watermap.replace(".tif", ".shp")

    cmd = f"gdal_polygonize.py {binary_water_map} {shp_w_map}"
    os.system(cmd)

    if not os.path.exists(shp_w_map):
        raise OSError(f"{shp_w_map} not exists")
    return shp_w_map


def compute_area_from_water_body(dam_info, shp_w_map):
    """Compute intersection between insider and water bodies."""
    gdf_w_map = gpd.GeoDataFrame().from_file(shp_w_map)
    gdf_w_map["area"] = gdf_w_map.geometry.area
    gdf_filtered = gdf_w_map[gdf_w_map.DN == 1]

    gdf_dam = gpd.GeoDataFrame().from_file(dam_info)

    gdf_insider = gdf_dam[gdf_dam.name == "Insider"]
    # convert point to polygon
    gdf_insider.geometry = gdf_insider.geometry.buffer(1)
    result = gdf_filtered.sjoin(gdf_insider)
    area = np.array(result.area)
    if len(area) != 1:
        logging.info("There is no matching water body. Check your data")
        return 0
    return area[0]


def compute_area_from_database_geom(database, damname, shp_wmap):
    """From database geojson compute area of water bodies."""
    gdf_db = gpd.GeoDataFrame().from_file(database)
    gdf_db = gdf_db[gdf_db.DAM_NAME == damname]

    gdf_water_map = gpd.GeoDataFrame().from_file(shp_wmap)
    gdf_water_map = gdf_water_map[gdf_water_map.DN == 1]
    gdf_db = gdf_db.to_crs(gdf_water_map.crs)
    result = gdf_water_map.overlay(gdf_db, how="intersection")
    result["area"] = result.geometry.area
    area = np.array(result.area)
    res = 0
    for val in area:
        res += val
    return res
