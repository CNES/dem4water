"""This module provides tools to manipulate the watermap."""

import logging
import os

import geopandas as gpd
import numpy as np
import rasterio as rio
#import otbApplication as otb
from dem4water.tools.save_raster import save_image

def create_water_mask(watermap, water_thres=0.15):
    """Find the water body."""
    binary_wmap = watermap.replace(".tif", "_binary.tif")
    raster_watermap=rio.open(watermap)
    raster_wmap = raster_watermap.read()
    raster_wmap_profile = raster_watermap.profile
    binary_watermap=np.where(raster_wmap > water_thres, 1,0)
    
    save_image(binary_watermap, raster_wmap_profile, binary_wmap)
    shp_wmap = watermap.replace(".tif", ".shp")

    cmd = f"gdal_polygonize.py {binary_wmap} {shp_wmap}"
    os.system(cmd)

    if not os.path.exists(shp_wmap):
        raise OSError(f"{shp_wmap} not exists")
    return shp_wmap


def get_largest_water_body(shp_wmap, watermap):
    """Find the largest area intersecting the point."""
    gdf = gpd.GeoDataFrame().from_file(shp_wmap)
    gdf["area"] = gdf.geometry.area
    # Get only water bodies
    gdf_filtered = gdf[gdf.DN == 1]
    gdf_out = gdf_filtered[gdf_filtered.area == max(gdf_filtered.area)]
    water_body = watermap.replace(".tif", "_water_body.shp")
    gdf_out.to_file(water_body)

    # Rasterize
    water_body_im = water_body.replace(".shp", ".tif")
    rasterize = otb.Registry.CreateApplication("Rasterization")
    rasterize.SetParameterString("in", water_body)
    rasterize.SetParameterString("im", watermap)
    rasterize.SetParameterString("out", water_body_im)
    rasterize.SetParameterString("mode", "binary")
    rasterize.SetParameterString("mode.binary.foreground", "1")
    rasterize.ExecuteAndWriteOutput()
    return water_body_im


def compute_area_from_water_body(daminfo, shp_wmap):
    """Compute intersection between insider and water bodies."""
    gdf_wmap = gpd.GeoDataFrame().from_file(shp_wmap)
    gdf_wmap["area"] = gdf_wmap.geometry.area
    gdf_filtered = gdf_wmap[gdf_wmap.DN == 1]

    gdf_dam = gpd.GeoDataFrame().from_file(daminfo)
    gdf_insider = gdf_dam[gdf_dam.name == "Insider"]
    # convert point to polygon
    gdf_insider.geometry = gdf_insider.geometry.buffer(1)
    result = gdf_filtered.sjoin(gdf_insider)
    area = np.array(result.area)
    if len(area) != 1:
        logging.info("There is no matching waterbody. Check your data")
        return 0
    return area[0]


def compute_area_from_database_geom(database, damname, shp_wmap):
    """From database geojson compute area of waterbodies."""
    gdf_db = gpd.GeoDataFrame().from_file(database)
    gdf_db = gdf_db[gdf_db.DAM_NAME == damname]

    gdf_wmap = gpd.GeoDataFrame().from_file(shp_wmap)
    gdf_wmap = gdf_wmap[gdf_wmap.DN == 1]
    gdf_db = gdf_db.to_crs(gdf_wmap.crs)
    # result = gdf_wmap.sjoin(gdf_db)
    result = gdf_wmap.overlay(gdf_db, how="intersection")
    result["area"] = result.geometry.area
    area = np.array(result.area)
    res = 0
    for val in area:
        res += val
    return res
