#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""This module define the study area from inputs."""
import logging
import os
import sys
import urllib.request
from time import perf_counter

import geopandas as gpd
import numpy as np
import rasterio as rio
from affine import Affine
from bmi_topography import Topography
from osgeo import gdal
from pyproj import Transformer
from rasterio.coords import BoundingBox

from dem4water.tools.save_raster import save_image


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


def download_cop30(minx, maxx, miny, maxy, output_path):
    """."""
    params = Topography.DEFAULT.copy()
    # TODO: check if needed to be in 4326
    params = {
        "dem_type": "COP30",
        "south": miny,
        "north": maxy,
        "west": minx,
        "east": maxx,
        "output_format": "GTiff",
        "cache_dir": output_path,
    }
    dem = Topography(**params)
    dem.fetch()
    # TODO: reproject in meters resolution if needed


def download_raster(layer, format, crs, bbox, width, height, outpath):
    link = (
        f"https://data.geopf.fr/wms-r?LAYERS={layer}&FORMAT={format}&SERVICE=WMS&"
        f"VERSION=1.3.0&REQUEST=GetMap&STYLES=&CRS={crs}&BBOX={bbox}"
        f"&WIDTH={width}&HEIGHT={height}"
    )
    print(f"Request:\n{link}\n")
    urllib.request.urlretrieve(link, outpath)


def translate(outfile, sourcefile, crs, bounds):
    srcDS = gdal.Open(sourcefile)
    gdal.Translate(outfile, srcDS, outputSRS=f"EPSG:{crs}", outputBounds=bounds)


def download_rgealti(minx, maxx, miny, maxy, epsg, out_file):
    """."""
    p1 = (minx, miny)
    p2 = (minx, maxy)
    p3 = (maxx, miny)
    width = int(np.ceil(compute_distance(p1, p2)))
    height = int(np.ceil(compute_distance(p1, p3)))
    bbox = f"{minx},{miny},{maxx},{maxy}"
    download_raster(
        "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES",
        "image/tiff",
        f"EPSG:{epsg}",
        bbox,
        width,
        height,
        out_file,
    )
    outtranslate = out_file.replace(".tif", "_translate.tif")

    translate(outtranslate, out_file, epsg, [minx, maxy, maxx, miny])


def extract_from_vrt(to_crop_file, minx, maxx, miny, maxy, t_epsg, out_file):
    """."""
    with rio.open(to_crop_file) as crop:
        transformer = Transformer.from_crs(t_epsg, crop.crs, always_xy=True)
        x_proj, y_proj = transformer.transform([minx, maxx], [miny, maxy])
        bounds_crop = BoundingBox(
            left=np.min(x_proj),
            bottom=np.min(y_proj),
            right=np.max(x_proj),
            top=np.max(y_proj),
        )

        # Crop the data
        window_crop = rio.windows.from_bounds(
            bounds_crop.left,
            bounds_crop.bottom,
            bounds_crop.right,
            bounds_crop.top,
            crop.transform,
        )
        crop_array = crop.read(window=window_crop)
        # Manage profile to write result in file
        profile = crop.profile
        dst_transform = Affine(
            crop.res[0], 0.0, bounds_crop.left, 0.0, -crop.res[1], bounds_crop.top
        )
        profile.update(
            {
                "transform": dst_transform,
                "height": crop_array.shape[1],
                "width": crop_array.shape[2],
                "driver": "GTiff",
            }
        )
        save_image(crop_array, profile, out_file)


def area_mapping(
    dam_database,
    dam_name,
    retrieve_mode,
    epsg,
    out_dir,
    t_epsg,
    to_crop_file=None,
    dam_name_col="DAM_NAME",
    buffer_roi=1000,
    debug=False,
):
    """."""
    t1_start = perf_counter()
    # Silence VRT related error (bad magic number)
    gdal.PushErrorHandler("CPLQuietErrorHandler")
    logging_format = (
        "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
    )
    if debug is True:
        logging.basicConfig(
            stream=sys.stdout, level=logging.DEBUG, format=logging_format
        )
    else:
        logging.basicConfig(
            stream=sys.stdout, level=logging.INFO, format=logging_format
        )
    logging.info(f"Starting area_mapping.py for {dam_name}")
    gdf_dam = gpd.read_file(dam_database)
    gdf_dam = gdf_dam.to_crs(epsg)
    gdf_dam = gdf_dam[gdf_dam[dam_name_col] == dam_name]
    # gdf_dam.to_file(
    #     "/home/btardy/Documents/activites/France2030/Dem_RGE/Alesani.geojson"
    # )

    gdf_dam.geometry = gdf_dam.geometry.buffer(buffer_roi)
    # gdf_dam.to_file(
    #     "/home/btardy/Documents/activites/France2030/Dem_RGE/Test_buffer_area.geojson"
    # )
    # # bbox = gdf_dam.total_bounds
    minx, miny, maxx, maxy = gdf_dam.total_bounds
    out_file = os.path.join(out_dir, f"dem_extract_{dam_name}.tif")
    if retrieve_mode == "RGE_ALTI":
        download_rgealti(minx, maxx, miny, maxy, epsg, out_file)
    elif retrieve_mode == "HPC":
        extract_from_vrt(to_crop_file, minx, maxx, miny, maxy, t_epsg, out_file)
    elif retrieve_mode == "":
        download_cop30(minx, maxx, miny, maxy, out_dir)
    else:
        raise ValueError(f"{retrieve_mode} is not a valid choice ")

    t1_stop = perf_counter()
    logging.info(f"Elapsed time: {t1_stop} s {t1_start} s")

    logging.info(f"Elapsed time during the whole program in s : {t1_stop-t1_start} s")


# area_mapping(
#     "/home/btardy/Documents/activites/WATER/Barrages-15m-France_INPE-V0_340barrages_20221002/392-retenues-pourLoiZSV_V4_sans_tampon.geojson",
#     "Alesani",
#     "RGE_ALTI",
#     2154,
#     "/home/btardy/Documents/activites/France2030/Test_area_mapping"
# )
