#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""This module define the study area from inputs."""


import argparse
import logging
import os
import sys
from time import perf_counter
import glob
import urllib.request
import shutil
import numpy as np
import otbApplication as otb
from osgeo import gdal, ogr, osr
from bmi_topography import Topography
from dem4water.tools.utils import distance


def area_mapping(
    infile,
    dam_id,
    id_db,
    out_dem,
    out_wmap,
    output_download_path,
    radius=None,
    debug=False,
):
    """Extract dem and watermap according the in-situ information provided in the DB.

    Retrieve dam coordinate
    """
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
    logging.info("Starting area_mapping.py")

    driver = ogr.GetDriverByName("GeoJSON")
    data_source = driver.Open(infile, 0)
    layer = data_source.GetLayer()
    clat = 0
    clon = 0
    dam_name = ""
    # dam_path = ""
    dam_404 = True

    for feature in layer:
        if str(int(feature.GetField(str(id_db)))) == str(dam_id):
            # Compute radius
            if radius is None:
                geom = feature.GetGeometryRef()
                bbox = geom.GetEnvelope()
                radius = distance(bbox[2], bbox[0], bbox[3], bbox[1])
                logging.info(f"=> RADIUS : {radius}")

            dam_404 = False
            logging.debug(feature.GetField("DAM_NAME"))
            dam_name = feature.GetField("DAM_NAME")
            # dam_path = dam_name.replace(" ", "-")
            clat = float(feature.GetField("LAT_DD"))
            clon = float(feature.GetField("LONG_DD"))
            logging.debug(f"Height of dam in meters: {feature.GetField('DAM_LVL_M')}")
            break
    layer.ResetReading()

    # Download DEM
    long_radius = abs(bbox[2] - bbox[3])
    lat_radius = abs(bbox[0] - bbox[1])
    params = Topography.DEFAULT.copy()
    params = {
        "dem_type": "COP30",
        "south": bbox[2] - long_radius,
        "north": bbox[3] + long_radius,
        "west": bbox[0] - lat_radius,
        "east": bbox[1] + lat_radius,
        "output_format": "GTiff",
        "cache_dir": output_download_path,
    }
    boulder = Topography(**params)
    boulder.fetch()
    dem = glob.glob(os.path.join(output_download_path, "COP30*"))[0]

    if dam_404 is True:
        logging.error("404 - Dam Not Found: {dam_id} is not present in {infile}")
    calt = float(
        os.popen(f"gdallocationinfo -valonly -wgs84 {dem} {clon} {clat}").read()
    )
    logging.info(
        f"Currently processing: {dam_name} (ID: {dam_id}) [Lat: {clat}, Lon: {clon}, Alt: {calt}]"
    )

    # Download occurence
    DATASET_NAME = "occurrence"
    lg = abs(bbox[0])
    lt = abs(bbox[2])

    long = int(lg // 10 * 10) + 10
    lat = int(lt // 10 * 10) + 10

    if bbox[0] < 0 and bbox[1] < 0:
        long = str(long) + "W"
    elif bbox[0] > 0 and bbox[1] > 0:
        long = str(long) + "E"
    if bbox[2] < 0 and bbox[3] < 0:
        lat = str(lat) + "S"
    if bbox[2] > 0 and bbox[3] > 0:
        lat = str(lat) + "N"

    filename = DATASET_NAME + "_" + str(long) + "_" + str(lat) + "v1_4_2021.tif"
    url = os.path.join(
        "https://storage.googleapis.com/global-surface-water/downloads2021",
        DATASET_NAME,
        filename,
    )
    code = urllib.request.urlopen(url).getcode()
    if code != 404:
        print("Downloading " + url + ")")
        urllib.request.urlretrieve(url, os.path.join(output_download_path, filename))
    else:
        print(url + " not found")

    watermap = glob.glob(os.path.join(output_download_path, "occurrence*"))[0]

    watermap_reproject = watermap.replace(".tif", "_reproject.tif")
    dst_crs = "EPSG:32630"
    src_ds = gdal.Open(watermap)

    largeur = src_ds.RasterXSize
    hauteur = src_ds.RasterYSize
    gdal.Warp(watermap_reproject, src_ds, dstSRS=dst_crs, width=largeur, height=hauteur)

    src = osr.SpatialReference()

    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    src.ImportFromEPSG(4326)
    dataset = gdal.Open(watermap_reproject, gdal.GA_ReadOnly)
    dst = osr.SpatialReference(wkt=dataset.GetProjection())
    coord_trans = osr.CoordinateTransformation(src, dst)
    point = ogr.Geometry(ogr.wkbPoint)
    # point.AddPoint(clat, clon)
    point.AddPoint(clon, clat)
    point.Transform(coord_trans)
    logging.debug(f"Coordinates: {point.GetX()} - {point.GetY()}")

    extw = otb.Registry.CreateApplication("ExtractROI")
    extw.SetParameterString("in", watermap_reproject)
    extw.SetParameterString("mode", "radius")
    extw.SetParameterString("mode.radius.unitr", "phy")
    extw.SetParameterFloat("mode.radius.r", float(radius))
    extw.SetParameterString("mode.radius.unitc", "phy")
    extw.SetParameterFloat("mode.radius.cx", point.GetX())
    extw.SetParameterFloat("mode.radius.cy", point.GetY())
    extw.SetParameterString("out", out_wmap)
    extw.ExecuteAndWriteOutput()

    app = otb.Registry.CreateApplication("Superimpose")
    app.SetParameterString("inr", out_wmap)
    app.SetParameterString("inm", dem)
    app.SetParameterString("out", out_dem)
    app.ExecuteAndWriteOutput()

    # Search dam bottom
    extw_bt = otb.Registry.CreateApplication("ExtractROI")
    extw_bt.SetParameterString("in", out_wmap)
    extw_bt.SetParameterString("mode", "radius")
    extw_bt.SetParameterString("mode.radius.unitr", "phy")
    extw_bt.SetParameterFloat("mode.radius.r", 500)
    extw_bt.SetParameterString("mode.radius.unitc", "phy")
    extw_bt.SetParameterFloat("mode.radius.cx", point.GetX())
    extw_bt.SetParameterFloat("mode.radius.cy", point.GetY())
    extw_bt.Execute()

    extd_bt = otb.Registry.CreateApplication("Superimpose")
    extd_bt.SetParameterInputImage("inr", extw_bt.GetParameterOutputImage("out"))
    extd_bt.SetParameterString("inm", out_dem)
    extd_bt.Execute()

    bandmath = otb.Registry.CreateApplication("BandMath")
    bandmath.AddImageToParameterInputImageList(
        "il", extw_bt.GetParameterOutputImage("out")
    )
    bandmath.AddImageToParameterInputImageList(
        "il", extd_bt.GetParameterOutputImage("out")
    )
    bandmath.SetParameterString("exp", "( im1b1  > 0.50 ) ? im2b1 : " + str(calt))
    bandmath.Execute()

    np_surf = bandmath.GetImageAsNumpyArray("out")
    bt_alt = np.amin(np_surf)
    """
    if os.path.isdir(output_download_path):
        shutil.rmtree(output_download_path)
    """
    logging.info(f"Bottom Alt: {bt_alt}")
    t1_stop = perf_counter()
    logging.info(f"Elapsed time: {t1_stop} s {t1_start} s")

    logging.info(f"Elapsed time during the whole program in s : {t1_stop-t1_start} s")


def area_mapping_args():
    """Define area mapping parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--infile", help="Input file")
    parser.add_argument("--id", help="Dam ID")
    parser.add_argument("--id_db", help="Dam id field in database")
    parser.add_argument("-w", "--watermap", help="Input water map file")
    parser.add_argument("-d", "--dem", help="Input DEM")
    parser.add_argument("--out_dem", help="Extracted dem")
    parser.add_argument("--out_wmap", help="Extracted wmap")
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
    return parser


def main():
    """Cli function to launch area mapping."""
    parser = area_mapping_args()
    args = parser.parse_args()
    area_mapping(
        args.infile,
        args.id,
        args.id_db,
        args.radius,
        args.out_dem,
        args.out_wmap,
        args.output_download_path,
        debug=args.debug,
    )


if __name__ == "__main__":
    sys.exit(main())
