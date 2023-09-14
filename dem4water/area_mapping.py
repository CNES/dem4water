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
from osgeo import gdal, ogr, osr
import rasterio as rio
from dem4water.tools.utils import distance
from dem4water.tools.save_raster import save_image
from dem4water.tools.extract_roi import ExtractROIParam,  extract_roi
from dem4water.tools.superimpose import SuperimposeParam, superimpose
from bmi_topography import Topography
from dem4water.tools.utils import distance
import rasterio as rio



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
        "south": bbox[2] - 2*long_radius,
        "north": bbox[3] + 2*long_radius,
        "west": bbox[0] - 2*lat_radius,
        "east": bbox[1] + 2*lat_radius,
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
    gdal.Warp(watermap_reproject, src_ds, dstSRS=dst_crs, width=largeur, height=hauteur, resampleAlg=gdal.GRA_Cubic)
   
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

    extract_roi_parameters_extw=ExtractROIParam(
        mode="radius",
        mode_radius_r= float(radius),   
        mode_radius_unitr="phy", 
        mode_radius_unitc="phy", 
        mode_radius_cx=point.GetX(),
        mode_radius_cy= point.GetY(), 
        dtype='float')
    
  
    extw, profile_etw = extract_roi(rio.open(watermap_reproject),extract_roi_parameters_extw)
    save_image(extw, profile_etw, out_wmap)
    superimpose_app=SuperimposeParam(interpolator ="bco",dtype= 'float')
    app, profile_app=superimpose(rio.open(dem), rio.open(out_wmap), superimpose_app)
    save_image( app, profile_app, out_dem)

    # Search dam bottom
    extract_roi_parameter_extw_bt=ExtractROIParam(
        mode="radius",
        mode_radius_r= 500,   
        mode_radius_unitr="phy", 
        mode_radius_unitc="phy", 
        mode_radius_cx=point.GetX(),
        mode_radius_cy= point.GetY(),
        dtype='float')
    
    extw_bt, profile_extw_bt = extract_roi(rio.open(out_wmap), extract_roi_parameter_extw_bt)
   
    superimpose_extd_bt=SuperimposeParam(interpolator ="bco",dtype= 'float')
 
    extd_bt,profile_extd_bt =superimpose(rio.open(out_dem), extw_bt, superimpose_extd_bt, profile_extw_bt)
    extd_bt=extd_bt[:, :extw_bt.shape[1], :extw_bt.shape[2]]
    np_surf=np.where(extw_bt > 0.50, extd_bt, str(calt))
    np_surf= np_surf.astype('float')
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
