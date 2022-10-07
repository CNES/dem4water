#!/usr/bin/env python3
"""
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CS Group
:copyright: 2020 CS Group. All rights reserved.
:license: see LICENSE file
:created: 2020
"""


import argparse
import logging
import os
import sys

import numpy as np
import otbApplication as otb
from osgeo import gdal, ogr, osr

from utils import distance
from time import perf_counter

def main(arguments):
    """area_mapping.py
    Retrieve dam coordinate
    """
    t1_start = perf_counter()
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--infile", help="Input file")
    parser.add_argument("--id", help="Dam ID")
    parser.add_argument("--id_db", help="Dam id field in database")
    parser.add_argument("-w", "--watermap", help="Input water map file")
    parser.add_argument("-d", "--dem", help="Input DEM")
    parser.add_argument("-r", "--radius", help="Extract radius (m)", default=2000)
    parser.add_argument("-o", "--out", help="Output directory")
    parser.add_argument("--debug", action="store_false", help="Activate Debug Mode")
    args = parser.parse_args(arguments)

    # Silence VRT related error (bad magic number)
    gdal.PushErrorHandler("CPLQuietErrorHandler")

    logging_format = (
        "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
    )
    if args.debug is True:
        logging.basicConfig(
            stream=sys.stdout, level=logging.DEBUG, format=logging_format
        )
    else:
        logging.basicConfig(
            stream=sys.stdout, level=logging.INFO, format=logging_format
        )
    logging.info("Starting area_mapping.py")

    driver = ogr.GetDriverByName("GeoJSON")
    dataSource = driver.Open(args.infile, 0)
    layer = dataSource.GetLayer()
    clat = 0
    clon = 0
    dam_name = ""
    dam_path = ""
    dam_404 = True
    radius = args.radius

    for feature in layer:

        if str(int(feature.GetField(str(args.id_db)))) == str(args.id):
            # Compute radius
            if radius == "":
                geom = feature.GetGeometryRef()
                bbox = geom.GetEnvelope()
                radius = distance(bbox[2], bbox[0], bbox[3], bbox[1])
                logging.info("=> RADIUS : {}".format(radius))

            dam_404 = False
            logging.debug(feature.GetField("DAM_NAME"))
            dam_name = feature.GetField("DAM_NAME")
            dam_path = dam_name.replace(" ", "-")
            clat = float(feature.GetField("LAT_DD"))
            clon = float(feature.GetField("LONG_DD"))
            logging.debug(
                "Height of dam in meters: " + str(feature.GetField("DAM_LVL_M"))
            )
            break
    layer.ResetReading()

    if dam_404 is True:
        logging.error(
            "404 - Dam Not Found: " + str(args.id) + " is not present in " + args.infile
        )
    # logging.info("Currently processing: "+ dam_name +"(ID:"+str(args.id)+") [Lat: "+ str(clat) +", Lon: "+ str(clon) "]")
    calt = float(
        os.popen(
            "gdallocationinfo -valonly -wgs84 %s %s %s" % (args.dem, clon, clat)
        ).read()
    )
    logging.info(
        "Currently processing: "
        + dam_name
        + "(ID:"
        + str(args.id)
        + ") [Lat: "
        + str(clat)
        + ", Lon: "
        + str(clon)
        + ", Alt: "
        + str(calt)
        + "]"
    )

    src = osr.SpatialReference()

    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    src.ImportFromEPSG(4326)
    ds = gdal.Open(args.watermap, gdal.GA_ReadOnly)
    dst = osr.SpatialReference(wkt=ds.GetProjection())
    ct = osr.CoordinateTransformation(src, dst)
    point = ogr.Geometry(ogr.wkbPoint)
    # point.AddPoint(clat, clon)
    point.AddPoint(clon, clat)
    point.Transform(ct)
    logging.debug("Coordinates: " + str(point.GetX()) + " - " + str(point.GetY()))

    extw = otb.Registry.CreateApplication("ExtractROI")
    extw.SetParameterString("in", args.watermap)
    extw.SetParameterString("mode", "radius")
    extw.SetParameterString("mode.radius.unitr", "phy")
    extw.SetParameterFloat("mode.radius.r", float(radius))
    extw.SetParameterString("mode.radius.unitc", "phy")
    extw.SetParameterFloat("mode.radius.cx", point.GetX())
    extw.SetParameterFloat("mode.radius.cy", point.GetY())
    extw.SetParameterString(
        "out", os.path.join(args.out, "wmap_extract-" + dam_path + ".tif")
    )
    extw.ExecuteAndWriteOutput()

    #  extd = otb.Registry.CreateApplication("ExtractROI")
    #  extd.SetParameterString("in", args.dem)
    #  extd.SetParameterString("mode","fit")
    #  extd.SetParameterString("mode.fit.im", os.path.join(args.out, "wmap_extract-"+dam_path+".tif"))
    #  extd.SetParameterString("out", os.path.join(args.out, "dem_extract-"+dam_path+".tif"))
    #  extd.ExecuteAndWriteOutput()

    app = otb.Registry.CreateApplication("Superimpose")
    app.SetParameterString(
        "inr", os.path.join(args.out, "wmap_extract-" + dam_path + ".tif")
    )
    app.SetParameterString("inm", args.dem)
    app.SetParameterString(
        "out", os.path.join(args.out, "dem_extract-" + dam_path + ".tif")
    )
    app.ExecuteAndWriteOutput()

    # Search dam bottom
    extw_bt = otb.Registry.CreateApplication("ExtractROI")
    extw_bt.SetParameterString(
        "in", os.path.join(args.out, "wmap_extract-" + dam_path + ".tif")
    )
    extw_bt.SetParameterString("mode", "radius")
    extw_bt.SetParameterString("mode.radius.unitr", "phy")
    extw_bt.SetParameterFloat("mode.radius.r", 500)
    extw_bt.SetParameterString("mode.radius.unitc", "phy")
    extw_bt.SetParameterFloat("mode.radius.cx", point.GetX())
    extw_bt.SetParameterFloat("mode.radius.cy", point.GetY())
    extw_bt.Execute()

    extd_bt = otb.Registry.CreateApplication("Superimpose")
    extd_bt.SetParameterInputImage("inr", extw_bt.GetParameterOutputImage("out"))
    extd_bt.SetParameterString(
        "inm", os.path.join(args.out, "dem_extract-" + dam_path + ".tif")
    )
    extd_bt.Execute()

    bm = otb.Registry.CreateApplication("BandMath")
    bm.AddImageToParameterInputImageList("il", extw_bt.GetParameterOutputImage("out"))
    bm.AddImageToParameterInputImageList("il", extd_bt.GetParameterOutputImage("out"))
    bm.SetParameterString("exp", "( im1b1  > 0.50 ) ? im2b1 : " + str(calt))
    bm.Execute()

    np_surf = bm.GetImageAsNumpyArray("out")
    bt_alt = np.amin(np_surf)
    logging.info("Bottom Alt: " + str(bt_alt))
    t1_stop = perf_counter()
    logging.info("Elapsed time:", t1_stop, 's', t1_start, 's')
 
    logging.info("Elapsed time during the whole program in s :",
       t1_stop-t1_start, 's')
    # Profiling:
    if args.debug is True:
        for r in range(200, 1001, 50):
            extw_l = otb.Registry.CreateApplication("ExtractROI")
            extw_l.SetParameterString(
                "in", os.path.join(args.out, "wmap_extract-" + dam_path + ".tif")
            )
            extw_l.SetParameterString("mode", "radius")
            extw_l.SetParameterString("mode.radius.unitr", "phy")
            extw_l.SetParameterFloat("mode.radius.r", r)
            extw_l.SetParameterString("mode.radius.unitc", "phy")
            extw_l.SetParameterFloat("mode.radius.cx", point.GetX())
            extw_l.SetParameterFloat("mode.radius.cy", point.GetY())
            extw_l.Execute()

            extd_l = otb.Registry.CreateApplication("Superimpose")
            extd_l.SetParameterInputImage("inr", extw_l.GetParameterOutputImage("out"))
            extd_l.SetParameterString(
                "inm", os.path.join(args.out, "dem_extract-" + dam_path + ".tif")
            )
            extd_l.Execute()

            np_extdl = extd_l.GetImageAsNumpyArray("out")
            extdl_alt = np.amin(np_extdl)

            bml = otb.Registry.CreateApplication("BandMath")
            bml.AddImageToParameterInputImageList(
                "il", extw_l.GetParameterOutputImage("out")
            )
            bml.AddImageToParameterInputImageList(
                "il", extd_l.GetParameterOutputImage("out")
            )
            bml.SetParameterString("exp", "( im1b1  > 0.50 ) ? im2b1 : " + str(calt))
            bml.Execute()

            np_bml = bml.GetImageAsNumpyArray("out")
            bml_alt = np.amin(np_bml)
            logging.info(
                "@radius= "
                + str(r)
                + "m: Local min (dem)= "
                + str(extdl_alt)
                + "m / Local min (dem+w>.5)= "
                + str(bml_alt)
                + "m"
            )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
