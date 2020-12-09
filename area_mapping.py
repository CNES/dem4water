#!/usr/bin/env python3
'''
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CS Group
:copyright: 2020 CS Group. All rights reserved.
:license: see LICENSE file
:created: 2020
'''

""" area_mapping.py
Retrieve dam coordinate
"""

import os
import sys
import argparse
from osgeo import ogr
from osgeo import osr
from osgeo import gdal
import numpy as np
import otbApplication as otb


def main(arguments):
    '''Entrypoint'''

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i',
                        '--infile',
                        help="Input file")
    parser.add_argument('-n',
                        '--name',
                        help="Dam Name")
    parser.add_argument('-w',
                        '--watermap',
                        help="Input water map file")
    parser.add_argument('-d',
                        '--dem',
                        help="Input DEM")
    parser.add_argument('-r',
                        '--radius',
                        help="Extract radius (m)",
                        default=2000)
    parser.add_argument('-o',
                        '--out',
                        help="Output directory")

    args = parser.parse_args(arguments)
    #  print(args)

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(args.infile, 0)
    layer = dataSource.GetLayer()

    clat = 0
    clon = 0
    for feature in layer:
        #  print(feature.GetField("Nom du bar"))
        if (feature.GetField("Nom du bar") == args.name):
            print(feature.GetField("Nom du bar"))
            geom = feature.GetGeometryRef()
            clat = float(geom.Centroid().ExportToWkt().split('(')[1].split(' ')[1].split(')')[0])
            clon = float(geom.Centroid().ExportToWkt().split('(')[1].split(' ')[0])
            #  print("Lat: " + geom.Centroid().ExportToWkt().split('(')[1].split(' ')[1].split(')')[0])
            #  print("Lon: " + geom.Centroid().ExportToWkt().split('(')[1].split(' ')[0])
            break
    layer.ResetReading()

    calt = float(os.popen('gdallocationinfo -valonly -wgs84 %s %s %s' % (args.dem, clon, clat)).read())
    print("Lat: " + str(clat))
    print("Lon: " + str(clon))
    print("Alt: " + str(calt))


    src = osr.SpatialReference();
    src.ImportFromEPSG(4326)
    ds = gdal.Open(args.watermap, gdal.GA_ReadOnly);
    dst = osr.SpatialReference(wkt=ds.GetProjection());
    ct = osr.CoordinateTransformation(src, dst);
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(clat, clon)
    point.Transform(ct)
    print("Coordinates: " + str(point.GetX())+" - "+str(point.GetY()))

    extw = otb.Registry.CreateApplication("ExtractROI")
    extw.SetParameterString("in", args.watermap)
    extw.SetParameterString("mode","radius")
    extw.SetParameterString("mode.radius.unitr", "phy")
    extw.SetParameterFloat("mode.radius.r", float(args.radius))
    extw.SetParameterString("mode.radius.unitc", "phy")
    extw.SetParameterFloat("mode.radius.cx", point.GetX())
    extw.SetParameterFloat("mode.radius.cy", point.GetY())
    extw.SetParameterString("out", os.path.join(args.out, "wmap_extract-"+args.name+".tif"))
    extw.ExecuteAndWriteOutput()

    #  extd = otb.Registry.CreateApplication("ExtractROI")
    #  extd.SetParameterString("in", args.dem)
    #  extd.SetParameterString("mode","fit")
    #  extd.SetParameterString("mode.fit.im", os.path.join(args.out, "wmap_extract-"+args.name+".tif"))
    #  extd.SetParameterString("out", os.path.join(args.out, "dem_extract-"+args.name+".tif"))
    #  extd.ExecuteAndWriteOutput()

    app = otb.Registry.CreateApplication("Superimpose")
    app.SetParameterString("inr", os.path.join(args.out, "wmap_extract-"+args.name+".tif"))
    app.SetParameterString("inm", args.dem)
    app.SetParameterString("out", os.path.join(args.out, "dem_extract-"+args.name+".tif"))
    app.ExecuteAndWriteOutput()


    # Search dam bottom
    extw_bt = otb.Registry.CreateApplication("ExtractROI")
    extw_bt.SetParameterString("in", args.watermap)
    extw_bt.SetParameterString("mode","radius")
    extw_bt.SetParameterString("mode.radius.unitr", "phy")
    extw_bt.SetParameterFloat("mode.radius.r", 500)
    extw_bt.SetParameterString("mode.radius.unitc", "phy")
    extw_bt.SetParameterFloat("mode.radius.cx", point.GetX())
    extw_bt.SetParameterFloat("mode.radius.cy", point.GetY())
    extw_bt.Execute()

    extd_bt = otb.Registry.CreateApplication("Superimpose")
    extd_bt.SetParameterInputImage("inr", extw_bt.GetParameterOutputImage("out"))
    extd_bt.SetParameterString("inm", args.dem)
    extd_bt.Execute()

    bm = otb.Registry.CreateApplication("BandMath")
    bm.AddImageToParameterInputImageList("il",extw_bt.GetParameterOutputImage("out"));
    bm.AddImageToParameterInputImageList("il",extd_bt.GetParameterOutputImage("out"));
    bm.SetParameterString("exp", "( im1b1  > 0.50 ) ? im2b1 : "+str(calt))
    bm.Execute()

    np_surf = bm.GetImageAsNumpyArray('out')
    bt_alt = np.amin(np_surf)
    print("Bottom Alt: " + str(bt_alt))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
