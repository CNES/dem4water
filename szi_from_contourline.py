#!/usr/bin/env python3
'''
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CS Group
:copyright: 2020 CS Group. All rights reserved.
:license: see LICENSE file
:created: 2020
'''

""" szi_from_contourline.py
Measure, for a given range of Z_i, the theorical water surface associated for the dem.
"""

import os
import sys
import logging
import argparse
from osgeo import ogr
from osgeo import osr
from osgeo import gdal
import numpy as np
import otbApplication as otb
import matplotlib.pyplot as plt

def nderiv(y,x):
    "Différence finie, dérivée de la fonction f"
    n = len(y)
    d = np.zeros(n,'d') # virgule flottante à double précision (double)
    # différences de part et d'autre
    # centrées sur les points intérieurs
    for i in range(1,n-1):
        d[i] = (y[i+1]-y[i-1])/(x[i+1]-x[i-1])
    # différences sur un seul côté pour les extrémités
    d[0] = (y[1]-y[0])/(x[1]-x[0])
    d[n-1] = (y[n-1]-y[n-2])/(x[n-1]-x[n-2])
    return d


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
                        type=int,
                        default=500,
                        help="Extract radius (m)")
    parser.add_argument('-s',
                        '--step',
                        type=int,
                        default=5,
                        help="Altitude sampling step")
    parser.add_argument('-o',
                        '--out',
                        help="Output directory")
    parser.add_argument('--debug',
                        action='store_true',
                        help='Activate Debug Mode')
    args = parser.parse_args(arguments)

    logging_format = '%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s'
    if (args.debug is True):
        logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format=logging_format)
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, format=logging_format)
    logging.info("Starting szi_from_contourline.py")

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(args.infile, 0)
    layer = dataSource.GetLayer()

    clat = 0
    clon = 0
    calt = 0

    for feature in layer:
        logging.debug(feature.GetField("Nom du bar"))
        if (feature.GetField("Nom du bar") == args.name):
            geom = feature.GetGeometryRef()
            clat = float(geom.Centroid().ExportToWkt().split('(')[1].split(' ')[1].split(')')[0])
            clon = float(geom.Centroid().ExportToWkt().split('(')[1].split(' ')[0])
            break
    layer.ResetReading()

    calt = float(os.popen('gdallocationinfo -valonly -wgs84 %s %s %s' % (args.dem, clon, clat)).read())
    logging.info("Currently processing: "+ args.name +" [Lat: "+ str(clat) +", Lon: "+ str(clon) +", Alt: "+ str(calt) +"]")


    src = osr.SpatialReference();
    src.ImportFromEPSG(4326)
    ds = gdal.Open(args.dem, gdal.GA_ReadOnly);
    dst = osr.SpatialReference(wkt=ds.GetProjection());
    ct = osr.CoordinateTransformation(src, dst);
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(clat, clon)
    point.Transform(ct)
    logging.debug("Coordinates: " + str(point.GetX())+" - "+str(point.GetY()))

    rad_l = []
    alt_l = []
    #TODO: @param 500
    for r in range(args.radius, 1, -1*args.step):
        ext_l = otb.Registry.CreateApplication("ExtractROI")
        ext_l.SetParameterString("in", args.dem)
        ext_l.SetParameterString("mode","radius")
        ext_l.SetParameterString("mode.radius.unitr", "phy")
        ext_l.SetParameterFloat("mode.radius.r", r)
        ext_l.SetParameterString("mode.radius.unitc", "phy")
        ext_l.SetParameterFloat("mode.radius.cx", point.GetX())
        ext_l.SetParameterFloat("mode.radius.cy", point.GetY())
        ext_l.Execute()

        #  extd_l = otb.Registry.CreateApplication("Superimpose")
        #  extd_l.SetParameterInputImage("inr", extw_l.GetParameterOutputImage("out"))
        #  extd_l.SetParameterString("inm", os.path.join(args.out, "dem_extract-"+args.name+".tif"))
        #  extd_l.Execute()

        np_ext_l = ext_l.GetImageAsNumpyArray('out')
        ext_l_alt = np.amin(np_ext_l)

        rad_l.append(r)
        alt_l.append(ext_l_alt)

        #  logging.debug("@radius= "
                      #  + str(r)
                      #  +"m: local min = "
                      #  + str(ext_l_alt))

    d = nderiv(alt_l, rad_l)
    #  d_2 = nderiv(d, rad_l)

    #  plt.plot(rad_l, alt_l, label='alt(rad)')
    #  plt.plot(rad_l, d*100, label='dérivée')
    #  plt.plot(rad_l, d_2*100, label='Dérivée seconde')
    #  plt.savefig(os.path.join(args.out, "plot.png"))

    found_pdb = False
    rad_pdb = 0
    alt_pdb = 0

    #TODO: @param 0.1
    for i_r, i_a, i_d in zip(rad_l, alt_l, d):
        #  print(i_r, i_a, i_d)
        if (abs(i_d) > 0.15):
            found_pdb = True
            rad_pdb = i_r
            alt_pdb = i_a
            #  logging.info("@radius= "+ str(rad_pdb)
                        #  +"m: local min = "
                        #  + str(alt_pdb))
            break

    if found_pdb is True:
        logging.info("@radius= "+ str(rad_pdb)
                    +"m: local min = "
                    + str(alt_pdb))
    else:
        logging.error("404 - PDB not Found")

    # Retrieve PDB coordinates
    ext = otb.Registry.CreateApplication("ExtractROI")
    ext.SetParameterString("in", args.dem)
    ext.SetParameterString("out", os.path.join(args.out, "dem_pdb.tif"))
    ext.SetParameterString("mode","radius")
    ext.SetParameterString("mode.radius.unitr", "phy")
    ext.SetParameterFloat("mode.radius.r", rad_pdb)
    ext.SetParameterString("mode.radius.unitc", "phy")
    ext.SetParameterFloat("mode.radius.cx", point.GetX())
    ext.SetParameterFloat("mode.radius.cy", point.GetY())
    ext.Execute()

    np_ext = ext.GetImageAsNumpyArray('out')
    ext_alt = np.amin(np_ext)
    print(ext_alt)
    #  indices = np.where(np_ext == [alt_pdb])
    indices = np.where(np_ext == [ext_alt])
    coordinates = zip(indices[0], indices[1])
    print(len(indices[0]))
    #  print(len(coordinates))
    print(indices[0])
    print(indices[1])

    ext.ExecuteAndWriteOutput()
    ds = gdal.Open(os.path.join(args.out, "dem_pdb.tif"))
    xoffset, px_w, rot1, yoffset, px_h, rot2 = ds.GetGeoTransform()

    print(ds.GetGeoTransform())

    posX = px_w * indices[1] + rot1 * indices[0] + xoffset
    posY = rot2 * indices[1] + px_h * indices[0] + yoffset
    #  posX = px_w * indices[0] + rot1 * indices[1] + xoffset
    #  posY = rot2 * indices[0] + px_h * indices[1] + yoffset

    # shift to the center of the pixel
    posX += px_w / 2.0
    posY += px_h / 2.0

    logging.debug("Coordinates: " + str(posX)+" - "+str(posY))

    # get CRS from dataset
    crs = osr.SpatialReference()
    crs.ImportFromWkt(ds.GetProjectionRef())
    # create lat/long crs with WGS84 datum
    crsGeo = osr.SpatialReference()
    crsGeo.ImportFromEPSG(4326) # 4326 is the EPSG id of lat/long crs
    t = osr.CoordinateTransformation(crs, crsGeo)
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(posX, posY)
    point.Transform(t)
    logging.debug("Coordinates: " + str(point.GetX())+" - "+str(point.GetY()))
    #  (pdblat, pdblong) = t.TransformPoint(posX, posY)

    pdbalt = float(os.popen('gdallocationinfo -valonly -wgs84 %s %s %s' % (args.dem, pdblon, pdblat)).read())
    logging.info("Currently processing: "+ args.name +" [pdbLat: "+ str(pdblat) +", Lon: "+ str(pdblon) +", Alt: "+ str(pdbalt) +"]")



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
