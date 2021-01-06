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
from PIL import Image
from osgeo import ogr
from osgeo import osr
from osgeo import gdal
import numpy as np
import otbApplication as otb
import matplotlib.pyplot as plt
import itertools
from math import sqrt,ceil


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


def pixel(dx,dy,ds):
    #  print(ds.GetGeoTransform())
    px = ds.GetGeoTransform()[0]
    py = ds.GetGeoTransform()[3]
    rx = ds.GetGeoTransform()[1]
    ry = ds.GetGeoTransform()[5]
    x = (dx - px) / rx
    y = (dy - py) / ry
    x -= .5
    y -= .5

    return round(x), round(y)

def coord(x,y,ds):
    #  print(ds.GetGeoTransform())
    xoffset, px_w, rot1, yoffset, rot2, px_h = ds.GetGeoTransform()
    dx = px_w * x + rot1 * y + xoffset
    dy = rot2 * x + px_h * y + yoffset
    dx += px_w / 2.0
    dy += px_h / 2.0

    return dx,dy

def points_in_circle(circle, arr):
    "A generator to return all points whose indices are within given circle."
    i0,j0,r = circle
    def intceil(x):
        return int(ceil(x))
    for i in xrange(intceil(i0-r),intceil(i0+r)):
        ri = sqrt(r**2-(i-i0)**2)
        for j in xrange(intceil(j0-ri),intceil(j0+ri)):
            yield arr[i][j]


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
    parser.add_argument('-t',
                        '--tmp',
                        required=True,
                        help="Temporary directory")
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


    geo = osr.SpatialReference();
    geo.ImportFromEPSG(4326)
    ds = gdal.Open(args.dem, gdal.GA_ReadOnly);
    carto = osr.SpatialReference(wkt=ds.GetProjection());
    ct = osr.CoordinateTransformation(geo, carto);
    dam = ogr.Geometry(ogr.wkbPoint)
    dam.AddPoint(clat, clon)
    dam.Transform(ct)
    logging.debug("Coordinates Carto: " + str(dam.GetX())+" - "+str(dam.GetY()))

    if os.path.exists(os.path.join(args.out, "map.json")):
        # It does exist, so remove it
        os.remove(os.path.join(args.out, "map.json"))
    drv = ogr.GetDriverByName( 'GeoJSON' )
    dst_ds = drv.CreateDataSource( os.path.join(args.out, "map.json"))
    dst_layer = dst_ds.CreateLayer('', srs=carto , \
                                   geom_type=ogr.wkbPoint )
    field_defn=ogr.FieldDefn( 'name', ogr.OFTString )
    dst_layer.CreateField( field_defn )
    wkt = "POINT ( %f %f )" % ( float(dam.GetX()), float(dam.GetY()) )
    feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
    p = ogr.CreateGeometryFromWkt( wkt )
    feat.SetGeometryDirectly( p )
    feat.SetField ( "name", "Dam" )
    dst_layer.CreateFeature( feat )
    feat.Destroy()

    rad_l = []
    alt_l = []
    for r in range(args.radius, 1, -1*args.step):
        ext_l = otb.Registry.CreateApplication("ExtractROI")
        ext_l.SetParameterString("in", args.dem)
        ext_l.SetParameterString("mode","radius")
        ext_l.SetParameterString("mode.radius.unitr", "phy")
        ext_l.SetParameterFloat("mode.radius.r", r)
        ext_l.SetParameterString("mode.radius.unitc", "phy")
        ext_l.SetParameterFloat("mode.radius.cx", dam.GetX())
        ext_l.SetParameterFloat("mode.radius.cy", dam.GetY())
        ext_l.Execute()

        np_ext_l = ext_l.GetImageAsNumpyArray('out')
        ext_l_alt = np.amin(np_ext_l)

        rad_l.append(r)
        alt_l.append(ext_l_alt)

        #  logging.debug("@radius= "
                      #  + str(r)
                      #  +"m: local min = "
                      #  + str(ext_l_alt))

    d = nderiv(alt_l, rad_l)

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
            break

    if found_pdb is True:
        logging.debug("@radius= "+ str(rad_pdb)
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
    ext.SetParameterFloat("mode.radius.cx", dam.GetX())
    ext.SetParameterFloat("mode.radius.cy", dam.GetY())
    ext.ExecuteAndWriteOutput()

    np_ext = ext.GetImageAsNumpyArray('out')
    indices = np.where(np_ext == [alt_pdb])
    #TODO: if multiple pdb detected
    if (len(indices[0]) > 1):
        logging.warning("Absolute minimum is not unique on the current area!["
                        + str(len(indices[0]))
                        + "]")
    #  print("Indices Length: "+str(len(indices[0])))
    #  print("X: "+str(indices[1][0]))
    #  print("Y: "+str(indices[0][0]))

    ds = gdal.Open(os.path.join(args.out, "dem_pdb.tif"))
    #  xoffset, px_w, rot1, yoffset, px_h, rot2 = ds.GetGeoTransform()
    #  xoffset, px_w, rot1, yoffset, rot2, px_h = ds.GetGeoTransform()
    #  print(ds.GetGeoTransform())

    #  posX = px_w * indices[1][0] + rot1 * indices[0][0] + xoffset
    #  posY = rot2 * indices[1][0] + px_h * indices[0][0] + yoffset

    # shift to the center of the pixel
    #  posX += px_w / 2.0
    #  posY += px_h / 2.0

    posX, posY = coord(indices[1][0], indices[0][0], ds)
    pX, pY = pixel(posX, posY, ds)

    wkt = "POINT ( %f %f )" % ( float(posX), float(posY) )
    feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
    p = ogr.CreateGeometryFromWkt( wkt )
    feat.SetGeometryDirectly( p )
    feat.SetField ( "name", "PDB" )
    dst_layer.CreateFeature( feat )
    feat.Destroy()

    # get CRS from dataset
    crs = osr.SpatialReference()
    crs.ImportFromWkt(ds.GetProjectionRef())
    # create lat/long crs with WGS84 datum
    crsGeo = osr.SpatialReference()
    crsGeo.ImportFromEPSG(4326) # 4326 is the EPSG id of lat/long crs
    t = osr.CoordinateTransformation(crs, crsGeo)
    pdb = ogr.Geometry(ogr.wkbPoint)
    pdb.AddPoint(float(posX), float(posY))
    pdb.Transform(t)
    pdblat = pdb.GetX()
    pdblon = pdb.GetY()

    pdbalt = float(os.popen('gdallocationinfo -valonly -wgs84 %s %s %s' % (args.dem, pdblon, pdblat)).read())
    logging.debug("Coordinates (pixel): " + str(pX)+" - "+str(pY))
    logging.debug("Coordinates (carto): " + str(posX)+" - "+str(posY))
    logging.debug("Coordinates (latlon): " + str(pdb.GetX())+" - "+str(pdb.GetY()))
    logging.info("PDB detected: "+ args.name +" [pdbLat: "+ str(pdblat) +", pdbLon: "+ str(pdblon) +", pdbAlt: "+ str(pdbalt) +"]")

    # Search for Dam Line
    ext_r = otb.Registry.CreateApplication("ExtractROI")
    ext_r.SetParameterString("in", args.dem)
    ext_r.SetParameterString("out", os.path.join(args.tmp, "extract@"+str(args.radius)+"mFromDam.tif"))
    ext_r.SetParameterString("mode","radius")
    ext_r.SetParameterString("mode.radius.unitr", "phy")
    ext_r.SetParameterFloat("mode.radius.r", args.radius)
    ext_r.SetParameterString("mode.radius.unitc", "phy")
    ext_r.SetParameterFloat("mode.radius.cx", dam.GetX())
    ext_r.SetParameterFloat("mode.radius.cy", dam.GetY())
    ext_r.ExecuteAndWriteOutput()

    r_ds = gdal.Open(os.path.join(args.tmp, "extract@"+str(args.radius)+"mFromDam.tif"), gdal.GA_ReadOnly);
    r_x, r_y = pixel(dam.GetX(), dam.GetY(), r_ds)

    np_r = ext_r.GetImageAsNumpyArray('out')
    (image_size_x, image_size_y) = np_r.shape
    # Disk definition:
    (center_x, center_y) = (r_x, r_y)
    radius = 10

    x_grid, y_grid = np.meshgrid(np.arange(image_size_x),
                                 np.arange(image_size_y))

    prev = 0
    prevprev = 0
    prevpoint1 = ogr.Geometry(ogr.wkbPoint)
    prevpoint1.AddPoint(dam.GetX(), dam.GetY())
    prevpoint2 = ogr.Geometry(ogr.wkbPoint)
    prevpoint2.AddPoint(dam.GetX(), dam.GetY())
    multiline = ogr.Geometry(ogr.wkbMultiLineString)
    for radius in range(5, 100, 2):
        # Array of booleans with the disk shape
        #  disk_out = (((x_grid-center_x)**2 + (y_grid-center_y)**2) <= radius**2) & (((x_grid-center_x)**2 + (y_grid-center_y)**2) <= (radius-1)**2)
        circle = np.logical_and(((x_grid-center_x)**2 + (y_grid-center_y)**2) <= radius**2, ((x_grid-center_x)**2 + (y_grid-center_y)**2) > (radius-2)**2)
        disk_out = ((x_grid-center_x)**2 + (y_grid-center_y)**2) <= radius**2
        disk_in = ((x_grid-center_x)**2 + (y_grid-center_y)**2) <= (radius-1)**2

        points_on_circle = np_r[circle]

        #  print(len(points_on_circle))
        if (len(points_on_circle)>prevprev):
            prevprev = prev
            prev = len(points_on_circle)

            #  print(np.amax(points_on_circle))
            # You can also use masked arrays:
            #  circle_2 = np.ma.masked_array(disk_out, ~disk_in)
            #  print(np.amax(circle_2))
            #  print(np.amax(points_on_circle))
            masked_tmp = np.where(disk_out, np_r, 0)
            masked = np.where(~disk_in, masked_tmp, 0)
            im_masked = Image.fromarray(masked)
            #  im_masked.save(os.path.join(args.tmp, "circle@"+str(radius)+".tif"))

            #  print(np.amax(im_masked))
            #  print(np.amax(im_masked))

            #  l_ds = gdal.Open(os.path.join(args.out, "circle@"+str(radius)+".tif"))

            l_indices = np.where(masked == [np.amax(masked)])
            l_posX, l_posY = coord(l_indices[1][0], l_indices[0][0], r_ds)
            l_pX, l_pY = pixel(l_posX, l_posY, r_ds)

            currpoint = ogr.Geometry(ogr.wkbPoint)
            currpoint.AddPoint(l_posX, l_posY)

            #TODO: if multiple max detected
            if (len(l_indices[0]) > 1):
                logging.warning("Absolute maximum is not unique on the current circle!["
                                + str(len(l_indices[0]))
                                + "]")
            #  print("Indices Length: "+str(len(l_indices[0])))
            #  print("X: "+str(l_indices[1][0]))
            #  print("Y: "+str(l_indices[0][0]))

            # Add Circle absolute max to json
            l_wkt = "POINT ( %f %f )" % ( float(l_posX), float(l_posY) )
            l_feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
            l_p = ogr.CreateGeometryFromWkt( l_wkt )
            l_feat.SetGeometryDirectly( l_p )
            l_feat.SetField ( "name", str(radius)+"/1" )
            dst_layer.CreateFeature( l_feat )
            l_feat.Destroy()

            #  logging.debug("First Point - Coordinates (pixel): " + str(l_pX)+" - "+str(l_pY))
            #  logging.debug("First Point - Coordinates (carto): " + str(l_posX)+" - "+str(l_posY))

            # Search for opposite relative max: mask half of the circle
            done_half = ((x_grid-l_indices[1][0])**2 + (y_grid-l_indices[0][0])**2) <= (radius*1.7)**2
            half_masked = np.where(~done_half, masked, 0)
            im_hmasked = Image.fromarray(half_masked)
            #  im_hmasked.save(os.path.join(args.tmp, "half_circle@"+str(radius)+".tif"))

            l_indices = np.where(half_masked == [np.amax(half_masked)])
            l_posX, l_posY = coord(l_indices[1][0], l_indices[0][0], r_ds)
            l_pX, l_pY = pixel(l_posX, l_posY, r_ds)

            nextpoint = ogr.Geometry(ogr.wkbPoint)
            nextpoint.AddPoint(l_posX, l_posY)

            #TODO: if multiple max detected
            if (len(l_indices[0]) > 1):
                logging.warning("Absolute maximum is not unique on the current circle!["
                                + str(len(l_indices[0]))
                                + "]")
            #  print("Indices Length: "+str(len(l_indices[0])))
            #  print("X: "+str(l_indices[1][0]))
            #  print("Y: "+str(l_indices[0][0]))

            # Add Circle absolute max to json
            l_wkt = "POINT ( %f %f )" % ( float(l_posX), float(l_posY) )
            l_feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
            l_p = ogr.CreateGeometryFromWkt( l_wkt )
            l_feat.SetGeometryDirectly( l_p )
            l_feat.SetField ( "name", str(radius)+"/1" )
            dst_layer.CreateFeature( l_feat )
            l_feat.Destroy()

            #  logging.debug("Second Point - Coordinates (pixel): " + str(l_pX)+" - "+str(l_pY))
            #  logging.debug("Second Point - Coordinates (carto): " + str(l_posX)+" - "+str(l_posY))

            distance1 = currpoint.Distance(prevpoint1)
            distance2 = currpoint.Distance(prevpoint2)
            #  logging.debug("distance1: " + str(distance1))
            #  logging.debug("distance2: " + str(distance2))

            if (distance1 <= distance2):
                logging.debug("First")
                line1 = ogr.Geometry(ogr.wkbLineString)
                line1.AddPoint(prevpoint1.GetX(),prevpoint1.GetY())
                line1.AddPoint(currpoint.GetX(),currpoint.GetY())
                multiline.AddGeometry(line1)
                prevpoint1 = currpoint

                line2 = ogr.Geometry(ogr.wkbLineString)
                line2.AddPoint(prevpoint2.GetX(),prevpoint2.GetY())
                line2.AddPoint(nextpoint.GetX(),nextpoint.GetY())
                multiline.AddGeometry(line2)
                prevpoint2 = nextpoint

            else:
                logging.debug("Second")
                line1 = ogr.Geometry(ogr.wkbLineString)
                line1.AddPoint(prevpoint2.GetX(),prevpoint2.GetY())
                line1.AddPoint(currpoint.GetX(),currpoint.GetY())
                multiline.AddGeometry(line1)
                prevpoint2 = currpoint

                line2 = ogr.Geometry(ogr.wkbLineString)
                line2.AddPoint(prevpoint1.GetX(),prevpoint1.GetY())
                line2.AddPoint(nextpoint.GetX(),nextpoint.GetY())
                multiline.AddGeometry(line2)
                prevpoint1 = nextpoint


            # dict2wkbMultiLineString
            #  multiline = ogr.Geometry(ogr.wkbMultiLineString)
            #  for i in itertools.combinations(pointDict.values(), 2):
            #      point1 = ogr.Geometry(ogr.wkbPoint)
            #      point1.AddPoint(i[0][0],i[0][1])
            #      point2 = ogr.Geometry(ogr.wkbPoint)
            #      point2.AddPoint(i[1][0],i[1][1])
            #
            #      distance = point1.Distance(point2)
            #
            #      if distance < maxDistance:
            #          line = ogr.Geometry(ogr.wkbLineString)
            #          line.AddPoint(i[0][0],i[0][1])
            #          line.AddPoint(i[1][0],i[1][1])
            #          multiline.AddGeometry(line)


        else:
            logging.debug("Search area outside of ROI")
            break

        #  from PIL import Image
        #  im = Image.fromarray(A)
        #  im.save("your_file.jpeg")

    shpDriver = ogr.GetDriverByName( 'GeoJSON' )
    if os.path.exists(os.path.join(args.out, "line.json")):
        shpDriver.DeleteDataSource(os.path.join(args.out, "line.json"))
    outDataSource = shpDriver.CreateDataSource(os.path.join(args.out, "line.json"))
    outLayer = outDataSource.CreateLayer('', srs=carto , \
                                         geom_type=ogr.wkbMultiLineString )
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(multiline)
    outLayer.CreateFeature(outFeature)

    #  for r in range(1, args.radius, args.step):
    #      ext_l = otb.Registry.CreateApplication("ExtractROI")
    #      ext_l.SetParameterString("in", args.dem)
    #      ext_l.SetParameterString("out", os.path.join(args.tmp, "extract@"+str(r)+"mFromDam.tif"))
    #      ext_l.SetParameterString("mode","radius")
    #      ext_l.SetParameterString("mode.radius.unitr", "phy")
    #      ext_l.SetParameterFloat("mode.radius.r", r)
    #      ext_l.SetParameterString("mode.radius.unitc", "phy")
    #      ext_l.SetParameterFloat("mode.radius.cx", dam.GetX())
    #      ext_l.SetParameterFloat("mode.radius.cy", dam.GetY())
    #      ext_l.ExecuteAndWriteOutput()
    #
    #      # Dam location in pixel indicies
    #      l_ds = gdal.Open(os.path.join(args.tmp, "extract@"+str(r)+"mFromDam.tif"), gdal.GA_ReadOnly);
    #      l_x, l_y = pixel(dam.GetX(), dam.GetY(), l_ds)
    #      logging.debug("Coordinates : [" + str(dam.GetX())+", "+str(dam.GetY())
    #                    +"] --> ["+str(l_x)+", "+str(l_y)+"]")
    #
    #      np_ext_l = ext_l.GetImageAsNumpyArray('out')
    #      ext_l_max = np.amax(np_ext_l)
    #      #  logging.debug("@radius= "+ str(r)
    #                  #  +"m: local max = "
    #                  #  + str(ext_l_max))


    # close json
    dst_ds = None


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
