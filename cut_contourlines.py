#!/usr/bin/env python3
'''
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CS Group
:copyright: 2020 CS Group. All rights reserved.
:license: see LICENSE file
:created: 2020
'''

""" cut_contourliness.py
Cut contour lines based on the cutline to estimate the virtual water surface.
"""

import os
import sys
import logging
import argparse
import json
import shapely
import shapely.wkt
from shapely.geometry import shape
from shapely.geometry import Point, Polygon
from shapely.ops import split
from osgeo import ogr
from osgeo import osr
from osgeo import gdal


def main(arguments):
    '''Entrypoint'''

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i',
                        '--info',
                        help="daminfo.json file")
    parser.add_argument('-c',
                        '--cut',
                        help="cutline.json file")
    parser.add_argument('-l',
                        '--level',
                        help="contournline.json file")
    parser.add_argument('-d',
                        '--dem',
                        help="Input DEM")
    parser.add_argument('-t',
                        '--tmp',
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
    logging.info("Starting cut_contourlines.py")

    # Silence Mathplotlib related debug messages (font matching)
    logging.getLogger('matplotlib').setLevel(logging.ERROR)

    ds = gdal.Open(args.dem, gdal.GA_ReadOnly);
    carto = osr.SpatialReference(wkt=ds.GetProjection());

    drv = ogr.GetDriverByName( 'GeoJSON' )
    if os.path.exists("splitted.json"):
        os.remove("splitted.json")
    dst_ds = drv.CreateDataSource("splitted.json")
    #  if os.path.exists(os.path.join(args.out, dam_path +"_daminfo.json")):
        #  os.remove(os.path.join(args.out, dam_path +"_daminfo.json"))
    #  dst_ds = drv.CreateDataSource( os.path.join(args.out, dam_path +"_daminfo.json"))
    dst_layer = dst_ds.CreateLayer('', srs=carto , \
                                   geom_type=ogr.wkbPolygon )
    field_defn_id=ogr.FieldDefn( 'ID', ogr.OFTString )
    field_defn=ogr.FieldDefn( 'level', ogr.OFTString )
    dst_layer.CreateField( field_defn_id )
    dst_layer.CreateField( field_defn )

    # load GeoJSON file containing info
    with open(args.info) as i:
        jsi = json.load(i)

    #  dam=""
    #  pdb=""
    for feature in jsi['features']:
        if feature['properties']['name'] == 'Dam':
            dam = shape(feature['geometry'])

        if feature['properties']['name'] == 'PDB':
            pdb = shape(feature['geometry'])

    # load GeoJSON file containing cutline
    with open(args.cut) as c:
        jsc = json.load(c)

    for feature in jsc['features']:
        lines = shape(feature['geometry'])
        line = shapely.ops.linemerge(lines)
        print(line)

    # load GeoJSON file containing contour lines
    with open(args.level) as l:
        jsl = json.load(l)

    r_id = 1
    for feature in jsl['features']:
        print('--')
        print(feature['properties']['ID'])
        #  print('--')
        #  print('--')
        #  print(feature)
        level = shape(feature['geometry'])
        results = split(level, line)
        max_area = 0
        for p in results:
            #  print(p.contains(dam))
            if not p.contains(pdb) and (p.area > max_area):
                print("Contains PDB: " + str(p.contains(pdb)))
                #  if p.touches(dam):
                    #  print("Contains DAM: " + str(p.contains(dam)))
                max_area = p.area
                print("Area: " + str(p.area) + "m2")
                r_feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
                r_p = ogr.CreateGeometryFromWkt( p.wkt )
                r_feat.SetGeometryDirectly( r_p )
                r_feat.SetField ( "ID", str(r_id) )
                r_feat.SetField ( "level", feature['properties']['ID'] )

        dst_layer.CreateFeature( r_feat )
        r_feat.Destroy()
        r_id =r_id + 1


    print(dam.wkt)
    print(pdb.wkt)
    #  print(results)

    #  driver = ogr.GetDriverByName("ESRI Shapefile")
    #  dataSource = driver.Open(args.level, 0)
    #  layer = dataSource.GetLayer()
    #
    #  for feature in layer:
    #      geom = feature.GetGeometryRef()
    #      #  print(geom.Centroid().ExportToWkt())
    #      P = shapely.wkt.loads(geom.ExportToWkt())
    #      #  print(P)




        #  polygon = shape(feature['geometry'])

    # load GeoJSON file containing levels
    #  with open(args.level) as f:
        #  levels = json.load(f)

    #  levels = ogr.Open(args.level)
    #  shape = levels.GetLayer(0)
    #  for feature in layer:
            #  print feature.GetField("STATE_NAME")
    #first feature of the shapefile
    #  for feature in shape:
    #  feature = shape.GetFeature(0)
    #  first = feature.ExportToJson()
    #  print(first) # (GeoJSON format)


    #  for feat in first['geometry']:
        #  polygon = Polygon(feat)
        #  print(type(polygon))
#

    #  from shapely.geometry import shape
    #  shp_geom = Polygon(first['geometry'])
    #  print(shp_geom)
    #  for levelfeature in js['features']:
        #  polygon = shape(feature['geometry'])
        #  if polygon.contains(point):
            #  print 'Found containing polygon:', feature



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
