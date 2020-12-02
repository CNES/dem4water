#!/usr/bin/env python3

""" area_mapping.py
Retrieve dam coordinate
"""

import os
import sys
import argparse
from osgeo import ogr

def main(arguments):
    '''Entrypoint'''

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i',
                        '--infile',
                        help="Input file")
    parser.add_argument('-d',
                        '--dam',
                        help="Dam Name")
    parser.add_argument('-o',
                        '--outfile',
                        help="Output file",
                        default=sys.stdout,
                        type=argparse.FileType('w'))

    args = parser.parse_args(arguments)
    #  print(args)

    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(args.infile, 0)
    layer = dataSource.GetLayer()

    for feature in layer:
        #  print(feature.GetField("Nom du bar"))
        if (feature.GetField("Nom du bar") == args.dam):
            print(feature.GetField("Nom du bar"))
            geom = feature.GetGeometryRef()
            print(geom.Centroid().ExportToWkt())
            #  print(geom.Centroid())
            #  print("Lat: " + geom.Centroid().ExportToWkt()[1])
            #  print("Lon: " + geom.Centroid().ExportToWkt()[0])
    layer.ResetReading()



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
