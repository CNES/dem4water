#!/usr/bin/env python3
'''
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CS Group
:copyright: 2020 CS Group. All rights reserved.
:license: see LICENSE file
:created: 2020
'''

""" szi_from_watermap.py
Measure, for a given range of Z_i, the water surface associated.
"""

import os
import sys
import math
import logging
import argparse
import numpy as np
import otbApplication as otb
from osgeo import gdal


def main(arguments):
    '''Entrypoint'''

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-w',
                        '--watermap',
                        required=True,
                        help="Input water map file")
    parser.add_argument('-d',
                        '--dem',
                        required=True,
                        help="Input DEM")
    parser.add_argument('--zmin',
                        type=float,
                        required=True,
                        help="Bottom altitude")
    parser.add_argument('--zmax',
                        type=float,
                        required=True,
                        help="Starting altitude")
    parser.add_argument('-t',
                        '--tmp',
                        required=True,
                        help="Temporary directory")
    parser.add_argument('-s',
                        '--step',
                        type=int,
                        default=5,
                        help="Altitude sampling step")
    parser.add_argument('-o',
                        '--outfile',
                        help="Output file")
    parser.add_argument('--loglevel',
                        help="Log Level",
                        default="INFO")

    args = parser.parse_args(arguments)

    numeric_level = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.loglevel)
    logging_format = '%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s'
    logging.basicConfig(stream=sys.stdout, level=numeric_level, format=logging_format)
    logging.info("Starting szi_from_watermap.py")

    raster = gdal.Open(args.watermap)
    gt =raster.GetGeoTransform()
    pixelSizeX = gt[1]
    pixelSizeY =-gt[5]
    logging.debug("pixelSizeX: " + str(pixelSizeX) +" - pixelSizeY: " + str(pixelSizeY))

    zi = []
    S_zi = []

    last_result = -1
    msk = otb.Registry.CreateApplication("BandMath")
    seg = otb.Registry.CreateApplication("Segmentation")
    for i in range(math.ceil(args.zmax), math.floor(args.zmin), -1 * args.step):

        msk.SetParameterStringList("il", [args.watermap, args.dem])
        msk.SetParameterString("out", os.path.join(args.tmp, "mask_S_Z-"+str(i)+".tif"))
        msk.SetParameterOutputImagePixelType("out", otb.ImagePixelType_uint8)
        msk.SetParameterString("exp", "( im1b1  > 0.10 ) ? ((im2b1 < "+str(i)+") ? 1 : 0 ) : 0")
        #  msk.ExecuteAndWriteOutput()
        msk.Execute()

        seg.SetParameterInputImage("in", msk.GetParameterOutputImage("out"))
        seg.SetParameterString("mode","raster")
        seg.SetParameterString("mode.raster.out", os.path.join(args.tmp, "tmplbmap_S_Z-"+str(i)+".tif"))
        seg.SetParameterOutputImagePixelType("mode.raster.out", 3)
        seg.SetParameterString("filter","cc")
        seg.SetParameterString("filter.cc.expr", "(distance < 2) * (p1b1 > 0)")
        #  seg.ExecuteAndWriteOutput()
        seg.Execute()

        bm = otb.Registry.CreateApplication("BandMath")
        bm.AddImageToParameterInputImageList("il",msk.GetParameterOutputImage("out"));
        bm.AddImageToParameterInputImageList("il",seg.GetParameterOutputImage("mode.raster.out"));
        bm.SetParameterString("out", os.path.join(args.tmp, "lbmap_S_Z-"+str(i)+".tif"))
        bm.SetParameterString("exp", "( im1b1 * im2b1 )")
        #  bm.ExecuteAndWriteOutput()
        bm.Execute()

        np_surf = bm.GetImageAsNumpyArray('out')
        max_label = np.amax(np_surf)
        logging.debug("pixelSizeX: " + str(pixelSizeX) +" - pixelSizeY: " + str(pixelSizeY))
        logging.debug("max_label: "+str(max_label))
        logging.debug("max_label: "+str(int(max_label)))

        best_surf = 0
        best_label = 0
        for l in range(1, int(max_label+1)):
            results_surf = np.count_nonzero((np_surf == [l]))
            if (results_surf > best_surf):
                logging.debug("new best_surf: "+str(results_surf)+" [label: "+str(l)+"]")
                best_label = l
                best_surf = results_surf

        surf_m2 = best_surf*pixelSizeX*pixelSizeY
        logging.info("@"+ str(i) +"m: "+ str(surf_m2)+" m2.")
        logging.debug("@" + str(i) + "m: best_surf: "+str(best_surf)+" pixels "+"[label: "+str(best_label)+"]")
        logging.debug("--> "+str(surf_m2)+" m2 ")
        logging.debug("--> "+str(surf_m2*0.000001)+" km2 ")

        if ( best_surf < last_result/2.0 ):
            logging.info("Estimated surface dropped --> Stop here." )
            break
        else:
            last_result = best_surf
            zi.append(i)
            S_zi.append(surf_m2)

    zi.append(args.zmin)
    S_zi.append(0.0)
    data = np.column_stack((zi, S_zi))
    np.savetxt(args.outfile, data)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
