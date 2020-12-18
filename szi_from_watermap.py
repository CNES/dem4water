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
    parser.add_argument('--debug',
                        action='store_true',
                        help='Activate Debug Mode')

    args = parser.parse_args(arguments)
    #  print(args)

    logging_format = '%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s'
    if (args.debug is True):
        logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format=logging_format)
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, format=logging_format)
    logging.info("Starting szi_from_watermap.py")

    raster = gdal.Open(args.watermap)
    gt =raster.GetGeoTransform()
    pixelSizeX = gt[1]
    pixelSizeY =-gt[5]
    logging.debug("pixelSizeX: " + str(pixelSizeX) +" - pixelSizeY: " + str(pixelSizeY))

    zi = []
    S_zi = []

    zi_up = []
    S_zi_up= []

    #  wshed = otb.Registry.CreateApplication("Segmentation")
    #  wshed.SetParameterString("in", args.dem)
    #  wshed.SetParameterString("mode","raster")
    #  wshed.SetParameterString("mode.raster.out", os.path.join(args.tmp, "wshed.tif"))
    #  wshed.SetParameterOutputImagePixelType("mode.raster.out", 3)
    #  wshed.SetParameterString("filter","watershed")
    #  wshed.SetParameterFloat("filter.watershed.threshold", 0.01)
    #  wshed.SetParameterFloat("filter.watershed.level", 0.25)
    #  wshed.ExecuteAndWriteOutput()


    msk = otb.Registry.CreateApplication("BandMath")
    seg = otb.Registry.CreateApplication("Segmentation")
    bm = otb.Registry.CreateApplication("BandMath")

    # First mask @zmax
    msk.SetParameterStringList("il", [args.watermap, args.dem])
    msk.SetParameterString("out", os.path.join(args.tmp, "mask_@"+str(args.zmax)+"m.tif"))
    msk.SetParameterOutputImagePixelType("out", otb.ImagePixelType_uint8)
    msk.SetParameterString("exp", "( im1b1  > 0.05 ) ? ((im2b1 < "+str(args.zmax)+") ? 1 : 0 ) : 0")
    #  msk.SetParameterString("exp", "( im1b1  > 0.10 ) ? ((im2b1 < "+str(args.zmax)+") ? 1 : 0 ) : 0")
    #  msk.SetParameterString("exp", "( im1b1  > 0.10 ) ? 1 : 0 ")
    if os.path.exists(os.path.join(args.tmp, "mask_@"+str(args.zmax)+"m.tif")) is False:
        msk.ExecuteAndWriteOutput()
    #  else:
        #  msk.Execute()
        #  seg.SetParameterInputImage("in", msk.GetParameterOutputImage("out"))

    mdil = otb.Registry.CreateApplication("BinaryMorphologicalOperation")
    #  mdil.SetParameterInputImage("in", kp.GetParameterOutputImage("out"))
    mdil.SetParameterString("in", os.path.join(args.tmp, "mask_@"+str(args.zmax)+"m.tif"))
    mdil.SetParameterString("out", os.path.join(args.tmp, "maskdil_@"+str(args.zmax)+"m.tif"))
    mdil.SetParameterInt("channel", 1)
    mdil.SetParameterInt("xradius", 5)
    mdil.SetParameterInt("yradius", 5)
    mdil.SetParameterString("filter","dilate")
    mdil.Execute()

    mbm = otb.Registry.CreateApplication("BandMath")
    mbm.SetParameterStringList("il", [args.dem])
    mbm.AddImageToParameterInputImageList("il",mdil.GetParameterOutputImage("out"));
    mbm.SetParameterString("out", os.path.join(args.tmp, "lbmap_@"+str(args.zmax)+"m.tif"))
    mbm.SetParameterString("exp", "( (im1b1 < "+str(args.zmax)+") ? im1b1 * im2b1 : 0 )")
    mbm.Execute()

    wshed = otb.Registry.CreateApplication("Segmentation")
    wshed.SetParameterInputImage("in", mbm.GetParameterOutputImage("out"))
    #  wshed.SetParameterString("in", os.path.join(args.tmp, "mask_@"+str(args.zmax)+"m.tif"))
    wshed.SetParameterString("mode","raster")
    wshed.SetParameterString("mode.raster.out", os.path.join(args.tmp, "wshed.tif"))
    wshed.SetParameterOutputImagePixelType("mode.raster.out", 3)
    wshed.SetParameterString("filter","watershed")
    wshed.SetParameterFloat("filter.watershed.threshold", 0.01)
    wshed.SetParameterFloat("filter.watershed.level", 0.25)
    wshed.ExecuteAndWriteOutput()

    seg.SetParameterString("in", os.path.join(args.tmp, "mask_@"+str(args.zmax)+"m.tif"))
    seg.SetParameterString("mode","raster")
    seg.SetParameterString("mode.raster.out", os.path.join(args.tmp, "tmplbmap_@"+str(args.zmax)+"m.tif"))
    seg.SetParameterOutputImagePixelType("mode.raster.out", 3)
    seg.SetParameterString("filter","cc")
    seg.SetParameterString("filter.cc.expr", "(distance < 20) * (p1b1 > 0)")
    if os.path.exists(os.path.join(args.tmp, "tmplbmap_@"+str(args.zmax)+"m.tif")) is False:
        seg.ExecuteAndWriteOutput()
    #  else:
        #  seg.Execute()

    #  #  bm.AddImageToParameterInputImageList("il",msk.GetParameterOutputImage("out"));
    #  #  bm.AddImageToParameterInputImageList("il",seg.GetParameterOutputImage("mode.raster.out"));
    #  bm.SetParameterStringList("il", [os.path.join(args.tmp, "mask_@"+str(args.zmax)+"m.tif"), os.path.join(args.tmp, "tmplbmap_@"+str(args.zmax)+"m.tif")])
    #  bm.SetParameterString("out", os.path.join(args.tmp, "lbmap_@"+str(args.zmax)+"m.tif"))
    #  bm.SetParameterString("exp", "( im1b1 * im2b1 )")
    #  if os.path.exists(os.path.join(args.tmp, "lbmap_@"+str(args.zmax)+"m.tif")) is False:
    #      bm.ExecuteAndWriteOutput()
    #  else:
    #      bm.Execute()

    best_surf = 0
    best_label = 0
    #  np_surf = bm.GetImageAsNumpyArray('out')
    np_surf = wshed.GetImageAsNumpyArray('mode.raster.out')
    max_label = np.amax(np_surf)
    for l in range(2, int(max_label+1)):
        results_surf = np.count_nonzero((np_surf == [l]))
        if (results_surf > best_surf):
            logging.debug("new best_surf: "+str(results_surf)+" [label: "+str(l)+"]")
            best_label = l
            best_surf = results_surf

    surf_m2 = best_surf*pixelSizeX*pixelSizeY
    logging.info("@"+ str(args.zmax) +"m: "+ str(surf_m2)+" m2.")
    logging.debug("@" + str(args.zmax) + "m: best_surf: "+str(best_surf)+" pixels "+"[label: "+str(best_label)+"]")
    logging.debug("--> "+str(surf_m2)+" m2 ")
    logging.debug("--> "+str(surf_m2*0.000001)+" km2 ")

    zi.append(args.zmax)
    S_zi.append(surf_m2)

    #  alt_on_surf = otb.Registry.CreateApplication("BandMath")
    #  alt_on_surf.SetParameterStringList("il", [args.dem, os.path.join(args.tmp, "wshed.tif")])
    #  alt_on_surf.AddImageToParameterInputImageList("il",bm.GetParameterOutputImage("out"));
    #  #  alt_on_surf.SetParameterString("out", os.path.join(args.tmp, "lbmap_S_Z-"+str(args.zmax)+".tif"))
    #  alt_on_surf.SetParameterString("exp", "( im2b1 == "+ str(best_label) +") ? im1b1 : 0")
    #  alt_on_surf.Execute()
    #
    #  np_alt_on_surf = alt_on_surf.GetImageAsNumpyArray('out')
    #  max_alt = np.amax(np_alt_on_surf)
    #  logging.info("Water surface max altitude: "+ str(max_alt) +"m.")

    kp = otb.Registry.CreateApplication("BandMath")
    #  kp.AddImageToParameterInputImageList("il", bm.GetParameterOutputImage("out"));
    kp.AddImageToParameterInputImageList("il", wshed.GetParameterOutputImage("mode.raster.out"));
    kp.SetParameterOutputImagePixelType("out", otb.ImagePixelType_uint8)
    #  kp.SetParameterString("out", os.path.join(args.tmp, "kp_@"+str(max_alt)+"m.tif"))
    kp.SetParameterString("out", os.path.join(args.tmp, "kp_@"+str(args.zmax)+"m.tif"))
    kp.SetParameterString("exp", "( im1b1 == "+ str(best_label) +" ) ? 1 : 0")
    kp.Execute()

    #  dil = otb.Registry.CreateApplication("BinaryMorphologicalOperation")
    #  dil.SetParameterInputImage("in", kp.GetParameterOutputImage("out"))
    #  dil.SetParameterString("out", os.path.join(args.tmp, "maskdil_@"+str(max_alt)+"m.tif"))
    #  dil.SetParameterInt("channel", 1)
    #  dil.SetParameterInt("xradius", 25)
    #  dil.SetParameterInt("yradius", 25)
    #  dil.SetParameterString("filter","dilate")
    #  dil.Execute()
    #
    #  sth = otb.Registry.CreateApplication("BandMath")
    #  sth.SetParameterStringList("il", [args.dem])
    #  sth.AddImageToParameterInputImageList("il",dil.GetParameterOutputImage("out"));
    #  sth.AddImageToParameterInputImageList("il",kp.GetParameterOutputImage("out"));
    #  sth.SetParameterString("out", os.path.join(args.tmp, "maskfil_@"+str(max_alt)+"m.tif"))
    #  sth.SetParameterString("exp", "( im2b1 * (im1b1 >= "+str(max_alt)+") * (im1b1 < "+str(args.zmax)+") ) ? 1 : im3b1 ")
    #  sth.Execute()

    #  np_surf_corr = sth.GetImageAsNumpyArray('out')
    #  results_surf_corr = np.count_nonzero((np_surf_corr == [1]))
    #  surf_m2 = results_surf_corr*pixelSizeX*pixelSizeY
    #  logging.info("@"+ str(args.zmax) +"m: "+ str(surf_m2)+" m2.")
    #  logging.debug("@" + str(args.zmax) + "m: corrected_surf: "+str(results_surf_corr)+" pixels")
    #  logging.debug("--> "+str(surf_m2)+" m2 ")
    #  logging.debug("--> "+str(surf_m2*0.000001)+" km2 ")

    #  bm.ExecuteAndWriteOutput()
    #  sth.ExecuteAndWriteOutput()
    kp.ExecuteAndWriteOutput()
    #  dil.ExecuteAndWriteOutput()

    last_result = -1
    last_valid_alt = math.ceil(args.zmax)
    #  for i in range(math.ceil(args.zmax - args.step), math.floor(args.zmin), -1 * args.step):
    for i in range(math.ceil(args.zmax), math.floor(args.zmin), -1 * args.step):

        msk = otb.Registry.CreateApplication("BandMath")
        #  msk.SetParameterStringList("il", [args.watermap, args.dem])
        msk.SetParameterStringList("il", [args.dem, os.path.join(args.tmp, "kp_@"+str(args.zmax)+"m.tif")])
        #  msk.SetParameterStringList("il", [args.dem])
        #  msk.AddImageToParameterInputImageList("il",kp.GetParameterOutputImage("out"));
        msk.SetParameterString("out", os.path.join(args.tmp, "mask_S_Z-"+str(i)+".tif"))
        msk.SetParameterOutputImagePixelType("out", otb.ImagePixelType_uint8)
        msk.SetParameterString("exp", "(im1b1 < "+str(i)+") ? im2b1 : 0 ")
        #  msk.SetParameterString("exp", "(im1b1 < "+str(i)+") ? im1b1 * im2b1 : 0 ")
        msk.ExecuteAndWriteOutput()
        #  msk.Execute()

        #  seg = otb.Registry.CreateApplication("Segmentation")
        #  seg.SetParameterInputImage("in", msk.GetParameterOutputImage("out"))
        #  seg.SetParameterString("mode","raster")
        #  seg.SetParameterString("mode.raster.out", os.path.join(args.tmp, "tmplbmap_S_Z-"+str(i)+".tif"))
        #  seg.SetParameterOutputImagePixelType("mode.raster.out", 3)
        #  seg.SetParameterString("filter","cc")
        #  seg.SetParameterString("filter.cc.expr", "(distance < 2) * (p1b1 > 0)")
        #  seg.Execute()

        #  bm = otb.Registry.CreateApplication("BandMath")
        #  bm.AddImageToParameterInputImageList("il",msk.GetParameterOutputImage("out"));
        #  bm.AddImageToParameterInputImageList("il",seg.GetParameterOutputImage("mode.raster.out"));
        #  bm.SetParameterString("out", os.path.join(args.tmp, "lbmap_S_Z-"+str(i)+".tif"))
        #  bm.SetParameterString("exp", "( im1b1 * im2b1 )")
        #  bm.Execute()

        #  np_surf = bm.GetImageAsNumpyArray('out')
        #  max_label = np.amax(np_surf)
        #  logging.debug("pixelSizeX: " + str(pixelSizeX) +" - pixelSizeY: " + str(pixelSizeY))
        #  logging.debug("max_label: "+str(max_label))
        #  logging.debug("max_label: "+str(int(max_label)))

        #  best_surf = 0
        #  best_label = 0
        #  for l in range(1, int(max_label+1)):
            #  results_surf = np.count_nonzero((np_surf == [l]))
            #  if (results_surf > best_surf):
                #  logging.debug("new best_surf: "+str(results_surf)+" [label: "+str(l)+"]")
                #  best_label = l
                #  best_surf = results_surf

        np_surf = msk.GetImageAsNumpyArray('out')
        best_surf = np.count_nonzero((np_surf == [1]))

        surf_m2 = best_surf*pixelSizeX*pixelSizeY
        logging.info("@"+ str(i) +"m: "+ str(surf_m2)+" m2.")
        #  logging.debug("@" + str(i) + "m: best_surf: "+str(best_surf)+" pixels "+"[label: "+str(best_label)+"]")
        #  logging.debug("--> "+str(surf_m2)+" m2 ")
        logging.debug("--> "+str(surf_m2*0.000001)+" km2 ")

        if ( best_surf < last_result/2.0 ):
            logging.info("Estimated surface dropped --> Stop here." )
            break
        else:
            last_result = best_surf
            last_valid_alt = i
            zi.append(i)
            S_zi.append(surf_m2)

        if (args.debug is True):
            msk.ExecuteAndWriteOutput()
            #  seg.ExecuteAndWriteOutput()
            #  bm.ExecuteAndWriteOutput()


    # Ascending Search
    init = True
    prev_alt = last_valid_alt
    for i in range(last_valid_alt, math.floor(args.zmax + 10*args.step), args.step):
        adil = otb.Registry.CreateApplication("BinaryMorphologicalOperation")
        #  adil.SetParameterString("in", os.path.join(args.tmp, "kp_@"+str(args.zmax)+"m.tif"))
        if (init is True):
            adil.SetParameterString("in", os.path.join(args.tmp, "mask_S_Z-"+str(prev_alt)+".tif"))
        else:
            adil.SetParameterString("in", os.path.join(args.tmp, "up_mask_S_Z-"+str(prev_alt)+".tif"))

        adil.SetParameterString("out", os.path.join(args.tmp, "up_maskdil_S_Z-"+str(i)+".tif"))
        adil.SetParameterOutputImagePixelType("out", otb.ImagePixelType_uint8)
        adil.SetParameterInt("channel", 1)
        adil.SetParameterInt("xradius", 20)
        adil.SetParameterInt("yradius", 20)
        adil.SetParameterString("filter","dilate")
        adil.ExecuteAndWriteOutput()
        #  adil.Execute()

        amsk = otb.Registry.CreateApplication("BandMath")
        if (init is True):
            amsk.SetParameterStringList("il", [args.dem,
                                               os.path.join(args.tmp, "mask_S_Z-"+str(prev_alt)+".tif"),
                                               os.path.join(args.tmp, "up_maskdil_S_Z-"+str(i)+".tif")])
            init = False
        else:
            amsk.SetParameterStringList("il", [args.dem,
                                               #  os.path.join(args.tmp, "up_mask_S_Z-"+str(prev_alt)+".tif"),
                                               os.path.join(args.tmp, "up_masktmp_S_Z-"+str(prev_alt)+".tif"),
                                               os.path.join(args.tmp, "up_maskdil_S_Z-"+str(i)+".tif")])
        #  amsk.AddImageToParameterInputImageList("il",kp.GetParameterOutputImage("out"));
        amsk.SetParameterString("out", os.path.join(args.tmp, "up_masktmp_S_Z-"+str(i)+".tif"))
        amsk.SetParameterOutputImagePixelType("out", otb.ImagePixelType_uint8)
        amsk.SetParameterString("exp", "(im2b1 == 1) ? im2b1 : (((im3b1 == 1)*(im1b1 > "+str(prev_alt)+")*(im1b1 <= "+str(i)+") ) ? 1 : 0 ) ")
        amsk.ExecuteAndWriteOutput()
        #  msk.Execute()

        aclo = otb.Registry.CreateApplication("BinaryMorphologicalOperation")
        aclo.SetParameterString("in", os.path.join(args.tmp, "up_masktmp_S_Z-"+str(i)+".tif"))
        aclo.SetParameterString("out", os.path.join(args.tmp, "up_mask_S_Z-"+str(i)+".tif"))
        aclo.SetParameterOutputImagePixelType("out", otb.ImagePixelType_uint8)
        aclo.SetParameterInt("channel", 1)
        aclo.SetParameterInt("xradius", 1)
        aclo.SetParameterInt("yradius", 1)
        aclo.SetParameterString("filter","erode")
        aclo.ExecuteAndWriteOutput()

        np_surf = aclo.GetImageAsNumpyArray('out')
        best_surf = np.count_nonzero((np_surf == [1]))

        surf_m2 = best_surf*pixelSizeX*pixelSizeY
        logging.info("@"+ str(i) +"m: "+ str(surf_m2)+" m2.")
        #  logging.debug("@" + str(i) + "m: best_surf: "+str(best_surf)+" pixels "+"[label: "+str(best_label)+"]")
        #  logging.debug("--> "+str(surf_m2)+" m2 ")
        logging.debug("--> "+str(surf_m2*0.000001)+" km2 ")

        zi_up.append(i)
        S_zi_up.append(surf_m2)

        prev_alt = i


    zi.append(args.zmin)
    S_zi.append(0.0)
    data = np.column_stack((zi, S_zi))
    np.savetxt(args.outfile, data)

    zi_up.append(args.zmin)
    S_zi_up.append(0.0)
    data_up = np.column_stack((zi_up, S_zi_up))
    np.savetxt(args.outfile+"_up.dat", data_up)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
