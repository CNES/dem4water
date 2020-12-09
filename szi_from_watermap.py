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
import argparse
import numpy as np
import otbApplication as otb

def main(arguments):
    '''Entrypoint'''

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-w',
                        '--watermap',
                        help="Input water map file")
    parser.add_argument('-d',
                        '--dem',
                        help="Input DEM")
    parser.add_argument('-z',
                        '--zmax',
                        type=int,
                        default=185,
                        help="Starting altitude")
    parser.add_argument('-s',
                        '--step',
                        type=int,
                        default=5,
                        help="Altitude sampling step")
    parser.add_argument('-t',
                        '--tmp',
                        help="Temporary directory")
    parser.add_argument('-o',
                        '--outfile',
                        help="Output file",
                        default=sys.stdout,
                        type=argparse.FileType('w'))

    args = parser.parse_args(arguments)
    #  print(args)

    last_result = -1
    msk = otb.Registry.CreateApplication("BandMath")
    seg = otb.Registry.CreateApplication("Segmentation")
    for i in range(args.zmax, 0, -1 * args.step):

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

        #  np_surf = msk.GetImageAsNumpyArray('out')
        #  np_surf = seg.GetImageAsNumpyArray('mode.raster.out')
        np_surf = bm.GetImageAsNumpyArray('out')
        max_label = np.amax(np_surf)
        #  print("max_label: "+str(max_label))
        #  print("max_label: "+str(int(max_label)))
        #  min_label = np.amin(np_surf)
        best_surf = 0
        best_label = 0
        for l in range(1, int(max_label+1)):
            results_surf = np.count_nonzero((np_surf == [l]))
            if (results_surf > best_surf):
                #  print("new best_surf: "+str(results_surf)+" [label: "+str(l)+"]")
                best_label = l
                best_surf = results_surf


        print("@" + str(i) + "m: best_surf: "+str(best_surf)+" [label: "+str(best_label)+"]")

        if ( best_surf < last_result/2.0 ):
            print("Estimated surface dropped --> Stop here." )
            break
        else:
            last_result = best_surf



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
