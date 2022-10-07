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
from time import perf_counter

def main(arguments):
    """cutline_score.py
    Qualify cutline quality by computing intersection with wmap
    """
    t1_start = perf_counter()
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--infile", help="Input json file")
    parser.add_argument("-w", "--watermap", help="Input water map file")
    parser.add_argument("-o", "--out", help="Output directory")
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
    args = parser.parse_args(arguments)

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
    logging.info("Starting cutline_score.py")

    # Rasterize cutline on wmap
    Rasterization = otb.Registry.CreateApplication("Rasterization")
    Rasterization.SetParameterString("in", args.infile)
    Rasterization.SetParameterString("im", args.watermap)
    Rasterization.SetParameterString("out", os.path.join(args.out, "cutline_mask.tif"))
    Rasterization.SetParameterFloat("background", 0)
    Rasterization.SetParameterString("mode", "binary")
    Rasterization.SetParameterFloat("mode.binary.foreground", 1)
    Rasterization.Execute()

    # Compute nb pixel in cutline
    np_mask = Rasterization.GetImageAsNumpyArray("out")
    np_mask_sum = np.sum(np_mask)
    logging.info("np_mask_sum= " + str(np_mask_sum))

    # Compute sum of occurence in wmap on cutline
    bm = otb.Registry.CreateApplication("BandMath")
    bm.SetParameterStringList("il", [args.watermap])
    bm.AddImageToParameterInputImageList(
        "il", Rasterization.GetParameterOutputImage("out")
    )
    bm.SetParameterString("out", os.path.join(args.out, "cutline_wmap.tif"))
    bm.SetParameterString("exp", "( (im2b1 == 1) ? im1b1 : 0 )")
    bm.Execute()

    np_wmask = bm.GetImageAsNumpyArray("out")
    np_wmask_sum = np.sum(np_wmask)
    logging.info("np_wmask_sum= " + str(np_wmask_sum))

    # Normalize
    score = np_wmask_sum / np_mask_sum

    # Output result
    logging.info("Score= " + str(score))

    if args.debug is True:
        Rasterization.ExecuteAndWriteOutput()
        bm.ExecuteAndWriteOutput()

    # Cleanup tmp files
    t1_stop = perf_counter()
    logger.info("Elapsed time:", t1_stop, 's', t1_start, 's')
 
    logger.info("Elapsed time during the whole program in s :",
       t1_stop-t1_start, 's')

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
