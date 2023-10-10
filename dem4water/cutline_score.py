#!/usr/bin/env python3
"""This module provides function to compute a validation score for cutlines."""


import argparse
import logging
import os
import sys
from time import perf_counter

import numpy as np


def cutline_score(infile, watermap, out, debug=False):
    """Qualify cutline quality by computing intersection with wmap."""
    t1_start = perf_counter()

    logging_format = (
        "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
    )
    if debug is True:
        logging.basicConfig(
            stream=sys.stdout, level=logging.DEBUG, format=logging_format
        )
    else:
        logging.basicConfig(
            stream=sys.stdout, level=logging.INFO, format=logging_format
        )
    logging.info("Starting cutline_score.py")

    # Rasterize cutline on wmap
    rasterization = otb.Registry.CreateApplication("Rasterization")
    rasterization.SetParameterString("in", infile)
    rasterization.SetParameterString("im", watermap)
    rasterization.SetParameterString("out", os.path.join(out, "cutline_mask.tif"))
    rasterization.SetParameterFloat("background", 0)
    rasterization.SetParameterString("mode", "binary")
    rasterization.SetParameterFloat("mode.binary.foreground", 1)
    rasterization.Execute()

    # Compute nb pixel in cutline
    np_mask = rasterization.GetImageAsNumpyArray("out")
    np_mask_sum = np.sum(np_mask)
    logging.info(f"np_mask_sum= {np_mask_sum}")

    # Compute sum of occurence in wmap on cutline
    bandmath = otb.Registry.CreateApplication("BandMath")
    bandmath.SetParameterStringList("il", [watermap])
    bandmath.AddImageToParameterInputImageList(
        "il", rasterization.GetParameterOutputImage("out")
    )
    bandmath.SetParameterString("out", os.path.join(out, "cutline_wmap.tif"))
    bandmath.SetParameterString("exp", "( (im2b1 == 1) ? im1b1 : 0 )")
    bandmath.Execute()

    np_wmask = bandmath.GetImageAsNumpyArray("out")
    np_wmask_sum = np.sum(np_wmask)
    logging.info(f"np_wmask_sum= {np_wmask_sum}")

    # Normalize
    score = np_wmask_sum / np_mask_sum

    # Output result
    logging.info(f"Score= {score}")

    if debug is True:
        rasterization.ExecuteAndWriteOutput()
        bandmath.ExecuteAndWriteOutput()

    # Cleanup tmp files
    t1_stop = perf_counter()
    logging.info(f"Elapsed time: {t1_stop}s {t1_start}s")

    logging.info(f"Elapsed time during the whole program in s : {t1_stop-t1_start}s")


def cutline_score_parameters():
    """Define cutline_score.py parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--infile", help="Input json file")
    parser.add_argument("-w", "--watermap", help="Input water map file")
    parser.add_argument("-o", "--out", help="Output directory")
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
    return parser


def main():
    """Cli function for cutline_score."""
    parser = cutline_score_parameters()
    args = parser.parse_args()
    cutline_score(args.infile, args.watermap, args.out, args.debug)


if __name__ == "__main__":
    sys.exit(main())
