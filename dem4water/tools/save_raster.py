#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Save Image"""
import argparse
import sys

import rasterio as rio


def save_image(raster, profile, output_path):
    print(profile)
    with rio.open(output_path, "w", **profile) as dataset:
        dataset.write(raster)


def main():
    """Define parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-raster", "--raster", help="Raster to save", required=True)
    parser.add_argument("-profile", "--profile", help="Profile's raster", required=True)

    parser.add_argument(
        "-output_path",
        "--output_path",
        help="Image output path",
        required=True,
    )
    args = parser.parse_args()
    save_image(args.raster, args.profile, args.output_path)


if __name__ == "__main__":
    sys.exit(main())
