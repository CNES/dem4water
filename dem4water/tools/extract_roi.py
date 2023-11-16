#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract ROI from an image."""
import argparse
import sys
from dataclasses import dataclass
from math import ceil, floor

import numpy as np
import rasterio as rio
from affine import Affine
from pyproj import Transformer
from rasterio.coords import BoundingBox
from rasterio.windows import Window

from dem4water.tools.save_raster import save_image

DTYPE = {
    "uint8": rio.uint8,
    "uint16": rio.uint16,
    "uint32": rio.uint32,
    "float": rio.float32,
}


@dataclass
class ExtractROIParam:
    """Class providing parameters for extract roi functions."""

    mode: str
    mode_radius_cx: int
    mode_radius_cy: int
    mode_radius_r: int = 0
    mode_radius_unitr: str = "phy"
    mode_radius_unitc: str = "phy"
    dtype: str = "float"


def coord_phys_to_pixel(
    in_raster: str,
    extractroi_parameters: ExtractROIParam,
):
    """
    Transform coordinates to index.

    Parameters
    ----------
    in_raster:str
    """
    x_pixel, y_pixel = in_raster.index(
        extractroi_parameters.mode_radius_cx, extractroi_parameters.mode_radius_cy
    )
    return x_pixel, y_pixel


def extract_roi(
    in_raster: rio.io.DatasetReader,
    extractroi_parameters: ExtractROIParam,
):
    """
    Extract ROI from a raster.

    Parameters
    ----------
    in_raster:rio.io.DatasetReader

    ExtractROI_parameters: ExtractROIParam,
    """
    if extractroi_parameters.mode == "radius":
        if (
            extractroi_parameters.mode_radius_unitr == "phy"
            and extractroi_parameters.mode_radius_unitc == "phy"
        ):
            resolution = in_raster.res[0]
            dist_radius = extractroi_parameters.mode_radius_r / resolution
            x_pixel, y_pixel = coord_phys_to_pixel(in_raster, extractroi_parameters)

            width = (dist_radius) * 2 + 1
            height = (dist_radius) * 2 + 1
            if (dist_radius % 2) == 0:
                col_off, row_off = (
                    x_pixel - dist_radius,
                    y_pixel - dist_radius,
                )

            else:
                col_off, row_off = (
                    ceil(x_pixel - dist_radius),
                    floor(y_pixel - dist_radius),
                )

            window = Window(row_off, col_off, width, height)
            data = in_raster.read(window=window)
            transform = rio.windows.transform(window, in_raster.transform)

            profile = in_raster.profile
            profile.update(
                {
                    "height": height,
                    "width": width,
                    "transform": transform,
                    "driver": "GTiff",
                }
            )
            return data, profile


def compute_roi_from_ref(ref_file, in_file, out_file):
    with rio.open(ref_file) as ref:
        with rio.open(in_file) as crop:
            # Define crop area from ref bounds
            bounds_ref = ref.bounds
            ll_ref = (bounds_ref.left, bounds_ref.bottom)
            lr_ref = (bounds_ref.right, bounds_ref.bottom)
            ul_ref = (bounds_ref.left, bounds_ref.top)
            ur_ref = (bounds_ref.right, bounds_ref.top)
            x_ref = [corner[0] for corner in [ul_ref, ur_ref, ll_ref, lr_ref]]
            y_ref = [corner[1] for corner in [ul_ref, ur_ref, ll_ref, lr_ref]]
            # Ensure the bounds are in correct projection
            transformer = Transformer.from_crs(ref.crs, crop.crs, always_xy=True)
            x_proj, y_proj = transformer.transform(x_ref, y_ref)

            bounds_crop = BoundingBox(
                left=np.min(x_proj),
                bottom=np.min(y_proj),
                right=np.max(x_proj),
                top=np.max(y_proj),
            )
            # Crop the data
            window_crop = rio.windows.from_bounds(
                bounds_crop.left,
                bounds_crop.bottom,
                bounds_crop.right,
                bounds_crop.top,
                crop.transform,
            )
            crop_array = crop.read(window=window_crop)
            # Manage profile to write result in file
            profile = crop.profile
            dst_transform = Affine(
                crop.res[0], 0.0, bounds_crop.left, 0.0, -crop.res[1], bounds_crop.top
            )
            profile.update(
                {
                    "transform": dst_transform,
                    "height": crop_array.shape[1],
                    "width": crop_array.shape[2],
                    "driver": "GTiff",
                }
            )
            save_image(crop_array, profile, out_file)


def main():
    """Define parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-in", "--imref", help="Image to be processed", required=True)
    parser.add_argument("-out", "--outim", help="Output Image", required=True)

    parser.add_argument(
        "-mode",
        "--mode",
        help="Extraction mode",
        choices=["standard", "fit", "extent", "radius"],
        required=True,
    )
    parser.add_argument(
        "-mode.radius.r", "--mode_radius_r", help="Radius", type=int, default=0
    )
    parser.add_argument(
        "-mode.radius.unitr",
        "--mode_radius_unitr",
        help="Radius unit",
        choices=["pxl", "phy"],
    )
    parser.add_argument(
        "-mode.radius.cx",
        "--mode_radius_cx",
        help="X coordinate of the center",
        type=int,
        default=0,
    )
    parser.add_argument(
        "-mode.radius.cy",
        "--mode_radius_cy",
        help="Y coordinate of the center",
        type=int,
        default=0,
    )
    parser.add_argument(
        "-mode.radius.unitc",
        "--mode_radius_unitc",
        help="Center unit",
        choices=["pxl", "phy"],
    )
    parser.add_argument(
        "-dtype",
        "--dtype",
        help="Data type of output raster",
        choices=["uint8", "uint16", "uint32", "float"],
        default="uint8",
    )
    args = parser.parse_args()
    params = ExtractROIParam(
        args.mode,
        args.mode_radius_r,
        args.mode_radius_unitr,
        args.mode_radius_cx,
        args.mode_radius_cy,
        args.mode_radius_unitc,
        DTYPE[args.dtype],
    )
    with rio.open(args.imref) as imref:
        data, profile = extract_roi(imref, params)
        save_image(data, profile, args.out)


if __name__ == "__main__":
    sys.exit(main())
