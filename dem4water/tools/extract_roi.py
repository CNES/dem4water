#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract ROI from an image"""
import argparse
import sys
from dataclasses import dataclass
import rasterio
from rasterio.windows import Window
from math import ceil, floor

DTYPE = {
    "uint8": rasterio.uint8,
    "uint16": rasterio.uint16,
    "uint32": rasterio.uint32,
    "float": rasterio.float32,
}


@dataclass
class ExtractROIParam:
    """Class providing parameters for extract roi functions."""

    mode: str
    mode_radius_r: int
    mode_radius_unitr: str
    mode_radius_cx: int
    mode_radius_cy: int
    mode_radius_unitc: str
    dtype: str


def coord_phys_to_pixel(
    in_raster: str,
    ExtractROI_parameters: ExtractROIParam,
):
    """
    Transform coordinates to index

    Parameters
    ----------

    in_raster:str
    """
    with rasterio.open(in_raster) as src:
        x_pixel, y_pixel = src.index(
            ExtractROI_parameters.mode_radius_cx, ExtractROI_parameters.mode_radius_cy
        )
        return x_pixel, y_pixel


def extract_roi(
    in_raster: str,
    out_raster: str,
    ExtractROI_parameters: ExtractROIParam,
):
    """
    Extract ROI from a raster

    Parameters
    ----------

    in_raster:str

    out_raster:str
    """
    if ExtractROI_parameters.mode == "radius":
        if (
            ExtractROI_parameters.mode_radius_unitr == "phy"
            and ExtractROI_parameters.mode_radius_unitc == "phy"
        ):
            with rasterio.open(in_raster) as src:
                resolution = src.res[0]
                dist_radius = ExtractROI_parameters.mode_radius_r / resolution

                x_pixel, y_pixel = coord_phys_to_pixel(in_raster, ExtractROI_parameters)

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
                data = src.read(window=window)
                transform = rasterio.windows.transform(window, src.transform)

                profile = src.profile
                profile.update(
                    {"height": height, "width": width, "transform": transform}
                )

                with rasterio.open(
                    out_raster, "w", dtypes=ExtractROI_parameters.dtype, **profile
                ) as dst:
                    dst.write(data)


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
    extract_roi(args.imref, args.outim, params)


if __name__ == "__main__":
    sys.exit(main())
