#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Superimpose."""
import argparse
import sys
from dataclasses import dataclass
from math import ceil
from typing import Optional, Union

import numpy as np
import rasterio as rio
from rasterio.warp import Resampling, reproject

from dem4water.tools.save_raster import save_image

DTYPE = {
    "uint8": rio.uint8,
    "uint16": rio.uint16,
    "uint32": rio.uint32,
    "float": rio.float32,
}


@dataclass
class SuperimposeParam:
    """Class providing parameters for superimpose function."""

    interpolator: str
    dtype: str


def get_interpolator(superimpose_parameters: SuperimposeParam):
    """Convert str name to rio object."""
    if superimpose_parameters.interpolator == "nn":
        resampling = Resampling.nearest
    if superimpose_parameters.interpolator == "bco":
        resampling = Resampling.cubic
    if superimpose_parameters.interpolator == "linear":
        resampling = Resampling.bilinear
    return resampling


def superimpose_ndarray(
    input_image,
    input_image_profile,
    input_ref_profile,
    superimpose_parameters: SuperimposeParam,
) -> rio.io.DatasetReader:
    """
    Poject an image into the geometry of another one.

    Parameters
    ----------
    input_image_path:str
    image_ref:str
    output_image_path:str
    """
    resampling = get_interpolator(superimpose_parameters)
    output_image= np.ones(np.shape(input_image), DTYPE[superimpose_parameters.dtype])
    result, _ = reproject(
        source=input_image,
        destination=output_image,
        src_transform=input_image_profile["transform"],
        src_crs=input_image_profile["crs"],
        dst_crs=input_ref_profile["crs"],
        dst_transform=input_ref_profile["transform"],
        resampling=resampling,
        src_nodata=input_image_profile["nodata"],
        dst_nodata=input_image_profile["nodata"]
    )
    height = ceil(input_ref_profile["height"])
    width = ceil(input_ref_profile["width"])

    result = result[:, 0:height, 0:width]
    profile = input_ref_profile
    profile.update(
        {
            "count": input_image_profile["count"],
            "dtype": DTYPE[superimpose_parameters.dtype],
            "driver": "GTiff",
            "nodata": input_image_profile["nodata"]
        }
    )
    return result, profile


def superimpose(
    input_image: Union[np.ndarray, rio.DatasetReader],
    image_ref: Union[np.ndarray, rio.DatasetReader],
    superimpose_parameters: SuperimposeParam,
    input_ref_profile: Optional[rio._base.DatasetBase] = None,
    input_image_profile: Optional[rio._base.DatasetBase] = None,
) -> rio.io.DatasetReader:
    """Reproject a image according a reference.

    Parameters
    ----------
    input_image
    """
    if isinstance(input_image, rio.DatasetReader):
        raster_input = input_image.read()
        raster_profile_input = input_image.profile
    elif isinstance(input_image, np.ndarray):
        raster_input = input_image
        raster_profile_input = input_image_profile
    else:
        raise ValueError(
            f"{type(input_image)} not handled. Only ndarray or rio.DatasetReader"
        )
    if isinstance(image_ref, rio.DatasetReader):
        raster_profile_ref = image_ref.profile
    elif isinstance(image_ref, np.ndarray):
        raster_profile_ref = input_ref_profile
    else:
        raise ValueError(
            f"{type(image_ref)} not handled. Only ndarray or rio.DatasetReader"
        )
    data, profile = superimpose_ndarray(
        raster_input,
        raster_profile_input,
        raster_profile_ref,
        superimpose_parameters,
    )

    return data, profile


def main():
    """Define parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-inr", "--imref", help="The input reference image", required=True
    )
    parser.add_argument(
        "-inm",
        "--im",
        help="The image to reproject into the geometry of the reference input.",
        required=True,
    )
    parser.add_argument("-out", "--outim", help="Output Image", required=True)
    parser.add_argument(
        "-interpolator",
        "--interpolator",
        help="Extraction mode",
        choices=["nn", "bco", "linear"],
        required=True,
    )
    parser.add_argument(
        "-dtype",
        "--dtype",
        help="Data type of output raster",
        choices=["uint8", "uint16", "uint32", "float"],
        default="uint8",
    )
    args = parser.parse_args()
    params = SuperimposeParam(
        args.interpolator,
        DTYPE[args.dtype],
    )
    with rio.open(args.im) as image:
        with rio.open(args.imref) as reference:
            data, profile = superimpose(image, reference, params)
            save_image(data, profile, args.outim)


if __name__ == "__main__":
    sys.exit(main())
