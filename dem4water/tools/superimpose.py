#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Superimpose"""
import argparse
import sys
from dataclasses import dataclass
import rasterio as rio
from rasterio.warp import reproject, Resampling
import numpy as np
from math import ceil, floor
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


def get_interpolator(Superimpose_parameters: SuperimposeParam):
    if Superimpose_parameters.interpolator == "nn":
        resampling = Resampling.nearest
    if Superimpose_parameters.interpolator == "bco":
        resampling = Resampling.cubic
    if Superimpose_parameters.interpolator == "linear":
        resampling = Resampling.bilinear
    return resampling


def superimpose_ndarray(
    input_image,
    input_image_profile,
    image_ref,
    input_ref_profile,
    Superimpose_parameters: SuperimposeParam,
)-> rio.io.DatasetReader:
    """
    Poject an image into the geometry of another one

    Parameters
    ----------

    input_image_path:str
    image_ref:str
    output_image_path:str
    """
   
    resampling = get_interpolator(Superimpose_parameters)
    output_image=input_image.astype(Superimpose_parameters.dtype)
    result, _ = reproject(
        source=input_image,
        destination=output_image,
        src_transform=input_image_profile['transform'],
        src_crs=input_image_profile['crs'],
        dst_crs=input_ref_profile['crs'],
        dst_transform=input_ref_profile['transform'],
        resampling=resampling,
        )
    h= ceil(input_ref_profile['height'])
    w= ceil(input_ref_profile['width'])
    

    result = result[:, 0 : h, 0 :w]
    profile = input_ref_profile
    profile.update(
        {"count": input_image_profile['count'], "dtype":DTYPE[Superimpose_parameters.dtype],  "driver" :"GTiff"}
            )
    return result, profile

def superimpose(input_image, image_ref, Superimpose_parameters: SuperimposeParam,  input_ref_profile=None, input_image_profile =None)-> rio.io.DatasetReader:
    
    if not isinstance(input_image, np.ndarray):

        raster_input =input_image.read()
        raster_profile_input=input_image.profile
    else :
        raster_input =input_image
        raster_profile_input=input_image_profile

    if not isinstance(image_ref, np.ndarray):
        raster_ref =image_ref.read()
        raster_profile_ref=image_ref.profile
    else :
        raster_ref =image_ref
        raster_profile_ref=input_ref_profile
    data, profile =superimpose_ndarray(
        raster_input,
        raster_profile_input,
        raster_ref,
        raster_profile_ref,
        Superimpose_parameters
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
    superimpose(rio.open(args.im), rio.open(args.imref), args.outim, params)


if __name__ == "__main__":
    sys.exit(main())
