#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract ROI from an image"""
import argparse
import sys
from dataclasses import dataclass
import rasterio as rio
from rasterio.warp import reproject, Resampling


DTYPE = {
    "uint8": rio.uint8,
    "uint16": rio.uint16,
    "uint32": rio.uint32,
    "float": rio.float32,
}


@dataclass
class SuperimposeParam:
    """Class providing parameters for superimpose function."""

    interpolator: str ="bco"
    dtype: str ="uint32"


def get_interpolator(Superimpose_parameters: SuperimposeParam):
    if Superimpose_parameters.interpolator == "nn":
        resampling = Resampling.nearest
    if Superimpose_parameters.interpolator == "bco":
        resampling = Resampling.cubic
    if Superimpose_parameters.interpolator == "linear":
        resampling = Resampling.bilinear
    return resampling


def superimpose_save_raster(raster_out:rio.io.DatasetReader, output_path:str, Superimpose_parameters: SuperimposeParam):
    with rio.open(
        output_path, 
        "w",
        driver="GTiff",
        height=raster_out.height,
        width=raster_out.width,
        count=raster_out.count,
        dtype=Superimpose_parameters.dtype,
        crs=raster_out.crs,
        transform=raster_out.transform,
        ) as dataset:
            dataset.write(raster_out.read())
    
def superimpose(
    input_image:rio.io.DatasetReader,
    image_ref:rio.io.DatasetReader,
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
    result, _ = reproject(
        source=input_image.read(),
        destination=input_image.read(),
        src_transform=input_image.transform,
        src_crs=input_image.crs,
        dst_crs=image_ref.crs,
        dst_transform=image_ref.transform,
        resampling=resampling,
        )
    result = result[:, 0 : image_ref.height, 0 : image_ref.width]
           
    with rio.MemoryFile() as memfile:
        with memfile.open(
            driver="GTiff",
            height=image_ref.height,
            width=image_ref.width,
            count=input_image.count,
            dtype=result.dtype,
            crs=image_ref.crs,
            transform=image_ref.transform,
        ) as dataset:
            dataset.write(result)
        dataset_reader = memfile.open()

        return dataset_reader


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
    superimpose(args.im, args.imref, args.outim, params)


if __name__ == "__main__":
    sys.exit(main())
