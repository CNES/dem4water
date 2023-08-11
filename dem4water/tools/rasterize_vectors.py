#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rasterize a vector file according a reference raster."""
import argparse
import sys
from dataclasses import dataclass
from typing import Optional

import geopandas as gpd
import rasterio as rio
from rasterio import features
from rasterio.enums import MergeAlg

from dem4water.tools.save_raster import save_image

DTYPE = {
    "uint8": rio.uint8,
    "uint16": rio.uint16,
    "uint32": rio.uint32,
    "float": rio.float32,
}


@dataclass
class RasterizarionParams:
    """Class providing parameters for rasterization functions."""

    mode: str
    binary_foreground_value: int
    background_value: int
    column_field: str = None
    dtype: str = "float"


def prepare_vector(in_vector: str, ref_proj: int) -> gpd.GeoDataFrame:
    """Open a vector file and change projection.

    Parameters
    ----------
    in_vector: str
        Any file with geometry supported by geopandas
    ref_proj: int
        The EPSG target code
    Returns
    -------
    vector: gpd.GeoDataFrame
        The reprojected geodataframe
    """
    # Read in vector

    vector = gpd.read_file(in_vector)
    vector = vector.to_crs(ref_proj)
    return vector


def binary_rasterize(
    in_vector: str,
    raster: rio.io.DatasetReader,
    rasterize_params: RasterizarionParams,
) -> None:
    """Rasterize a shapefile using binary value and a reference raster.

    Parameters
    ----------
    in_vector:

    in_raster:

    binary_foreground_value:

    """
    ref_proj = raster.crs.to_epsg()
    if isinstance(in_vector, str):
        vector = prepare_vector(in_vector, ref_proj)
    else:
        # If geodataframe
        vector = in_vector
    # Convert dataframe into a list for rasterize input format
    # no pair value means a binary rasterization
    geom = list(vector.geometry)
    # Rasterize vector using the shape and coordinate system of the raster

    vector_rasterized = features.rasterize(
        geom,
        out_shape=raster.shape,
        fill=rasterize_params.background_value,
        out=None,
        transform=raster.transform,
        all_touched=False,
        default_value=rasterize_params.binary_foreground_value,
        dtype=rasterize_params.dtype,
    )
    profile = raster.profile
    profile.update({"dtype": DTYPE[rasterize_params.dtype], "driver": "GTiff"})

    return vector_rasterized, profile


def attribute_rasterize(
    in_vector: str,
    raster: rio.io.DatasetReader,
    rasterize_params: RasterizarionParams,
) -> None:
    """Rasterize a vector file according a field columns and a reference raster.

    Parameters
    ----------
    in_vector: str
        the input vector file
    in_raster: str
        the reference image
    column_attribute: str

    """
    ref_proj = raster.crs.to_epsg()
    vector = prepare_vector(in_vector, ref_proj)
    # Convert dataframe into a list of pair (geometry, value)
    # where value is the attribute to burn
    geom_value = (
        (geom, int(value))
        for geom, value in zip(vector.geometry, vector[rasterize_params.column_field])
    )
    # Rasterize vector using the shape and transform of the raster
    vector_rasterized = features.rasterize(
        geom_value,
        out_shape=raster.shape,
        transform=raster.transform,
        all_touched=True,
        fill=rasterize_params.background_value,  # background value
        merge_alg=MergeAlg.replace,
        dtype=None,
    )
    profile = raster.profile
    profile.update({"dtype": DTYPE[rasterize_params.dtype], "driver": "GTiff"})

    return vector_rasterized, profile


def rasterize(
    in_vector: str,
    raster: rio.DatasetReader,
    rasterization_params: RasterizarionParams,
    out_raster: Optional[str] = None,
) -> None:
    """Rasterize a vector file according a reference image.

    Parameters
    ----------
    in_vector:
        the input vector file
    in_raster:
        the reference image
    rasterization_params: RasterizarionParams
        dataclass containing all parameters for rasterization
    out_raster:
        the output file

    """
    if rasterization_params.mode == "attribute":
        vector_rasterized, profile = attribute_rasterize(
            in_vector, raster, rasterization_params
        )
    elif rasterization_params.mode == "binary":
        vector_rasterized, profile = binary_rasterize(
            in_vector,
            raster,
            rasterization_params,
        )
    else:
        raise ValueError(
            f"{rasterization_params.mode} is not supported."
            " Only 'attribute' or 'binary' are allowed"
        )
    if out_raster is not None:
        save_image(vector_rasterized, profile, out_raster)
    return vector_rasterized, profile


def main():
    """Define parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-iv", "--invector", help="Input vector file", required=True)
    parser.add_argument("-imr", "--imref", help="Reference raster file", required=True)
    parser.add_argument("-out", "--oraster", help="Output raster file", required=True)

    parser.add_argument(
        "-mode",
        "--mode",
        help="Rasterization mode",
        choices=["attribute", "binary"],
        required=True,
    )

    parser.add_argument(
        "-attribute.field",
        "--field",
        help="Field name for attribute rasterization",
        default=None,
    )
    parser.add_argument(
        "-background", "--background", help="Background value (no data)", default=0
    )
    parser.add_argument(
        "-binary_foreground", "--foreground", help="Binary foreground value", default=1
    )
    parser.add_argument(
        "-dtype",
        "--dtype",
        help="Data type of output raster",
        choices=["uint8", "uint16", "uint32", "float"],
        default="uint8",
    )
    args = parser.parse_args()
    params = RasterizarionParams(
        args.mode, args.foreground, args.background, args.field, DTYPE[args.dtype]
    )
    with rio.open(args.imref) as imref:
        rasterize(args.invector, imref, params, args.oraster)


if __name__ == "__main__":
    sys.exit(main())
