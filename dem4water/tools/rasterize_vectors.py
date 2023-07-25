#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rasterize a vector file according a reference raster."""
import argparse
import sys
from dataclasses import dataclass

import geopandas as gpd
import numpy as np
import rasterio
from rasterio import features
from rasterio.enums import MergeAlg

DTYPE = {
    "uint8": rasterio.uint8,
    "uint16": rasterio.uint16,
    "uint32": rasterio.uint32,
    "float": rasterio.float32,
}


@dataclass
class RasterizarionParams:
    """Class providing parameters for rasterization functions."""

    mode: str
    binary_foreground_value: int
    background_value: int
    column_field: str
    dtype: str


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


def save_rasterization(
    out_raster: str,
    raster: rasterio.io.DatasetReader,
    vector_rasterized: np.ndarray,
    dtype,
) -> None:
    """Save the result of rasterization into a raster.

    Parameters
    ----------
    out_raster: str
        Any file with geometry supported by geopandas
    raster: rasterio.io.DatasetReader
        The rasterio raster dataset, providing all georeferences
    vector_rasterized: np.ndarray
        The numpy array containing the rasterized vector data
    """
    with rasterio.open(
        out_raster,
        "w",
        driver="GTiff",
        crs=raster.crs,
        transform=raster.transform,
        dtype=dtype,
        count=1,
        width=raster.width,
        height=raster.height,
    ) as dst:
        dst.write(vector_rasterized, indexes=1)


def binary_rasterize(
    in_vector: str,
    in_raster: str,
    out_raster: str,
    rasterize_params: RasterizarionParams,
) -> None:
    """Rasterize a shapefile using binary value and a reference raster.

    Parameters
    ----------
    in_vector:

    in_raster:

    out_raster:

    binary_foreground_value:

    """
    raster = rasterio.open(in_raster)
    ref_proj = raster.crs.to_epsg()
    vector = prepare_vector(in_vector, ref_proj)
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
    save_rasterization(out_raster, raster, vector_rasterized, rasterize_params.dtype)


def attribute_rasterize(
    in_vector: str,
    in_raster: str,
    out_raster: str,
    rasterize_params: RasterizarionParams,
) -> None:
    """Rasterize a vector file according a field columns and a reference raster.

    Parameters
    ----------
    in_vector: str

    in_raster: str

    out_raster: str

    column_attribute: str

    """
    raster = rasterio.open(in_raster)
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
    save_rasterization(out_raster, raster, vector_rasterized, rasterize_params.dtype)


def rasterize(
    in_vector, in_raster, out_raster, rasterization_params: RasterizarionParams
) -> None:
    """Rasterize a vector file according a reference image.

    Parameters
    ----------
    in_vector:
        the input vector file
    in_raster:
        the reference image
    out_raster:
        the output file
    mode:
        the rasterization mode ('attribute' or 'binary')
    user_type:
        the type of out_raster image
    column_attribute:
        the field to rasterize in 'attribute' mode
    binary_foreground_value:
        the foreground value for 'binary' mode
    """
    if rasterization_params.mode == "attribute":
        attribute_rasterize(in_vector, in_raster, out_raster, rasterization_params)
    elif rasterization_params.mode == "binary":
        binary_rasterize(
            in_vector,
            in_raster,
            out_raster,
            rasterization_params,
        )
    else:
        raise ValueError(
            f"{rasterization_params.mode} is not supported."
            " Only 'attribute' or 'binary' are allowed"
        )


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
    rasterize(args.invector, args.imref, args.oraster, params)


if __name__ == "__main__":
    sys.exit(main())
