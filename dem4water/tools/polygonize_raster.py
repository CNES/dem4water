#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Polygonize a raster file."""
import os
from dataclasses import dataclass

import rasterio
from osgeo import gdal, ogr, osr


@dataclass
class PolygonizeParams:
    """Class providing parameters for rasterization functions."""

    layer_name: str
    field_name: str
    driver: str
    overwrite: bool


def polygonize(in_raster: str, out_vector: str, params_poly: PolygonizeParams):
    """Polygonize a raster in the raster projection.

    Parameters
    ----------
    in_raster: str
        the input raster file to polygonize
    out_vector: str
        the vector file to create
    params_poly: PolygonizeParams
        the polygonize parameters
    """
    #  get raster datasource
    src_ds = gdal.Open(in_raster)
    srcband = src_ds.GetRasterBand(1)
    dst_layername = params_poly.layer_name
    drv = ogr.GetDriverByName(params_poly.driver)
    if os.path.exists(out_vector):
        if params_poly.overwrite:
            print(f"{out_vector} was overwrited")
            os.remove(out_vector)
        else:
            return None
    dst_ds = drv.CreateDataSource(out_vector)
    with rasterio.open(in_raster) as src:
        epsg = src.crs.to_epsg()

    sp_ref = osr.SpatialReference()
    sp_ref.SetFromUserInput(f"EPSG:{epsg}")

    dst_layer = dst_ds.CreateLayer(dst_layername, srs=sp_ref)

    fld = ogr.FieldDefn(params_poly.field_name, ogr.OFTInteger)
    dst_layer.CreateField(fld)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex(params_poly.field_name)

    gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)
    del src_ds
    del dst_ds
    return None
