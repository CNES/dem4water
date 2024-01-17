#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Polygonize a raster file."""
import geopandas as gpd
import rasterio
from rasterio import features
from shapely.geometry import shape


def polygonize(in_raster, out_vector):
    with rasterio.open(in_raster) as src:
        array = src.read().astype(rasterio.float32)
        no_data = src.nodata
        mask = array != no_data

        shapes = features.shapes(array, mask=mask, transform=src.transform)
        values = []
        geometry = []
        for shapedict, value in shapes:
            if value != 0:
                values.append(value)
                geometry.append(shape(shapedict))
        gdf = gpd.GeoDataFrame(
            {"DN": values, "geometry": geometry}, crs=src.crs.to_epsg()
        )
        gdf.to_file(out_vector)
