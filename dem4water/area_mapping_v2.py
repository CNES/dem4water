#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""This module define the study area from inputs."""
import argparse
import glob
import logging
import os
import sys
from time import perf_counter
from typing import Optional

import geopandas as gpd
import numpy as np
import rasterio as rio
from affine import Affine
from bmi_topography import Topography
from osgeo import gdal
from pyproj import Transformer
from rasterio import windows
from rasterio.coords import BoundingBox

from dem4water.tools.save_raster import save_image

logger = logging.getLogger("area_mapping_v2")
log = logging.getLogger()
log.setLevel(logging.ERROR)
logging.getLogger("geopandas").setLevel(logging.WARNING)

logging.getLogger("rasterio").setLevel(logging.WARNING)

logging.getLogger("fiona").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


def download_cop30(
    min_x: float, max_x: float, min_y: float, max_y: float, output_path: str
) -> str:
    """Download the Copernicus 30 meters DEM."""
    logger.info("Default DEM mode - automatic download: COPERNICUS-30m")

    params = {
        "dem_type": "COP30",
        "south": min_y,
        "north": max_y,
        "west": min_x,
        "east": max_x,
        "output_format": "GTiff",
        "cache_dir": output_path,
    }
    dem = Topography(**params)
    dem.fetch()
    dem_file = glob.glob(os.path.join(output_path, "COP30*"))[0]
    return dem_file


def extract_from_vrt(
    to_crop_file: Optional[str],
    min_x: float,
    max_x: float,
    min_y: float,
    max_y: float,
    t_epsg: int,
    out_file: str,
) -> str:
    """Extract a raster according a bounding box from a raster or vrt file a change proj.

    Parameters
    ----------
    to_crop_file:
        The file to crop
    min_x:
        Min x of the bounding box
    max_x:
        Max x of the bounding box
    min_y:
        Min y of the bounding box
    max_y:
        Max y of the bounding box
    t_epsg:
        The target projection
    out_file:
        The output file name
    """
    with rio.open(to_crop_file) as crop:
        transformer = Transformer.from_crs(t_epsg, crop.crs, always_xy=True)
        x_proj, y_proj = transformer.transform([min_x, max_x], [min_y, max_y])
        bounds_crop = BoundingBox(
            left=np.min(x_proj),
            bottom=np.min(y_proj),
            right=np.max(x_proj),
            top=np.max(y_proj),
        )

        # Crop the data
        window_crop = windows.from_bounds(
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
                "nodata": crop.profile["nodata"],
                "transform": dst_transform,
                "height": crop_array.shape[1],
                "width": crop_array.shape[2],
                "driver": "GTiff",
            }
        )

        save_image(crop_array, profile, out_file.replace(".tif", "_proj_ori.tif"))
        return out_file.replace(".tif", "_proj_ori.tif")


def get_dam(
    gdf_dam_db: gpd.GeoDataFrame,
    dam_name: Optional[str] = None,
    dam_name_col: Optional[str] = None,
    dam_id: Optional[int] = None,
    dam_id_col: Optional[str] = None,
) -> gpd.GeoDataFrame:
    if dam_name is not None:
        return gdf_dam_db[gdf_dam_db[dam_name_col] == dam_name]
    elif dam_id is not None:
        return gdf_dam_db[gdf_dam_db[dam_id_col] == dam_id]
    else:
        raise ValueError(f"{dam_name} or {dam_id} must be provided")


def area_mapping(
    dam_name: str,
    dam_id: int,
    dam_database: str,
    out_dir: str,
    dem: str,
    retrieve_mode: str = "local",
    epsg: int = 2154,
    target_resolution: int = 5,
    dam_name_col: str = "DAM_NAME",
    dam_id_col: Optional[str] = "ID_DB",
    buffer_roi: int = 1000,
    debug: bool = False,
) -> None:
    """Extract the DEM centered on the water body according a buffer.

    Parameters
    ----------
    dam_database:
        the geojson base containing dams
    dam_name:
        the dam to be processed
    retrieve_mode:
        choose between local or online download
    epsg:
        the working projection, must be meter compatible
    out_dir:
        the folder to store extracts
    target_resolution:float
        the resolution needed for DEM
    dem:
        if local, provide the DEM file path
    dam_name_col:
        the columns in dam_database containing the name field
    buffer_roi:
        the size to kept around the reservoir. Must be in meter and consistent with epsg
    dam_id:
        the dam id
    dam_id_col:
        the column name for id in dataframe
    debug:
        enable debug mode for more logs
    """
    t1_start = perf_counter()
    # Silence VRT related error (bad magic number)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    gdal.PushErrorHandler("CPLQuietErrorHandler")
    logging_format = (
        "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
    )
    if debug is True:
        logging.basicConfig(
            stream=sys.stdout, level=logging.DEBUG, format=logging_format
        )
    else:
        logging.basicConfig(
            stream=sys.stdout, level=logging.INFO, format=logging_format
        )
    logger.info(f"Starting area_mapping.py for {dam_name}")
    gdf_dam = gpd.read_file(dam_database)
    gdf_dam = gdf_dam.to_crs(epsg)
    gdf_dam = get_dam(gdf_dam, dam_name, dam_name_col, dam_id, dam_id_col)
    gdf_dam.to_file(os.path.join(out_dir, f"DB_{dam_name.replace(' ','-')}.geojson"))

    gdf_dam.geometry = gdf_dam.geometry.buffer(buffer_roi)

    min_x, min_y, max_x, max_y = gdf_dam.total_bounds
    out_file = os.path.join(out_dir, f"dem_extract_{dam_name.replace(' ', '-')}.tif")
    if retrieve_mode == "local":
        dem_roi = extract_from_vrt(dem, min_x, max_x, min_y, max_y, epsg, out_file)
    elif retrieve_mode == "cop30":
        gdf_dam_wgs84 = gdf_dam.to_crs(4326)
        min_x, min_y, max_x, max_y = gdf_dam_wgs84.total_bounds
        dem_roi = download_cop30(min_x, max_x, min_y, max_y, out_dir)
    else:
        raise ValueError(f"{retrieve_mode} is not a valid choice ")
    gdal.Warp(
        out_file,
        dem_roi,
        dstSRS=f"EPSG:{epsg}",
        xRes=target_resolution,
        yRes=-target_resolution,
        resampleAlg=gdal.GRA_CubicSpline,
        dstNodata=-10000,
    )
    t1_stop = perf_counter()
    logger.info(f"Elapsed time: {t1_stop} s {t1_start} s")

    logger.info(f"Elapsed time during the whole program in s : {t1_stop-t1_start} s")


def area_mapping_args() -> argparse.ArgumentParser:
    """Define area mapping parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--dam_database", help="Input file")
    parser.add_argument(
        "--retrieve_mode",
        help="The retrieve mode",
        choices=["local", "cop30"],
        default="local",
    )
    parser.add_argument("--dam_id", help="Dam ID")
    parser.add_argument(
        "--dam_id_col", help="Dam id field in database", default="ID_DB"
    )
    parser.add_argument("--dam_name", help="Dam name")
    parser.add_argument("--dam_name_col", help="Dam column field", default="DAM_NAME")
    parser.add_argument("-d", "--dem", help="Input DEM")
    parser.add_argument("--out_dir", help="Folder to store extracted dem")
    parser.add_argument(
        "--buffer_roi", help="Buffer area around the water body", default=1000
    )
    parser.add_argument("--epsg", help="The target projection as integer", default=2154)
    parser.add_argument(
        "--target_resolution",
        help="The target resolution (supposed to be squared",
        default=10,
    )
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
    return parser


def main() -> None:
    """Cli function to launch area mapping."""
    parser = area_mapping_args()
    args = parser.parse_args()
    area_mapping(
        dam_database=args.dam_database,
        dam_name=args.dam_name,
        retrieve_mode=args.retrieve_mode,
        epsg=args.epsg,
        out_dir=args.out_dir,
        target_resolution=args.target_resolution,
        dem=args.dem,
        dam_name_col=args.dam_name_col,
        buffer_roi=args.buffer_roi,
        dam_id=args.dam_id,
        dam_id_col=args.dam_id_col,
        debug=args.debug,
    )


if __name__ == "__main__":
    sys.exit(main())
