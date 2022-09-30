#!/usr/bin/python
#  -*- coding: utf-8 -*-
"""
:author: Benjamin Tardy
:organization: CNES
:copyright: 2022 CNES. All rights reserved.
:license: see LICENSE file
:created: 2022
"""
import argparse
import logging
import sys

import geopandas as gpd

try:
    import pygeos
except:
    raise ModuleNotFoundError("You must install pygeos to use sjoin_nearest")

if __name__ == "__main__":
    """
    Usage : python generate_list_from_DB.py dams_file dam_id dam_name ouput_list
    """

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("poly_file", type=str, help="Water polygons")
    parser.add_argument("point_file", type=str, help="DAM information file")  # ID_SWOT
    parser.add_argument("output_file", type=str, help="Fused output file")

    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    print(args.poly_file)
    gdf_geom = gpd.GeoDataFrame().from_file(args.poly_file)
    # Il faut obligatoirement être dans une projection en metre pour utiliser
    # sjoin_nearest sinon les résultats ne sont pas fiables cf doc et warning
    gdf_points = gpd.GeoDataFrame().from_file(args.point_file).to_crs("EPSG:2154")
    out_gdf = gpd.sjoin_nearest(gdf_geom, gdf_points, how="left")
    out_gdf = out_gdf.to_crs("EPSG:4326")
    out_gdf.to_file(args.output_file)
