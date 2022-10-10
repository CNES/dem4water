#!/usr/bin/python
#  -*- coding: utf-8 -*-
"""
:author: Gaël Nicolas
:organization: CNES
:copyright: 2021 CNES. All rights reserved.
:license: see LICENSE file
:created: 2021
"""

import argparse
import json
import logging
import sys

if __name__ == "__main__":
    """
    Usage : python generate_list_from_DB.py dams_file dam_id dam_name ouput_list
    """

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "dams_file", type=str, help="Dams database json file"
    )  # /work/OT/siaa/Work/SCO_StockWater/dams_database/Andalousie/andalousie_selected_max_extent_v7.geojson
    parser.add_argument("dam_id", type=str, help="DAM id field")  # ID_SWOT
    parser.add_argument("dam_name", type=str, help="DAM name field")  # DAM_NAME
    parser.add_argument("ouput_list", type=str, help="Dams list")

    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    # Read anw write id and dams names
    with open(args.dams_file, "r") as json_file:
        reader = json.load(json_file)
        crs = reader["crs"]["properties"]["name"].split(":")[-1]
        if not crs == "CRS84":
            raise ValueError(f"CRS must be in WGS84/EPSG:4326, {crs} found.")
        with open(args.ouput_list, "w") as output_file:

            for feature in reader["features"]:
                output_file.write(
                    "%s,%s \n"
                    % (
                        int(feature["properties"][args.dam_id]),
                        feature["properties"][args.dam_name],
                    )
                )
