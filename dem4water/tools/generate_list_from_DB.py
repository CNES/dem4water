#!/usr/bin/python
#  -*- coding: utf-8 -*-
"""This tools parse the database file and provide a list of DAM ready to process."""

import argparse
import json
import logging
import sys


def create_dam_list_from_db(dams_file, dam_id, dam_name, output_list):
    """Parse geojson file and return the dam valid for processing."""
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    count_correct = 0
    count_incorrect = 0
    dam_to_process = {}
    # Read anw write id and dams names
    with open(dams_file, "r", encoding="utf-8") as json_file:
        reader = json.load(json_file)
        crs = reader["crs"]["properties"]["name"].split(":")[-1]
        if not crs == "CRS84":
            raise ValueError(f"CRS must be in WGS84/EPSG:4326, {crs} found.")
        with open(output_list, "w", encoding="utf-8") as output_file:

            for feature in reader["features"]:
                # DAM_LVL_M is hard coded in dem4water
                dam_lvl_m = feature["properties"]["DAM_LVL_M"]
                name_dam = feature["properties"][dam_name]
                if dam_lvl_m is None:
                    print(
                        f"{name_dam} has a empty value for DAM_LVL_M. "
                        "Impossible to process this dam."
                    )
                    count_incorrect += 1
                    continue
                try:
                    float(dam_lvl_m)
                    id_dam = int(feature["properties"][dam_id])

                    output_file.write(f"{id_dam},{name_dam} \n")
                    dam_to_process[name_dam] = id_dam
                    count_correct += 1
                except ValueError:
                    print(
                        f"DAM_LVL_M for {name_dam} is not a float compatible value. "
                        f"{dam_lvl_m} detected. Please check your DB"
                    )
                    count_incorrect += 1

    print(f"{count_correct + count_incorrect} dams found in {dams_file}")
    print(f"Ready to launch : {count_correct} dams")
    print(f"{count_incorrect} dams need corrections")
    return dam_to_process


if __name__ == "__main__":
    desc = (
        "Usage : python generate_list_from_DB.py dams_file dam_id dam_name ouput_list."
    )

    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("dams_file", type=str, help="Dams database json file")
    parser.add_argument("dam_id", type=str, help="DAM id field")  # ID_SWOT
    parser.add_argument("dam_name", type=str, help="DAM name field")  # DAM_NAME
    parser.add_argument("output_list", type=str, help="Dams list")

    args = parser.parse_args()
    create_dam_list_from_db(
        args.dams_file, args.dam_id, args.dam_name, args.output_list
    )
