"""Console script for dem4water."""
import argparse
import json
import logging
import os
import sys

from dem4water.area_mapping import area_mapping
from dem4water.cut_contourlines import cut_countourlines
from dem4water.cutline_score import cutline_score
from dem4water.find_pdb_and_cutline import find_pdb_and_cutline
from dem4water.szi_to_model import szi_to_model
from dem4water.val_report import val_report


def launch_full_process(input_config_json):
    """Console script for dem4water."""
    # parser = argparse.ArgumentParser()
    # parser.add_argument("_", nargs="*")
    # args = parser.parse_args()

    with open(input_config_json, encoding="utf-8") as in_config:
        config = json.load(in_config)

    extract_watermap = config["area_mapping"]["out_wmap"]
    extract_dem = config["area_mapping"]["out_dem"]

    if os.path.exists(extract_watermap) and os.path.exists(extract_dem):
        logging.info(
            f"{extract_watermap} already exists.\n"
            f"{extract_dem} already exists.\n"
            ".Skipping area_mapping"
        )
    else:
        area_mapping(**config["area_mapping"])

    find_pdb_and_cutline(**config["find_pdb_and_cutline"])
    cutline_score(**config["cutline_score"])
    cut_countourlines(**config["cut_contourlines"])
    szi_to_model(**config["szi_to_model"])
    if "val_report" in config:
        val_report(**config["val_report"])
    return 0


def process_parameters():
    """."""
    # CLI
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")

    # create sub-parser
    sub_parsers = parser.add_subparsers(
        title="Modes",
        description="Select a mode",
        dest="mode",
        required=True,
    )
    # mode campaign
    # 1 VRT DEM WMAP 1 BD
    parser_camp = sub_parsers.add_parser(
        "campaign", help="1- Campaign mode, runs the model estimations."
    )
    parser_camp.add_argument(
        "-json_campaign", help="Configuration file for a complete campaign"
    )

    # mode autovalidation
    # lancer les 40  fichiers json andalousie & occitanie

    # mode DAM unique
    # dem4water --json agly.json


def main():
    """."""


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
