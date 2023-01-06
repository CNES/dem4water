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
from dem4water.tools.generate_dam_json_config import write_json
from dem4water.val_report import val_report


def launch_pbs(conf, log_out, log_err, cpu=12, ram=60, h_wall=1, m_wall=0):
    """Submit a job to pbs."""
    # Export system variables simulating loading modules and venv
    pbs_file = (
        f"#!/usr/bin/env bash\n"
        f"#PBS -l select=1:ncpus={cpu}:mem={ram}000MB:os=rh7\n"
        f"#PBS -l walltime={h_wall:02d}:{m_wall:02d}:0\n\n"
        f"#PBS -e {log_err}\n"
        f"#PBS -o {log_out}\n"
        "\nmodule purge\n"
        f"export PYTHONPATH={os.environ.get('PYTHONPATH')}\n"
        f"export PATH={os.environ.get('PATH')}\n"
        f"export LD_LIBRARY_PATH={os.environ.get('LD_LIBRARY_PATH')}\n"
        f"export OTB_APPLICATION_PATH={os.environ.get('OTB_APPLICATION_PATH')}\n"
        f"export GDAL_DATA={os.environ.get('GDAL_DATA')}\n"
        f"export GEOTIFF_CSV={os.environ.get('GEOTIFF_CSV')}\n"
        f"export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={cpu}\n\n"
        f"dem4water single -dam_json {conf}"
    )
    out_file = log_out.replace(".log", ".pbs")
    with open(out_file, "w", encoding="utf-8") as ofile:
        ofile.write(pbs_file)
    os.system(f"qsub {out_file}")


def launch_campaign(json_campaign, scheduler):
    """Launch on PBS or local."""
    config_list = write_json(json_campaign)
    # config_list = [config_list[0]]
    if scheduler == "local":
        for conf in config_list:
            launch_full_process(conf)
    else:
        for conf in config_list:
            with open(conf, encoding="utf-8") as in_config:
                config = json.load(in_config)
                log_out = config["chain"]["log_out"]
                log_err = config["chain"]["log_err"]
                launch_pbs(conf, log_out, log_err)


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
        "-json_campaign", help="Configuration file for a complete campaign", required=True
    )
    parser_camp.add_argument("-scheduler", help="Local or PBS", default="PBS")
    # mode autovalidation
    # lancer les 40  fichiers json andalousie & occitanie

    # mode DAM unique
    # dem4water --json agly.json
    parser_single = sub_parsers.add_parser(
        "single", help="2- Single mode, process one dam"
    )
    parser_single.add_argument("-dam_json", help="Configuration for an unique dam", required=True)
    parser_single.add_argument("-scheduler", help="Local or PBS", default="PBS")

    return parser


def main():
    """."""
    parser = process_parameters()
    args = parser.parse_args()

    if args.mode == "campaign":
        launch_campaign(args.json_campaign, args.scheduler)
    elif args.mode == "single":
        launch_full_process(args.dam_json)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
