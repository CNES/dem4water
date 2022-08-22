#!/usr/bin/env python3
"""
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CNES
:copyright: 2022 CNES. All rights reserved.
:license: see LICENSE file
:created: 2022
"""

import argparse
import json
import logging
import os
import pathlib
import subprocess
import sys


def get_current_git_rev():
    return (
        subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
        .decode("ascii")
        .strip()
    )


def run_campaign(args):
    """Run the baseline campaign."""
    # Check output directory:
    if pathlib.Path(args.out).is_dir() is False:
        logging.error("Output directory " + args.out + " does not exist.")
        raise RuntimeError("Output directory " + args.out + " does not exist.")

    # Retrieve baseline sites info
    for site in args.sites:
        p_site = pathlib.Path(site)
        if p_site.is_dir() is False:
            logging.error("Baseline site: " + site + "Not Found.")
        else:
            logging.info("Included baseline site: " + p_site.stem)
            with open(pathlib.Path(p_site, p_site.stem + ".cfg"), "r") as cfg:
                site_cfg = json.load(cfg)
            logging.debug("\n" + json.dumps(site_cfg, indent=4, sort_keys=True))

            # with open(pathlib.Path(p_site, p_site.stem + ".lst"), "r") as lst:
            #     dam_list = lst.readlines()
            # logging.debug(dam_list)

            os.system(
                str(
                    str(pathlib.Path(args.exec).absolute())
                    + " "
                    + str(pathlib.Path(p_site, p_site.stem + ".lst"))
                    + " "
                    + str(pathlib.Path(p_site, p_site.stem + ".geojson"))
                    + " "
                    + site_cfg["dem_path"]
                    + " "
                    + str(pathlib.Path(p_site, p_site.stem + "_ref.json"))
                    + " "
                    + site_cfg["wmap_path"]
                    + " "
                    + str(pathlib.Path(args.exec).absolute().parent)
                    + " "
                    + args.out
                )
            )


def run_report(args):
    """Collect and compile performance metrics from an existing campaign."""
    logging.error("NOT IMPLEMENTED.")


def run_dashboard(args):
    """Generate a trend dashboard from an existing reports."""
    logging.error("NOT IMPLEMENTED.")


def main(arguments):
    """gen_report.py
    Entry point to generate performance report of the current revision
    """
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

    # Campaign sub-command
    parser_camp = sub_parsers.add_parser(
        "campaign", help="Campaign mode, runs the model estimations."
    )
    parser_camp.add_argument(
        "--sites",
        nargs="+",
        default=["perf/data/occitanie/", "perf/data/andalousie/"],
        help="Paths to mega-site direcories",
    )
    parser_camp.add_argument(
        "--exec",
        default="run_processors.py",
        help="Paths run_processors.py script",
    )
    parser_camp.add_argument("-o", "--out", help="Output directory", required=True)

    # Report sub-command
    parser_rep = sub_parsers.add_parser(
        "report",
        help="Report mode, collect and compile validation report scores of an existing campaign.",
    )
    parser_rep.add_argument(
        "--sites",
        nargs="+",
        default=["perf/data/occitanie/", "perf/data/andalousie/"],
        help="Paths to mega-site direcories",
    )
    parser_rep.add_argument(
        "-i",
        "--in",
        help="Directory containing the results of a campaign",
        required=True,
    )
    parser_rep.add_argument(
        "-o", "--out", help="Output directory", default="perf/reports/"
    )

    # Dashboard sub-command
    parser_dash = sub_parsers.add_parser(
        "dashboard",
        help="Dashboard mode, generate a comprehensive performance trend dashboard.",
    )
    parser_dash.add_argument(
        "--sites",
        nargs="+",
        default=["perf/data/occitanie/", "perf/data/andalousie/"],
        help="Paths to mega-site direcories",
    )
    parser_dash.add_argument(
        "-i",
        "--in",
        help="Directory containing the reports",
        default="perf/reports/",
    )
    parser_dash.add_argument(
        "-o", "--out", help="Output directory", default="perf/dashboards/"
    )

    args = parser.parse_args(arguments)

    # Setup Logger
    logging_format = (
        "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
    )
    if args.debug is True:
        logging.basicConfig(
            stream=sys.stdout, level=logging.DEBUG, format=logging_format
        )
    else:
        logging.basicConfig(
            stream=sys.stdout, level=logging.INFO, format=logging_format
        )
    # Expose revision
    logging.info("Performance Estimation of revision " + get_current_git_rev())

    if args.mode == "campaign":
        run_campaign(args)
    elif args.mode == "report":
        run_report(args)
    elif args.mode == "dashboard":
        run_dashboard(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
