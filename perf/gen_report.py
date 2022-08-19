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


def main(arguments):
    """gen_report.py
    Entry point to generate performance report of the current revision
    """
    # CLI
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--sites",
        nargs="+",
        default=["perf/data/occitanie/", "perf/data/andalousie/"],
        help="Paths to mega-site direcories",
    )
    parser.add_argument(
        "--exec",
        default="run_processors.py",
        help="Paths run_processors.py script",
    )
    parser.add_argument("-o", "--out", help="Output directory", default="perf/reports/")
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
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

    # Check output directory:
    if pathlib.Path(args.out).is_dir() is False:
        logging.error("Output directory " + args.out + "does not exist.")
        raise RuntimeError("Output directory " + args.out + "does not exist.")

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
                    + args.out
                )
            )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
