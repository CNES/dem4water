#!/usr/bin/env python3
"""
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CNES
:copyright: 2022 CNES. All rights reserved.
:license: see LICENSE file
:created: 2022
"""

import argparse
import csv
import json
import logging
import os
import pathlib
import subprocess
import sys
from datetime import datetime


def get_current_git_rev():
    return (
        subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
        .decode("ascii")
        .strip()
    )


def get_score(json, dam, measure):
    score = "☓"
    for record in json:
        if str(record["ID"]) == str(dam):
            if str(record[measure[0]][measure[1]]["mean"]) == "NaN":
                score = str(record[measure[0]][measure[1]]["mean"])
            else:
                score = "%.4f" % float(record[measure[0]][measure[1]]["mean"])

    return score


def get_error(json, dam, measure):
    score = "☓"
    for record in json:
        if str(record["ID"]) == str(dam):
            if str(record[measure[0]][measure[1]]) == "NaN":
                score = str(record[measure[0]][measure[1]])
            else:
                score = "%.4f" % float(record[measure[0]][measure[1]])

    return score


def get_trend(line):
    trend = "⊙"
    if line.split(" | ")[1] == "☓" or line.split(" | ")[1] == "NaN":
        trend = "☓"
    elif line.split(" | ")[2] == "☓" or line.split(" | ")[2] == "NaN":
        trend = "⊙"
    else:
        try:
            if abs(float(line.split(" | ")[1])) > abs(float(line.split(" | ")[2])):
                trend = '<span style="color:red">↘</span>'
            elif abs(float(line.split(" | ")[1])) < abs(float(line.split(" | ")[2])):
                trend = '<span style="color:green">↗</span>'
            elif abs(float(line.split(" | ")[1])) == abs(float(line.split(" | ")[2])):
                trend = '<span style="color:blue">↔</span>'
        except ValueError:
            logging.debug("N/A")
            trend = "☓"

    return get_trophy(line) + " " + trend


def get_trophy(line):
    trophy = '<span style="color:yellow">🏅</span>'
    if line.split(" | ")[1] == "☓" or line.split(" | ")[1] == "NaN":
        trophy = " "
    else:
        best = abs(float(line.split(" | ")[1]))
        for measure in line.split(" | ")[2:]:
            if measure == "☓" or measure == "NaN" or measure == "\n":
                continue
            else:
                # logging.debug("measure:<" + measure + ">")
                if abs(float(measure)) < best:
                    trophy = " "
                    break

    return trophy


def get_table(measure, sites, rep_files, func):
    # Report Table Title
    table = "\n## " + measure[0] + " " + measure[1] + " quality trend\n"
    # Report Table Header
    table += "| Dam ID | Dam Name | Trend | "
    sep_line = "|:------|:------|:------"
    for rep_file in rep_files:
        table += rep_file.stem + " | "
        sep_line += "|:------"

    table += "\n"
    sep_line += ":|\n"
    table += sep_line

    for site in sites:
        p_site = pathlib.Path(site)
        if p_site.is_dir() is False:
            logging.error("Baseline site: " + site + "Not Found.")
        else:
            with open(pathlib.Path(p_site, p_site.stem + ".lst"), "r") as f:
                reader = csv.reader(f)
                dams = {rows[0]: rows[1] for rows in reader}
            logging.debug(
                "Baseline site: " + p_site.stem + " [" + str(len(dams)) + " sites]"
            )

        for dam in dams:
            beg_line = "| " + dam + " | " + dams[dam] + " | "
            line = " | "
            for rep_file in rep_files:
                with open(rep_file, "r") as f:
                    records = json.load(f)
                line += func(records, dam, measure) + " | "
                # line += get_error(records, dam, error) + " | "

            line += "\n"
            trend = get_trend(line)
            table += beg_line + trend + line
        # logging.debug(table)

    return table


def run_campaign(args):
    """Run the baseline campaign."""
    logging.info("Starting baseline campaign execution.")
    # Check output directory:
    if pathlib.Path(args.outdir).is_dir() is False:
        logging.error("Output directory " + args.outdir + " does not exist.")
        raise RuntimeError("Output directory " + args.outdir + " does not exist.")

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
                    + site_cfg["wmap_path"]
                    + " "
                    + str(pathlib.Path(args.exec).absolute().parent)
                    + " "
                    + str(pathlib.Path(args.outdir, get_current_git_rev()).absolute())
                    + " --ref_model "
                    + str(pathlib.Path(p_site, p_site.stem + "_ref.json"))
                )
            )
        with open(
            pathlib.Path(args.outdir, get_current_git_rev(), "version.txt"), "w+"
        ) as vf:
            vf.write(get_current_git_rev())


def run_report(args):
    """Collect and compile performance metrics from an existing campaign."""
    logging.info("Starting campaign report generation.")
    try:
        with open(pathlib.Path(args.indir, "version.txt")) as vf:
            version = vf.readlines()[0]
    except FileNotFoundError:
        logging.error(
            "No version file, Directory does not seems to contain campaign results."
        )
        raise RuntimeError(
            "No version file, Directory "
            + args.indir
            + " does not seems to contain campaign results."
        )
    logging.info("Generating report for revision " + version + ".")

    rep_files = sorted(pathlib.Path(args.indir).glob("**/*_report.json"))
    logging.debug(rep_files)

    report = []

    for rep in rep_files:
        with open(rep, "rb") as infile:
            report.append(json.load(infile))

    logging.debug("\n" + json.dumps(report, indent=4, sort_keys=True))

    with open(
        pathlib.Path(
            args.outdir, datetime.now().strftime("%Y%m%d") + "_" + version + ".json"
        ),
        "w",
    ) as outfile:
        json.dump(report, outfile, indent=4)


def run_dashboard(args):
    """Generate a trend dashboard from an existing reports."""
    # ↗ Better
    # ↘ Worst
    # ↔ Same
    # ⊙ New
    # ☓ Missing
    # 🏆 Best

    measures = [
        ["S(z)_quality", "glob"],
        ["V(S)_quality", "glob"],
        ["Vr(S)_quality", "glob"],
    ]

    errors = [["Dam_bottom_estimation", "Dam_bottom_error"]]

    rep_files = sorted(pathlib.Path(args.indir).glob("**/*.json"), reverse=True)
    logging.debug(rep_files)
    rev = rep_files[0].stem.split("_")[1]
    page = "# Reference Campaign Dashboard for revision " + rev

    for measure in measures:
        page += get_table(measure, args.sites, rep_files, get_score)

    page += get_table(errors[0], args.sites, rep_files, get_error)

    if (
        pathlib.Path(args.outdir, datetime.now().strftime("%Y%m%d") + ".md").exists()
        is False
    ):
        outfile = pathlib.Path(args.outdir, datetime.now().strftime("%Y%m%d") + ".md")
    else:
        outfile = pathlib.Path(
            args.outdir, datetime.now().strftime("%Y%m%d_%H%M%S") + ".md"
        )

    with open(outfile, "w") as outfile:
        outfile.write(page)

    logging.info(
        "Up-to-date dashboard relative to revision "
        + rev
        + " can be found at "
        + str((outfile.name))
        + "."
    )


def run_full(args):
    """Run the whole benchmark end-to-end."""
    logging.error("End-to-end benchmark - NOT YET IMPLEMENTED.")
    # NOTE: IMHO running steps through CI,
    # including pushing resulting markdown to wiki would be better


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
        "campaign", help="1- Campaign mode, runs the model estimations."
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
    parser_camp.add_argument("-o", "--outdir", help="Output directory", required=True)

    # Report sub-command
    parser_rep = sub_parsers.add_parser(
        "report",
        help="2- Report mode, collect and compile validation report scores of an existing campaign.",
    )
    parser_rep.add_argument(
        "--sites",
        nargs="+",
        default=["perf/data/occitanie/", "perf/data/andalousie/"],
        help="Paths to mega-site direcories",
    )
    parser_rep.add_argument(
        "-i",
        "--indir",
        help="Directory containing the results of a campaign",
        required=True,
    )
    parser_rep.add_argument(
        "-o", "--outdir", help="Output directory", default="perf/reports/"
    )

    # Dashboard sub-command
    parser_dash = sub_parsers.add_parser(
        "dashboard",
        help="3- Dashboard mode, generate a comprehensive performance trend dashboard.",
    )
    parser_dash.add_argument(
        "--sites",
        nargs="+",
        default=["perf/data/occitanie/", "perf/data/andalousie/"],
        help="Paths to mega-site direcories",
    )
    parser_dash.add_argument(
        "-i",
        "--indir",
        help="Directory containing the reports",
        default="perf/reports/",
    )
    parser_dash.add_argument(
        "-o", "--outdir", help="Output directory", default="perf/dashboards/"
    )

    # Full sub-command
    parser_full = sub_parsers.add_parser(
        "end-2-end",
        help="End-2-end benchmark mode, run the whole process.",
    )
    parser_full.add_argument(
        "--sites",
        nargs="+",
        default=["perf/data/occitanie/", "perf/data/andalousie/"],
        help="Paths to mega-site directories",
    )
    parser_full.add_argument(
        "--exec",
        default="run_processors.py",
        help="Paths run_processors.py script",
    )
    parser_full.add_argument(
        "-o", "--outdir", help="Campaign output directory", required=True
    )
    parser_full.add_argument(
        "-r", "--repdir", help="Report directory", default="perf/reports/"
    )
    parser_full.add_argument(
        "-d", "--dashdir", help="Dashboard directory", default="perf/dashboards/"
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

    if args.mode == "campaign":
        # Expose revision
        logging.info("Currently using revision " + get_current_git_rev())
        run_campaign(args)
    elif args.mode == "report":
        run_report(args)
    elif args.mode == "dashboard":
        run_dashboard(args)
    elif args.mode == "end-2-end":
        # Expose revision
        logging.info("Currently using revision " + get_current_git_rev())
        run_full(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
