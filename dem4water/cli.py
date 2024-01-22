"""Console script for dem4water."""
import argparse
import json
import logging
import os
import subprocess
import sys

from dem4water.area_mapping_v2 import area_mapping
from dem4water.cut_contourlines import cut_countourlines
from dem4water.find_cutline_and_pdb import find_cutline_and_pdb
from dem4water.find_pdb_and_cutline import find_pdb_and_cutline
from dem4water.szi_to_model import szi_to_model
from dem4water.tools.generate_dam_json_config import write_json
from dem4water.val_report import val_report


def run_command(args):
    """"""
    output = subprocess.run(args, capture_output=True)
    print("###############")
    print("Return code:", output.returncode)
    # use decode function to convert to string
    print("Output:", output.stdout.decode("utf-8"))


def get_current_git_rev():
    """Get the current revision number from git sources."""
    git_folder = os.path.join(os.path.dirname(__file__), "..", ".git")

    return (
        subprocess.check_output(
            ["git", f"--git-dir={git_folder}", "rev-parse", "--short", "HEAD"]
        )
        .decode("ascii")
        .strip()
    )


def launch_pbs(conf, log_out, log_err, cpu=12, ram=60, h_wall=1, m_wall=0):
    """Submit a job to pbs."""
    # Export system variables simulating loading modules and venv
    pbs_file = (
        f"#!/usr/bin/env bash\n"
        f"#PBS -l select=1:ncpus={cpu}:mem={ram}000MB:os=rh7\n"
        f"#PBS -l walltime={int(h_wall):02d}:{int(m_wall):02d}:0\n\n"
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
        f"dem4water single -dam_json {conf} -scheduler_type local"
    )
    out_file = log_out.replace(".log", ".pbs")
    with open(out_file, "w", encoding="utf-8") as ofile:
        ofile.write(pbs_file)
    run_command(["qsub", out_file])


def launch_slurm(
    conf,
    log_out,
    log_err,
    cpu=12,
    ram=60,
    h_wall=1,
    m_wall=0,
    account="campus",
    dam_name=None,
):
    """Submit a job to pbs."""
    # Export system variables simulating loading modules and venv
    if dam_name is None:
        name = "dem4water"
    else:
        name = f"d4w_{dam_name}"
    pbs_file = (
        "#!/bin/bash\n"
        f"#SBATCH --job-name {name}"
        "#SBATCH -N 1\n"
        "#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task={cpu}\n"
        f"#SBATCH --mem={ram}gb\n"
        f"#SBATCH --time={int(h_wall):02d}:{int(m_wall):02d}:00\n"
        f"#SBATCH --error={log_err}\n"
        f"#SBATCH --output={log_out}\n"
        f"#SBATCH --account={account}\n"
        "\nmodule purge\n"
        f"export PYTHONPATH={os.environ.get('PYTHONPATH')}\n"
        f"export PATH={os.environ.get('PATH')}\n"
        f"export LD_LIBRARY_PATH={os.environ.get('LD_LIBRARY_PATH')}\n"
        "export OTB_APPLICATION_PATH="
        f"{os.environ.get('OTB_APPLICATION_PATH')}\n"
        f"export GDAL_DATA={os.environ.get('GDAL_DATA')}\n"
        f"export GEOTIFF_CSV={os.environ.get('GEOTIFF_CSV')}\n"
        f"export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={cpu}\n\n"
        f"dem4water single -dam_json {conf} -scheduler_type local"
    )
    out_file = log_out.replace(".log", ".slurm")
    with open(out_file, "w", encoding="utf-8") as ofile:
        ofile.write(pbs_file)
    print("Submit slurm job ", out_file)
    run_command(["sbatch", out_file])


def launch_campaign(
    json_campaign,
    scheduler,
    walltime_hour,
    walltime_minutes,
    ram,
    cpu,
    input_force_list,
):
    """Launch on PBS or local."""
    config_list = write_json(json_campaign, input_force_list)
    # config_list = [config_list[0]]
    if scheduler == "local":
        for conf in config_list:
            launch_full_process(conf)
    else:
        for conf in config_list:
            name = conf.split("/")[-1].split("_")[-1].split(".")[0]
            with open(conf, encoding="utf-8") as in_config:
                config = json.load(in_config)
                log_out = config["chain"]["log_out"]
                log_err = config["chain"]["log_err"]
                if scheduler == "PBS":
                    launch_pbs(
                        conf,
                        log_out,
                        log_err,
                        h_wall=walltime_hour,
                        m_wall=walltime_minutes,
                        ram=ram,
                        cpu=cpu,
                    )
                elif scheduler == "Slurm":
                    launch_slurm(
                        conf,
                        log_out,
                        log_err,
                        h_wall=walltime_hour,
                        m_wall=walltime_minutes,
                        ram=ram,
                        cpu=cpu,
                        dam_name=name,
                    )


def launch_reference_validation_campaign(
    targets,
    output_folder,
    campaign_name,
    scheduler,
    walltime_hour,
    walltime_minutes,
    ram,
    cpu,
    only_ref,
):
    """Launch andalousie or occitanie reference campaign."""
    config_list = []
    git_folder = os.path.dirname(__file__)
    json_campaign_andalousie = os.path.join(
        git_folder, "..", "data", "andalousie", "campaign_andalousie_params.json"
    )
    json_campaign_occitanie = os.path.join(
        git_folder, "..", "data", "occitanie", "campaign_occitanie_params.json"
    )
    if "andalousie" in targets:
        if not os.path.exists(json_campaign_andalousie):
            print(f"Trying to access {json_campaign_andalousie} but not found.")
        else:
            print(f"Using {json_campaign_andalousie}")
            config_list += write_json(
                json_campaign_andalousie,
                output_folder,
                campaign_name,
                ref_only=only_ref,
            )
    if "occitanie" in targets:
        if not os.path.exists(json_campaign_occitanie):
            print(f"Trying to access {json_campaign_occitanie} but not found.")
        else:
            print(f"Using {json_campaign_occitanie}")
            config_list += write_json(
                json_campaign_occitanie,
                output_folder,
                campaign_name,
                concat=True,
                ref_only=only_ref,
            )

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
                launch_pbs(
                    conf,
                    log_out,
                    log_err,
                    h_wall=walltime_hour,
                    m_wall=walltime_minutes,
                    ram=ram,
                    cpu=cpu,
                )


def launch_single(conf, scheduler, walltime_hour, walltime_minutes, ram, cpu):
    """Launch a single dam on PBS or local mode."""
    if scheduler == "local":
        launch_full_process(conf)
    else:
        name = conf.split("/")[-1].split("_")[-1].split(".")[0]
        with open(conf, encoding="utf-8") as in_config:
            config = json.load(in_config)
            log_out = config["chain"]["log_out"]
            log_err = config["chain"]["log_err"]
            launch_slurm(
                conf,
                log_out,
                log_err,
                h_wall=walltime_hour,
                m_wall=walltime_minutes,
                ram=ram,
                cpu=cpu,
                dam_name=name,
            )


def launch_full_process(input_config_json):
    """Console script for dem4water."""
    with open(input_config_json, encoding="utf-8") as in_config:
        config = json.load(in_config)

    # extract_watermap = config["area_mapping"]["out_wmap"]
    extract_dem = config["area_mapping"]["out_dem"]
    algo = config["area_mapping"]["mode"]
    # if (os.path.exists(extract_watermap) and
    if os.path.exists(extract_dem):
        logging.info(
            # f"{extract_watermap} already exists.\n"
            f"{extract_dem} already exists.\n"
            ".Skipping area_mapping"
        )
    else:
        area_mapping(**config["area_mapping"])
    if algo == "GDP":
        find_cutline_and_pdb(**config["find_cutline_and_pdb"])
    else:
        find_pdb_and_cutline(**config["find_pdb_and_cutline"])
    # cutline_score(**config["cutline_score"])
    cut_countourlines(**config["cut_contourlines"])
    szi_to_model(**config["szi_to_model"])
    if "val_report" in config:
        val_report(**config["val_report"])
    return 0


def process_parameters():
    """Define all parameters."""
    # CLI
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
    parser.add_argument("--walltime_hour", help="Job duration for PBS: hour", default=1)
    parser.add_argument(
        "--walltime_minutes", help="Job duration for PBS: minutes", default=0
    )
    parser.add_argument("--cpu", help="PBS number of CPU ressource", default=12)
    parser.add_argument("--ram", help="PBS ram value (in GB)", default=60)
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
        "-json_campaign",
        help="Configuration file for a complete campaign",
        required=True,
    )
    parser_camp.add_argument(
        "-scheduler_type",
        help="Local or PBS",
        default="Slurm",
        choices=["local", "Slurm"],
    )
    parser_camp.add_argument(
        "-input_force_list",
        help="Text file containing id_dam,dam_name to force the processing of these dams.",
        default=None,
    )
    # mode autovalidation
    # lancer les 40  fichiers json andalousie & occitanie
    parser_ref = sub_parsers.add_parser(
        "camp_ref",
        help="3- Launch the pre-designed Andalousie and /or Occitanie campaign.",
    )
    parser_ref.add_argument(
        "-name", help="name the output folder", default=get_current_git_rev()
    )
    parser_ref.add_argument(
        "-output_folder", help="path to store outputs", required=True
    )
    parser_ref.add_argument(
        "-sites",
        nargs="+",
        default=["occitanie", "andalousie"],
        help="Paths to mega-site direcories",
    )
    parser_ref.add_argument(
        "-only_ref", action="store_true", help="Launch only dams with a ref."
    )
    parser_ref.add_argument(
        "-scheduler_type",
        help="Local or PBS",
        default="Slurm",
        choices=["local", "PBS", "Slurm"],
    )
    # mode DAM unique
    # dem4water --json agly.json
    parser_single = sub_parsers.add_parser(
        "single", help="2- Single mode, process one dam"
    )
    parser_single.add_argument(
        "-dam_json", help="Configuration for an unique dam", required=True
    )
    parser_single.add_argument(
        "-scheduler_type",
        help="Local or PBS",
        default="Slurm",
        choices=["local", "PBS", "Slurm"],
    )

    return parser


def main():
    """."""
    parser = process_parameters()
    args = parser.parse_args()

    if args.mode == "campaign":
        launch_campaign(
            args.json_campaign,
            args.scheduler_type,
            args.walltime_hour,
            args.walltime_minutes,
            args.ram,
            args.cpu,
            args.input_force_list,
        )
    elif args.mode == "single":
        launch_single(
            args.dam_json,
            args.scheduler_type,
            args.walltime_hour,
            args.walltime_minutes,
            args.ram,
            args.cpu,
        )
    elif args.mode == "camp_ref":
        launch_reference_validation_campaign(
            args.sites,
            args.output_folder,
            args.name,
            args.scheduler_type,
            args.walltime_hour,
            args.walltime_minutes,
            args.ram,
            args.cpu,
            args.only_ref,
        )


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
