#!/usr/bin/python3
#  -*- coding: utf-8 -*-
"""
:author: Gaël Nicolas.

:organization: CS SI
:copyright: 2021 CS. All rights reserved.
:license: see LICENSE file
:created: 2021
"""

import argparse
import glob
import logging
import os
import shutil
import subprocess
import time
from datetime import datetime

# import sys


def format_walltime(wall_h, wall_m):
    wall_h = int(wall_h)
    wall_m = int(wall_m)
    if wall_h < 0:
        raise ValueError(
            f"Hours for walltime cannot be negative : {wall_h} provided as parameter"
        )
    if wall_m < 0 or wall_m > 59:
        raise ValueError(f"Minutes must be between 0 and 59: {wall_m} provided")
    return f"{wall_h:02d}:{wall_m:02d}:00"


def save_previous_run(path, dam_name):
    """Copy old file to time named folder."""
    # find model.json
    path_model = os.path.join(path, "camp", dam_name, f"{dam_name}_model.json")
    file_of_interest = [
        "*_daminfo.json",
        "*_cutline.json",
        "*.png",
        "*model.json",
        "tmp/*.png",
        f"../../log/*{dam_name}*.log",
    ]

    if not os.path.exists(path_model):
        logging.info(f"Model {path_model} not found. No saving done.")
    else:
        # get the last modification time
        ti_m = os.path.getmtime(path_model)
        m_ti = time.ctime(ti_m)
        t_obj = time.strptime(m_ti)
        t_stamp = time.strftime("%Y%m%dT%H%M%S", t_obj)
        save_path = os.path.join(path, "camp", dam_name, f"{dam_name}_{t_stamp}")
        mk_dir(save_path)
        for pred in file_of_interest:
            files = glob.glob(os.path.join(path, "camp", dam_name, pred))
            for infile in files:
                shutil.copy(infile, save_path)
                logging.info(f"{infile} copied to {save_path}")


def find_corrected_input(path, dam_name, opt_path=None):
    """Look for manually corrected files."""
    dam_info_old = []
    dam_info_new = []
    cutline_old = []
    cutline_new = []
    additionnal_params = ""
    if opt_path is not None:

        dam_info_new = glob.glob(os.path.join(opt_path, f"{dam_name}_daminfo*.json"))
        cutline_new = glob.glob(os.path.join(opt_path, f"{dam_name}_cutline*.*json"))
        if len(dam_info_new) > 1:
            raise ValueError(
                f"More than one file found for dam info for {dam_name} dam."
            )
        if len(cutline_new) > 1:
            raise ValueError(
                f"More than one file found for cutline for {dam_name} dam."
            )

    else:
        dam_info_old = glob.glob(
            os.path.join(path, "camp", dam_name, f"{dam_name}_daminfo_custom.json")
        )
        # Look for json or geojson files
        cutline_old = glob.glob(
            os.path.join(path, "extracts", dam_name, f"{dam_name}_cutline_custom.*json")
        )
        if len(dam_info_old) > 1:
            raise ValueError(
                f"More than one file found for dam info for {dam_name} dam."
            )
        if len(cutline_old) > 1:
            raise ValueError(
                f"More than one file found for cutline for {dam_name} dam."
            )
    if dam_info_new:
        in_daminfo = os.path.join(
            path, "camp", dam_name, f"{dam_name}_daminfo_custom.json"
        )
        os.makedirs(os.path.dirname(in_daminfo), exist_ok=True)
        shutil.copy(dam_info_new[0], in_daminfo)
        logging.info(f"{dam_info_new[0]} copied to {in_daminfo}")
        additionnal_params += f",INFO_DAM={in_daminfo}"
    elif dam_info_old:
        additionnal_params += f",INFO_DAM={dam_info_old[0]}"
    else:
        logging.info(f"{dam_name} No custom daminfo file found.")

    if cutline_new:
        ext = cutline_new[0].split(".")[-1]
        in_cut = os.path.join(
            path, "extracts", dam_name, f"{dam_name}_cutline_custom.{ext}"
        )
        os.makedirs(os.path.dirname(in_cut), exist_ok=True)
        shutil.copy(cutline_new[0], in_cut)
        logging.info(f"{cutline_new[0]} copied to {in_cut}")

        additionnal_params += f",CUTLINE={in_cut}"
    elif cutline_old:
        additionnal_params += f",CUTLINE={cutline_old[0]}"
    else:
        logging.info(f"{dam_name} No custom cutline found")
    return additionnal_params


def run_processing(cmd_list, stdoutfile, stderrfile, title="", nb_procs="1"):
    """Run qsub processings.

    :param cmd_list : the commands to run
    :type cmd_list : str
    :param stdoutfile : name of stdoutfile
    :type stdoutfile : str
    :param stderrfile : name of stderrfile
    :type stderrfile : str
    :param title : optional title
    :type title : str
    :param nb_procs :number of processors
    :type nb_procs : str
    """
    nb_cmd = len(cmd_list)
    start = time.time()
    pids = []
    stdout_param = open(stdoutfile, "a")
    stderr_param = open(stderrfile, "a")

    logging.info(f"Running :  {title} {cmd_list}")

    while len(cmd_list) > 0 or len(pids) > 0:
        if (len(pids) < int(nb_procs)) and (len(cmd_list) > 0):
            pids.append(
                [
                    subprocess.Popen(
                        cmd_list[0],
                        stdout=stdout_param,
                        stderr=stderr_param,
                        shell=True,
                    ),
                    cmd_list[0],
                ]
            )
            cmd_list.remove(cmd_list[0])

        for i, pid in enumerate(pids):
            status = pid[0].poll()

            if status is not None and status != 0:
                print("!! ERROR in pid #" + str(i) + " id=" + str(pid[0]))
                print(pid[1])
                del pids[i]
                break
            if status == 0:
                del pids[i]
                print(
                    title
                    + "... "
                    + str(int((nb_cmd - (len(cmd_list) + len(pids))) * 100.0 / nb_cmd))
                    + "%"
                )
                time.sleep(0.2)
                break
    end = time.time()

    print(str(title + " done, elapsed time : " + str(end - start)))


def mk_dir(path):
    """Test existance and create a directorie.

    :param path: SW input directory
    :type path : str
    """
    if os.path.exists(path):
        logging.warning(f"!! Directory {path} already exist")
    else:
        os.makedirs(path, mode=0o755)


if __name__ == "__main__":
    """
    Usage : python run_processors.py dams_list dams_db dem_path ref_model wmap_path chain_dir out_dir
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("dams_list", type=str, help="Dams list")
    parser.add_argument("dams_db", type=str, help="Dams database path")
    parser.add_argument("dem_path", type=str, help="DEM path")
    parser.add_argument("wmap_path", type=str, help="surfwater map path")
    parser.add_argument("chain_dir", type=str, help="dem4water chain directory")
    parser.add_argument("out_dir", type=str, help="HSV directory")
    parser.add_argument(
        "--ref_model", type=str, help="Reference model path", default=None
    )
    parser.add_argument("--id_field", type=str, help="DAM ID column", default="ID_SWOT")
    parser.add_argument("--radius", type=str, help="PDB radius search", default=None)
    parser.add_argument(
        "--elev_off", type=int, help="Offset added to dam elevation", default=60
    )
    parser.add_argument("--walltime_h", type=int, help="Walltime hours", default=30)
    parser.add_argument("--walltime_m", type=int, help="Walltime minutes", default=0)
    parser.add_argument("--ncpu", type=int, help="Number of cpu", default=12)
    parser.add_argument("--ram", type=int, help="Ram in Mb", default=60000)

    parser.add_argument(
        "--correct_folder",
        type=str,
        help="Path to folder containing user edited files",
        default=None,
    )
    args = parser.parse_args()

    log_dir = os.path.join(args.out_dir, "log")
    logging_format = (
        "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
    )

    # Create output directories
    mk_dir(args.out_dir)
    mk_dir(log_dir)
    logging.basicConfig(
        level=logging.INFO,
        filename=os.path.join(log_dir, "run_processors.log"),
        filemode="a",
        format=logging_format,
    )
    walltime = format_walltime(args.walltime_h, args.walltime_m)
    # Dams name and id
    dams_dict = {}
    with open(args.dams_list, "r", encoding="utf-8") as file_name:
        dams_to_process = file_name.readlines()
    for dam in dams_to_process:
        dams_dict[dam.split(",")[0]] = dam.split(",")[1].rstrip()

    # ======================#

    # Processor compute_hsv #
    # ======================#
    current_time = datetime.now()
    date_time = current_time.strftime("%Y%m%dT%Hh%Mm")
    all_cmd = []
    for cle in dams_dict.keys():

        dame_name = dams_dict[cle].replace(" ", "-")

        # search for custom files
        add_params = find_corrected_input(args.out_dir, dame_name, args.correct_folder)
        if add_params == "":
            if os.path.exists(
                os.path.join(args.out_dir, "camp", dame_name, f"{dame_name}_model.json")
            ):
                logging.info(f"!! {dame_name} already processed. Skip !")
                continue
        else:
            # if exists save previous results
            save_previous_run(args.out_dir, dame_name)

        if args.ref_model is not None:
            add_params += f",REF_MODEL={args.ref_model}"
        if args.radius is not None:
            add_params += f",RADIUS={args.radius}"
        if args.elev_off is not None:
            add_params += f",ELEV_OFF_DAM={args.elev_off}"

        cmd_compute_hsv = []
        cmd_compute_hsv.append(
            str(
                "qsub -W umask=117"
                + " -v WD="
                + args.chain_dir
                + ",DAM="
                + dame_name
                + ",DAM_ID="
                + cle
                + ",ID_FIELD="
                + args.id_field
                + ",DB_PATH="
                + args.dams_db
                + ",DEM_PATH="
                + args.dem_path
                + ",WMAP_PATH="
                + args.wmap_path
                + ",ROOT_DIR="
                + args.out_dir
                + add_params
                + f" -l walltime={walltime}"
                + f" -l select=1:ncpus={args.ncpu}:mem={args.ram}MB:os=rh7"
                + " -o "
                + os.path.join(log_dir, f"{dame_name}_{cle}_out.log")
                + " -e "
                + os.path.join(log_dir, f"{dame_name}_{cle}_err.log")
                + " compute_hsv.pbs"
            )
        )
        all_cmd += cmd_compute_hsv

        run_processing(
            cmd_compute_hsv,
            os.path.join(log_dir, "qsub_dem4water_out.log"),
            os.path.join(log_dir, "qsub_dem4water_err.log"),
            title="qsub_dem4water",
        )
    with open(
        os.path.join(args.out_dir, f"command_list_{date_time}.txt"),
        "a",
        encoding="utf-8",
    ) as out_file:
        out_file.write("\n".join(all_cmd))
