#!/usr/bin/python
#  -*- coding: utf-8 -*-
'''
:author: Gaël Nicolas
:organization: CS SI
:copyright: 2021 CS. All rights reserved.
:license: see LICENSE file
:created: 2021
'''

import os
import sys
import argparse
import logging
import time
import subprocess



def run_processing(cmd_list, stdoutfile, stderrfile, title="", nb_procs="1"):
    """ Run qsub processings.

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
    stdout_param = open(stdoutfile, 'a')
    stderr_param = open(stderrfile, 'a')

    logging.info("Running :  {} {}".format(title, cmd_list))

    while len(cmd_list) > 0 or len(pids) > 0:
        if (len(pids) < int(nb_procs)) and (len(cmd_list) > 0):
            pids.append([subprocess.Popen(cmd_list[0],
                                          stdout=stdout_param,
                                          stderr=stderr_param,
                                          shell=True),
                         cmd_list[0]])
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
                print(title + "... " + str(int((nb_cmd - (len(cmd_list) \
                                     + len(pids))) * 100. / nb_cmd)) + "%")
                time.sleep(0.2)
                break
    end = time.time()

    print(str(title + " done, elapsed time : " + str(end - start)))


def mk_dir(path):
    ''' Test existance and create a directorie.

     :param path: SW input directory
     :type path : str
    '''
    if os.path.exists(path):
        logging.warning("!! Directory {} already exist".format(path))
    else:
        os.makedirs(path, mode=0o755)



if __name__ == "__main__":
    """
    Usage : python run_processors.py dams_list dams_db dem_path wmap_path chain_dir out_dir
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('dams_list',type=str, help='Dams list')
    parser.add_argument('dams_db',type=str, help='Dams database path')
    parser.add_argument('dem_path',type=str, help='DEM path')
    parser.add_argument('wmap_path',type=str, help='surfwater map path')
    parser.add_argument('chain_dir',type=str, help='dem4water chain directory')
    parser.add_argument('out_dir',type=str, help='HSV directory')

    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    # Create output directories
    mk_dir(args.out_dir)
    log_dir = os.path.join(args.out_dir, 'log')
    mk_dir(log_dir)

    # Dams name and id
    dams_dict={}
    with open(args.dams_list, 'r') as file_name:
        dams_to_process = file_name.readlines()
    for dam in dams_to_process:
            dams_dict[dam.split(',')[0]] = dam.split(',')[1].rstrip()

    # ======================#
    # Processor compute_hsv #
    # ======================#

    for cle in dams_dict.keys():

        dame_name = dams_dict[cle].replace(" ", "_")
        cmd_compute_hsv = []
        cmd_compute_hsv.append(str("qsub -W umask=117"
                             + " -v WD=" + args.chain_dir
                             + ",DAM=" + dame_name
                             + ",DAM_ID="+ cle
                             + ",ID_FIELD="+ "ID_SWOT"
                             + ",DB_PATH="+ args.dams_db
                             + ",DEM_PATH="+ args.dem_path
                             + ",WMAP_PATH="+ args.wmap_path
                             + ",ROOT_DIR="+ args.out_dir
                             + " -l walltime=20:00:00"
                             + " -l select=1:ncpus=12:mem=60000MB:os=rh7"
                             + " -o " + os.path.join(log_dir, "compute_hsv_" + cle + "_out.log")
                             + " -e " + os.path.join(log_dir, "compute_hsv_" + cle + "_err.log")
                             + " compute_hsv.pbs"))

        run_processing(cmd_compute_hsv,
                       os.path.join(log_dir, "qsub_dem4water_out.log"),
                       os.path.join(log_dir, "qsub_dem4water_err.log"),
                       title = "qsub_dem4water")

