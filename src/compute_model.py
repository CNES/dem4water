#!/usr/bin/env python3
"""
:author: Benjamin Tardy <benjamin.tardy@csgroup.eu>
:organization: CS Group
:copyright: 2022 CNES. All rights reserved.
:license: see LICENSE file
:created: 2022
"""
# import argparse
import json
import logging
import math

# import os
import sys
from statistics import median

import numpy as np

import water_body as wb


def select_szi(args, Zi=None, S_Zi=None):
    if Zi is None and S_Zi is None:
        if args.custom_szi is not None:
            print("Dat file used : ", args.custom_szi)
            infile = args.custom_szi
        else:
            infile = args.infile
        data = np.loadtxt(infile)

        if data.size <= 2:
            logging.error(f"Not enought S(Zi) data inside file {infile}")
            sys.exit("Error")

        Zi = data[:, 0]
        S_Zi = data[:, 1]

    # remove outliers / virtual surface overflow
    stop_i = 0
    break_found = False
    prev = S_Zi[1]
    for z, sz in zip(Zi, S_Zi):
        if sz != 0:
            ratio = prev / sz
        else:
            ratio = 1
        #  print(str(ratio))
        # TODO: @parameters max ratio
        if ratio < 4:
            prev = sz
            stop_i = stop_i + 1
        else:
            break_found = True
            break

    if break_found is True:
        logging.debug(
            "Dropping S_ZI after index "
            + str(stop_i)
            + " with a delta ratio of "
            + str(ratio)
            + "."
        )
        Zi = Zi[stop_i:]
        S_Zi = S_Zi[stop_i:]
    else:
        logging.debug("No outliers detected, keeping all S_ZI data.")
    return Zi, S_Zi


def filter_szi(args, damname, max_elev, min_elev):
    if args.custom_szi is not None:
        print("Dat file used : ", args.custom_szi)
        infile = args.custom_szi
    else:
        infile = args.infile
    data = np.loadtxt(infile)
    if data.size <= 2:
        logging.error("Not enought S(Zi) data inside file " + infile)
        sys.exit("Error")

    z_i = data[:, 0]
    s_zi = data[:, 1]
    shp_wmap = wb.create_water_mask(args.watermap, 0.05)
    # water_body_area = wb.compute_area_from_water_body(args.daminfo, shp_wmap)
    water_body_area = wb.compute_area_from_database_geom(
        args.database, damname, shp_wmap
    )
    logging.info(f"water body area: {water_body_area}")
    thres_wb = (water_body_area * 15) / 100
    # zi[0] is the PDB
    zi_min = z_i[1]
    zi_max = z_i[-1]
    if zi_max > max_elev:
        logging.info("Too high contour detected filter S_ZI data")
    else:
        logging.info("Contour seems correct for max bound. Process")
    if zi_min < min_elev:
        logging.info("Too low contour detected filter S_ZI data")
    else:
        logging.info("Contour seems correct for min bound. Process")

    filter_zi_out = []
    filter_szi_out = []
    for val_zi, val_szi in zip(z_i[:-1], s_zi[:-1]):
        if val_szi > thres_wb:
            filter_zi_out.append(val_zi)
            filter_szi_out.append(val_szi)
        else:
            print(
                f"Szi {val_szi} for altitude {val_zi} is too small. Check the cutline."
            )
    filter_szi_out.append(s_zi[-1])
    filter_zi_out.append(z_i[-1])
    return filter_zi_out, filter_szi_out


def get_info_dam(daminfo):
    """Extract information from daminfo name."""
    with open(daminfo) as info:
        jsi = json.load(info)
        for feature in jsi["features"]:
            if feature["properties"]["name"] == "Dam":
                logging.debug(feature)
                damname = feature["properties"]["damname"]
                damelev = feature["properties"]["elev"]
                dam_id = feature["properties"]["ID"]
                return damname, damelev, dam_id
        raise ValueError(f"Field Dam not found in {daminfo}")
    raise OSError(f"Unable to open {daminfo}")


def select_lower_szi(z_i, sz_i, damelev, max_offset):
    """Select all valid point under the damelev to find law."""
    filtered_szi = [szi for zi_, szi in zip(z_i, sz_i) if zi_ < damelev + max_offset]
    filtered_zi = [zi_ for zi_, szi in zip(z_i, sz_i) if zi_ < damelev + max_offset]
    if len(filtered_szi) < 11:
        if len(z_i) > 11:
            filtered_szi = sz_i[:11]
            filtered_zi = z_i[:11]
        else:
            filtered_szi = sz_i[:]
            filtered_zi = z_i[:]
    return filtered_zi, filtered_szi


def compute_model(z_i, s_zi, z_0, s_z0):
    """Compute model from points."""
    logging.info(f"Model computed using {len(s_zi)} values.")
    poly = np.polyfit(z_i[1:], s_zi[1:], 1, rcond=None, full=False, w=None, cov=False)
    beta = poly[0] * (median([z_i]) - z_0) / (median(s_zi) - s_z0)
    alpha = poly[0] * (math.pow(median(z_i) - z_0, 1 - beta)) / beta

    mae_sum = 0
    for z, sz in zip(z_i, s_zi):
        s = s_zi[0] + alpha * math.pow((z - z_0), beta)
        mae_sum += abs(sz - s)
    mae = mae_sum / len(z_i)
    logging.info(f"MAE computed: {mae}")
    return alpha, beta, mae, poly


# def main(arguments):  # noqa: C901  #FIXME: Function is too complex
#     """
#     Prototype scrip allowing to derive a HSV model from a set of S(Z_i) values.

#     The first value of the set should be S_0 = S(Z_0) = 0 with Z_0 the altitude of dam bottom

#     The output
#     - the list of parameters [alpha;beta]
#     - a quality measurment of how well the model fit the S(Z_i) values
#     - additionnaly the plot of S(Z) and V(S)
#     """
#     parser = argparse.ArgumentParser(
#         description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
#     )

#     parser.add_argument("-i", "--infile", help="Input file")
#     parser.add_argument("-d", "--daminfo", help="daminfo.json file")
#     parser.add_argument("-wa", "--watermap", help="Water map product")
#     parser.add_argument(
#         "-db", "--database", help="The database geojson file with geometry"
#     )
#     parser.add_argument(
#         "-w", "--winsize", type=int, default=11, help="S(Zi) used for model estimation."
#     )
#     parser.add_argument(
#         "-z",
#         "--zmaxoffset",
#         type=int,
#         default=30,
#         help="Elevation offset on top of dam elevation used for ending optimal model search",
#     )
#     parser.add_argument(
#         "--zminoffset",
#         type=int,
#         default=10,
#         help="Elevation offset from dam elevation used for starting optimal model search",
#     )
#     list_of_mode = ["absolute", "first", "hybrid"]
#     parser.add_argument("-m", "--maemode", default="absolute", choices=list_of_mode)
#     parser.add_argument(
#         "--dslopethresh",
#         default=1000,
#         help="Threshold used to identify slope breaks in hybrid model selection",
#     )
#     parser.add_argument(
#         "--custom_szi", type=str, default=None, help="Custom SZi.dat file"
#     )
#     parser.add_argument("-o", "--outfile", help="Output file")
#     parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")

#     args = parser.parse_args(arguments)
#     logging_format = (
#         "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
#     )
#     if args.debug is True:
#         logging.basicConfig(
#             stream=sys.stdout, level=logging.DEBUG, format=logging_format
#         )
#     else:
#         logging.basicConfig(
#             stream=sys.stdout, level=logging.INFO, format=logging_format
#         )
#     logging.info("Starting szi_to_model.py")

#     damname, damelev, dam_id = get_info_dam(args.daminfo)
#     z_i, s_zi = filter_szi(args, damname, 100000, 0)
#     z_i, s_zi = select_szi(args, z_i, s_zi)
#     logging.debug(f"Number of S_Zi used for compute model: {len(s_zi)}")
#     z_i = z_i[::-1]
#     s_zi = s_zi[::-1]

#     logging.debug("Zi: ")
#     logging.debug(z_i[:])
#     logging.debug("S_Zi: ")
#     logging.debug(s_zi[:])
#     logging.debug("Zi_max: " + str(z_i[-1]) + " - S(Zi_max): " + str(s_zi[-1]))
#     logging.info("Z0: " + str(z_i[0]) + " - S(Z0): " + str(s_zi[0]))

#     select_z_i, select_s_zi = select_lower_szi(z_i, s_zi, damelev, max_offset)
#     alpha, beta = compute_model(select_z_i, select_s_zi)
#     model_json = {
#         "ID": dam_id,
#         "Name": damname,
#         "Elevation": damelev,
#         "Model": {
#             "Z0": z_i[0],
#             "S0": s_zi[0],
#             "V0": 0.0,
#             "alpha": alpha,
#             "beta": beta,
#         },
#     }

#     with open(os.path.splitext(args.outfile)[0] + ".json", "w") as write_file:
#         json.dump(model_json, write_file)
