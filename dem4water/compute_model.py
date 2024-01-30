#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module providing tools for computing models and filter points."""
import json
import logging
import math
import sys
from statistics import median

import numpy as np

import dem4water.water_body as wb


def remove_jump_szi(infile, z_i=None, s_zi=None, jump_ratio=4):
    """Use a ratio to remove first jump found.

    PARAMETERS
    ---------
    args: argparse object
    z_i: optional list of altitude
    s_zi: optional list of surface
    jump_ratio: ratio between current and previous surface.
    """
    if z_i is None and s_zi is None:
        print(infile)
        data = np.loadtxt(infile)

        if data.size <= 2:
            logging.error(f"Not enough S(Zi) data inside file {infile}")
            sys.exit("Error")

        z_i = data[:, 0]
        s_zi = data[:, 1]

    # remove outliers / virtual surface overflow
    stop_i = 0
    break_found = False
    prev = s_zi[1]
    for s_z in s_zi:
        if s_z != 0:
            ratio = prev / s_z
        else:
            ratio = 1
        if ratio < jump_ratio:
            prev = s_z
            stop_i = stop_i + 1
        else:
            break_found = True
            break

    if break_found is True:
        logging.debug(
            f"Dropping S_ZI after index {stop_i} with a delta ratio of {ratio}."
        )
        z_i = z_i[stop_i:]
        s_zi = s_zi[stop_i:]
    else:
        logging.debug("No outliers detected, keeping all S_ZI data.")
    return z_i, s_zi


def filter_szi(
    infile,
    database,
    watermap,
    damname,
    max_elev,
    min_elev,
    water_map_thres=0.05,
    area_threshold=15,
):
    """Filter szi looking at the watermap."""
    data = np.loadtxt(infile)
    if data.size <= 2:
        logging.error(f"Not enough S(Zi) data inside file {infile}")
        sys.exit("Error")

    z_i = data[:, 0]
    s_zi = data[:, 1]
    shp_wmap = wb.create_water_mask(watermap, water_map_thres)
    # water_body_area = wb.compute_area_from_water_body(daminfo, shp_wmap)
    water_body_area = wb.compute_area_from_database_geom(database, damname, shp_wmap)
    logging.info(f"water body area: {water_body_area}")
    thres_wb = (water_body_area * area_threshold) / 100
    logging.info(f"Minimal surface allowed : {thres_wb}")
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
            logging.info(
                f"Szi {val_szi} for altitude {val_zi} is too small. Check the cutline."
            )
    filter_szi_out.append(s_zi[-1])
    filter_zi_out.append(z_i[-1])
    return filter_zi_out, filter_szi_out


def get_info_dam(daminfo):
    """Extract information from daminfo name."""
    with open(daminfo, encoding="utf-8") as info:
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


def select_lower_szi(z_i, sz_i, damelev, max_offset, winsize):
    """Select all valid point under the damelev to find law."""
    filtered_szi = [
        szi for zi_, szi in zip(z_i, sz_i) if zi_ < float(damelev) + float(max_offset)
    ]
    filtered_zi = [
        zi_ for zi_, szi in zip(z_i, sz_i) if zi_ < float(damelev) + float(max_offset)
    ]
    # Add 1 to winsize as the first point is Z0, and it must not be used as a valid point
    # to compute model
    if len(filtered_szi) < winsize + 1:
        if len(z_i) > winsize + 1:
            filtered_szi = sz_i[: winsize + 1]
            filtered_zi = z_i[: winsize + 1]
        else:
            filtered_szi = sz_i[:]
            filtered_zi = z_i[:]
    else:
        filtered_szi = filtered_szi[: winsize + 1]
        filtered_zi = filtered_zi[: winsize + 1]
    return filtered_zi, filtered_szi


def compute_model(z_i, s_zi, z_0, s_z0):
    """Compute model from points."""
    logging.info(f"Model computed using {len(s_zi)} values.")
    poly = np.polyfit(z_i, s_zi, 1, rcond=None, full=False, w=None, cov=False)
    beta = poly[0] * (median(z_i) - z_0) / (median(s_zi) - s_z0)
    print(poly[0], z_i, z_0, beta)
    alpha = poly[0] * (math.pow(median(z_i) - z_0, 1 - beta)) / beta

    mae_sum = 0
    for alt, s_z in zip(z_i, s_zi):
        surf = s_z0 + alpha * math.pow((alt - z_0), beta)
        mae_sum += abs(s_z - surf)
    mae = mae_sum / len(z_i)
    logging.info(f"MAE computed: {mae}")
    return alpha, beta, mae, poly
