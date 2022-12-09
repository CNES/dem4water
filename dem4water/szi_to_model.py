#!/usr/bin/env python3
"""
From Szi compute model.

:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CS Group
:copyright: 2020 CNES. All rights reserved.
:license: see LICENSE file
:created: 2020
"""
import argparse
import json
import logging
import math
import os
import sys
from statistics import median
from time import perf_counter

# import matplotlib.gridspec as gridspec  # noqa: F401
# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker
import numpy as np


def select_szi(args, z_i=None, s_zi=None):
    """Apply a ratio filtering on SZi.

    If Zi and S_Zi are None, the dat file is loaded.
    """
    if z_i is None and s_zi is None:
        if args.custom_szi is not None:
            print("Dat file used : ", args.custom_szi)
            infile = args.custom_szi
        else:
            infile = args.infile
        data = np.loadtxt(infile)

        if data.size <= 2:
            logging.error(f"Not enought S(Zi) data inside file {infile}")
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
        #  print(str(ratio))
        # TODO: @parameters max ratio
        if ratio < 4:
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


def filter_szi(args, damname, max_elev, min_elev):
    """Use area criterion to filter szi."""
    if args.custom_szi is not None:
        print("Dat file used : ", args.custom_szi)
        infile = args.custom_szi
    else:
        infile = args.infile
    data = np.loadtxt(infile)
    if data.size <= 2:
        logging.error(f"Not enought S(Zi) data inside file {infile}")
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


def main(arguments):  # noqa: C901  #FIXME: Function is too complex
    """
    Prototype scrip allowing to derive a HSV model from a set of S(Z_i) values.

    The first value of the set should be S_0 = S(Z_0) = 0 with Z_0 the altitude of dam bottom

    The output
    - the list of parameters [alpha;beta]
    - a quality measurment of how well the model fit the S(Z_i) values
    - additionnaly the plot of S(Z) and V(S)
    """
    t1_start = perf_counter()
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-i", "--infile", help="Input file")
    parser.add_argument("-d", "--daminfo", help="daminfo.json file")
    parser.add_argument("-wa", "--watermap", help="Water map product")
    parser.add_argument(
        "-db", "--database", help="The database geojson file with geometry"
    )
    parser.add_argument(
        "-w", "--winsize", type=int, default=11, help="S(Zi) used for model estimation."
    )
    parser.add_argument(
        "-z",
        "--zmaxoffset",
        type=int,
        default=30,
        help="Elevation offset on top of dam elevation used for ending optimal model search",
    )
    parser.add_argument(
        "--zminoffset",
        type=int,
        default=10,
        help="Elevation offset from dam elevation used for starting optimal model search",
    )
    list_of_mode = ["absolute", "first", "hybrid"]
    parser.add_argument("-m", "--maemode", default="absolute", choices=list_of_mode)
    parser.add_argument(
        "--dslopethresh",
        default=1000,
        help="Threshold used to identify slope breaks in hybrid model selection",
    )
    parser.add_argument(
        "--custom_szi", type=str, default=None, help="Custom SZi.dat file"
    )
    parser.add_argument("-o", "--outfile", help="Output file")
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")

    # Silence Mathplotlib related debug messages (font matching)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)

    args = parser.parse_args(arguments)

    #  print(args)

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
    logging.info("Starting szi_to_model.py")

    # load GeoJSON file containing info
    with open(args.daminfo) as i:
        jsi = json.load(i)
    for feature in jsi["features"]:
        if feature["properties"]["name"] == "Dam":
            logging.debug(feature)
            damname = feature["properties"]["damname"]
            damelev = feature["properties"]["elev"]
            dam_id = feature["properties"]["ID"]
        # if feature["properties"]["name"] == "PDB":
        #     pdbelev = feature["properties"]["elev"]
    # shp_wmap = wb.create_water_mask(args.watermap, 0.05)
    # water_body_area = wb.compute_area_from_database_geom(args.database, damname, shp_wmap)
    # wm_thres = (water_body_area * 15)/100
    z_i, s_zi = filter_szi(args, damname, 100000, 0)
    z_i, s_zi = select_szi(args, z_i, s_zi)
    logging.debug(f"Number of S_Zi used for compute model: {len(s_zi)}")
    z_i = z_i[::-1]
    s_zi = s_zi[::-1]
    print("Filtered : ", z_i, s_zi)
    logging.debug("z_i: ")
    logging.debug(z_i[:])
    logging.debug("s_zi: ")
    logging.debug(s_zi[:])
    logging.debug("Zi_max: " + str(z_i[-1]) + " - S(Zi_max): " + str(s_zi[-1]))
    logging.info("Z0: " + str(z_i[0]) + " - S(Z0): " + str(s_zi[0]))

    # find start_i
    start_i = 0
    for z in z_i:
        if z > float(damelev) - args.zminoffset:
            logging.debug(
                "start_i = :" + str(start_i) + " - search stopped at Zi = " + str(z)
            )
            break
        else:
            start_i = start_i + 1

    i = start_i
    best = -10000
    best_i = i
    best_P = 0
    best_alpha = 0
    best_beta = 0
    l_i = []
    l_z = []
    l_sz = []
    l_P = []
    l_slope = []
    l_mae = []
    l_alpha = []
    l_beta = []

    # Shortcut if just enough data
    data_shortage = False
    # autorise d'avoir un seul point...
    if (i + args.winsize) >= (len(z_i) - 1):
        data_shortage = True
        logging.warning("Just enough data! Compute unique model.")
        alpha, beta, mae, poly = cm.compute_model(z_i, s_zi, z_i[0], s_zi[0])
        logging.debug(
            "Zrange ["
            + str(z_i[0])
            + "; "
            + str(z_i[-1])
            + "]"
            + " --> alpha= "
            + str(alpha)
            + " - beta= "
            + str(beta)
            + " with a local mae of: "
            + str(mae)
            + " m2"
        )
        logging.debug(
            "i: "
            + str(i)
            + " - Slope= "
            + str(poly[0])
            + " - z_med= "
            + str(median(z_i))
            + " - Sz_med= "
            + str(median(s_zi))
        )

        # Select MEA to be used:
        best = mae
        best_P = poly
        best_alpha = alpha
        best_beta = beta

    # Enough data but not in specified distance to the dam
    # (maybe estimated dam elevation is false, maybe offsets are to strict)
    if median(z_i[i : i + args.winsize]) >= args.zmaxoffset + float(damelev):
        data_shortage = True
        logging.error("zmaxoffset to restrictive, no S(Zi) data within search range.")
        logging.error("Maybe dam elevation estimation is not correct.")
        logging.warning(
            "Computing unique model with the first available S(Zi) dataset."
        )
        alpha, beta, mae, poly = cm.compute_model(
            z_i[i : i + args.winsize], s_zi[i : i + args.winsize], z_i[0], s_zi[0]
        )

        logging.debug(
            "i: "
            + str(i)
            + " - Zrange ["
            + str(z_i[i])
            + "; "
            + str(z_i[i + args.winsize])
            + "]"
            + " --> alpha= "
            + str(alpha)
            + " - beta= "
            + str(beta)
            + " with a local mae of: "
            + str(mae)
            + " m2"
        )
        logging.debug(
            "i: "
            + str(i)
            + " - Slope= "
            + str(poly[0])
            + " - z_med= "
            + str(median(z_i[i : i + args.winsize]))
            + " - Sz_med= "
            + str(median(s_zi[i : i + args.winsize]))
        )

        best = mae
        best_i = i
        best_P = poly
        best_alpha = alpha
        best_beta = beta

    while ((i + args.winsize) < (len(z_i) - 1)) and (
        median(z_i[i : i + args.winsize]) < args.zmaxoffset + float(damelev)
    ):
        logging.debug(
            "len(z_i[i:i+args.winsize]): " + str(len(z_i[i : i + args.winsize]))
        )
        #  logging.debug(Zi[i:i+args.winsize])
        # p = np.polyfit(
        #     Zi[i : i + args.winsize],
        #     s_zi[i : i + args.winsize],
        #     1,
        #     rcond=None,
        #     full=False,
        #     w=None,
        #     cov=False,
        # )
        # #  P = np.poly1d(p)

        # beta = (
        #     p[0]
        #     * (median(Zi[i : i + args.winsize]) - Zi[0])
        #     / (median(s_zi[i : i + args.winsize]) - s_zi[0])
        # )
        # alpha = (
        #     p[0] * (math.pow(median(Zi[i : i + args.winsize]) - Zi[0], 1 - beta)) / beta
        # )

        # mae_sum = 0
        # for z, sz in zip(Zi, s_zi):
        #     s = s_zi[0] + alpha * math.pow((z - Zi[0]), beta)
        #     mae_sum += abs(sz - s)
        # # glo_mae = mae_sum / len(Zi)

        # mae_sum = 0
        # for z, sz in zip(Zi[i : i + args.winsize], s_zi[i : i + args.winsize]):
        #     s = s_zi[0] + alpha * math.pow((z - Zi[0]), beta)
        #     mae_sum += abs(sz - s)
        # loc_mae = mae_sum / len(Zi[i : i + args.winsize])
        alpha, beta, loc_mae, poly = cm.compute_model(
            z_i[i : i + args.winsize], s_zi[i : i + args.winsize], z_i[0], s_zi[0]
        )
        logging.debug(
            "i: "
            + str(i)
            + " - Zrange ["
            + str(z_i[i])
            + "; "
            + str(z_i[i + args.winsize])
            + "]"
            + " --> alpha= "
            + str(alpha)
            + " - beta= "
            + str(beta)
            #  + " with a glogal mae of: " +str(glo_mae) +" m2"
            + " with a local mae of: "
            + str(loc_mae)
            + " m2"
        )
        logging.debug(
            "i: "
            + str(i)
            + " - Slope= "
            + str(poly[0])
            + " - z_med= "
            + str(median(z_i[i : i + args.winsize]))
            + " - Sz_med= "
            + str(median(s_zi[i : i + args.winsize]))
        )

        # Select MEA to be used:
        #  mae = glo_mae
        mae = loc_mae
        l_i.append(i)
        l_z.append(median(z_i[i : i + args.winsize]))
        l_sz.append(median(s_zi[i : i + args.winsize]))
        l_P.append(poly)
        l_slope.append(poly[0])
        l_mae.append(mae)
        l_alpha.append(alpha)
        l_beta.append(beta)

        if (best == -10000) or (mae < best):
            best = mae
            best_i = i
            best_P = poly
            best_alpha = alpha
            best_beta = beta
            #  logging.debug("i: "+ str(i)
            #  +" - alpha= " +str(alpha)
            #  +" - beta= "  +str(beta)
            #  +" - mae= " +str(mae))

        i = i + 1

    # For testing
    abs_i = best_i
    # abs_P = best_P
    abs_mae = best
    abs_beta = best_beta
    abs_alpha = best_alpha

    logging.debug("Best Abs")
    logging.debug(
        "i: "
        + str(best_i)
        + " - alpha= "
        + str(best_alpha)
        + " - beta= "
        + str(best_beta)
        + " - mae= "
        + str(best)
    )
    logging.info(
        damname
        + ": S(Z) = "
        + format(s_zi[0], ".2F")
        + " + "
        + format(best_alpha, ".3E")
        + " * ( Z - "
        + format(z_i[0], ".2F")
        + " ) ^ "
        + format(best_beta, ".3E")
    )

    found_first = False
    if args.maemode == "first":
        logging.debug("Reanalizing local mae to find the first local minimum.")
        logging.debug("len(l_mae): " + str(len(l_mae)))
        if len(l_i) > 4:
            if l_mae[0] < l_mae[1] and l_mae[0] < l_mae[2]:

                found_first = True
                logging.debug(
                    "First local minimum found at "
                    + str(l_z[0])
                    + " (i= "
                    + str(l_i[0])
                    + ")."
                )
                best_i = l_i[0]
                best_P = l_P[0]
                best = l_mae[0]
                best_beta = l_beta[0]
                best_alpha = l_alpha[0]
                logging.debug(
                    "i: "
                    + str(best_i)
                    + " - alpha= "
                    + str(best_alpha)
                    + " - beta= "
                    + str(best_beta)
                    + " - mae= "
                    + str(best)
                )

            elif l_mae[1] < l_mae[0] and l_mae[1] < l_mae[2] and l_mae[1] < l_mae[3]:

                found_first = True
                logging.debug(
                    "First local minimum found at "
                    + str(l_z[1])
                    + " (i= "
                    + str(l_i[1])
                    + ")."
                )
                best_i = l_i[1]
                best_P = l_P[1]
                best = l_mae[1]
                best_beta = l_beta[1]
                best_alpha = l_alpha[1]
                logging.debug(
                    "i: "
                    + str(best_i)
                    + " - alpha= "
                    + str(best_alpha)
                    + " - beta= "
                    + str(best_beta)
                    + " - mae= "
                    + str(best)
                )

            else:
                x = range(0, len(l_i) - 1)
                logging.debug(x)
                for j in x[2:-2]:
                    if (
                        l_mae[j] < l_mae[j - 2]
                        and l_mae[j] < l_mae[j - 1]
                        and l_mae[j] < l_mae[j + 1]
                        and l_mae[j] < l_mae[j + 2]
                    ):

                        found_first = True
                        logging.debug(
                            "First local minimum found at "
                            + str(l_z[j])
                            + " (i= "
                            + str(l_i[j])
                            + ")."
                        )
                        best_i = l_i[j]
                        best_P = l_P[j]
                        best = l_mae[j]
                        best_beta = l_beta[j]
                        best_alpha = l_alpha[j]
                        logging.debug(
                            "i: "
                            + str(best_i)
                            + " - alpha= "
                            + str(best_alpha)
                            + " - beta= "
                            + str(best_beta)
                            + " - mae= "
                            + str(best)
                        )
                        break

            if found_first is False:
                logging.info(
                    "Reanalizing local mae to find the first local minimum --> FAILLED."
                )
        else:
            logging.debug("Reanalizing Impossible, not enough local mae data!")

        #  prev_mae_id=0
        #  next_mae_id=4
        #  #  for e in l_mae[2:]:
        #  for e in l_i[2:]:
        #      if (l_mae[e] < l_mae[e-2]) and (l_mae[e] < l_mae[e+2]):
        #          found_first = True
        #          logging.debug("@i= "+str(e)+": "+str(l_mae[e-2])+" > "+str(l_mae[e])+" < "+str(l_mae[e+2]))
        #          logging.debug("First local minimum found at "+ str(l_z[prev_mae_id+2])
        #                        +" (i= "+str(l_i[prev_mae_id+2])+").")
        #          logging.debug("First local minimum found at "+ str(l_z[e])
        #                        +" (i= "+str(e)+").")
        #          #  best_i = l_i[prev_mae_id+2]
        #          best_i = e+start_i-2
        #          best_P = l_P[best_i]
        #          #  best = e
        #          best = l_mae[e]
        #          best_beta = l_beta[best_i]
        #          best_alpha = l_alpha[best_i]
        #          logging.debug("i: "+ str(best_i)
        #                        +" - alpha= " +str(best_alpha)
        #                        +" - beta= "  +str(best_beta)
        #                        +" - mae= " +str(best))
        #          break
        #      else:
        #          prev_mae_id=prev_mae_id+1
        #          next_mae_id=next_mae_id+1
        #  if found_first is False:
        #      logging.debug("Reanalizing local mae to find the first local minimum --> FAILLED.")

    elif args.maemode == "hybrid" and data_shortage is False:
        # first find absolute min at z > damelev > damelev+offset
        logging.debug("Reanalizing local mae to find the optimal local minimum.")
        best = -10000
        x = range(0, len(l_i) - 1)
        for j in x:
            #  logging.debug("j: "+ str(j)
            #  +" - l_z[j]: "+ str(l_z[j]))
            if (
                l_z[j] <= args.zmaxoffset + float(damelev)
                and l_z[j] >= float(damelev)
                and ((best == -10000) or (l_mae[j] < best))
            ):

                best_j = j
                best_z = l_z[j]
                best_i = l_i[j]
                best_P = l_P[j]
                best = l_mae[j]
                best_beta = l_beta[j]
                best_alpha = l_alpha[j]

                logging.debug(
                    "j: "
                    + str(best_j)
                    + " - alpha= "
                    + str(best_alpha)
                    + " - beta= "
                    + str(best_beta)
                    + " - mae= "
                    + str(best)
                    + "@"
                    + str(best_z)
                    + "m"
                )

        logging.info(
            "Hybrid pass #1 => i: "
            + str(best_i)
            + " - alpha= "
            + str(best_alpha)
            + " - beta= "
            + str(best_beta)
            + " - mae= "
            + str(best)
            + "@"
            + str(best_z)
            + "m"
        )

        # then refine with a local lmae minima for smaller Zs
        # and until the slope "breaks"
        x = range(0, best_j - 1)
        ds = np.diff(l_slope) / np.diff(l_z)
        for k in x[::-1]:
            if (
                l_z[k] >= float(damelev) - args.zminoffset
                and abs(ds[k]) < args.dslopethresh
            ):  # and
                #  (l_mae[k] < best)):
                logging.debug(
                    "k: %s - lmae[k]: %s - l_z[k]: %s - ds[k]: %s",
                    str(k),
                    str(l_mae[k]),
                    str(l_z[k]),
                    str(ds[k]),
                )
                best_z = l_z[k]
                best_i = l_i[k]
                best_P = l_P[k]
                best = l_mae[k]
                best_beta = l_beta[k]
                best_alpha = l_alpha[k]
            else:
                break

        logging.info(
            "Hybrid pass #2 => i: "
            + str(best_i)
            + " - alpha= "
            + str(best_alpha)
            + " - beta= "
            + str(best_beta)
            + " - mae= "
            + str(best)
            + "@"
            + str(best_z)
            + "m"
        )

        logging.info("Model updated using hybrid optimal model selection.")
        logging.info(
            "alpha= "
            + str(best_alpha)
            + " - beta= "
            + str(best_beta)
            + " (Computed on the range ["
            + str(z_i[best_i])
            + "; "
            + str(z_i[best_i + args.winsize])
            + "])"
            + " (i= "
            + str(best_i)
            + ")."
        )
        logging.info(
            damname
            + ": S(Z) = "
            + format(s_zi[0], ".2F")
            + " + "
            + format(best_alpha, ".3E")
            + " * ( Z - "
            + format(z_i[0], ".2F")
            + " ) ^ "
            + format(best_beta, ".3E")
        )

    if found_first is True:
        logging.info("Model updated using first LMAE minimum.")
        logging.info(
            "alpha= "
            + str(best_alpha)
            + " - beta= "
            + str(best_beta)
            + " (Computed on the range ["
            + str(z_i[best_i])
            + "; "
            + str(z_i[best_i + args.winsize])
            + "])"
            + " (i= "
            + str(best_i)
            + ")."
        )
        logging.info(
            damname
            + ": S(Z) = "
            + format(s_zi[0], ".2F")
            + " + "
            + format(best_alpha, ".3E")
            + " * ( Z - "
            + format(z_i[0], ".2F")
            + " ) ^ "
            + format(best_beta, ".3E")
        )

    model_json = {
        "ID": dam_id,
        "Name": damname,
        "Elevation": damelev,
        "Model": {
            "Z0": z_i[0],
            "S0": s_zi[0],
            "V0": 0.0,
            "alpha": best_alpha,
            "beta": best_beta,
        },
    }

    with open(os.path.splitext(args.outfile)[0] + ".json", "w") as write_file:
        json.dump(model_json, write_file)

    z = range(int(z_i[0]) + 1, int(z_i[-1]))
    mod_sz = []
    abs_sz = []
    for h in z:
        val_s = s_zi[0] + best_alpha * math.pow((h - z_i[0]), best_beta)
        mod_sz.append(val_s)
        val_s = s_zi[0] + abs_alpha * math.pow((h - z_i[0]), abs_beta)
        abs_sz.append(val_s)

    # reg_abs = np.poly1d(abs_P)
    # reg_best = np.poly1d(best_P)

    # Moldel Plot
    pl.plot_model(
        s_zi[best_i : best_i + args.winsize],
        z_i[best_i : best_i + args.winsize],
        best_alpha,
        best_beta,
        damname,
        args.outfile,
    )
    # fig, ax = plt.subplots()
    # ax.plot(z, sz, "r--", label="S(z)")
    # ax.plot(Zi[0], s_zi[0], "ro")
    # ax.scatter(Zi[1:], s_zi[1:], marker=".", label="S(Zi)")
    # ax.scatter(
    #     median(Zi[best_i : best_i + args.winsize]),
    #     median(s_zi[best_i : best_i + args.winsize]),
    #     color="#ff7f0e",
    #     marker="p",
    #     label="Selected S(Zi)",
    # )
    # #  ax.plot(z, reg_abs(z), 'b-')
    # #  ax.plot(z, reg_best(z), 'g-')
    # ax.grid(b=True, which="major", linestyle="-")
    # ax.grid(b=False, which="minor", linestyle="--")
    # ax.set(
    #     xlabel="Virtual Water Surface Elevation (m)",
    #     ylabel="Virtual Water Surface (ha)",
    # )
    # plt.title(
    #     damname
    #     + ": S(Z) = "
    #     + format(s_zi[0], ".2F")
    #     + " + "
    #     + format(best_alpha, ".3E")
    #     + " * ( Z - "
    #     + format(Zi[0], ".2F")
    #     + " ) ^ "
    #     + format(best_beta, ".3E"),
    #     fontsize=10,
    # )
    # ax.set_xlim(z[0] - 10, z[-1] + 10)
    # # Trick to display in Ha
    # ticks_m2 = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 10000.0))
    # ax.yaxis.set_major_formatter(ticks_m2)
    # plt.minorticks_on()
    # plt.legend(prop={"size": 6}, loc="upper left")
    # plt.savefig(args.outfile, dpi=300)

    # Superimpose local MAE to model
    #  ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    #  ax2.set_ylabel('Local MAE')  # we already handled the x-label with ax1
    #  ax2.set_yscale('log')
    #  ax2.yaxis.set_major_formatter(ticks_m2)
    #  ax2.plot(l_z, l_mae,
    #  marker='.',
    #  linestyle='dashed')
    #  ax2.plot(median(Zi[best_i:best_i+args.winsize]), best, 'p', color='#ff7f0e')
    #  plt.savefig(os.path.splitext(args.outfile)[0]+"_imposed.png", dpi=300)

    # Plot Local MAE
    # TODO : call function plot_slope()
    pl.plot_slope(
        z_i[best_i : best_i + args.winsize],
        l_mae,
        damelev,
        z_i[abs_i : abs_i + args.winsize],
        abs_mae,
        best,
        damname,
        l_z,
        l_slope,
        best_P,
        os.path.splitext(args.outfile)[0] + "_slope.png",
    )
    # ms_fig = plt.figure(dpi=300)
    # ms_fig.subplots_adjust(hspace=0)
    # ms_gs = ms_fig.add_gridspec(3, 1, height_ratios=[1, 1, 1])

    # maeax = plt.subplot(ms_gs[0])
    # #  maefig, maeax = plt.subplots()
    # maeax.axvline(x=float(damelev), ls="--", lw=1, color="teal", label="Dam elevation")
    # maeax.plot(l_z, l_mae, marker=".", color="purple", linestyle="dashed", label="MAE")
    # maeax.plot(
    #     median(Zi[best_i : best_i + args.winsize]),
    #     best,
    #     "p",
    #     color="#ff7f0e",
    #     label="Selected min(MAE)",
    # )
    # maeax.plot(
    #     median(Zi[abs_i : abs_i + args.winsize]),
    #     abs_mae,
    #     "+",
    #     color="black",
    #     label="Absolute min(MAE)",
    # )
    # maeax.grid(b=True, which="major", linestyle="-")
    # maeax.grid(b=False, which="minor", linestyle="--")
    # maeax.set(xlabel="Virtual Water Surface Elevation (m)", ylabel="Local MAE (ha)")
    # maeax.set_xlim(z[0], z[-1] + 10)
    # maeax.set_yscale("log")
    # maeax.label_outer()
    # # Trick to display in Ha
    # maeax.yaxis.set_major_formatter(ticks_m2)
    # plt.minorticks_on()
    # plt.title(damname + ": Local Maximum Absolute Error")
    # plt.legend(prop={"size": 6})
    # #  if (args.debug is True):
    # #  plt.savefig(os.path.splitext(args.outfile)[0]+"_mae.png", dpi=300)

    # # Plot slope
    # #  slpfig, slpax = plt.subplots()
    # slpax = plt.subplot(ms_gs[1])
    # slpax.axvline(x=float(damelev), ls="--", lw=1, color="teal", label="Dam elevation")
    # slpax.plot(
    #     l_z, l_slope, marker=".", color="purple", linestyle="dashed", label="Slope"
    # )
    # slpax.plot(
    #     median(Zi[best_i : best_i + args.winsize]),
    #     best_P[0],
    #     "p",
    #     color="#ff7f0e",
    #     label="Selected min(MAE)",
    # )
    # slpax.grid(b=True, which="major", linestyle="-")
    # slpax.grid(b=False, which="minor", linestyle="--")
    # slpax.set(xlabel="Virtual Water Surface Elevation (m)", ylabel="Local Slope")
    # slpax.set_xlim(z[0], z[-1] + 10)
    # #  slpax.set_yscale('log')
    # # Trick to display in Ha
    # #  slpax.yaxis.set_major_formatter(ticks_m2)
    # plt.minorticks_on()
    # #  plt.title(damname+": Local Slope")
    # #  plt.legend(prop={'size': 6})
    # dslpax = plt.subplot(ms_gs[2])
    # ds = np.diff(l_slope) / np.diff(l_z)
    # dz = (np.array(l_z)[:-1] + np.array(l_z)[1:]) / 2
    # dslpax.axvline(x=float(damelev), ls="--", lw=1, color="teal", label="Dam elevation")
    # dslpax.plot(dz, ds, marker=".", color="blue", label="dSlope")
    # dslpax.grid(b=True, which="major", linestyle="-")
    # dslpax.grid(b=False, which="minor", linestyle="--")
    # dslpax.set(
    #     xlabel="Virtual Water Surface Elevation (m)", ylabel="Local Slope Derivative"
    # )
    # dslpax.set_xlim(z[0], z[-1] + 10)
    # #  slpax.set_yscale('log')
    # # Trick to display in Ha
    # #  dslpax.yaxis.set_major_formatter(ticks_m2)
    # plt.minorticks_on()

    # if args.debug is True:
    #     plt.savefig(os.path.splitext(args.outfile)[0] + "_slope.png", dpi=300)

    # Combined Local MAE / model plot
    pl.plot_model_combo(
        z_i,
        s_zi,
        z_i[best_i : best_i + args.winsize],
        s_zi[best_i : best_i + args.winsize],
        damelev,
        data_shortage,
        alpha,
        beta,
        damname,
        abs_sz,
        mod_sz,
        l_z,
        l_mae,
        abs_mae,
        best,
        os.path.splitext(args.outfile)[0] + "_combo.png",
    )
    # fig = plt.figure(dpi=300)
    # fig.subplots_adjust(hspace=0)
    # gs = fig.add_gridspec(2, 1, height_ratios=[4, 1])

    # maeax0 = plt.subplot(gs[0])
    # maeax0.axvline(x=float(damelev), ls=":", lw=2, color="teal", label="Dam Elevation")
    # maeax0.plot(z, abs_sz, "g--", label="S(z) abs_model")
    # maeax0.plot(z, sz, "r--", label="S(z)")
    # maeax0.plot(Zi[0], s_zi[0], "ro")
    # maeax0.scatter(Zi, s_zi, marker=".", label="S(Z_i)")
    # maeax0.plot(
    #     Zi[best_i : best_i + args.winsize],
    #     s_zi[best_i : best_i + args.winsize],
    #     marker="o",
    #     ls=":",
    #     lw=1,
    #     color="yellow",
    #     markerfacecolor="blue",
    #     label="Selected S(Z_i)",
    # )
    # maeax0.plot(
    #     median(Zi[best_i : best_i + args.winsize]),
    #     median(s_zi[best_i : best_i + args.winsize]),
    #     color="#ff7f0e",
    #     marker="p",
    #     markerfacecolor="#ff7f0e",
    # )
    # maeax0.grid(b=True, which="major", linestyle="-")
    # maeax0.grid(b=False, which="minor", linestyle="--")
    # maeax0.set(
    #     xlabel="Virtual Water Surface Elevation (m)",
    #     ylabel="Virtual Water Surface (ha)",
    # )
    # maeax0.label_outer()
    # #  maeax0.set_xlim(z[0]-10, z[-1]+10)
    # if data_shortage is False:
    #     maeax0.set_xlim(z[0] - 5, median(Zi[best_i : best_i + args.winsize]) + 30)
    #     maeax0.set_ylim(-10, s_zi[best_i + args.winsize] * 1.1)
    # # Trick to display in Ha
    # maeax0.yaxis.set_major_formatter(ticks_m2)
    # plt.minorticks_on()
    # plt.legend(prop={"size": 6}, loc="upper left")

    # maeax1 = plt.subplot(gs[1])
    # maeax1.axvline(x=float(damelev), ls=":", lw=2, color="teal", label="Dam Elevation")
    # maeax1.plot(l_z, l_mae, color="purple", marker=".", linestyle="dashed", label="MAE")
    # maeax1.plot(
    #     median(Zi[best_i : best_i + args.winsize]),
    #     best,
    #     "p",
    #     color="#ff7f0e",
    #     label="Selected MAE",
    # )
    # maeax1.plot(
    #     median(Zi[abs_i : abs_i + args.winsize]),
    #     abs_mae,
    #     "+",
    #     color="black",
    #     label="Absolute min(MAE)",
    # )
    # maeax1.grid(b=True, which="major", linestyle="-")
    # maeax1.grid(b=True, which="minor", linestyle="--")
    # maeax1.set(xlabel="Virtual Water Surface Elevation (m)", ylabel="Local MAE (ha)")
    # #  maeax1.set_xlim(z[0]-10, z[-1]+10)
    # if data_shortage is False:
    #     maeax1.set_xlim(z[0] - 5, median(Zi[best_i : best_i + args.winsize]) + 30)
    # maeax1.set_yscale("log")
    # # Trick to display in Ha
    # maeax1.yaxis.set_major_formatter(ticks_m2)
    # plt.minorticks_on()
    # plt.legend(prop={"size": 6}, loc="lower left")

    # plt.suptitle(
    #     damname
    #     + ": S(Z) = "
    #     + format(s_zi[0], ".2F")
    #     + " + "
    #     + format(best_alpha, ".3E")
    #     + " * ( Z - "
    #     + format(Zi[0], ".2F")
    #     + " ) ^ "
    #     + format(best_beta, ".3E"),
    #     fontsize=10,
    # )
    # plt.savefig(os.path.splitext(args.outfile)[0] + "_combo.png", dpi=300)

    # V(S)
    pl.plot_vs(
        z_i,
        s_zi,
        damelev,
        best_alpha,
        best_beta,
        damname,
        os.path.splitext(args.outfile)[0] + "_VS.png",
    )
    # S_dam = s_zi[0] + best_alpha * math.pow((float(damelev) - Zi[0]), best_beta)
    # s = range(0, math.ceil(S_dam) + 10000)
    # v0 = 0.0
    # vs = []
    # for a in s:
    #     v = v0 + math.pow(((a - s_zi[0]) / best_alpha), 1 / best_beta) * (
    #         s_zi[0] + ((a - s_zi[0]) / (best_beta + 1.0))
    #     )
    #     vs.append(v)

    # # V(S) Moldel Plot
    # figv, axv = plt.subplots()
    # axv.axvline(x=S_dam, ls=":", lw=2, color="teal", label="Max Dam Surface")
    # axv.plot(s, vs, "r-", label="V(S)")
    # axv.grid(b=True, which="major", linestyle="-")
    # axv.grid(b=False, which="minor", linestyle="--")
    # axv.set(
    #     xlabel="Estimated Water Surface (ha)", ylabel="Modelized Water Volume (hm3)"
    # )
    # plt.title(
    #     damname + " : V=f(S) ",
    #     #  +": S(Z) = "+format(s_zi[0], '.2F')+" + "+format(best_alpha, '.3E')
    #     #  +" * ( Z - "+format(Zi[0], '.2F')+" ) ^ "+ format(best_beta, '.3E'),
    #     fontsize=10,
    # )
    # # Trick to display in Ha
    # axv.xaxis.set_major_formatter(ticks_m2)
    # ticks_m3 = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000000.0))
    # axv.yaxis.set_major_formatter(ticks_m3)
    # plt.minorticks_on()
    # plt.legend(prop={"size": 6}, loc="upper left")
    # plt.savefig(os.path.splitext(args.outfile)[0] + "_VS.png", dpi=300)
    t1_stop = perf_counter()
    logging.info(f"Elapsed time: {t1_stop}s, {t1_start}s")

    logging.info(f"Elapsed time during the whole program in s :{t1_stop-t1_start}s")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
