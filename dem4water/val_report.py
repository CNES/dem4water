#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""This modules provides functions to compare a model to a reference law."""
import argparse
import json
import logging
import math
import os
import statistics
import sys

from dem4water import plot_lib as pl


def secured_mean(vals):
    """Compute mean or return NaN."""
    try:
        result = statistics.mean(vals)
    except statistics.StatisticsError:
        return "NaN"

    return result


def secured_stdev(vals):
    """Compute stdev or return NaN."""
    try:
        result = statistics.stdev(vals)
    except statistics.StatisticsError:
        return "NaN"

    return result


def val_report(infile, outfile, reffile, debug):
    """Compare a model to a reference file."""
    # Silence Mathplotlib related debug messages (font matching)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)

    logging_format = (
        "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
    )
    if debug is True:
        logging.basicConfig(
            stream=sys.stdout, level=logging.DEBUG, format=logging_format
        )
    else:
        logging.basicConfig(
            stream=sys.stdout, level=logging.INFO, format=logging_format
        )
    logging.info("Starting val_report.py")

    with open(infile, encoding="utf-8") as model_in:
        model = json.load(model_in)

    with open(reffile, encoding="utf-8") as ref_in:
        ref_db = json.load(ref_in)

    print("\ninfile =", infile)
    print("reffile =", reffile)

    damid = model["ID"]
    damname = model["Name"]
    damelev = model["Elevation"]
    z_0 = model["Model"]["Z0"]
    s_0 = model["Model"]["S0"]
    v_0 = model["Model"]["V0"]
    alpha = model["Model"]["alpha"]
    beta = model["Model"]["beta"]

    print("\ndamelev =", damelev)
    print("Z0 =", z_0)
    print("S0 =", s_0)
    print("V0 =", v_0)
    print("alpha =", alpha)
    print("beta =", beta)

    logging.info(
        "Model for "
        + damname
        + ": S(Z) = "
        + format(s_0, ".2F")
        + " + "
        + format(alpha, ".3E")
        + " * ( Z - "
        + format(z_0, ".2F")
        + " ) ^ "
        + format(beta, ".3E")
    )

    if str(model["ID"]) not in ref_db:
        logging.error(
            f"No reference model available for {model['ID']} in reference DB."
            " Aborting report generation."
        )
        sys.exit(
            "No reference model available for " + str(model["ID"]) + " in reference DB."
        )

    ref_z0 = ref_db[str(model["ID"])]["Model"]["Z0"]
    ref_s0 = ref_db[str(model["ID"])]["Model"]["S0"]  # original values in m2 in DB
    ref_v0 = ref_db[str(model["ID"])]["Model"]["V0"]  # original values in m3 in DB
    ref_alpha = ref_db[str(model["ID"])]["Model"]["alpha"]
    ref_beta = ref_db[str(model["ID"])]["Model"]["beta"]
    ref_zmax = ref_db[str(model["ID"])]["Model"]["Zmax"]
    ref_zmin = ref_db[str(model["ID"])]["Model"]["Zmin"]
    ref_z25 = ref_zmin + 0.25 * (ref_zmax - ref_zmin)
    ref_z75 = ref_zmin + 0.75 * (ref_zmax - ref_zmin)

    print("\nref_Z0 =", ref_z0)
    print("ref_S0 =", ref_s0)
    print("ref_V0 =", ref_v0)
    print("ref_alpha =", ref_alpha)
    print("ref_beta =", ref_beta)

    logging.info(
        "Reference for "
        + damname
        + ": S(Z) = "
        + format(ref_s0, ".2F")
        + " + "
        + format(ref_alpha, ".3E")
        + " * ( Z - "
        + format(ref_z0, ".2F")
        + " ) ^ "
        + format(ref_beta, ".3E")
    )

    print("\ndamelev  =", damelev)
    print("ref_Zmax =", ref_zmax)

    # Figures:
    z_min = max(int(float(z_0)), int(float(ref_z0)))
    alt = range(z_min, int(float(damelev) * 1.1))

    #  z_ref = range(int(float(ref_z_0)), int(float(damelev)*1.1))
    #  s_m_Zmax = S0 + alpha * math.pow((float(damelev) - z_0), beta)
    s_m_zmax = s_0 + alpha * math.pow((float(ref_zmax) - z_0), beta)

    s_r_zmax_m2 = ref_s0 + ref_alpha * math.pow((float(ref_zmax) - ref_z0), ref_beta)
    s_r_zmin_m2 = ref_s0 + ref_alpha * math.pow((float(ref_zmin) - ref_z0), ref_beta)
    s_r_z25_m2 = ref_s0 + ref_alpha * math.pow((float(ref_z25) - ref_z0), ref_beta)
    s_r_z75_m2 = ref_s0 + ref_alpha * math.pow((float(ref_z75) - ref_z0), ref_beta)

    surf = range(0, math.ceil(s_r_zmax_m2) + 10000, 10000)

    sz_model_scatter = []
    sz_ref_scatter = []
    szg = []
    szh = []
    szm = []
    szl = []
    for elev in alt:
        s_m = s_0 + alpha * math.pow((elev - int(z_0)), beta)
        sz_model_scatter.append(s_m)
        s_r = ref_s0 + ref_alpha * math.pow((elev - int(ref_z0)), ref_beta)
        sz_ref_scatter.append(s_r)  # m2

        if s_r != 0 and elev >= ref_zmin:
            szg.append((s_r - s_m) / (s_r))
            if elev >= ref_z75:
                szh.append((s_r - s_m) / (s_r))
            elif elev >= ref_z25:
                szm.append((s_r - s_m) / (s_r))
            else:
                szl.append((s_r - s_m) / (s_r))

    print("\n== s_m =", s_m)
    print("== s_r =", s_r)

    print("\n== Zmax model =", s_m_zmax)
    print("== Zmax ref   =", s_r_zmax_m2)

    vs_model_scatter = []
    vs_ref_scatter = []

    v_m_zmax = v_0 + math.pow(((s_r_zmax_m2 - s_0) / alpha), 1 / beta) * (
        s_0 + ((s_r_zmax_m2 - s_0) / (beta + 1.0))
    )
    v_r_zmax = ref_v0 + math.pow(((s_r_zmax_m2 - ref_s0) / ref_alpha), 1 / ref_beta) * (
        ref_s0 + ((s_r_zmax_m2 - ref_s0) / (ref_beta + 1.0))
    )

    tx_model_scatter = []
    tx_ref_scatter = []
    vsg = []
    vsh = []
    vsm = []
    vsl = []
    tsg = []
    tsh = []
    tsm = []
    tsl = []
    for area in surf:
        v_m = v_0 + math.pow(((area - s_0) / alpha), 1 / beta) * (
            s_0 + ((area - s_0) / (beta + 1.0))
        )
        vs_model_scatter.append(v_m)
        tx_model_scatter.append(v_m / v_m_zmax)

        v_r = ref_v0 + math.pow(((area - ref_s0) / ref_alpha), 1 / ref_beta) * (
            ref_s0 + ((area - ref_s0) / (ref_beta + 1.0))
        )
        vs_ref_scatter.append(v_r)
        tx_ref_scatter.append(v_r / (v_r_zmax))

        if v_r != 0 and area >= s_r_zmin_m2:
            vsg.append((v_r - v_m) / (v_r))
            tsg.append(((v_r / v_r_zmax) - (v_m / v_m_zmax)) / (v_r / v_r_zmax))
            if area >= s_r_z75_m2:
                vsh.append((v_r - v_m) / (v_r))
                tsh.append(((v_r / v_r_zmax) - (v_m / v_m_zmax)) / (v_r / v_r_zmax))
            if area >= s_r_z25_m2:
                vsm.append((v_r - v_m) / (v_r))
                tsm.append(((v_r / v_r_zmax) - (v_m / v_m_zmax)) / (v_r / v_r_zmax))
            else:
                vsl.append((v_r - v_m) / (v_r))
                tsl.append(((v_r / v_r_zmax) - (v_m / v_m_zmax)) / (v_r / v_r_zmax))

    print("\n== v_m_Zmax =", v_m_zmax)
    print("== v_r_Zmax =", v_r_zmax)

    pl.plot_report_sz(
        alt,
        sz_ref_scatter,
        sz_model_scatter,
        ref_zmax,
        damelev,
        damname,
        os.path.splitext(outfile)[0] + "_Sz.png",
    )

    pl.plot_report_vs(
        surf,
        vs_ref_scatter,
        vs_model_scatter,
        s_r_zmax_m2,
        damname,
        os.path.splitext(outfile)[0] + "_Vs.png",
    )

    pl.plot_report_volume_rate(
        surf,
        tx_ref_scatter,
        tx_model_scatter,
        s_r_zmax_m2,
        damname,
        os.path.splitext(outfile)[0] + "_VolumeRate.png",
    )

    results_json = {
        "ID": damid,
        "Name": damname,
        "Zmin": ref_zmin,
        "Z_25": ref_z25,
        "Z_75": ref_z75,
        "Zmax": ref_zmax,
        "Smin": s_r_zmin_m2,
        "S_25": s_r_z25_m2,
        "S_75": s_r_z75_m2,
        "Smax": s_r_zmax_m2,
        "S(z)_quality": {
            "glob": {"mean": secured_mean(szg), "stdev": secured_stdev(szg)},
            "high": {"mean": secured_mean(szh), "stdev": secured_stdev(szh)},
            "mid": {"mean": secured_mean(szm), "stdev": secured_stdev(szm)},
            "low": {"mean": secured_mean(szl), "stdev": secured_stdev(szl)},
        },
        "V(S)_quality": {
            "glob": {"mean": secured_mean(vsg), "stdev": secured_stdev(vsg)},
            "high": {"mean": secured_mean(vsh), "stdev": secured_stdev(vsh)},
            "mid": {"mean": secured_mean(vsm), "stdev": secured_stdev(vsm)},
            "low": {"mean": secured_mean(vsl), "stdev": secured_stdev(vsl)},
        },
        "Vr(S)_quality": {
            "glob": {"mean": secured_mean(tsg), "stdev": secured_stdev(tsg)},
            "high": {"mean": secured_mean(tsh), "stdev": secured_stdev(tsh)},
            "mid": {"mean": secured_mean(tsm), "stdev": secured_stdev(tsm)},
            "low": {"mean": secured_mean(tsl), "stdev": secured_stdev(tsl)},
        },
        "Dam_bottom_estimation": {
            "Z0_mod": z_0,
            "Z0_ref": ref_z0,
            "Dam_bottom_error": abs(z_0 - ref_z0),
        },
    }

    with open(
        os.path.splitext(outfile)[0] + ".json", "w", encoding="utf-8"
    ) as write_file:
        json.dump(results_json, write_file, indent=4)


def val_report_parameters():
    """Define val_report.py parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-i", "--infile", help="dam_model.json file in SI units")
    parser.add_argument("-r", "--reffile", help="validation_DB.json file in SI units")
    parser.add_argument("-o", "--outfile", help="Report.png file")
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
    return parser


def main():
    """Cli function to val_report."""
    parser = val_report_parameters()
    args = parser.parse_args()
    val_report(args.infile, args.outfile, args.reffile, args.debug)


if __name__ == "__main__":
    sys.exit(main())
