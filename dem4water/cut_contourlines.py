#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""This module extract virtual surface using DEM and cutline."""


import argparse
import json
import logging
import os
import sys
from time import perf_counter

import numpy as np
import shapely
import shapely.wkt
from osgeo import gdal, ogr, osr
from shapely.geometry import shape
from shapely.ops import polygonize, split, unary_union

from dem4water.plot_lib import plot_szi_points


def load_info_file(info, cartotogeo, dem):
    """Load info from daminfo file."""
    # load GeoJSON file containing info
    with open(info, encoding="utf-8") as i:
        jsi = json.load(i)
    for feature in jsi["features"]:
        if feature["properties"]["name"] == "Dam":
            logging.debug(feature)
            # dam = shape(feature["geometry"])
            dam_elev = float(feature["properties"]["elev"])
            damname = feature["properties"]["damname"]
            # damname = damname.replace("-","_")
            dam_path = damname.replace(" ", "-")
            # dam_path  = dam_path.replace("-","_")
        if feature["properties"]["name"] == "PDB":
            logging.debug(feature)
            pdbin = shape(feature["geometry"])

            pdb = ogr.Geometry(ogr.wkbPoint)
            pdb.AddPoint(float(pdbin.x), float(pdbin.y))
            pdb.Transform(cartotogeo)
            pdblat = pdb.GetX()
            pdblon = pdb.GetY()
            pdb_elev = float(
                os.popen(
                    f'gdallocationinfo -valonly -wgs84 "{dem}" {pdblon} {pdblat}'
                ).read()
            )

            logging.debug(f"Coordinates (carto): {pdbin.x} - {pdbin.y}")
            logging.debug(f"Coordinates (latlon): {pdblat} - {pdblon})")
            logging.info(
                "PDB detected: "
                f" [pdbLat: {pdblat}, pdbLon: {pdblon}, pdbAlt: {pdb_elev}]"
            )

        if feature["properties"]["name"] == "Insider":
            logging.debug(feature)
            in_w = shape(feature["geometry"])

    try:
        in_w
    except NameError:
        logging.error(
            f"Point inside water body for dam {damname} "
            f"is not present in {info}. Can not process."
        )
    return damname, dam_path, dam_elev, pdb_elev, in_w


def manage_cutline(jsc, lines, out, debug):
    """Manage cutline format."""
    if lines.geom_type == "MultiLineString":
        outcoords = [list(i.coords) for i in lines]
        line = shapely.geometry.LineString(
            [i for sublist in outcoords for i in sublist]
        )
        line = shapely.geometry.LineString(line.coords)
    else:
        line = lines

    # Working around looping LineString

    if not line.is_simple:
        logging.debug(line.length)

        poly = list(polygonize(unary_union(line)))
        logging.debug(f"Loop detected in cutline. Looping complexity: {len(poly)}")
        # This is not correct. If the line is growing on y-axis it will produce nothing good
        line = shapely.geometry.LineString(sorted(line.coords))
        logging.debug(line.length)

        if line.is_simple:
            logging.info("A loop was detected in cutline. It has been simplified.")

    if debug is True:
        dbg_simplified_cutline = shapely.geometry.mapping(line)
        dbg_simplified_cutline["crs"] = jsc["crs"]
        with open(
            os.path.join(out, "simplified_cutline.geojson"), "w", encoding="utf-8"
        ) as outfile:
            json.dump(dbg_simplified_cutline, outfile)

    return line


def generate_countourlines(
    cache, dam_path, elevsampling, dam_elev, elevoffset, dem, pdb_elev, tmp
):
    """Generate countourlines using gdal."""
    logging.debug("No contour line provided, generating to cache.")
    # Generate contour lines from DEM
    logging.debug(
        f"cache: {cache} - dam_path: {dam_path} - args.elevsampling: {elevsampling}"
    )
    contourline_fname = os.path.join(
        cache,
        dam_path + f"_contourlines@{elevsampling}m.json",
    )
    logging.debug(f"contourline_fname: {contourline_fname}")
    elev_margin = 3 * elevsampling
    target_elev = dam_elev + elevoffset

    if os.path.exists(contourline_fname):
        os.remove(contourline_fname)

    logging.info("gen_contourline_polygons.sh parameters: ")
    logging.info(f"dem: {dem} ")
    logging.info(f"start elev: {str(int(pdb_elev - elev_margin))}")
    logging.info(f"elevsampling: {elevsampling} ")
    logging.info(f"end elev: {str(int(target_elev + elev_margin))} ")
    logging.info(f"output file: {contourline_fname} ")
    logging.info(f"TMPDIR: {tmp} ")
    start_elev = int(pdb_elev - elev_margin)
    end_elev = int(target_elev + elev_margin)

    if start_elev > end_elev:
        raise ValueError(
            f"Start elevation {start_elev} is upper than target_elev {end_elev}"
        )

    # os.system(
    #     './gen_contourline_polygons.sh "%s" "%s" "%s" "%s" "%s" "%s"'
    #     % (
    #         dem,
    #         str(int(pdb_elev - elev_margin)),
    #         str(elevsampling),
    #         str(int(target_elev + elev_margin)),
    #         contourline_fname,
    #         tmp,
    #     )
    # )
    # path for auxillary script
    script_path = os.path.dirname(__file__)
    os.system(
        f"{script_path}/gen_contourline_polygons.sh {dem} {int(pdb_elev - elev_margin)} "
        f"{elevsampling} {int(target_elev + elev_margin)} {contourline_fname} {tmp}"
    )
    return contourline_fname
    # with open(contourline_fname, encoding="utf-8") as lvl:
    #     jsl = json.load(lvl)
    # return jsl


def cut_countourlines(
    info, dem, cutline, level, elevoffset, elevsampling, cache, tmp, out, debug=False
):
    """Cut contour lines based on the cutline to estimate the virtual water surface."""
    t1_start = perf_counter()
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
    logging.info("Starting cut_contourlines.py")

    # Silence Mathplotlib related debug messages (font matching)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)

    geo = osr.SpatialReference()
    geo.ImportFromEPSG(4326)

    dataset = gdal.Open(dem, gdal.GA_ReadOnly)
    carto = osr.SpatialReference(wkt=dataset.GetProjection())

    cartotogeo = osr.CoordinateTransformation(carto, geo)

    damname, dam_path, dam_elev, pdb_elev, in_w = load_info_file(info, cartotogeo, dem)

    drv = ogr.GetDriverByName("GeoJSON")
    if os.path.exists(os.path.join(out, damname + "_vSurfaces.json")):
        os.remove(os.path.join(out, damname + "_vSurfaces.json"))
    dst_ds = drv.CreateDataSource(os.path.join(out, damname + "_vSurfaces.json"))
    dst_layer = dst_ds.CreateLayer("", srs=carto, geom_type=ogr.wkbPolygon)
    field_defn_id = ogr.FieldDefn("ID", ogr.OFTString)
    field_defn = ogr.FieldDefn("level", ogr.OFTString)
    dst_layer.CreateField(field_defn_id)
    dst_layer.CreateField(field_defn)

    # load GeoJSON file containing cutline
    with open(cutline, "r", encoding="utf-8") as cutline_in:
        jsc = json.load(cutline_in)
    for feature in jsc["features"]:
        lines = shape(feature["geometry"])
        logging.debug(len(lines.geoms))
        lines = shapely.ops.linemerge(lines)

    # Fixing linemerge not merging every part of MultiLineString
    line = manage_cutline(jsc, lines, out, debug)

    if level is None:
        level = generate_countourlines(
            cache,
            dam_path,
            elevsampling,
            dam_elev,
            elevoffset,
            dem,
            pdb_elev,
            tmp,
        )
    # else:
    # If provided, load GeoJSON file containing contour lines
    logging.debug("Using provided contour line file.")
    with open(level, "r", encoding="utf-8") as lvl:
        jsl = json.load(lvl)

    r_id = 1
    r_elev = []
    r_area = []
    for feature in jsl["features"]:
        level = shape(feature["geometry"])
        results = split(level, line)
        found = False
        max_area = -10000
        max_elev = -10000
        # Combien d'éléments dans results?
        # l'indentation parait etrange.
        # a revoir
        for poly in results.geoms:
            if poly.contains(in_w):
                max_area = poly.area
                max_elev = float(feature["properties"]["ID"])
                found = True
                logging.info(f"Elevation: {max_elev}m - Area: {poly.area} m2")
                r_feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
                r_p = ogr.CreateGeometryFromWkt(poly.wkt)
                r_feat.SetGeometryDirectly(r_p)
                r_feat.SetField("ID", str(r_id))
                r_feat.SetField("level", feature["properties"]["ID"])

        if found is True:
            dst_layer.CreateFeature(r_feat)
            r_feat.Destroy()
            r_elev.append(max_elev)
            r_area.append(max_area)
            r_id = r_id + 1
        else:
            logging.debug("No relevant polygon found for Elevation {max_elev} m")

    logging.debug(f"Identified levels: {r_id}")

    r_elev.append(pdb_elev)
    r_area.append(0.0)

    # make sure data is sorted by elevation
    r_elev = np.sort(r_elev)[::-1]
    r_area = np.sort(r_area)[::-1]

    plot_szi_points(
        r_elev, r_area, pdb_elev, damname, os.path.join(out, damname + "_SZi.png")
    )

    data = np.column_stack((r_elev, r_area))
    np.savetxt(os.path.join(out, damname + "_SZi.dat"), data)
    t1_stop = perf_counter()
    logging.info(f"Elapsed time: {t1_stop}s {t1_start}s")

    logging.info(f"Elapsed time during the whole program in s : {t1_stop-t1_start}s")


def cut_countourlines_ars():
    """Define cut_countourlines parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--info", help="daminfo.json file")
    parser.add_argument("-d", "--dem", help="Input DEM")
    parser.add_argument("-c", "--cut", help="cutline.json file")
    parser.add_argument("-l", "--level", help="contourline.json file")
    parser.add_argument(
        "--elevoffset",
        type=float,
        default=50,
        help="Elevation offset target for the cutline wrt the estimated dam elevation",
    )
    parser.add_argument(
        "--elevsampling",
        type=int,
        default=1,
        help="Elevation sampling step for contour lines generation.",
    )
    parser.add_argument(
        "--cache",
        help="Cache directory to store <DAM>_contourlines@*m.json files.",
    )
    parser.add_argument("-t", "--tmp", help="Temporary directory")
    parser.add_argument("-o", "--out", help="Output directory")
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
    return parser


def main():
    """Cli function to launch cut countourlines."""
    parser = cut_countourlines_ars()
    args = parser.parse_args()
    cut_countourlines(
        args.info,
        args.dem,
        args.cut,
        args.level,
        args.elevoffset,
        args.elevsampling,
        args.cache,
        args.tmp,
        args.out,
        args.debug,
    )


if __name__ == "__main__":
    sys.exit(main())
