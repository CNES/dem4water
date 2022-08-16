#!/usr/bin/env python3
"""
:author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
:organization: CS Group
:copyright: 2020 CS Group. All rights reserved.
:license: see LICENSE file
:created: 2020
"""


import argparse
import json
import logging
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import shapely
import shapely.wkt
from osgeo import gdal, ogr, osr
from shapely.geometry import shape
from shapely.ops import polygonize, split, unary_union


def main(arguments):  # noqa: C901  #FIXME: Function is too complex
    """cut_contourliness.py
    Cut contour lines based on the cutline to estimate the virtual water surface.
    """

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
    args = parser.parse_args(arguments)

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
    logging.info("Starting cut_contourlines.py")

    # Silence Mathplotlib related debug messages (font matching)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)

    geo = osr.SpatialReference()
    geo.ImportFromEPSG(4326)

    ds = gdal.Open(args.dem, gdal.GA_ReadOnly)
    carto = osr.SpatialReference(wkt=ds.GetProjection())

    cartotogeo = osr.CoordinateTransformation(carto, geo)

    # load GeoJSON file containing info
    with open(args.info) as i:
        jsi = json.load(i)
    for feature in jsi["features"]:
        if feature["properties"]["name"] == "Dam":
            logging.debug(feature)
            # dam = shape(feature["geometry"])
            dam_elev = float(feature["properties"]["elev"])
            damname = feature["properties"]["damname"]
            dam_path = damname.replace(" ", "_")

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
                    'gdallocationinfo -valonly -wgs84 "%s" %s %s'
                    % (args.dem, pdblon, pdblat)
                ).read()
            )

            logging.debug("Coordinates (carto): " + str(pdbin.x) + " - " + str(pdbin.y))
            logging.debug("Coordinates (latlon): " + str(pdblat) + " - " + str(pdblon))
            logging.info(
                "PDB detected: "
                + " [pdbLat: "
                + str(pdblat)
                + ", pdbLon: "
                + str(pdblon)
                + ", pdbAlt: "
                + str(pdb_elev)
                + "]"
            )

        if feature["properties"]["name"] == "Insider":
            logging.debug(feature)
            in_w = shape(feature["geometry"])

    try:
        in_w
    except NameError:
        logging.error(
            "Point inside water body for dam "
            + damname
            + " is not present in "
            + args.info
            + ". Can not process."
        )

    drv = ogr.GetDriverByName("GeoJSON")
    if os.path.exists(os.path.join(args.out, damname + "_vSurfaces.json")):
        os.remove(os.path.join(args.out, damname + "_vSurfaces.json"))
    dst_ds = drv.CreateDataSource(os.path.join(args.out, damname + "_vSurfaces.json"))
    dst_layer = dst_ds.CreateLayer("", srs=carto, geom_type=ogr.wkbPolygon)
    field_defn_id = ogr.FieldDefn("ID", ogr.OFTString)
    field_defn = ogr.FieldDefn("level", ogr.OFTString)
    dst_layer.CreateField(field_defn_id)
    dst_layer.CreateField(field_defn)

    # load GeoJSON file containing cutline
    with open(args.cut) as c:
        jsc = json.load(c)
    for feature in jsc["features"]:
        lines = shape(feature["geometry"])
        logging.debug(len(lines.geoms))
        lines = shapely.ops.linemerge(lines)

    # Fixing linemerge not merging every part of MultiLineString
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
        logging.debug("Loop detected in cutline. Looping complexity: " + str(len(poly)))

        line = shapely.geometry.LineString(sorted(line.coords))
        logging.debug(line.length)

        if line.is_simple:
            logging.info("A loop was detected in cutline. It has been simplified.")

    if args.debug is True:
        dbg_simplified_cutline = shapely.geometry.mapping(line)
        dbg_simplified_cutline["crs"] = jsc["crs"]
        with open(os.path.join(args.out, "simplified_cutline.geojson"), "w") as outfile:
            json.dump(dbg_simplified_cutline, outfile)

    if args.level is None:
        logging.debug("No contour line provided, generating to cache.")
        # Generate contour lines from DEM
        logging.debug(
            "args.cache: "
            + args.cache
            + " - dam_path: "
            + dam_path
            + " - args.elevsampling: "
            + args.elevsampling
        )
        contourline_fname = os.path.join(
            args.cache,
            dam_path + "_contourlines@" + str(args.elevsampling) + "m.json",
        )
        logging.debug("contourline_fname: " + contourline_fname)
        elev_margin = 3 * args.elevsampling
        target_elev = dam_elev + args.elevoffset

        if os.path.exists(contourline_fname):
            os.remove(contourline_fname)

        os.system(
            './gen_contourline_polygons.sh "%s" "%s" "%s" "%s" "%s" "%s"'
            % (
                args.dem,
                str(int(pdb_elev - elev_margin)),
                str(args.elevsampling),
                str(int(target_elev + elev_margin)),
                contourline_fname,
                args.tmp,
            )
        )

        with open(contourline_fname) as lvl:
            jsl = json.load(lvl)

    else:
        # If provided, load GeoJSON file containing contour lines
        logging.debug("Using provided contour line file.")
        with open(args.level) as lvl:
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
        for p in results:
            if p.contains(in_w):
                max_area = p.area
                max_elev = float(feature["properties"]["ID"])
                found = True
                logging.info(
                    "Elevation: "
                    + str(feature["properties"]["ID"])
                    + " m - Area: "
                    + str(p.area)
                    + " m2"
                )
                r_feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
                r_p = ogr.CreateGeometryFromWkt(p.wkt)
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
            logging.debug(
                "No relevant polygon found for Elevation "
                + str(feature["properties"]["ID"])
                + " m"
            )

    logging.debug("Identified levels: " + str(r_id))

    r_elev.append(pdb_elev)
    r_area.append(0.0)

    # make sure data is sorted by elevation
    r_elev = np.sort(r_elev)[::-1]
    r_area = np.sort(r_area)[::-1]

    fig, ax = plt.subplots()
    # Trick to display in Ha
    ticks_y = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 10000.0))
    ax.yaxis.set_major_formatter(ticks_y)
    ax.plot(
        r_elev[:-1],
        r_area[:-1],
        color="#ff7f0e",
        marker="o",
        linestyle="dashed",
        markerfacecolor="blue",
    )
    ax.set_ylim(bottom=0.0)
    ax.set(
        xlabel="Virtual Water Surface Elevation (m)",
        ylabel="Virtual Water Surface (ha)",
    )
    ax.plot(pdb_elev, 0.0, "ro")
    ax.grid(b=True, which="major", linestyle="-")
    ax.grid(b=False, which="minor", linestyle="--")
    plt.minorticks_on()
    plt.title(damname + ": S(Z_i)")
    fig.savefig(os.path.join(args.out, damname + "_SZi.png"))

    data = np.column_stack((r_elev, r_area))
    np.savetxt(os.path.join(args.out, damname + "_SZi.dat"), data)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
