#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""This module extract virtual surface using DEM and cutline."""
import argparse
import json
import logging
import os
import sys
from time import perf_counter

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio as rio
import shapely
import shapely.wkt
from osgeo import gdal, ogr, osr
from shapely.geometry import shape
from shapely.ops import polygonize, split, unary_union

from dem4water.plot_lib import plot_szi_points

logger = logging.getLogger("cut_contourlines")
log = logging.getLogger()
log.setLevel(logging.ERROR)
logging.getLogger("geopandas").setLevel(logging.WARNING)
logging.getLogger("rasterio").setLevel(logging.WARNING)
logging.getLogger("fiona").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.ERROR)


def extract_dam_info(dam_info):
    gdf = gpd.read_file(dam_info)
    insider_info = gdf.loc[gdf.name == "Insider"]
    in_w = insider_info.geometry.values[0]
    dam_name = str(insider_info.damname.values[0])
    dam_path = dam_name.replace(" ", "-")
    d_temp = dict(zip(gdf.name, gdf.elev))
    dam_elev = d_temp["Dam"]
    pdb_elev = d_temp["PDB"]
    print(dam_name, in_w, dam_elev, pdb_elev)
    return dam_name, dam_path, dam_elev, pdb_elev, in_w


def load_info_file(info, cartotogeo, dem):
    """Load info from daminfo file."""
    # load GeoJSON file containing info
    with open(info, encoding="utf-8") as i:
        jsi = json.load(i)
    for feature in jsi["features"]:
        if feature["properties"]["name"] == "Dam":
            logger.debug(feature)
            # dam = shape(feature["geometry"])
            dam_elev = float(feature["properties"]["elev"])
            damname = feature["properties"]["damname"]
            # damname = damname.replace("-","_")
            dam_path = damname.replace(" ", "-")
            # dam_path  = dam_path.replace("-","_")
        if feature["properties"]["name"] == "PDB":
            logger.debug(feature)
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

            logger.debug(f"Coordinates (carto): {pdbin.x} - {pdbin.y}")
            logger.debug(f"Coordinates (latlon): {pdblat} - {pdblon})")
            logger.info(
                "PDB detected: "
                f" [pdbLat: {pdblat}, pdbLon: {pdblon}, pdbAlt: {pdb_elev}]"
            )

        if feature["properties"]["name"] == "Insider":
            logger.debug(feature)
            in_w = shape(feature["geometry"])

    try:
        in_w
    except NameError:
        logger.error(
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
        logger.debug(line.length)

        poly = list(polygonize(unary_union(line)))
        logger.debug(f"Loop detected in cutline. Looping complexity: {len(poly)}")
        # This is not correct. If the line is growing on y-axis it will produce nothing good
        line = shapely.geometry.LineString(sorted(line.coords))
        logger.debug(line.length)

        if line.is_simple:
            logger.info("A loop was detected in cutline. It has been simplified.")

    if debug is True:
        dbg_simplified_cutline = shapely.geometry.mapping(line)
        dbg_simplified_cutline["crs"] = jsc["crs"]
        with open(
            os.path.join(out, "simplified_cutline.geojson"), "w", encoding="utf-8"
        ) as outfile:
            json.dump(dbg_simplified_cutline, outfile)

    return line


def ensure_elev_in_dem(dem, start_elev, target_elev, pdb_alt):
    """

    Parameters
    ----------
    dem
    start_elev
    target_elev

    Returns
    -------

    """
    with rio.open(dem) as dem_raster:
        dem_array = dem_raster.read()
        min_elev = np.nanmin(dem_array)
        max_elev = np.nanmax(dem_array)
        if start_elev < min_elev:
            logger.info(
                f"Targeted start elev was {start_elev} "
                f"but {min_elev} found in dem. New bound defined"
            )
            start_elev = min_elev
        if start_elev < pdb_alt:
            logger.info("The start elev must no be lower than the PDB.")
            start_elev = pdb_alt
        if target_elev > max_elev:
            logger.info(
                f"Targeted end_elev was {target_elev} "
                f"but {max_elev} found in dem. New bound defined"
            )
            target_elev = max_elev
        return start_elev, target_elev


def create_contour_lines(
    dem, elev_sampling, start_elev, end_elev, tmp_path, level_file
):
    ds = gdal.Open(dem)
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    list_gdf = []
    list_temp_file_to_remove = []
    for elev in range(start_elev, end_elev, elev_sampling):
        output_shp = os.path.join(tmp_path, f"contour_{elev}.shp")
        ogr_ds = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(output_shp)
        ogr_lyr = ogr_ds.CreateLayer("contour", geom_type=ogr.wkbMultiPolygon, srs=proj)
        field_defn = ogr.FieldDefn("ID", ogr.OFTInteger)
        ogr_lyr.CreateField(field_defn)
        field_defn = ogr.FieldDefn("elevMin", ogr.OFTReal)
        ogr_lyr.CreateField(field_defn)
        field_defn = ogr.FieldDefn("level", ogr.OFTReal)
        ogr_lyr.CreateField(field_defn)

        gdal.ContourGenerateEx(
            ds.GetRasterBand(1),
            ogr_lyr,
            options=[
                f"FIXED_LEVELS={elev}",
                "ID_FIELD=0",
                "ELEV_FIELD_MIN=1",
                "ELEV_FIELD_MAX=2",
                "POLYGONIZE=TRUE",
                f"NODATA={ds.GetRasterBand(1).GetNoDataValue()}",
            ],
        )
        ogr_ds = None
        del ogr_ds
        # Remove small surface not related with the biggest
        gdf = gpd.read_file(output_shp)
        gdf = gdf.loc[gdf.level == elev]
        gdf = gdf.explode(ignore_index=True)
        # Ensure no multipolygon
        gdf = gdf.loc[gdf.area == np.max(gdf.area)]
        gdf.to_file(output_shp)

        list_temp_file_to_remove.append(output_shp)
        list_gdf.append(gdf)
    # reverse to store them in a visual convenience

    list_gdf.reverse()
    gdf_final = pd.concat(list_gdf)
    gdf_final.to_file(level_file)
    for file_to_rm in list_temp_file_to_remove:
        os.remove(file_to_rm)


def generate_countourlines(
    cache, dam_path, elev_sampling, dam_elev, elevoffset, dem, pdb_elev, tmp
):
    """Generate countourlines using gdal."""
    logger.debug("No contour line provided, generating to cache.")
    # Generate contour lines from DEM
    logger.debug(
        f"cache: {cache} - dam_path: {dam_path} - args.elevsampling: {elev_sampling}"
    )
    level_file = os.path.join(
        cache,
        f"{dam_path}_contourlines@{elev_sampling}m.geojson",
    )
    logger.debug(f"contourline_fname: {level_file}")

    # TODO: set to 0 as it increase the surface over the cutline
    elev_margin = 0  # 3 * elevsampling
    target_elev = dam_elev + elevoffset

    if os.path.exists(level_file):
        os.remove(level_file)
    start_elev = int(pdb_elev - elev_margin)
    end_elev = int(target_elev + elev_margin)
    start_elev, end_elev = ensure_elev_in_dem(dem, start_elev, end_elev, pdb_elev)
    logger.info("gen_contourline_polygons.sh parameters: ")
    logger.info(f"dem: {dem} ")
    logger.info(f"start elev: {start_elev}")
    logger.info(f"elevsampling: {elev_sampling} ")
    logger.info(f"end elev: {end_elev} ")
    logger.info(f"output file: {level_file} ")
    logger.info(f"TMPDIR: {tmp} ")

    if start_elev > end_elev:
        raise ValueError(
            f"Start elevation {start_elev} is upper than target_elev {end_elev}"
        )
    create_contour_lines(dem, elev_sampling, start_elev, end_elev, cache, level_file)
    # path for auxillary script
    # script_path = os.path.dirname(__file__)
    # os.system(
    #     f"{script_path}/gen_contourline_polygons.sh {dem} {int(pdb_elev - elev_margin)} "
    #     f"{elevsampling} {int(target_elev + elev_margin)} {contourline_fname} {tmp}"
    # )
    return level_file


def cut_countourlines(
    info,
    dem,
    cutline,
    level,
    elevoffset,
    elevsampling,
    cache,
    tmp,
    out,
    mode,
    debug=False,
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
    logger.info("Starting cut_contourlines.py")

    # Silence Mathplotlib related debug messages (font matching)

    geo = osr.SpatialReference()
    geo.ImportFromEPSG(4326)

    dataset = gdal.Open(dem, gdal.GA_ReadOnly)
    carto = osr.SpatialReference(wkt=dataset.GetProjection())

    # cartotogeo = osr.CoordinateTransformation(carto, geo)

    damname, dam_path, dam_elev, pdb_elev, in_w = extract_dam_info(info)
    # load_info_file(info, cartotogeo, dem)

    drv = ogr.GetDriverByName("GeoJSON")
    if os.path.exists(os.path.join(out, damname + "_vSurfaces.geojson")):
        os.remove(os.path.join(out, damname + "_vSurfaces.geojson"))
    dst_ds = drv.CreateDataSource(os.path.join(out, damname + "_vSurfaces.geojson"))
    dst_layer = dst_ds.CreateLayer("", srs=carto, geom_type=ogr.wkbPolygon)
    field_defn_id = ogr.FieldDefn("ID", ogr.OFTString)
    field_defn = ogr.FieldDefn("level", ogr.OFTString)
    dst_layer.CreateField(field_defn_id)
    dst_layer.CreateField(field_defn)
    if mode == "GDP":
        gdf_cutline = gpd.read_file(cutline)
        line = gdf_cutline.geometry.values[0]
    else:
        # load GeoJSON file containing cutline
        with open(cutline, "r", encoding="utf-8") as cutline_in:
            jsc = json.load(cutline_in)
        for feature in jsc["features"]:
            lines = shape(feature["geometry"])
            logger.debug(len(lines.geoms))
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

    # If provided, load GeoJSON file containing contour lines
    logger.debug("Using provided contour line file.")
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
        for poly in results.geoms:
            if poly.contains(in_w):
                max_area = poly.area
                max_elev = float(feature["properties"]["level"])
                found = True
                logger.info(f"Elevation: {max_elev}m - Area: {poly.area} m2")
                r_feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
                r_p = ogr.CreateGeometryFromWkt(poly.wkt)
                r_feat.SetGeometryDirectly(r_p)
                r_feat.SetField("ID", str(r_id))
                r_feat.SetField("level", max_elev)

        if found is True:
            dst_layer.CreateFeature(r_feat)
            r_feat.Destroy()
            r_elev.append(max_elev)
            r_area.append(max_area)
            r_id = r_id + 1
        else:
            logger.debug(f"No relevant polygon found for Elevation {max_elev} m")

    logger.debug(f"Identified levels: {r_id}")

    r_elev.append(pdb_elev)
    r_area.append(0.0)

    # make sure data is sorted by elevation
    # r_elev = np.sort(r_elev)[::-1]
    # r_area = np.sort(r_area)[::-1]

    plot_szi_points(
        r_elev, r_area, pdb_elev, damname, os.path.join(out, damname + "_SZi.png")
    )

    data = np.column_stack((r_elev, r_area))
    np.savetxt(os.path.join(out, damname + "_SZi.dat"), data)
    t1_stop = perf_counter()
    logger.info(f"Elapsed time: {t1_stop}s {t1_start}s")

    logger.info(f"Elapsed time during the whole program in s : {t1_stop-t1_start}s")


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
