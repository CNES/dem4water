#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""This module search the PDB and the cutline."""
import argparse
import logging
import math
import os
import sys
from itertools import groupby

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from shapely.geometry import LineString, MultiPoint, Point

from dem4water.tools.compute_grandient_dot_product import compute_gradient_product
from dem4water.tools.cutlines_tools import find_extent_to_line
from dem4water.tools.polygonize_raster import polygonize
from dem4water.tools.rasterize_vectors import RasterizarionParams, rasterize
from dem4water.tools.remove_holes_in_shapes import close_holes

logger = logging.getLogger("find_cutline_and_pdb")
log = logging.getLogger()
log.setLevel(logging.ERROR)
logging.getLogger("geopandas").setLevel(logging.WARNING)

logging.getLogger("rasterio").setLevel(logging.WARNING)

logging.getLogger("fiona").setLevel(logging.WARNING)
logger = logging.getLogger("find_cutline_and_pdb")


# ########################################
# Database utils
# ########################################
def remove_mutlipolygon(in_vector, epsg, buffer_size=10):
    """Apply a buffer in both sign to close small gap between polygons."""
    gdf = gpd.GeoDataFrame().from_file(in_vector)
    # print(epsg)
    gdf = gdf.to_crs(epsg)
    gdf["geometry"] = gdf.geometry.buffer(buffer_size)
    gdf["geometry"] = gdf.geometry.buffer(-buffer_size)
    gdf = gdf.explode(ignore_index=True)
    return gdf


def preprocess_water_body(in_vector, epsg, simplify=False, buff_simplify=200):
    """Simplify the water body to GDP."""
    # Fuse multipolygon
    gdf_wb = remove_mutlipolygon(in_vector, epsg, 10)
    # Remove holes
    gdf_wb.geometry = gdf_wb.geometry.apply(lambda p: close_holes(p))
    if simplify:
        gdf_wb.geometry = gdf_wb.geometry.buffer(buff_simplify)
        gdf_wb.geometry = gdf_wb.geometry.buffer(-buff_simplify)
    if len(gdf_wb.index) > 1:
        logger.error("WARNING: multiple water body detected, process only the larger")
        gdf_wb = gdf_wb[gdf_wb.geometry.area == max(gdf_wb.geometry.area)]
    return gdf_wb


# ##############################################
# GDP utils
# ##############################################
def clear_polygonize(in_vector, epsg):
    """Remove all 0 polygons from a vector file."""
    gdf = gpd.GeoDataFrame().from_file(in_vector)
    gdf = gdf.to_crs(epsg)
    if gdf.empty:
        return None
    gdf = gdf.loc[gdf["DN"] == 1]
    # input(gdf.crs)
    return gdf


def merge_close_gdp(gdf_gdp, max_dist):
    """Merge gdp entity if they are close enough."""
    # find intersection according the max distance
    gdf_gdp.geometry = gdf_gdp.geometry.buffer(max_dist)
    # then fuse
    gdf_gdp = gdf_gdp.dissolve(by="DN")
    # remove the buffer
    gdf_gdp.geometry = gdf_gdp.geometry.buffer(-max_dist)
    # then separe into unique object
    gdf_gdp = gdf_gdp.explode(ignore_index=True)
    # dissolve remove the DN column
    gdf_gdp["gdp_unique_id"] = range(len(gdf_gdp.index))
    return gdf_gdp


def filter_by_convex_hull(gdf, id_col):
    """Filter insider polygon from decreasing area."""
    gdf_convex_hull = gdf.copy()
    gdf_convex_hull.geometry = gdf_convex_hull.geometry.convex_hull
    gdf_convex_hull["area"] = gdf_convex_hull.geometry.area
    gdf_convex_hull = gdf_convex_hull.sort_values("area", ascending=False)
    gdf_convex_hull = gdf_convex_hull.reset_index(drop=True)
    id_line = list(gdf_convex_hull[id_col])
    view_id = []
    not_to_see = []
    for line in id_line:
        if line not in not_to_see:
            df_work = gdf_convex_hull.loc[gdf_convex_hull[id_col] == line]
            intersect = gpd.sjoin(gdf, df_work, how="inner")
            view_id.append(line)
            not_to_see += list(intersect[f"{id_col}_left"])
    gdf_filtered = gdf[gdf[id_col].isin(view_id)]
    return gdf_filtered


# ####################################################
# Format conversion
# ####################################################
def oversampling_polygon_boundary(gdf_wb, max_dist):
    """Ensure that each point of boundary is lower than the max_dist parameter.

    Parameters
    ----------
    gdf_wb:
         geodataframe containing only one polygon
    max_dist:
        the maximum distance in meter allowed between two points
    """
    # ensure projection is in meters
    gdf_w = gdf_wb.to_crs(2154)
    list_of_points = [
        Point(x[0], x[1]) for x in gdf_w.geometry.values[0].exterior.coords
    ]

    new_list_of_points = []
    while len(list_of_points) >= 2:
        ori_point = list_of_points.pop(0)
        new_list_of_points.append(ori_point)

        target_point = list_of_points.pop(0)
        line = LineString([ori_point, target_point])
        if max_dist < line.length:
            # Too spaced generate point
            new_point = line.interpolate(max_dist)
            list_of_points = [new_point, target_point] + list_of_points
        else:
            # Close enough
            list_of_points = [target_point] + list_of_points
    # gdf_w.geometry = [LineString(new_list_of_points)]
    # gdf_w = gdf_w.to_crs(gdf_wb.crs)
    gdf_f = gpd.GeoDataFrame(
        {"id_point": range(len(new_list_of_points))},
        geometry=new_list_of_points,
        crs=2154,
    )
    gdf_f = gdf_f.to_crs(gdf_wb.crs)
    return gdf_f


# ########################################################
# Points extraction
# ########################################################
def extract_points_from_raster(raster, poly_dataframe):
    """Extract all points inside a polygon."""
    geoms = poly_dataframe.geometry.values
    with rasterio.open(raster) as src:
        out_image, out_transform = mask(src, geoms, crop=True)
        # reference the pixel centre
        # transf = out_transform * Affine.translation(0.5, 0.5)
        transformer = rasterio.transform.AffineTransformer(out_transform)

        no_data = src.nodata
        if no_data is None:
            raise ValueError("The DEM extracted as no nodata value set. Error")
        data = out_image[0, :, :]
        row, col = np.where(data != no_data)
        values = np.extract(data != no_data, data)
        df_p = pd.DataFrame({"col": col, "row": row, "val": values})
        df_p["x"] = df_p.apply(
            lambda row_: transformer.xy(row_.row, row_.col, offset="center")[0],
            axis=1,
        )
        df_p["y"] = df_p.apply(
            lambda row_: transformer.xy(row_.row, row_.col, offset="center")[1],
            axis=1,
        )
        gdf = gpd.GeoDataFrame(
            df_p,
            geometry=df_p.apply(lambda row_: Point(row_["x"], row_["y"]), axis=1),
            crs=poly_dataframe.crs,
        )
        return gdf


# #######################################################
# PDB
# ######################################################


def find_pdb(poly_dataframe, raster):
    """Find pdb point."""
    gdf_pdb = extract_points_from_raster(raster, poly_dataframe)
    min_dem = gdf_pdb["val"].idxmin()
    value_min = gdf_pdb.loc[min_dem]
    pdb = Point(value_min.x, value_min.y)
    return pdb, gdf_pdb["val"].min()


def find_insider(wb_poly, poly_gdp):
    """Find insider point."""
    poly = poly_gdp.copy()
    poly.geometry = poly.geometry.buffer(30)
    inter = wb_poly.overlay(poly)
    if inter.empty:
        return None
    insider = inter.geometry.representative_point().values[0]
    return insider


def find_dam(wb_poly, pdb_point, insider_point):
    """Find the dam."""
    print(pdb_point, insider_point)
    line = LineString([pdb_point, insider_point])
    poly_wb = wb_poly.geometry.iloc[0]
    inter = line.intersection(poly_wb.boundary)
    # input(inter)
    if isinstance(inter, MultiPoint):
        logger.info("Multipoint dam")
        return None
    return inter


# ##################################################
# Approach using angles
# ##################################################
def get_angle(a, b, c):
    """Provide the angle A between the points A,B,C."""
    ab = math.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]))
    ac = math.sqrt((a[0] - c[0]) * (a[0] - c[0]) + (a[1] - c[1]) * (a[1] - c[1]))
    cb = math.sqrt((c[0] - b[0]) * (c[0] - b[0]) + (c[1] - b[1]) * (c[1] - b[1]))
    ang = math.acos((ac * ac + ab * ab - cb * cb) / (2 * ac * ab))
    return np.degrees(ang)


def find_base_line_using_segments(
    gdf_wb,
    gdf_gdp,
    ident,
    index_line,
    work_dir,
    buffer_size_max=50,
    step_buff=5,
    angle_threshold=130,
):
    """."""
    gdf_simpl = gdf_wb.copy()
    # 1. Generate segments
    # Reduce the number of angles to try to find direction changes
    gdf_simpl.geometry = gdf_simpl.geometry.simplify(1)
    coords_points = gdf_simpl.geometry.values[0].exterior.coords
    simple_coords = list(coords_points)[:-1]
    list_a = simple_coords[:]
    list_b = simple_coords[1:] + [simple_coords[0]]
    list_c = [simple_coords[-1]] + simple_coords[:-1]
    list_rupt = []
    list_angle = []
    for a, b, c in zip(list_a, list_b, list_c):
        list_angle.append(get_angle(a, b, c))
        if get_angle(a, b, c) < angle_threshold:
            list_rupt.append(a)

    coords_points_ori = list(gdf_wb.geometry.values[0].exterior.coords)
    segments = []
    seg = []
    for point in coords_points_ori:
        seg.append(point)
        if point in list_rupt:
            segments.append(seg)
            seg = []
    ori_lines = []
    for seg in segments:
        if len(seg) > 2:
            ori_lines.append(LineString(seg))

    gdf = gpd.GeoDataFrame(
        {"segments": range(len(ori_lines))},
        geometry=ori_lines,
        crs=gdf_wb.crs,
    )
    gdf.to_file(os.path.join(work_dir, "segments.geojson"))
    for buff in range(0, buffer_size_max + step_buff, step_buff):
        gdf_gdp_buff = gdf_gdp.copy()
        gdf_gdp_buff.geometry = gdf_gdp_buff.geometry.buffer(buff)
        inter = gpd.sjoin(gdf, gdf_gdp_buff, how="inner")
        # handle the case of multiple segment intersecting the same GDP
        if inter.empty:
            logger.info("Search base for cutline failed. No intersection found")
            logger.info(
                f"No intersection for buffer size {buff}. "
                f"Try higher (max {buffer_size_max})"
            )
        else:
            logger.info(f"Intersection found for buffer {buff}. Process")
            break
        if buffer_size_max == buff:
            return None, ident
    if not len(inter.index) == 1:
        logger.info("More than one segment intersect the same GDP. Try to fuse them")

    lines = list(inter.geometry.values)
    id_seg = list(inter.segments.values)
    coords_new_seg = []
    group = (n - i for i, n in enumerate(id_seg))
    seg_group = [g for _, (*g,) in groupby(lines, lambda _: next(group))]
    if len(seg_group) == 2:
        print("2 pair of segments")
        seg1 = []
        for s in seg_group[0]:
            seg1 += s.coords
        seg2 = []
        for s in seg_group[1]:
            seg2 += s.coords
        if id_seg[0] == 0 and id_seg[-1] == len(ori_lines) - 1:
            print("Join by 0")
            coords_new_seg = seg2 + seg1
        else:
            coords_new_seg = seg1 + seg2
    else:
        logger.info("Only two segments found.")
        if id_seg[0] == 0 and id_seg[-1] == len(ori_lines) - 1:
            coords_new_seg = list(lines[-1].coords) + list(lines[0].coords)
            for seg in lines[1:-1]:
                coords_new_seg += seg.coords
        else:
            for seg in lines:
                coords_new_seg += seg.coords

    df_coords = pd.DataFrame()
    df_coords["x_cut"] = [x[0] for x in coords_new_seg]
    df_coords["y_cut"] = [x[1] for x in coords_new_seg]
    df_coords["gdf_unique_id"] = index_line
    df_coords["gdf_sub_line"] = ident
    ident += 1
    inter = gpd.GeoDataFrame(
        df_coords,
        geometry=gpd.points_from_xy(df_coords["x_cut"], df_coords["y_cut"]),
        crs=gdf_wb.crs,
    )
    return inter, ident


# #######################################################
# Main function
# #######################################################
def find_cutline_and_pdb(
    database_file,
    dem_raster,
    work_dir,
    gdp_buffer_size,
    radius_search_size,
    maximum_alt,
    id_db,
    dam_name,
    elevoffset,
    debug=False,
):
    """Search the PDB and cutline from the database and the DEM.

    Parameters
    ----------
    database_file:
    dem_raster:
    """
    logger_format = (
        "%(asctime)s - %(filename)s:%(lineno)s - %(levelname)s - %(message)s"
    )
    if debug:
        logging.basicConfig(
            stream=sys.stdout, level=logging.DEBUG, format=logger_format
        )
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, format=logger_format)
    logger.info("Starting search cutline and PDB using GDP")
    maximum_alt = maximum_alt + elevoffset
    # 1. Manage the particularity of database:
    # - Fuse multipolygon
    # - Remove holes inside water body
    # - Smooth waterbody (optional)
    # - Rasterize the waterbody
    with rasterio.open(dem_raster) as dem:
        epsg = dem.crs.to_epsg()
        gdf_wb = preprocess_water_body(database_file, epsg)
        gdf_wb.to_file(os.path.join(work_dir, "bd_clean.geojson"))
        params_raster = RasterizarionParams(
            mode="binary",
            binary_foreground_value=1,
            background_value=0,
            column_field=None,
            dtype="uint8",
            nodata=0,
        )
        waterbody_bin = os.path.join(work_dir, "waterbody_bin.tif")
        rasterize(gdf_wb, dem, params_raster, waterbody_bin)
        logger.info("Ending the water body preparation")
        # 2. Compute the gradient dot product using dem and cleaned water body
        # - The output is a raster then vectorize it over the whole area
        # - Remove the 0 polygons which means no GDP found
        gdp_raster = os.path.join(work_dir, "gdp_raster.tif")
        compute_gradient_product(waterbody_bin, dem_raster, gdp_raster)
        gdp_vector = os.path.join(work_dir, "gdp_vector.geojson")
        polygonize(gdp_raster, gdp_vector)
        gdf_gdp = clear_polygonize(gdp_vector, epsg)
        if gdf_gdp is None:
            logger.info("ERROR: when computing GDP, no slope found")
            logger.info("ERROR: Stopping the chain")
            return None
        logger.info("GDP have run successfuly")
        # 3. Fuse close GDP
        # 5. TODO: filter by area to remove small GDP
        gdf_gdp = gdf_gdp[gdf_gdp.geometry.area > 500]
        gdf_gdp = gdf_gdp.reset_index(drop=True)
        gdf_gdp.to_file(
            os.path.join(work_dir, "gdp_fusion_nearest_remove_small.geojson")
        )
        gdf_gdp = merge_close_gdp(gdf_gdp, gdp_buffer_size)
        if gdf_gdp.empty:
            logger.info("ERROR: Empty dataframe after merging GDP.")
            return None
        gdf_gdp.to_file(os.path.join(work_dir, "gdp_fusion_nearest.geojson"))
        logger.info("Fusion ended. Remove very small object")
        # 4. Filter by convex hull to remove all insider GDP
        gdf_gdp = filter_by_convex_hull(gdf_gdp, "gdp_unique_id")
        gdf_gdp.to_file(os.path.join(work_dir, "gdp_fusion_nearest_convex.geojson"))

        logger.info("End filtering. Look for a PDB")
        # 6 Convert the water body into a sett of point with regular sampling
        list_ident = []
        list_line = []
        ident = 0
        list_pdb = []
        list_alt_pdb = []
        for index in gdf_gdp.index:
            row = gdf_gdp.loc[[index]]
            pdb, alt = find_pdb(row, dem_raster)
            list_pdb.append(pdb)
            if alt >= -100:
                list_alt_pdb.append(alt)
            else:
                list_alt_pdb.append(np.nan)
        print(list_alt_pdb)
        index_min = np.nanargmin(list_alt_pdb)
        pdb_point = list_pdb[index_min]
        logger.info(f"PDB {pdb_point} found at altitude {list_alt_pdb[index_min]}")
        for index in [index_min]:  # gdf_gdp.index:
            row = gdf_gdp.loc[[index]]
            row.to_file(os.path.join(work_dir, "search.geojson"))
            # 7. For each GDP find a PDB
            # pdb = find_pdb(row, dem_raster)
            # 8. For each GDP found a insider point
            logger.info("Look for an insider point")
            insider = find_insider(gdf_wb, row)
            logger.info(f"Insider : {insider}")
            if insider is None:
                continue
            # Find dam
            logging.info("Try to find the DAM")
            dam_point = find_dam(gdf_wb, pdb_point, insider)
            logger.info(f"DAM: {dam_point}")
            # 9. For each GDP draw baselines
            gdf_line, ident = find_base_line_using_segments(
                gdf_wb, row, ident, index, work_dir, angle_threshold=150
            )
            logger.info("Base line found. looking for extent")
            # print("aft", ident)
            if gdf_line is None:
                continue
            # 10. For each baseline find extents

            coords_new_line = find_extent_to_line(
                gdf_line, dem_raster, gdf_wb, radius_search_size, maximum_alt, work_dir
            )
            logger.info("Extents found.")
            list_line.append(LineString(coords_new_line))
            list_ident.append(ident)

        # 11. Merge all cutlines into one file
        if not list_line:
            return None
        logger.info("Fuse all files")
        gdf_final = gpd.GeoDataFrame(
            {"ident_line": list_ident}, geometry=list_line, crs=gdf_wb.crs
        )
        dam_name = dam_name.replace(" ", "-")
        logging.info("Write cutline")
        out_cutline = os.path.join(work_dir, f"{dam_name}_cutline.geojson")
        gdf_final.to_file(out_cutline)
        # Dam info
        print([list_alt_pdb[index_min], maximum_alt - elevoffset, ""])
        gdf_dam = gpd.GeoDataFrame(
            {
                "name": ["PDB", "Dam", "Insider"],
                "ID": [id_db, id_db, id_db],
                "elev": [
                    int(list_alt_pdb[index_min]),
                    int(maximum_alt) - elevoffset,
                    0,
                ],
                "damname": [
                    dam_name.replace(" ", "-"),
                    dam_name.replace(" ", "-"),
                    dam_name.replace(" ", "-"),
                ],
            },
            geometry=[pdb_point, dam_point, insider],
            crs=epsg,
        )
        gdf_dam.to_file(os.path.join(work_dir, f"{dam_name}_daminfo.json"))
        return out_cutline


def find_cutline_and_pdb_args():
    """Define find_cutline_and_pdb parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--database_dam", help="GeoJSON database file")
    parser.add_argument("--dem_raster", help="Input DEM raster")
    parser.add_argument("--work_dir", help="Working directory")
    parser.add_argument("--gdp_buffer_size", help="contourline.json file", default=50)
    parser.add_argument(
        "--radius_search_size",
        type=int,
        default=5,
        help="Radius of the search area for new point in line",
    )
    parser.add_argument(
        "--maximum_alt",
        type=int,
        help="Maximum altitude to reach before ending line.",
    )
    parser.add_argument(
        "--id_db",
        help="DAM id in database",
    )
    parser.add_argument("--elevoffset", help="Offset for maximum altitude", default=20)
    parser.add_argument("--dam_name", help="DAM name in DB")
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
    return parser


def main():
    """Cli function for find_cutline_and_pdb."""
    parser = find_cutline_and_pdb_args()
    args = parser.parse_args()
    find_cutline_and_pdb(
        database_file=args.database_dam,
        dem_raster=args.dem_raster,
        work_dir=args.work_dir,
        gdp_buffer_size=args.gdp_buffer_size,
        radius_search_size=args.radius_search_size,
        maximum_alt=args.maximum_alt,
        id_db=args.id_db,
        dam_name=args.dam_name,
        elevoffset=args.elevoffset,
        debug=args.debug,
    )


if __name__ == "__main__":
    sys.exit(main())
