#!/usr/bin/env python3
"""This module provides function to find the PDB and the cutline."""


import argparse
import json
import logging
import os
import sys
from math import ceil, sqrt
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np

import rasterio as rio
from osgeo import gdal, ogr, osr
from shapely.geometry import shape

from dem4water.tools.extract_roi import ExtractROIParam, extract_roi
from dem4water.tools.save_raster import save_image
from dem4water.tools.superimpose import SuperimposeParam, superimpose

from dem4water.tools.utils import distance

# from PIL import Image


def nderiv(y, x):
    "Différence finie, dérivée de la fonction f."
    n = len(y)
    d = np.zeros(n, "d")  # virgule flottante à double précision (double)
    # différences de part et d'autre
    # centrées sur les points intérieurs
    for i in range(1, n - 1):
        d[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
    # différences sur un seul côté pour les extrémités
    d[0] = (y[1] - y[0]) / (x[1] - x[0])
    d[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])
    return d


def pixel(dx, dy, ds):
    #  print(ds.GetGeoTransform())
    px = ds.GetGeoTransform()[0]
    py = ds.GetGeoTransform()[3]
    rx = ds.GetGeoTransform()[1]
    ry = ds.GetGeoTransform()[5]
    x = (dx - px) / rx
    y = (dy - py) / ry
    x -= 0.5
    y -= 0.5

    return round(x), round(y)


def coord(x, y, ds):
    #  print(ds.GetGeoTransform())
    xoffset, px_w, rot1, yoffset, rot2, px_h = ds.GetGeoTransform()
    dx = px_w * x + rot1 * y + xoffset
    dy = rot2 * x + px_h * y + yoffset
    dx += px_w / 2.0
    dy += px_h / 2.0

    return dx, dy


def points_in_circle(circle, arr):
    "A generator to return all points whose indices are within given circle."
    i0, j0, r = circle

    def intceil(x):
        return int(ceil(x))

    for i in range(intceil(i0 - r), intceil(i0 + r)):
        ri = sqrt(r**2 - (i - i0) ** 2)
        for j in range(intceil(j0 - ri), intceil(j0 + ri)):
            yield arr[i][j]


def find_pdb_and_cutline(
    infile,
    dam_id,
    id_db,
    watermap,
    dem,
    pdbstep,
    pdbradius,
    mradius,
    maxdist,
    elevoffset,
    elevsampling,
    info,
    tmp,
    out,
    radius=None,
    debug=False,
):  # noqa: C901  #FIXME: Function is too complex
    """Find the PDB and create the cutline."""
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
    logging.info("Starting szi_from_contourline.py")

    # Silence Mathplotlib related debug messages (font matching)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)

    # Init global srs:
    geo = osr.SpatialReference()
    geo.ImportFromEPSG(4326)

    ds = gdal.Open(dem, gdal.GA_ReadOnly)
    carto = osr.SpatialReference(wkt=ds.GetProjection())

    geotocarto = osr.CoordinateTransformation(geo, carto)
    cartotogeo = osr.CoordinateTransformation(carto, geo)

    # Start processing by qualifying the Dam
    driver = ogr.GetDriverByName("GeoJSON")
    dataSource = driver.Open(infile, 0)
    layer = dataSource.GetLayer()

    clat = 0
    clon = 0
    calt = 0
    clat_in = 0
    clon_in = 0
    dam_name = ""
    dam_path = ""
    dam_404 = True
    calt_from_DB = False
    # radius = radius

    for feature in layer:
        if str(int(feature.GetField(str(id_db)))) == str(dam_id):
            # Compute radius
            if radius is None:
                geom = feature.GetGeometryRef()
                bbox = geom.GetEnvelope()
                radius = distance(bbox[2], bbox[0], bbox[3], bbox[1])
                logging.info(f"=> RADIUS : {radius}")

            dam_404 = False
            logging.debug(feature.GetField("DAM_NAME"))
            dam_name = feature.GetField("DAM_NAME")
            # dam_id = int(feature.GetField(str(id_db)))
            dam_path = dam_name.replace(" ", "-")
            # dam_path = dam_path.replace("-","_")
            clat = float(feature.GetField("LAT_DD"))
            clon = float(feature.GetField("LONG_DD"))
            # extraire l'altitude directement depuis le MNT

            if bool(feature.GetField("DAM_LVL_M")):
                calt = float(feature.GetField("DAM_LVL_M"))
                calt_from_DB = True
            else:
                logging.warning(
                    "Altitude for dam "
                    + dam_name
                    + "(ID:"
                    + str(id)
                    + ") is not present in "
                    + infile
                    + "."
                )
            if bool(feature.GetField("LAT_WW")):
                clat_in = float(feature.GetField("LAT_WW"))
                clon_in = float(feature.GetField("LONG_WW"))
            else:
                logging.error(
                    "Point inside water body for dam "
                    + dam_name
                    + "(ID:"
                    + str(id)
                    + ") is not present in "
                    + infile
                    + ". Can not process."
                )
            break
    layer.ResetReading()

    if dam_404 is True:
        logging.error(
            "404 - Dam Not Found: "
            + dam_name
            + "(ID:"
            + str(id)
            + ") is not present in "
            + infile
        )

    # Init output vector data files:
    shpDriver = ogr.GetDriverByName("GeoJSON")
    if os.path.exists(os.path.join(out, dam_path + "_cutline.geojson")):
        shpDriver.DeleteDataSource(os.path.join(out, dam_path + "_cutline.geojson"))
    outDataSource = shpDriver.CreateDataSource(
        os.path.join(out, dam_path + "_cutline.geojson")
    )
    outLayer = outDataSource.CreateLayer(
        "", srs=carto, geom_type=ogr.wkbMultiLineString
    )

    drv_dbg = ogr.GetDriverByName("GeoJSON")
    if os.path.exists(os.path.join(tmp, dam_path + "_cutline_points.geojson")):
        os.remove(os.path.join(tmp, dam_path + "_cutline_points.geojson"))
    dbg_ds = drv_dbg.CreateDataSource(
        os.path.join(tmp, dam_path + "_cutline_points.geojson")
    )
    dbg_layer = dbg_ds.CreateLayer("", srs=carto, geom_type=ogr.wkbPoint)
    dbg_field_defn = ogr.FieldDefn("name", ogr.OFTString)
    dbg_layer.CreateField(dbg_field_defn)

    logging.info(
        "Currently processing: "
        + dam_name
        + "(id:"
        + str(id)
        + ") [Lat: "
        + str(clat)
        + ", Lon: "
        + str(clon)
        + "]"
    )

    if calt_from_DB is True:
        logging.info("Alt from DB: " + str(calt))
    else:
        # altcmd = (
        #     'gdallocationinfo -valonly -wgs84 "'
        #     + dem
        #     + '" '
        #     + str(clon)b
        #     + " "
        #     + str(clat)
        # )
        calt = float(
            os.popen(
                'gdallocationinfo -valonly -wgs84 "%s" %s %s' % (dem, clon, clat)
            ).read()
        )
        logging.info("Alt from DEM: " + str(calt))

    logging.info(
        "Currently processing: "
        + dam_name
        + "(id:"
        + str(id)
        + ") [Lat: "
        + str(clat)
        + ", Lon: "
        + str(clon)
        + ", Alt: "
        + str(calt)
        + "]"
    )

    dam = ogr.Geometry(ogr.wkbPoint)
    dam.AddPoint(clat, clon)
    dam.Transform(geotocarto)
    logging.debug("Coordinates Carto: " + str(dam.GetX()) + " - " + str(dam.GetY()))

    # Dam elevation from watermap
    if calt_from_DB is True:
        logging.info("Alt Extracted from DB.")
        bml_alt = calt
    else:
        logging.info("Alt Estimated from Watermap.")
        extract_roi_parameters_ext = ExtractROIParam(
            mode="radius",
            mode_radius_r=pdbradius,
            mode_radius_unitr="phy",
            mode_radius_unitc="phy",
            mode_radius_cx=dam.GetX(),
            mode_radius_cy=dam.GetY(),
            dtype="float",
        )

        ext, profile_ext = extract_roi(rio.open(dem), extract_roi_parameters_ext)
        save_image(ext, profile_ext, os.path.join(out, "dem_pdb.tif"))

        superimpose_extw = SuperimposeParam(interpolator="bco", dtype="float")
        extw, profile_extw = superimpose(
            rio.open(watermap), ext, superimpose_extw, profile_ext, None
        )
        bml = np.where(extw > 0.05, ext, 0)
        np_bml = bml.reshape(bml.shape[1], bml.shape[2])


        bml_alt = np.amax(np_bml)

    targetelev = bml_alt + elevoffset
    if calt_from_DB is True:
        logging.info("Extracted dam elevation= " + str(bml_alt) + "m")
    else:
        logging.info("Estimated dam elevation= " + str(bml_alt) + "m")
    logging.info("Target Elevation for cutline search= " + str(targetelev) + "m")

    # BEGIN:
    # Can be conflincting with other launcher mode
    # Set info to default is enough?
    if info is None:
        drv = ogr.GetDriverByName("GeoJSON")
        if os.path.exists(os.path.join(out, dam_path + "_daminfo.json")):
            os.remove(os.path.join(out, dam_path + "_daminfo.json"))
        dst_ds = drv.CreateDataSource(os.path.join(out, dam_path + "_daminfo.json"))
        dst_layer = dst_ds.CreateLayer("", srs=carto, geom_type=ogr.wkbPoint)
        field_defn = ogr.FieldDefn("name", ogr.OFTString)
        field_defnID = ogr.FieldDefn("ID", ogr.OFTInteger64)
        field_defnelev = ogr.FieldDefn("elev", ogr.OFTString)
        field_defndam = ogr.FieldDefn("damname", ogr.OFTString)
        dst_layer.CreateField(field_defn)
        dst_layer.CreateField(field_defnID)
        dst_layer.CreateField(field_defnelev)
        dst_layer.CreateField(field_defndam)

        rad_l = []
        alt_l = []
        for r in range(pdbradius, 1, -1 * pdbstep):
            extract_roi_parameters_ext_l = ExtractROIParam(
                mode="radius",
                mode_radius_r=r,
                mode_radius_unitr="phy",
                mode_radius_unitc="phy",
                mode_radius_cx=dam.GetX(),
                mode_radius_cy=dam.GetY(),
            )

            ext_l, profile_ext_l = extract_roi(
                rio.open(dem), extract_roi_parameters_ext_l
            )

            np_ext_l = ext_l.reshape(ext_l.shape[1], ext_l.shape[2])
            ext_l_alt = np.amin(np_ext_l)

            rad_l.append(r)
            alt_l.append(ext_l_alt)

        d = nderiv(alt_l, rad_l)

        fig, axs = plt.subplots(2)
        axs[0].plot(rad_l, alt_l, "r")
        axs[0].set(xlabel="Search Area to the Dam (m)", ylabel="Minimum Elevation")
        axs[0].label_outer()
        axs[1].plot(rad_l, abs(d), "b")
        axs[1].set(xlabel="Search Area to the Dam (m)", ylabel="d(Minimum Elevation)")
        axs[1].label_outer()

        found_pdb = False
        rad_pdb = 0
        alt_pdb = 0

        # TODO: @param 0.15
        for i_r, i_a, i_d in zip(rad_l, alt_l, d):
            #  print(i_r, i_a, i_d)
            if abs(i_d) > 0.15:
                found_pdb = True
                rad_pdb = i_r
                alt_pdb = i_a
                fig.suptitle("PDB profile (1st pass detection)")
                axs[0].plot(rad_pdb, alt_pdb, "x")
                axs[1].plot(rad_pdb, abs(i_d), "x")
                break

        if found_pdb is True:
            logging.debug("@radius= " + str(rad_pdb) + "m: local min = " + str(alt_pdb))
        else:
            for i_r, i_a, i_d in zip(rad_l, alt_l, d):
                #  print(i_r, i_a, i_d)
                if abs(i_d) > 0.10:
                    found_pdb = True
                    rad_pdb = i_r
                    alt_pdb = i_a
                    fig.suptitle("PDB profile (2nd pass detection)")
                    axs[0].plot(rad_pdb, alt_pdb, "x")
                    axs[1].plot(rad_pdb, abs(i_d), "x")
                    logging.warning(
                        "PDB found during 2nd pass - It may not be reliable"
                    )
                    logging.debug(
                        "@radius= " + str(rad_pdb) + "m: local min = " + str(alt_pdb)
                    )
                    break
                else:
                    logging.debug(
                        "@radius= "
                        + str(rad_pdb)
                        + "m: elev = "
                        + str(alt_pdb)
                        + "m: delev = "
                        + str(i_d)
                    )

        fig.savefig(os.path.join(tmp, "pdb_profile.png"))

        if found_pdb is False:
            logging.error("404 - PDB not Found")
            sys.exit("PDB search failed for dam " + dam_name + ". Aborting.")

        # Retrieve PDB coordinates
        extract_roi_parameters_ext = ExtractROIParam(
            mode="radius",
            mode_radius_r=rad_pdb,
            mode_radius_unitr="phy",
            mode_radius_unitc="phy",
            mode_radius_cx=dam.GetX(),
            mode_radius_cy=dam.GetY(),
        )

        ext, profile_ext = extract_roi(rio.open(dem), extract_roi_parameters_ext)
        save_image(ext, profile_ext, os.path.join(out, "dem_pdb.tif"))
        np_ext = ext.reshape(ext.shape[1], ext.shape[2])

        indices = np.where(np_ext == [alt_pdb])
        # TODO: if multiple pdb detected
        if len(indices[0]) > 1:
            logging.warning(
                "Absolute minimum is not unique on the current area!["
                + str(len(indices[0]))
                + "]"
            )


        ds = gdal.Open(os.path.join(out, "dem_pdb.tif"))


        posX, posY = coord(indices[1][0], indices[0][0], ds)
        pX, pY = pixel(posX, posY, ds)

        pdb = ogr.Geometry(ogr.wkbPoint)
        pdb.AddPoint(float(posX), float(posY))
        pdb.Transform(cartotogeo)
        pdblat = pdb.GetX()
        pdblon = pdb.GetY()

        pdbalt = float(
            os.popen(
                'gdallocationinfo -valonly -wgs84 "%s" %s %s' % (dem, pdblon, pdblat)
            ).read()
        )
        logging.debug("Coordinates (pixel): " + str(pX) + " - " + str(pY))
        logging.debug("Coordinates (carto): " + str(posX) + " - " + str(posY))
        logging.debug(
            "Coordinates (latlon): " + str(pdb.GetX()) + " - " + str(pdb.GetY())
        )
        logging.info(
            "PDB detected: "
            + dam_name
            + "(id:"
            + str(id)
            + ") [pdbLat: "
            + str(pdblat)
            + ", pdbLon: "
            + str(pdblon)
            + ", pdbAlt: "
            + str(pdbalt)
            + "]"
        )

        wkt = "POINT ( %f %f )" % (float(posX), float(posY))
        feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
        p = ogr.CreateGeometryFromWkt(wkt)
        feat.SetGeometryDirectly(p)
        feat.SetField("name", "PDB")
        feat.SetField("elev", pdbalt)

        feat.SetField("damname", dam_path)
        feat.SetField("ID", dam_id)
        dst_layer.CreateFeature(feat)
        feat.Destroy()

        # Point on dam
        wkt = "POINT ( %f %f )" % (float(dam.GetX()), float(dam.GetY()))
        feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
        p = ogr.CreateGeometryFromWkt(wkt)
        feat.SetGeometryDirectly(p)
        feat.SetField("name", "Dam")
        feat.SetField("elev", str(bml_alt))
        feat.SetField("damname", dam_path)
        feat.SetField("ID", dam_id)
        dst_layer.CreateFeature(feat)
        feat.Destroy()

        # Point inside water body from DB
        in_w = ogr.Geometry(ogr.wkbPoint)
        in_w.AddPoint(clat_in, clon_in)
        in_w.Transform(geotocarto)
        logging.debug(
            "Coordinates Carto Point Inside Water Body: "
            + str(in_w.GetX())
            + " - "
            + str(in_w.GetY())
        )

        wkt = "POINT ( %f %f )" % (float(in_w.GetX()), float(in_w.GetY()))
        feat = ogr.Feature(feature_def=dst_layer.GetLayerDefn())
        p = ogr.CreateGeometryFromWkt(wkt)
        feat.SetGeometryDirectly(p)
        feat.SetField("name", "Insider")
        feat.SetField("elev", "")
        feat.SetField("damname", dam_path)
        feat.SetField("ID", dam_id)
        dst_layer.CreateFeature(feat)
        feat.Destroy()
    # END
    else:
        # User-defined daminfo.json was provided
        logging.info("User-defined daminfo.json has been provided.")
        with open(info) as i:
            jsi = json.load(i)
        for feature in jsi["features"]:
            if feature["properties"]["name"] == "Dam":
                logging.debug(feature)
                bml_alt = float(feature["properties"]["elev"])
                logging.info("Dam elevation extracted from user provided daminfo.json")
                logging.info("Extracted dam elevation= " + str(bml_alt) + "m")
                logging.info(
                    "Target Elevation for cutline search= " + str(targetelev) + "m"
                )

            if feature["properties"]["name"] == "PDB":
                logging.debug(feature)
                pdbin = shape(feature["geometry"])

                pdb = ogr.Geometry(ogr.wkbPoint)
                pdb.AddPoint(float(pdbin.x), float(pdbin.y))
                pdb.Transform(cartotogeo)
                pdblat = pdb.GetX()
                pdblon = pdb.GetY()
                pdbalt = float(
                    os.popen(
                        'gdallocationinfo -valonly -wgs84 "%s" %s %s'
                        % (dem, pdblon, pdblat)
                    ).read()
                )

                logging.debug(
                    "Coordinates (carto): " + str(pdbin.x) + " - " + str(pdbin.y)
                )
                logging.debug(
                    "Coordinates (latlon): " + str(pdblat) + " - " + str(pdblon)
                )
                logging.info(
                    "PDB extracted from user provided daminfo.json: "
                    + " [pdbLat: "
                    + str(pdblat)
                    + ", pdbLon: "
                    + str(pdblon)
                    + ", pdbAlt: "
                    + str(pdbalt)
                    + "]"
                )

    # Search for Dam Line
    # TODO: fix confusion between radius and pdbradius needed here
    extract_roi_parameters_ext_r = ExtractROIParam(
        mode="radius",
        mode_radius_r=radius,
        mode_radius_unitr="phy",
        mode_radius_unitc="phy",
        mode_radius_cx=dam.GetX(),
        mode_radius_cy=dam.GetY(),
    )

    ext_r, profile_ext_r = extract_roi(rio.open(dem), extract_roi_parameters_ext_r)
    save_image(
        ext_r,
        profile_ext_r,
        os.path.join(tmp, "extract@" + str(radius) + "mFromDam.tif"),
    )

    r_ds = gdal.Open(
        os.path.join(tmp, "extract@" + str(radius) + "mFromDam.tif"),
        gdal.GA_ReadOnly,
    )
    r_x, r_y = pixel(dam.GetX(), dam.GetY(), r_ds)


    np_r = ext_r.reshape(ext_r.shape[1], ext_r.shape[2])

    (image_size_y, image_size_x) = np_r.shape
    # Disk definition:
    (center_x, center_y) = (r_x, r_y)
    radius = 10

    x_grid, y_grid = np.meshgrid(np.arange(image_size_x), np.arange(image_size_y))

    prev = 0
    prevprev = 0
    prevpoint1 = ogr.Geometry(ogr.wkbPoint)
    prevpoint1.AddPoint(dam.GetX(), dam.GetY())
    prev1alt = 0
    prevpoint2 = ogr.Geometry(ogr.wkbPoint)
    prevpoint2.AddPoint(dam.GetX(), dam.GetY())
    multiline = ogr.Geometry(ogr.wkbMultiLineString)
    prev2alt = 0
    # TODO: @param
    step_lc = 2
    lc_first_it = True
    stop_side_1 = False
    stop_side_2 = False
    # NB: range is define pixel-wise
    # TODO: @param 500 -> pdbradius
    for radius in range(5, 500, step_lc):
        # Array of booleans with the disk shape
        circle = np.logical_and(
            ((x_grid - center_x) ** 2 + (y_grid - center_y) ** 2) <= radius**2,
            ((x_grid - center_x) ** 2 + (y_grid - center_y) ** 2) > (radius - 2) ** 2,
        )
        disk_out = ((x_grid - center_x) ** 2 + (y_grid - center_y) ** 2) <= radius**2
        disk_in = ((x_grid - center_x) ** 2 + (y_grid - center_y) ** 2) <= (
            radius - 1
        ) ** 2
        points_on_circle = np_r[circle]

        if len(points_on_circle) > prevprev:
            prevprev = prev
            prev = len(points_on_circle)

            masked_tmp = np.where(disk_out, np_r, 0)
            masked = np.where(~disk_in, masked_tmp, 0)
            # im_masked = Image.fromarray(masked)
            # im_masked.save(os.path.join(tmp, "circle@"+str(radius)+".tif"))

            l_indices = np.where(masked == [np.amax(masked)])
            l_posX, l_posY = coord(l_indices[1][0], l_indices[0][0], r_ds)
            l_pX, l_pY = pixel(l_posX, l_posY, r_ds)

            currpoint = ogr.Geometry(ogr.wkbPoint)
            currpoint.AddPoint(l_posX, l_posY)

            # TODO: if multiple max detected
            if len(l_indices[0]) > 1:
                logging.warning(
                    "Absolute maximum is not unique on the current circle!["
                    + str(len(l_indices[0]))
                    + "]"
                )

            # Add Circle absolute max to json
            l_wkt = "POINT ( %f %f )" % (float(l_posX), float(l_posY))
            l_feat = ogr.Feature(feature_def=dbg_layer.GetLayerDefn())
            l_p = ogr.CreateGeometryFromWkt(l_wkt)
            l_feat.SetGeometryDirectly(l_p)
            l_feat.SetField("name", str(radius) + "/1")
            dbg_layer.CreateFeature(l_feat)
            l_feat.Destroy()

            #  logging.debug("First Point - Coordinates (pixel): " + str(l_pX)+" - "+str(l_pY))
            #  logging.debug("First Point - Coordinates (carto): " + str(l_posX)+" - "+str(l_posY))

            distance1 = currpoint.Distance(prevpoint1) / r_ds.GetGeoTransform()[1]
            distance2 = currpoint.Distance(prevpoint2) / r_ds.GetGeoTransform()[1]
            if (
                (distance1 > radius * maxdist)
                and (distance2 > radius * maxdist)
                and (lc_first_it is False)
            ):
                logging.debug(
                    "New detected point "
                    + "(first point of the current iteration) "
                    + "too distant from previous one!"
                )
                #  logging.debug("radius    (px): " + str(radius))
                #  logging.debug("maxdist   (px): " + str(radius*maxdist))
                #  logging.debug("distance1 (px): " + str(distance1))
                #  logging.debug("distance2 (px): " + str(distance2))
                # Search a local maxima:
                px, py = pixel(prevpoint1.GetX(), prevpoint1.GetY(), r_ds)
                if distance1 > distance2:
                    px, py = pixel(prevpoint2.GetX(), prevpoint2.GetY(), r_ds)
                force_local = ((x_grid - px) ** 2 + (y_grid - py) ** 2) <= (
                    2 * step_lc
                ) ** 2

                force_local_masked = np.where(force_local, masked, 0)
                # im_fmasked = Image.fromarray(force_local_masked)
                f_indices = np.where(
                    force_local_masked == [np.amax(force_local_masked)]
                )
                f_posX, f_posY = coord(f_indices[1][0], f_indices[0][0], r_ds)
                f_pX, f_pY = pixel(f_posX, f_posY, r_ds)

                # TODO: if multiple max detected
                if len(f_indices[0]) > 1:
                    logging.warning(
                        "Absolute maximum is not unique on the current local mask!["
                        + str(len(f_indices[0]))
                        + "]"
                    )
                # TODO: find a better way to detect this problematic case:
                if len(f_indices[0]) < 42:  # Detecting when masked area is only zeros
                    currpoint = ogr.Geometry(ogr.wkbPoint)
                    currpoint.AddPoint(f_posX, f_posY)
                    # Add local absolute max to json
                    l_wkt = "POINT ( %f %f )" % (float(f_posX), float(f_posY))
                    l_feat = ogr.Feature(feature_def=dbg_layer.GetLayerDefn())
                    l_p = ogr.CreateGeometryFromWkt(l_wkt)
                    l_feat.SetGeometryDirectly(l_p)
                    l_feat.SetField("name", str(radius) + "/1/alt")
                    dbg_layer.CreateFeature(l_feat)
                    l_feat.Destroy()
                    logging.debug(
                        "New detected point "
                        + "(first point of the current iteration) "
                        + "corrected by local restricted search!"
                    )
                else:
                    logging.warning(
                        "Too distant detected point "
                        + "(first point of the current iteration) "
                        + "NOT corrected by local restricted search!"
                    )

            # Search for opposite relative max: mask half of the circle
            done_half = (
                (x_grid - l_indices[1][0]) ** 2 + (y_grid - l_indices[0][0]) ** 2
            ) <= (radius * mradius) ** 2
            half_masked = np.where(~done_half, masked, 0)
            # im_hmasked = Image.fromarray(half_masked)
            # im_hmasked.save(os.path.join(tmp, "half_circle@"+str(radius)+".tif"))

            l_indices = np.where(half_masked == [np.amax(half_masked)])
            l_posX, l_posY = coord(l_indices[1][0], l_indices[0][0], r_ds)
            l_pX, l_pY = pixel(l_posX, l_posY, r_ds)

            nextpoint = ogr.Geometry(ogr.wkbPoint)
            nextpoint.AddPoint(l_posX, l_posY)

            # TODO: if multiple max detected
            if len(l_indices[0]) > 1:
                logging.warning(
                    "Absolute maximum is not unique on the current circle!["
                    + str(len(l_indices[0]))
                    + "]"
                )

            # Add Circle absolute max to json
            l_wkt = "POINT ( %f %f )" % (float(l_posX), float(l_posY))
            l_feat = ogr.Feature(feature_def=dbg_layer.GetLayerDefn())
            l_p = ogr.CreateGeometryFromWkt(l_wkt)
            l_feat.SetGeometryDirectly(l_p)
            l_feat.SetField("name", str(radius) + "/2")
            dbg_layer.CreateFeature(l_feat)
            l_feat.Destroy()

            #  logging.debug("Second Point - Coordinates (pixel): " + str(l_pX)+" - "+str(l_pY))
            #  logging.debug("Second Point - Coordinates (carto): " + str(l_posX)+" - "+str(l_posY))

            distance3 = nextpoint.Distance(prevpoint1) / r_ds.GetGeoTransform()[1]
            distance4 = nextpoint.Distance(prevpoint2) / r_ds.GetGeoTransform()[1]
            if (
                (distance3 > radius * maxdist)
                and (distance4 > radius * maxdist)
                and (lc_first_it is False)
            ):
                logging.debug(
                    "New detected point "
                    + "(second point of the current iteration) "
                    + "too distant from previous one!"
                )
                #  logging.debug("radius    (px): " + str(radius))
                #  logging.debug("maxdist   (px): " + str(radius*maxdist))
                #  logging.debug("distance3 (px): " + str(distance3))
                #  logging.debug("distance4 (px): " + str(distance4))
                # Search a local maxima:
                px, py = pixel(prevpoint1.GetX(), prevpoint1.GetY(), r_ds)
                if distance3 > distance4:
                    px, py = pixel(prevpoint2.GetX(), prevpoint2.GetY(), r_ds)
                force_local = ((x_grid - px) ** 2 + (y_grid - py) ** 2) <= (
                    2 * step_lc
                ) ** 2

                force_local_masked = np.where(force_local, half_masked, 0)
                # im_fmasked = Image.fromarray(force_local_masked)
                f_indices = np.where(
                    force_local_masked == [np.amax(force_local_masked)]
                )
                f_posX, f_posY = coord(f_indices[1][0], f_indices[0][0], r_ds)
                f_pX, f_pY = pixel(f_posX, f_posY, r_ds)

                # TODO: if multiple max detected
                if len(f_indices[0]) > 1:
                    logging.warning(
                        "Absolute maximum is not unique on the current local mask!["
                        + str(len(f_indices[0]))
                        + "]"
                    )
                # TODO: find a better way to detect this problematic case:
                if len(f_indices[0]) < 42:  # Detecting when masked area is only zeros
                    nextpoint = ogr.Geometry(ogr.wkbPoint)
                    nextpoint.AddPoint(f_posX, f_posY)
                    # Add local absolute max to json
                    f_wkt = "POINT ( %f %f )" % (float(f_posX), float(f_posY))
                    f_feat = ogr.Feature(feature_def=dbg_layer.GetLayerDefn())
                    f_p = ogr.CreateGeometryFromWkt(f_wkt)
                    f_feat.SetGeometryDirectly(f_p)
                    f_feat.SetField("name", str(radius) + "/2/alt")
                    dbg_layer.CreateFeature(f_feat)
                    f_feat.Destroy()
                    logging.debug(
                        "New detected point "
                        + "(second point of the current iteration) "
                        + "corrected by local restricted search!"
                    )
                else:
                    logging.warning(
                        "Too distant detected point "
                        + "(second point of the current iteration) "
                        + "NOT corrected by local restricted search!"
                    )

            # Check if under max elev
            prev1geo = ogr.Geometry(ogr.wkbPoint)
            prev1geo.AddPoint(float(prevpoint1.GetX()), float(prevpoint1.GetY()))
            prev1geo.Transform(cartotogeo)
            prev1alt = float(
                os.popen(
                    'gdallocationinfo -valonly -wgs84 "%s" %s %s'
                    % (dem, prev1geo.GetY(), prev1geo.GetX())
                ).read()
            )

            if (prev1alt > targetelev) and (stop_side_1 is False):
                stop_side_1 = True
                logging.info(
                    "Stop cutline search on side #1 "
                    + "[prevalt: "
                    + str(prev1alt)
                    + " ; targeted elevation: "
                    + str(targetelev)
                    + "]"
                )

            prev2geo = ogr.Geometry(ogr.wkbPoint)
            prev2geo.AddPoint(prevpoint2.GetX(), prevpoint2.GetY())
            prev2geo.Transform(cartotogeo)
            prev2alt = float(
                os.popen(
                    'gdallocationinfo -valonly -wgs84 "%s" %s %s'
                    % (dem, str(prev2geo.GetY()), str(prev2geo.GetX()))
                ).read()
            )

            if (prev2alt > targetelev) and (stop_side_2 is False):
                stop_side_2 = True
                logging.info(
                    "Stop cutline search on side #2 "
                    + "[prevalt: "
                    + str(prev2alt)
                    + " ; targeted elevation: "
                    + str(targetelev)
                    + "]"
                )

            if (stop_side_1 is True) and (stop_side_2 is True):
                logging.info(
                    "Cutline over target elevation ("
                    + str(targetelev)
                    + " m) on both sides - Stopping Search."
                )
                break

            if distance1 <= distance2:
                line1 = ogr.Geometry(ogr.wkbLineString)
                line1.AddPoint(prevpoint1.GetX(), prevpoint1.GetY())
                line1.AddPoint(currpoint.GetX(), currpoint.GetY())
                if stop_side_1 is False:
                    multiline.AddGeometry(line1)
                prevpoint1 = currpoint

                line2 = ogr.Geometry(ogr.wkbLineString)
                line2.AddPoint(prevpoint2.GetX(), prevpoint2.GetY())
                line2.AddPoint(nextpoint.GetX(), nextpoint.GetY())
                if stop_side_2 is False:
                    multiline.AddGeometry(line2)
                prevpoint2 = nextpoint

            else:
                #  logging.debug("Second")
                line1 = ogr.Geometry(ogr.wkbLineString)
                line1.AddPoint(prevpoint2.GetX(), prevpoint2.GetY())
                line1.AddPoint(currpoint.GetX(), currpoint.GetY())
                if stop_side_2 is False:
                    multiline.AddGeometry(line1)
                prevpoint2 = currpoint

                line2 = ogr.Geometry(ogr.wkbLineString)
                line2.AddPoint(prevpoint1.GetX(), prevpoint1.GetY())
                line2.AddPoint(nextpoint.GetX(), nextpoint.GetY())
                if stop_side_1 is False:
                    multiline.AddGeometry(line2)
                prevpoint1 = nextpoint

            # Activate Correction after first iteration
            lc_first_it = False

        else:
            logging.debug(
                "Search area outside of ROI @radius= " + str(radius) + " pixels"
            )
            break

    if (stop_side_1 is False) or (stop_side_2 is False):
        logging.warning(
            "Target elevation for cutline extremities ("
            + str(targetelev)
            + " m) not reached! ["
            + str(prev1alt)
            + " m; "
            + str(prev2alt)
            + " m]"
        )

    # Export line to line.json
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(multiline)
    outLayer.CreateFeature(outFeature)
    t1_stop = perf_counter()
    logging.info(f"Elapsed time:{t1_stop}s, {t1_start}s")

    logging.info(f"Elapsed time during the whole program in s : {t1_stop-t1_start}s")
    print("End SZI_FROM_CONTOURLINE")


def find_pdb_and_cutline_parameters():
    """Define parameters."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--infile", help="Input file")
    parser.add_argument("--id", help="Dam ID")
    parser.add_argument("--id_db", help="Dam id field in database")
    parser.add_argument("-w", "--watermap", help="Input water map file")
    parser.add_argument("-d", "--dem", help="Input DEM")
    parser.add_argument("-r", "--radius", help="Extract radius (m)", default=None)
    parser.add_argument(
        "-s", "--pdbstep", type=int, default=5, help="Sampling step for pdb search"
    )
    parser.add_argument(
        "-p", "--pdbradius", type=int, default=500, help="PDB Search Radius"
    )
    parser.add_argument(
        "--mradius",
        type=float,
        default=1.7,
        help="Masking radius for second point search wrt internal loop radius",
    )
    parser.add_argument(
        "--maxdist",
        type=float,
        default=0.3,
        help="Max distance between two concecutive points in the cutline wrt internal loop radius",
    )
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
    parser.add_argument("--info", help="Optional user-defined daminfo.json file")

    parser.add_argument("-t", "--tmp", help="Temporary directory")
    parser.add_argument("-o", "--out", help="Output directory")
    parser.add_argument("--debug", action="store_true", help="Activate Debug Mode")
    return parser


def main():
    """Cli for find_pdb_and_cutline."""
    parser = find_pdb_and_cutline_parameters()
    args = parser.parse_args()
    find_pdb_and_cutline(
        args.infile,
        args.id,
        args.id_db,
        args.watermap,
        args.dem,
        args.pdbstep,
        args.pdbradius,
        args.mradius,
        args.maxdist,
        args.elevoffset,
        args.elevsampling,
        args.info,
        args.tmp,
        args.out,
        args.radius,
        args.debug,
    )


if __name__ == "__main__":
    sys.exit(main())
