#!/usr/bin/python
#  -*- coding: utf-8 -*-
"""
:author: Benjamin TARDY
:organization: CNES
:copyright: 2021 CNES. All rights reserved.
:license: see LICENSE file
:created: 2022
"""

"""
Avant d'utiliser ce script:
- assurer que le csv n'est pas corrompu
- assurer que le csv a bien été enregistré en utf-8
- si besoin convertir le fichier excel en csv
- assurer les projections, si besoin éditer le script
- les colonnes sont censées être fixées pour dem4water
"""
import argparse
import unicodedata

import geopandas as gpd
import pandas as pd


def normalize_list_str(list_str):
    out_list = []
    for name in list_str:
        if name is None:
            out_list.append(None)
        elif not type(name) == str:
            out_list.append(None)
        else:
            tmp_name = unicodedata.normalize("NFKD", name)
            tmp_name = tmp_name.encode("ascii", "ignore")
            tmp_name = tmp_name.decode("utf-8")
            tmp_name = tmp_name.replace("(", "").replace(")", "").replace("'", " ")
            out_list.append(tmp_name)

    return out_list


def normalize_list_float(list_float):
    out_list = []
    for value in list_float:
        if value is None:
            out_list.append(None)
        else:
            if type(value) == str:
                value = value.replace(",", ".")
                out_list.append(value)
            elif type(value) == float:
                out_list.append(value)
            else:
                try:
                    value = float(value)
                    out_list.append(value)
                except ValueError:
                    # The value is unexpected replace by None for futher process
                    out_list.append(None)


if __name__ == "__main__":
    """
    Usage : python generate_list_from_DB.py dams_file dam_id dam_name ouput_list
    """

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("dams_file", type=str, help="Dams database shp file")
    parser.add_argument("dams_csv", type=str, help="DAM info csv file")
    parser.add_argument("reservoirs_file", type=str, help="Reservoirs pôlygon files")
    parser.add_argument("output_file", type=str, help="Corrected DB")

    args = parser.parse_args()

    list_of_fields_str = ["RES_NAME", "DAM_NAME", "MAIN_USE", "NEAR_CITY"]
    list_of_fields_float = ["DAM_LVL_M"]
    gdf = gpd.GeoDataFrame().from_file(args.dams_file)
    df = pd.read_csv(args.dams_csv)
    gdf2 = gpd.GeoDataFrame().from_file(args.reservoirs_file)

    for field in list_of_fields_str:
        list_values = list(df[field])
        corrected_values = normalize_list_str(list_values)
        gdf[field] = corrected_values

    for field in list_of_fields_float:
        list_values = list(gdf[field])
        corrected_values = normalize_list_float(list_values)
        gdf[field] = corrected_values

    gdf = gdf.to_crs("EPSG:2154")
    out_gdf = gpd.sjoin_nearest(gdf2, gdf, how="left")
    out_gdf = out_gdf.drop(["index_right"], axis="columns")
    out_gdf = out_gdf.to_crs("EPSG:4326")
    out_gdf.to_file(args.output_file)
