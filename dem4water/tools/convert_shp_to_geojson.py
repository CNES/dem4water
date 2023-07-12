#!/usr/bin/python
#  -*- coding: utf-8 -*-
"""
Merge a csv and a geometry file to create a geojson database.

Avant d'utiliser ce script:
- assurer que le csv n'est pas corrompu
- assurer que le csv a bien été enregistré en utf-8
- si besoin convertir le fichier excel en csv
- assurer les projections, si besoin éditer le script
- les colonnes sont censées être fixées pour dem4water
"""
import argparse
import sys
import unicodedata

import geopandas as gpd
import numpy as np
import pandas as pd


def normalize_list_str(list_str):
    """Remove special letters from names."""
    out_list = []
    for name in list_str:
        if name is None:
            out_list.append(None)
        elif not isinstance(name, str):
            out_list.append(None)
        else:
            tmp_name = unicodedata.normalize("NFKD", name)
            tmp_name = tmp_name.encode("ascii", "ignore")
            tmp_name = tmp_name.decode("utf-8")
            tmp_name = tmp_name.replace("(", "").replace(")", "").replace("'", " ")
            out_list.append(tmp_name)

    return out_list


def normalize_list_float(list_float):
    """Ensure float are correctly encoded."""
    out_list = []
    for value in list_float:
        if value is None:
            out_list.append(np.nan)
        else:
            if isinstance(value, str):
                value = value.replace(",", ".")
                # out_list.append(value)
            if isinstance(value, float):
                out_list.append(value)
            else:
                try:
                    value = float(value)
                    out_list.append(value)
                except ValueError:
                    # The value is unexpected replace by None for futher process
                    out_list.append(np.nan)
    return out_list


def main():
    """Launch function."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("dams_file", type=str, help="Dams database shp file")
    parser.add_argument("dams_csv", type=str, help="DAM info csv file")
    parser.add_argument("reservoirs_file", type=str, help="Reservoirs pôlygon files")
    parser.add_argument("output_file", type=str, help="Corrected DB")

    args = parser.parse_args()

    list_of_fields_str = ["RES_NAME", "DAM_NAME", "MAIN_USE", "NEAR_CITY"]
    list_of_fields_float = ["DAM_LVL_M", "AREA_HA"]
    gdf = gpd.GeoDataFrame().from_file(args.dams_file)
    df_csv = pd.read_csv(args.dams_csv)
    gdf2 = gpd.GeoDataFrame().from_file(args.reservoirs_file)

    for field in list_of_fields_str:
        list_values = list(df_csv[field])
        corrected_values = normalize_list_str(list_values)
        gdf[field] = corrected_values

    for field in list_of_fields_float:
        list_values = list(gdf[field])
        corrected_values = normalize_list_float(list_values)
        gdf[field] = corrected_values
        # input(gdf[field])
    gdf = gdf.to_crs("EPSG:2154")
    out_gdf = gpd.sjoin_nearest(gdf2, gdf, how="left")
    out_gdf = out_gdf.drop(["index_right"], axis="columns")
    out_gdf = out_gdf.to_crs("EPSG:4326")
    out_gdf.to_file(args.output_file)


if __name__ == "__main__":

    sys.exit(main())
