#!/usr/bin/python
#  -*- coding: utf-8 -*-
'''
:author: Gaël Nicolas
:organization: CNES
:copyright: 2021 CNES. All rights reserved.
:license: see LICENSE file
:created: 2021
'''

import os
import sys
import csv
import argparse
import logging


if __name__ == "__main__":
    """
    Usage : python generate_list_from_DB.py dams_csv id_field id_name ouput_list
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('dams_csv',type=str, help='Dams database csv file') # /work/OT/siaa/Work/SCO_StockWater/dams_database/Andalousie/andalousie_selected_max_extent_v7.csv
    parser.add_argument('id_field',type=str, help='id_field list') # ID_SWOT
    parser.add_argument('id_name',type=str, help='id_name list')   # DAM_NAME
    parser.add_argument('ouput_list',type=str, help='Dams list')

    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    # Read anw write id and dams names
    with open(args.dams_csv, 'r') as csv_file:
        reader = csv.DictReader(csv_file, delimiter =";")
        with open(args.ouput_list, 'w') as output_file:

            for row in reader:
                output_file.write("%s,%s \n" % (row[args.id_field], row[args.id_name]))

