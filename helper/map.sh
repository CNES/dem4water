#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2020 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2020
# ----

SRC_DIR="/home/ad/briciera/dem4water/dem4water"
ROOT_DIR="/home/ad/briciera/scratch/HSV"
DAM=${1:-"Agly"}
RADIUS=7500


[ -d "$ROOT_DIR/${DAM}_${RADIUS}/tmp" ] || mkdir -p "$ROOT_DIR/${DAM}_${RADIUS}/tmp"

cd $SRC_DIR || exit
export OTB_MAX_RAM_HINT=4000
export OTB_LOGGER_LEVEL=CRITICAL
module load otb/7.2-python3.7.2


#Lat: 42.7461
#Lon: 2.58972
# gdallocationinfo -valonly -wgs84 ../data/dem/dem.vrt 2.58972 42.7461

# Agly
# Lat: 42.7461
# Lon: 2.58972
# Alt: 184.119995117188
# Coordinates: 466421.36695527623 - 4732701.863316036
# Bottom Alt: 123.55


python3 area_mapping.py \
  --infile   "../data/synth_names.shp" \
  --watermap "../data/T31TDH/all_cumul.tif" \
  --dem      "../data/dem/dem.vrt" \
  --name     "$DAM" \
  --radius   "$RADIUS" \
  --out      "$ROOT_DIR/${DAM}_${RADIUS}" \
  --debug

gdal_contour -p -i 10.0 \
  "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-Agly.tif" \
  "$ROOT_DIR/${DAM}_${RADIUS}/contour.shp"
