#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2020 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2020
# ----

cd ~/dem4water/dem4water
ROOT_DIR="/home/ad/briciera/scratch/HSV"
DAM=${1:-"Agly"}

[ -d "$ROOT_DIR/$DAM/tmp" ] || mkdir -p "$ROOT_DIR/$DAM/tmp"

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


if [ -f "$ROOT_DIR/$DAM/dem_extract-$DAM.tif" ] \
  && [ -f "$ROOT_DIR/$DAM/wmap_extract-$DAM.tif" ] ; then
  echo "Extracts already available --> Skipping area_mapping."
else
  python3 area_mapping.py \
    --infile   "../data/synth_names.shp" \
    --watermap "../data/T31TDH/all_cumul.tif" \
    --dem      "../data/dem/dem.vrt" \
    --name     "$DAM" \
    --radius   2000 \
    --out      "$ROOT_DIR/$DAM" # --loglevel "DEBUG"
fi


python3 szi_from_watermap.py \
  --watermap "$ROOT_DIR/$DAM/wmap_extract-$DAM.tif" \
  --dem      "$ROOT_DIR/$DAM/dem_extract-$DAM.tif" \
  --tmp      "$ROOT_DIR/$DAM/tmp" \
  --zmin     123.55 \
  --zmax     184.11 \
  --step     5 \
  --outfile  "$ROOT_DIR/$DAM/szi_wmap.txt" # --loglevel "DEBUG"
