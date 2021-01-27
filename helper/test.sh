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
DAMNAME=${1:-"Agly"}
DAM=${DAMNAME// /_}
RADIUS=${2:-7500}


module load otb/7.2-python3.7.2
[ -z "$OTB_MAX_RAM_HINT" ] && export OTB_MAX_RAM_HINT=4000
# [ -z "$OTB_LOGGER_LEVEL" ] && export OTB_LOGGER_LEVEL=CRITICAL
export OTB_LOGGER_LEVEL=CRITICAL
echo "OTB_MAX_RAM_HINT: $OTB_MAX_RAM_HINT"
echo "OTB_LOGGER_LEVEL: $OTB_LOGGER_LEVEL"


cd $SRC_DIR
[ -d "$ROOT_DIR/${DAM}_${RADIUS}/tmp" ] || mkdir -p "$ROOT_DIR/${DAM}_${RADIUS}/tmp"


#Lat: 42.7461
#Lon: 2.58972
# gdallocationinfo -valonly -wgs84 ../data/dem/dem.vrt 2.58972 42.7461

# Agly
# Lat: 42.7461
# Lon: 2.58972
# Alt: 184.119995117188
# Coordinates: 466421.36695527623 - 4732701.863316036
# Bottom Alt: 123.55


if [ -f "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif" ] \
  && [ -f "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" ] ; then
  echo "Extracts already available --> Skipping area_mapping."
else
  python3 area_mapping.py \
    --name     "${DAMNAME}" \
    --infile   "../data/DB_Barrages_Fixed_v3/DB_Barrages_Fixed.shp" \
    --watermap "../data/wmap/wmap.vrt" \
    --dem      "../data/dem/dem.vrt" \
    --radius   5000 \
    --out      "$ROOT_DIR/${DAM}_${RADIUS}" \
    # --debug

fi

python3 szi_from_contourline.py \
  --name       "${DAMNAME}" \
  --infile     "../data/DB_Barrages_Fixed_v3/DB_Barrages_Fixed.shp" \
  --watermap   "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
  --dem        "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
  --radius     5000 \
  --pdbstep    5 \
  --pdbradius  500 \
  --elevoffset 100 \
  --tmp        "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
  --out        "$ROOT_DIR/${DAM}_${RADIUS}" \
  --debug

exit

if [ -f "$ROOT_DIR/${DAM}_${RADIUS}/contour.shp" ] ; then
  echo "Contours already available --> Skipping gdal_contour."
else
  gdal_contour -p -i 10.0 \
    "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-Agly.tif" \
    "$ROOT_DIR/${DAM}_${RADIUS}/contour.shp"
fi

if [ -f "$ROOT_DIR/${DAM}_${RADIUS}/szi_wmap.dat" ] ; then
  echo "S(Zi)s already available --> Skipping szi_from_watermap."
else
  python3 szi_from_watermap.py \
    --watermap "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
    --dem      "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif" \
    --tmp      "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
    --zmin     123.55 \
    --zmax     184.11 \
    --step     2 \
    --outfile  "$ROOT_DIR/${DAM}_${RADIUS}/szi_wmap.dat" \
    --debug

fi


if [ -f "$ROOT_DIR/${DAM}_${RADIUS}/szi_wmap.dat" ] ; then
  python3 szi_to_model.py \
    --infile   "$ROOT_DIR/${DAM}_${RADIUS}/szi_wmap.dat" \
    --outfile  "$ROOT_DIR/${DAM}_${RADIUS}/szi_wmap.png" \
    --debug

fi
