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
EXTR_DIR="/home/ad/briciera/scratch/HSV/Extracts"
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

[ -d "$EXTR_DIR/${DAM}_${RADIUS}" ] || mkdir -p "$EXTR_DIR/${DAM}_${RADIUS}"
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


if [ -f "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif" ] \
  && [ -f "$EXTR_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" ] ; then
  echo "Extracts already available --> Skipping area_mapping."
else
  python3 area_mapping.py \
    --name     "${DAMNAME}" \
    --infile   "/home/ad/briciera/dem4water/dem4water/data/DB_v5_KL.geojson" \
    --watermap "../data/wmap/wmap.vrt" \
    --dem      "../data/dem/dem.vrt" \
    --radius   "${RADIUS}" \
    --out      "$EXTR_DIR/${DAM}_${RADIUS}" \
    --debug

fi

python3 szi_from_contourline.py \
  --name       "${DAMNAME}" \
  --infile     "/home/ad/briciera/dem4water/dem4water/data/DB_v5_KL.geojson" \
  --watermap   "$EXTR_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
  --dem        "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
  --radius     "${RADIUS}" \
  --pdbstep    5 \
  --pdbradius  500 \
  --elevoffset 60 \
  --elevsampling 5 \
  --tmp        "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
  --out        "$ROOT_DIR/${DAM}_${RADIUS}" \
  --debug

python3 cutline_score.py \
  --watermap "$EXTR_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
  --infile   "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
  --out      "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
  --debug

python3 cut_contourlines.py --debug \
  --name      "${DAMNAME}" \
  --dem       "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
  --info      "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
  --cut       "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
  --level     "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_contourlines@5m.json" \
  --out       "$ROOT_DIR/${DAM}_${RADIUS}"

python3 szi_to_model.py \
  --infile   "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_SZi.dat" \
  --daminfo  "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
  --outfile  "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_model.png" \
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
