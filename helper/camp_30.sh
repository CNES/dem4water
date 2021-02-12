#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2021 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2021
# ----

SRC_DIR="/home/ad/briciera/dem4water/dem4water"
DB_PATH="../data/DB_Barrages_Fixed_v3/DB_Barrages_Fixed.shp"
DEM_PATH="../data/dem/dem.vrt"
WMAP_PATH="../data/wmap/wmap_30.vrt"
ROOT_DIR="/home/ad/briciera/scratch/HSV/camp_20210212"
EXTR_DIR="/home/ad/briciera/scratch/HSV/Extracts"
RADIUS=${1:-10000}

module load otb/7.2-python3.7.2
[ -z "$OTB_MAX_RAM_HINT" ] && export OTB_MAX_RAM_HINT=4000
# [ -z "$OTB_LOGGER_LEVEL" ] && export OTB_LOGGER_LEVEL=CRITICAL
export OTB_LOGGER_LEVEL=CRITICAL
echo "OTB_MAX_RAM_HINT: $OTB_MAX_RAM_HINT"
echo "OTB_LOGGER_LEVEL: $OTB_LOGGER_LEVEL"

cd $SRC_DIR

declare -a StringArray=('Arrêt-Darré' 'Astarac' 'Gimone'    'Lac d'\''Oô'
                        'Puydarrieux' 'Miélan'  'Portillon' 'Louet'
                        'Cap de Long' 'Boues (Sere-Rustaing)')

for DAMNAME in "${StringArray[@]}"; do

  DAM=${DAMNAME// /_}

  [ -d "$EXTR_DIR/${DAM}_${RADIUS}" ] || mkdir -p "$EXTR_DIR/${DAM}_${RADIUS}"
  [ -d "$ROOT_DIR/${DAM}_${RADIUS}/tmp" ] || mkdir -p "$ROOT_DIR/${DAM}_${RADIUS}/tmp"
  [ -d "$ROOT_DIR/${DAM}_${RADIUS}/log" ] || mkdir -p "$ROOT_DIR/${DAM}_${RADIUS}/log"

  if [ -f "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif" ] \
    && [ -f "$EXTR_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" ] ; then
    echo "Extracts already available --> Skipping area_mapping."
  else
    # python3 area_mapping.py \
    python3 area_mapping.py --debug \
      --name     "${DAMNAME}" \
      --infile   "${DB_PATH}" \
      --watermap "${WMAP_PATH}" \
      --dem      "${DEM_PATH}" \
      --radius   "$RADIUS" \
      --out      "$EXTR_DIR/${DAM}_${RADIUS}" \
      2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/area_mapping.log"
  fi

  # python3 szi_from_contourline.py \
  python3 szi_from_contourline.py --debug \
    --name         "${DAMNAME}" \
    --infile       "${DB_PATH}" \
    --watermap     "$EXTR_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
    --dem          "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
    --radius       "$RADIUS" \
    --pdbstep      5 \
    --pdbradius    500 \
    --elevsampling 1 \
    --elevoffset   60 \
    --tmp          "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
    --out          "$ROOT_DIR/${DAM}_${RADIUS}" \
    2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/szi_from_contourline.log"

  python3 cutline_score.py \
    --watermap "$EXTR_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
    --infile   "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
    --out      "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
    --debug    \
    2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/cutline_score.log"

  # python3 cut_contourlines.py \
  python3 cut_contourlines.py --debug \
    --name     "${DAMNAME}" \
    --dem      "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
    --info     "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
    --cut      "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
    --level    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_contourlines@1m.json" \
    --fpoints  "../data/DB_in_points/Retenues-TETIS.shp" \
    --out      "$ROOT_DIR/${DAM}_${RADIUS}" \
    2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/cut_contourlines.log"

  # python3 szi_to_model.py \
  python3 szi_to_model.py --debug \
    --infile   "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_SZi.dat" \
    --outfile  "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_model.png" \
    2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/szi_to_model.log"

done
