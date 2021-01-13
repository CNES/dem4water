#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2021 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2021
# ----

SRC_DIR="/home/ad/briciera/dem4water/dem4water"
ROOT_DIR="/home/ad/briciera/scratch/HSV/camp"
# DAM=${1:-"Agly"}
RADIUS=${1:-5000}

module load otb/7.2-python3.7.2
[ -z "$OTB_MAX_RAM_HINT" ] && export OTB_MAX_RAM_HINT=4000
# [ -z "$OTB_LOGGER_LEVEL" ] && export OTB_LOGGER_LEVEL=CRITICAL
export OTB_LOGGER_LEVEL=CRITICAL
echo "OTB_MAX_RAM_HINT: $OTB_MAX_RAM_HINT"
echo "OTB_LOGGER_LEVEL: $OTB_LOGGER_LEVEL"

cd $SRC_DIR

declare -a StringArray=("Agly" "Cavayère" "Montbel" "Vinca" "Puyvalador" "Matemale" "Bouillouses" "Naguilhes")

for DAM in ${StringArray[@]}; do

  [ -d "$ROOT_DIR/${DAM}_${RADIUS}/tmp" ] || mkdir -p "$ROOT_DIR/${DAM}_${RADIUS}/tmp"

  if [ -f "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif" ] \
    && [ -f "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" ] ; then
    echo "Extracts already available --> Skipping area_mapping."
  else
    python3 area_mapping.py \
      --name     "$DAM" \
      --infile   "../data/synth_names.shp" \
      --watermap "../data/T31TDH/all_cumul.tif" \
      --dem      "../data/dem/dem.vrt" \
      --radius   "$RADIUS" \
      --out      "$ROOT_DIR/${DAM}_${RADIUS}" 2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/area_mapping.log"
      # --debug

  fi

  python3 szi_from_contourline.py \
    --name       "$DAM" \
    --infile     "../data/synth_names.shp" \
    --watermap   "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
    --dem        "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
    --radius     "$RADIUS" \
    --pdbstep    5 \
    --pdbradius  500 \
    --tmp        "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
    --out        "$ROOT_DIR/${DAM}_${RADIUS}" \
    --debug 2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/szi_from_contourline.log"

  done
