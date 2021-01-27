#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2021 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2021
# ----

SRC_DIR="/home/ad/briciera/dem4water/dem4water"
ROOT_DIR="/home/ad/briciera/scratch/HSV/camp_20210126"
DB_PATH="../data/DB_Barrages_Fixed_v3/DB_Barrages_Fixed.shp"
DEM_PATH="../data/dem/dem.vrt"
WMAP_PATH="../data/wmap/wmap_30.vrt"
RADIUS=${1:-5000}

module load otb/7.2-python3.7.2
[ -z "$OTB_MAX_RAM_HINT" ] && export OTB_MAX_RAM_HINT=4000
# [ -z "$OTB_LOGGER_LEVEL" ] && export OTB_LOGGER_LEVEL=CRITICAL
export OTB_LOGGER_LEVEL=CRITICAL
echo "OTB_MAX_RAM_HINT: $OTB_MAX_RAM_HINT"
echo "OTB_LOGGER_LEVEL: $OTB_LOGGER_LEVEL"

cd $SRC_DIR

declare -a StringArray=('Agly'        'Astarac'       'Aussoue'      'Grande Patures'
                        'Balerme'     'Cammazes'      'Filhet'       'Pla de Soulcem'
                        'Galaube'     'Ganguise'      'Izourt'       'Raschas'
                        'Laparan'     'Laprade'       'Matemale'     'Gouyre'
                        'Salagou'     'Montbel'       'Olivettes'    'Monts d'\''Orb (Avène)'
                        'Cap de Long' 'Pareloup'      'Tordre'       'Vinca'
                        'Puylaurent'  'Saint Ferréol' 'Saint géraud' 'Sainte Peyres'
                        'Charpal'     'Puyvalador'   'Villeneuve de la raho')

for DAMNAME in "${StringArray[@]}"; do

  DAM=${DAMNAME// /_}

  [ -d "$ROOT_DIR/${DAM}_${RADIUS}/tmp" ] || mkdir -p "$ROOT_DIR/${DAM}_${RADIUS}/tmp"
  [ -d "$ROOT_DIR/${DAM}_${RADIUS}/log" ] || mkdir -p "$ROOT_DIR/${DAM}_${RADIUS}/log"

  if [ -f "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif" ] \
    && [ -f "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" ] ; then
    echo "Extracts already available --> Skipping area_mapping."
  else
    # python3 area_mapping.py \
    python3 area_mapping.py --debug \
      --name     "${DAMNAME}" \
      --infile   "${DB_PATH}" \
      --watermap "${WMAP_PATH}" \
      --dem      "${DEM_PATH}" \
      --radius   "$RADIUS" \
      --out      "$ROOT_DIR/${DAM}_${RADIUS}" \
      2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/area_mapping.log"
  fi

  # python3 szi_from_contourline.py \
  python3 szi_from_contourline.py --debug \
    --name         "${DAMNAME}" \
    --infile       "${DB_PATH}" \
    --watermap     "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
    --dem          "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
    --radius       "$RADIUS" \
    --pdbstep      5 \
    --pdbradius    500 \
    --elevsampling 5 \
    --elevoffset   60 \
    --tmp          "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
    --out          "$ROOT_DIR/${DAM}_${RADIUS}" \
    2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/szi_from_contourline.log"

  python3 cutline_score.py \
    --infile   "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
    --watermap "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
    --out      "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
    --debug    \
    2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/cutline_score.log"

done
