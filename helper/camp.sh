#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2021 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2021
# ----

SRC_DIR="/home/ad/briciera/dem4water/dem4water"
DB_PATH="${SRC_DIR}/data/DB_v5_KL.geojson"
GT_PATH="${SRC_DIR}/data/validation_DB.json"
DEM_PATH="../data/dem/dem.vrt"
WMAP_PATH="../data/wmap/wmap.vrt"
ROOT_DIR=${2:-"/home/ad/briciera/scratch/HSV/camp_20210315"}
EXTR_DIR="/home/ad/briciera/scratch/HSV/Extracts"
RADIUS=${1:-10000}

module load otb/7.2-python3.7.2
[ -z "$OTB_MAX_RAM_HINT" ] && export OTB_MAX_RAM_HINT=4000
# [ -z "$OTB_LOGGER_LEVEL" ] && export OTB_LOGGER_LEVEL=CRITICAL
export OTB_LOGGER_LEVEL=CRITICAL
echo "OTB_MAX_RAM_HINT: $OTB_MAX_RAM_HINT"
echo "OTB_LOGGER_LEVEL: $OTB_LOGGER_LEVEL"

cd $SRC_DIR

declare -a DamDict=(   [3]='Agly'        [109]='Charpal'        [25]='Aussoue'      [213]='Grande Patures'
                      [32]='Balerme'      [83]='Cammazes'      [177]='Filhet'       [367]='Pla de Soulcem'
                     [190]='Galaube'     [193]='Ganguise'      [231]='Izourt'       [398]='Raschas'
                     [253]='Laparan'     [254]='Laprade'       [294]='Matemale'     [210]='Gouyre'
                     [463]='Salagou'     [320]='Montbel'       [336]='Olivettes'    [323]='Monts d'\''Orb (Avène)'
                     [392]='Puyvalador'  [348]='Pareloup'      [500]='Tordre'       [540]='Vinca'
                     [391]='Puylaurent'  [438]='Saint Ferréol' [440]='Saint Géraud' [462]='Sainte Peyres' )


for DAMID in "${!DamDict[@]}"; do

  DAM=${DamDict[$DAMID]// /_}

  [ -d "$EXTR_DIR/${DAM}_${RADIUS}" ] || mkdir -p "$EXTR_DIR/${DAM}_${RADIUS}"
  [ -d "$ROOT_DIR/${DAM}_${RADIUS}/tmp" ] || mkdir -p "$ROOT_DIR/${DAM}_${RADIUS}/tmp"
  [ -d "$ROOT_DIR/${DAM}_${RADIUS}/log" ] || mkdir -p "$ROOT_DIR/${DAM}_${RADIUS}/log"

  if [ -f "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif" ] \
    && [ -f "$EXTR_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" ] ; then
    echo "Extracts already available --> Skipping area_mapping."
  else
    # python3 area_mapping.py \
    python3 area_mapping.py --debug \
      --id       "${DAMID}" \
      --infile   "${DB_PATH}" \
      --watermap "${WMAP_PATH}" \
      --dem      "${DEM_PATH}" \
      --radius   "$RADIUS" \
      --out      "$EXTR_DIR/${DAM}_${RADIUS}" \
      2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/area_mapping.log"
  fi

  # python3 szi_from_contourline.py \
  python3 szi_from_contourline.py --debug \
    --id           "${DAMID}" \
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

  # python3 cutline_score.py \
  python3 cutline_score.py --debug \
    --watermap "$EXTR_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
    --infile   "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
    --out      "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
    2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/cutline_score.log"

  # python3 cut_contourlines.py \
  python3 cut_contourlines.py --debug \
    --dem      "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
    --info     "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
    --cut      "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
    --level    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_contourlines@1m.json" \
    --out      "$ROOT_DIR/${DAM}_${RADIUS}" \
    2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/log/cut_contourlines.log"

  # python3 szi_to_model.py \
  python3 szi_to_model.py --debug \
    --infile     "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_SZi.dat" \
    --daminfo    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
    --zminoffset 10 \
    --zmaxoffset 30 \
    --maemode    "first" \
    --outfile    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_model.png" \
    2>&1 | tee   "$ROOT_DIR/${DAM}_${RADIUS}/log/szi_to_model.log"

  python3 val_report.py --debug \
    -i "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_model.json" \
    -r "${GT_PATH}" \
    -o "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_report.json" \
    2>&1 | tee   "$ROOT_DIR/${DAM}_${RADIUS}/log/val_report.log"

done
