#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2021 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2021
# ----

DB_PATH=$1
DEM_PATH=$2
WMAP_PATH=$3
ROOT_DIR=$4

SRC_DIR="/home/mp/nicolasg/dem4water"
CAMP_DIR="${ROOT_DIR}/camp"
EXTR_DIR="${ROOT_DIR}/extracts"
RADIUS=""


module purge
module load otb/7.2-python3.7.2
[ -z "$OTB_MAX_RAM_HINT" ] && export OTB_MAX_RAM_HINT=4000
# [ -z "$OTB_LOGGER_LEVEL" ] && export OTB_LOGGER_LEVEL=CRITICAL
export OTB_LOGGER_LEVEL=CRITICAL
echo "OTB_MAX_RAM_HINT: $OTB_MAX_RAM_HINT"
echo "OTB_LOGGER_LEVEL: $OTB_LOGGER_LEVEL"

cd $SRC_DIR

declare -a DamDict=( [7052]='Andevalo'       [2917]='La Vinuela'        [2888]='Jose Toran'               
                     [7056]='Los Melonares'  [2915]='Beznar'            [2873]='San Rafael de Navallana'
                     [7054]='La Brena II'    [2911]='Los Bermejales'    [2869]='Yeguas'
                     [7053]='Arenoso'        [2909]='Puebla de Cazalla' [2867]='Puente Nuevo'
                     [2925]='Guadarranque'   [2904]='Iznajar'           [2861]='Tranco de Beas'
                     [2923]='La Concepcion'  [2900]='Piedras'           [2868]='Giribaile'
                     [2919]='Zahara'         [2894]='Vadomojon'         [2864]='Guadalen'
                                                                        [2863]='Encinarejo' )              
                              

#declare -a DamDict=( [2861]='Tranco de Beas' )
#declare -a DamDict=( [2925]='Guadarranque' )


for DAMID in "${!DamDict[@]}"; do

  DAM=${DamDict[$DAMID]// /_}

  echo ""
  echo "DAMID: $DAMID"
  echo "DB_PATH: $DB_PATH"
  echo "WMAP_PATH: $WMAP_PATH"
  echo "DEM_PATH: $DEM_PATH"
  echo "RADIUS: $RADIUS"
  echo "EXTR_DIR: $EXTR_DIR/${DAM}"
  echo "CAMP_DIR: $CAMP_DIR/${DAM}/log/area_mapping.log"
  echo ""

  [ -d "$EXTR_DIR/${DAM}" ] || mkdir -p "$EXTR_DIR/${DAM}"
  [ -d "$CAMP_DIR/${DAM}/tmp" ] || mkdir -p "$CAMP_DIR/${DAM}/tmp"
  [ -d "$CAMP_DIR/${DAM}/log" ] || mkdir -p "$CAMP_DIR/${DAM}/log"

  if [ -f "$EXTR_DIR/${DAM}/dem_extract-$DAM.tif" ] \
    && [ -f "$EXTR_DIR/${DAM}/wmap_extract-$DAM.tif" ] ; then
    echo "Extracts already available --> Skipping area_mapping."
  else
    echo "== AREA_MAPPING =="
    # python3 area_mapping.py \
    python3 area_mapping.py --debug \
      --id       "${DAMID}" \
      --infile   "${DB_PATH}" \
      --watermap "${WMAP_PATH}" \
      --dem      "${DEM_PATH}" \
      --radius   "$RADIUS" \
      --out      "$EXTR_DIR/${DAM}" \
     2>&1 | tee "$CAMP_DIR/${DAM}/log/area_mapping.log"
  fi

  # python3 szi_from_contourline.py \
  echo "== SZI_FROM_CONTOURLINE =="
  python3 szi_from_contourline.py --debug \
    --id           "${DAMID}" \
    --infile       "${DB_PATH}" \
    --watermap     "$EXTR_DIR/${DAM}/wmap_extract-$DAM.tif" \
    --dem          "$EXTR_DIR/${DAM}/dem_extract-$DAM.tif"  \
    --radius       "$RADIUS" \
    --pdbstep      5 \
    --pdbradius    500 \
    --elevsampling 1 \
    --elevoffset   60 \
    --tmp          "$CAMP_DIR/${DAM}/tmp" \
    --out          "$CAMP_DIR/${DAM}" \
    2>&1 | tee "$CAMP_DIR/${DAM}/log/szi_from_contourline.log"

    echo "== CUTLINE_SCORE =="
  # python3 cutline_score.py \
  python3 cutline_score.py --debug \
    --watermap "$EXTR_DIR/${DAM}/wmap_extract-$DAM.tif" \
    --infile   "$CAMP_DIR/${DAM}/${DAM}_cutline.json" \
    --out      "$CAMP_DIR/${DAM}/tmp" \
    2>&1 | tee "$CAMP_DIR/${DAM}/log/cutline_score.log"

    echo "== CUT_CONTOURLINES =="
  # python3 cut_contourlines.py \
  python3 cut_contourlines.py --debug \
    --dem      "$EXTR_DIR/${DAM}/dem_extract-$DAM.tif"  \
    --info     "$CAMP_DIR/${DAM}/${DAM}_daminfo.json" \
    --cut      "$CAMP_DIR/${DAM}/${DAM}_cutline.json" \
    --level    "$CAMP_DIR/${DAM}/${DAM}_contourlines@1m.json" \
    --out      "$CAMP_DIR/${DAM}" \
    2>&1 | tee "$CAMP_DIR/${DAM}/log/cut_contourlines.log"

    echo "== SZI_TO_MODEL =="
  # python3 szi_to_model.py \
  python3 szi_to_model.py --debug \
    --infile     "$CAMP_DIR/${DAM}/${DAM}_SZi.dat" \
    --daminfo    "$CAMP_DIR/${DAM}/${DAM}_daminfo.json" \
    --zminoffset 10 \
    --zmaxoffset 30 \
    --maemode    "first" \
    --outfile    "$CAMP_DIR/${DAM}/${DAM}_model.png" \
    2>&1 | tee   "$CAMP_DIR/${DAM}/log/szi_to_model.log"

  #python3 val_report.py --debug \
  #  -i "$CAMP_DIR/${DAM}/${DAM}_model.json" \
  #  -r "${GT_PATH}" \
  #  -o "$CAMP_DIR/${DAM}/${DAM}_report.json" \
  #  2>&1 | tee   "$CAMP_DIR/${DAM}/log/val_report.log"

done
