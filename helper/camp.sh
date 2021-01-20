#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2021 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2021
# ----

SRC_DIR="/home/ad/briciera/dem4water/dem4water"
ROOT_DIR="/home/ad/briciera/scratch/HSV/camp_20210120"
# DAM=${1:-"Agly"}
RADIUS=${1:-5000}

module load otb/7.2-python3.7.2
[ -z "$OTB_MAX_RAM_HINT" ] && export OTB_MAX_RAM_HINT=4000
# [ -z "$OTB_LOGGER_LEVEL" ] && export OTB_LOGGER_LEVEL=CRITICAL
export OTB_LOGGER_LEVEL=CRITICAL
echo "OTB_MAX_RAM_HINT: $OTB_MAX_RAM_HINT"
echo "OTB_LOGGER_LEVEL: $OTB_LOGGER_LEVEL"

cd $SRC_DIR

# declare -a StringArray=("Agly" "Cavayère" "Montbel" "Vinca" "Puyvalador"
                        # "Matemale" "Bouillouses" "Naguilhes" "Ganguise"
                        # "Aussoue" "Gimone" "Marcaoue" "Pessoulens"
                        # "Comberouger" "Malause" "Tordre" "Fontbouysse"
                        # "Pinet" "Bally" "Araing")

declare -a StringArray=("Agly"    "Astarac"  "Aussoue"             "Grandes Patures"
                        "Balerme" "Cammazes" "Filhet"              "Pla de Soulcem"
                        "Galaube" "Ganguise" "Izourt"              "Raschas"
                        "Laparan" "Laprade"  "Matemale"            "Gouyre"
                        "Salagou" "Montbel"  "Monts d'Orb (Avène)" "Olivettes")
                        # "Cap de Long" "Pareloup" "Portillon" "Puydarrieux"
                        # "Puylaurent" "Saint Ferréol" "Saint géraud" "Sainte Peyres"
                        # "Tordre" "Vinca" "Boues - Serres-Rustaing" "Miélan"
                        # "Charpal" "Villeneuve de la raho" "Gouyre" "Lac d'Oô"
                        # "Puyvalador" "Arrêt-Darré" "Louet" "Gimone"

for DAM in ${StringArray[@]}; do

  [ -d "$ROOT_DIR/${DAM}_${RADIUS}/tmp" ] || mkdir -p "$ROOT_DIR/${DAM}_${RADIUS}/tmp"

  if [ -f "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif" ] \
    && [ -f "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" ] ; then
    echo "Extracts already available --> Skipping area_mapping."
  else
    python3 area_mapping.py \
      --name     "$DAM" \
      --infile   "../data/synth_names.shp" \
      --watermap "../data/wmap/wmap.vrt" \
      --dem      "../data/dem/dem.vrt" \
      --radius   "$RADIUS" \
      --out      "$ROOT_DIR/${DAM}_${RADIUS}" 2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/area_mapping.log"
      # --debug
      # --watermap "../data/T31TDH/all_cumul.tif" \

  fi

  python3 szi_from_contourline.py \
    --name       "$DAM" \
    --infile     "../data/synth_names.shp" \
    --watermap   "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
    --dem        "$ROOT_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
    --radius     "$RADIUS" \
    --pdbstep    5 \
    --pdbradius  500 \
    --elevsampling 1 \
    --elevoffset 60 \
    --tmp        "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
    --out        "$ROOT_DIR/${DAM}_${RADIUS}"  2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/szi_from_contourline.log"
    # --out        "$ROOT_DIR/${DAM}_${RADIUS}" --debug 2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/szi_from_contourline.log"

  python3 cutline_score.py \
    --infile   "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
    --watermap "$ROOT_DIR/${DAM}_${RADIUS}/wmap_extract-$DAM.tif" \
    --out      "$ROOT_DIR/${DAM}_${RADIUS}/tmp" \
    --debug 2>&1 | tee "$ROOT_DIR/${DAM}_${RADIUS}/cutline_score.log"

  done
