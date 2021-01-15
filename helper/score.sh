#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2021 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2021
# ----


SRC_DIR="/home/ad/briciera/dem4water/dem4water"
ROOT_DIR="/home/ad/briciera/scratch/HSV"


cd $SRC_DIR
export OTB_MAX_RAM_HINT=4000
export OTB_LOGGER_LEVEL=CRITICAL
module load otb/7.2-python3.7.2


python3 cutline_score.py \
  --infile   "$1" \
  --watermap "$2" \
  --out      "$PWD" \
  --debug
