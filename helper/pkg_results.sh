#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2021 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2021
# ----

CAMP_DIR=${1}
EXPO_DIR=${2}

# Extract JSONs
cd "${CAMP_DIR}" || exit

# find . \
  # -name '*.png' \
  # | cpio -pdm "${EXPO_DIR}"

# find . \
  # -name '*@*m.json' \
  # | cpio -pdm "${EXPO_DIR}"

# find . \
  # -name '*_cutline_points.json' \
  # | cpio -pdm "${EXPO_DIR}"

find . \
  -name    '*_daminfo.json' \
  -o -name '*_cutline.json' \
  -o -name '*_vSurfaces.json' \
  -o -name '*_model.json' \
  -o -name '*_report.json' \
  -o -name '*.png' \
  -o -name '*.dat' \
  -o -name '*.log' \
  | cpio -pdm "${EXPO_DIR}"

# Archive
cd "${EXPO_DIR}" || exit
zip -v -r "${EXPO_DIR}.zip" "."

# Clean the mess
cd "$HOME" || exit
rm -rf "${EXPO_DIR}"
