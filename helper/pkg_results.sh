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

cd $CAMP_DIR
find . \
  -name    '*_daminfo.json' \
  -o -name '*_cutline.json' \
  -o -name '*@5m.json' \
  -o -name '*.png' \
  -o -name '*.log' \
  | cpio -pdm "${EXPO_DIR}"

zip "${EXPO_DIR}.zip" "${EXPO_DIR}"
