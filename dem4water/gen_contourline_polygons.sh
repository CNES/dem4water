#!/usr/bin/env bash
# ----
# :author: Aurélien Bricier <aurelien.bricier@csgroup.eu>
# :organization: CS Group
# :copyright: 2021 CS Group. All rights reserved.
# :license: see LICENSE file
# :created: 2021
# ----

TMP_DIR=${6}

for i in $(seq "$2" "$3" "$4"); do
  gdal_contour -p -amax ID -fl "$i" "$1" "${TMP_DIR}/__TMP__poly_${i}.json"
  ogr2ogr -where ID="$i" "${TMP_DIR}/__TMP__poly_${i}_FIL.json" "${TMP_DIR}/__TMP__poly_${i}.json"
  rm "${TMP_DIR}/__TMP__poly_${i}.json"
done

consolidated_file="${TMP_DIR}/__TMP__consolidated.shp"
if [ -f "$consolidated_file" ]; then
  rm "$consolidated_file"
fi

for i in $(find "${TMP_DIR}" -name '__TMP__poly_*_FIL.json' | sort -r ); do
  if [ ! -f "$consolidated_file" ]; then
    # echo "new "
    # first file - create the consolidated output file
    ogr2ogr -f "ESRI Shapefile" "$consolidated_file" "$i"
  else
    # echo "update"
    # update the output file with new file content
    ogr2ogr -f "ESRI Shapefile" -update -append "$consolidated_file" "$i"
  fi
done

ogr2ogr "${5}" "${consolidated_file}"

rm "${TMP_DIR}"/__TMP__poly_*.json
rm "${TMP_DIR}"/__TMP__consolidated.*
