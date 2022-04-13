# dem4water

## Description

The aim of dem4water is to estimate water surface elevation (Z), surface area (S) and volume (V) relationships of reservoirs. It uses as input a digital elevation model (DEM) acquired after the reservoir construction and a water occurrence map (watermap).

## Contributing

The project rely on [pre-commit](https://pre-commit.com/) to apply a given set of coding rules and formating to insure a certain level of code quality.

```bash
cd dem4water.git
python3 -m pip install --user pre-commit
pre-commit install
```

Pre-commit checks will then be performed at each commit to insure code formating and quality.

## Execution

The main script used to run the chain is [camp.sh](helper/camp.sh), gathering all steps together to process a list of dams. Using the script requires adapting the input definition section (first part of camp.sh file) as well as the dam list.
The script [camp.sh](helper/camp.sh) run a whole campaign and is in charge of setting the right environment and dependencies.
Each application embbed a documentation that can be accessed using --help.

It's also possible to use qsub in order to parallelize by dam process :
```
python3 run_processors.py dams_list dams_db dem_path wmap_path out_dir chain_dir

```
where "dam_list" is a csv like file containing "dam_id,dam_name" on each line.


## Examples

The following examples are directly extracted from [camp.sh](helper/camp.sh)

### Step 1 - area_mapping

This application extracts area specific data from the input watermap and DEM. Then, it uses the water map and concentric cercles on the DEM to identify upstream/downstream areas (= orthogonal to flow direction) and the dam bottom location. Radius defines the cropping window size around the dam location used for the following steps and it is expressed in meters.

``` sh
python3 area_mapping.py --debug \
  --id       "${DAMID}" \
  --infile   "${DB_PATH}" \
  --watermap "${WMAP_PATH}" \
  --dem      "${DEM_PATH}" \
  --radius   "$RADIUS" \
  --out      "$EXTR_DIR/${DAM}_${RADIUS}"
```

### Step 2 - szi_from_contourline

This application uses the dam bottom location to derive Z0 as well as the dam cut line (imaginary line cutting the valley on the dam position). It also extracts the contour lines from the DEM.

```sh
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
  --out          "$ROOT_DIR/${DAM}_${RADIUS}"
```

### Step 3 - cut_contourlines

Using the dam cutting line as well as the contour lines, this application estimates the virtual surfaces and the associated S(Zi).

```sh
python3 cut_contourlines.py --debug \
  --dem      "$EXTR_DIR/${DAM}_${RADIUS}/dem_extract-$DAM.tif"  \
  --info     "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
  --cut      "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_cutline.json" \
  --level    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_contourlines@1m.json" \
  --out      "$ROOT_DIR/${DAM}_${RADIUS}"
```

### Step 4 - szi_to_model

This application computes model parameters of the bathymetric relationship based on the best S(Zi) subset.

```sh
python3 szi_to_model.py --debug \
  --infile     "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_SZi.dat" \
  --daminfo    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_daminfo.json" \
  --zminoffset 10 \
  --zmaxoffset 30 \
  --maemode    "first" \
  --outfile    "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_model.png"
```

### Step 5 - val_report

This application compares the estimated model to a ground truth model to evaluate its quality.

```sh
python3 val_report.py --debug \
  -i "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_model.json" \
  -r "${GT_PATH}" \
  -o "$ROOT_DIR/${DAM}_${RADIUS}/${DAM}_report.json"
```
