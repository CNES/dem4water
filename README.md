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

## Data preparation

### Harmonize DAM_NAME

Before trying to use dem4water, ensure that the name in your database are all correct.
Some cases can cause issues during exectution:

- Special symbols : for instance "é,è,à,ç", but also "'"
- Non unicode characters which can be interpreted like "@,!,#" and interfere with bash


### Fuse shapefile to geojson

The script `convert_shp_to_geojson.py` fuse two shapefiles into one geojson database.

The inputs are:
- A shapefile containing polygones for water bodies
- A shapefile containing points for DAM informations

It outputs one geojson file in WGS84 projection.

Example:

```sh
python3 convert_shp_to_geojson.py INPE-V0-79-retenues.shp INPE-V0-79-barrages.shp INPE-V0-79-barrages.geojson
```

This script requires the installation of `pygeos` to use the `sjoin_nearest` geopandas function, required to handle the case of dam info not intersect exactly the water bodies.

## Execution

The main script used to run the chain is [compute_hsv.pbs](compute_hsv.pbs), gathering all steps together to process an unique dam.
```sh
qsub -v WD=\$PWD,DAM=dame_name,DAM_ID=dam_id,ID_FIELD=id_field,DB_PATH=database,DEM_PATH=dem,REF_MODEL=ref_model,WMAP_PATH=surfwater_map,ROOT_DIR=hsv_directory
                      compute_hsv.pbs
```

with:
- DAM: the dam name
- ID_FIELD: the name of field containing DAM ID in the database
- DB_PATH: the geojson database
- DEM_PATH, WMAP_PATH: path to the DEM and the watermap vrt
- ROOT_DIR: output path to store results


It's also possible to use [run_processors.py](run_processors.py) in order to process a whole campaign and parallelize the execution by dam process :

```sh
python3 run_processors.py dams_list dams_db dem_path wmap_path chain_dir out_dir [--id_field] [--ref_model] [--radius] [--elev_off] [--correct_folder]
```

where "dam_list" is a csv like file containing "dam_id,dam_name" on each line.

## Examples

The following examples are directly extracted from [compute_hsv.pbs](compute_hsv.pbs)
Each application have a dedicated help, which is available using `--help` parameter.

### Step 1 - area_mapping

This application extracts area specific data from the input watermap and DEM. Then, it uses the water map and concentric cercles on the DEM to identify upstream/downstream areas (= orthogonal to flow direction) and the dam bottom location. Radius defines the cropping window size around the dam location used for the following steps and it is expressed in meters.

```sh
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

## Performance report generator

Based on a selected reference campaign (Occitania & Andalousia), the goal is to provide a simple interface to run the reference campaign, collect the results and display them in a markdown dashboard. All the tool options can be explored using:

```sh
python3 perf/gen_report.py --help
python3 perf/gen_report.py campaign --help
python3 perf/gen_report.py report --help
python3 perf/gen_report.py dashboard --help
```

### Run the campaign

The following command will run the reference campaign and store the output in a subdirectory using the git revision hash.

```sh
python3 perf/gen_report.py campaign --outdir /tmp/campaign
```

### Collect the results

The following command will collect and creat a new report in /perf/reports/ with the current date and the revision used by the campaign runner.

```sh
python3 perf/gen_report.py report --indir /tmp/campaign/494a81a/
```

### Generate the updated dashboard

The following command will retrieve all the reports in /perf/reports and generate a markdown dashboard displaying the performance history in perf/dashboard/

```sh
python3 perf/gen_report.py dashboard
```

At the moment, the resulting dashboard looks like the following table:

| Dam ID     | Dam Name                | 20220823_f1a959e | 20220823_e107d70 | 20220822_a31d185 | 20220820_f1a959e |
| :--------- | :---------------------- | :--------------- | :--------------- | :--------------- | :--------------: |
| 2320028893 | Pareloup                | 0.0866           | 0.0866           | 0.0866           |      0.0866      |
| 2320028933 | Saint Geraud            | ☓                | ☓                | ☓                |        ☓         |
| 2160030183 | Avene                   | ☓                | ☓                | ☓                |        ☓         |
| 2160026973 | Salagou                 | -7.2406          | ☓                | ☓                |     -7.2406      |
| 2320033043 | Saints Peyres           | -28.2802         | -26.4857         | -26.4857         |     -28.2802     |
| 2320030823 | Cammazes                | 0.0619           | 0.0619           | 0.0619           |      0.0619      |
| 2320039133 | Astarac                 | 0.1639           | 0.1813           | 0.1813           |      0.1639      |
| 2320038733 | La Gimone               | 0.1938           | 0.1938           | 0.1938           |      0.1938      |
| 2320030233 | Montbel                 | -6.3207          | -6.0946          | -6.0946          |     -6.3207      |
| 2160030553 | Agly                    | 0.0310           | 0.0273           | 0.0273           |      0.0310      |
| 2320031293 | Laparan                 | -0.0354          | -0.0354          | -0.0354          |     -0.0354      |
| 2320031303 | Pla de Soulcem          | -0.0383          | -0.0383          | -0.0383          |     -0.0383      |
| 2160030123 | Vinca                   | -6.4771          | -6.2043          | -6.2043          |     -6.4771      |
| 2160029943 | Puyvalador              | -0.8011          | -0.5677          | -0.5677          |     -0.8011      |
| 2160029873 | Villeneuve la Raho      | NaN              | NaN              | NaN              |       NaN        |
| 2160028013 | Matemale                | 0.3337           | 0.3337           | 0.3337           |      0.3337      |
| 2160004602 | Casasola                | -0.0233          | ☓                | ☓                |     -0.0233      |
| 2160004532 | Limonero                | -0.0211          | ☓                | ☓                |     -0.0211      |
| 2310018492 | Arcos de la Frontera    | ☓                | ☓                | ☓                |        ☓         |
| 2310018272 | Almodovar               | -12.3807         | ☓                | ☓                |     -12.3807     |
| 2160004623 | Guadarranque            | 0.0163           | ☓                | ☓                |      0.0163      |
| 2160004663 | La Concepcion           | 0.0090           | ☓                | ☓                |      0.0090      |
| 2310018393 | Zahara                  | 0.0593           | ☓                | ☓                |      0.0593      |
| 2160004413 | La Vinuela              | -0.0111          | ☓                | ☓                |     -0.0111      |
| 2160004183 | Beznar                  | 0.0365           | ☓                | ☓                |      0.0365      |
| 2310022903 | Los Bermejales          | 0.0982           | ☓                | ☓                |      0.0982      |
| 2310024113 | Puebla de Cazalla       | 0.0563           | ☓                | ☓                |      0.0563      |
| 2310024073 | Iznajar                 | 0.0240           | ☓                | ☓                |      0.0240      |
| 2310020663 | Piedras                 | ☓                | ☓                | ☓                |        ☓         |
| 2310023943 | Vadomojon               | -0.0260          | ☓                | ☓                |     -0.0260      |
| 2310016612 | Andevalo                | ☓                | ☓                | ☓                |        ☓         |
| 2310023523 | Jose Toran              | -0.0750          | ☓                | ☓                |     -0.0750      |
| 2310027322 | Los Melonares           | NaN              | ☓                | ☓                |       NaN        |
| 2310024253 | La Brena II             | -0.0159          | ☓                | ☓                |     -0.0159      |
| 2310000173 | San Rafael de Navallana | -0.0637          | ☓                | ☓                |     -0.0637      |
| 2310022062 | Arenoso                 | ☓                | ☓                | ☓                |        ☓         |
| 2310022743 | Yeguas                  | -0.0438          | ☓                | ☓                |     -0.0438      |
| 2310027683 | Puente Nuevo            | 0.0992           | ☓                | ☓                |      0.0992      |
| 2310020933 | Tranco de Beas          | 0.1796           | ☓                | ☓                |      0.1796      |
| 2310021843 | Encinarejo              | ☓                | ☓                | ☓                |        ☓         |
| 2310022963 | Giribaile               | -0.2160          | ☓                | ☓                |     -0.2160      |
| 2310022023 | Guadalen                | -547824.0700     | ☓                | ☓                |   -547824.0700   |
| 2160003673 | Cuevas de Almanzora     | 0.1080           | ☓                | ☓                |      0.1080      |
| 2160004373 | Guadalhorce             | ☓                | ☓                | ☓                |        ☓         |
| 2160004383 | Guadalteba              | 0.6970           | ☓                | ☓                |      0.6970      |
| 2160004403 | Conde de Guadalhorce    | -0.0575          | ☓                | ☓                |     -0.0575      |
| 2160004253 | Rules                   | 0.0083           | ☓                | ☓                |      0.0083      |
| 2310018403 | Bornos                  | 0.0686           | ☓                | ☓                |      0.0686      |
| 2310018233 | Los Hurones             | 0.0996           | ☓                | ☓                |      0.0996      |
| 2310020153 | Guadalcacin 2           | -0.1066          | ☓                | ☓                |     -0.1066      |
| 2310018473 | Barbate                 | -0.0077          | ☓                | ☓                |     -0.0077      |
| 2310020223 | Celemin                 | 0.0189           | ☓                | ☓                |      0.0189      |
| 2160004643 | Charco Redondo          | -0.0706          | ☓                | ☓                |     -0.0706      |

### Use this tool for a custom study campaign

#### Prepare the data

To use this example replace all test-site occurences by the corresponding name site (for example andalousie or occitanie to use the git data).

1. Create a folder named test_site
2. In this folder create 4 files:
   - test-site.lst : output of generate_list_from_DB.py
   - test-site.geojson: the DAM database
   - test-site.cfg: a config file containing two entries, the watermap and the dem vrts
   - test-site_ref.json : the reference data for all models to compare with

3. Create an output folder in the location you want: /home/dev/dem4water_tests_campaign
   and the sub folder campaign, reports and dashboards

#### Prepare working env

If not done clone the dem4water git repository.
Put your sources on the correct branch
Ensure you are using python3

#### Launch the campaigns and produce dashboard

1. First campaign all parameters by default

```sh
python perf/gen_report.py campaign --outdir /home/dev/dem4water_tests_campaign/campaign --name all_default
```

2. Second campaign with a new value for the elevation offset

```sh
python perf/gen_report.py campaign --outdir /home/dev/dem4water_tests_campaign/campaign --name elev_15 --elev_off 15
```

3. Third campaign with new values for elevation and radius

```sh
python perf/gen_report.py campaign --outdir /home/dev/dem4water_tests_campaign/campaign --name elev_15_rad_100 --elev_off 15 --radius 100
```

These three campaigns used differents parameters but same source code, then it is possible to launch all three in same time. If you want to test algorithm modifications, be careful to not modify the source code while a campaign is running.

After several minutes (or hours), the campaigns should be done.

Then produce the reports for the campaign with:

```h
python perf/gen_report.py report -i /home/dev/dem4water_tests_campaign/campaign/all_default -o /home/dev/dem4water_tests_campaign/reports/
```

Update the -i argument to produce all the reports

Then produce the dashboard:

```sh
python perf/gen_report.py dashboard -i /home/dev/dem4water_tests_campaign/reports/ -o /home/dev/dem4water_tests_campaign/dashboards/
```