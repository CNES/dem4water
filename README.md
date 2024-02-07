# dem4water

## Description

The aim of dem4water is to estimate water surface elevation (Z), surface area (S) and volume (V) relationships of
reservoirs. It uses as input a digital elevation model (DEM) acquired after the reservoir construction and a water
occurrence map (watermap).

## Installation

### Install a working version of GDAL

The chain version is set to GDAL 3.4.3 according to the TREX cluster.
For a recent version, update the setup.py file before installing the source.

### Install the source project

Once GDAL is installed, you can clone this repository then proceed to the installation:

```bash
# On TREX cluster load python which provide GDAL
# module load python/3.8.4

# then install dem4water in a dedicated virtual env
python3 -m pip install --user virtualenv
virtualenv -p `which python3` <YOUR_VENV_PATH>
<YOUR_VENV_PATH>/bin/pip install -e .

```

The pip installation create the entry_point `dem4water` which will be used in the following. It refers to the source
file `dem4water/cli.py`.

### Test your installation

Considering a fresh terminal, first load python then your virtual env. Always in this order.

```bash
# On TREX load the correct version of python
# module load python/3.8.4
source <YOUR_VENV_PATH>/bin/activate

# then try to access entry point
dem4water -h

```

## Usage

### General information

dem4water steps can be divided in 5 steps, with an optional end step to validate a model from reference. Each step is
usable in command line, and tools are provided to launch a full processing of several dams.

The algorithm workflow is the following:

1. `area_mapping.py` : extract the DEM and the watermap around the targeted dam.
2. `find_pdb_and_cutline.py` : automatic methods to find the PDB (Pied De Barrage @Santiago en anglais ?) and the site
   cutline.
3. `cutline_score.py` : provide a score of the founded cutline, which indicate if the cutline is often in the water body
   of the dam.
4. `cut_contourlines.py` : extract the contourlines which represent the valley containing the dam and the water body.
5. `szi_to_model.py` : compute the model from the previous results.
6. (opt) `val_report.py` : compare the model to a reference if provided.

Every script named in the previous list can be launched after loading the virtual env using `python3 XXXX.py -h`. There
are many parameters for some function, this is why help scripts are provided.

### Mode campaign

This mode allows to launch all dam from a database file. All the parameters are handled using json files.

A template is provided in this git folder under `example` folder. But it can be automatically generated using the
command `python3 dem4water/tools/generate_default_configuration.py -o /YOUR_OUTPUT_PATH` which create a
file: `/YOUR_OUTPUT_PATH/campaign_template_file.json`.

A parameter `mode` is available, it allows to choose between the `GDP` or `standard` mode.

You can edit every parameters in this file. If you not fill a parameter, for instance `reference` let it to `null`
value.

Once the parameters are filled, you can simply launch the entry point using this json file:

```bash
dem4water campaign -json_campaign /YOUR_OUTPUT_PATH/campaign_template_file.json -scheduler_type Slurm
```

The parameter `scheduler_type` is set by default to `PBS` which create a PBS file and send the `qsub` command for each
dam to process.
The other allowed value is `local`, then each dam is processed sequentially.

At the end, you can find the output in the folder defined by `output_path` is the json file. Each dam is stored
as `output_path/camp/dam_name`.
The others folders `extracts/dam_name`and `log` contain the dem and watermap extract, and the PBS logs.

#### Advanced parameters for campaign mode

- `input_force_list`: allow the user to provide the file generate by `dem4water/tools/generate_list_from_DB.py` as an
  input. The main usage is to remove some dam from the processing list.

### Mode single

The entry point allow to process a single dam too.

```bash
dem4water single -dam_json params_dam_name.json -scheduler_type PBS
```

The `params_dam_name.json` is created by the campaign mode as it contains informations dedicated to the dam, like the
ID, the name etc.

### Mode autovalidation

This mode allow to launch the test dataset provided to the git folder.

### Use PBS launcher

It is possible to choose Slurm resources values for all modes. There is 4 arguments and they must be placed before the
mode parameter.

- `-ram` : the amount of RAM per job, in GB (default = 60)
- `-cpu` : the amount of CPU per job (default = 12)
- `-walltime_hour`: indicate the hour for walltime (default = 1)
- `-walltime_minutes` : indicate the minutes for walltime (default = 0)

Exemple for the campaign mode:

```bash
dem4water --ram 50 --walltime_hour 4 campaign -json_campaign /YOUR_OUTPUT_PATH/campaign_template_file.json -scheduler_type PBS
```

### Prepare your run

You need to fill a configuration json file to launch the chain.

An example is provided in notebook/demo_cfg.json

You have to update the `output_path`

If you have a DEM already downloaded, fill the `dem` field with the path and set `retrieve_mode` to `local`.

Note that you need only one file or create a vrt

If you want to automatically download the DEM, let `dem` to `null` and set the parameter `retrieve_mode` to
`cop30`

DEM can be automatically downloaded on (https://pypi.org/project/bmi-topography/)

Other DEM sources will be added in future updates.

## Try it on Google Colab

- Create a Gmail account.
- The dataset contains a dam database and a reference set to produce plots at the end
- Get your APIKEY : https://opentopography.org/blog/introducing-api-keys-access-opentopography-global-datasets
- Copy the sources or clone them into your own drive
- Follow the notebook in "demo_colab.ipynb"

## Try it on a local notebook

- Clone the source
- Get the opentopography api key if you need to download the DEM
- If you already have a DEM, ensure that you have only one file or create a vrt
- Look at the notebook "demo_local.ipynb" and follow the instructions
- Note that you need to handle the installation of GDAL

## Database preparation

Dem4water provides two execution mode for searching the cutline.
One algorithm compute the gradient dot product over the DEM to find some information.
The other approach use more expert information to compute the cutline.

### Mode Gradient Dot Product (GDP)

For the GDP mode, the following inputs must be provided in the geojson database

- Lake geometry: the geometry field of the geojson
- The dam name: a string naming the dam
- An unique identifier for each dam in integer

More information can be stored in the geojson and will be ignored during processing.

### Mode classic

For the classic mode some information are mandatory in the input geojson and the columns names are imposed

- Dam localisation:
    - LONG_DD : the longitude in 4326 projection
    - LAT_DD : the latitude in 4326 projection
- Insider water point: a point inside the reservoir
    - LONG_WW : the longitude in 4326
    - LAT_WW : the latitude in 4326
- Dam elevation: the altitude in meters of the dam structure
    - DAM_LVL_M : float or integer value
- Dam name:
    - DAM_NAME: a str value
- Unique identifier:
    - name free for example ID_DB: a str or int value

### Custom corrections

The chain try is best to provide the dam location, the PDB and the cutline.

In some case, it is not possible to find correctly them or the reservoir is very complex.
Then you need to manually provide correction to the daminfo.json and/or the cutline.geojson file.

The best way to provide the corrections is to create a folder containing the corrected files.
Once done, update the configuration file and especially the `custom_files` field with the path to the folder.

Then launch the campaign. Messages will inform you that your corrections are found and used in the several steps of the
workflow.

## Contributing

See Contributing.md for all the information.
