u'i# dem4water

## Description

The aim of dem4water is to estimate water surface elevation (Z), surface area (S) and volume (V) relationships of reservoirs. It uses as input a digital elevation model (DEM) acquired after the reservoir construction and a water occurrence map (watermap).

## Installation

### Install a working version of GDAL

The chain version is set to GDAL 3.4.3 according to the TREX cluster.
For a recent version, update the setup.py file before installing the source.

### Install the source project

Once GDAL is installed, you can clone this repository then proceed to the installation:

```bash
# On TREX cluster load python which provide GDAL
# module load otb7.4/python3.8.4

# then install dem4water in a dedicated virtual env
python3 -m pip install --user virtualenv
virtualenv -p `which python3` <YOUR_VENV_PATH>
<YOUR_VENV_PATH>/bin/pip install -e .

```

The pip installation create the entry_point `dem4water` which will be used in the following. It refers to the source file `dem4water/cli.py`.

### Test your installation

Considering a fresh terminal, first load OTB the your virtual env. Always in this order.

```bash
# Source the otb modules
source /XXXX/otbenv.profile
# or module load otb7.4/python3.7.2
source <YOUR_VENV_PATH>/bin/activate

# then try to access entry point
dem4water -h

```

## Usage

### General information

dem4water steps can be divided in 5 steps, with an optional end step to validate a model from reference. Each step is usable in command line, and tools are provided to launch a full processing of several dams.

The algorithm workflow is the following:

1. `area_mapping.py` : extract the DEM and the watermap around the targeted dam.
2. `find_pdb_and_cutline.py` : automatic methods to find the PDB (Pied De Barrage @Santiago en anglais ?) and the site cutline.
3. `cutline_score.py` : provide a score of the founded cutline, which indicate if the cutline is often in the water body of the dam.
4. `cut_contourlines.py` : extract the contourlines which represent the valley containing the dam and the water body.
5. `szi_to_model.py` : compute the model from the previous results.
6. (opt) `val_report.py` : compare the model to a reference if provided.

Every script named in the previous list can be launched after loading the virtual env using `python3 XXXX.py -h`. There are many parameters for some function, this is why help scripts are provided.

### Mode campaign

This mode allows to launch all dam from a database file. All the parameters are handled using json files.

A template is provided in this git folder under `example` folder. But it can be automatically generated using the command `python3 dem4water/tools/generate_default_configuration.py -o /YOUR_OUTPUT_PATH` which create a file: `/YOUR_OUTPUT_PATH/campaign_template_file.json`.

You can edit every parameters in this file. If you not fill a parameter, for instance `reference` let it to `null` value.

Once the parameters are filled, you can simply launch the entry point using this json file:

```bash
dem4water campaign -json_campaign /YOUR_OUTPUT_PATH/campaign_template_file.json -scheduler_type PBS
```

The parameter `scheduler_type` is set by default to `PBS` which create a PBS file and send the `qsub` command for each dam to process.
The other allowed value is `local`, then each dam is processed sequentially.

At the end, you can find the output in the folder defined by `output_path` is the json file. Each dam is stored as `output_path/camp/dam_name`.
The others folders `extracts/dam_name` and `log` contain the dem and watermap extract, and the PBS logs.

#### Advanced parameters for campaign mode

- `input_force_list`: allow the user to provide the file generate by `dem4water/tools/generate_list_from_DB.py` as an input. The main usage is to remove some dam from the processing list.

### Mode single

The entry point allow to process a single dam too.

```bash
dem4water single -dam_json params_dam_name.json -scheduler_type PBS
```

The `params_dam_name.json` is created by the campaign mode as it contains informations dedicated to the dam, like the ID, the name etc.

### Mode autovalidation

This mode allow to launch the test dataset provided to the git folder.

### Use PBS launcher

It is possible to choose PBS ressources values for all modes. There is 4 arguments and they must be placed before the mode parameter.

- `-ram` : the amount of RAM per job, in GB (default = 60)
- `-cpu` : the amount of CPU per job (default = 12)
- `-walltime_hour`: indicate the hour for walltime (default = 1)
- `-walltime_minutes` : indicate the minutes for walltime (default = 0)

Exemple for the campaign mode:

```bash
dem4water --ram 50 --walltime_hour 4 campaign -json_campaign /YOUR_OUTPUT_PATH/campaign_template_file.json -scheduler_type PBS
```
### Automatic download DEM/Watermap 

DEM and Watermap can be automatically downloaded on (https://pypi.org/project/bmi-topography/) and (https://global-surface-water.appspot.com/download).

In /YOUR_OUTPUT_PATH/campaign_template_file.json :

```bash
{
    "campaign": {
        "output_path": "/YOUR_OUTPUT_PATH/",
        "watermap": null,
        "dem": null,
        "database": "/../data/andalousie/andalousie.geojson",
        "id_dam_column": "ID_SWOT",
        "dam_name_column": "DAM_NAME",
        "reference": "/../data/andalousie/andalousie_ref.json",
        "customs_files": null
    },
 .....

}
```

## Installation on Google Collab

- Create a gmail account. 
- Get your APIKEY : https://opentopography.org/blog/introducing-api-keys-access-opentopography-global-datasets

### Install the source project

First you need to download code in google drive :

```bash
# then install dem4water in a dedicated virtual env in google Collab :

!pip install virtualenv
!virtualenv dem4water
!source dem4water/bin/activate

# code in google drive 
from google.colab import drive
drive.mount('/content/drive')

!pip install -e /content/drive/MyDrive/Dem4Water/script/dem4water/

 #Install Gdal/OGR, rasterio, bmi-topography
!apt-get install gdal-bin -y


#export apikey
!export OPENTOPOGRAPHY_API_KEY='YOUR API KEY'

#file .sh Permission denied 
!chmod 755  /content/drive/MyDrive/Dem4Water/script/dem4water/dem4water/gen_contourline_polygons.sh

```

### Test your installation

```bash
!dem4water -h

```

## Usage


### Mode campaign

Be careful with paths in campaign_andalousie_params.json. No need wmap and dem parameters.
For keeping results, output folder must be in Mydrive

You can launch with :

```bash
!dem4water campaign -json_campaign /~/campaign_andalousie_params.json -scheduler_type local

```

See Contributing.md for all the information.

## Contributing

See Contributing.md for all the information.
