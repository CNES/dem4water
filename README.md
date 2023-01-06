# dem4water

## Description

The aim of dem4water is to estimate water surface elevation (Z), surface area (S) and volume (V) relationships of reservoirs. It uses as input a digital elevation model (DEM) acquired after the reservoir construction and a water occurrence map (watermap).

## Installation

### OTB is mandatory

dem4water is a OTB based processing chain. You must have a correct OTB installation before starting any installation process.

You can use your own compiled version or use [otb binaries](https://www.orfeo-toolbox.org/download/). If you are using an HPC (PBS only) cluster where OTB is installed you can load it with "module load otb". See compatibility section to know the list of tested OTB versions.

### Install the source project

Once OTB is installed, you can clone this repository then proceed to the installation:
```bash
# Source the otb modules
source /XXXX/otbenv.profile
# or module load otb7.4/python3.7.2

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

You can edit every parameters in this file. If you not fill a parameter, for instance `reference` let it to `null`.



### Mode single

### Mode autovalidation





## Contributing

See Contributing.md for all the information.