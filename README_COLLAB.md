# dem4water

## Installation

- Create a gmail account. 
- Get your APIKEY : https://opentopography.org/blog/introducing-api-keys-access-opentopography-global-datasets

### Install the source project

First you need to download code in google drive :
(git branch : integration_functions_lib)

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
!pip install rasterio
!pip install bmi-topography

#export apikey
!export OPENTOPOGRAPHY_API_KEY='YOUR API KEY'

#file .sh Permission denied 
!chmod 755  /content/drive/MyDrive/Dem4Water/script/dem4water/dem4water/gen_contourline_polygons.sh

```

### Test your installation

```bash

# then try to access entry point
dem4water -h

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
