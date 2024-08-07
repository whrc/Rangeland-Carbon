# Rangeland-Carbon
Repository for the Rangeland Carbon Tracking Model (RCTM). The RCTM uses remote sensing fusion coupled with a lightweight process-based model to quantify carbon stocks in vegetation and soil, and carbon fluxes (GPP, RECO, NEE). The model is currently configured to run at a five-day interval and 30-m spatial resolution.

Code for running the full process through pipelines is contained in the ```RCTM/``` directory. Scripts used for model developement are located in ```dev/```

The general process for running a new region goes as follows:
1. Download Landsat and MODIS imagery from Google Earth Engine
3. Download RCTM covariates from Google Earth Engine
3. Download landcover from Google Earth Engine
4. Preprocessing MODIS imagery
5. Run STARFM: Fuses 30-m Landsat and 250-m MODIS to generate a 30-m image for each day with a MODIS acquisition
6. QA/QC STARFM images, merge STARFM with covariates, and gap-fill as needed to be input into RCTM
7. Run the RCTM model
8. Postprocessing, analysis, uploading assets for Google Earth Engine - based web tool.

# Usage
Major processes have been bundled into three pipelines. At this stage of development, the pipelines are the preferable way to run the full process:

## Setup
### Installing the conda environent
The conda environment required to run RCTM is stored in conda_env/rctm.yml

To install the environment: 

```
conda env create -f rctm.yml
```

### Setting environment variables
Make sure to set the PYTHONPATH environment variable to the Rangeland-Carbon directory

```
export PYTHONPATH='{...}/Rangeland-Carbon'
```

We will define another environment variable 'RCTMPATH' that the model requires. This is the path to the RCTM directory.

```
export RCTMPATH='{...}/Rangeland-Carbon/RCTM'
```

### STARFM setup
The Spatial and Temporal Adaptive Reflectance Fusion Model (STARFM; Gao et al., 2006; Gao et al., 2015) is used to fuse 30-m spatial resolution Landsat imagery with daily 500-m MODIS imagery. We use STARFM to produce a consistent time series of 30-m surface reflectance imagery at a 5-day time step. The setup is straightforward.

To set up the source, navigate to Rangeland-Carbon/RCTM/remote_sensing/starfm_source/ and run ```$ make```.

STARFM source code was obtained from (https://www.ars.usda.gov/research/software/download/?softwareid=432&modecode=80-42-05-10). STARFM source code in this repository is nearly identical to code from the link, except for a minor bug fix in StarFM_util.c

References
 
Gao, F., Masek, J., Schwaller M. and Hall, F. On the blending of the Landsat and MODIS surface reflectance: predict daily Landsat surface reflectance. IEEE Transactions on Geoscience and Remote Sensing. 44 (8): 2207-2218. 2006.

Gao, F., Hilker, T., Zhu, X., Anderson, M. A., Masek, J., Wang, P. and Yang, Y. Fusing Landsat and MODIS data for vegetation monitoring, IEEE Geoscience and Remote Sensing Magazine. 3(3): 47-60. doi: 10.1109/MGRS.2015.2434351. 2015.

### Config file

A config file is required to run the model. Two examples of config files are located in examples/example_configs. **test_config_minimal.yaml** contains all required config variables and their definitions. **test_config.yaml** contains all possible variables that can be defined through the config file.

## 1. GEE_pipeline.py
   - handles downloading of landsat and modis imagery for input into starfm algorithm, covariate data downloads, and landcover downloads
   - only required inputs are a configuration file and a region of interest in .geojson format
   - can handle a .geojson file that contains multiple regions for downloading i.e. if you are working with a large area split into tiles
```
from RCTM.pipelines.GEE_pipeline import GEEPipeline

#initialize downloads pipeline with geometry
pipe = GEEPipeline(config_filename='/home/amullen/Rangeland-Carbon/examples/example_configs/test_config.yaml')
pipe.init_pipeline_from_geometry() #populates rows for each geometry in workflow status file to dictate pocessing
```
```init_pipeline_from_geometry()``` will attempt to upload a local .geojson to Earth Engine and populate a local status file to track progress. You can also specify in the config file to use a geometry asset
that already exists within Earth Engine. The following steps require the geometry to be uploaded as an EE asset, but there is currently no listener implemented to automaticly trigger downstream processes upon completion,
so you just need to wait a few minutes for your upload to complete.

Then call ```pipe.update_status()``` to update asset upload statuses.

Once the geometry asset(s) are uploaded, the following processes can be run sequentially.

```
pipe.download_landsat_imagery() #downloads landsat imagery for regions in workflow status file
pipe.download_modis_imagery() #uses landsat image dates from landsat status file to download every five days of modis + landsat dates
pipe.download_covariates() #uses modis status file to download covariates for each modis date
pipe.download_landcover()#downloads landcover for regions in workflow status file
```

The processing and downloading of imagery can take a while to complete, and you can update the statuses of individual image requests using ```pipe.update_status()```

Once all imagery, covariates, and landcover datasets are downloaded you may proceed.

## 2. RCTM_preprocessing_pipeline.py
   - implements a local MODIS image time series smoothing algorithm (optional)
   - runs the STARFM algorithm to fuse 30 m landsat imagery and daily MODIS imagery
   - merges STARFM NDVI and covariate datasets for ingestion into RCTM
   - generates spinup dataset based on specified temporal extent
   - fuses RAP and NLCD landcover products to determine spatial plant-functional types (PFTs) for RCTM
   - based on the generated PFT maps and parameter file, creates spatial set of parameters for running RCTM
```
from RCTM.pipelines.RCTM_preprocessing_pipeline import RCTMPrePipeline

pipe=RCTMPrePipeline(config_filename='/home/amullen/Rangeland-Carbon/examples/example_configs/test_config.yaml')
pipe.smooth_modis()
pipe.starfm()
pipe.starfm_postprocessing()
pipe.gen_RCTM_params()
```

## 3. RCTM_model_pipeline.py
   - Runs RCTM based on config file
   - Can run spinup in either point or spatial mode
   - Transient modeling only supports spatial runs currently
   - Can run with PFT maps or force a single PFT
```
from RCTM.pipelines.RCTM_model_pipeline import RCTMPipeline

pipe=RCTMPipeline(config_filename='/home/amullen/Rangeland-Carbon/examples/example_configs/test_config.yaml')
pipe.run_RCTM()
```
## Citing this work

Xia, Y., Sanderman, J., Watts, J.D. et al. (in review). Coupling Remote Sensing with a Process Model for the Simulation of Rangeland Carbon Dynamics. 
