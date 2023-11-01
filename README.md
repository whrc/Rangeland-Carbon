# Rangeland-Carbon
Repository for the Rangeland Carbon project. The repository is split into two components.
1. Remote Sensing
2. Modeling

'remote-sensing' contains the code to download MODIS and Landsat imagery and fuse them together with the STARFM algorithm. This directory also handles covariate and landcover downloads from Google Earth Engine, as well as all processing required to generate RCTM input data.

'modeling' contains code related to the RCTM model and parameters.

The general process for running a new region goes as follows:
1. Download Landsat and MODIS imagery from Google Earth Engine
   - requirements: roi uploaded as Earth Engine asset
   - GEE_starfm.py
   - batch_GEE.py
3. Download RCTM covariates from Google Earth Engine
  - covariates are curently downloaded corresponding to dates of downloaded MODIS imagery
  - requirements: path to modis imagery, roi uploaded as Earth Engine asset
  - GEE_covariates.py
  - batch_get_covariates.py
3. Download landcover from Google Earth Engine
  - requirements: roi uploaded as Earth Engine asset
  - GEE_landcover.py
4. Preprocessing MODIS imagery
  - Modis imagery is smoothed based on a rolling linear regression
  - starfm_preprocessing.py
5. Run STARFM: Fuses 30-m Landsat and 250-m MODIS to generate a 30-m image for each day with a MODIS acquisition
  - requirements: modis directory, landsat directory, output starfm directory
  - starfm.py
  - Steps 4 and 5 can automatically be run sequentially and in a batch environment using batch_local.py
6. QA/QC STARFM images, merge STARFM with covariates, and gap-fill as needed to be input into RCTM
  - starfm_postprocessing.py
    - This file generates gap-filled input and spinup data for the RCTM model in NetCDF format, and generates spatial parameters based on landcover.
7. Run the RCTM model
  - RCTM.py
