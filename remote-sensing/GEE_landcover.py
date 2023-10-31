import ee
import time 
import numpy as np
from functools import partial
import starfm_preprocessing as psi
from google.cloud import storage
from datetime import datetime as dt
import pandas as pd
import sys
sys.path.insert(1, '../utils')
import utils


utils.authorize()
#utils.print_root_assets()
bucket_name = 'rangelands'

covariate_mapping = {
  'meteorology': 'NASA/ORNL/DAYMET_V4', #daymet
  'soil_texture': 'projects/rangelands-explo-1571664594580/assets/Covariates/Soils/SoilGrid100',
  'soil_temperature': 'projects/rangelands-explo-1571664594580/assets/Covariates/Soils/NLDAST', #NLDAST
  'soil_moisture': 'projects/rangelands-explo-1571664594580/assets/Covariates/Soils/NLDAS', #NLDAS
  'shortwave': 'NASA/NLDAS/FORA0125_H002',
  'NLCD_landcover': 'USGS/NLCD_RELEASES/2019_REL/NLCD',
  'RAP_landcover': 'projects/rap-data-365417/assets/vegetation-cover-v3'
}


def export_landcover(bucket_name, roi_asset_path, out_directory):

  roi=ee.FeatureCollection(roi_asset_path)
  
  NLCD = ee.ImageCollection(covariate_mapping['NLCD_landcover']).filterBounds(roi.geometry()).filterDate('2019-01-01', '2020-01-01').first().clip(roi)
  RAP = ee.ImageCollection(covariate_mapping['RAP_landcover']).filterBounds(roi.geometry()).filterDate('2019-01-01', '2020-01-01').first().clip(roi)
  print('exporting NLCD')
  task = ee.batch.Export.image.toCloudStorage(
                  image=NLCD,
                  scale=30,
                  crs='EPSG:4326',
                  region=roi.geometry(),
                  bucket=bucket_name,
                  fileNamePrefix='{}NLCD_2019'.format(out_directory))
  task.start()
    
  print('exporting RAP')              
  task = ee.batch.Export.image.toCloudStorage(
                  image=RAP,
                  scale=30,
                  crs='EPSG:4326',
                  region=roi.geometry(),
                  bucket=bucket_name,
                  fileNamePrefix='{}RAP_2019'.format(out_directory))
                  
  task.start()
  
  return

roi_file = '/home/amullen/Rangeland-Carbon/res/site_footprints/sites.txt'
#landcover_dir = '/HLD/G1/landcover/'
#export_landcover(bucket_name, roi_asset_path, landcover_dir)
with open(roi_file) as f:
  
  rois = [line.rstrip('\n') for line in f]
  sites = [roi.split('/')[-1] for roi in rois]
  
  for i, roi, in enumerate(rois):
    print(sites[i])
    export_landcover(bucket_name, roi, 'Ameriflux_sites/{}_starfm/landcover/'.format(sites[i]))
