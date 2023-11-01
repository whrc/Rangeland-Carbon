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
import argparse

utils.authorize()
bucket_name = 'rangelands'

covariate_mapping = {
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

if __name__ == "__main__":
  parser=argparse.ArgumentParser()
  parser.add_argument("--roi_asset_path", help="path to roi Earth Engine Asset")
  parser.add_argument("--landcover_out_dir", help="path to export landcover")
  parser.add_argument("--bucket_name", help="bucket name")

  args=parser.parse_args()
  
  export_landcover(args.bucket_name, args.roi_asset_path, args.landcover_out_dir) 
  
#python GEE_landcover.py  --roi_asset_path="projects/rangelands-explo-1571664594580/assets/Shapefiles/HB_shp" --landcover_out_dir="Ranch_Runs/HB/landcover/" --bucket_name="rangelands"

#roi_file = '/home/amullen/Rangeland-Carbon/res/site_footprints/sites.txt'
#landcover_dir = '/HLD/G1/landcover/'
#export_landcover(bucket_name, roi_asset_path, landcover_dir)
#with open(roi_file) as f:
#  
#  rois = [line.rstrip('\n') for line in f]
#  sites = [roi.split('/')[-1] for roi in rois]
  
#  for i, roi, in enumerate(rois):
#    print(sites[i])
#    export_landcover(bucket_name, roi, 'Ameriflux_sites/{}_starfm/landcover/'.format(sites[i]))
