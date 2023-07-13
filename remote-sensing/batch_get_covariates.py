import ee
import utils
import time 
import numpy as np
from functools import partial
import preprocess_starfm_imagery as psi
import GEE_fetch_starfm_imagery as fsi
from google.cloud import storage
from datetime import datetime as dt
import sys
import argparse
import os
import pandas as pd
import get_covariates

bad_sites=['Kon']
start_date = '2002-01-01'
end_date = '2023-01-01'

daymet_collection_path = 'NASA/ORNL/DAYMET_V4'
  
if __name__ == "__main__":

  parser=argparse.ArgumentParser()
  parser.add_argument("--init", action='store_true')
  parser.add_argument("--covariate_status_file", help="directory with output modis imagery")
  parser.add_argument("--local_status_file", help="directory with output modis imagery")
  parser.add_argument("--bucket_name", help="bucket name")

  args=parser.parse_args()
  
  storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/remote-sensing/gee_key.json')
  bucket = storage_client.get_bucket(args.bucket_name)
  
  df_process_stats = pd.read_csv(args.local_status_file)
  
  if args.init:
  
    df_covariates = pd.DataFrame(columns=['roi', 'daymet'])
    df_covariates.to_csv(args.covariate_status_file, index=False)
    
  df_covariates = pd.read_csv(args.covariate_status_file)

  for index, row in df_process_stats.iterrows():
      print(row['roi'])
      
      if df_covariates['roi'].str.contains(row['roi']).any():
        if isinstance(df_covariates.loc[df_covariates['roi']==row['roi'], 'daymet'].values[0], str):
          print('already done')
          continue
      
      print('starting {}'.format(row['roi']))
      
      if row['roi'] in bad_sites:
        continue
      
      print('getting daymet')
      
      roi_asset_path = 'projects/rangelands-explo-1571664594580/assets/Ameriflux_RS/{}/{}'.format(row['roi'], row['roi'])
      sfm_in_dir = 'Ameriflux_sites/{}_starfm/starfm_test/'.format(row['roi'])
      daymet_out_dir = 'Ameriflux_sites/{}_starfm/daymet/'.format(row['roi'])
      
      daymet_collection=get_covariates.get_daymet(start_date, end_date, roi_asset_path, daymet_collection_path)
      get_covariates.export_daymet(daymet_collection, args.bucket_name, roi_asset_path, sfm_in_dir, daymet_out_dir)
      daymet_datetime = dt.now() 
      
      row = pd.DataFrame({'roi':[row['roi']], 'daymet':[daymet_datetime]})
      df_covariates=pd.concat([df_covariates, row], ignore_index = True)
      df_covariates.to_csv(args.covariate_status_file, index=False)
      
# python batch_get_covariates.py --init --covariate_status_file='ameriflux_covariate_status.csv' --local_status_file='ameriflux_processing_status_local.csv' --bucket_name='rangelands'