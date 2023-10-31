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
  
if __name__ == "__main__":

  parser=argparse.ArgumentParser()
  parser.add_argument("--init", action='store_true')
  parser.add_argument("--covariate_status_file", help="status file to track covariate downloads")
  parser.add_argument("--local_status_file", help="local status file to track which rois have been processed")
  parser.add_argument("--bucket_name", help="bucket name")
  parser.add_argument("--prefix", help="path to site subdirectories")
  parser.add_argument("--start_date", help="start date (YYYY-mm-dd)")
  parser.add_argument("--end_date", help="end date (YYYY-mm-dd)")

  args=parser.parse_args()
  
  storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/remote-sensing/gee_key.json')
  bucket = storage_client.get_bucket(args.bucket_name)
  
  df_process_stats = pd.read_csv(args.local_status_file)
  
  if args.init:
  
    df_covariates = pd.DataFrame(columns=['roi', 'daymet'])
    df_covariates.to_csv(args.covariate_status_file, index=False)
    
  df_covariates = pd.read_csv(args.covariate_status_file)

  for index, row in df_process_stats.iterrows():
      site = row['roi'].split('/')[-1]
      print(row['roi'])
      
      if df_covariates['roi'].str.contains(row['roi']).any():
        if isinstance(df_covariates.loc[df_covariates['roi']==row['roi'], 'daymet'].values[0], str):
          print('already done')
          continue
      
      print('starting {}'.format(row['roi']))
      
      if row['roi'] in bad_sites:
        continue
      
      print('getting covariates')
      prefix=args.prefix
      out_dir = f'{prefix}/{site}/covariates/'
      modis_dir = f'{prefix}/{site}/modis/'
      roi_asset_path = row['roi']
      covariate_tuple = get_covariates.get_datasets(get_covariates.covariate_mapping, roi_asset_path, args.start_date, args.end_date)
      get_covariates.export_daymet(covariate_tuple, args.bucket_name, roi_asset_path, modis_dir, out_dir)
      daymet_datetime = dt.now() 
      
      row = pd.DataFrame({'roi':[row['roi']], 'daymet':[daymet_datetime]})
      df_covariates=pd.concat([df_covariates, row], ignore_index = True)
      df_covariates.to_csv(args.covariate_status_file, index=False)
      
# python batch_get_covariates.py --init --covariate_status_file='res/status_files/HLD_grid_covariates.csv' --local_status_file='res/status_files/HLD_grid_local.csv' --bucket_name='rangelands' --prefix='Ranch_Runs/HLD' --start_date='2002-01-01' --end_date='2023-01-01'