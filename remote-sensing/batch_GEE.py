import ee
import time 
import numpy as np
from functools import partial
import starfm_preprocessing as psi
import GEE_starfm as fsi
from google.cloud import storage
from datetime import datetime as dt
import sys
import argparse
import os
import pandas as pd
import sys
sys.path.insert(1, '../utils')
import utils

start_date = '2002-01-01'
end_date = '2023-01-01'

LS8_collection = 'LANDSAT/LC08/C02/T1_L2' #do no change
LS7_collection = 'LANDSAT/LE07/C02/T1_L2' #do not change
LS5_collection = 'LANDSAT/LT05/C02/T1_L2' #do not change
MODIS_collection = 'MODIS/061/MCD43A4' #do not change

bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')

bad_sites=['Kon']

def initialize_batch_with_rois(roi_file, status_filename):
  with open(roi_file) as f:
  
    sites = [line.rstrip('\n') for line in f]
    #print('sites: {}'.format(sites))
  geometries=[]
  
  for site in sites:
    if site in bad_sites:
        geometries.append('')
        continue
    roi_asset_path = site
    roi=ee.FeatureCollection(roi_asset_path)
    geometries.append(roi.geometry().getInfo())
  
  df_process_status=pd.DataFrame(columns=['roi', 'geometry', 'landsat_total', 'landsat_success', 'landsat_failed', 'modis_total', 'modis_success', 'modis_failed', 'imagery_downloaded'])
  df_process_status['roi']=sites
  df_process_status['geometry']=geometries
  print(df_process_status)
  
  df_process_status.to_csv(status_filename, index=False)
  
  return

def poll_for_success(task_ids):
    
  success=0
  failed=0
    
  for l_id in task_ids:
    status_string=''
      
    while status_string!='COMPLETED' and status_string!='FAILED':
      status=ee.data.getTaskStatus(l_id)
      status_string=status[0]['state']
      time.sleep(1)
      
        
    if status_string=='COMPLETED':
      success+=1
       
    if status_string=='FAILED':
      print(status)
      failed+=1
      
  return success, failed

def batch(status_filename, prefix):

  df_process_stats=pd.read_csv(status_filename)
  
  for index, row in df_process_stats.iterrows():
  
    if isinstance(row['imagery_downloaded'], str):
      continue
      
    print('starting {}'.format(row['roi']))

    if row['roi'] in bad_sites:
      continue
    
    site = row['roi'].split('/')[-1]
    #get roi asset path  
    roi_asset_path = row['roi']
    
    #get landsat
    landsat_out_dir = f'{prefix}/{site}/landsat_v2/'
    landsat = fsi.get_landsat(LS8_collection, LS7_collection, LS5_collection, roi_asset_path, start_date, end_date)
    landsat_dates, landsat_task_ids = fsi.export_landsat_collection(landsat, roi_asset_path,bucket_name, landsat_out_dir, overwrite=False)
    df_process_stats.loc[index,'landsat_total'] = len(landsat_task_ids)
     
    #get modis   
    modis_out_dir = f'{prefix}/{site}/modis/'
    modis = fsi.get_MODIS(MODIS_collection, roi_asset_path, start_date, end_date)
    modis_task_ids = fsi.export_modis_collection(modis, roi_asset_path, bucket_name, modis_out_dir, landsat_dates=landsat_dates, overwrite=False)
    df_process_stats.loc[index,'modis_total'] = len(modis_task_ids) 
    
    print('waiting for landsat')
    landsat_succes, landsat_failed = poll_for_success(landsat_task_ids)
    df_process_stats.loc[index,'landsat_success'] = landsat_succes
    df_process_stats.loc[index,'landsat_failed'] = landsat_failed
    print('success')
    
    print('waiting for modis')
    modis_success, modis_failed = poll_for_success(modis_task_ids)
    df_process_stats.loc[index,'modis_success'] = modis_success
    df_process_stats.loc[index,'modis_failed'] = modis_failed
    print('success')
    
    df_process_stats.loc[index,'imagery_downloaded'] = dt.now() 
    df_process_stats.to_csv(status_filename, index=False)
    
    
if __name__ == "__main__":
  parser=argparse.ArgumentParser()
  parser.add_argument("--init", action='store_true', help='Use flag if running the batch task for the first time. Creates new status file.')
  parser.add_argument("--footprint_list", help="List of footprints to download imagery for")
  parser.add_argument("--status_file", help="GEE task status file")
  parser.add_argument("--prefix", help="path to export imagery subdirectories")

  args=parser.parse_args()
  
  if args.init:
    initialize_batch_with_rois(args.footprint_list, args.status_file)
     
  batch(args.status_file, args.prefix)
  