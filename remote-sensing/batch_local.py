import starfm_preprocessing as psi
import pandas as pd
from google.cloud import storage
import numpy as np
import time 
from datetime import datetime as dt
import starfm
import argparse
import os

bad_sites=['Kon']
    
os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = '/home/amullen/Rangeland-Carbon/res/gee_key.json'

if __name__ == "__main__":

  parser=argparse.ArgumentParser()
  parser.add_argument("--init", action='store_true')
  parser.add_argument("--gee_status_file", help="gee imagery download status")
  parser.add_argument("--local_status_file", help="status file for local processes")
  parser.add_argument("--bucket_name", help="bucket name")
  parser.add_argument("--in_dir", help="directory to save smoothed modis imagery")

  args=parser.parse_args()
  
  storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
  bucket = storage_client.get_bucket(args.bucket_name)
  
  df_process_stats = pd.read_csv(args.gee_status_file)
  
  # if initializing new local processing, create status file with appropriate fields
  if args.init:
    df_local = pd.DataFrame(columns=['roi', 'modis_smoothed', 'starfm'])
    df_local.to_csv(args.local_status_file)
    
  df_local = pd.read_csv(args.local_status_file)

  for index, row in df_process_stats.iterrows():
  
      print('starting {}'.format(row['roi']))
      
      if row['roi'] in bad_sites:
        continue
        
      imagery_downloaded=False
      
      print('waiting for imagery')
      while imagery_downloaded==False:
        
        df_process_stats = pd.read_csv(args.gee_status_file)
        imagery_downloaded = isinstance(df_process_stats.loc[index, 'imagery_downloaded'],str) #checks GEE status file to see if string occupies 'imagery downloaded fiels'
        time.sleep(1)
      
      #ignore roi if already processed
      if row['roi'] in df_local['roi'].to_list():
        print('this roi is done')
        continue
      
      print('smoothing modis')
      modis_in_dir = os.path.join(args.in_dir, row['roi'].split('/')[-1], 'modis/')
      modis_smooth_dir = os.path.join(args.in_dir, row['roi'].split('/')[-1], 'modis_smooth_v2')
      psi.smooth_modis_col(modis_in_dir, modis_smooth_dir, bucket)
      modis_smoothed_datetime = dt.now()
      
      #modis_in_dir = os.path.join(args.in_dir, row['roi'].split('/')[-1] + '/modis/')
      #modis_smooth_dir = os.path.join(args.in_dir, row['roi'].split('/')[-1] + '/modis_smooth')
      #psi.smooth_modis_col(modis_in_dir, modis_smooth_dir, bucket)
      #modis_smoothed_datetime = dt.now() 
      
      
      print('running starfm')
      landsat_dir = os.path.join(args.in_dir, row['roi'].split('/')[-1], 'landsat_v2')
      sfm_out_dir = os.path.join(args.in_dir, row['roi'].split('/')[-1], 'starfm_v2')
      
      #landsat_dir = os.path.join(args.in_dir, row['roi'].split('/')[-1] + '/landsat')
      #sfm_out_dir = os.path.join(args.in_dir, row['roi'].split('/')[-1] + '/starfm')
      
      runlist=starfm.get_runlist(modis_smooth_dir, landsat_dir, sfm_out_dir, args.bucket_name)
      print(runlist['outname'])
      starfm.run_starfm(runlist, args.bucket_name)
      starfm_datetime = dt.now()
      
      row = pd.DataFrame({'roi':[row['roi']], 'modis_smoothed':[modis_smoothed_datetime], 'starfm':[starfm_datetime]})
      df_local=pd.concat([df_local, row], ignore_index = True)
      df_local.to_csv(args.local_status_file, index=False)

#python batch_local.py --init=True --gee_status_file='ameriflux_processing_status.csv' --local_status_file='ameriflux_processing_status_local.csv' --bucket_name='rangelands'

#python batch_local.py --init --gee_status_file='res/status_files/HLD_grid_v2.csv' --local_status_file='res/status_files/HLD_grid_v2_local.csv' --bucket_name='rangelands' --in_dir='Ranch_Runs/HLD'
