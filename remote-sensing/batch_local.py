import preprocess_starfm_imagery as psi
import pandas as pd
from google.cloud import storage
import numpy as np
import time 
from datetime import datetime as dt
import starfm
import argparse

bad_sites=['Kon']
    

if __name__ == "__main__":

  parser=argparse.ArgumentParser()
  parser.add_argument("--init", action='store_true')
  parser.add_argument("--gee_status_file", help="directory with input modis imagery")
  parser.add_argument("--local_status_file", help="directory with output modis imagery")
  parser.add_argument("--bucket_name", help="bucket name")
  #parser.add_argument("--modis_in_dir", help="directory with input modis imagery")
  #parser.add_argument("--modis_smooth_dir", help="directory with output modis imagery")
  #parser.add_argument("--landsat_dir", help="directory with input landsat imagery")
  #parser.add_argument("--sfm_out_dir", help="directory to output starfm imagery")

  args=parser.parse_args()
  
  storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/remote-sensing/gee_key.json')
  bucket = storage_client.get_bucket(args.bucket_name)
  
  df_process_stats = pd.read_csv(args.gee_status_file)
  
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
        imagery_downloaded = isinstance(df_process_stats.loc[index, 'imagery_downloaded'],str)
        time.sleep(1)
      
      if row['roi'] in df_local['roi'].to_list():
        print('this roi is done')
        continue
      
      print('smoothing modis')
      modis_in_dir = 'Ameriflux_sites/{}_starfm/modis_test/'.format(row['roi'])
      modis_smooth_dir = 'Ameriflux_sites/{}_starfm/modis_test_smooth_5d/'.format(row['roi'])
      psi.smooth_modis_col(modis_in_dir, modis_smooth_dir, bucket)
      modis_smoothed_datetime = dt.now() 
      
      
      print('running starfm')
      landsat_dir = 'Ameriflux_sites/{}_starfm/landsat_test/'.format(row['roi'])
      sfm_out_dir = 'Ameriflux_sites/{}_starfm/starfm_test_no_smooth/'.format(row['roi']) 
      
      runlist=starfm.get_runlist(modis_smooth_dir, landsat_dir, sfm_out_dir, args.bucket_name)
      starfm.run_starfm(runlist, args.bucket_name)
      starfm_datetime = dt.now()
      
      row = pd.DataFrame({'roi':[row['roi']], 'modis_smoothed':[modis_smoothed_datetime], 'starfm':[starfm_datetime]})
      df_local=pd.concat([df_local, row], ignore_index = True)
      df_local.to_csv(args.local_status_file)

#python batch_local.py --init=True --gee_status_file='ameriflux_processing_status.csv' --local_status_file='ameriflux_processing_status_local.csv' --bucket_name='rangelands'

#python batch_local.py --gee_status_file='ameriflux_processing_status.csv' --local_status_file='ameriflux_processing_status_local.csv' --bucket_name='rangelands'
