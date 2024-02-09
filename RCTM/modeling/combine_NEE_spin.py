from google.cloud import storage
import os
import pandas as pd
import utils

bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/remote-sensing/gee_key.json')
bucket = storage_client.get_bucket(bucket_name)

SPIN_dir='Ameriflux/SPIN_AVG_ALL/'
NEE_dir='Ameriflux/NEEINS_AVG_ALL/'


NEE_dfs = []
#loop through site data, calculate reference values, append to df, save to new csv
for filename in utils.gs_listdir(NEE_dir, bucket):
  if filename.endswith('.csv'):
    site_df = pd.read_csv('gs://' + bucket_name + '/' + filename)
  
    site_df['system:index'] = pd.to_datetime(site_df['system:index'], format='%Y_%m_%d')
  
    site_df['site'] = [filename[:3]]*len(site_df)
  
    NEE_dfs.append(site_df)

NEE_df = pd.concat(NEE_dfs)

spin_dfs = []
#loop through site data, calculate reference values, append to df, save to new csv
for filename in utils.gs_listdir(SPIN_dir, bucket):
  if filename.endswith('.csv'):
    site_df = pd.read_csv('gs://' + bucket_name + '/' + filename)
  
    site_df['system:index'] = pd.to_datetime(site_df['system:index'], format='%Y_%m_%d')
  
    site_df['site'] = [filename[:-4]]*len(site_df)
  
    spin_dfs.append(site_df)

spin_df = pd.concat(spin_dfs)

print(NEE_df['site'].unique())
print(spin_df['site'].unique())
print(NEE_df.columns)
print(spin_df.columns)

NEE_df.to_csv('')
spin_df.to_csv('')
