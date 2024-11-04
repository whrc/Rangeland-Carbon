import rioxarray as rxr
import xarray as xr
from google.cloud import storage
import ee
import sys
sys.path.insert(1, '../utils')
import utils
from google.auth.transport.requests import AuthorizedSession
from ee import oauth
import os
from google.auth import compute_engine
import json
import numpy as np
from uuid import uuid4
import time as t
from multiprocessing import Pool
import random
import pandas as pd
from utils import asset_exists, gs_listdir

def ds_to_cog_gcloud(ds, out_path, path_to_temp, bucket):
  #save as COG, Use tile dimensions of 256x256, include power of 2 overviews, LZW, compressor's zlevel to the maximum of 9
  #Use a predictor of 2. Sometimes a predictor of 3 will do better
  num = random.randint(0, 10000000000)
  ds.rio.to_raster(os.path.join(path_to_temp, f'temp_gtiff_{num}.tif'), driver="COG", nodata=0)
  utils.gs_write_blob(os.path.join(path_to_temp, f'temp_gtiff_{num}.tif'), out_path, bucket)
  os.remove(os.path.join(path_to_temp, f'temp_gtiff_{num}.tif'))
  
  return
       
def gen_request_body(path_to_cog, time):
  #path_to_cog: e.g. 'gs://ee-docs-demos/COG_demo.tif'
  #time: time string
  request = {
    'type': 'IMAGE',
    'gcs_location': {
      'uris': [path_to_cog] 
    },
    'properties': {
    },
    'startTime': time
  }
  
  return request
  
def send_request(project_folder, asset_id, request, session):
  # project_folder: Earth Engine enabled Cloud Project. e.g. 'your-project'
  # asset_id: A folder (or ImageCollection) name and the new asset name. e.g. 'cog-collection/your-cog-asset'
  
  url = 'https://earthengine.googleapis.com/v1alpha/projects/{}/assets?assetId={}'.format(project_folder, asset_id)
  
  response = session.post(
    url = url,
    data = json.dumps(request)
  )
  
  return json.loads(response.content)

def build_manifest_for_tiles(cogs, asset):
  #storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
  #bucket_name='rangelands'
  #bucket = storage_client.get_bucket(bucket_name)
  #  {
  #"name": "projects/earthengine-legacy/assets/users/username/some_folder/some_asset_id",
  #"tilesets": [
  #  {
  #    "sources": [
  #      {
  #        "uris": [
  #          "gs://bucket/N30W22.tif"
  #        ]
  #      },
  #      {
  #        "uris": [
  #          "gs://bucket/N31W22.tif"
  #        ]
  #      }
  #    ]
  #  }
  #]
  #}
  params = {
              'name': asset,
              'tilesets': [
                          {'sources': []}
                          ]
              }
  date=''
  for cog in cogs:
    #print(cog)
    #print('/'.join(cog.split('/')[3:])
    #print(storage.Blob(bucket=bucket, name=cog.split('/'][3:).exists(storage_client)))
    params['tilesets'][0]['sources'].append({'uris': [cog]})
    
    date = cog.split('/')[-1].split('.')[0]+'T00:00:00.000000000Z'
  print(date)
  params['start_time'] = date
  
  print(params)
  
  request_id = ee.data.newTaskId()[0]
  task = ee.data.startIngestion(request_id=request_id, params=params)
  print(task)
  
  while True:
    status = ee.data.getTaskStatus(task['id'])
    print(status)
    
    if status[0]['state'] == 'FAILED' or status[0]['state'] == 'SUCCESS'  or status[0]['state'] == 'SUCCEEDED':
      break
    t.sleep(5)
  
  #manifest to json
  return

def netCDF_to_cog(path_to_netcdf, path_to_cog_dir, 
                  bucket, path_to_temp, image_prefix = '', 
                  bands=None, date_range = None,
                  overwrite=False):

  if isinstance(path_to_netcdf, str):
    ds = rxr.open_rasterio(path_to_netcdf)
  
  if not bands is None:
    ds = ds[bands]

  for i in range(0, len(ds.indexes['time'])):

    time = ds.indexes['time'][i]

    if not date_range is None:  
      if pd.to_datetime(str(time).split(' ')[0]) < pd.to_datetime(date_range[0]) or pd.to_datetime(str(time).split(' ')[0]) > pd.to_datetime(date_range[1]):
        continue
    
    cog_path = os.path.join(path_to_cog_dir, '{}_{}.tif'.format(image_prefix, str(time).split(' ')[0]))
    print(cog_path)

    if overwrite==False:
      blob = bucket.blob(cog_path)

      if blob.exists():
        print('exists, skipping!\n')
        continue
    
    time_start = str(time).split(' ')[0]+'T00:00:00.000000000Z'
    print(time_start)
    
    ds_date = ds.sel(time = ds['time'] == time).squeeze()
    ds_date = ds_date.assign_attrs({'system:time_start': time_start})
    
    ds_to_cog_gcloud(ds_date, cog_path, path_to_temp, bucket)
  
def cog_to_asset(project_folder, short_asset_path, 
                 path_to_cog_dir, bucket, image_prefix = '',date_range = None, overwrite_images=False):

  #create folder if if necessarry
  if asset_exists(f'projects/{project_folder}/assets/' + '/'.join(short_asset_path.split('/')[:-1])):
    print(f'folder exists')
  
  else:
    ee.data.createAsset({'type': 'Folder'}, f'projects/{project_folder}/assets/' + '/'.join(short_asset_path.split('/')[:-1]))
  
  #create image collection if necesarry
  if asset_exists(f'projects/{project_folder}/assets/' + short_asset_path):
    print(f'image collection exists')

  else:
    ee.data.createAsset({'type': 'ImageCollection'}, f'projects/{project_folder}/assets/' + short_asset_path)
  
  blob_list = gs_listdir(path_to_cog_dir, bucket)

  for path in blob_list:

    time = path.split('/')[-1].split('_')[-1].split('.')[0]

    if not date_range is None:
      if pd.to_datetime(time) < pd.to_datetime(date_range[0]) or pd.to_datetime(time) > pd.to_datetime(date_range[1]):
        continue
    
    asset_path= short_asset_path + '/' + image_prefix + '_'+ time
    
    print(str(time).split(' ')[0]+'T00:00:00.000000000Z')
      
    if asset_exists(f'projects/{project_folder}/assets/' + asset_path):
      print(f'asset: {asset_path} already exists!')
        
      if overwrite_images == False:
          continue
          
      else:
        ee.data.deleteAsset(f'projects/{project_folder}/assets/' + asset_path)
    print(asset_path)   
    request = gen_request_body(os.path.join(f'gs://', bucket.name, path), 
                               str(time).split(' ')[0]+'T00:00:00.000000000Z')

    result = send_request(project_folder, asset_path, request, session)
    t.sleep(5)

def netCDF_to_EE_asset(path_to_nc = None, bands = None, path_to_cog_dir = None, 
                       project_folder = None, bucket_name = None, export_to_cog = False, 
                       save_to_image_collection = False, short_asset_path = None, 
                       date_range = None, overwrite_images = False):
  #open results net cdf for fluxes
  if isinstance(path_to_nc, str):
    ds = rxr.open_rasterio(path_to_nc)
    ds = ds[bands]
  else:
    ds = path_to_nc[bands]
  
  if save_to_image_collection:
    if asset_exists(f'projects/{project_folder}/assets/' + '/'.join(short_asset_path.split('/')[:-1])):
        print(f'folder exists')
        #ee.data.deleteAsset(f'projects/{project_folder}/assets/' + '/'.join(short_asset_path.split('/')[:-1]))
  
    else:
      ee.data.createAsset({'type': 'Folder'}, f'projects/{project_folder}/assets/' + '/'.join(short_asset_path.split('/')[:-1]))
  
  if save_to_image_collection:
    if asset_exists(f'projects/{project_folder}/assets/' + short_asset_path):
        print(f'image collection exists')
    else:
      ee.data.createAsset({'type': 'ImageCollection'}, f'projects/{project_folder}/assets/' + short_asset_path)
  
  for i in range(0, len(ds.indexes['time'])):

    time = ds.indexes['time'][i]
    
    if pd.to_datetime(str(time).split(' ')[0]) < pd.to_datetime(date_range[0]) or pd.to_datetime(str(time).split(' ')[0]) > pd.to_datetime(date_range[1]):
      continue
    
    cog_path = os.path.join(path_to_cog_dir, '{}.tif'.format(str(time).split(' ')[0]))
    print(cog_path)
    asset_path= short_asset_path + '/' + str(time).split(' ')[0]
    
    print(str(time).split(' ')[0]+'T00:00:00.000000000Z')
    ds_date = ds.sel(time = ds['time'] == time).squeeze()
    
    if export_to_cog:
      print('export_to_cog')
      ds_to_cog_gcloud(ds_date, cog_path, path_to_temp, bucket)
        
    if save_to_image_collection:
      
      if asset_exists(f'projects/{project_folder}/assets/' + asset_path):
        print(f'asset: {asset_path} already exists!')
        
        if overwrite_images == False:
          continue
          
        else:
          ee.data.deleteAsset(f'projects/{project_folder}/assets/' + asset_path)
        
      request = gen_request_body(f'gs://{bucket_name}/' + cog_path, str(time).split(' ')[0]+'T00:00:00.000000000Z')
      print(request)
      result = send_request(project_folder, asset_path, request, session)
      
      t.sleep(5)
        
  ds = None
    
#Authenticate EE
 
#service_account = 'rctm-07dea6fa718fdd0921f124263@rangelands-explo-1571664594580.iam.gserviceaccount.com'
credentials = compute_engine.Credentials()
ee.Initialize(credentials, project='rangelands-explo-1571664594580')
session = AuthorizedSession(credentials)

storage_client = storage.Client.from_service_account_json('/home/amullen/res/gee_key.json')
path_to_temp = '/home/amullen/temp/'
project_folder = 'rangelands-explo-1571664594580'
bucket_name='rangelands'
bucket = storage_client.get_bucket(bucket_name)

path_to_fluxes = f'gs://rangelands/NY-FieldPlots/RCTM_output/transient/flux_hist.nc'
path_to_C_stocks = f'gs://rangelands/NY-FieldPlots/RCTM_output/transient/C_stock_hist.nc'

path_to_cog_dir_fluxes = f'NY-FieldPlots/RCTM_output/transient/cog/flux_hist'
path_to_cog_dir_C_stocks = f'NY-FieldPlots/RCTM_output/transient/cog/C_stock_hist'
  
short_asset_path_fluxes = f'NY-FieldPlots/fluxes'
short_asset_path_C_stocks = f'NY-FieldPlots/C_stocks'

print('fluxes to COG')

#netCDF_to_cog(path_to_fluxes, path_to_cog_dir_fluxes, 
#                  bucket, path_to_temp, image_prefix = 'fluxes')

#netCDF_to_cog(path_to_C_stocks, path_to_cog_dir_C_stocks, 
#                  bucket, path_to_temp, image_prefix = 'C_stocks')

print('COG to asset')

cog_to_asset(project_folder, short_asset_path_fluxes, path_to_cog_dir_fluxes, 
             bucket, image_prefix = 'fluxes')

cog_to_asset(project_folder, short_asset_path_C_stocks, path_to_cog_dir_C_stocks, 
             bucket, image_prefix = 'C_stocks')








