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

def export_ds_to_cog_gcloud(ds, out_path, path_to_temp, bucket):
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

def netCDF_to_EE_asset(path_to_nc = None, bands = None, path_to_cog_dir = None, project_folder = None, bucket_name = None, export_to_cog = False, save_to_image_collection = False, short_asset_path = None):
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
    cog_path = os.path.join(path_to_cog_dir, '{}.tif'.format(str(time).split(' ')[0]))
    print(cog_path)
    asset_path= short_asset_path + '/' + str(time).split(' ')[0]
    
    print(str(time).split(' ')[0]+'T00:00:00.000000000Z')
    ds_date = ds.sel(time = ds['time'] == time).squeeze()
    
    if export_to_cog:
      print('export_to_cog')
      export_ds_to_cog_gcloud(ds_date, cog_path, path_to_temp, bucket)
        
    if save_to_image_collection:
      
      if asset_exists(f'projects/{project_folder}/assets/' + asset_path):
        print(f'asset: {asset_path} already exists! skipping')
        #ee.data.deleteAsset(f'projects/{project_folder}/assets/' + asset_path)
      
      else:
        
        request = gen_request_body(f'gs://{bucket_name}/' + cog_path, str(time).split(' ')[0]+'T00:00:00.000000000Z')
        result = send_request(project_folder, asset_path, request, session)
        t.sleep(5)
        
  ds = None
    
#Authenticate EE
 
#service_account = 'rctm-07dea6fa718fdd0921f124263@rangelands-explo-1571664594580.iam.gserviceaccount.com'
credentials = compute_engine.Credentials()
ee.Initialize(credentials, project='rangelands-explo-1571664594580')
session = AuthorizedSession(credentials)

storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
path_to_temp = '/home/amullen/temp/'
project_folder = 'rangelands-explo-1571664594580'
site = 'HLD'
bucket_name='rangelands'
bucket = storage_client.get_bucket(bucket_name)


roi_file='/home/amullen/Rangeland-Carbon/res/site_footprints/HLD_tiles.txt'

#cogs need to point to 'HLD/tile/date'
cogs_fluxes=[]
cogs_C_stocks=[]

#asset needs to point to 'HLD/date'
assets_fluxes=[]
assets_C_stocks=[]

with open(roi_file) as f:
  tiles = [line.rstrip('\n') for line in f]
  
  for tile in tiles:
    site='/'.join(tile.split('/')[-2:])
    #open results net cdf for fluxes
    if site not in ['HLD/G1', 'HLD/G2', 'HLD/G3', 'HLD/G4', 'HLD/G5', 'HLD/G6', 'HLD/F1', 'HLD/F2', 'HLD/F3', 'HLD/F4', 'HLD/F5', 'HLD/F6', 
                    'HLD/E1', 'HLD/E2', 'HLD/E3', 'HLD/E4', 'HLD/E5', 'HLD/E6', 'HLD/D2', 'HLD/D3', 'HLD/D4', 'HLD/D5', 'HLD/C2', 'HLD/C3', 'HLD/C4', 'HLD/C5', 'HLD/C6']:
      path_to_fluxes = f'gs://rangelands/Ranch_Runs/{site}/results/flux_hist_weighted.nc'
      path_to_C_stocks = f'gs://rangelands/Ranch_Runs/{site}/results/C_stock_hist_weighted.nc'
      #C_stocks=rxr.open_rasterio(path_to_C_stocks)
      #C_stocks=C_stocks[['POC', 'HOC']].astype(np.float64)
      #C_stocks['TSOC'] = C_stocks['POC'] + C_stocks['HOC'] + 1000
      #C_stocks['TSOC'].attrs = C_stocks['POC'].attrs
      
      path_to_cog_dir_fluxes = f'Ranch_Runs/{site}/results/cog/fluxes'
      path_to_cog_dir_C_stocks = f'Ranch_Runs/{site}/results/cog/C_stocks'
      
      cogs_fluxes.append('gs://rangelands/' + path_to_cog_dir_fluxes)
      cogs_C_stocks.append('gs://rangelands/' + path_to_cog_dir_C_stocks)
  
      
      short_asset_path_fluxes = f'Ranches/{site}/fluxes'
      short_asset_path_C_stocks = f'Ranches/{site}/C_stocks'
      
      assets_fluxes.append(f'projects/{project_folder}/assets/' + short_asset_path_fluxes)
      assets_C_stocks.append(f'projects/{project_folder}/assets/' + short_asset_path_C_stocks)
      
      #netCDF_to_EE_asset(path_to_fluxes, ['GPP', 'NEE', 'Rh'], path_to_cog_dir_fluxes, project_folder, bucket_name, save_to_image_collection = False, short_asset_path = short_asset_path_fluxes)
      netCDF_to_EE_asset(path_to_nc = path_to_fluxes, bands = ['GPP', 'NEE', 'Rh'], path_to_cog_dir = path_to_cog_dir_fluxes, project_folder = project_folder, bucket_name = bucket_name, export_to_cog = False, save_to_image_collection = True, short_asset_path = short_asset_path_fluxes)
      
      #netCDF_to_EE_asset(C_stocks, ['POC', 'HOC', 'TSOC'], path_to_cog_dir_C_stocks, project_folder, bucket_name, save_to_image_collection = False, short_asset_path = short_asset_path_C_stocks)
      netCDF_to_EE_asset(path_to_nc = path_to_C_stocks, bands = ['POC', 'HOC'], path_to_cog_dir = path_to_cog_dir_C_stocks, project_folder = project_folder, bucket_name = bucket_name, export_to_cog = False, save_to_image_collection = True, short_asset_path = short_asset_path_C_stocks)
      
      C_stocks = None

#if asset_exists(f'projects/{project_folder}/assets/' + 'Ranches/HLD/fluxes'):
#  print(f'image collection exists')
        
#else:
#  ee.data.createAsset({'type': 'ImageCollection'}, f'projects/{project_folder}/assets/' + 'Ranches/HLD/fluxes')
  
#if asset_exists(f'projects/{project_folder}/assets/' + 'Ranches/HLD/C_stocks'):
#  print(f'image collection exists')
        
#else:
#  ee.data.createAsset({'type': 'ImageCollection'}, f'projects/{project_folder}/assets/' + 'Ranches/HLD/C_stocks')
         
#dates_fluxes = ''

#flux_dates = [d.split('/')[-1].split('.')[0] for d in utils.gs_listdir('Ranch_Runs/HLD/G1/results/cog/fluxes', bucket)]
#C_stock_dates = [d.split('/')[-1].split('.')[0] for d in utils.gs_listdir('Ranch_Runs/HLD/G1/results/cog/C_stocks', bucket)]

#print(cogs_fluxes)
#for i, d in enumerate(flux_dates):
#
#  cogs = [c + '/' + d + '.tif' for c in cogs_fluxes]
#  asset = f'projects/{project_folder}/assets/' + 'Ranches/HLD/fluxes/' + d
#  if asset_exists(asset):
#    print(f'{asset} exists')
    
#  else:
#    build_manifest_for_tiles(cogs, asset)
#  break
  
  #cogs_C_stocks = 
#print(cogs_C_stocks)

#print(assets_fluxes)
#print(assets_C_stocks)
#pool = Pool(maxtasksperchild=1)                      # Create a multiprocessing Pool
#pool.starmap(netCDF_to_EE_asset, args_fluxes)
#pool.close()
#pool.join()






