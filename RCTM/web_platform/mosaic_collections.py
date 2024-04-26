import sys
sys.path.insert(1, '../utils')
import utils
from google.auth.transport.requests import AuthorizedSession
import os
from google.auth import compute_engine
from utils import asset_exists
import ee


## this file takes multiple Earth Engine Image Collections and mosaics them by date
#authorization
credentials = compute_engine.Credentials()
ee.Initialize(credentials, project='rangelands-explo-1571664594580')
session = AuthorizedSession(credentials)

#ee project folder
project_folder = 'rangelands-explo-1571664594580'

#parent folder containing subdirectories that contain the image collections
path_to_parent_folder = f'projects/{project_folder}/assets/Ranches/HLD/'

#list of dictionaries containing subdirectories
subfolder_dicts = ee.data.listAssets(path_to_parent_folder)

#subdirectories
subfolders = [d['id'] for d in subfolder_dicts['assets']]

#dates to mosaic images by
dates = [d['id'].split('/')[-1] for d in ee.data.listAssets(os.path.join(subfolders[0], 'C_stocks'))['assets']]

#roi, for uploading mosaics as EE assets
roi=ee.FeatureCollection('projects/rangelands-explo-1571664594580/assets/Ranches/HLD_boundary')

#if asset_exists(f'projects/{project_folder}/assets/Ranches/HLD_C_stocks'):
#  ee.data.deleteAsset(f'projects/{project_folder}/assets/Ranches/HLD_C_stocks')
  
#if asset_exists(f'projects/{project_folder}/assets/Ranches/HLD_fluxes'):
#  ee.data.deleteAsset(f'projects/{project_folder}/assets/Ranches/HLD_fluxes')

#ee.data.createAsset({'type': 'ImageCollection'}, f'projects/{project_folder}/assets/Ranches/HLD_C_stocks')
#ee.data.createAsset({'type': 'ImageCollection'}, f'projects/{project_folder}/assets/Ranches/HLD_fluxes')

#iterate dates
for i, date in enumerate(dates):

  if asset_exists(f'projects/{project_folder}/assets/Ranches/HLD_C_stocks/{date}'):
    print(f'projects/{project_folder}/assets/Ranches/HLD_C_stocks/{date} exists, skipping!')
    
  else:
    print(f'mosaicing for projects/{project_folder}/assets/Ranches/HLD_C_stocks/{date}')

    #get paths to image collections, each subdirectory contains a 'C_stocks' and 'fluxes' ImageCollection
    C_stock_paths = [os.path.join(s, 'C_stocks', date) for s in subfolders if asset_exists(os.path.join(s, 'C_stocks', date))]
    
    #build new image collections based on date, containing images corresponding to date from each ImageCollection
    C_stock_col_date = ee.ImageCollection.fromImages([ee.Image(i).clipToBoundsAndScale(roi.geometry(), scale=30) for i in C_stock_paths])
    
    #mosaic new image collection
    C_stock_mosaic = C_stock_col_date.mosaic()
    C_stock_mosaic = C_stock_mosaic.clip(roi.geometry())
    C_stock_mosaic = C_stock_mosaic.set({'system:time_start': ee.Date(date).millis()})
    
    #get projection
    C_stock_projection = C_stock_mosaic.select('TSOC').projection().getInfo()
    
    #export C stocks 
    task = ee.batch.Export.image.toAsset(
      image=C_stock_mosaic,
      description=f'C_stock_mosaic_{date}',
      assetId=f'projects/{project_folder}/assets/Ranches/HLD_C_stocks/{date}',
      crs=C_stock_projection['crs'],
      scale=30,
      region = roi.geometry(),
      maxPixels = 1e13
    )
    task.start()
    
  if asset_exists(f'projects/{project_folder}/assets/Ranches/HLD_fluxes/{date}'):
    print(f'projects/{project_folder}/assets/Ranches/HLD_fluxes/{date} exists, skipping!')
  
  else:
    print(f'mosaicing for projects/{project_folder}/assets/Ranches/HLD_fluxes/{date}')

    #get paths to image collections, each subdirectory contains a 'C_stocks' and 'fluxes' ImageCollection
    flux_paths = [os.path.join(s, 'fluxes', date) for s in subfolders if asset_exists(os.path.join(s, 'fluxes', date))]
    
    #build new image collections based on date, containing images corresponding to date from each ImageCollection
    flux_col_date = ee.ImageCollection.fromImages([ee.Image(i).clipToBoundsAndScale(roi.geometry(), scale=30) for i in flux_paths])
    
    #mosaic new image collection
    flux_moasic = flux_col_date.mosaic()
    flux_moasic = flux_moasic.clip(roi.geometry())
    flux_moasic = flux_moasic.set({'system:time_start': ee.Date(date).millis()})
    
    #get projection
    flux_projection = flux_moasic.select('GPP').projection().getInfo()
    
    #export fluxes
    task = ee.batch.Export.image.toAsset(
      image=flux_moasic,
      description=f'flux_mosaic_{date}',
      assetId=f'projects/{project_folder}/assets/Ranches/HLD_fluxes/{date}',
      crs=flux_projection['crs'],
      scale=30,
      region = roi.geometry(),
      maxPixels = 1e13
      )
    task.start()
  
  