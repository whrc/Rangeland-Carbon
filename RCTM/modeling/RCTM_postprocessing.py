import rioxarray as rxr
import xarray as xr
from google.cloud import storage
import os
import sys
sys.path.insert(1, '../utils')
import utils
from datetime import datetime

bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
bucket = storage_client.get_bucket(bucket_name)
path_to_temp = '/home/amullen/temp/'

#RAP_in_dir = 'gs://rangelands/Ranch_Runs/HB/landcover/RAP_2019.tif'
#C_stocks_grass_path = 'gs://rangelands/Ranch_Runs/HB/results/HB_C_stocks_hist_grassland.nc'
#covariates_grass_path = 'gs://rangelands/Ranch_Runs/HB/results/HB_covariate_hist_grassland.nc'
#C_stocks_shrub_path = 'gs://rangelands/Ranch_Runs/HB/results/HB_C_stocks_hist_grass-shrub.nc'
#covariates_shrub_path = 'gs://rangelands/Ranch_Runs/HB/results/HB_covariate_hist_grass-shrub.nc'
#C_stocks_tree_path = 'gs://rangelands/Ranch_Runs/HB/results/HB_C_stocks_hist_grass-tree.nc'
#covariates_tree_path = 'gs://rangelands/Ranch_Runs/HB/results/HB_covariate_hist_grass-tree.nc'
#covariate_save_path = 'Ranch_Runs/HB/results/HB_covariate_hist_weighted.nc'
#C_stock_save_path = 'Ranch_Runs/HB/results/HB_C_stock_hist_weighted.nc'

def average_end_members(grass_path, shrub_path, tree_path, RAP_in_dir, save_path, columns=None):
  
  if columns is None:
    print('processing all data variables')
    print('loading grass')
    grass = rxr.open_rasterio(grass_path, masked=True)
    print('loading shrub')
    shrub = rxr.open_rasterio(shrub_path, masked=True)
    print('loading tree')
    tree = rxr.open_rasterio(tree_path, masked=True)
  
  else:
    print('loading grass')
    grass = rxr.open_rasterio(grass_path, masked=True)[columns]
    print('loading shrub')
    shrub = rxr.open_rasterio(shrub_path, masked=True)[columns]
    print('loading tree')
    tree = rxr.open_rasterio(tree_path, masked=True)[columns]
  
  print('loading landcover')
  RAP_landcover = rxr.open_rasterio(RAP_in_dir, masked=True)
  RAP_landcover = RAP_landcover.rio.reproject_match(grass)
  RAP_landcover.values[0] = (RAP_landcover.values[0] + RAP_landcover.values[3]) / (RAP_landcover.values[0] + RAP_landcover.values[3] + RAP_landcover.values[4] + RAP_landcover.values[5]) #grass ratio
  RAP_landcover.values[4] = (RAP_landcover.values[4]) / (RAP_landcover.values[0] + RAP_landcover.values[3] + RAP_landcover.values[4] + RAP_landcover.values[5]) #shrub ratio
  RAP_landcover.values[5] = (RAP_landcover.values[5]) / (RAP_landcover.values[0] + RAP_landcover.values[3] + RAP_landcover.values[4] + RAP_landcover.values[5]) #tree ratio
  
  print('averaging')
  weighted = (grass*RAP_landcover.values[0] + shrub*RAP_landcover.values[4] + tree*RAP_landcover.values[5])
  #print(weighted)
  
  print('saving')
  weighted.rio.write_crs(grass.rio.crs.to_wkt(), inplace=True)
  weighted.rio.write_transform(grass.rio.transform(), inplace=True)
  weighted.rio.nodata=0
  weighted.to_netcdf(os.path.join(path_to_temp, 'temp_average_write.nc'))
  utils.gs_write_blob(os.path.join(path_to_temp, 'temp_average_write.nc'), save_path, bucket)
  
  return

roi_file='/home/amullen/Rangeland-Carbon/res/site_footprints/HLD_tiles.txt'

with open(roi_file) as f:
  paths = [line.rstrip('\n') for line in f]
  
  for path in paths:
    
    site='/'.join(path.split('/')[-2:])
    
    if site in ['HLD/G1', 'HLD/G2', 'HLD/G3', 'HLD/G4', 'HLD/G5', 'HLD/G6', 'HLD/G4', 'HLD/F1', 'HLD/F2', 'HLD/F3', 'HLD/F4', 'HLD/F5', 'HLD/F6', 
                'HLD/E1', 'HLD/E2', 'HLD/E3', 'HLD/E4', 'HLD/E5', 'HLD/E6', 'HLD/D2', 'HLD/D3', 'HLD/D4','HLD/D5', 'HLD/C2', 'HLD/C3',
                'HLD/C4', 'HLD/C5', 'HLD/C6', 'HLD/B3', 'HLD/B4', 'HLD/B5', 'HLD/B6', 'HLD/A3', 'HLD/A4']:
      continue
      
    print(site)
    
    RAP_in_dir = f'gs://rangelands/Ranch_Runs/{site}/landcover/RAP_2019.tif'
    
    C_stocks_grass_path = f'gs://rangelands/Ranch_Runs/{site}/results/C_stocks_hist_grassland.nc'
    covariates_grass_path = f'gs://rangelands/Ranch_Runs/{site}/results/flux_hist_grassland.nc'
    
    C_stocks_shrub_path = f'gs://rangelands/Ranch_Runs/{site}/results/C_stocks_hist_grass-shrub.nc'
    covariates_shrub_path = f'gs://rangelands/Ranch_Runs/{site}/results/flux_hist_grass-shrub.nc'
    
    C_stocks_tree_path = f'gs://rangelands/Ranch_Runs/{site}/results/C_stocks_hist_grass-tree.nc'
    covariates_tree_path = f'gs://rangelands/Ranch_Runs/{site}/results/flux_hist_grass-tree.nc'
    
    covariate_save_path = f'Ranch_Runs/{site}/results/flux_hist_weighted.nc'
    C_stock_save_path = f'Ranch_Runs/{site}/results/C_stock_hist_weighted.nc'
    
      
    print(datetime.now())
    average_end_members(covariates_grass_path, covariates_shrub_path, covariates_tree_path, RAP_in_dir, covariate_save_path, columns=['GPP', 'NPP', 'Ra', 'Rh', 'NEE'])
    average_end_members(C_stocks_grass_path, C_stocks_shrub_path, C_stocks_tree_path, RAP_in_dir, C_stock_save_path)
    print(datetime.now())