import rioxarray as rxr
import xarray as xr
from google.cloud import storage
import os
import sys
sys.path.insert(1, '../utils')
import utils
from datetime import datetime
import numpy as np

bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/res/gee_key.json')
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

def average_end_members(grass_path, shrub_path, tree_path, RAP_in_dir, save_path, bucket, columns=None):
  """
  RCTM v1 calibrated parameters can produce visible discontinuities between pixel outputs from different pfts. 
  The discontinuities are most visible between grass-tree and other pfts. To alleviate this it is often advantageous 
  to run the model for all PFTs, and create a weighted average of outputs using the RAP landcover fractions. This function creates a
  weighted averaging of results from end members based on RAP cover fractions.
  
  Args:
  grass_path (str): path to grassland netcdf
  shrub_path (str): path to grass-shrub netcdf 
  tree_path (str): path to grass-tree netcdf 
  RAP_in_dir (str): path to RAP landcover (.tif)
  save_path (str): path to save result netcdf on google cloud bucket
  columns list(str): None or subset of data variable names to process

  Returns:
     None
  
  """
  
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
  
def spatial_filter(netcdf_file, save_path, bucket, sd_threshold=10):

  """
  Spatial filtering for netCDF files with dimensions (time, y, x). Could also work on 2d data with
  dimensions (y,x). Saves filtered dataset to Google Cloud.
  
  Args:
  netcdf_file (str): path to netcdf
  save_path (str): path to save result netcdf on google cloud bucket
  bucket (google.cloud.storage.bucket): bucket to save output in
  sd_threshold (numeric): None or subset of data variable names to process

  Returns:
     None
  """
  ds = rxr.open_rasterio(netcdf_file, masked=True)
  
  median = ds.median(dim=['x', 'y'])
  std = ds.std(dim=['x', 'y'])
  
  ds = ds.where(abs(ds - median) < std * sd_threshold)
  ds = ds.interpolate_na(dim=['y', 'x'], method="linear")
  
  utils.gs_write_blob(os.path.join(path_to_temp, 'temp_average_write.nc'), save_path, bucket)
  utils.image_average_variables(ds, ['GPP', 'NEE', 'NPP', 'Ra', 'Rh'], time_index = 'time', plot_dir='/home/amullen/figs')
  ds.to_zarr('gs://rangelands/Ranch_Runs/MCK/results/flux_hist_weighted_filtered.zarr')

  
#spatial_filter('gs://rangelands/Ranch_Runs/MCK/results/C_stock_hist_weighted.nc', 'Ranch_Runs/MCK/results/C_stock_hist_weighted_filtered.nc', bucket)
spatial_filter('gs://rangelands/Ranch_Runs/MCK/results/flux_hist_weighted.nc', 'Ranch_Runs/MCK/results/flux_hist_weighted_filtered.nc', bucket)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  