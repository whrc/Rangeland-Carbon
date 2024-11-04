import xarray as xr
import rioxarray as rxr
import numpy as np
from google.cloud import storage
import sys
sys.path.insert(1, '../modeling')
sys.path.insert(1, '../utils')
import utils
import RCTM_utils as rctmu
import os
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import fsspec
import pandas as pd
import xarrayMannKendall as mk
from osgeo import gdal, gdalconst
import time
from dask.distributed import Client

def export_ds_to_geotiff_gcloud(ds, out_path, path_to_temp):

  ds.rio.to_raster(os.path.join(path_to_temp, 'temp_gtiff_write.tif'))
  utils.gs_write_blob(os.path.join(path_to_temp, 'temp_gtiff_write.tif'), out_path, bucket)
  os.remove(os.path.join(path_to_temp, 'temp_gtiff_write.tif'))
  
  return

def export_cumulative_flux(covariates, out_path, path_to_temp, days = 5):
  
  crs = covariates.rio.crs.to_wkt()
  transform = covariates.rio.transform()
  
  covariates = covariates*days
  covariates = covariates.sum(dim='time')
  
  covariates = covariates.fillna(0)
  covariates.rio.nodata=0
  covariates.rio.write_crs(crs, inplace=True)
  covariates.rio.write_transform(transform, inplace=True)
  
  #save to temp, write to google storage
  export_ds_to_geotiff_gcloud(covariates, out_path, path_to_temp)
  
  return covariates
  
def export_average(ds, out_path, path_to_temp):

  crs = ds.rio.crs.to_wkt()
  transform = ds.rio.transform()
  
  ds = ds.mean(dim='time')
  
  ds = ds.fillna(0)
  
  ds.rio.write_crs(crs, inplace=True)
  ds.rio.write_transform(transform, inplace=True)
  ds.rio.write_nodata(0, inplace=True) 
  #save to temp, write to google storage
  export_ds_to_geotiff_gcloud(ds, out_path, path_to_temp)
  
  return ds
  
def spatial_stats_by_value(ds, pixel_area, ancillary_ds, value, variables):
  
  dates = ds.indexes['time'].to_datetimeindex()
  
  masked_ds  = ds.where(ancillary_ds[0]==value)
  average = masked_ds.mean(dim=['x', 'y'])
  sum = masked_ds.sum(dim=['x', 'y'])
  std = masked_ds.std(dim=['x', 'y'])
  
  df = pd.DataFrame({'Date': dates})
  df['value'] = value
  
  for variable in variables:
    df[variable+'_mean_gC/m2'] = average[variable].values #gC/m2
    df[variable+'_sum_gC'] = sum[variable].values * pixel_area #gC
    df[variable+'_std_gC'] = std[variable].values * pixel_area #gC
  
  return df
  
def export_change_analysis(ds, out_path, path_to_temp, out_p_val, out_signif):
  
  crs = ds.rio.crs.to_wkt()
  transform = ds.rio.transform()
  
  yearly = ds.groupby('time.year').sum()
  
  MK_class = mk.Mann_Kendall_test(yearly, coords_name = {'time':'year','y':'y','x':'x'})
  
  change_analysis = MK_class.compute(progress_bar=True)
  change_analysis.rio.write_crs(crs, inplace=True)
  change_analysis.rio.write_transform(transform, inplace=True)
  change_analysis.rio.nodata=0
  
  export_ds_to_geotiff_gcloud(change_analysis['trend'], out_path, path_to_temp)
  export_ds_to_geotiff_gcloud(change_analysis['p'], out_p_val, path_to_temp)
  export_ds_to_geotiff_gcloud(change_analysis['signif'], out_signif, path_to_temp)
  
  return change_analysis




################################################ C_stocks #######################################################################


def spatial_analysis_C_stocks(C_stocks_path, path_to_landcover, ranch):

  C_stocks = rxr.open_rasterio(C_stocks_path, masked=True)
  landcover = rxr.open_rasterio(path_to_landcover)
  
  landcover_res = landcover.rio.reproject_match(C_stocks)
  landcover_metric = landcover_res.rio.reproject('EPSG:5070')
  gdal_transform = landcover_metric.rio.transform().to_gdal()
  pixel_area_m2 = np.abs(gdal_transform[1]) * np.abs(gdal_transform[5])
  
  C_stocks=C_stocks.sel(time = (C_stocks['time.year'] != 2023) & (C_stocks['time.year'] != 2022))
    
  C_stocks['TSOC'] = C_stocks['POC'] + C_stocks['HOC'] + 1000 #in gC/m2
    
  #Change in C stocks
  C_stocks_2022 = export_average(C_stocks['TSOC'].sel(time=(C_stocks['TSOC']['time.year'] == 2021)), f'Ranch_Runs/{ranch}/analysis/C_stock_average_2021.tif', path_to_temp)
  C_stocks_2003 = export_average(C_stocks['TSOC'].sel(time=(C_stocks['TSOC']['time.year'] == 2003)), f'Ranch_Runs/{ranch}/analysis/C_stock_average_2003.tif', path_to_temp)
    
  export_ds_to_geotiff_gcloud(C_stocks_2022 - C_stocks_2003, f'Ranch_Runs/{ranch}/analysis/C_stock_average_change_2021-2003.tif', path_to_temp)
  
  print('averaging/plotting C stocks')
  C_stock_site_dfs=[]
  
  print('averaging grass')
  C_stock_grass_df = spatial_stats_by_value(C_stocks, pixel_area_m2, landcover_res, 1, ['TSOC'])
  C_stock_site_dfs.append(C_stock_grass_df)
      
  print('averaging shrubs')
  C_stock_shrub_df = spatial_stats_by_value(C_stocks, pixel_area_m2, landcover_res, 2, ['TSOC'])
  C_stock_site_dfs.append(C_stock_shrub_df)
      
  print('averaging trees')
  C_stock_tree_df = spatial_stats_by_value(C_stocks, pixel_area_m2, landcover_res, 3, ['TSOC'])
  C_stock_site_dfs.append(C_stock_tree_df)
  
  C_stock_site_dfs=pd.concat(C_stock_site_dfs)
  C_stock_site_dfs.to_csv(f'gs://rangelands/Ranch_Runs/{ranch}/analysis/C_stock_timeseries.csv')

##################################################################################################################################

#Fluxes
def spatial_analysis_fluxes(covariates_path, path_to_landcover, ranch):
  covariates = rxr.open_rasterio(covariates_path, masked=True)
  covariates= covariates.sel(time = (covariates['time.year'] != 2023) & (covariates['time.year'] != 2022))
  
  landcover = rxr.open_rasterio(path_to_landcover)
  landcover_res = landcover.rio.reproject_match(covariates)
  landcover_metric = landcover_res.rio.reproject('EPSG:5070')
  gdal_transform = landcover_metric.rio.transform().to_gdal()
  pixel_area_m2 = np.abs(gdal_transform[1]) * np.abs(gdal_transform[5])
  
  #Change Analysis
  GPP_change = export_change_analysis(covariates['GPP']*5, out_path=f'Ranch_Runs/{ranch}/analysis/GPP_change_analysis.tif', path_to_temp = path_to_temp, out_p_val = f'Ranch_Runs/{ranch}/analysis/GPP_change_analysis_p_value.tif',
                                        out_signif=f'Ranch_Runs/{ranch}/analysis/GPP_change_analysis_signif.tif')
                                        
  #Mean annual GPP/Reco
  covariates['Reco'] = covariates['Ra'] + covariates['Rh']
  GPP_RECO_means = (covariates[['GPP', 'Reco']]*5).groupby('time.year').sum().mean(dim='year')
  export_ds_to_geotiff_gcloud(GPP_RECO_means['GPP'], out_path=f'Ranch_Runs/{ranch}/analysis/mean_annual_GPP.tif', path_to_temp = path_to_temp)
  export_ds_to_geotiff_gcloud(GPP_RECO_means['Reco'], out_path=f'Ranch_Runs/{ranch}/analysis/mean_annual_Reco.tif', path_to_temp = path_to_temp)
  
  #covariates = utils.xr_dataset_to_data_array(covariates)
  print('averaging/plotting fluxes')
  site_dfs=[]
  
  print('averaging grass')
  grass_df = spatial_stats_by_value(covariates, pixel_area_m2, landcover_res, 1, ['GPP', 'Reco', 'NEE'])
  site_dfs.append(grass_df)
      
  print('averaging shrubs')
  shrub_df = spatial_stats_by_value(covariates, pixel_area_m2, landcover_res, 2, ['GPP', 'Reco', 'NEE'])
  site_dfs.append(shrub_df)
      
  print('averaging trees')
  tree_df = spatial_stats_by_value(covariates, pixel_area_m2, landcover_res, 3, ['GPP', 'Reco', 'NEE'])
  site_dfs.append(tree_df)
  
  site_dfs=pd.concat(site_dfs)
  site_dfs.to_csv(f'gs://rangelands/Ranch_Runs/{ranch}/analysis/flux_timeseries.csv')

def spatial_average_dask(ds, vars):
  #calculates spatial average

  ds = ds[vars]
  mean = ds.mean(dim=['x', 'y'])
  
  return mean
#utils.image_average_variables(covariates, ['GPP', 'Reco', 'Ra', 'NEE', 'NPP'], plot_dir='output/RCTM/transient')


####################################################################################################################
#roi_file='/home/amullen/Rangeland-Carbon/res/site_footprints/HLD_tiles.txt'

#with open(roi_file) as f:
#  paths = [line.rstrip('\n') for line in f]
  
#  for path in paths:
    
#    ranch='/'.join(path.split('/')[-2:])
#    print(ranch)
#    if ranch in ['...']:
#      continue
   
#C_stocks_path = f'gs://rangelands/Ranch_Runs/{ranch}/results/C_stock_hist_weighted.nc'
#covariates_path = f'gs://rangelands/Ranch_Runs/{ranch}/results/flux_hist_weighted.nc'
#path_to_landcover = 'gs://' + bucket.name + f'/Ranch_Runs/{ranch}/landcover/fused_landcover.tif'
#    
#spatial_analysis_C_stocks(C_stocks_path, path_to_landcover, ranch)
#spatial_analysis_fluxes(covariates_path, path_to_landcover, ranch)

#google storage etc.
if __name__ == '__main__':
  
  bucket_name='rangelands'
  storage_client = storage.Client.from_service_account_json('/home/amullen/res/gee_key.json')
  bucket = storage_client.get_bucket(bucket_name)
  path_to_temp = '/home/amullen/temp/'

  #initiate timer
  start_time = time.time()

  with Client(n_workers=15, threads_per_worker=1, memory_limit='1.5GB') as client:
  #if True:
    #open C stock file
    C_stocks_path = f'gs://rangelands/Ranch_Runs/MCK/results/C_stock_hist_weighted.nc'
    C_stocks = rxr.open_rasterio(C_stocks_path, chunks='auto')
    
    #3 workers
    #auto: 1028s
    #time = 100: DNC
    #time = 500:
    #time = 1000: 
    #no chunks: 437s

    #get doy for each datetime index in file
    doys = np.array(C_stocks.indexes['time'].dayofyear)

    #create index to bin individual timestamps into 5-day intervals
    doy_bins_lower = (np.floor((doys-1)/5).astype(np.int16)*5)+1

    #assign coordinates to doy bins
    C_stocks = C_stocks.assign_coords(doy_bins_lower=('time', doy_bins_lower))
    C_stocks = C_stocks.swap_dims({'time':'doy_bins_lower'})

    #mean = C_stocks.groupby('doy_bins_lower').mean(dim=['x', 'y'])
    #std = C_stocks.groupby('doy_bins_lower').std(dim=['x', 'y'])
    mean = C_stocks.mean(dim=['x', 'y'])
    std = C_stocks.std(dim=['x', 'y'])
    sum = C_stocks.sum(dim=['x', 'y'])

    mean_df = mean.to_dataframe()
    std_df = std.to_dataframe()
    sum_df = sum.to_dataframe()
    mean_df.to_csv('MCK_C_stock_means.csv')
    std_df.to_csv('MCK_C_stock_stds.csv')
    sum_df.to_csv('MCK_C_stock_sums.csv')
    
  #end timer
  end_time = time.time()
  elapsed_time = end_time - start_time
  print(f"Elapsed time: {elapsed_time:.6f} seconds")
