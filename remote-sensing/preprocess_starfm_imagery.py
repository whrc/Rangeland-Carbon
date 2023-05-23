from google.cloud import storage
import rioxarray as rxr
import os
import utils
import pandas as pd
import numpy as np
from itertools import chain
from scipy.stats import t
import xarray as xr
from matplotlib import pyplot as plt
import seaborn as sns


storage_client = storage.Client.from_service_account_json('gee_key.json')

bucket='rangelands'
in_landsat_dir = ''
out_landsat_dir = ''
in_modis_dir = 'Ameriflux_sites/Rws_starfm/modis_test/'
out_modis_dir = 'Ameriflux_sites/Rws_starfm/modis_test_smooth/'
path_to_temp = '/home/amullen/temp/'

bucket = storage_client.get_bucket(bucket)

#########
# skip landsat images with all zero
# gap-fill landsat
# 

########

########
#fliter modis collection

def get_modis_date_df(in_modis_dir, bucket):
  """gets dataframe of MODIS images with date in Y_m_d format and milliseconds, image name, and full path

  Args:
    in_modis_dir (string): path to google storage directory to list (should not include bucket name, should not start with forward slash)
    bucket (storage.bucket): google cloud bucket

  Returns:
     pandas dataframe
  """ 
  im_paths = utils.gs_listdir(in_modis_dir, bucket)
  im_paths = [k for k in im_paths if '.tif' in k]
  
  df_im_paths = pd.DataFrame({'im_path':im_paths})
  df_im_paths['im_name'] = [i[-1] for i in df_im_paths['im_path'].str.split('/').to_list()]
  df_im_paths['im_date'] = pd.to_datetime(df_im_paths['im_path'].str[-14:-4], format='%Y_%m_%d')
  df_im_paths['millis'] = df_im_paths['im_date'].astype(np.int64) / int(1e6)
  df_im_paths = df_im_paths.set_index(df_im_paths['im_date'])
  df_im_paths = df_im_paths.sort_index()
  
  return df_im_paths
  
def get_landsat_date_df(in_landsat_dir, bucket):
  """gets dataframe of MODIS images with date in Y_m_d format and milliseconds, image name, and full path

  Args:
    in_modis_dir (string): path to google storage directory to list (should not include bucket name, should not start with forward slash)
    bucket (storage.bucket): google cloud bucket

  Returns:
     pandas dataframe
  """ 
  im_paths = utils.gs_listdir(in_landsat_dir, bucket)
  im_paths = [k for k in im_paths if '.tif' in k]
  
  df_im_paths = pd.DataFrame({'im_path':im_paths})
  df_im_paths['im_name'] = [i[-1] for i in df_im_paths['im_path'].str.split('/').to_list()]
  df_im_paths['im_date'] = pd.to_datetime(df_im_paths['im_path'].str[-12:-4], format='%Y%m%d')
  df_im_paths['millis'] = df_im_paths['im_date'].astype(np.int64) / int(1e6)
  df_im_paths = df_im_paths.set_index(df_im_paths['im_date'])
  df_im_paths = df_im_paths.sort_index()
  
  return df_im_paths
 
def lag_linregress_3D(x, y, lagx=0, lagy=0):
  """
  Input: Two xr.Datarrays of any dimensions with the first dim being time. 
  Thus the input data could be a 1D time series, or for example, have three 
  dimensions (time,lat,lon). 
  Datasets can be provided in any order, but note that the regression slope 
  and intercept will be calculated for y with respect to x.
  Output: Covariance, correlation, regression slope and intercept, p-value, 
  and standard error on regression between the two datasets along their 
  aligned time dimension.  
  Lag values can be assigned to either of the data, with lagx shifting x, and
  lagy shifting y, with the specified lag amount. 
  """ 
  #1. Ensure that the data are properly alinged to each other. 

  x,y = xr.align(x,y)

  #2. Add lag information if any, and shift the data accordingly
  if lagx!=0:
  
      # If x lags y by 1, x must be shifted 1 step backwards. 
      # But as the 'zero-th' value is nonexistant, xr assigns it as invalid 
      # (nan). Hence it needs to be dropped
      x   = x.shift(time = -lagx).dropna(dim='time')
  
      # Next important step is to re-align the two datasets so that y adjusts
      # to the changed coordinates of x
      x,y = xr.align(x,y)
  
  if lagy!=0:
      y   = y.shift(time = -lagy).dropna(dim='time')
      x,y = xr.align(x,y)
  
  #3. Compute data length, mean and standard deviation along time axis: 
  n = y.notnull().sum(dim='time')
  xmean = x.mean(axis=0)
  ymean = y.mean(axis=0)
  xstd  = x.std(axis=0)
  ystd  = y.std(axis=0)

  #4. Compute covariance along time axis
  cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)
  
  #5. Compute correlation along time axis
  cor   = cov/(xstd*ystd)
  
  #6. Compute regression slope and intercept:
  slope     = cov/(xstd**2)
  intercept = ymean - xmean*slope  
  
  #7. Compute P-value and standard error
  #Compute t-statistics
  tstats = cor*np.sqrt(n-2)/np.sqrt(1-cor**2)
  stderr = slope/tstats

  #pval   = t.sf(tstats, n-2)*2
  #pval   = xr.DataArray(pval, dims=cor.dims, coords=cor.coords)
  
  return cov,cor,slope,intercept,stderr 

def smooth_modis_col(in_modis_dir, bucket, windowsize='20d', min_periods=1):
  
  """smooths MODIS collection using linear regression over a moving window

  Args:
    in_modis_dir (string): path to google storage directory containing only modis images (should not include bucket name, should not start with forward slash)
    bucket (storage.bucket): google cloud bucket
    windowsize (string): time period for moving window
    min_periods (int): minimum periods within the moving window for the window to be valid

  Returns:
     pandas none
  """ 
  
  df_im_paths = get_modis_date_df(in_modis_dir, bucket)

  #build rolling window for image names
  rolling_dfs = list(df_im_paths.rolling(windowsize, center=True))
  df_im_paths = df_im_paths.reset_index(drop=True)

  #create variables to hold image data for moving window
  temp_im_dict = {}
  temp_time_dict = {}
  
  blue=[]
  blue_smooth=[]
  green=[]
  green_smooth=[]
  red=[]
  red_smooth=[]
  nir=[]
  nir_smooth=[]
  
  for index, row in df_im_paths.iterrows():
    
    #target image is image being smoothed
    tgt_image = rxr.open_rasterio('gs://' + bucket.name + '/' + row['im_path'])
    tgt_image = tgt_image.astype(np.float32)
    tgt_image.values[tgt_image.values==0] = np.nan
    tgt_image_name = row['im_path'].split('/')[-1]

    #create array to hold 'x' for prediction after linear regression. 'x' in this case is time (milliseconds)
    tgt_millis = np.full(np.shape(tgt_image[0:5]), row['millis']) 
    
    #on each iteration, determine which files are being opened and which files are being closed
    rolling_df = rolling_dfs[index]

    if index == 0:
      open_rows = rolling_df
      temp_images = list(open_rows['im_path'])
      close_rows = pd.DataFrame(columns=open_rows.columns)
      
    else:
      open_rows = rolling_df[~rolling_df['im_path'].isin(temp_images)]
      temp_images = temp_images + list(open_rows['im_path'])
      close_rows = rolling_dfs[index-1][~rolling_dfs[index-1]['im_path'].isin(rolling_df['im_path'])]
      temp_images = [x for x in temp_images if x not in list(close_rows['im_path'])]
    
    #open all rows that need to be open, add to image dict
    for open_row_index, open_row in open_rows.iterrows():
      if open_row['im_path'] not in temp_im_dict:
      
        temp_im_dict[open_row['im_path']] = rxr.open_rasterio('gs://' + bucket.name + '/' + open_row['im_path'])[0:5]
        
        if temp_im_dict[open_row['im_path']].rio.transform() != tgt_image.rio.transform():
          
          print('reprojecting image')
          temp_im_dict[open_row['im_path']] = temp_im_dict[open_row['im_path']].rio.reproject_match(tgt_image)
          
        temp_im_dict[open_row['im_path']] = temp_im_dict[open_row['im_path']].astype(np.float32)
        temp_im_dict[open_row['im_path']].values[temp_im_dict[open_row['im_path']].values==0]=np.nan
        time = np.full(np.shape(temp_im_dict[open_row['im_path']]), open_row['millis'])
 
        temp_time_dict[open_row['im_path']] = xr.DataArray(data=time, coords=temp_im_dict[open_row['im_path']].coords)
        
    #close datasets that need to be closed, remove from dict
    for close_row_index, close_row in close_rows.iterrows():
      temp_im_dict[close_row['im_path']].close()
      temp_time_dict[close_row['im_path']].close()
      temp_im_dict.pop(close_row['im_path'])
      temp_time_dict.pop(close_row['im_path'])
    
    #stack images over new dimension "time" and run linear regression
    image_stack = xr.concat(list(temp_im_dict.values()), dim='time')
    time_stack = xr.concat(list(temp_time_dict.values()), dim='time')
    cov,cor,slope,intercept,stderr = lag_linregress_3D(time_stack, image_stack)
    
    #predict on target image
    for i in range(0,5):
    
      tgt_image.values[i] = tgt_millis[i]*slope[i]+intercept[i]
      tgt_image.values[i][np.isnan(tgt_image.values[i])] = 0
    
    if np.nanmax(tgt_image.values[0])==0:
      continue

    tgt_image.rio.write_nodata(0, inplace=True)

    tgt_image.rio.to_raster(os.path.join(path_to_temp, 'temp_gs_write.tif'), dtype=np.uint16, driver='GTiff')
    utils.gs_write_blob(os.path.join(path_to_temp, 'temp_gs_write.tif'), os.path.join(out_modis_dir, tgt_image_name), bucket)
      
  return


#smooth_modis_col(in_modis_dir, bucket)