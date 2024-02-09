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
from skimage.filters import rank
from skimage.morphology import disk
import argparse
import preprocess_starfm_imagery as psi
import os
import fsspec

path_to_temp = '/home/amullen/temp/'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/remote-sensing/gee_key.json')

bucket_name='rangelands'
path_to_footprint_list = 'res/sites_v2_rest.txt'
bucket = storage_client.get_bucket(bucket_name)

def gen_covariates(starfm_in_dir, daymet_in_dir, out_dir, out_name, bucket_name):
  """reads in starfm imagery from site, calculates ndvi and evi, saves to netcdf

  Args:
    starfm_in_dir (string): path to google storage directory containing starfm images (should not include bucket name, should not start with forward slash)
    daymet_in_dir (string): path to google storage directory containing daymet images (should not include bucket name, should not start with forward slash)
    out_dir (string): path to google storage directory to output covariate files (should not include bucket name, should not start with forward slash)
    bucket_name (string): google storage bucket name

  Returns:
     None
  """ 
  
  bucket = storage_client.get_bucket(bucket_name)
  starfm_df = psi.get_sfm_date_df(starfm_in_dir, bucket)
  daymet_df= psi.get_landsat_date_df(daymet_in_dir, bucket)
  
  dates=[]
  ndvis=[]
  daymets=[]
  
  i=0
  for index, row in starfm_df.iterrows():
    
    if len(daymet_df.loc[daymet_df['im_date']==row['im_date']])>0:
  
      starfm=rxr.open_rasterio('gs://' + bucket.name + '/' + row['im_path'], masked=True)
      starfm.values=starfm.values/10000
      
      ndvi=(starfm.values[1]-starfm.values[0])/(starfm.values[1]+starfm.values[0])
      ndvi[np.isnan(ndvi)] = 0
      
      daymet_index=daymet_df.loc[daymet_df['im_date']==row['im_date']].index[0]
      daymet=rxr.open_rasterio('gs://' + bucket.name + '/' + daymet_df.loc[daymet_index, 'im_path'], masked=True)
      daymet=daymet.rio.reproject_match(starfm)
      daymet.values[0][np.isnan(daymet.values[0])] = 0
      
      dates.append(row['im_date'])
      ndvis.append(ndvi)
      daymets.append(daymet.values[0])
    i+=1
    

  covariate_DataArray = xr.Dataset(data_vars={
                                           'ndvi':(('time', 'y', 'x'),ndvis),
                                           'daymet': (('time', 'y', 'x'),daymets)},
                                  coords={'time':dates,
                                          'y':starfm.y.values,
                                          'x':starfm.x.values 
                                          })
  covariate_DataArray.rio.write_crs(starfm.rio.crs.to_wkt(), inplace=True)
  covariate_DataArray.rio.write_transform(starfm.rio.transform(), inplace=True)
  covariate_DataArray.rio.nodata=0
  covariate_DataArray.to_netcdf(os.path.join(path_to_temp, 'temp_covariates_write.nc'))
  utils.gs_write_blob(os.path.join(path_to_temp, 'temp_covariates_write.nc'), os.path.join(out_dir, out_name), bucket)
  
  return
  
def load_dataset(filename, engine="scipy", *args, **kwargs) -> xr.Dataset:
    """Load a NetCDF dataset from local file system or cloud bucket."""
    with fsspec.open(filename, mode="rb") as file:
        dataset = xr.load_dataset(file, engine=engine, *args, **kwargs)
    return dataset
      
def daymet_qaqc(covariate_nc, bucket_name):
  """QA/QCs evi and ndvi time series by filtering vegetation indices in imagery where daymet temp is <=0 and index value-rolling median value < rolling median * 2(sd). Temporally interpolates na values with nearest temporal    neighbor

  Args:
    covariate_nc (string): path to .nc file containing covariates (currently supports variables 'evi', 'ndvi', and 'daymet')
    bucket_name (string): name of google storage bucket

  Returns:
     xarray.Dataset with QA/QC'ed vegetation indices
  """ 
  bucket = storage_client.get_bucket(bucket_name)
  print('gs://' + bucket_name + '/' + covariate_nc)
  #ds=load_dataset('gs://' + bucket_name + '/' + covariate_nc)
  ds=rxr.open_rasterio('gs://' + bucket_name + '/' + covariate_nc)
  
  print('masking with daymet')
  
  
  ds['rolling_std_ndvi'] = ds['ndvi'].rolling(time=365, center=True, min_periods=1).std()
  ds['rolling_median_ndvi'] = ds['ndvi'].rolling(time=14, center=True, min_periods=1).mean()
  
  ds['qa'] = xr.where((ds['ndvi'] < 0) | (ds['daymet'] < 0) | (abs(ds['ndvi'] - ds['rolling_median_ndvi']) > ds['rolling_median_ndvi']+ds['rolling_std_ndvi']*2), 1, 0)
  
  ds['ndvi'] = ds['ndvi'].where((ds['ndvi'] > 0) & (ds['daymet'] > 0) & (abs(ds['ndvi'] - ds['rolling_median_ndvi']) < ds['rolling_median_ndvi']+ds['rolling_std_ndvi']*2))
  
  ds['ndvi'] =  ds['ndvi'].interpolate_na(dim='time', method='linear', fill_value="extrapolate")
  
  ds['ndvi'] = ds['ndvi'].where((ds['ndvi'] > 0) & (ds['daymet'] > 0) & (abs(ds['ndvi'] - ds['rolling_median_ndvi']) < ds['rolling_median_ndvi']+ds['rolling_std_ndvi']*2))
  
  ds['ndvi'] =  ds['ndvi'].interpolate_na(dim='time', method='nearest', fill_value="extrapolate")
  
  return ds

def image_average_variables(ds, variable_list, plot=None):
  """Take spatial average of xarray Dataset for variables in variable_list. Returns pandas dataframe

  Args:
    ds (xarray.Dataset): dataset with variables to average and dimensions x, y, and time
    variable_list (list of strings): variable names to spatially average
    plot (string): full path to save timeseries plot to, will only plot if this argument is set

  Returns:
     pandas.DataFrame with a column for time, and spatial average of each variable in variable_list
  """ 
  site_dates = pd.to_datetime(pd.DatetimeIndex(ds.indexes['time'].to_datetimeindex()))
  
  ds['ndvi'] = ds['ndvi'].where(ds['ndvi']!=0)
  
  mean_indices=ds[variable_list].mean(dim=['y', 'x']).to_pandas()
  mean_indices = mean_indices.sort_values(by='time', ascending=True)

  if plot!=None:
    sns.lineplot(x=site_dates, y=mean_indices['ndvi'].values, label='ndvi')
    plt.ylabel('index')
    print('saving plot')
    plt.savefig(os.path.join(path_to_temp, 'temp_plot.jpg'), dpi=300)
    utils.gs_write_blob(os.path.join(path_to_temp, 'temp_plot.jpg'), plot, bucket)
    plt.show()
    plt.clf()

  return mean_indices



bad_sites=['Kon']
with open(path_to_footprint_list) as f:
  
  sites = [line.rstrip('\n') for line in f]
  for site in sites:
    if site != '':
      print(site)
      
      site = site.split('/')[-1]    
      if site in bad_sites:
          
          continue
      starfm_in_dir='Ameriflux_sites/{}_starfm/starfm_test_smooth_v2/'.format(site)
      daymet_in_dir='Ameriflux_sites/{}_starfm/daymet/'.format(site)
      out_dir='Ameriflux_sites/{}_starfm/covariates_v2/'.format(site)
      out_name='{}_covariates.nc'.format(site)
      out_csv = '{}_indices_v2.csv'.format(site)
      plot_name='{}_index_plot_QA_spatial_rolling_v2.jpg'.format(site)
      
      gen_covariates(starfm_in_dir, daymet_in_dir, out_dir, out_name, bucket_name)
      ds = daymet_qaqc(os.path.join(out_dir, out_name), bucket_name)
      df = image_average_variables(ds, ['ndvi', 'daymet', 'qa'], plot=out_dir+plot_name)
  
      df.to_csv('gs://' + bucket_name + '/' + out_dir + out_csv)

  