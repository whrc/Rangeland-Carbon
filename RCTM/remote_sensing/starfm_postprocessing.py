from google.cloud import storage
import rioxarray as rxr
import os
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
import os
import fsspec
import sys
from RCTM.utils.utils import gs_write_blob, get_sfm_date_df, get_landsat_date_df
import yaml

############################################ starfm_postprocessing.py #############################################################
# This file contains the processes necesarry to package and postprocess starfm imagery and covariates to be fed into the RCTM model.
# 
# There are 4 major processes that occur
# 1. Fetch STARFM and covariate images, and place into a single dataset with a consistent timestep (currently 5 days). 
#    - function(s): gen_covariates()
# 2. QA/QC the STARFM retrievals based on NDVI and DAYMET air temperature, filter as needed, and add QA band
#    - function(s): daymet_qaqc()
# 3. Temporal gap-filling of NDVI and covariates, currently a simple nearest-neighbor approximation 
#    - function(s): daymet_qaqc()
# 4. Generate spinup dataset
#    - function(s): aggregate_for_spinup()
# 5. Generate an image of spatially distributed RCTM parameters based on plant-functional type that is determined by landcover. Currently based on NLCD class.
#    - function(s): get_spatial_RCTM_params()
# 
# Additional useful functions:
#  - get_ndvi_nc(): builds netcdf time series with only NDVI calculated from STARFM. Useful for debugging image fusion methods.
#  - image_average_variables(): spatially averages one or multiple variables in an XArray dataset and saves to a pandas dataframe. 
#    Useful for assessing covariate time series and running the model in point mode.
#
#Example usage:
# The following call will run all 5 processes sequentially
#
#python starfm_postprocessing.py --starfm_in_dir='Ranch_Runs/HV/starfm/' --covariate_in_dir='Ranch_Runs/HV/covariates/' --out_dir='Ranch_Runs/HV/covariates_nc/' --NLCD_in_dir='Ranch_Runs/HV/landcover/NLCD_2019.tif' --RAP_in_dir='Ranch_Runs/HV/landcover/RAP_2019.tif' --param_file='/home/amullen/Rangeland-Carbon/modeling/RCTM_params.yaml' --bucket_name='rangelands'
#
###################################################################################################################################


#path_to_temp = '/home/amullen/temp/'
#storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')

#bucket_name='rangelands'
#path_to_footprint_list = '/home/amullen/Rangeland-Carbon/res/site_footprints/HLR_tiles_subregion.txt'
#bucket = storage_client.get_bucket(bucket_name)

  
def get_ndvi_nc(starfm_in_dir, out_dir, out_name, bucket_name):
  """Depracated
      reads in starfm imagery from site, calculates ndvi, saves time series to netcdf

  Args:
    starfm_in_dir (string): path to google storage directory containing starfm images (should not include bucket name, should not start with forward slash)
    daymet_in_dir (string): path to google storage directory containing daymet images (should not include bucket name, should not start with forward slash)
    out_dir (string): path to google storage directory to output covariate files (should not include bucket name, should not start with forward slash)
    bucket_name (string): google storage bucket name

  Returns:
     None
  """ 
  
  bucket = storage_client.get_bucket(bucket_name)
  starfm_df = get_sfm_date_df(starfm_in_dir, bucket).reset_index(drop=True)
  
  ref_im = rxr.open_rasterio('gs://' + bucket.name + '/' + starfm_df.loc[0, 'im_path'], masked=True)
  print(ref_im)
  dates=[]
  ndvis=[]
  
  #create 5 day date range starting on 2002-01-01
  date_range = pd.DataFrame({'im_date': pd.date_range(start='2002-01-01', end='2022-12-31', freq='5d').to_series()}).reset_index(drop=True)
  starfm_df = date_range.merge(starfm_df, how='left', on='im_date')
  
  #iterate over dates to generate covariate timeseries
  for index, row in starfm_df.iterrows():
  
    dates.append(row['im_date'])
    print(row['im_path'])
    #STARFM results
    if not pd.isnull(row['im_path']):
      starfm=rxr.open_rasterio('gs://' + bucket.name + '/' + row['im_path'], masked=True)
      starfm=starfm.rio.reproject_match(ref_im)
      
      ndvi=(starfm.values[1]-starfm.values[0])/(starfm.values[1]+starfm.values[0])
      ndvis.append(ndvi)
      
    else:
      empty_arr=np.empty(np.shape(ref_im.values[0]))
      empty_arr[:] = np.nan
      ndvis.append(empty_arr)
      
  covariate_DataArray = xr.Dataset(data_vars={'ndvi':(('time', 'y', 'x'),ndvis)},
                                  coords={'time':dates,
                                          'y':starfm.y.values,
                                          'x':starfm.x.values 
                                          })
  
  covariate_DataArray.rio.write_crs(starfm.rio.crs.to_wkt(), inplace=True)
  covariate_DataArray.rio.write_transform(starfm.rio.transform(), inplace=True)
  covariate_DataArray.rio.nodata=0
  covariate_DataArray.to_netcdf(os.path.join(path_to_temp, 'temp_covariates_write.nc'))
  gs_write_blob(os.path.join(path_to_temp, 'temp_covariates_write.nc'), os.path.join(out_dir, out_name), bucket)
  
  return covariate_DataArray
  
def get_spatial_RCTM_params(NLCD_in_dir, RAP_in_dir, params, fused_landcover_outname, param_out_name, bucket, path_to_temp, param_type='starfm'):
  """Generates n-band geotiff holding the spatial parameters for RCTM, each band is a parameter. 
     Parameterization is based on plant functional type which are currently determined using NLCD.

  Args:
    landcover_in_dir (string): path to google storage directory containing starfm images (should not include bucket name, should not start with forward slash)
    param_file (string): path to parameter file (.yaml) that holds pft-specific RCTM parameters.
    param_out_name (string): path to output parameter geotiff on google storage
    bucket_name (string): google storage bucket name
    param_type (string): which parameter set to reference. Currently the RCTM_params.yaml file only contains starfm params.

  Returns:
     None
  """ 
  #open landcover
  NLCD_landcover = rxr.open_rasterio('gs://' + bucket.name + '/' + NLCD_in_dir, masked=True)
  RAP_landcover = rxr.open_rasterio('gs://' + bucket.name + '/' + RAP_in_dir, masked=True)
  RAP_landcover = RAP_landcover.rio.reproject_match(NLCD_landcover)
  #RAP_landcover = utils.xr_dataset_to_data_array(RAP_landcover)
  
  #get keys (param names), doesn't matter which pft since they all have the same parameters
  keys = params['RCTM_params']['grassland'][param_type].keys()
  
  #create numpy array of nans with proper shape, shape = (params, y, x)
  spatial_params = np.zeros((len(keys), np.shape(NLCD_landcover.values)[1], np.shape(NLCD_landcover.values)[2]))
  spatial_params[:] = np.nan
  
  landcover = np.zeros(np.shape(NLCD_landcover.values[0]))
  
  #RAP Thresholds
  RAP_thresholds = {
                    'grass-tree': 30,
                    'grass-shrub': 30,
                    'grassland': 30
                    }
                    
  #NLCD pft mapping
  NLCD_LUT = {
            '41': 'grass-tree',
            '42': 'grass-tree',
            '43': 'grass-tree',
            '51': 'grass-shrub',
            '52': 'grass-shrub',
            '71': 'grassland',
            '72': 'grassland',
            '81': 'grassland',
            '82': 'grassland',
            '90': 'grass-shrub',
            '95': 'grassland'
            }
            
  NLCD_grass = [71, 72, 81, 82, 90, 95]
  NLCD_shrub = [51, 52, 90]
  NLCD_tree = [41, 42, 43]
  
  landcover[np.isin(NLCD_landcover.values[0], NLCD_grass)] = 1
  landcover[np.isin(NLCD_landcover.values[0], NLCD_shrub)] = 2
  landcover[np.isin(NLCD_landcover.values[0], NLCD_tree)] = 3
  
  landcover[RAP_landcover.values[0]>=RAP_thresholds['grassland']] = 1
  landcover[RAP_landcover.values[3]>=RAP_thresholds['grassland']] = 1
  landcover[RAP_landcover.values[4]>=RAP_thresholds['grass-shrub']] = 2
  landcover[RAP_landcover.values[5]>=RAP_thresholds['grass-tree']] = 3
  
  #loop through params and pfts setting spatial params
  for i, key in enumerate(keys):
    spatial_params[i][landcover==1] = params['RCTM_params']['grassland'][param_type][key]
    spatial_params[i][landcover==2] = params['RCTM_params']['grass-shrub'][param_type][key]
    spatial_params[i][landcover==3] = params['RCTM_params']['grass-tree'][param_type][key]
  
  #create DataArray for new landcover
  landcover_DataArray = xr.DataArray(data = landcover,
                                      dims = ['y', 'x'],
                                      coords={'y':NLCD_landcover.y.values,
                                              'x':NLCD_landcover.x.values 
                                              },
                                      attrs=dict(description="landcover"))
                                      
  #write spatial information and set nodata value to 0                                 
  landcover_DataArray.rio.write_crs(NLCD_landcover.rio.crs.to_wkt(), inplace=True)
  landcover_DataArray.rio.write_transform(NLCD_landcover.rio.transform(), inplace=True)
  landcover_DataArray.rio.write_nodata(0, inplace=True)
  
  #save to temp, write to google storage
  landcover_DataArray.rio.to_raster(os.path.join(path_to_temp, 'temp_landcover_write.tif'))
  gs_write_blob(os.path.join(path_to_temp, 'temp_landcover_write.tif'), fused_landcover_outname, bucket)
  
  #create DataArray with spatial parameters
  param_DataArray = xr.DataArray(data = spatial_params,
                                      dims = ['band', 'y', 'x'],
                                      coords={'y':NLCD_landcover.y.values,
                                              'x':NLCD_landcover.x.values 
                                              },
                                      attrs=dict(description="landcover", long_name=list(keys)))
                                      
  #write spatial information and set nodata value to 0                                 
  param_DataArray.rio.write_crs(NLCD_landcover.rio.crs.to_wkt(), inplace=True)
  param_DataArray.rio.write_transform(NLCD_landcover.rio.transform(), inplace=True)
  param_DataArray.rio.write_nodata(0, inplace=True)
  
  #save to temp, write to google storage
  param_DataArray.rio.to_raster(os.path.join(path_to_temp, 'temp_spatial_params_write.tif'))
  gs_write_blob(os.path.join(path_to_temp, 'temp_spatial_params_write.tif'), param_out_name, bucket)
  
def load_dataset(filename, engine="scipy", *args, **kwargs) -> xr.Dataset:
    """Load a NetCDF dataset from local file system or cloud bucket."""
    with fsspec.open(filename, mode="rb") as file:
        dataset = xr.load_dataset(file, engine=engine, *args, **kwargs)
    return dataset
      
def daymet_qaqc(ds, bucket):
  """QA/QCs evi and ndvi time series by filtering vegetation indices in imagery where daymet temp is <=0 and index value-rolling median value < rolling median * 2(sd). Temporally interpolates na values with nearest temporal neighbor

  Args:
    covariate_nc (string): path to .nc file containing covariates (currently supports variables 'evi', 'ndvi', and 'daymet')
    bucket_name (string): name of google storage bucket

  Returns:
     xarray.Dataset with QA/QC'ed vegetation indices, QA band, gap-filled
  """ 
  if not (isinstance(ds, xr.core.dataarray.DataArray) | isinstance(ds, xr.core.dataarray.Dataset)):
    #ds=load_dataset('gs://' + bucket_name + '/' + covariate_nc)
    ds=rxr.open_rasterio('gs://' + bucket.name + '/' + ds)
  
  print('masking with daymet')
  
  print('rolling statistics')
  ds['rolling_std_ndvi'] = ds['ndvi'].rolling(time=365, center=True, min_periods=1).std()
  ds['rolling_median_ndvi'] = ds['ndvi'].rolling(time=14, center=True, min_periods=1).mean()
  
  print('qa band from rolling statistics')
  
  ds['qa'] = xr.where(ds['ndvi'] < 0, 1, 0)
  ds['qa'] = xr.where(ds['tavg'] < 0, 1, 0)
  ds['qa'] = xr.where((abs(ds['ndvi'] - ds['rolling_median_ndvi']) > ds['rolling_median_ndvi']+ds['rolling_std_ndvi']*2), 1, 0)
  
  print('checking NDVI validity')
  ds['ndvi'] = ds['ndvi'].where((ds['ndvi'] > 0) & (ds['tavg'] > 0) & (abs(ds['ndvi'] - ds['rolling_median_ndvi']) < ds['rolling_median_ndvi']+ds['rolling_std_ndvi']*2))
  
  ds['ndvi'] =  ds['ndvi'].interpolate_na(dim='time', method='linear', fill_value="extrapolate")
  
  ds['ndvi'] = ds['ndvi'].where((ds['ndvi'] > 0) & (ds['tavg'] > 0) & (abs(ds['ndvi'] - ds['rolling_median_ndvi']) < ds['rolling_median_ndvi']+ds['rolling_std_ndvi']*2))
  
  ds['ndvi'] =  ds['ndvi'].interpolate_na(dim='time', method='nearest', fill_value="extrapolate")
  
  print('interpolating')
  ds =  ds.interpolate_na(dim='time', method='nearest', fill_value="extrapolate")
  
  return ds

def gen_covariates(starfm_in_dir, daymet_in_dir,  out_dir, out_name, bucket, start_date, end_date, path_to_temp, gap_fill = True):
  """reads in starfm imagery from site, calculates ndvi and evi, saves to netcdf

  Args:
    starfm_in_dir (string): path to google storage directory containing starfm images (should not include bucket name, should not start with forward slash)
    daymet_in_dir (string): path to google storage directory containing daymet images (should not include bucket name, should not start with forward slash)
    out_dir (string): path to google storage directory to output covariate files (should not include bucket name, should not start with forward slash)
    bucket_name (string): google storage bucket name

  Returns:
     None
  """ 
  
  starfm_df = get_sfm_date_df(starfm_in_dir, bucket).reset_index(drop=True)
  daymet_df= get_landsat_date_df(daymet_in_dir, bucket)
  
  ref_im = rxr.open_rasterio('gs://' + bucket.name + '/' + starfm_df.loc[0, 'im_path'], masked=True)
  
  dates=[]
  ndvis=[]
  daymets=[]
  
  #create 5 day date range starting on 2002-01-01
  date_range = pd.DataFrame({'im_date': pd.date_range(start=start_date, end=end_date, freq='5d').to_series()}).reset_index(drop=True)
  starfm_df = date_range.merge(starfm_df, how='left', on='im_date')
  
  #iterate over dates to generate covariate timeseries
  for index, row in starfm_df.iterrows():
  
    dates.append(row['im_date'])
    print(row['im_path'])
    print(row['im_date'])
    #STARFM results
    if not pd.isnull(row['im_path']):
      starfm=rxr.open_rasterio('gs://' + bucket.name + '/' + row['im_path'], masked=True)
      starfm=starfm.rio.reproject_match(ref_im)
      starfm.values=starfm.values/10000
      
      ndvi=(starfm.values[1]-starfm.values[0])/(starfm.values[1]+starfm.values[0])
      ndvis.append(ndvi)
      
    else:
      fill = np.empty(np.shape(ref_im.values[0]))
      fill[:]=np.nan
      ndvis.append(fill)
    
    if len(daymet_df.loc[daymet_df['im_date']==row['im_date']])>0:
      
      daymet_index=daymet_df.loc[daymet_df['im_date']==row['im_date']].index[0]
      daymet=rxr.open_rasterio('gs://' + bucket.name + '/' + daymet_df.loc[daymet_index, 'im_path'], masked=True)
  
      if len(daymet.attrs['long_name'])!=11:
        print('missing covariate')
        cov_fill = np.empty((11, np.shape(ref_im.values[0])[0],np.shape(ref_im.values[0])[1]))
        
        cov_fill[:] = np.nan
        daymets.append(cov_fill)
        continue
        
      daymet=daymet.rio.reproject_match(ref_im)
      daymets.append(daymet)
      
    else:
      print('missing image')
      cov_fill = np.empty((11, np.shape(ref_im.values[0])[0],np.shape(ref_im.values[0])[1]))
      cov_fill[:] = np.nan
      daymets.append(cov_fill)
    
  print('building dataset')
  daymets = np.array(daymets)

  covariate_DataArray = xr.Dataset(data_vars={
                                           'ndvi':(('time', 'y', 'x'),ndvis),
                                           'srad': (('time', 'y', 'x'),daymets[:,0,:,:]),
                                           'vpd': (('time', 'y', 'x'),daymets[:,2,:,:]),
                                           'tsoil': (('time', 'y', 'x'),daymets[:,3,:,:]),
                                           'sm1': (('time', 'y', 'x'),daymets[:,4,:,:]),
                                           'sm2': (('time', 'y', 'x'),daymets[:,5,:,:]),
                                           'shortwave_radition': (('time', 'y', 'x'),daymets[:,6,:,:]),
                                           'tavg': (('time', 'y', 'x'),daymets[:,7,:,:]),
                                           'tmin': (('time', 'y', 'x'),daymets[:,8,:,:]),
                                           'prcp': (('time', 'y', 'x'),daymets[:,9,:,:]),
                                           'clay': (('time', 'y', 'x'),daymets[:,10,:,:])},
                                  coords={'time':dates,
                                          'y':starfm.y.values,
                                          'x':starfm.x.values 
                                          })
  
  if gap_fill==True:
    covariate_DataArray = daymet_qaqc(covariate_DataArray, bucket)
  
  covariate_DataArray.rio.write_crs(starfm.rio.crs.to_wkt(), inplace=True)
  covariate_DataArray.rio.write_transform(starfm.rio.transform(), inplace=True)
  covariate_DataArray.rio.nodata=0
  covariate_DataArray.to_netcdf(os.path.join(path_to_temp, 'temp_covariates_write.nc'))
  gs_write_blob(os.path.join(path_to_temp, 'temp_covariates_write.nc'), os.path.join(out_dir, out_name), bucket)
  
  return covariate_DataArray

def image_average_variables(ds, variable_list, bucket, path_to_temp, plot_dir=None):
  """Take spatial average of xarray Dataset for variables in variable_list. Returns pandas dataframe

  Args:
    ds (xarray.Dataset): dataset with variables to average and dimensions x, y, and time
    variable_list (list of strings): variable names to spatially average
    plot (string): full path to save timeseries plot to, will only plot if this argument is set

  Returns:
     pandas.DataFrame with a column for time, and spatial average of each variable in variable_list
  """
  
  if 'time' in ds.indexes:
    site_dates = pd.to_datetime(pd.DatetimeIndex(ds.indexes['time']))
    mean_indices=ds[variable_list].mean(dim=['y', 'x']).to_pandas()
    mean_indices = mean_indices.sort_values(by='time', ascending=True)
    
  elif 'doy_bins_lower' in ds.indexes:
    site_dates = ds.indexes['doy_bins_lower']
    mean_indices=ds[variable_list].mean(dim=['y', 'x']).to_pandas()
    mean_indices = mean_indices.sort_values(by='doy_bins_lower', ascending=True)
    
  else: return
  
  for variable in variable_list:
    if plot_dir!=None:
      fig, ax = plt.subplots()
      sns.lineplot(x=site_dates, y=mean_indices[variable].values, label=variable)
      plt.ylabel(variable)
      print('saving plot')
      plt.savefig(os.path.join(path_to_temp, 'temp_plot.jpg'), dpi=300)
      gs_write_blob(os.path.join(path_to_temp, 'temp_plot.jpg'), plot_dir+variable+'.jpg', bucket)
      plt.show()
      plt.clf()

  return mean_indices

def aggregate_for_spinup(ds, out_dir, out_name, date_min, date_max, bucket, path_to_temp, period = 5):
  
  """Average data in time slice temporally for RCTM spinup dataset

  Args:
    ds (xarray.Dataset): dataset with variables to average and dimensions x, y, and time
    out_dir (string): Google storage directory to write results to, does not include bucket name or 'gs://' prefix. e.g. Ranch_Runs/HLD/covariates_nc/
    out_name (string): name for file with extension (should be a .nc file)
    date_min (string): datetime to begin spinup aggregation, formatted YYYY-MM-DD
    date_max (string): datetime to end spinup aggregation, formatted YYYY-MM-DD
    bucket_name (string): google storage bucket name
    period (int): number of days per spinup time period, used to bin individual image dates

  Returns:
     xarray.Dataset and saves the spinup file to Google Cloud Storage
  """
  
  if not (isinstance(ds, xr.core.dataarray.DataArray) | isinstance(ds, xr.core.dataarray.Dataset)):
    ds=rxr.open_rasterio('gs://' + bucket.name + '/' + ds)
  
  ds=ds.sel(time=slice(date_min, date_max))
  
  orig_crs = ds.rio.crs.to_wkt()
  orig_transform = ds.rio.transform()
  
  #get day of year from datetime index
  doys = np.array(ds.indexes['time'].dayofyear)
  doy_bins_lower = (np.floor(doys/period).astype(np.int16)*period)+1
  
  ds = ds.assign_coords(doy_bins_lower=('time', doy_bins_lower))
  ds = ds.swap_dims({'time':'doy_bins_lower'})
  ds = ds.groupby('doy_bins_lower').mean()
  
  ds.rio.write_crs(orig_crs, inplace=True)
  ds.rio.write_transform(orig_transform, inplace=True)
  ds.rio.nodata=0
  
  ds.to_netcdf(os.path.join(path_to_temp, 'temp_covariates_spin_write.nc'))
  gs_write_blob(os.path.join(path_to_temp, 'temp_covariates_spin_write.nc'), os.path.join(out_dir, out_name), bucket)

  return ds

if __name__ == "__main__":

  parser=argparse.ArgumentParser()
  
  parser.add_argument("--starfm_in_dir", help="directory to output STARFM")
  parser.add_argument("--covariate_in_dir", help="directory to output STARFM")
  parser.add_argument("--out_dir", help="directory to output STARFM")
  parser.add_argument("--NLCD_in_dir", help="directory to output STARFM")
  parser.add_argument("--RAP_in_dir", help="directory to output STARFM")
  parser.add_argument("--param_file", help="directory to output STARFM")
  parser.add_argument("--bucket_name", help="bucket name")

  args=parser.parse_args()
  
  covariate_nc_outname = 'covariates.nc'
  covariate_csv_outname = 'covariates.csv'
  spin_nc_outname = 'covariates_spin.nc'
  spin_csv_outname = 'covariates_spin.csv'
  param_outname = args.out_dir + 'params.tif'
  
  
  #gen covariates
  ds = gen_covariates(args.starfm_in_dir, args.covariate_in_dir,  args.out_dir, covariate_nc_outname, args.bucket_name, '2002-01-01', '2023-01-01', gap_fill = True)
  #average indices
  df = image_average_variables(ds, ['ndvi','srad','vpd','tsoil','sm1','sm2','shortwave_radition','tavg','tmin','prcp','clay'], plot_dir=os.path.join(args.out_dir, 'transient_figs/'))
  df.to_csv('gs://' + args.bucket_name + '/' + args.out_dir + covariate_csv_outname)
  #gen spinup
  spin_ds = aggregate_for_spinup(ds, args.out_dir, spin_nc_outname, '2002-01-01', '2005-12-31', args.bucket_name,  period = 5)
  df = image_average_variables(spin_ds, ['ndvi','srad','vpd','tsoil','sm1','sm2','shortwave_radition','tavg','tmin','prcp','clay'], plot_dir=os.path.join(args.out_dir, 'spin_figs/'))
  df.to_csv('gs://' + args.bucket_name + '/' + args.out_dir + spin_csv_outname)
  #get spatial params
  get_spatial_RCTM_params(args.NLCD_in_dir, args.RAP_in_dir, args.param_file, param_outname, args.bucket_name, param_type='starfm')



  