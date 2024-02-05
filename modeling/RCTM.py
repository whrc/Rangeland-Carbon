import RCTM_utils as rctmu
import rioxarray as rxr
import xarray as xr
from google.cloud import storage
import sys
sys.path.insert(1, '../utils')
import utils
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import argparse
import os
import yaml
import time
#google storage init
bucket_name='rangelands'

storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
bucket = storage_client.get_bucket(bucket_name)
path_to_temp = '/home/amullen/temp/'

def InitCarbonModel(params, covariates, days=5):
  
  #initialize params for carbon modeling before a run
  #
  
  Tmult = np.exp((-(1/(covariates['tsoil']-params['beta2'])) + (1/params['beta1'])) * params['beta0'])
  #Tmult = np.exp(params['beta0']*((1/params['beta1'])-(1/(covariates['tsoil']-params['beta2']))))
  Wmult = rctmu.get_mult(params['MIN_SMsf'], params['MAX_SMsf'], covariates['sm1'])
  Smult = 0.7+0.3*(np.exp(-0.08*covariates['clay']))
  
  fsdethg = xr.where(covariates['sm1'] < params['MIN_SMsd'], params['h_SMsd'], params['l_SMsd'])
  
  fsdethg = xr.where((covariates['sm1'] > params['MIN_SMsd']) & (covariates['sm1'] < params['MAX_SMsd']), 
                      params['h_SMsd']-(covariates['sm1']-params['MIN_SMsd'])/(params['MAX_SMsd']-params['MIN_SMsd'])*(params['h_SMsd']-params['l_SMsd']), fsdethg)
    
  fsdeths = params['snsd']
  fallrt = params['frt']
  
  rdrg = xr.where(covariates['sm2'] < params['MIN_SMrd'], params['h_SMrd'], params['l_SMrd'])
  rdrg = xr.where((covariates['sm2'] > params['MIN_SMrd']) & (covariates['sm2'] < params['MAX_SMrd']), 
                      params['h_SMrd']-(covariates['sm2']-params['MIN_SMrd'])/(params['MAX_SMrd']-params['MIN_SMrd'])*(params['h_SMrd']-params['l_SMrd']), rdrg)
  
  rdrs = params['rdt']
  
  rdr = xr.where(covariates['tsoil'] > rdrs, rdrg, 0)
  
  
  
  #covariates['GPP'] = covariates['GPP'] * days #total gpp for period
  
  covariates['NPP'] = covariates['GPP'] * (1-params['F_Ra']) #subtract 'F_Ra' (fraction autotrophic respiration) to get NPP
  
  return Tmult, Wmult, Smult, fsdethg, fsdeths, fallrt, rdrg, rdrs, rdr, covariates
  
def get_time_index(indices):
  """Determines the name of the index used for the time dimension. This function is necesarry because covariate generation scripts use the index 'doy_bins_lower'
     for the time dimension in spinup files and 'time' for the time dimension in transient/contemporary covariate files

  Args:
    indices (xarray.DataArray.indexes or xarray.Dataset.indexes): xarray DataArray or Dataset indexes
  
  Returns:
     (string) time index name
  """ 
  
  if 'time' in indices: 
    return 'time'
  
  elif 'doy_bins_lower' in indices: 
    return 'doy_bins_lower'
  
  else: 
    print('uncertain time index')
    return ''
  
def CarbonModel(C_stocks, Tmult, Wmult, Smult, fsdethg, fsdeths, fallrt, rdrg, rdrs, rdr, covariates, params, days=5, spinup=False, point_mode=False):
  
  #run the carbon model
  C_stock_log = {
                'time': [],
                'AGB' : [],
                'BGB' : [],
                'AGL' : [],
                'BGL' : [],
                'POC' : [],
                'HOC' : []
  }

  #Tmult, Wmult, Smult, fsdethg, fsdeths, fallrt, rdrg, rdrs, rdr, covariates = InitCarbonModel(params, covariates)
  
  covariates['GPP'] = covariates['GPP'] * days
  covariates['NPP'] = covariates['NPP'] * days
  
  #determine time index based on input covariate file
  time_index = get_time_index(covariates.indexes)
  
  #create empty xarray data variable for heterotrophic respiration
  if point_mode:
    covariates['Rh'] = (covariates.indexes, np.zeros(np.shape(covariates['GPP'])))
  else:
    covariates['Rh'] = ([time_index, 'y', 'x'], np.zeros(np.shape(covariates['GPP'])))
    
  covariates['Rh'][:] = np.nan
  
  #iterate time steps
  for i in range(0, len(covariates.indexes[time_index])):
  
    if spinup==False:
      print(covariates.indexes[time_index][i])
    
    #get day-of-year based on time dimension
    if time_index == 'time':
      doy = pd.DatetimeIndex([pd.to_datetime(str(covariates.indexes['time'][i]))]).dayofyear.values[0]
      C_stock_log['time'].append(covariates.indexes['time'][i])
    
    if time_index == 'doy_bins_lower':
      doy = covariates.indexes[time_index][i]
      C_stock_log['time'].append(covariates.indexes[time_index][i])
    
    #if day-of-year before 230th, fsdethg equals function of soil moisture, else fsdeth equals 'snsd' param 
    if doy < 230:
      fsdeth = fsdethg[i]
    else: 
      fsdeth = fsdeths  
    
    #change in litter from previous timestep
    ddt_C_AGL = -1*(Tmult[i] * Wmult[i] * params['k_AGL'] * C_stocks['AGL'] * days) + C_stocks['AGB'] * fsdeth * fallrt #aboveground litter
    ddt_C_BGL = C_stocks['BGB'] * rdr[i] - (Tmult[i] * Wmult[i] * params['k_BGL'] * C_stocks['BGL'] * days) #belowground litter
    
    #update aboveground C stocks 
    C_stocks['AGB'] = C_stocks['AGB'] * (1 - fsdeth*fallrt) + covariates['NPP'][i] * (1.0 - params['F_Root']) #subtract root fraction from NPP
    C_stocks['BGB'] = C_stocks['BGB'] * (1 - rdr[i]) + covariates['NPP'][i] * params['F_Root'] #root fraction parameterization to get belowground biomass  
    
    #change in particulate and humic carbon from previous timestep                                    
    ddt_C_POC = (Tmult[i] * Wmult[i] * params['k_AGL'] * days * params['CSE_AGL'] * C_stocks['AGL'] + Tmult[i] * Wmult[i] * params['k_BGL'] * days * params['CSE_BGL'] * C_stocks['BGL']) - Tmult[i] * Wmult[i] * Smult[i] * params['k_POC'] * days * C_stocks['POC']   
    ddt_C_HOC = Tmult[i] * Wmult[i] * Smult[i] * params['k_POC'] * days * C_stocks['POC'] * params['CSE_POC'] - Tmult[i] * Wmult[i] * Smult[i] * params['k_HOC'] * days * C_stocks['HOC'] * (1 - params['CSE_HOC'])
    
    #update litter and soil carbon  
    C_stocks['AGL'] = (ddt_C_AGL) + C_stocks['AGL']
    C_stocks['BGL'] = (ddt_C_BGL) + C_stocks['BGL']
    C_stocks['POC'] = (ddt_C_POC) + C_stocks['POC']
    C_stocks['HOC'] = (ddt_C_HOC) + C_stocks['HOC']
    
    #total soil organic carbon  
    TSOC = C_stocks['POC'] + C_stocks['HOC'] + 1000
    
    #calculate heterotrophic respiration from litter and soil organic carbon
    R_AGL = C_stocks['AGL'] * (Tmult[i] * Wmult[i] * params['k_AGL'] * days * (1-params['CSE_AGL']))
    R_BGL = C_stocks['BGL'] * (Tmult[i] * Wmult[i] * params['k_BGL'] * days * (1-params['CSE_BGL']))
    R_POC = C_stocks['POC'] * (Tmult[i] * Wmult[i] * Smult[i] * params['k_POC'] * days * (1-params['CSE_POC']))
    R_HOC = C_stocks['HOC'] * (Tmult[i] * Wmult[i] * Smult[i] * params['k_HOC'] * days * (1-params['CSE_HOC']))
    
    #sum heterotrophic respiration components
    covariates['Rh'][i] = (-1) * (R_AGL + R_BGL + R_POC + R_HOC)
    
    C_stock_log['AGB'].append(C_stocks['AGB'])
    C_stock_log['BGB'].append(C_stocks['BGB'])
    C_stock_log['AGL'].append(C_stocks['AGL'])
    C_stock_log['BGL'].append(C_stocks['BGL'])
    C_stock_log['POC'].append(C_stocks['POC'])
    C_stock_log['HOC'].append(C_stocks['HOC'])
    
  #derive autotrophic respiration from NPP and GPP
  covariates['Ra'] = covariates['NPP'] - covariates['GPP']
  
  #Net ecosystem exchange
  covariates['NEE'] = covariates['GPP'] + (covariates['Rh'] + covariates['Ra'])
  
  #carbon flux components currently cumulative for timestep, so we will normalize by the timestep to get flux/day  
  covariates['NEE'] = covariates['NEE']/days
  covariates['GPP'] = covariates['GPP']/days
  covariates['Rh'] = covariates['Rh']/days
  covariates['Ra'] = covariates['Ra']/days
  covariates['NPP'] = covariates['NPP']/days
  
  return C_stocks, C_stock_log, covariates

def spinup_RCTM(C_stocks, covariates, params, spinup_years=0, spin_log_period=5, point_mode=False, save_checkpoint=False, save_path=''):
  
  covariates['ndvi_sr'] = rctmu.calc_sr(covariates['ndvi'])
  covariates['fpar'] = rctmu.calc_fPAR(covariates['ndvi'], params['ndvi_02'], params['ndvi_98'], covariates['ndvi_sr'], params['ndvi_sr_02'], params['ndvi_sr_98'], params['fPAR_min'], params['fPAR_max'])
  covariates['GPP'] = rctmu.calcGPP(params, covariates['tsoil'], covariates['sm2'], covariates['vpd'], covariates['shortwave_radition'], covariates['fpar'])
  
  #Change unit to match with C modeling
  covariates['sm1'] =  covariates['sm1'] * 100
  covariates['sm2'] =  covariates['sm2'] * 100
  covariates['tsoil'] = covariates['tsoil'] + 273.15
  
  AGB_stock_hist = []
  AGL_stock_hist = []
  BGB_stock_hist = []
  BGL_stock_hist = []
  POC_stock_hist = []
  HOC_stock_hist = []
  
  hist_years = []
  
  Tmult, Wmult, Smult, fsdethg, fsdeths, fallrt, rdrg, rdrs, rdr, covariates = InitCarbonModel(params, covariates)
  
  #log output to track spinup
  for spinup_year in range(0,spinup_years):
    
    if (spinup_year % spin_log_period == 0) & (spinup_year!=0):
        mean_AGB = np.nanmedian(C_stocks['AGB']).round(2)
        sd_AGB = np.nanstd(C_stocks['AGB']).round(2)
        AGB_stock_hist.append(mean_AGB)
        
        mean_AGL = np.nanmedian(C_stocks['AGL']).round(2)
        sd_AGL = np.nanstd(C_stocks['AGL']).round(2)
        AGL_stock_hist.append(mean_AGL)
        
        mean_BGB = np.nanmedian(C_stocks['BGB']).round(2)
        sd_BGB = np.nanstd(C_stocks['BGB']).round(2)
        BGB_stock_hist.append(mean_BGB)
        
        mean_BGL = np.nanmedian(C_stocks['BGL']).round(2)
        sd_BGL = np.nanstd(C_stocks['BGL']).round(2)
        BGL_stock_hist.append(mean_BGL)
        
        mean_POC = np.nanmedian(C_stocks['POC']).round(2)
        sd_POC = np.nanstd(C_stocks['POC']).round(2)
        POC_stock_hist.append(mean_POC)
        
        mean_HOC = np.nanmedian(C_stocks['HOC']).round(2)
        sd_HOC = np.nanstd(C_stocks['HOC']).round(2)
        HOC_stock_hist.append(mean_HOC)
        
        hist_years.append(spinup_year)
        print(f'year {spinup_year}/{spinup_years}:, AGB: {mean_AGB} ({sd_AGB}), BGB: {mean_BGB} ({sd_BGB}), AGL: {mean_AGL} ({sd_AGL}), BGL: {mean_BGL} ({sd_BGL}), POC: {mean_POC} ({sd_POC}), HOC: {mean_HOC} ({sd_HOC})')
        if save_checkpoint:
          spinup_C_stocks = xr.Dataset(data_vars={
                                               'AGB':(('y', 'x'),C_stocks['AGB'].values),
                                               'BGB': (('y', 'x'),C_stocks['BGB'].values),
                                               'AGL': (('y', 'x'),C_stocks['AGL'].values),
                                               'BGL': (('y', 'x'),C_stocks['BGL'].values),
                                               'POC': (('y', 'x'),C_stocks['POC'].values),
                                               'HOC': (('y', 'x'),C_stocks['HOC'].values)
                                              },
                                      coords={
                                              'y':covariates.y.values,
                                              'x':covariates.x.values 
                                              })
      
          spinup_C_stocks.rio.write_crs(covariates.rio.crs.to_wkt(), inplace=True)
          spinup_C_stocks.rio.write_transform(covariates.rio.transform(), inplace=True)
          spinup_C_stocks.rio.nodata=0
          
          ts = '_'.join(str(time.time()).split('.'))
          spinup_C_stocks.rio.to_raster(os.path.join(path_to_temp, f'temp_C_stocks_write_{ts}.tif'))
          utils.gs_write_blob(os.path.join(path_to_temp, f'temp_C_stocks_write_{ts}.tif'), save_path, bucket)
          os.remove(os.path.join(path_to_temp, f'temp_C_stocks_write_{ts}.tif'))
          
    C_stocks, C_stock_log, covariates = CarbonModel(C_stocks, Tmult, Wmult, Smult, fsdethg, fsdeths, fallrt, rdrg, rdrs, rdr, covariates, params, days=5, spinup=True, point_mode=point_mode)
    
        
  C_stock_hist = {'year': hist_years, 'AGB': AGB_stock_hist, 'AGL': AGL_stock_hist, 'BGB': BGB_stock_hist, 'BGL': BGL_stock_hist, 'POC': POC_stock_hist, 'HOC': HOC_stock_hist}
        
  return covariates, C_stocks, C_stock_hist

def RCTM(C_stocks, covariates, params, spinup_years=0, spin_log_period=1):
  print('Running RCTM')
  print('calculating GPP')
  covariates['ndvi_sr'] = rctmu.calc_sr(covariates['ndvi'])
  covariates['fpar'] = rctmu.calc_fPAR(covariates['ndvi'], params['ndvi_02'], params['ndvi_98'], covariates['ndvi_sr'], params['ndvi_sr_02'], params['ndvi_sr_98'], params['fPAR_min'], params['fPAR_max'])
  covariates['GPP'] = rctmu.calcGPP(params, covariates['tsoil'], covariates['sm2'], covariates['vpd'], covariates['shortwave_radition'], covariates['fpar'])
  
  
  
  #Change unit to match with C modeling
  covariates['sm1'] =  covariates['sm1'] * 100
  covariates['sm2'] =  covariates['sm2'] * 100
  covariates['tsoil'] = covariates['tsoil'] + 273.15
  
  #init carbon model
  Tmult, Wmult, Smult, fsdethg, fsdeths, fallrt, rdrg, rdrs, rdr, covariates = InitCarbonModel(params, covariates)
  
  print('NEE model')
  C_stocks, C_stock_log, covariates = CarbonModel(C_stocks, Tmult, Wmult, Smult, fsdethg, fsdeths, fallrt, rdrg, rdrs, rdr, covariates, params, days=5, spinup=False)   
  
  return covariates, C_stocks, C_stock_log

def mask_variables(dataset):
    # Get the list of variable names in the dataset
    variable_names = list(dataset.data_vars.keys())

    # Create a mask for each variable based on NaN values in any other variable
    for var_name in variable_names:
        # Create a mask for NaN values in the entire dataset (excluding the current variable)
        nan_mask = np.isnan(dataset.drop_vars(var_name).to_array().values).any(axis=0)
        # Apply the mask to the current variable
        dataset[var_name] = dataset[var_name].where(~nan_mask)

    return dataset

def main(force_pft, point_mode, spin_years, run_transient, init_C_stocks_with_image, ranch,
       path_to_RCTM_params, path_to_RCTM_spatial_params, path_to_spin_covariates_point, 
       path_to_spin_covariates_spatial, C_stock_spin_out_path_point, C_stock_spin_out_path,
       C_stock_inits_yaml, C_stock_init_image, spatial_spin_fig_path, transient_covariate_path,
       transient_C_stock_hist, transient_flux_hist):
  print(ranch)
  #force_pft='grass-tree'
  #point_mode=False
  #spin_years=0
  #init_C_stocks_with_image=True
  #ranch='VDZ'
  
  
  #C_stock_inits_yaml = '/home/amullen/Rangeland-Carbon/modeling/Ranch_C_Stock_Inits.yaml'
  #path_to_RCTM_params = '/home/amullen/Rangeland-Carbon/modeling/RCTM_params.yaml'
  #path_to_spin_covariates_point = f'gs://rangelands/Ranch_Runs/{ranch}/covariates_nc/covariates_spin.csv'
  #path_to_RCTM_spatial_params = f'gs://rangelands/Ranch_Runs/{ranch}/covariates_nc/params.tif'
  #path_to_spin_covariates_spatial = f'gs://rangelands/Ranch_Runs/{ranch}/covariates_nc/covariates_spin.nc'
  #C_stock_spin_out_path_point = f'Ranch_Runs/{ranch}/covariates_nc/C_stocks_{force_pft}.yaml'
  #C_stock_spin_out_path = f'Ranch_Runs/{ranch}/covariates_nc/C_stocks_{force_pft}.tif'
  #C_stock_init_image = 'gs://rangelands/' + C_stock_spin_out_path
  #spatial_spin_fig_path = f'output/RCTM/{ranch}_spinup_stocks.jpg'
  #transient_covariate_path = f'gs://rangelands/Ranch_Runs/{ranch}/covariates_nc/covariates.nc'
  #transient_C_stock_hist = f'Ranch_Runs/{ranch}/results/C_stocks_hist_{force_pft}.nc'
  #transient_flux_hist = f'Ranch_Runs/{ranch}/results/flux_hist_{force_pft}.nc'
  
  path_to_C_stock_inits = C_stock_inits_yaml

  ######################################################################################
  #if args.covariate_in_dir.endswith('.csv'):
  #  point_mode=True
  #else:
  #  point_mode=False
  blob=bucket.blob(C_stock_inits_yaml)
  with blob.open("r") as file:
      init_C_stocks = yaml.safe_load(file)
              
  if point_mode:
    print('running model in point mode')
    params=[]
    with open(path_to_RCTM_params, 'r') as file:
      params = yaml.safe_load(file)
    
    #get keys (param names), doesn't matter which pft since they all have the same parameters
    params = params['RCTM_params'][force_pft]['starfm']
    C_stocks = {'AGB' : init_C_stocks['AGB'],
                'BGB' : init_C_stocks['BGB'],
                'AGL' : init_C_stocks['AGL'],
                'BGL' : init_C_stocks['BGL'],
                'POC' : init_C_stocks['POC'],
                'HOC' : init_C_stocks['HOC']}
    
    #read spinup covariates  
    spin_covariates = pd.read_csv(path_to_spin_covariates_point)
    spin_covariates = spin_covariates.set_index('doy_bins_lower')
    spin_covariates = spin_covariates.to_xarray()
    
  else:
    print('running model in spatial mode')
    #Define parameters - these are calibrated parameters for RCTM v1.0
    params = rxr.open_rasterio(path_to_RCTM_spatial_params, masked=True)
    params = rctmu.xr_dataset_to_data_array(params)
    
    if force_pft!='':
      print(f'forcing params to {force_pft}:')
      #open point param yaml file
      params_point=[]
      
      with open(path_to_RCTM_params, 'r') as file:
        params_point = yaml.safe_load(file)
        
      params_point = params_point['RCTM_params'][force_pft]['starfm']
      
      for p in params_point.keys():
        print(f'\t{p}: {params_point[p]}')
        params[p].values = np.where(~np.isnan(params[p].values), params_point[p], np.nan)
    
    #read spinup covariates
    spin_covariates = rxr.open_rasterio(path_to_spin_covariates_spatial, masked=True)
    
    spin_covariates = mask_variables(spin_covariates) #this line consumes a lot of memory!!!
    
    #Define initial C pools
    if init_C_stocks_with_image:
      print('initializing C stocks from image')
      C_stocks = rxr.open_rasterio(C_stock_init_image, masked=True)
      C_stocks = rctmu.xr_dataset_to_data_array(C_stocks)
      
    else:
      print('initializing C stocks from yaml file')
      C_stocks = {'AGB' : np.full(np.shape(spin_covariates['ndvi'].values[0]), init_C_stocks['AGB']),
                  'BGB' : np.full(np.shape(spin_covariates['ndvi'].values[0]), init_C_stocks['BGB']),
                  'AGL' : np.full(np.shape(spin_covariates['ndvi'].values[0]), init_C_stocks['AGL']),
                  'BGL' : np.full(np.shape(spin_covariates['ndvi'].values[0]), init_C_stocks['BGL']),
                  'POC' : np.full(np.shape(spin_covariates['ndvi'].values[0]), init_C_stocks['POC']),
                  'HOC' : np.full(np.shape(spin_covariates['ndvi'].values[0]), init_C_stocks['HOC'])
                  }
                  
    if params.rio.transform()!=spin_covariates.rio.transform():
      params = params.rio.reproject_match(spin_covariates)
        
  if spin_years > 0: 
    
    if point_mode:
      #run spinup in point mode
      spin_covariates, C_stocks, C_stock_hist = spinup_RCTM(C_stocks, spin_covariates, params, spinup_years=spin_years, point_mode=point_mode)
      
      #values is C_stocks dict are xarray DataArrays, so we need to extrack the numerical  values
      for key in C_stocks.keys():
        C_stocks[key] = float(C_stocks[key].values.mean())
      print(C_stocks)
      
      #write final C stocks to yaml
      ts = '_'.join(str(time.time()).split('.'))
      file=open(os.path.join(path_to_temp, f'temp_C_stocks_write_{force_pft}_{ts}.yaml'),'w')
      yaml.dump(C_stocks,file)
      file.close()
      utils.gs_write_blob(os.path.join(path_to_temp, f'temp_C_stocks_write_{force_pft}_{ts}.yaml'), C_stock_spin_out_path_point, bucket)
      os.remove(os.path.join(path_to_temp, f'temp_C_stocks_write_{force_pft}_{ts}.yaml'))
      
    else:
      spin_covariates, C_stocks, C_stock_hist = spinup_RCTM(C_stocks, spin_covariates, params, spinup_years=spin_years, point_mode=point_mode, save_checkpoint = True, save_path = C_stock_spin_out_path)
      
      #Save C stocks
      spinup_C_stocks = xr.Dataset(data_vars={
                                               'AGB':(('y', 'x'),C_stocks['AGB'].values),
                                               'BGB': (('y', 'x'),C_stocks['BGB'].values),
                                               'AGL': (('y', 'x'),C_stocks['AGL'].values),
                                               'BGL': (('y', 'x'),C_stocks['BGL'].values),
                                               'POC': (('y', 'x'),C_stocks['POC'].values),
                                               'HOC': (('y', 'x'),C_stocks['HOC'].values)
                                              },
                                      coords={
                                              'y':spin_covariates.y.values,
                                              'x':spin_covariates.x.values 
                                              })
      
      spinup_C_stocks.rio.write_crs(spin_covariates.rio.crs.to_wkt(), inplace=True)
      spinup_C_stocks.rio.write_transform(spin_covariates.rio.transform(), inplace=True)
      spinup_C_stocks.rio.nodata=0
      
      ts = '_'.join(str(time.time()).split('.'))
      spinup_C_stocks.rio.to_raster(os.path.join(path_to_temp, f'temp_C_stocks_write_{ts}.tif'))
      utils.gs_write_blob(os.path.join(path_to_temp, f'temp_C_stocks_write_{ts}.tif'), C_stock_spin_out_path, bucket)
      os.remove(os.path.join(path_to_temp, f'temp_C_stocks_write_{ts}.tif'))
      
      #Plot averaged time series
      #fig, axes = plt.subplots(6,1, figsize=(5, 18))
      
      #sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['AGB'], ax=axes[0], label = 'AGB')
      #sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['AGL'], ax=axes[1], label = 'AGL')
      #sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['BGB'], ax=axes[2], label = 'BGB')
      #sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['BGL'], ax=axes[3], label = 'BGL')
      #sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['POC'], ax=axes[4], label = 'POC')
      #sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['HOC'], ax=axes[5], label = 'HOC')
      
      #plt.savefig(spatial_spin_fig_path, dpi=300)
      
      #utils.image_average_variables(spin_covariates, ['GPP'], time_index = 'doy_bins_lower', plot_dir='output/RCTM')
  
  if point_mode or not run_transient:
    pass
  
  else:
    print('reading in covariates') 
    covariates = rxr.open_rasterio(transient_covariate_path, masked=True)
    
    crs = covariates.rio.crs.to_wkt()
    transform = covariates.rio.transform()
    #C_stocks = rxr.open_rasterio(C_stock_spin_out_path, masked=True)
    #C_stocks = rctmu.xr_dataset_to_data_array(C_stocks)
    
    covariates, C_stocks, C_stock_log = RCTM(C_stocks, covariates, params)
    
    cstock_DataArray = xr.Dataset(data_vars={
                                             'AGB':(('time', 'y', 'x'), C_stock_log['AGB']),
                                             'BGB': (('time', 'y', 'x'), C_stock_log['BGB']),
                                             'AGL': (('time', 'y', 'x'), C_stock_log['AGL']),
                                             'BGL': (('time', 'y', 'x'), C_stock_log['BGL']),
                                             'POC': (('time', 'y', 'x'), C_stock_log['POC']),
                                             'HOC': (('time', 'y', 'x'), C_stock_log['HOC'])},
                                    coords={'time':covariates.indexes['time'].values,
                                            'y':covariates.y.values,
                                            'x':covariates.x.values 
                                            })
    
    cstock_DataArray.rio.write_crs(crs, inplace=True)
    cstock_DataArray.rio.write_transform(transform, inplace=True)
    cstock_DataArray.rio.nodata=0
    
    ts = '_'.join(str(time.time()).split('.'))
    cstock_DataArray.to_netcdf(os.path.join(path_to_temp, f'temp_C_stocks_write_{ts}.nc'))
    utils.gs_write_blob(os.path.join(path_to_temp, f'temp_C_stocks_write_{ts}.nc'), transient_C_stock_hist, bucket)
    os.remove(os.path.join(path_to_temp, f'temp_C_stocks_write_{ts}.nc'))
    cstock_DataArray = None
    C_stocks = None
    C_stock_log = None
    
    covariates.rio.write_crs(crs, inplace=True)
    covariates.rio.write_transform(transform, inplace=True)
    covariates.rio.nodata=0
    
    ts = '_'.join(str(time.time()).split('.'))
    covariates[['GPP', 'Rh', 'Ra', 'NEE', 'NPP']].to_netcdf(os.path.join(path_to_temp, f'temp_fluxes_write_{ts}.nc'))
    utils.gs_write_blob(os.path.join(path_to_temp, f'temp_fluxes_write_{ts}.nc'), transient_flux_hist, bucket)
    os.remove(os.path.join(path_to_temp, f'temp_fluxes_write_{ts}.nc'))
    
      
    #utils.image_average_variables(covariates, ['GPP', 'Rh', 'Ra', 'NEE', 'NPP'], plot_dir='output/RCTM/transient')

if __name__ == '__main__':
  parser=argparse.ArgumentParser()
  
  parser.add_argument("--force_pft", help="directory to output STARFM", const='')
  parser.add_argument("--point_mode", help="directory to output STARFM", const=False)
  parser.add_argument("--spin_years", help="directory to output STARFM", const=0)
  parser.add_argument("--run_transient", help="directory to output STARFM", const=True)
  parser.add_argument("--init_C_stocks_with_image", help="directory to output STARFM", const=False)
  parser.add_argument("--ranch", help="directory to output STARFM", const='')
  
  parser.add_argument("--path_to_RCTM_params", help="directory to output STARFM", const='')
  parser.add_argument("--path_to_RCTM_spatial_params", help="directory to output STARFM", const='')
  parser.add_argument("--path_to_spin_covariates_point", help="directory to output STARFM", const='')
  
  parser.add_argument("--path_to_spin_covariates_spatial", help="directory to output STARFM", const='')
  
  parser.add_argument("--C_stock_spin_out_path_point", help="directory to output STARFM", const='')
  parser.add_argument("--C_stock_spin_out_path", help="directory to output STARFM", const='')
  
  parser.add_argument("--C_stock_inits_yaml", help="directory to output STARFM", const='')
  parser.add_argument("--C_stock_init_image", help="directory to output STARFM", const='')
  
  parser.add_argument("--spatial_spin_fig_path", help="directory to output STARFM", const='')
  
  parser.add_argument("--transient_covariate_path", help="bucket name", const='')
  parser.add_argument("--transient_C_stock_hist", help="bucket name", const='')
  parser.add_argument("--transient_flux_hist", help="bucket name", const='')

  args=parser.parse_args()
  
  main(args.force_pft, args.point_mode, args.spin_years, args.args.run_transient, args.init_C_stocks_with_image, args.ranch,
       args.path_to_RCTM_params, args.path_to_RCTM_spatial_params, args.path_to_spin_covariates_point, 
       args.path_to_spin_covariates_spatial, args.C_stock_spin_out_path_point, args.C_stock_spin_out_path,
       args.C_stock_inits_yaml, args.C_stock_init_image, args.spatial_spin_fig_path, args.transient_covariate_path,
       args.transient_C_stock_hist, args.transient_flux_hist)