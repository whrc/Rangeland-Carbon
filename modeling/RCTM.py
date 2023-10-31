import RCTM_utils as rctmu
import rioxarray as rxr
import xarray as xr
from google.cloud import storage
import sys
sys.path.insert(1, '../utils')
import utils
import numpy as np
import pandas as pd

#google storage init
bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')

def InitCarbonModel(params, covariates, days=5):
  #initialize params for carbon modeling before a run
  #
  
  Tmult= np.exp((-(1/(covariates['tsoil']-params['beta2'])) + (1/params['beta1'])) * params['beta0'])
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
                      params['h_SMrd']-(covariates['sm2']-params['MIN_SMrd'])/(params['MAX_SMrd']-params['MIN_SMrd'])*(params['h_SMrd']-params['l_SMrd']), fsdethg)
  
  rdrs = params['rdt']
  
  rdr = xr.where(covariates['tsoil'] > rdrs, rdrg, 0)
      
  covariates['GPP'] = covariates['GPP'] * days #total gpp for period
      
  covariates['NPP'] = covariates['GPP'] * (1-params['F_Ra']) #subtract 'F_Ra' (fraction autotrophic respiration) to get NPP
  
  return Tmult, Wmult, Smult, fsdethg, fsdeths, fallrt, rdrg, rdrs, rdr, covariates
  

def CarbonModel(C_stocks, covariates, params, days=5, spinup=False):
  
  #run the carbon model
  
  Tmult, Wmult, Smult, fsdethg, fsdeths, fallrt, rdrg, rdrs, rdr, covariates = InitCarbonModel(params, covariates)
  
  covariates['Rh'] = (['time', 'y', 'x'], np.zeros(np.shape(covariates['GPP'])))
  covariates['Rh'][:] = np.nan
  
  for i in range(0, len(covariates.indexes['time'])):
    print(covariates.indexes['time'][i])
    doy = pd.DatetimeIndex([pd.to_datetime(str(covariates.indexes['time'][i]))]).dayofyear.values[0]
    if doy < 230:
      fsdeth = fsdethg[i]
    else: 
      fsdeth = fsdeths  #if before 46th period, fsdeth equals 'snsd' param, else fsdethg equals function of soil moisture
      
    if spinup and i==0:
      C_stocks['AGB'] = C_stocks['AGB'] * (1 - fsdeth*fallrt) + covariates['NPP'][i] * (1.0 - params['F_Root']) #subtract root fraction from NPP
      C_stocks['BGB'] = C_stocks['BGB'] * (1 - rdr[i]) + covariates['NPP'][i] * params['F_Root'] #root fraction parameterization to get belowground biomass
      
      ddt_C_AGL = C_stocks['AGB'] * fsdeth * fallrt - (Tmult[i] * Wmult[i] * params['k_AGL'] * days * C_stocks['AGL']) #aboveground litter
      ddt_C_BGL = C_stocks['BGB'] * rdr[i] - (Tmult[i] * Wmult[i] * params['k_BGL'] * days * C_stocks['BGL']) #belowground litter
        
    else:
      ddt_C_AGL = -1*(Tmult[i] * Wmult[i] * params['k_AGL'] * C_stocks['AGL'] * days) + C_stocks['AGB'] * fsdeth * fallrt #aboveground litter
      ddt_C_BGL = C_stocks['BGB'] * rdr[i] - (Tmult[i] * Wmult[i] * params['k_BGL'] * C_stocks['BGL'] * days) #belowground litter
      
      C_stocks['AGB'] = C_stocks['AGB'] * (1 - fsdeth*fallrt) + covariates['NPP'][i] * (1.0 - params['F_Root']) #subtract root fraction from NPP
      C_stocks['BGB'] = C_stocks['BGB'] * (1 - rdr[i]) + covariates['NPP'][i] * params['F_Root'] #root fraction parameterization to get belowground biomass
      
                                          
    ddt_C_POC = (Tmult[i] * Wmult[i] * params['k_AGL'] * days * params['CSE_AGL'] * C_stocks['AGL'] + Tmult[i] * Wmult[i] * params['k_BGL'] * days * params['CSE_BGL'] * C_stocks['BGL']) - Tmult[i] * Wmult[i] * Smult[i] * params['k_POC'] * days * C_stocks['POC']
      
    ddt_C_HOC = Tmult[i] * Wmult[i] * Smult[i] * params['k_POC'] * days * C_stocks['POC'] * params['CSE_POC'] - Tmult[i] * Wmult[i] * Smult[i] * params['k_HOC'] * days * C_stocks['HOC'] * (1 - params['CSE_HOC'])
      
    C_stocks['AGL'] = (ddt_C_AGL) + C_stocks['AGL']
    C_stocks['BGL'] = (ddt_C_BGL) + C_stocks['BGL']
    C_stocks['POC'] = (ddt_C_POC) + C_stocks['POC']
    C_stocks['HOC'] = (ddt_C_HOC) + C_stocks['HOC']
      
    TSOC = C_stocks['POC'] + C_stocks['HOC'] + 1000
    
    R_AGL = C_stocks['AGL'] * (Tmult[i] * Wmult[i] * params['k_AGL'] * days * (1-params['CSE_AGL']))
    R_BGL = C_stocks['BGL'] * (Tmult[i] * Wmult[i] * params['k_BGL'] * days * (1-params['CSE_BGL']))
    R_POC = C_stocks['POC'] * (Tmult[i] * Wmult[i] * Smult[i] * params['k_POC'] * days * (1-params['CSE_POC']))
    R_HOC = C_stocks['HOC'] * (Tmult[i] * Wmult[i] * Smult[i] * params['k_HOC'] * days * (1-params['CSE_HOC']))
      
    covariates['Rh'][i].values = (-1) * (R_AGL + R_BGL + R_POC + R_HOC) 
  
  covariates['Ra'] = covariates['NPP'] - covariates['GPP'] 
  covariates['NEE'] = covariates['GPP'] + (covariates['Rh'] + covariates['Ra'])
    
  covariates['NEE'] = covariates['NEE']/days
  covariates['GPP'] = covariates['GPP']/days
  covariates['Rh'] = covariates['Rh']/days
  covariates['Ra'] = covariates['Ra']/days
  covariates['NPP'] = covariates['NPP']/days
    
  return C_stocks, covariates

def RCTM(C_stocks, covariates, params, spinup_years=0, spin_log_period=1):
  
  covariates['ndvi_sr'] = rctmu.calc_sr(covariates['ndvi'])
  covariates['fpar'] = rctmu.calc_fPAR(covariates['ndvi'], params['ndvi_02'], params['ndvi_98'], covariates['ndvi_sr'], params['ndvi_sr_02'], params['ndvi_sr_98'], params['fPAR_min'], params['fPAR_max'])
  covariates['GPP'] = rctmu.calcGPP(params, covariates['tsoil'], covariates['sm2'], covariates['vpd'], covariates['shortwave_radition'], covariates['fpar'])
  
  #Change unit to match with C modeling
  covariates['sm1'] =  covariates['sm1'] * 100
  covariates['sm2'] =  covariates['sm2'] * 100
  covariates['tsoil'] = covariates['tsoil'] + 273.15
  
  C_stock_hist = []
  
  for spinup_year in range(0,spinup_years):
    C_stocks, covariates = CarbonModel(C_stocks, covariates, params, days=5, spinup=False)
    if spinup_years > 0:
      if spinup_year % spin_log_period == 0:
        C_stock_hist.append(np.nanmean(C_stocks['AGB']))
        
  return covariates, C_stocks, C_stock_hist

#Define parameters - these are calibrated parameters for RCTM v1.0
param_path = 'gs://rangelands/Ameriflux_sites/A32_starfm/landcover/params_NLCD_2019.tif'
params = rxr.open_rasterio(param_path, masked=True)
params = rctmu.xr_dataset_to_data_array(params)

spin_cov_path = 'gs://rangelands/Ameriflux_sites/A32_starfm/covariates_v2/A32_covariates_site_test.nc'
covariates = rxr.open_rasterio(spin_cov_path, masked=True)

#Define initial C pools
C_stocks = {'AGB' : np.full(np.shape(covariates['ndvi'].values[0]), 200),
            'BGB' : np.full(np.shape(covariates['ndvi'].values[0]), 700),
            'AGL' : np.full(np.shape(covariates['ndvi'].values[0]), 200),
            'BGL' : np.full(np.shape(covariates['ndvi'].values[0]), 700),
            'POC' : np.full(np.shape(covariates['ndvi'].values[0]), 1000),
            'HOC' : np.full(np.shape(covariates['ndvi'].values[0]), 5000)
            } 

covariates, C_stocks, C_stock_hist = RCTM(C_stocks, covariates, params, spinup_years=5)

print(C_stock_hist)

utils.image_average_variables(covariates, ['GPP'], plot_dir='output/RCTM')
