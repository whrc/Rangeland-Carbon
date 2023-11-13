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
  
  if 'time' in covariates.indexes: 
    time_index = 'time'
  elif 'doy_bins_lower' in covariates.indexes: 
    time_index = 'doy_bins_lower'
  else: 
    print('uncertain time index')
    return
    
  covariates['Rh'] = ([time_index, 'y', 'x'], np.zeros(np.shape(covariates['GPP'])))
  covariates['Rh'][:] = np.nan
  
  for i in range(0, len(covariates.indexes[time_index])):
    if spinup==False:
      print(covariates.indexes[time_index][i])
    
    if time_index == 'time':
      doy = pd.DatetimeIndex([pd.to_datetime(str(covariates.indexes['time'][i]))]).dayofyear.values[0]
    
    if time_index == 'doy_bins_lower':
      doy = covariates.indexes[time_index][i]
      
    if doy < 230:
      fsdeth = fsdethg[i]
    else: 
      fsdeth = fsdeths  #if before 46th period, fsdeth equals 'snsd' param, else fsdethg equals function of soil moisture
      
    #if spinup and i==0:
    
    #  print(np.nanmean(C_stocks['AGB']))
    #  C_stocks['AGB'] = C_stocks['AGB'] * (1 - fsdeth*fallrt) + covariates['NPP'][i] * (1.0 - params['F_Root']) #subtract root fraction from NPP
    #  C_stocks['BGB'] = C_stocks['BGB'] * (1 - rdr[i]) + covariates['NPP'][i] * params['F_Root'] #root fraction parameterization to get belowground biomass
    #  print(np.nanmean(C_stocks['AGB']))
      
    #  ddt_C_AGL = C_stocks['AGB'] * fsdeth * fallrt - (Tmult[i] * Wmult[i] * params['k_AGL'] * days * C_stocks['AGL']) #aboveground litter
    #  ddt_C_BGL = C_stocks['BGB'] * rdr[i] - (Tmult[i] * Wmult[i] * params['k_BGL'] * days * C_stocks['BGL']) #belowground litter
        
    
    ddt_C_AGL = -1*(Tmult[i] * Wmult[i] * params['k_AGL'] * C_stocks['AGL'] * days) + C_stocks['AGB'] * fsdeth * fallrt #aboveground litter
    ddt_C_BGL = C_stocks['BGB'] * rdr[i] - (Tmult[i] * Wmult[i] * params['k_BGL'] * C_stocks['BGL'] * days) #belowground litter
    #print(np.nanmean(C_stocks['BGB']))  
    C_stocks['AGB'] = C_stocks['AGB'] * (1 - fsdeth*fallrt) + covariates['NPP'][i] * (1.0 - params['F_Root']) #subtract root fraction from NPP
    C_stocks['BGB'] = C_stocks['BGB'] * (1 - rdr[i]) + covariates['NPP'][i] * params['F_Root'] #root fraction parameterization to get belowground biomass
    #print(np.nanmean(C_stocks['BGB']))  
                                          
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

    covariates['Rh'][i] = (-1) * (R_AGL + R_BGL + R_POC + R_HOC)
  
  covariates['Ra'] = covariates['NPP'] - covariates['GPP'] 
  covariates['NEE'] = covariates['GPP'] + (covariates['Rh'] + covariates['Ra'])
    
  covariates['NEE'] = covariates['NEE']/days
  covariates['GPP'] = covariates['GPP']/days
  covariates['Rh'] = covariates['Rh']/days
  covariates['Ra'] = covariates['Ra']/days
  covariates['NPP'] = covariates['NPP']/days
    
  return C_stocks, covariates

def spinup_RCTM(C_stocks, covariates, params, spinup_years=0, spin_log_period=5):
  
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
  
  for spinup_year in range(0,spinup_years):
    
    if (spinup_year % spin_log_period == 0) & (spinup_year!=0):
        mean_AGB = np.nanmean(C_stocks['AGB'])
        AGB_stock_hist.append(mean_AGB)
        
        mean_AGL = np.nanmean(C_stocks['AGL'])
        AGL_stock_hist.append(mean_AGL)
        
        mean_BGB = np.nanmean(C_stocks['BGB'])
        BGB_stock_hist.append(mean_BGB)
        
        mean_BGL = np.nanmean(C_stocks['BGL'])
        BGL_stock_hist.append(mean_BGL)
        
        mean_POC = np.nanmean(C_stocks['POC'])
        POC_stock_hist.append(mean_POC)
        
        mean_HOC = np.nanmean(C_stocks['HOC'])
        HOC_stock_hist.append(mean_HOC)
        
        hist_years.append(spinup_year)
        print(f'year {spinup_year}:, AGB: {mean_AGB}, AGL: {mean_AGL}, BGB: {mean_BGB}, BGL: {mean_BGL}, POC: {mean_POC}, HOC: {mean_HOC}')
        
    C_stocks, covariates = CarbonModel(C_stocks, covariates, params, days=5, spinup=True)
        
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
  
  print('NEE model')
  C_stocks, covariates = CarbonModel(C_stocks, covariates, params, days=5, spinup=False)   
  
  return covariates, C_stocks

if __name__ == "__main_off__":

  parser=argparse.ArgumentParser()
  
  parser.add_argument("--starfm_in_dir", help="directory to output STARFM")
  parser.add_argument("--covariate_in_dir", help="directory to output STARFM")
  parser.add_argument("--out_dir", help="directory to output STARFM")
  parser.add_argument("--covariate_nc_outname", help="directory to output STARFM")
  parser.add_argument("--covariate_csv_outname", help="directory to output STARFM")
  parser.add_argument("--spin_nc_outname", help="directory to output STARFM")
  parser.add_argument("--landcover_in_dir", help="directory to output STARFM")
  parser.add_argument("--param_file", help="directory to output STARFM")
  parser.add_argument("--param_outname", help="directory to output STARFM")
  parser.add_argument("--bucket_name", help="bucket name")

  args=parser.parse_args()
  
  #Define parameters - these are calibrated parameters for RCTM v1.0
  param_path = 'gs://rangelands/Ranch_Runs/HB/covariates_nc/HB_params.tif'
  params = rxr.open_rasterio(param_path, masked=True)
  params = rctmu.xr_dataset_to_data_array(params)
  
  spin_cov_path = 'gs://rangelands/Ranch_Runs/HB/covariates_nc/HB_covariates_spin.nc'
  spin_covariates = rxr.open_rasterio(spin_cov_path, masked=True)
  
  
  #Define initial C pools
  C_stocks = {'AGB' : np.full(np.shape(spin_covariates['ndvi'].values[0]), 200),
              'BGB' : np.full(np.shape(spin_covariates['ndvi'].values[0]), 700),
              'AGL' : np.full(np.shape(spin_covariates['ndvi'].values[0]), 200),
              'BGL' : np.full(np.shape(spin_covariates['ndvi'].values[0]), 700),
              'POC' : np.full(np.shape(spin_covariates['ndvi'].values[0]), 1000),
              'HOC' : np.full(np.shape(spin_covariates['ndvi'].values[0]), 5000)
              } 
            

  spin_covariates, C_stocks, C_stock_hist = spinup_RCTM(C_stocks, spin_covariates, params, spinup_years=100)
  
  fig, axes = plt.subplots(6,1, figsize=(5, 18))
  
  sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['AGB'], ax=axes[0], label = 'AGB')
  sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['AGL'], ax=axes[1], label = 'AGL')
  sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['BGB'], ax=axes[2], label = 'BGB')
  sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['BGL'], ax=axes[3], label = 'BGL')
  sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['POC'], ax=axes[4], label = 'POC')
  sns.lineplot(x=C_stock_hist['year'], y=C_stock_hist['HOC'], ax=axes[5], label = 'HOC')
  
  plt.savefig('output/RCTM/spinup_stocks_iter_2.jpg', dpi=300)
  
  utils.image_average_variables(spin_covariates, ['GPP'], time_index = 'doy_bins_lower', plot_dir='output/RCTM')
  
#Define parameters - these are calibrated parameters for RCTM v1.0
print('reading in params')
param_path = 'gs://rangelands/Ranch_Runs/HB/covariates_nc/HB_params.tif'
params = rxr.open_rasterio(param_path, masked=True)
params = rctmu.xr_dataset_to_data_array(params)

print('reading in covariates') 
cov_path = 'gs://rangelands/Ranch_Runs/HB/covariates_nc/HB_covariates.nc'
covariates = rxr.open_rasterio(cov_path, masked=True)

C_stocks = {'AGB' : np.full(np.shape(covariates['ndvi'].values[0]), 4.18),
            'BGB' : np.full(np.shape(covariates['ndvi'].values[0]), 30.8),
            'AGL' : np.full(np.shape(covariates['ndvi'].values[0]), 755),
            'BGL' : np.full(np.shape(covariates['ndvi'].values[0]), 1885),
            'POC' : np.full(np.shape(covariates['ndvi'].values[0]), 3473),
            'HOC' : np.full(np.shape(covariates['ndvi'].values[0]), 5161)
            }

covariates, C_stocks = RCTM(C_stocks, covariates, params, spinup_years=100)
  
utils.image_average_variables(covariates, ['GPP', 'Rh', 'Ra', 'NEE', 'NPP'], plot_dir='output/RCTM/transient')

