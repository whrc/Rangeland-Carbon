from google.cloud import storage
import os
import utils
import pandas as pd
import numpy as np
from res import RCTM_params
from sklearn.metrics import mean_squared_error
from scipy.optimize import minimize 
from sklearn.metrics import r2_score
import os
import fsspec
from matplotlib import pyplot as plt
import seaborn as sns

path_to_temp = '/home/amullen/temp/'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/remote-sensing/gee_key.json')

bucket_name='rangelands'
path_to_footprint_list = 'res/sites.txt'
bucket = storage_client.get_bucket(bucket_name)

site_pfts=RCTM_params.ameriflux_site_pfts

#MODIS params for grasslands
# Updated with posteriors - 3
mod_refPars_grasslands = []
mod_refPars_grasslands.append([0.493, 0.45, 1]) #LUEmax
mod_refPars_grasslands.append([-21.219, -65, 2]) #MIN_Tmn
mod_refPars_grasslands.append([10.697, 5, 25]) #MAX_Tmx
mod_refPars_grasslands.append([0.460, 0, 5]) #MIN_VPD
mod_refPars_grasslands.append([0.713, 0, 30]) #MAX_VPD
mod_refPars_grasslands.append([0.021, 0, 0.03]) #MIN_SMrz
mod_refPars_grasslands.append([0.514, 0.5, 1]) #MAX_SMrz
mod_refPars_grasslands.append([3.25, 0.1, 5]) #error-sd
mod_refPars_grasslands=np.array(mod_refPars_grasslands)

#STARFM params for grasslands
# Updated with posteriors - 3
sfm_refPars_grasslands = []
sfm_refPars_grasslands.append([1.01, 1, 2.2]) #LUEmax
sfm_refPars_grasslands.append([-22.32, -65, -5]) #MIN_Tmn
sfm_refPars_grasslands.append([17.419, 10, 35]) #MAX_Tmx
sfm_refPars_grasslands.append([0.335, 0, 2]) #MIN_VPD
sfm_refPars_grasslands.append([2.128, 2, 50]) #MAX_VPD
sfm_refPars_grasslands.append([0.068, 0, 0.3]) #MIN_SMrz
sfm_refPars_grasslands.append([0.994, 0.3, 1]) #MAX_SMrz
sfm_refPars_grasslands.append([3.25, 0.1, 5]) #error-sd
sfm_refPars_grasslands=np.array(sfm_refPars_grasslands)

#STARFM params for grass-shrub
# Updated with posteriors - 3
sfm_refPars_grass_shrub = []
sfm_refPars_grass_shrub.append([1.194, 0.95, 1.5]) #LUEmax
sfm_refPars_grass_shrub.append([-0.918, -15, 5]) #MIN_Tmn
sfm_refPars_grass_shrub.append([10.308, 5, 20]) #MAX_Tmx
sfm_refPars_grass_shrub.append([0.222, 0, 1.8]) #MIN_VPD
sfm_refPars_grass_shrub.append([1.932, 1.8, 30]) #MAX_VPD
sfm_refPars_grass_shrub.append([0.002, 0, 0.05]) #MIN_SMrz
sfm_refPars_grass_shrub.append([0.919, 0.7, 1]) #MAX_SMrz
sfm_refPars_grass_shrub.append([1.64, 0.1, 5]) #error-sd
sfm_refPars_grass_shrub=np.array(sfm_refPars_grass_shrub)

#STARFM params for trees
# Updated with posteriors - 3
sfm_refPars_tree = []
sfm_refPars_tree.append([1.031, 1, 1.6]) #LUEmax
sfm_refPars_tree.append([-31.049, -60, 0]) #MIN_Tmn
sfm_refPars_tree.append([38.676, 25, 40]) #MAX_Tmx
sfm_refPars_tree.append([1.329, 0, 1.9]) #MIN_VPD
sfm_refPars_tree.append([2.049, 1.9, 20]) #MAX_VPD
sfm_refPars_tree.append([0.039, 0, 0.25]) #MIN_SMrz
sfm_refPars_tree.append([0.976, 0.7, 1]) #MAX_SMrz
sfm_refPars_tree.append([1.64, 0.1, 5]) #error-sd
sfm_refPars_tree=np.array(sfm_refPars_tree)

path_to_covariates = 'gs://rangelands/RCTM_GPP_calibration/final_extracted_GPP_covariates_all.csv'
path_to_all_GPP = 'gs://rangelands/RCTM_GPP_calibration/final_measured_GPP.csv'

paths_to_GPP=['/content/gdrive/Shareddrives/Rangeland Carbon/Ongoing projects/Model calibration and validation/RCTM-Calibration/GPP_Cov_Extraction_v4/data.T.csv',
              '/content/gdrive/Shareddrives/Rangeland Carbon/Ongoing projects/Model calibration and validation/RCTM-Calibration/GPP_Cov_Extraction_v4/data.S.csv',
              '/content/gdrive/Shareddrives/Rangeland Carbon/Ongoing projects/Model calibration and validation/RCTM-Calibration/GPP_Cov_Extraction_v4/data.H.csv',
              '/content/gdrive/Shareddrives/Rangeland Carbon/Ongoing projects/Model calibration and validation/RCTM-Calibration/GPP_Cov_Extraction_v4/data.G.csv']
              
starfm_covariates_df = pd.read_csv(path_to_covariates)
starfm_covariates_df['date'] = pd.to_datetime(starfm_covariates_df['date'])
starfm_covariates_df['Site'] = starfm_covariates_df['Site'].str[:-4]
starfm_covariates_df = starfm_covariates_df.rename(columns={"Site": "site", "date": "time"})
starfm_covariates_df = starfm_covariates_df.drop(columns=['evi'])

def calc_sr(index):
  return (1 + index)/(1-index)
  
def calc_bias(observed, predicted):
  return np.nansum((predicted-observed))/len(predicted)
  
def calc_fPAR(this_evi, min_evi, max_evi, this_sr, min_sr, max_sr, max_fpar, min_fpar):
  fpar_EVI = (((this_evi - min_evi)*(max_fpar - min_fpar))/(max_evi - min_evi)) + min_fpar
  
  fpar_SR = (((this_sr - min_sr)*(max_fpar - min_fpar))/(max_sr - min_sr)) + min_fpar
  
  fpar = (fpar_SR + fpar_EVI)/2
  fpar[fpar < 0] = 0
  
  return fpar
  
def get_mult(min, max, in_val):
  out = np.zeros(len(in_val))
  out[in_val>max] = 1 #if val>max, set to 1
  out[(in_val>min) & (in_val<max)] = (in_val[(in_val>min) & (in_val<max)]-min)/(max-min) #if val in between min and max, scale between min and max
  return out
  
def calcGPP(params, tsoil_in, sm_in, vpd_in, SW_IN, fPAR):
  LUEmax = params[0, 0]
  MIN_Tmn = params[1, 0]
  MAX_Tmx = params[2, 0]
  MIN_VPD = params[3, 0]
  MAX_VPD = params[4, 0]
  MIN_SMrz = params[5, 0]
  MAX_SMrz = params[6, 0]
  tmult = get_mult(MIN_Tmn, MAX_Tmx, tsoil_in)
  smult = get_mult(MIN_SMrz, MAX_SMrz, sm_in)
  wmult = get_mult(MIN_VPD, MAX_VPD, vpd_in)

  LUE = LUEmax * tmult * smult * wmult
  GPP = LUE * SW_IN * fPAR * 0.45

  return(GPP)

def get_site_pft(site, site_pfts):
  pft = None
  for key in site_pfts:
      if site in site_pfts[key]:
        pft = key
        
  if pft==None:
    print(f'{site} pft not found')
  return pft
  
##### Read in vegetation index csvs as pandas dataframes #######
df_sites=[]

bad_sites=['Kon']

with open(path_to_footprint_list) as f:
  
  sites = [line.rstrip('\n') for line in f]
  for site in sites:
    print(site)
    if site in bad_sites:
        continue
    
    #read in site csv as pandas dataframe with evi and ndvi columns
    csv_dir = f'Ameriflux_sites/{site}_starfm/covariates/{site}_indices.csv'
    df=pd.read_csv('gs://' + bucket_name + '/' + csv_dir)
    df['time'] = pd.to_datetime(df['time'])
    df['Year'] = df['time'].dt.year
    df['evi_sr'] = calc_sr(df['evi'])
    df['ndvi_sr'] = calc_sr(df['ndvi'])
    df['site'] = site
    
    #calculate yearly reference values for min and max evi/ndvi
    df_groupby_yearly = df.groupby(by=['Year']).agg({'ndvi':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)],
                                                          'ndvi_sr':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)],
                                                          'evi':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)],
                                                          'evi_sr':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)]})
                                                
    df_groupby_yearly = df_groupby_yearly.reset_index()
    df_groupby_yearly.columns = ['Year', 'ndvi_02_yr', 'ndvi_98_yr', 'ndvi_sr_02_yr', 'ndvi_sr_98_yr',
                                 'evi_02_yr', 'evi_98_yr', 'evi_sr_02_yr', 'evi_sr_98_yr']
    ##### special case sites have combined imagery #####
    df = pd.merge(df, df_groupby_yearly, on='Year').reset_index()

    if site=='Rwe_Rms':
      rwe=df
      rwe['site']= 'Rwe'
      rwe['pft'] = get_site_pft('Rwe', site_pfts)
      rms=df.copy()
      rms['site']= 'Rms'
      rms['pft'] = get_site_pft('Rms', site_pfts)
      df_sites.append(rwe)
      df_sites.append(rms)
      
    elif site=='ARbc':
      ARb=df
      ARb['site']= 'ARb'
      ARb['pft'] = get_site_pft('ARb', site_pfts)
      ARc=df.copy()
      ARc['site']= 'ARc'
      ARc['pft'] = get_site_pft('ARc', site_pfts)
      df_sites.append(ARb)
      df_sites.append(ARc)
      
    elif site=='Hn2_3':
      Hn2=df
      Hn2['site']= 'Hn2'
      Hn2['pft'] = get_site_pft('Hn2', site_pfts)
      Hn3=df.copy()
      Hn3['site']= 'Hn3'
      Hn3['pft'] = get_site_pft('Hn3', site_pfts)
      df_sites.append(Hn2)
      df_sites.append(Hn3)
      
    elif site=='KM2_3':
      KM2=df
      KM2['site']= 'KM2'
      KM2['pft'] = get_site_pft('KM2', site_pfts)
      KM3=df.copy()
      KM3['site']= 'KM3'
      KM3['pft'] = get_site_pft('KM3', site_pfts)
      df_sites.append(KM2)
      df_sites.append(KM3)
      
    elif site=='LS1_2':
      LS1=df
      LS1['site']= 'LS1'
      LS1['pft'] = get_site_pft('LS1', site_pfts)
      LS2=df.copy()
      LS2['site']= 'LS2'
      LS2['pft'] = get_site_pft('LS2', site_pfts)
      df_sites.append(LS1)
      df_sites.append(LS2)
      
    elif site=='CZ1_xSJ':
      df['site']= 'xSJ'
      df['pft'] = get_site_pft('xSJ', site_pfts)
      df_sites.append(df)
      
    elif site=='Fpe':
      df['site']= 'FPe'
      df['pft'] = get_site_pft('FPe', site_pfts)
      df_sites.append(df)
      
    else:
      df['pft'] = get_site_pft(site, site_pfts) 
      #merge monthly and yearly reference values back to original df
      df_sites.append(df)
      
      
    
  df_sites=pd.concat(df_sites, ignore_index=True)
  print(df_sites.columns)
df_sites.to_csv('ameriflux_starfm_20230712.csv')  
#read in GPP
GPP_df = pd.read_csv(path_to_all_GPP)
GPP_df['Date'] = pd.to_datetime(GPP_df['Date'])
GPP_df['Site'] = GPP_df['Site'].str[:-4]
GPP_df = GPP_df.rename(columns={"Site": "site", "Date": "time"})

#merge both covariates and GPP to starfm
STARFM_df = pd.merge(df_sites, starfm_covariates_df, how='inner', on=['site', 'time'])
STARFM_df = pd.merge(STARFM_df, GPP_df, how='inner', on=['site', 'time'])

#fig, axes = plt.subplots(3, 1, figsize = (12, 15))
#sns.lineplot(data=GPP_df[GPP_df['time']<='01-01-2005'], x='time', y='GPP_measured', label='Measured', ax=axes[0])
#plt.savefig('GPP_timeseries.jpg')
#STARFM_df.to_csv('starfm_df.csv')


def get_fPAR_ref_yearly(merged_dataframe, ref_pars, method, iterations_per_year=2000):
  
  #min_dist = np.random.uniform(0, 0.15, iterations_per_year) 
  #max_dist = np.random.uniform(0.15, 1, iterations_per_year)
  
  min_dist = np.random.uniform(-0.3, 0.5, iterations_per_year) 
  max_dist = np.random.uniform(0.05, 1.3, iterations_per_year)

  min_swap_vals = min_dist[min_dist>max_dist]
  max_swap_vals = max_dist[max_dist<min_dist]

  min_dist[min_dist>max_dist] = max_swap_vals
  max_dist[max_dist==min_dist] = min_swap_vals

  
  yearly_results = {'Year':[], 'fPAR_min': [], 'fPAR_max': [], 'rmse':[], 'r2': []}

  for year in pd.unique(merged_dataframe['Year']):
    df_merged=merged_dataframe[(merged_dataframe['Year']==year)]
    rmses=[]
    r2s=[]

    for i in range(0, len(min_dist)):
      if method=='EVI':
        STARFM_fPAR=calc_fPAR(df_merged['evi'], 
                            df_merged['evi_02_yr'], 
                            df_merged['evi_98_yr'],
                            df_merged['evi_sr'], 
                            df_merged['evi_sr_02_yr'], 
                            df_merged['evi_sr_98_yr'], 
                            max_dist[i], 
                            min_dist[i])
                             
      if method=='NDVI':
        STARFM_fPAR=calc_fPAR(
                            df_merged['ndvi'], 
                            df_merged['ndvi_02_yr'], 
                            df_merged['ndvi_98_yr'],  
                            df_merged['ndvi_sr'], 
                            df_merged['ndvi_sr_02_yr'], 
                            df_merged['ndvi_sr_98_yr'], 
                            max_dist[i], 
                            min_dist[i])
      STARFM_GPP = calcGPP(ref_pars, 
                        df_merged['TS'], 
                        df_merged['SWC2'], 
                        df_merged['VPD'], 
                        df_merged['SW_IN_NLDAS'], 
                        STARFM_fPAR)
      
      r2 = r2_score(df_merged['GPP_measured'], STARFM_GPP)
      rmse = mean_squared_error(df_merged['GPP_measured'], STARFM_GPP, squared = False)
      rmses.append(rmse)
      r2s.append(r2)

    minimum_rmse_index= np.argmin(rmses)
    maximum_r2_index= np.argmax(r2s)
    
    best_min_val = min_dist[minimum_rmse_index]
    best_max_val = max_dist[minimum_rmse_index]
    min_rmse = rmses[minimum_rmse_index]
    #MODIS_RMSE = mean_squared_error(df_merged['GPP_measured'], df_merged['MODIS_GPP'], squared = False)

    best_min_val_r2 = min_dist[maximum_r2_index]
    best_max_val_r2 = max_dist[maximum_r2_index]
    max_r2 = r2s[maximum_r2_index]
    #MODIS_R2 = r2_score(df_merged['GPP_measured'], df_merged['MODIS_GPP'])
    
    print('YEAR: {}'.format(year))
    print('best fPAR min value: {}\nbest fPAR max value: {}\nminimum RMSE: {}'.format(best_min_val, best_max_val, min_rmse))
    print('best fPAR min value: {}\nbest fPAR max value: {}\nmaximum r2: {}'.format(best_min_val_r2, best_max_val_r2, max_r2))
    print('=======================================================================================')
    
    yearly_results['Year'].append(year)
    yearly_results['fPAR_min'].append(np.nanmean([best_min_val,best_min_val_r2]))
    yearly_results['fPAR_max'].append(np.nanmean([best_max_val,best_max_val_r2]))
    yearly_results['rmse'].append(min_rmse)
    yearly_results['r2'].append(max_r2)
    
  return yearly_results

print(STARFM_df['pft'].unique())
print(STARFM_df[STARFM_df['pft']==None]['site'])  
STARFM_df = STARFM_df[STARFM_df['GPP_measured']>0]
STARFM_df = STARFM_df[STARFM_df['qa']==0]


method='EVI'
pft='grassland'
print(len(STARFM_df[STARFM_df['pft']=='grass']))
#opt_results = get_fPAR_ref_yearly(STARFM_df[STARFM_df['pft']=='grass-shrub'], sfm_refPars_grass_shrub, method=method, iterations_per_year=500)
opt_results = get_fPAR_ref_yearly(STARFM_df[STARFM_df['pft']=='grass'], mod_refPars_grasslands, method=method, iterations_per_year=100)

df_opt_results = pd.DataFrame(opt_results)
df_opt_results = df_opt_results.sort_values(by='Year')
df_STARFM_opt = pd.merge(STARFM_df, df_opt_results, on = ['Year'],how='inner')

if method=='EVI':
  df_STARFM_opt['STARFM_fPAR']=calc_fPAR(df_STARFM_opt['evi'], 
  df_STARFM_opt['evi_02_yr'], 
  df_STARFM_opt['evi_98_yr'],
  df_STARFM_opt['evi_sr'], 
  df_STARFM_opt['evi_sr_02_yr'], 
  df_STARFM_opt['evi_sr_98_yr'], 
  df_STARFM_opt['fPAR_max'], 
  df_STARFM_opt['fPAR_min'])

if method=='NDVI':
  df_STARFM_opt['STARFM_fPAR']=calc_fPAR(df_STARFM_opt['ndvi'], 
  df_STARFM_opt['ndvi_02_yr'], 
  df_STARFM_opt['ndvi_98_yr'], 
  df_STARFM_opt['ndvi_sr'], 
  df_STARFM_opt['ndvi_sr_02_yr'], 
  df_STARFM_opt['ndvi_sr_98_yr'], 
  df_STARFM_opt['fPAR_max'], 
  df_STARFM_opt['fPAR_min'])

df_STARFM_opt['STARFM_GPP'] = calcGPP(mod_refPars_grasslands,
                        df_STARFM_opt['TS'], 
                        df_STARFM_opt['SWC2'], 
                        df_STARFM_opt['VPD'], 
                        df_STARFM_opt['SW_IN_NLDAS'], 
                        df_STARFM_opt['STARFM_fPAR'])


STARFM_r2=r2_score(df_STARFM_opt['GPP_measured'], df_STARFM_opt['STARFM_GPP']).round(2)

STARFM_rmse = mean_squared_error(df_STARFM_opt['GPP_measured'], df_STARFM_opt['STARFM_GPP'], squared = False).round(2)

STARFM_bias = calc_bias(df_STARFM_opt['GPP_measured'], df_STARFM_opt['STARFM_GPP']).round(2)

fig, axes = plt.subplots(3, 1, figsize = (12, 15))

sns.lineplot(data=df_STARFM_opt[df_STARFM_opt['Year']<=2004], x='time', y='STARFM_GPP', label='STARFM', ax=axes[0])
sns.lineplot(data=df_STARFM_opt[df_STARFM_opt['Year']<=2004], x='time', y='GPP_measured', label='Measured', ax=axes[0])

sns.lineplot(data=df_STARFM_opt[(df_STARFM_opt['Year']>2004) & (df_STARFM_opt['Year']<=2008)], x='time', y='STARFM_GPP', label='STARFM', ax=axes[1])
sns.lineplot(data=df_STARFM_opt[(df_STARFM_opt['Year']>2004) & (df_STARFM_opt['Year']<=2008)], x='time', y='GPP_measured', label='Measured', ax=axes[1])

sns.lineplot(data=df_STARFM_opt[(df_STARFM_opt['Year']>2008) & (df_STARFM_opt['Year']<=2010)], x='time', y='STARFM_GPP', label='STARFM', ax=axes[2])
sns.lineplot(data=df_STARFM_opt[(df_STARFM_opt['Year']>2008) & (df_STARFM_opt['Year']<=2010)], x='time', y='GPP_measured', label='Measured', ax=axes[2])

axes[0].set_ylabel('GPP')
axes[1].set_ylabel('GPP')
axes[2].set_ylabel('GPP')

axes[0].set_title('STARFM fPAR calculated using {} and yearly optimized fPAR min and max for {} sites'.format(method, pft), color='grey')

axes[0].text(0.05, 0.9,"error against observed GPP",transform=axes[0].transAxes, fontsize=12, fontweight='bold')
axes[0].text(0.05, 0.85,"STARFM rmse: {}, STARFM r2: {}, STARFM bias: {}".format(STARFM_rmse, STARFM_r2, STARFM_bias),transform=axes[0].transAxes)

fig.tight_layout()

plt.savefig('GPP_{}_method_{}_sfm_cal.jpg'.format(method, pft), dpi=300)

print('STARFM r2: {}'.format(STARFM_r2))
#get_fPAR_ref_yearly(STARFM_df[STARFM_df['pft']=='tree'], sfm_refPars_grass_shrub, method='EVI', iterations_per_year=1000)
#get_fPAR_ref_yearly(STARFM_df[STARFM_df['pft']=='grass-shrub'], sfm_refPars_tree, method='EVI', iterations_per_year=1000)
#print(set(STARFM_df['site'])-set(GPP_df['site']))
#print(set(GPP_df['site'])-set(STARFM_df['site']))

