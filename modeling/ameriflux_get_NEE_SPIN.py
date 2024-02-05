import fsspec
from matplotlib import pyplot as plt
import seaborn as sns
import RCTM_utils as rctmu
import yaml
from sklearn.cluster import KMeans
from sklearn import preprocessing
import pickle
from sklearn.cluster import KMeans
from sklearn import preprocessing
import pickle
from google.cloud import storage
import RCTM_utils as rctmu
import pandas as pd
import numpy as np
from datetime import date, datetime, timedelta
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import sys
sys.path.insert(1, '../utils')
import utils

bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
bucket = storage_client.get_bucket(bucket_name)


path_to_temp = '/home/amullen/temp/'
path_to_footprint_list = '/home/amullen/Rangeland-Carbon/res/site_footprints/sites.txt'

params=[]
with open('RCTM_params.yaml', 'r') as file:
  params = yaml.safe_load(file)
  
ref_value_calc = 'all'

### read NEE/spin covariates              
SPIN_dir='Ameriflux/SPIN_AVG_ALL/'
NEE_dir='Ameriflux/NEEINS_AVG_ALL/'

NEE_meas_dir = 'gs://' + bucket_name + '/' + 'Ameriflux/final_measured_NEE'

NEE_meas_df = pd.read_csv(NEE_meas_dir, parse_dates=['Date'])
NEE_meas_df = NEE_meas_df.rename(columns={'Site': 'site', 'Date': 'time'})

NEE_meas_df['site'] = NEE_meas_df['site'].str[:3]

NEE_dfs = []

for filename in utils.gs_listdir(NEE_dir, bucket):
  if filename.endswith('.csv'):
    site_df = pd.read_csv('gs://' + bucket_name + '/' + filename)
    site_df['system:index'] = pd.to_datetime(site_df['system:index'], format='%Y_%m_%d')
    site_df=site_df.drop(columns=['date'])
    site_df=site_df.rename(columns={'system:index': 'time'})
    
  
    site_df['site'] = [filename.split('/')[-1][:3]]*len(site_df)
    
    NEE_dfs.append(site_df)

NEE_df = pd.concat(NEE_dfs)

spin_dfs = []
#loop through site data, calculate reference values, append to df, save to new csv
for filename in utils.gs_listdir(SPIN_dir, bucket):
  if filename.endswith('.csv'):
    site_df = pd.read_csv('gs://' + bucket_name + '/' + filename)
    site_df=site_df.drop(columns=['date'])
    site_df['system:index'] = pd.to_datetime(site_df['system:index'], format='%Y_%m_%d')
    site_df=site_df.rename(columns={'system:index': 'time'})
    
  
    site_df['site'] = [filename.split('/')[-1][:3]]*len(site_df)
    
    spin_dfs.append(site_df)

spin_df = pd.concat(spin_dfs)

##### Read in vegetation index csvs as pandas dataframes #######
df_sites = rctmu.read_in_sites_as_df(path_to_footprint_list, params, bucket_name, ref_value_calc=ref_value_calc)

df_sites['ndvi_02_yr'] = df_sites['ndvi'].quantile(0.02)
df_sites['ndvi_98_yr'] = df_sites['ndvi'].quantile(0.98)
df_sites['ndvi_sr_02_yr'] = df_sites['ndvi_sr'].quantile(0.02)
df_sites['ndvi_sr_98_yr'] = df_sites['ndvi_sr'].quantile(0.98) 

print('ndvi_02: {}'.format(df_sites['ndvi'].quantile(0.02)))
print('ndvi_98: {}'.format(df_sites['ndvi'].quantile(0.98)))
print('ndvi_sr_02: {}'.format(df_sites['ndvi_sr'].quantile(0.02)))
print('ndvi_sr_98: {}'.format(df_sites['ndvi_sr'].quantile(0.98)))
 
#create 5-day timestep from STARFM df
print('creating 5-day time series')
df_STARFM = []
for site in df_sites['site'].unique():
  
  df_site = df_sites.loc[df_sites['site']==site]
  df_site = df_site.sort_values(by='time')
  
  pft = df_site['pft'].iloc[0]
  first_date = df_site['time'].min()
  last_date = df_site['time'].max()

  df_site = df_site.sort_values(by='time')
  
  date_range = pd.DataFrame({'time': pd.date_range(start='2002-01-01', end='2023-01-01', freq='5d').to_series()})

  df_site = date_range.merge(df_site, how='left', on='time')
  
  curr = pd.to_datetime('2002-01-01')
  every_five = [curr]
  
  while curr <= pd.to_datetime('2023-01-01'):

        curr += timedelta(days=5)

        every_five.append(curr)
        
  df_site = df_site.loc[df_site['time'].isin(every_five)]
  df_site['site'] = site
  df_site['pft'] = pft
  
  
  df_STARFM.append(df_site)
  print(f'length of {site} time series: {len(df_site)}')
  
df_STARFM=pd.concat(df_STARFM, ignore_index=True)
 
#merge both covariates and GPP to starfm
print('merging NEE dataset with STARFM')
NEE_df = pd.merge(df_STARFM.loc[df_STARFM['site'].isin(NEE_df['site'].unique())], NEE_df, how='left', on=['site', 'time']) 

nee_covariates__grouped = NEE_df.groupby(by=['site']).agg({'VPD':np.nanmean})                                                
nee_covariates__grouped = nee_covariates__grouped.reset_index()

nee_covariates__grouped.columns = ['site', 'vpd_mean']
NEE_df = pd.merge(NEE_df, nee_covariates__grouped, on='site').reset_index()

#cluster
NEE_df['cluster'] = 0
#calc fPAR and GPP
print('calculating fPAR and GPP for dates with data (NEE)')
NEE_df['fPAR_min'] = [params['RCTM_params'][row['pft']]['starfm']['fPAR_min'][row['cluster']] for index, row in NEE_df.iterrows()]
NEE_df['fPAR_max'] = [params['RCTM_params'][row['pft']]['starfm']['fPAR_max'][row['cluster']] for index, row in NEE_df.iterrows()]

NEE_df['STARFM_fPAR']=rctmu.calc_fPAR(
                            NEE_df['ndvi'], 
                            NEE_df['ndvi_02_yr'], 
                            NEE_df['ndvi_98_yr'],  
                            NEE_df['ndvi_sr'], 
                            NEE_df['ndvi_sr_02_yr'], 
                            NEE_df['ndvi_sr_98_yr'], 
                            NEE_df['fPAR_min'], 
                            NEE_df['fPAR_max'])

for pft in NEE_df['pft'].unique():                                                       
  NEE_df.loc[NEE_df['pft']==pft, 'STARFM_GPP'] = rctmu.calcGPP(params['RCTM_params'][pft]['starfm'], 
                           NEE_df.loc[NEE_df['pft']==pft,'TS'], 
                           NEE_df.loc[NEE_df['pft']==pft,'SWC2'], 
                           NEE_df.loc[NEE_df['pft']==pft,'VPD'], 
                           NEE_df.loc[NEE_df['pft']==pft,'SW_IN_NLDAS'], 
                           NEE_df.loc[NEE_df['pft']==pft,'STARFM_fPAR'])

print('merging spin data')                
spin_df = pd.merge(df_STARFM.loc[df_STARFM['site'].isin(spin_df['site'].unique())], spin_df, how='left', on=['site', 'time'])

spin_covariates__grouped = spin_df.groupby(by=['site']).agg({'VPD':np.nanmean}) 
                                              
spin_covariates__grouped = spin_covariates__grouped.reset_index()

spin_covariates__grouped.columns = ['site', 'vpd_mean']
spin_df = pd.merge(spin_df, spin_covariates__grouped, on='site').reset_index()

#cluster
#scaler = preprocessing.StandardScaler().fit(spin_df['vpd_mean'].to_numpy().reshape(-1, 1))
#X_scaled = scaler.transform(spin_df['vpd_mean'].to_numpy().reshape(-1, 1))
#loaded_model = pickle.load(open('output/fPAR_optimization/clusters/kmeans.pickle', "rb"))
#spin_df['cluster'] = loaded_model.predict(X_scaled)
##spin_df = pd.merge(spin_df, cluster_df, on='site')

spin_df['cluster'] = 0
print('caluculating fPAR and GPP for dates with data (spin)')
spin_df['fPAR_min'] = [params['RCTM_params'][row['pft']]['starfm']['fPAR_min'][row['cluster']] for index, row in spin_df.iterrows()]
spin_df['fPAR_max'] = [params['RCTM_params'][row['pft']]['starfm']['fPAR_max'][row['cluster']] for index, row in spin_df.iterrows()]
spin_df['STARFM_fPAR']=rctmu.calc_fPAR(
                            spin_df['ndvi'], 
                            spin_df['ndvi_02_yr'], 
                            spin_df['ndvi_98_yr'],  
                            spin_df['ndvi_sr'], 
                            spin_df['ndvi_sr_02_yr'], 
                            spin_df['ndvi_sr_98_yr'], 
                            spin_df['fPAR_min'], 
                            spin_df['fPAR_max'])

for pft in spin_df['pft'].unique():                                                       
  spin_df.loc[spin_df['pft']==pft, 'STARFM_GPP'] = rctmu.calcGPP(params['RCTM_params'][pft]['starfm'], 
                           spin_df.loc[spin_df['pft']==pft,'TS'], 
                           spin_df.loc[spin_df['pft']==pft,'SWC2'], 
                           spin_df.loc[spin_df['pft']==pft,'VPD'], 
                           spin_df.loc[spin_df['pft']==pft,'SW_IN_NLDAS'], 
                           spin_df.loc[spin_df['pft']==pft,'STARFM_fPAR'])
NEE_df['STARFM_fPAR_min'] = np.nan
NEE_df['STARFM_fPAR_max'] = np.nan

for pft in NEE_df['pft'].unique():
  NEE_df.loc[NEE_df['pft']==pft, 'STARFM_fPAR_min'] = NEE_df.loc[NEE_df['pft']==pft, 'STARFM_fPAR'].quantile(0.02)
  NEE_df.loc[NEE_df['pft']==pft, 'STARFM_fPAR_max'] = NEE_df.loc[NEE_df['pft']==pft, 'STARFM_fPAR'].quantile(0.98)
  
model, scaler = utils.predict_column_with_neural_network(NEE_df, ['TS', 'SWC1', 'SWC2', 'VPD', 'SW_IN_NLDAS', 'STARFM_fPAR_min', 'STARFM_fPAR_max'], 'STARFM_fPAR')

print('fitting PCA to sites (NEE)') 
for site in NEE_df['site'].unique():
  print(NEE_df.loc[NEE_df['site']==site,'time'].min())
  #r2, model, pca = utils.principal_component_regression(NEE_df.loc[(NEE_df['site']==site)], ['TS', 'SWC2', 'VPD', 'SW_IN_NLDAS'], 'STARFM_fPAR')
  #print(r2)
  #interpolate covariates
  NEE_df = NEE_df.sort_values(by='time')
  NEE_df.loc[NEE_df['site']==site, 'Clay'] = NEE_df.loc[NEE_df['site']==site, 'Clay'].interpolate()
  NEE_df.loc[NEE_df['site']==site, 'TS'] = NEE_df.loc[NEE_df['site']==site, 'TS'].interpolate()
  NEE_df.loc[NEE_df['site']==site, 'SWC1'] = NEE_df.loc[NEE_df['site']==site, 'SWC1'].interpolate()
  NEE_df.loc[NEE_df['site']==site, 'SWC2'] = NEE_df.loc[NEE_df['site']==site, 'SWC2'].interpolate()
  NEE_df.loc[NEE_df['site']==site, 'VPD'] = NEE_df.loc[NEE_df['site']==site, 'VPD'].interpolate()
  NEE_df.loc[NEE_df['site']==site, 'SW_IN_NLDAS'] = NEE_df.loc[NEE_df['site']==site, 'SW_IN_NLDAS'].interpolate()
  
  #interpolate fPAR within 5 observations (25 days)
  NEE_df.loc[NEE_df['site']==site, 'STARFM_fPAR'] = NEE_df.loc[NEE_df['site']==site, 'STARFM_fPAR'].interpolate(limit=50)
  
  #predict fPAR using PCA if gap is longer than 25 days
  if len(NEE_df.loc[(NEE_df['site']==site) & (NEE_df['STARFM_fPAR'].isna()) & (~NEE_df['SWC2'].isna())])>0:
    #NEE_df.loc[(NEE_df['site']==site) & (NEE_df['STARFM_fPAR'].isna()) & (~NEE_df['SWC2'].isna()), 'STARFM_fPAR'] = model.predict(pca.transform(NEE_df.loc[(NEE_df['site']==site) & (NEE_df['STARFM_fPAR'].isna()) & (~NEE_df['SWC2'].isna()), ['TS', 'SWC2', 'VPD', 'SW_IN_NLDAS']]))
    NEE_df.loc[(NEE_df['site']==site) & (NEE_df['STARFM_fPAR'].isna()) & (~NEE_df['SWC2'].isna()), 'STARFM_fPAR'] = model.predict(scaler.transform(NEE_df.loc[(NEE_df['site']==site) & (NEE_df['STARFM_fPAR'].isna()) & (~NEE_df['SWC2'].isna()), ['TS', 'SWC1', 'SWC2', 'VPD', 'SW_IN_NLDAS', 'STARFM_fPAR_min', 'STARFM_fPAR_max']]))
  
  #set fPAR < 0 to 0  
  NEE_df.loc[(NEE_df['site']==site) & (NEE_df['STARFM_fPAR']<0), 'STARFM_fPAR'] = 0  

#filter nan
NEE_df = NEE_df.loc[~NEE_df['SWC2'].isna()]  

#filter nan
NEE_df.loc[NEE_df['STARFM_fPAR']>1, 'STARFM_fPAR'] = 1 

#merge NEE covariate dataframe with observed NEE
NEE_df = pd.merge(NEE_df, NEE_meas_df, on = ['site', 'time'], how='left')

#filter dates
NEE_df = NEE_df.loc[(NEE_df['time']>='2002-01-01') & (NEE_df['time']<'2022-05-01')]
  
print('calculating gap-filled GPP (NEE)')
#calculate GPP for each pft with gap-filled fPAR
for pft in NEE_df['pft'].unique():
  print('pft length: ' + str(len(NEE_df.loc[NEE_df['pft']==pft, 'STARFM_GPP'])))
  print('pft length GPP NA: ' + str(len(NEE_df.loc[(NEE_df['pft']==pft) & (NEE_df['STARFM_GPP'].isna()), 'STARFM_GPP'])))
  print('pft length fPAR NA: ' + str(len(NEE_df.loc[(NEE_df['pft']==pft) & (NEE_df['STARFM_fPAR'].isna()), 'STARFM_fPAR']))) 
                                                
  NEE_df.loc[NEE_df['pft']==pft, 'STARFM_GPP'] = rctmu.calcGPP(params['RCTM_params'][pft]['starfm'], 
                           NEE_df.loc[NEE_df['pft']==pft,'TS'], 
                           NEE_df.loc[NEE_df['pft']==pft,'SWC2'], 
                           NEE_df.loc[NEE_df['pft']==pft,'VPD'], 
                           NEE_df.loc[NEE_df['pft']==pft,'SW_IN_NLDAS'], 
                           NEE_df.loc[NEE_df['pft']==pft,'STARFM_fPAR'])
                           
  print('pft length GPP NA post calc: ' + str(len(NEE_df.loc[(NEE_df['pft']==pft) & (NEE_df['STARFM_GPP'].isna()), 'STARFM_GPP'])))
  
NEE_df = NEE_df.loc[~NEE_df['SWC2'].isna()] 
                           
for site in NEE_df['site'].unique():
  print(f'length of gap-filled {site} nan: ' + str(len(NEE_df.loc[(NEE_df['site']==site) & (NEE_df['STARFM_GPP'].isna())])))
  print(f'length of gap-filled {site} data without nan: ' + str(len(NEE_df.loc[(NEE_df['site']==site) & (~NEE_df['STARFM_GPP'].isna())])))
  
  print('length of spinup: {}'.format(len(NEE_df.loc[(NEE_df['site']==site) & (NEE_df['Year']>=2002) & (NEE_df['Year']<=2005)])))
  print(NEE_df.loc[(NEE_df['site']==site), 'Year'].min())
  
# save
NEE_df['Year']=NEE_df['time'].dt.year

NEE_df = NEE_df.sort_values(by=['site', 'time']) 
NEE_df = NEE_df.rename(columns={'site': 'Site', 'time':'Date'})
for pft in params['Ameriflux_site_pfts']:
  NEE_df.loc[NEE_df['Site'].isin(params['Ameriflux_site_pfts'][pft])].to_csv('output/NEE/Ameriflux_STARFM_NEE_{}.csv'.format(pft))

### generate spinup data
spin_df = NEE_df.loc[(NEE_df['Year']>=2002) & (NEE_df['Year']<=2005)]

spin_df = spin_df.drop(columns = ['.geo'])

spin_df.loc[spin_df['Date']>pd.to_datetime('2004-02-28'), 'Date'] = spin_df.loc[spin_df['Date']>pd.to_datetime('2004-02-28'), 'Date'] + pd.to_timedelta(1, unit = 'D')

spin_df['Month'] = spin_df['Date'].dt.month
spin_df['Day'] = spin_df['Date'].dt.day
spin_df = spin_df.drop(columns = ['Date'])

spin_df = spin_df.groupby(by=['Month', 'Day', 'pft', 'Site']).mean().reset_index()

spin_df['Date'] = pd.to_datetime(spin_df[['Year', 'Month', 'Day']])

for site in spin_df['Site'].unique():
  print(len(spin_df.loc[spin_df['Site']==site]))


spin_df = spin_df.sort_values(by=['Site', 'Date'])

for pft in params['Ameriflux_site_pfts']:

  spin_df.loc[spin_df['Site'].isin(params['Ameriflux_site_pfts'][pft])].to_csv('output/SPIN/Ameriflux_STARFM_spin_{}.csv'.format(pft))