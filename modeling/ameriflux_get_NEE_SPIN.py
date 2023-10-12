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
import utils
import RCTM_utils as RCTMU
import pandas as pd
import numpy as np
from datetime import date, datetime, timedelta
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/remote-sensing/gee_key.json')
bucket = storage_client.get_bucket(bucket_name)


path_to_temp = '/home/amullen/temp/'
path_to_footprint_list = 'res/sites.txt'

params=[]
with open('res/RCTM_params.yaml', 'r') as file:
  params = yaml.safe_load(file)

### read NEE/spin covariates              
SPIN_dir='Ameriflux/SPIN_AVG_ALL/'
NEE_dir='Ameriflux/NEEINS_AVG_ALL/'

NEE_meas_dir = 'gs://' + bucket_name + '/' + 'Ameriflux/final_measured_NEE'

NEE_meas_df = pd.read_csv(NEE_meas_dir, parse_dates=['Date'])
NEE_meas_df = NEE_meas_df.rename(columns={'Site': 'site', 'Date': 'time'})

NEE_meas_df['site'] = NEE_meas_df['site'].str[:3]

path_to_clusters='output/fPAR_optimization/clusters/STARFM_fPAR_optimization_vpd_clusters.csv'
cluster_df = pd.read_csv(path_to_clusters)
cluster_df['cluster'] = cluster_df['cluster'].astype(int)

NEE_dfs = []
#loop through site data, calculate reference values, append to df, save to new csv
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
print('reading in STARFM data')
df_sites=[]

bad_sites=['Kon']

with open(path_to_footprint_list) as f:
  sites = [line.rstrip('\n') for line in f]
  
  for site in sites:
    site=site.split('/')[-1]
    print(site)
    if site in bad_sites:
        continue
    
    #read in site csv as pandas dataframe with evi and ndvi columns
    csv_dir = f'Ameriflux_sites/{site}_starfm/covariates_v2/{site}_indices_v2.csv'
    #csv_dir = f'Ameriflux_sites/{site}_starfm/covariates_smooth_missing/{site}_indices.csv'
    df=pd.read_csv('gs://' + bucket_name + '/' + csv_dir)
    df['time'] = pd.to_datetime(df['time'])
    df['Year'] = df['time'].dt.year
    #df['evi_sr'] = rctmu.calc_sr(df['evi'])
    df['ndvi_sr'] = rctmu.calc_sr(df['ndvi'])
    df['site'] = site
    
    #calculate yearly reference values for min and max evi/ndvi
    df_groupby_yearly = df.groupby(by=['site']).agg({'ndvi':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)],
                                                          'ndvi_sr':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)]})
                                                          #'evi':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)],
                                                          #'evi_sr':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)]})

                                                
    df_groupby_yearly = df_groupby_yearly.reset_index()
    df_groupby_yearly.columns = ['site', 'ndvi_02_yr', 'ndvi_98_yr', 'ndvi_sr_02_yr', 'ndvi_sr_98_yr']
                                 #'evi_02_yr', 'evi_98_yr', 'evi_sr_02_yr', 'evi_sr_98_yr']
    ##### special case sites have combined imagery #####
    df = pd.merge(df, df_groupby_yearly, on='site').reset_index()

    if site=='Rwe_Rms':
      rwe=df
      rwe['site']= 'Rwe'
      rwe['pft'] = rctmu.get_site_pft('Rwe', params['Ameriflux_site_pfts'])
      rms=df.copy()
      rms['site']= 'Rms'
      rms['pft'] = rctmu.get_site_pft('Rms', params['Ameriflux_site_pfts'])
      df_sites.append(rwe)
      df_sites.append(rms)
      
    elif site=='ARbc':
      ARb=df
      ARb['site']= 'ARb'
      ARb['pft'] = rctmu.get_site_pft('ARb', params['Ameriflux_site_pfts'])
      ARc=df.copy()
      ARc['site']= 'ARc'
      ARc['pft'] = rctmu.get_site_pft('ARc', params['Ameriflux_site_pfts'])
      df_sites.append(ARb)
      df_sites.append(ARc)
      
    elif site=='Hn2_3':
      Hn2=df
      Hn2['site']= 'Hn2'
      Hn2['pft'] = rctmu.get_site_pft('Hn2', params['Ameriflux_site_pfts'])
      Hn3=df.copy()
      Hn3['site']= 'Hn3'
      Hn3['pft'] = rctmu.get_site_pft('Hn3', params['Ameriflux_site_pfts'])
      df_sites.append(Hn2)
      df_sites.append(Hn3)
      
    elif site=='KM2_3':
      KM2=df
      KM2['site']= 'KM2'
      KM2['pft'] = rctmu.get_site_pft('KM2', params['Ameriflux_site_pfts'])
      KM3=df.copy()
      KM3['site']= 'KM3'
      KM3['pft'] = rctmu.get_site_pft('KM3', params['Ameriflux_site_pfts'])
      df_sites.append(KM2)
      df_sites.append(KM3)
      
    elif site=='LS1_2':
      LS1=df
      LS1['site']= 'LS1'
      LS1['pft'] = rctmu.get_site_pft('LS1', params['Ameriflux_site_pfts'])
      LS2=df.copy()
      LS2['site']= 'LS2'
      LS2['pft'] = rctmu.get_site_pft('LS2', params['Ameriflux_site_pfts'])
      df_sites.append(LS1)
      df_sites.append(LS2)
      
    elif site=='CZ1_xSJ':
      df['site']= 'xSJ'
      df['pft'] = rctmu.get_site_pft('xSJ', params['Ameriflux_site_pfts'])
      df_sites.append(df)
      
    elif site=='Fpe':
      df['site']= 'FPe'
      df['pft'] = rctmu.get_site_pft('FPe', params['Ameriflux_site_pfts'])
      df_sites.append(df)
      
    else:
      df['pft'] = rctmu.get_site_pft(site, params['Ameriflux_site_pfts']) 
    
      df_sites.append(df)
    
      
  #merge monthly and yearly reference values back to original df
  df_sites=pd.concat(df_sites, ignore_index=True)

 
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
#scaler = preprocessing.StandardScaler().fit(NEE_df['vpd_mean'].to_numpy().reshape(-1, 1))
#X_scaled = scaler.transform(NEE_df['vpd_mean'].to_numpy().reshape(-1, 1))
#loaded_model = pickle.load(open('output/fPAR_optimization/clusters/kmeans.pickle', "rb"))
#NEE_df['cluster'] = loaded_model.predict(X_scaled)
NEE_df = pd.merge(NEE_df, cluster_df, on='site')
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
spin_df = pd.merge(spin_df, cluster_df, on='site')
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
  NEE_df.loc[NEE_df['site']==site, 'Clay'] = NEE_df.loc[NEE_df['site']==site, 'SWC2'].interpolate()
  NEE_df.loc[NEE_df['site']==site, 'TS'] = NEE_df.loc[NEE_df['site']==site, 'TS'].interpolate()
  NEE_df.loc[NEE_df['site']==site, 'SWC1'] = NEE_df.loc[NEE_df['site']==site, 'SWC2'].interpolate()
  NEE_df.loc[NEE_df['site']==site, 'SWC2'] = NEE_df.loc[NEE_df['site']==site, 'SWC2'].interpolate()
  NEE_df.loc[NEE_df['site']==site, 'VPD'] = NEE_df.loc[NEE_df['site']==site, 'VPD'].interpolate()
  NEE_df.loc[NEE_df['site']==site, 'SW_IN_NLDAS'] = NEE_df.loc[NEE_df['site']==site, 'SW_IN_NLDAS'].interpolate()
  
  #interpolate fPAR within 5 observations (25 days)
  NEE_df.loc[NEE_df['site']==site, 'STARFM_fPAR'] = NEE_df.loc[NEE_df['site']==site, 'STARFM_fPAR'].interpolate(limit=5)
  
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
  NEE_df.loc[NEE_df['Site'].isin(params['Ameriflux_site_pfts'][pft])].to_csv('output/NEE_SPIN/Ameriflux_STARFM_NEE_{}.csv'.format(pft))

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

  spin_df.loc[spin_df['Site'].isin(params['Ameriflux_site_pfts'][pft])].to_csv('output/NEE_SPIN/Ameriflux_STARFM_spin_{}.csv'.format(pft))