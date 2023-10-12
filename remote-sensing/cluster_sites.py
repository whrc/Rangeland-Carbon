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

bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/remote-sensing/gee_key.json')
bucket = storage_client.get_bucket(bucket_name)

path_to_clusters='output/fPAR_optimization/clusters/STARFM_fPAR_optimization_vpd_clusters.csv'
NEE_dir='Ameriflux/NEEINS_AVG_ALL/'

path_to_footprint_list = 'res/sites.txt'

params=[]
with open('res/RCTM_params.yaml', 'r') as file:
  params = yaml.safe_load(file)

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
    csv_dir = f'Ameriflux_sites/{site}_starfm/covariates_smooth_all_modis/{site}_indices.csv'
    #csv_dir = f'Ameriflux_sites/{site}_starfm/covariates_smooth_missing/{site}_indices.csv'
    df=pd.read_csv('gs://' + bucket_name + '/' + csv_dir)
    df['time'] = pd.to_datetime(df['time'])
    df['Year'] = df['time'].dt.year
    df['ndvi_sr'] = rctmu.calc_sr(df['ndvi'])
    df['site'] = site
    
    #calculate yearly reference values for min and max evi/ndvi
    df_groupby_yearly = df.groupby(by=['site']).agg({'ndvi':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)],
                                                          'ndvi_sr':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)]})

                                                
    df_groupby_yearly = df_groupby_yearly.reset_index()
    df_groupby_yearly.columns = ['site', 'ndvi_02_site', 'ndvi_98_site', 'ndvi_sr_02_site', 'ndvi_sr_98_site',
                                 ]
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
  
NEE_df = pd.merge(NEE_df, df_sites, how='left', on=['site', 'time'])

#get mean vpd per site
nee_covariates__grouped = NEE_df.groupby(by=['site']).agg({'VPD':np.nanmean,
                                                            'evi_x':np.nanmean,
                                                            'ndvi_02_site': np.nanmean,
                                                            'ndvi_98_site': np.nanmean,
                                                            'ndvi': np.nanmean}) 
#nee_covariates__grouped = NEE_df.groupby(by=['site']).agg({'evi':np.nanmean})                                               
nee_covariates__grouped = nee_covariates__grouped.reset_index()
nee_covariates__grouped.columns = ['site', 'vpd_mean', 'evi_mean', 'ndvi_02_site_mean', 'ndvi_98_site_mean', 'ndvi_mean']
NEE_df = pd.merge(NEE_df, nee_covariates__grouped, on='site').reset_index()

NEE_df = NEE_df.loc[~NEE_df['ndvi_02_site'].isna()]

#cluster
def cluster_by_attributes(df, attributes, n_clusters):

  scaler = preprocessing.StandardScaler().fit(df[attributes])
  X_scaled = scaler.transform(df[attributes])
  
  model = KMeans(n_clusters=n_clusters, n_init=10)
  model.fit(X_scaled)
  predicted_labels = model.predict(X_scaled)
  df['cluster'] = predicted_labels
  return df, model

print('clustering')
NEE_df, model = cluster_by_attributes(NEE_df, ['ndvi_02_site_mean', 'ndvi_98_site_mean'], n_clusters=4)
clusters = NEE_df[['site', 'cluster']].groupby(by='site').mean().reset_index()

pickle.dump(model, open('output/fPAR_optimization/clusters/kmeans.pickle', "wb"))
clusters.to_csv('output/fPAR_optimization/clusters/STARFM_fPAR_optimization_vpd_clusters.csv')