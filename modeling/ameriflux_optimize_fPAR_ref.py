from google.cloud import storage
import os
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from scipy.optimize import minimize 
from sklearn.metrics import r2_score
import os
import fsspec
from matplotlib import pyplot as plt
import seaborn as sns
import RCTM_utils as rctmu
import yaml
from sklearn.cluster import KMeans
from sklearn import preprocessing
import pickle
import sys
sys.path.insert(1, '../utils')
import utils

#Google Cloud Storage config
bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
bucket = storage_client.get_bucket(bucket_name)

path_to_temp = '/home/amullen/temp/' #temp directory for saving files to be uploaded to GCloud
path_to_footprint_list = '/home/amullen/Rangeland-Carbon/res/site_footprints/sites.txt' #footprints to process, also points to appropriate folders in GCloud

path_to_covariates = 'gs://rangelands/RCTM_GPP_calibration/final_extracted_GPP_covariates_all.csv' #GPP model covariates
path_to_all_GPP = 'gs://rangelands/RCTM_GPP_calibration/final_measured_GPP.csv' #measured GPP

#open model param YAML file
params=[]
with open('RCTM_params.yaml', 'r') as file:
  params = yaml.safe_load(file)
  
ref_value_calc = 'all'
iterations = 5000

### read GPP covariates              
starfm_covariates_df = pd.read_csv(path_to_covariates)
starfm_covariates_df['date'] = pd.to_datetime(starfm_covariates_df['date'])
starfm_covariates_df['Site'] = starfm_covariates_df['Site'].str[:-4]
starfm_covariates_df = starfm_covariates_df.rename(columns={"Site": "site", "date": "time"})

# average covariates for each site time series
starfm_covariates__grouped = starfm_covariates_df.groupby(by=['site']).agg({'TS':'mean',
                                                                            'SWC2':'mean',
                                                                            'VPD':'mean',
                                                                            'SW_IN_NLDAS':'mean'})                                                      
starfm_covariates__grouped = starfm_covariates__grouped.reset_index()
starfm_covariates__grouped.columns = ['site', 'ts_mean', 'swc_mean', 'vpd_mean', 'swin_mean']
starfm_covariates_df = pd.merge(starfm_covariates_df, starfm_covariates__grouped, on='site').reset_index()
  
##### Read in vegetation index csvs as pandas dataframes #######
df_sites = rctmu.read_in_sites_as_df(path_to_footprint_list, params, bucket_name, ref_value_calc=ref_value_calc)
 
#read in GPP
GPP_df = pd.read_csv(path_to_all_GPP)
GPP_df['Date'] = pd.to_datetime(GPP_df['Date'])
GPP_df['Site'] = GPP_df['Site'].str[:-4]
GPP_df = GPP_df.rename(columns={"Site": "site", "Date": "time"})

#merge both covariates and GPP to starfm
STARFM_df = pd.merge(df_sites, starfm_covariates_df, how='inner', on=['site', 'time'])
STARFM_df = pd.merge(STARFM_df, GPP_df, how='inner', on=['site', 'time'])

def get_global_fPAR_ref(merged_dataframe, params, iterations=2000, seed=2023):
  
  np.random.seed(seed)
  
  #fPAR min-max distributions to select from
  min_dist = np.random.uniform(-0.3, 0.5, iterations) 
  max_dist = np.random.uniform(0.05, 1.3, iterations)
  
  #if there is a min value > max value, swap their values in the sample arrays
  min_swap_vals = min_dist[min_dist>max_dist]
  max_swap_vals = max_dist[max_dist<min_dist]
  min_dist[min_dist>max_dist] = max_swap_vals
  max_dist[max_dist==min_dist] = min_swap_vals
  
  #dictionary to store results
  results = {'cluster': [], 'fPAR_min': [], 'fPAR_max': [], 'rmse':[], 'r2': []}

  rmses=[]
  r2s=[]

  for i in range(0, len(min_dist)):
    
    merged_dataframe['STARFM_fPAR'] = np.nan
    merged_dataframe['STARFM_GPP'] = np.nan
    
    pft_rmses=[]
    pft_r2s=[] 
    lens=[]
    for pft in merged_dataframe['pft'].unique(): 
      
      
      
      df_pft=merged_dataframe.loc[merged_dataframe['pft']==pft]
      #calculate fPAR with sampled fPAR min and max  
      merged_dataframe.loc[merged_dataframe['pft']==pft,'STARFM_fPAR']=rctmu.calc_fPAR(
                             df_pft['ndvi'], 
                             df_pft['ndvi_02_yr'], 
                             df_pft['ndvi_98_yr'],  
                             df_pft['ndvi_sr'], 
                             df_pft['ndvi_sr_02_yr'], 
                             df_pft['ndvi_sr_98_yr'], 
                             min_dist[i], 
                             max_dist[i])
        
      #calculate GPP with fPAR                      
      merged_dataframe.loc[merged_dataframe['pft']==pft, 'STARFM_GPP'] = rctmu.calcGPP(params['RCTM_params'][pft]['starfm'], 
                          df_pft['TS'], 
                          df_pft['SWC2'], 
                          df_pft['VPD'], 
                          df_pft['SW_IN_NLDAS'], 
                          merged_dataframe.loc[merged_dataframe['pft']==pft,'STARFM_fPAR'])
      
      #comparison statistics
      pft_rmses.append(mean_squared_error(merged_dataframe.loc[merged_dataframe['pft']==pft,'GPP_measured'], merged_dataframe.loc[merged_dataframe['pft']==pft,'STARFM_GPP'], squared = False))
      pft_r2s.append(r2_score(merged_dataframe.loc[merged_dataframe['pft']==pft,'GPP_measured'], merged_dataframe.loc[merged_dataframe['pft']==pft,'STARFM_GPP']))
      lens.append(len(df_pft))
      
    pft_rmses = np.array(pft_rmses)  
    pft_r2s = np.array(pft_r2s)
    lens = np.array(pft_rmses)  
    #comparison statistics
    r2=np.sum(pft_r2s*lens)/np.sum(lens)
    #r2 = r2_score(merged_dataframe['GPP_measured'], merged_dataframe['STARFM_GPP'])
    #rmse = mean_squared_error(merged_dataframe['GPP_measured'], merged_dataframe['STARFM_GPP'], squared = False)
    rmse=np.sum(pft_rmses*lens)/np.sum(lens)
    rmses.append(rmse)
    r2s.append(r2)
    
  #find best run based on r2 or rmse
  minimum_rmse_index= np.argmin(rmses)
  maximum_r2_index= np.argmax(r2s)
    
  best_min_val = min_dist[minimum_rmse_index]
  best_max_val = max_dist[minimum_rmse_index]
  min_rmse = rmses[minimum_rmse_index]

  best_min_val_r2 = min_dist[maximum_r2_index]
  best_max_val_r2 = max_dist[maximum_r2_index]
  max_r2 = r2s[maximum_r2_index]
  cluster=0
  print('Cluster: {}'.format(cluster))
  print('best fPAR min value: {}\nbest fPAR max value: {}\nminimum RMSE: {}'.format(best_min_val, best_max_val, min_rmse))
  print('best fPAR min value: {}\nbest fPAR max value: {}\nmaximum r2: {}'.format(best_min_val_r2, best_max_val_r2, max_r2))
  print('=======================================================================================')
    
  #append results to dictionary
  results['cluster'].append(0)
  results['fPAR_min'].append(np.nanmean([best_min_val,best_min_val_r2]))
  results['fPAR_max'].append(np.nanmean([best_max_val,best_max_val_r2]))
  results['rmse'].append(min_rmse)
  results['r2'].append(max_r2)
    
  return results

def get_fPAR_ref(merged_dataframe, ref_pars, iterations=2000, seed=2023):
  
  np.random.seed(seed)
  
  #fPAR min-max distributions to select from
  min_dist = np.random.uniform(-0.3, 0.5, iterations) 
  max_dist = np.random.uniform(0.05, 1.3, iterations)
  
  #if there is a min value > max value, swap their values in the sample arrays
  min_swap_vals = min_dist[min_dist>max_dist]
  max_swap_vals = max_dist[max_dist<min_dist]
  min_dist[min_dist>max_dist] = max_swap_vals
  max_dist[max_dist==min_dist] = min_swap_vals
  
  #dictionary to store results
  cluster_results = {'cluster':[], 'fPAR_min': [], 'fPAR_max': [], 'rmse':[], 'r2': []}

  for cluster in pd.unique(merged_dataframe['cluster'].sort_values()):
    df_merged=merged_dataframe[(merged_dataframe['cluster']==cluster)]
    rmses=[]
    r2s=[]

    for i in range(0, len(min_dist)):
      
      #calculate fPAR with sampled fPAR min and max  
      STARFM_fPAR=rctmu.calc_fPAR(
                            df_merged['ndvi'], 
                            df_merged['ndvi_02_yr'], 
                            df_merged['ndvi_98_yr'],  
                            df_merged['ndvi_sr'], 
                            df_merged['ndvi_sr_02_yr'], 
                            df_merged['ndvi_sr_98_yr'], 
                            min_dist[i], 
                            max_dist[i])
      
      #calculate GPP with fPAR                      
      STARFM_GPP = rctmu.calcGPP(ref_pars, 
                        df_merged['TS'], 
                        df_merged['SWC2'], 
                        df_merged['VPD'], 
                        df_merged['SW_IN_NLDAS'], 
                        STARFM_fPAR)
      
      #comparison statistics
      r2 = r2_score(df_merged['GPP_measured'], STARFM_GPP)
      rmse = mean_squared_error(df_merged['GPP_measured'], STARFM_GPP, squared = False)
      rmses.append(rmse)
      r2s.append(r2)
    
    #find best run based on r2 or rmse
    minimum_rmse_index= np.argmin(rmses)
    maximum_r2_index= np.argmax(r2s)
    
    best_min_val = min_dist[minimum_rmse_index]
    best_max_val = max_dist[minimum_rmse_index]
    min_rmse = rmses[minimum_rmse_index]

    best_min_val_r2 = min_dist[maximum_r2_index]
    best_max_val_r2 = max_dist[maximum_r2_index]
    max_r2 = r2s[maximum_r2_index]
    
    print('Cluster: {}'.format(cluster))
    print('best fPAR min value: {}\nbest fPAR max value: {}\nminimum RMSE: {}'.format(best_min_val, best_max_val, min_rmse))
    print('best fPAR min value: {}\nbest fPAR max value: {}\nmaximum r2: {}'.format(best_min_val_r2, best_max_val_r2, max_r2))
    print('=======================================================================================')
    
    #append results to dictionary
    cluster_results['cluster'].append(cluster)
    cluster_results['fPAR_min'].append(np.nanmean([best_min_val,best_min_val_r2]))
    cluster_results['fPAR_max'].append(np.nanmean([best_max_val,best_max_val_r2]))
    cluster_results['rmse'].append(min_rmse)
    cluster_results['r2'].append(max_r2)
    
  return cluster_results


######################## run optimization ###############################   
print(STARFM_df['site'].unique())
method='NDVI'
pfts=['grassland', 'grass-shrub', 'grass-tree']

STARFM_df['cluster']=0

STARFM_df = STARFM_df[~STARFM_df['ndvi'].isna()]
STARFM_df = STARFM_df[~STARFM_df['GPP_measured'].isna()]

#opt_results = get_global_fPAR_ref(STARFM_df, params)

for pft in pfts:

  sfm_params=params['RCTM_params'][pft]['starfm']
  
  print('optimizing')
  opt_results = get_fPAR_ref(STARFM_df[STARFM_df['pft']==pft], sfm_params,  iterations=iterations)
  df_opt_results = pd.DataFrame(opt_results)
  df_opt_results.to_csv(f'output/fPAR_optimization/data/ref_values/STARFM_fPAR_ref_values_{pft}.csv', index=False)
  
  df_STARFM_opt = pd.merge(STARFM_df[STARFM_df['pft']==pft], df_opt_results, on = ['cluster'],how='inner')
  
  print(df_STARFM_opt['fPAR_min'].unique())
  print(df_STARFM_opt['fPAR_max'].unique())

  df_STARFM_opt['STARFM_fPAR']=rctmu.calc_fPAR(df_STARFM_opt['ndvi'], 
                                              df_STARFM_opt['ndvi_02_yr'], 
                                              df_STARFM_opt['ndvi_98_yr'], 
                                              df_STARFM_opt['ndvi_sr'], 
                                              df_STARFM_opt['ndvi_sr_02_yr'], 
                                              df_STARFM_opt['ndvi_sr_98_yr'], 
                                              df_STARFM_opt['fPAR_min'], 
                                              df_STARFM_opt['fPAR_max'])
  
  df_STARFM_opt['STARFM_GPP'] = rctmu.calcGPP(sfm_params,
                          df_STARFM_opt['TS'], 
                          df_STARFM_opt['SWC2'], 
                          df_STARFM_opt['VPD'], 
                          df_STARFM_opt['SW_IN_NLDAS'], 
                          df_STARFM_opt['STARFM_fPAR'])
  if pft=='grassland':
    grass = df_STARFM_opt.loc[~df_STARFM_opt['site'].isin(params['Ameriflux_site_pfts']['hay-pasture'])]
    hay_pasture = df_STARFM_opt.loc[df_STARFM_opt['site'].isin(params['Ameriflux_site_pfts']['hay-pasture'])]
    grass.to_csv(f'output/fPAR_optimization/data/GPP/STARFM_fPAR_optimization_{method}_GPP_grassland.csv')
    hay_pasture.to_csv(f'output/fPAR_optimization/data/GPP/STARFM_fPAR_optimization_{method}_GPP_hay-pasture.csv')
    
  else:
    df_STARFM_opt.to_csv(f'output/fPAR_optimization/data/GPP/STARFM_fPAR_optimization_{method}_GPP_{pft}.csv')
  print(df_STARFM_opt['pft'].unique())
  STARFM_r2=r2_score(df_STARFM_opt['GPP_measured'], df_STARFM_opt['STARFM_GPP']).round(2)
  
  STARFM_rmse = mean_squared_error(df_STARFM_opt['GPP_measured'], df_STARFM_opt['STARFM_GPP'], squared = False).round(2)
  
  STARFM_bias = rctmu.calc_bias(df_STARFM_opt['GPP_measured'], df_STARFM_opt['STARFM_GPP']).round(2)
  
  ###### Plot Results ######
  for site in df_STARFM_opt['site'].unique():
    fig, axes = plt.subplots(figsize = (8, 5))
    sns.scatterplot(data=df_STARFM_opt[df_STARFM_opt['site']==site], x='time', y='STARFM_GPP', label='STARFM', color = 'black')
    sns.lineplot(data=df_STARFM_opt[df_STARFM_opt['site']==site], x='time', y='GPP_measured', label='Measured')
    fig.tight_layout()
  
    plt.savefig('output/fPAR_optimization/plots/GPP_comparisons/GPP_comp_{}.jpg'.format(site), dpi=300)
    plt.show()
  
  print('STARFM r2: {}'.format(STARFM_r2))
  
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
  
  plt.savefig('output/fPAR_optimization/plots/GPP_comparisons/GPP_{}_{}.jpg'.format(method, pft), dpi=300)

