import numpy as np
from google.cloud import storage
import sys
sys.path.insert(1, '../modeling')
sys.path.insert(1, '../utils')
import utils
import os
import pandas as pd

def write_list_to_file(lst, file_path):
  """
  Write each list from the input list of lists to a new line in a text file.

  Parameters:
  - lists (list): A list of lists containing numerical values.
  - file_path (str): The path to the text file.

  Returns:
  - None
  """
  if os.path.isfile(file_path):
    os.remove(file_path)
  with open(file_path, 'w') as file:
    for i, item in enumerate(lst):
      if i<len(lst)-1:
        file.write(item + '\n')
      else:
        file.write(item)
        
        
bucket_name='rangelands'

storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
bucket = storage_client.get_bucket(bucket_name)
path_to_temp = '/home/amullen/temp/'

paths_to_merge={'GPP_path': [], 'Reco_path': [], 'GPP_change_path': [], 'GPP_change_signif_path': [], 'SOC_change_path': [], 'landcover_path': [], 'landcover_path': [], 'flux_csv': [], 'C_stock_csv': []}
C_stock_dfs = []
flux_dfs = []

roi_file='/home/amullen/Rangeland-Carbon/res/site_footprints/HLD_tiles.txt'

with open(roi_file) as f:
  paths = [line.rstrip('\n') for line in f]
  
  for path in paths:
    
    ranch='/'.join(path.split('/')[-2:])
    print(ranch)
    if ranch in ['...']:
      continue
   
    paths_to_merge['GPP_path'].append(f'/vsigs/rangelands/Ranch_Runs/{ranch}/analysis/mean_annual_GPP.tif')
    paths_to_merge['Reco_path'].append(f'/vsigs/rangelands/Ranch_Runs/{ranch}/analysis/mean_annual_Reco.tif')
    paths_to_merge['GPP_change_path'].append(f'/vsigs/rangelands/Ranch_Runs/{ranch}/analysis/GPP_change_analysis.tif')
    paths_to_merge['GPP_change_signif_path'].append(f'/vsigs/rangelands/Ranch_Runs/{ranch}/analysis/GPP_change_analysis_signif.tif')
    paths_to_merge['SOC_change_path'].append(f'/vsigs/rangelands/Ranch_Runs/{ranch}/analysis/C_stock_average_change_2021-2003.tif')
    paths_to_merge['landcover_path'].append(f'/vsigs/rangelands/Ranch_Runs/{ranch}/landcover/fused_landcover.tif')
    
    flux_df = pd.read_csv(f'gs://rangelands/Ranch_Runs/{ranch}/analysis/flux_timeseries.csv')
    flux_df['ranch'] = ranch
    
    C_stock_df = pd.read_csv(f'gs://rangelands/Ranch_Runs/{ranch}/analysis/C_stock_timeseries.csv')
    C_stock_df['ranch'] = ranch
    
    flux_dfs.append(flux_df)
    C_stock_dfs.append(C_stock_df)
    
C_stock_dfs = pd.concat(C_stock_dfs)
flux_dfs = pd.concat(flux_dfs)
print(C_stock_dfs.columns)

print(flux_dfs.columns)    
C_stock_dfs = C_stock_dfs.groupby(by = ['Date', 'value']).agg({'TSOC_mean_gC/m2': 'mean',
                                                                'TSOC_sum_gC': 'sum',
                                                                'TSOC_std_gC': 'sum'
                                                                })
                                                                
flux_dfs = flux_dfs.groupby(by = ['Date', 'value']).agg({'GPP_mean_gC/m2': 'mean',
                                                          'GPP_sum_gC': 'sum',
                                                          'GPP_std_gC': 'sum',
                                                          'Reco_mean_gC/m2': 'mean',
                                                          'Reco_sum_gC': 'sum',
                                                          'Reco_std_gC': 'sum',
                                                          'NEE_mean_gC/m2': 'mean',
                                                          'NEE_sum_gC': 'sum',
                                                          'NEE_std_gC': 'sum'
                                                        })
flux_dfs.to_csv('flux_timeseries.csv')
C_stock_dfs.to_csv('C_stock_timeseries.csv')                                                                        

write_list_to_file(paths_to_merge['GPP_path'], 'GPP_file.txt')
write_list_to_file(paths_to_merge['Reco_path'], 'Reco_file.txt')
write_list_to_file(paths_to_merge['GPP_change_path'], 'GPP_change_file.txt')
write_list_to_file(paths_to_merge['GPP_change_signif_path'], 'GPP_change_signif_file.txt')
write_list_to_file(paths_to_merge['SOC_change_path'], 'SOC_change_file.txt')
write_list_to_file(paths_to_merge['landcover_path'], 'landcover_file.txt')

#gdal.merge('/vsigs/rangelands/Ranch_Runs/HLD/merged_results/mean_annual_GPP.tif', optfile='GPP_file.txt')
#gdal_merge.py -o 'mean_annual_GPP.tif' -n 0 --optfile GPP_file.txt
#gdal_merge.py -o 'mean_annual_Reco.tif' -n 0 --optfile Reco_file.txt
#gdal_merge.py -o 'GPP_change_analysis.tif' -n 0 --optfile GPP_change_file.txt
#gdal_merge.py -o 'GPP_change_analysis_signif.tif' -n 0 --optfile GPP_change_signif_file.txt
#gdal_merge.py -o 'C_stock_average_change_2021-2003.tif' -n 0 --optfile SOC_change_file.txt
#gdal_merge.py -o 'fused_landcover.tif' --optfile landcover_file.txt
