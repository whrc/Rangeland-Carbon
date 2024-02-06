from google.cloud import storage
import rioxarray as rxr
import os
import numpy as np
import pandas as pd
import multiprocessing
from dask import array as da
import xarray as xr
import seaborn as sns
from matplotlib import pyplot as plt

storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')

bucket_name='rangelands'
path_to_footprint_list = '/home/amullen/Rangeland-Carbon/res/site_footprints/HLR_tiles_subregion.txt'
bucket = storage_client.get_bucket(bucket_name)

path_to_landcover = 'gs://' + bucket.name + '/Ranch_Runs/HLR/landcover/NLCD_HLR_2019.tif'

#open landcover
landcover = rxr.open_rasterio(path_to_landcover)

datelist = pd.date_range('2001-01-01', '2023-01-01').strftime('%Y-%m-%d').tolist()

pft_dict = {'grass': [], 'shrub': [], 'tree': []}
index_dict = dict(zip(datelist, [pft_dict]*len(datelist)))

filepaths=[]


# for each tile
with open(path_to_footprint_list) as f:
  
  sites = [line.rstrip('\n').split('/')[-1] for line in f]
  
  site_dfs=[]
  
  for site in sites:
  
    print(site)
    ndvi_file = 'gs://' + bucket.name + f'/Ranch_Runs/HLR/{site}/covariates_nc/{site}_ndvi.nc'
    filepaths.append(ndvi_file)
    
    ndvi = rxr.open_rasterio(ndvi_file, masked=True)
    site_dates = ndvi.indexes['time'].to_datetimeindex()
    
    # reproject landcover
    print('reprojecting landcover')
    landcover_res = landcover.rio.reproject_match(ndvi)

    print('masking grass')
    grass_mask = ndvi.where((landcover_res[0]==71) | (landcover_res[0]==81) | (landcover_res[0]==82))
    print('averaging')
    grass_average = grass_mask.mean(dim=['x', 'y'])
    
    grass_df = pd.DataFrame({'site': [site]*len(site_dates), 
                             'pft': ['grass']*len(site_dates), 
                             'Date': site_dates, 
                             'NDVI': grass_average.values})
    site_dfs.append(grass_df)
    
    print('masking shrubs')
    shrub_mask = ndvi.where(landcover_res[0]==52)
    print('averaging')
    shrub_average = shrub_mask.mean(dim=['x', 'y'])
    shrub_df = pd.DataFrame({'site': [site]*len(site_dates), 
                             'pft': ['shrub']*len(site_dates), 
                             'Date': site_dates, 
                             'NDVI': shrub_average.values})
    site_dfs.append(shrub_df)
    
    
    print('masking trees')
    tree_mask = ndvi.where((landcover_res[0]==41) | (landcover_res[0]==42) | (landcover_res[0]==43))
    print('averaging')
    tree_average = tree_mask.mean(dim=['x', 'y'])
    tree_df = pd.DataFrame({'site': [site]*len(site_dates), 
                             'pft': ['tree']*len(site_dates), 
                             'Date': site_dates, 
                             'NDVI': tree_average.values})
    site_dfs.append(tree_df)
    
    
    
    #fig, ax = plt.subplots(figsize = (15, 5))
    
    #sns.lineplot(x=site_dates, y=grass_average.values, label='grass')
    #sns.lineplot(x=site_dates, y=shrub_average.values, label = 'shrub')
    #sns.lineplot(x=site_dates, y=tree_average.values, label = 'tree')
    #plt.xticks(rotation=45)
    
    #plt.savefig(f'output/HLR/HLR_ndvi_{site}.jpg')
    
  site_dfs=pd.concat(site_dfs, ignore_index=True)
  site_dfs.to_csv('output/HLR/HLR_ndvi.csv')
  
  fig, ax = plt.subplots(figsize = (20, 4))
  
  sns.lineplot(data = site_dfs, x='Date', y='NDVI', hue='pft')
  plt.xticks(rotation=45)
  fig.tight_layout()
  plt.savefig(f'output/HLR/HLR_ndvi.jpg', dpi=300)
  
    
    #for i, date in enumerate(site_dates):
    #  index_dict[date]['grass'].extend(ndvi_grass[i])
    #  index_dict[date]['shrub'].extend(ndvi_shrub[i])
    #  index_dict[date]['tree'].extend(ndvi_tree[i])
    #break
    
#average
#grass_means = []
#shrub_means = []
#tree_means = []

#for k in datelist:
#    print(k)
#    grass_mean = da.nanmean(da.from_array(index_dict[k]['grass'], chunks='auto')).compute()
#    shrub_mean = da.nanmean(da.from_array(index_dict[k]['shrub'], chunks='auto')).compute()
#    tree_mean = da.nanmean(da.from_array(index_dict[k]['tree'], chunks='auto')).compute()
    
    #grass_mean = distribute_mean_calc(index_dict[k]['grass'])
    #shrub_mean = distribute_mean_calc(index_dict[k]['shrub'])
    #tree_mean = distribute_mean_calc(index_dict[k]['tree'])
    
#    grass_means.append(grass_mean)
#    shrub_means.append(shrub_mean)
#    tree_means.append(tree_mean)
    
#df = pd.DataFrame({'Date': datelist, 'grass_ndvi': grass_means, 'shrub_ndvi': shrub_means, 'tree_ndvi': tree_means})
#df.to_csv('output/HLR/HLR_ndvi_summary.csv')