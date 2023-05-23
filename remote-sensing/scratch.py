import rioxarray as rxr
from google.cloud import storage
import rioxarray as rxr
import os
import utils
import pandas as pd
import numpy as np
from itertools import chain
from scipy.stats import t
import xarray as xr
from matplotlib import pyplot as plt
import seaborn as sns
import preprocess_starfm_imagery as psi


storage_client = storage.Client.from_service_account_json('gee_key.json')

bucket='rangelands'
modis_dir = 'Ameriflux_sites/Rws_starfm/modis_test_smooth/'
modis_smooth_dir = 'Ameriflux_sites/Rws_starfm/modis_smooth/'
landsat_dir = 'Ameriflux_sites/Rws_starfm/landsat_scenes/'
landsat_new_dir = 'Ameriflux_sites/Rws_starfm/landsat_test/'

bucket = storage_client.get_bucket(bucket)

modis_df = psi.get_modis_date_df(modis_dir, bucket)

modis_smooth_df = psi.get_modis_date_df(modis_smooth_dir, bucket)
modis_smooth_df['year'] = pd.DatetimeIndex(modis_smooth_df['im_date']).year
modis_smooth_df=modis_smooth_df[modis_smooth_df['year']==2002]

landsat_df = psi.get_landsat_date_df(landsat_dir, bucket)
landsat_df['year'] = pd.DatetimeIndex(landsat_df['im_date']).year
landsat_df=landsat_df[landsat_df['year']==2002]

landsat_new_df = psi.get_landsat_date_df(landsat_new_dir, bucket)
landsat_new_df['year'] = pd.DatetimeIndex(landsat_new_df['im_date']).year
landsat_new_df=landsat_new_df[landsat_new_df['year']==2002]

modis_df['evi'] = np.nan
modis_smooth_df['evi'] = np.nan
landsat_df['evi'] = np.nan
landsat_new_df['evi'] = np.nan

modis_df['ndvi'] = np.nan
modis_smooth_df['ndvi'] = np.nan
landsat_df['ndvi'] = np.nan
landsat_new_df['ndvi'] = np.nan

print(len(landsat_df))
print(len(landsat_new_df))

print('MODIS')
for index, row in modis_df.iterrows():
  im=rxr.open_rasterio('gs://' + bucket.name + '/' + row['im_path'])
  im=im.astype(np.float32)
  im.values=im.values/10000
  im.values[im.values==0] = np.nan

  evi = np.nanmean(2.5 * ((im.values[3]-im.values[2])/(im.values[3] + 6*im.values[2] - 7.5*im.values[0] + 1)))
  ndvi = np.nanmean((im.values[3]-im.values[2])/(im.values[3]+im.values[2]))
  
  modis_df.loc[index, 'evi']=evi
  modis_df.loc[index, 'ndvi']=ndvi

print('MODIS smoothed')
for index, row in modis_smooth_df.iterrows():
  im=rxr.open_rasterio('gs://' + bucket.name + '/' + row['im_path'])
  im=im.astype(np.float32)
  im.values=im.values/10000
  im.values[im.values==0] = np.nan
  
  #evi = np.nanmean(2.5 * ((im.values[3]-im.values[2])/(im.values[3] + 6*im.values[2] - 7.5*im.values[0] + 1)))
  #ndvi = np.nanmean((im.values[3]-im.values[2])/(im.values[3]+im.values[2]))
  
  evi = np.nanmean(2.5 * ((im.values[1]-im.values[0])/(im.values[1] + 6*im.values[0] - 7.5*im.values[2] + 1)))
  ndvi = np.nanmean((im.values[1]-im.values[0])/(im.values[1]+im.values[0]))
  
  modis_smooth_df.loc[index, 'evi']=evi
  modis_smooth_df.loc[index, 'ndvi']=ndvi

print('Landsat')
for index, row in landsat_df.iterrows():
  im=rxr.open_rasterio('gs://' + bucket.name + '/' + row['im_path'])
  im=im.astype(np.float32)
  im.values=im.values/10000
  im.values[im.values==0] = np.nan

  evi = np.nanmean(2.5 * ((im.values[3]-im.values[2])/(im.values[3] + 6*im.values[2] - 7.5*im.values[0] + 1)))
  ndvi = np.nanmean((im.values[3]-im.values[2])/(im.values[3]+im.values[2]))
  
  landsat_df.loc[index, 'evi']=evi
  landsat_df.loc[index, 'ndvi']=ndvi
  
print('landsat new')
for index, row in landsat_new_df.iterrows():
  im=rxr.open_rasterio('gs://' + bucket.name + '/' + row['im_path'])
  im=im.astype(np.float32)
  im.values=im.values/10000
  im.values[im.values==0] = np.nan

  evi = np.nanmean(2.5 * ((im.values[3]-im.values[2])/(im.values[3] + 6*im.values[2] - 7.5*im.values[0] + 1)))
  ndvi = np.nanmean((im.values[3]-im.values[2])/(im.values[3]+im.values[2]))
  
  landsat_new_df.loc[index, 'evi']=evi
  landsat_new_df.loc[index, 'ndvi']=ndvi

fig, ax = plt.subplots()
sns.lineplot(data=modis_df, x='im_date', y='evi', label='new modis smooth')
sns.lineplot(data=modis_smooth_df, x='im_date', y='evi', label='original modis smooth')
sns.lineplot(data=landsat_df, x='im_date', y='evi', label='original landsat')
sns.lineplot(data=landsat_new_df, x='im_date', y='evi', label='new landsat')

plt.ylim(0,0.2)
plt.show()
plt.savefig('sfm_in_imagery_evi_comp_c2.jpg', dpi=300)

fig, ax = plt.subplots()
sns.lineplot(data=modis_df, x='im_date', y='ndvi', label='new modis smooth')
sns.lineplot(data=modis_smooth_df, x='im_date', y='ndvi', label='original modis smooth')
sns.lineplot(data=landsat_df, x='im_date', y='ndvi', label='original landsat')
sns.lineplot(data=landsat_new_df, x='im_date', y='ndvi', label='new landsat')
plt.ylim(0,0.4)
plt.show()
plt.savefig('sfm_in_imagery_ndvi_comp_c2.jpg', dpi=300)
  

