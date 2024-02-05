from google.cloud import storage
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import sys
sys.path.insert(1, '../utils')
import utils
import os
import geopandas as gpd
import gcsfs
import numpy as np

bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
bucket = storage_client.get_bucket(bucket_name)
path_to_temp = '/home/amullen/temp/'

ranch='HLD'

#path_to_site_roi = f'gs://rangelands/Ranch_Runs/{ranch}/site_roi/Harts_Basin.shp'
#path_to_site_roi = f'gs://rangelands/Ranch_Runs/{ranch}/site_roi/Hanks_Valley.shp'
path_to_site_roi = f'gs://rangelands/Ranch_Runs/{ranch}/site_roi/MCK.shp'

#site_roi = gpd.read_file(path_to_site_roi)

#site_roi = site_roi.to_crs('EPSG:5070')
#area=site_roi.area
#area_ha = area / 10000


C_stock_hist = pd.read_csv(f'gs://rangelands/Ranch_Runs/{ranch}/analysis/C_stock_timeseries.csv', parse_dates=['Date'])
flux_hist = pd.read_csv(f'gs://rangelands/Ranch_Runs/{ranch}/analysis/flux_timeseries.csv', parse_dates=['Date'])

hist = pd.merge(flux_hist, C_stock_hist, on=['Date', 'value'])
hist['month']=hist['Date'].dt.month
hist['year'] = hist['Date'].dt.year

hist = hist.loc[hist['year']<2023]
#hist = hist.loc[hist['year']<2022]

hist['class'] = ''
hist.loc[hist['value']==1, 'class'] = 'grassland'
hist.loc[hist['value']==2, 'class'] = 'grass-shrub'
hist.loc[hist['value']==3, 'class'] = 'grass-tree'

hist = hist.rename(columns = {'GPP_mean_gC/m2': 'GPP_mean_gC/m2/d', 'NEE_mean_gC/m2': 'NEE_mean_gC/m2/d', 'Reco_mean_gC/m2': 'Reco_mean_gC/m2/d',
                              'GPP_sum_gC': 'GPP_sum_gC/d', 'NEE_sum_gC': 'NEE_sum_gC/d', 'Reco_sum_gC': 'Reco_sum_gC/d',
                              'GPP_std_gC': 'GPP_std_gC/d', 'NEE_std_gC': 'NEE_std_gC/d', 'Reco_std_gC': 'Reco_std_gC/d'})

#convert gC/m2 to tC/ha, gC to tC
hist['GPP_mean_tC/ha/d'] = hist['GPP_mean_gC/m2/d'] * 0.01
hist['NEE_mean_tC/ha/d'] = hist['NEE_mean_gC/m2/d'] * 0.01
hist['Reco_mean_tC/ha/d'] = hist['Reco_mean_gC/m2/d'] * 0.01
hist['GPP_sum_tC/d'] = hist['GPP_sum_gC/d'] / 1000000
hist['NEE_sum_tC/d'] = hist['NEE_sum_gC/d'] / 1000000
hist['Reco_sum_tC/d'] = hist['Reco_sum_gC/d'] / 1000000
hist['GPP_std_tC/d'] = hist['GPP_std_gC/d'] / 1000000
hist['NEE_std_tC/d'] = hist['NEE_std_gC/d'] / 1000000
hist['Reco_std_tC/d'] = hist['Reco_std_gC/d'] / 1000000

hist['TSOC_mean_tC/ha'] = hist['TSOC_mean_gC/m2'] * 0.01
hist['TSOC_sum_tC'] = hist['TSOC_sum_gC'] / 1000000
hist['TSOC_std_tC'] = hist['TSOC_std_gC'] / 1000000

hist['TSOC_sum_MtC'] = hist['TSOC_sum_tC'] / 1000000
hist['TSOC_std_MtC'] = hist['TSOC_std_tC'] / 1000000

pal = sns.color_palette(['#424235', '#7E7C52', '#CFA67C'])

#plotting daily flux time series
fig, ax = plt.subplots(figsize = (8, 2.25))

sns.lineplot(data = hist, x='Date', y='GPP_mean_gC/m2/d', color=pal[0],label='GPP', ax=ax)
sns.lineplot(data = hist, x='Date', y='NEE_mean_gC/m2/d', color=pal[1],label='NEE', ax=ax)
sns.lineplot(data = hist, x='Date', y='Reco_mean_gC/m2/d', color=pal[2],label='Reco', ax=ax)

ax.set_ylabel('Flux\n(gC/$m^2$/day)', fontsize=14)
ax.set_xlabel('')

ax.legend(title = '', frameon=False, fontsize=10, ncol=3)

fig.tight_layout()
plt.savefig(os.path.join(path_to_temp, 'temp_im_write.png'), dpi=300, transparent=True)
utils.gs_write_blob(os.path.join(path_to_temp, 'temp_im_write.png'), f'Ranch_Runs/{ranch}/analysis/5d_timeseries_summary.png', bucket)
plt.show()

#plotting yearly total C stocks
yrly_C_stocks = hist[['year', 'class', 'TSOC_sum_tC', 'TSOC_std_tC', 'TSOC_sum_MtC', 'TSOC_std_MtC']].groupby(by=['year', 'class']).mean().reset_index()
yrly_C_stocks['year'] = pd.to_datetime(yrly_C_stocks['year'], format='%Y')



#yrly_C_stocks['TSOC Mt C'] = (yrly_C_stocks['TSOC (tC/ha)'] * area_ha) / 1000000
def plot_sd_as_error(vector):
  return(-1*vector, vector)
  
fig, axes = plt.subplots(3,1,figsize = (3.8, 4.2), sharex=True)
#fig, axes = plt.subplots(2,1,figsize = (3.8, 3.2), sharex=True)

sns.lineplot(data = yrly_C_stocks.loc[yrly_C_stocks['class']=='grassland'], x='year', y='TSOC_sum_MtC', color=pal[0], label=None, ax=axes[0])
#axes[0].fill_between(x=yrly_C_stocks.loc[yrly_C_stocks['class']=='grassland', 'year'], 
#                    y1=yrly_C_stocks.loc[yrly_C_stocks['class']=='grassland', 'TSOC_sum_MtC'] - yrly_C_stocks.loc[yrly_C_stocks['class']=='grassland', 'TSOC_std_tC'],
#                    y2=yrly_C_stocks.loc[yrly_C_stocks['class']=='grassland', 'TSOC_sum_MtC'] + yrly_C_stocks.loc[yrly_C_stocks['class']=='grassland', 'TSOC_std_tC'])

sns.lineplot(data = yrly_C_stocks.loc[yrly_C_stocks['class']=='grass-shrub'], x='year', y='TSOC_sum_MtC', color=pal[0], label=None, ax=axes[1])

sns.lineplot(data = yrly_C_stocks.loc[yrly_C_stocks['class']=='grass-tree'], x='year', y='TSOC_sum_MtC', color=pal[0], label=None, ax=axes[2])

axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[2].spines['top'].set_visible(False)
axes[2].spines['right'].set_visible(False)

axes[0].set_xlabel('')
axes[0].yaxis.set_label_position("right")
axes[0].set_ylabel('grass', fontsize=12, fontweight='bold')

axes[1].set_xlabel('')
axes[1].yaxis.set_label_position("right")
axes[1].set_ylabel('grass-shrub', fontsize=12, fontweight='bold')

axes[2].set_xlabel('')
axes[2].yaxis.set_label_position("right")
axes[2].set_ylabel('grass-tree', fontsize=12, fontweight='bold')

axes[0].yaxis.set_tick_params(labelsize=14)
axes[1].yaxis.set_tick_params(labelsize=14)
axes[2].yaxis.set_tick_params(labelsize=14)

fig.supylabel('SOC (MtC)', fontsize=14)

plt.xticks(rotation=90, fontsize=14)

fig.tight_layout()
plt.savefig(os.path.join(path_to_temp, 'temp_im_write.png'), dpi=300, transparent=True)
utils.gs_write_blob(os.path.join(path_to_temp, 'temp_im_write.png'), f'Ranch_Runs/{ranch}/analysis/soil_carbon_summary.png', bucket)

#plotting yearly cumulative fluxes
pal = sns.color_palette(['#7E7C52', '#424235', '#CFA67C'])
#pal = sns.color_palette(['#7E7C52', '#CFA67C'])
hist[['GPP_mean_gC/m2/d', 'GPP_sum_tC/d', 'GPP_std_tC/d', 'NEE_mean_gC/m2/d', 'NEE_sum_tC/d', 'Reco_mean_gC/m2/d', 'Reco_sum_tC/d', 'Reco_std_tC/d', 'NEE_std_tC/d']] = 5 * hist[['GPP_mean_gC/m2/d', 'GPP_sum_tC/d', 'GPP_std_tC/d', 'NEE_mean_gC/m2/d', 'NEE_sum_tC/d', 'Reco_mean_gC/m2/d', 'Reco_sum_tC/d', 'Reco_std_tC/d', 'NEE_std_tC/d']]

yrly = hist[['year', 'class', 'GPP_mean_gC/m2/d', 'GPP_sum_tC/d', 'GPP_std_tC/d', 'NEE_mean_gC/m2/d', 'NEE_sum_tC/d', 'Reco_mean_gC/m2/d', 'Reco_sum_tC/d', 'Reco_std_tC/d', 'NEE_std_tC/d']].groupby(by=['year', 'class']).sum().reset_index()

#yrly = yrly.loc[yrly['class']!='grass-tree']

yrly['year'] = pd.to_datetime(yrly['year'], format='%Y')

fig, axes=plt.subplots(3,1,figsize = (6, 4.8), sharex=True)

sns.lineplot(data = yrly, x='year', y='GPP_mean_gC/m2/d', hue='class', ax=axes[0], palette=pal)
sns.lineplot(data = yrly, x='year', y='Reco_mean_gC/m2/d', hue='class', ax=axes[1], palette=pal)
sns.lineplot(data = yrly, x='year', y='NEE_mean_gC/m2/d', hue='class', ax=axes[2], palette=pal)

#axes[2].fill_between(x=yrly.loc[yrly['class']=='grassland', 'year'], 
#                    y1=yrly.loc[yrly['class']=='grassland', 'NEE_sum_tC/d'] - yrly.loc[yrly['class']=='grassland', 'NEE_std_tC/d'],
#                    y2=yrly.loc[yrly['class']=='grassland', 'NEE_sum_tC/d'] + yrly.loc[yrly['class']=='grassland', 'NEE_std_tC/d'])
                    
axes[2].axhline(y=0, xmin=0, xmax=1, c="grey", linewidth=0.5, linestyle='--', alpha=0.6)

#despine axes
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[2].spines['top'].set_visible(False)
axes[2].spines['right'].set_visible(False)

#set x and y labels
axes[0].set_xlabel('')
axes[0].set_ylabel('GPP\n(gC/m2/yr)', fontsize=14)
axes[1].set_xlabel('')
axes[1].set_ylabel('Reco\n(gC/m2/yr)', fontsize=14)
axes[2].set_xlabel('')
axes[2].set_ylabel('NEE\n(gC/m2/yr)', fontsize=14)

#legend props
axes[2].legend(title = '', frameon=True, fontsize=10, ncol=3, loc='lower center')
axes[1].get_legend().remove()
axes[0].get_legend().remove()

#ticks
plt.xticks(rotation=45, fontsize=14)
axes[0].tick_params(axis='y', labelsize=12)
axes[1].tick_params(axis='y', labelsize=12)
axes[2].tick_params(axis='y', labelsize=12)

fig.tight_layout()
plt.savefig(os.path.join(path_to_temp, 'temp_im_write.png'), dpi=300, transparent=True)
utils.gs_write_blob(os.path.join(path_to_temp, 'temp_im_write.png'), f'Ranch_Runs/{ranch}/analysis/yearly_fluxes.png', bucket)
plt.show()

yrly = pd.merge(yrly, yrly_C_stocks, on=['year', 'class'])
yrly[['GPP_mean_gC/m2/d', 'GPP_sum_tC/d', 'GPP_std_tC/d', 'NEE_mean_gC/m2/d', 'NEE_sum_tC/d', 'Reco_mean_gC/m2/d', 'Reco_sum_tC/d', 'Reco_std_tC/d', 'NEE_std_tC/d', 'TSOC_sum_tC', 'TSOC_std_tC']] = np.round(yrly[['GPP_mean_gC/m2/d', 'GPP_sum_tC/d', 'GPP_std_tC/d', 'NEE_mean_gC/m2/d', 'NEE_sum_tC/d', 'Reco_mean_gC/m2/d', 'Reco_sum_tC/d', 'Reco_std_tC/d', 'NEE_std_tC/d', 'TSOC_sum_tC', 'TSOC_std_tC']])

yrly = yrly[['year', 'class', 'GPP_sum_tC/d', 'Reco_sum_tC/d', 'NEE_sum_tC/d', 'TSOC_sum_tC']]
yrly.columns = ['Year', 'class', 'GPP (tC)', 'Reco (tC)', 'NEE (tC)', 'TSOC (tC)']
yrly['Year'] = yrly['Year'].dt.year
yrly = yrly.set_index('Year')

yrly_grass = yrly.loc[yrly['class']=='grassland'].T
yrly_grass['PFT'] = 'grassland'

yrly_shrub = yrly.loc[yrly['class']=='grass-shrub'].T
yrly_shrub['PFT'] = 'grass-shrub'

yrly_tree = yrly.loc[yrly['class']=='grass-tree'].T
yrly_tree['PFT'] = 'grass-tree'

merged = pd.concat([yrly_grass, yrly_shrub, yrly_tree])
merged.to_csv(f'gs://rangelands/Ranch_Runs/{ranch}/analysis/yearly_summary_stats.csv')

#fig, axes = plt.subplots(2,1,figsize = (3.5, 2.5), sharex=True)

#sns.lineplot(data = yrly_C_stocks, x='year', y='TSOC', color=pal[0],label=None, ax=axes[0])
#sns.lineplot(data = yrly_C_stocks, x='year', y=yrly_C_stocks['AGB'] + yrly_C_stocks['BGB'], color=pal[0],label=None, ax=axes[1])
#despine axes
#axes[0].spines['top'].set_visible(False)
#axes[0].spines['right'].set_visible(False)
#axes[1].spines['top'].set_visible(False)
#axes[1].spines['right'].set_visible(False)

#axes[0].set_xlabel('')
#axes[0].set_ylabel('Soil Carbon\n(gC/$m^2$)', fontsize=10)
#axes[1].set_xlabel('')
#axes[1].set_ylabel('Veg. Carbon\n(gC/$m^2$)', fontsize=10)

#plt.xticks(rotation=45, fontsize=12)

#fig.tight_layout()
#plt.savefig(os.path.join(path_to_temp, 'temp_im_write.png'), dpi=300, transparent=True)
#utils.gs_write_blob(os.path.join(path_to_temp, 'temp_im_write.png'), f'Ranch_Runs/{ranch}/analysis/carbon_summary.png', bucket)


