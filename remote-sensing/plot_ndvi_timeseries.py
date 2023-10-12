import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

path_to_ndvi_timeseries = 'output/HLR/HLR_ndvi.csv'

ndvi = pd.read_csv(path_to_ndvi_timeseries, parse_dates=['Date'])
ndvi['month']=ndvi['Date'].dt.month
ndvi['year'] = ndvi['Date'].dt.year

#pal = sns.color_palette(['#000000', '#242420', '#424235', '#5F5F44', '#7E7C52', '#A99965', '#CFA67C', '#EAAD98', '#FCC7C3', '#FFE6E6')
pal = sns.color_palette(['#424235', '#7E7C52', '#CFA67C'])
fig, ax = plt.subplots(figsize = (20, 3.5))
  
sns.lineplot(data = ndvi.loc[ndvi['pft']=='grass'], x='Date', y='NDVI', color=pal[0],label='grass')
sns.lineplot(data = ndvi.loc[ndvi['pft']=='tree'], x='Date', y='NDVI', color=pal[1],label='tree')
sns.lineplot(data = ndvi.loc[ndvi['pft']=='shrub'], x='Date', y='NDVI', color=pal[2],label='shrub')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.ylabel('NDVI', fontsize=16)
plt.xlabel('')
plt.legend(title = '', frameon=False, fontsize=14)
plt.xticks(rotation=45, fontsize=16)
plt.yticks(fontsize=16)
fig.tight_layout()
plt.savefig(f'output/HLR/HLR_ndvi_v2.jpg', dpi=300)



ndvi_growing = ndvi.loc[(ndvi['month']<=8) & (ndvi['month']>=5)]
ndvi_growing = ndvi_growing[['year', 'NDVI', 'pft']].groupby(by=['pft', 'year']).mean().reset_index()

ndvi_growing['year'] = pd.to_datetime(ndvi_growing['year'], format='%Y')

fig, ax = plt.subplots(figsize = (10, 4))
  
sns.lineplot(data = ndvi_growing.loc[ndvi_growing['pft']=='grass'], x='year', y='NDVI', color=pal[0],label=None)
sns.scatterplot(data = ndvi_growing.loc[ndvi_growing['pft']=='grass'], x='year', y='NDVI', color=pal[0],label='grass')

sns.lineplot(data = ndvi_growing.loc[ndvi_growing['pft']=='tree'], x='year', y='NDVI', color=pal[1],label=None)
sns.scatterplot(data = ndvi_growing.loc[ndvi_growing['pft']=='tree'], x='year', y='NDVI', color=pal[1],label='tree')

sns.lineplot(data = ndvi_growing.loc[ndvi_growing['pft']=='shrub'], x='year', y='NDVI', color=pal[2],label=None)
sns.scatterplot(data = ndvi_growing.loc[ndvi_growing['pft']=='shrub'], x='year', y='NDVI', color=pal[2],label='shrub')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.ylabel('NDVI', fontsize=16)
plt.xlabel('')
plt.legend(loc='upper right', title = '', frameon=False, fontsize=14)
plt.xticks(rotation=45, fontsize=16)
plt.yticks(fontsize=16)
fig.tight_layout()
plt.savefig(f'output/HLR/HLR_ndvi_growing.jpg', dpi=300)



