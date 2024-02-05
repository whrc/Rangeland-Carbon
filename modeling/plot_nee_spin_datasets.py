import geopandas as gpd
import yaml
import pandas as pd
import re
import seaborn as sns
from matplotlib import pyplot as plt
import os
from sklearn.metrics import r2_score
import sys
sys.path.insert(1, '../utils')
import utils

gpp_df=pd.read_csv('output/fPAR_optimization/data/GPP/STARFM_fPAR_optimization_NDVI_GPP_grassland.csv', parse_dates=['time'])

nee_df = pd.read_csv('/home/amullen/Rangeland-Carbon/remote-sensing/output/NEE_SPIN/NEE.g_before_starfm_20230901.csv', parse_dates=['Date'])
nee_df_Bkg = pd.read_csv('/home/amullen/Rangeland-Carbon/remote-sensing/output/NEE_SPIN/NEE.g_before_starfm_site.csv', parse_dates=['Date'])

nee_df_modis = pd.read_csv('/home/amullen/Rangeland-Carbon/remote-sensing/output/NEE_SPIN/NEE.g_before.csv', parse_dates=['Date'])

nee_df['Month'] = nee_df['Date'].dt.month
gpp_sites = gpp_df['site'].unique()

r2_gpp_sites = r2_score(nee_df.loc[nee_df['Site'].isin(gpp_sites),'NEE_measured'], nee_df.loc[nee_df['Site'].isin(gpp_sites),'NEE_p'])
r2_all_sites = r2_score(nee_df['NEE_measured'], nee_df['NEE_p'])

r2_gpp_sites_modis = r2_score(nee_df_modis.loc[nee_df_modis['Site'].isin(gpp_sites),'NEE_measured'], nee_df_modis.loc[nee_df_modis['Site'].isin(gpp_sites),'NEE_p'])
r2_all_sites_modis = r2_score(nee_df_modis['NEE_measured'], nee_df_modis['NEE_p'])

sites = []
r2s = []
r2s_modis = []
gpp_site = []

print(nee_df['Site'].unique())
print(nee_df_modis['Site'].unique())

def find_intersection(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    intersection = set1.intersection(set2)
    return list(intersection)

for site in find_intersection(nee_df['Site'].unique(), nee_df_modis['Site'].unique()):
  sites.append(site)
  site_r2 = r2_score(nee_df.loc[nee_df['Site']== site,'NEE_measured'], nee_df.loc[nee_df['Site']== site,'NEE_p'])
  site_r2_modis = r2_score(nee_df_modis.loc[nee_df_modis['Site']== site,'NEE_measured'], nee_df_modis.loc[nee_df_modis['Site']== site,'NEE_p'])
  is_gpp_site = site in gpp_sites
  
  r2s.append(site_r2)
  r2s_modis.append(site_r2_modis)
  gpp_site.append(is_gpp_site)
  
  fig, ax= plt.subplots(figsize=(8,5))
  sns.lineplot(data=nee_df.loc[nee_df['Site']==site], x='Date', y='NEE_p', label='STARFM', alpha=0.7)
  sns.scatterplot(data=nee_df.loc[nee_df['Site']==site], x='Date', y='NEE_measured', label='STARFM_NEE_measured', alpha=0.6, s=5)
    
  sns.lineplot(data=nee_df_modis.loc[nee_df_modis['Site']==site], x='Date', y='NEE_p', label='MODIS', alpha=0.7)
  sns.scatterplot(data=nee_df_modis.loc[nee_df_modis['Site']==site], x='Date', y='NEE_measured', label='MODIS_NEE_measured', alpha=0.6, s=5)
    
  plt.xticks(rotation = 90)
  fig.tight_layout()
  plt.savefig('output/NEE/plots/{}_NEE_comp.jpg'.format(site), dpi=300)
  plt.show()
  
df = pd.DataFrame({'Site': sites, 'r2': r2s, 'r2_modis': r2s_modis, 'GPP site': gpp_site})
  
print('r2 for gpp sites {}'.format(r2_gpp_sites))
print('r2 for all sites {}'.format(r2_all_sites))

print('modis r2 for gpp sites {}'.format(r2_gpp_sites_modis))
print('modis r2 for all sites {}'.format(r2_all_sites_modis))

fig, ax= plt.subplots(figsize=(8,5))
sns.lineplot(data=df, x='Site', y='r2', label='STARFM')
sns.lineplot(data=df, x='Site', y='r2_modis', label='MODIS', color = sns.color_palette()[2])
plt.xticks(rotation = 90)
fig.tight_layout()

plt.savefig('output/NEE/plots/NEE_r2_by_site.jpg', dpi=300)
plt.show()

months=[]
r2_months=[]
for month in nee_df['Month'].unique():
  months.append(month)
  site_r2 = r2_score(nee_df.loc[nee_df['Month']== month,'NEE_measured'], nee_df.loc[nee_df['Month']== month,'NEE_p'])
  
  r2_months.append(site_r2)

df_monthly = pd.DataFrame({'Month': months, 'r2_months': r2_months})

fig, ax= plt.subplots(figsize=(8,5))
sns.scatterplot(data=df_monthly, x='Month', y='r2_months')
plt.xticks(rotation = 90)
fig.tight_layout()
plt.savefig('output/NEE/plots/NEE_r2_by_month.jpg', dpi=300)
plt.show()

path_to_modis_in = '/home/amullen/Rangeland-Carbon/remote-sensing/output/modis_data_grasslands.csv'

modis_in = pd.read_csv(path_to_modis_in, parse_dates=['Date'])

path_to_starfm_in = 'output/NEE/Ameriflux_STARFM_NEE_grassland.csv'
starfm_in = pd.read_csv(path_to_starfm_in, parse_dates=['Date'])

for site in find_intersection(starfm_in['Site'].unique(), modis_in['Site'].unique()):

  fig, ax= plt.subplots(figsize=(8,5))
  sns.lineplot(data=starfm_in.loc[starfm_in['Site']==site], x='Date', y='STARFM_fPAR', label='STARFM', alpha=0.7)
  sns.lineplot(data=modis_in.loc[modis_in['Site']==site], x='Date', y='fPAR', label='MODIS', alpha=0.7)
  plt.xticks(rotation = 90)
  fig.tight_layout()
  plt.savefig('output/fPAR_optimization/plots/fPAR_comparisons/{}_fPAR_comp.jpg'.format(site), dpi=300)
  plt.show()
    
  #fig, ax= plt.subplots(figsize=(8,5))
  #sns.lineplot(data=starfm_in.loc[starfm_in['Site']==site], x='Date', y='STARFM_GPP', label='GPP', alpha=0.7)
  #plt.xlim(pd.to_datetime('2005-01-01'), pd.to_datetime('2009-01-01'))
  #plt.xticks(rotation = 90)
  #fig.tight_layout()
  #plt.savefig('plots/fPAR_comparisons/{}_GPP.jpg'.format(site), dpi=300)
  #plt.show()
    
path_to_starfm_spin = 'output/SPIN/Ameriflux_STARFM_spin_grassland.csv'
starfm_spin = pd.read_csv(path_to_starfm_spin, parse_dates=['Date'])

for site in find_intersection(starfm_spin['Site'].unique(), modis_in['Site'].unique()):
  fig, ax= plt.subplots(figsize=(8,5))
  sns.lineplot(data=starfm_spin.loc[starfm_spin['Site']==site], x='Date', y='STARFM_fPAR', label='STARFM', alpha=0.7)
  plt.xticks(rotation = 90)
  fig.tight_layout()
  plt.savefig('output/fPAR_optimization/plots/fPAR_SPIN_comparisons/{}_fPAR_comp.jpg'.format(site), dpi=300)
  plt.show()






