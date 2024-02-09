import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
import sys
sys.path.insert(1, '../utils')
import utils


utils.authorize()

date_range = pd.DataFrame({'im_date': pd.date_range(start='2002-01-01', end='2004-07-19', freq='4d').to_series()}).reset_index(drop=True)

  
path_to_covariates = 'gs://rangelands/Ameriflux/NEEINS_AVG_ALL/A32.csv'
path_to_covariate_test = 'gs://rangelands/Ameriflux_sites/A32_starfm/covariates_v2/A32_covariates_site_test.csv'

covariates_v1 = pd.read_csv(path_to_covariates, parse_dates=['date'])
covariates_v1 = date_range.merge(covariates_v1, how='left', left_on='im_date', right_on='date')

covariates_v2 = pd.read_csv(path_to_covariate_test, parse_dates=['time'])
covariates_v2 = date_range.merge(covariates_v2, how='left', left_on='im_date', right_on='time')

covariates_v1 = covariates_v1.loc[covariates_v1['date']<=covariates_v2['time'].max()]
covariates_v1 = covariates_v1.loc[covariates_v1['date']>=covariates_v2['time'].min()]

covariates_v2 = covariates_v2.loc[covariates_v2['time']<=covariates_v1['date'].max()]
covariates_v2 = covariates_v2.loc[covariates_v2['time']>=covariates_v1['date'].min()]

print(covariates_v1.columns)
print(covariates_v2.columns)

fig, ax=plt.subplots()
sns.lineplot(data=covariates_v1, x='date', y='TA', label='v1', alpha=0.7)
#sns.lineplot(data=covariates_v2, x='time', y='tavg', label='v2', alpha=0.7)
plt.savefig('output/covariate_site_test/tavg_4d.jpg')
plt.show()

fig, ax=plt.subplots()
sns.lineplot(data=covariates_v1, x='date', y='SWC1', label='v1', alpha=0.7)
#sns.lineplot(data=covariates_v2, x='time', y='sm1', label='v2', alpha=0.7)
plt.savefig('output/covariate_site_test/sm1_4d.jpg')
plt.show()

fig, ax=plt.subplots()
sns.lineplot(data=covariates_v1, x='date', y='SWC2', label='v1', alpha=0.7)
#sns.lineplot(data=covariates_v2, x='time', y='sm2', label='v2', alpha=0.7)
plt.savefig('output/covariate_site_test/sm2_4d.jpg')
plt.show()

fig, ax=plt.subplots()
sns.lineplot(data=covariates_v1, x='date', y='SW_IN_NLDAS', label='v1', alpha=0.7)
#sns.lineplot(data=covariates_v2, x='time', y='shortwave_radition', label='v2', alpha=0.7)
plt.savefig('output/covariate_site_test/sw_in_4d.jpg')
plt.show()

fig, ax=plt.subplots()
sns.lineplot(data=covariates_v1, x='date', y='TA_min', label='v1', alpha=0.7)
#sns.lineplot(data=covariates_v2, x='time', y='tmin', label='v2', alpha=0.7)
plt.savefig('output/covariate_site_test/tmin_4d.jpg')
plt.show()

fig, ax=plt.subplots()
sns.lineplot(data=covariates_v1, x='date', y='TA_min', label='v1', alpha=0.7)
#sns.lineplot(data=covariates_v2, x='time', y='tmin', label='v2', alpha=0.7)
plt.savefig('output/covariate_site_test/tmin_4d.jpg')
plt.show()

fig, ax=plt.subplots()
sns.lineplot(data=covariates_v1, x='date', y='TS', label='v1', alpha=0.7)
#sns.lineplot(data=covariates_v2, x='time', y='tsoil', label='v2', alpha=0.7)
plt.savefig('output/covariate_site_test/tsoil_4d.jpg')
plt.show()

fig, ax=plt.subplots()
sns.lineplot(data=covariates_v1, x='date', y='VPD', label='v1', alpha=0.7)
#sns.lineplot(data=covariates_v2, x='time', y='vpd', label='v2', alpha=0.7)
plt.savefig('output/covariate_site_test/vpd_4d.jpg')
plt.show()

fig, ax=plt.subplots()
sns.lineplot(data=covariates_v1, x='date', y='ppt', label='v1', alpha=0.7)
#sns.lineplot(data=covariates_v2, x='time', y='prcp', label='v2', alpha=0.7)
plt.savefig('output/covariate_site_test/prcp_4d.jpg')
plt.show()