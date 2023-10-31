import numpy as np
import pandas as pd
import xarray as xr

def xr_dataset_to_data_array(dataset):

  band_names = list(dataset.attrs['long_name'])
  band_name_dict = {}
  for i, band_name in enumerate(band_names):
    band_name_dict[i+1] = band_name 
  dataset = dataset.to_dataset(dim='band')
  dataset = dataset.rename_vars(name_dict=band_name_dict)
  return dataset
  
def calc_sr(index):
  """calculate simple ratio

  Args:
    index (numeric, or numeric array-like): index to calculate simple ratio

  Returns:
     simple ratio of same type as input
  """ 
  return (1 + index)/(1-index)
  
def calc_bias(observed, predicted):
  """calculate bias

  Args:
    observed (numeric, or numeric array-like): observation, or array of observations
    predicted (numeric, or numeric array-like): predictions, or array of predictions
  
  Returns:
     (numeric) bias
  """ 
  return np.nansum((predicted-observed))/len(predicted)
  
def calc_fPAR(this_index, min_index, max_index, this_sr, min_sr, max_sr, min_fpar, max_fpar):
  """calculate fPAR - fraction of photosynthetically active radiation absorbed by vegetation

  Args:
    this_index (numeric, or numeric array-like): observed evi or ndvi
    min_index (numeric, or numeric array-like): minimum evi/ndvi (reference value), calculated as 2nd percentile of all evi/ndvi values
    max_index (numeric, or numeric array-like): maximum evi/ndvi (reference value), calculated as 98th percentile of all evi/ndvi values
    this_sr (numeric, or numeric array-like): sr of evi/ndvi
    min_sr (numeric, or numeric array-like): minimum sr of evi/ndvi (reference value), calculated as 2nd percentile of all evi/ndvi sr values
    max_sr (numeric, or numeric array-like): maximum sr of evi/ndvi (reference value), calculated as 98th percentile of all evi/ndvi sr values
    min_fpar (numeric, or numeric array-like): reference value for minimum fpar, optimized using site observations
    max_fpar (numeric, or numeric array-like): reference value for maximum fpar, optimized using site observations
  
  Returns:
     (numeric) fPAR
  """ 
  
  fpar_EVI = (((this_index - min_index)*(max_fpar - min_fpar))/(max_index - min_index)) + min_fpar
  
  fpar_SR = (((this_sr - min_sr)*(max_fpar - min_fpar))/(max_sr - min_sr)) + min_fpar
  
  fpar = (fpar_SR + fpar_EVI)/2
  if isinstance(fpar, xr.core.dataarray.DataArray):
    fpar = xr.where(fpar >= 0, fpar, 0)
  else:
    fpar[fpar < 0] = 0
  
  return fpar
  
def get_mult(min, max, in_val):
  """scales input data between 0 and 1 based on minimum and maximum reference values, if value > maximum, scaled value=1, if value < minimum, scaled value = 0

  Args:
    in_val (numeric array-like): input value to scale
    min (numeric): minimum (reference value)
    max (numeric): maximum (reference value)
  
  Returns:
     (numeric) fPAR
  """ 
  if isinstance(in_val, xr.core.dataarray.DataArray):
  
    out = xr.where(in_val > max, 1, 0)
    out = xr.where((in_val > min) & (in_val < max), (in_val-min)/(max-min), out)
  
  else:
    out = np.zeros(len(in_val))
    out[in_val>max] = 1 #if val>max, set to 1
    out[(in_val>min) & (in_val<max)] = (in_val[(in_val>min) & (in_val<max)]-min)/(max-min) #if val in between min and max, scale between min and max
  return out
  
def calcGPP(params, tsoil_in, sm_in, vpd_in, SW_IN, fPAR):
  """calculate GPP using the fPAR method

  Args:
    params (dictionary): dict containing key-value pairs for parameters
    tsoil_in (numeric array-like): soil temperature (unit)
    sm_in (numeric array-like): soil moisture (unit)
    vpd_in (numeric array-like): vapor pressure defecit (unit)
    SW_IN (numeric array-like): Incoming shortwave radiation (w/m2)
    fPAR (numeric array-like): fPAR (unitless)
  
  Returns:
     (numeric array-like) GPP
  """
  LUEmax = params['LUEmax']
  MIN_Tmn = params['MIN_Tmn']
  MAX_Tmx = params['MAX_Tmx']
  MIN_VPD = params['MIN_VPD']
  MAX_VPD = params['MAX_VPD']
  MIN_SMrz = params['MIN_SMrz']
  MAX_SMrz = params['MAX_SMrz']
  tmult = get_mult(MIN_Tmn, MAX_Tmx, tsoil_in)
  smult = get_mult(MIN_SMrz, MAX_SMrz, sm_in)
  wmult = get_mult(MIN_VPD, MAX_VPD, vpd_in)
  LUE = tmult * smult * wmult * LUEmax
  GPP = LUE * SW_IN * fPAR * 0.45
  
  return GPP
  
def get_site_pft(site, site_pfts):
  """get plant functional type (PFT) corresponding to site, based on RCTM parameter file

  Args:
    site (string): site string representation
    site_pfts (dictionary): dict containing key-value pairs
  
  Returns:
     (string) PFT
  """ 
  pft = None
  for key in site_pfts:
      if site in site_pfts[key]:
        pft = key
        return pft
  if pft==None:
    print(f'{site} pft not found')
  
def read_in_sites_as_df(path_to_footprint_list, params, bucket_name, ref_value_calc='all'):
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
      df['ndvi_sr'] = calc_sr(df['ndvi'])
      df['site'] = site
      
      if ref_value_calc=='site':
        #calculate reference values for min and max evi/ndvi based on site
        df_ref_values = df.groupby(by=['site']).agg({'ndvi':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)],
                                                              'ndvi_sr':[lambda x: x.quantile(0.02), lambda x: x.quantile(0.98)]})                                    
        df_ref_values = df_ref_values.reset_index()
        df_ref_values.columns = ['site', 'ndvi_02', 'ndvi_98', 'ndvi_sr_02', 'ndvi_sr_98']
      
      ##### special case sites have combined imagery #####
      if ref_value_calc=='site':
        df = pd.merge(df, df_ref_values, on='site').reset_index()
  
      if site=='Rwe_Rms':
        rwe=df
        rwe['site']= 'Rwe'
        rwe['pft'] = get_site_pft('Rwe', params['Ameriflux_site_pfts'])
        rms=df.copy()
        rms['site']= 'Rms'
        rms['pft'] = get_site_pft('Rms', params['Ameriflux_site_pfts'])
        df_sites.append(rwe)
        df_sites.append(rms)
        
      elif site=='ARbc':
        ARb=df
        ARb['site']= 'ARb'
        ARb['pft'] = get_site_pft('ARb', params['Ameriflux_site_pfts'])
        ARc=df.copy()
        ARc['site']= 'ARc'
        ARc['pft'] = get_site_pft('ARc', params['Ameriflux_site_pfts'])
        df_sites.append(ARb)
        df_sites.append(ARc)
        
      elif site=='Hn2_3':
        Hn2=df
        Hn2['site']= 'Hn2'
        Hn2['pft'] = get_site_pft('Hn2', params['Ameriflux_site_pfts'])
        Hn3=df.copy()
        Hn3['site']= 'Hn3'
        Hn3['pft'] = get_site_pft('Hn3', params['Ameriflux_site_pfts'])
        df_sites.append(Hn2)
        df_sites.append(Hn3)
        
      elif site=='KM2_3':
        KM2=df
        KM2['site']= 'KM2'
        KM2['pft'] = get_site_pft('KM2', params['Ameriflux_site_pfts'])
        KM3=df.copy()
        KM3['site']= 'KM3'
        KM3['pft'] = get_site_pft('KM3', params['Ameriflux_site_pfts'])
        df_sites.append(KM2)
        df_sites.append(KM3)
        
      elif site=='LS1_2':
        LS1=df
        LS1['site']= 'LS1'
        LS1['pft'] = get_site_pft('LS1', params['Ameriflux_site_pfts'])
        LS2=df.copy()
        LS2['site']= 'LS2'
        LS2['pft'] = get_site_pft('LS2', params['Ameriflux_site_pfts'])
        df_sites.append(LS1)
        df_sites.append(LS2)
        
      elif site=='CZ1_xSJ':
        df['site']= 'xSJ'
        df['pft'] = get_site_pft('xSJ', params['Ameriflux_site_pfts'])
        df_sites.append(df)
        
      elif site=='Fpe':
        df['site']= 'FPe'
        df['pft'] = get_site_pft('FPe', params['Ameriflux_site_pfts'])
        df_sites.append(df)
        
      else:
        df['pft'] = get_site_pft(site, params['Ameriflux_site_pfts']) 
        
        df_sites.append(df)
        
    #merge monthly and yearly reference values back to original df  
    df_sites=pd.concat(df_sites, ignore_index=True)
    
    if ref_value_calc == 'all':
      # all data grouped reference value calculation
      df_sites['ndvi_02_yr'] = df_sites['ndvi'].quantile(0.02)
      df_sites['ndvi_98_yr'] = df_sites['ndvi'].quantile(0.98)
      df_sites['ndvi_sr_02_yr'] = df_sites['ndvi_sr'].quantile(0.02)
      df_sites['ndvi_sr_98_yr'] = df_sites['ndvi_sr'].quantile(0.98)
    
    if ref_value_calc == 'pft': 
      # pft based reference value caluculation
      for pft in df_sites['pft'].unique():
        df_sites.loc[df_sites['pft']==pft, 'ndvi_02_yr'] = df_sites.loc[df_sites['pft']==pft, 'ndvi'].quantile(0.02)
        df_sites.loc[df_sites['pft']==pft, 'ndvi_98_yr'] = df_sites.loc[df_sites['pft']==pft, 'ndvi'].quantile(0.98)
        df_sites.loc[df_sites['pft']==pft, 'ndvi_sr_02_yr'] = df_sites.loc[df_sites['pft']==pft, 'ndvi_sr'].quantile(0.02)
        df_sites.loc[df_sites['pft']==pft, 'ndvi_sr_98_yr'] = df_sites.loc[df_sites['pft']==pft, 'ndvi_sr'].quantile(0.98)

    df_sites = df_sites.loc[~df_sites['ndvi'].isna()]
    
    return df_sites