import numpy as np

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
  LUEmax = params['LUEmax'][0]
  MIN_Tmn = params['MIN_Tmn'][0]
  MAX_Tmx = params['MAX_Tmx'][0]
  MIN_VPD = params['MIN_VPD'][0]
  MAX_VPD = params['MAX_VPD'][0]
  MIN_SMrz = params['MIN_SMrz'][0]
  MAX_SMrz = params['MAX_SMrz'][0]
  tmult = get_mult(MIN_Tmn, MAX_Tmx, tsoil_in)
  smult = get_mult(MIN_SMrz, MAX_SMrz, sm_in)
  wmult = get_mult(MIN_VPD, MAX_VPD, vpd_in)

  LUE = LUEmax * tmult * smult * wmult
  GPP = LUE * SW_IN * fPAR * 0.45
  if np.isnan(np.sum(LUE)):
    
    print(f'LUEmax: {LUEmax[np.isnan(LUE)]}, tmult: {tmult[np.isnan(LUE)]}, smult: {smult[np.isnan(LUE)]}, wmult: {wmult[np.isnan(LUE)]}')
  if np.isnan(np.sum(GPP)):
    
    print(f'LUE: {LUE[np.isnan(GPP)]}, SW_IN: {SW_IN[np.isnan(GPP)]}, fPAR: {fPAR[np.isnan(GPP)]}')

  return(GPP)
  
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
  