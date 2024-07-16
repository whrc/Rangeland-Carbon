import sys
sys.path.append("..")
import numpy as np
import rasterio as rio
import rioxarray as rxr
import xarray as xr
import sys
from RCTM.utils.utils import get_landsat_date_df, get_modis_date_df, gs_write_blob
from google.cloud import storage
import subprocess
import ntpath
import pandas as pd
import os
import argparse


def get_runlist(in_modis_dir, in_landsat_dir, sfm_out_dir, bucket, modland_match_tol=0, modpred_match_tol=365):
  bucket_name = bucket.name

  #initial modpred_match_tol = 15 
  print('getting runlist')
  
  #landsat df
  df_landfiles = get_landsat_date_df(in_landsat_dir, bucket)
  df_landfiles['im_path'] = 'gs://' + bucket_name + '/' + df_landfiles['im_path']
  df_landfiles = df_landfiles.reset_index(drop=True)
  
  #df to hold modis images that match landsat image dates 
  df_modfiles = get_modis_date_df(in_modis_dir, bucket)
  df_modfiles['im_path'] = 'gs://' + bucket_name + '/' + df_modfiles['im_path']
  df_modfiles=df_modfiles.rename(columns={'im_path': 'modis_file'})
  df_modfiles = df_modfiles.reset_index(drop=True)
  df_modfiles['landsat_file'] = [None]*len(df_modfiles)
  df_modfiles['landsat_timediff'] = [None]*len(df_modfiles)
  
  #df to hold modis images for starfm prediction dates 
  df_modpred = get_modis_date_df(in_modis_dir, bucket)
  df_modpred['im_path'] = 'gs://' + bucket_name + '/' + df_modpred['im_path']
  df_modpred=df_modpred.rename(columns={'im_path': 'modpred'})
  df_modpred = df_modpred.reset_index(drop=True)
  
  #intersect modis images with landsat dates
  for index, row in df_modfiles.iterrows():
    timediffs = (df_landfiles['im_date']-row['im_date']).dt.days
    timediffs = timediffs.abs().sort_values()
    closest = timediffs.index[0]
    df_modfiles.loc[index,'landsat_file'] = df_landfiles.loc[closest,'im_path']
    df_modfiles.loc[index,'landsat_timediff'] = timediffs[closest]
  
  #filter matching imagery based on tolerance (minimum should always be zero so may replace this + previous block with intersect on date)
  df_modfiles=df_modfiles[df_modfiles['landsat_timediff']<=modland_match_tol]
  out_dir_series = pd.Series([sfm_out_dir]*len(df_modpred))
  df_modpred['outname'] = out_dir_series.str.cat('modpred_'+ df_modpred['im_date'].astype(str) + '.tif')
  
  df_modpred['landsat1'] = [None]*len(df_modpred)
  df_modpred['landsat2'] = [None]*len(df_modpred)
  df_modpred['modland1'] = [None]*len(df_modpred)
  df_modpred['modland2'] = [None]*len(df_modpred)
  df_modpred['modpred_timediff'] = [None]*len(df_modpred)
  
  for index, row in df_modpred.iterrows():
    timediffs = (df_modfiles['im_date']-row['im_date']).dt.days
    timediffs_before = timediffs.loc[lambda x : x < 0].sort_values()
    timediffs_after = timediffs.loc[lambda x : x > 0].sort_values()
    pair_indices=[]
    pair_diffs=[]
    #print(row['modpred'])
    
    if len(timediffs_before>0):
      before = timediffs_before.index[-1]
      if np.abs(timediffs[before])<modpred_match_tol:
        pair_indices.append(before)
        pair_diffs.append(timediffs[before])
        #print('before: {}, timediff: {}'.format(df_modfiles.loc[before,'landsat_file'], timediffs[before]))
    if len(timediffs_after>0):
      after = timediffs_after.index[0]
      if np.abs(timediffs[after])<modpred_match_tol:
        pair_indices.append(after)
        pair_diffs.append(timediffs[after])
        #print('after: {}, timediff: {}'.format(df_modfiles.loc[after,'landsat_file'], timediffs[after]))
    #print()
    #print()
    if pair_indices==[] or 0 in pair_diffs:
      continue
  
    if len(pair_indices)>=1:
      df_modpred.loc[index,'landsat1'] = df_modfiles.loc[pair_indices[0],'landsat_file']
      df_modpred.loc[index,'modland1'] = df_modfiles.loc[pair_indices[0],'modis_file']
      df_modpred.loc[index,'modpred_timediff'] = timediffs[pair_indices[0]]
    
    if len(pair_indices)==2:
      df_modpred.loc[index,'landsat2'] = df_modfiles.loc[pair_indices[1],'landsat_file']
      df_modpred.loc[index,'modland2'] = df_modfiles.loc[pair_indices[1],'modis_file']
    
  
  df_modpred = df_modpred[df_modpred['landsat1'].notnull()]
  df_modpred = df_modpred[df_modpred['modland1'].notnull()]

  print('done')
  return df_modpred
  
def select_bands(rxr_data, source_path):
    
    #rxr_data = rxr_data.sel(band=[1, 3, 4])
    #rxr_data.attrs['long_name'] = ["blue","red","nir"]
    rxr_data = rxr_data.sel(band=[1, 2])
    rxr_data.attrs['long_name'] = ['red','nir']
    
    return rxr_data

def resample_and_save(source_path, out_path, target_xds = None, landsat1=False):
  
  bands = ['red', 'nir']
  source = rxr.open_rasterio(source_path, chunks='auto')
  source= source.astype('uint16')
  source = select_bands(source, source_path)

  if landsat1==False:
    source = source.rio.reproject_match(target_xds)

  source_mask = source.where(source == 0, 1)

  for i, band in enumerate(bands):

    source_band=source.sel(band=source.band[i])#.astype('uint16')
    source_mask_band=source_mask.sel(band=source_mask.band[i])#.astype('uint16')

    source_band.attrs['long_name'] = [band]
    source_mask_band.attrs['long_name'] = [band]

    if np.nanmax(source_mask_band.values)==0:
      print('all zeros in source image')
      return False

    source_band.rio.to_raster(out_path[:-4]+'_{}.bin'.format(band), driver='ENVI', dtype='uint16')
    source_mask_band.rio.to_raster(out_path[:-4]+'_{}_mask.bin'.format(band), driver='ENVI', dtype='uint8')

  if landsat1==True:
    return source
  
  return True

def run_starfm(runlist, path_to_temp_dir, bucket, starfm_source, starfm_config):
    bands = ['red', 'nir']

    for index, row in runlist.iterrows():
      shape_x = 0
      shape_y = 0
      if row['modland1'] == row['modpred']:
        continue

      num_in_pairs=1
      print('opening landsat1')
      print(row['landsat1'])
      ls1 = resample_and_save(row['landsat1'], os.path.join(path_to_temp_dir, 'ls1_res.bin'), landsat1=True)
      
      if type(ls1)== type(False):
        subprocess.run(['rm', '-r', os.path.join(path_to_temp_dir,'*')])
        continue
      common_transform = ls1.rio.transform
      shape_x = ls1.shape[2]
      shape_y = ls1.shape[1]

      print('resampling...')
      print('mod1')
      print(row['modland1'])
      if resample_and_save(row['modland1'], os.path.join(path_to_temp_dir,'mod1_res.bin'), target_xds=ls1)==False:
        subprocess.run(['rm', '-r', os.path.join(path_to_temp_dir,'*')])
        continue
      
      print('modpred')
      print(row['modpred'])
      if resample_and_save(row['modpred'], os.path.join(path_to_temp_dir,'modpred_res.bin'), target_xds=ls1)==False:
        subprocess.run(['rm', '-r', os.path.join(path_to_temp_dir,'*')])
        continue

      ls2 = False
      mod2 = False

      if row['landsat2']!='None' and row['modland2']!=None:
        num_in_pairs = 2

        if resample_and_save(row['landsat2'], os.path.join(path_to_temp_dir,'ls2_res.bin'), target_xds=ls1)==False:
          num_in_pairs = 1
          continue
        if resample_and_save(row['modland2'], os.path.join(path_to_temp_dir,'mod2_res.bin'), target_xds=ls1)==False:
          num_in_pairs = 1
          continue
        
      
      #run starfm through bands
      for band in bands:

        numland = 'NUM_IN_PAIRS = 1'

        ls1_name = 'landsat1.bin landsat2.bin'
        ls1_new_name = os.path.join(path_to_temp_dir,'ls1_res_{}.bin'.format(band))

        ls1_mask_name = 'landsat1_mask.bin landsat2_mask.bin'
        ls1_mask_new_name = os.path.join(path_to_temp_dir,'ls1_res_{}_mask.bin'.format(band))

        mod1_name = 'modland1.bin modland2.bin'
        mod1_new_name = os.path.join(path_to_temp_dir,'mod1_res_{}.bin'.format(band))

        mod1_mask_name = 'modland1_mask.bin'
        mod1_mask_new_name = os.path.join(path_to_temp_dir,'mod1_res_{}_mask.bin'.format(band))

        modpred_name = 'modpred.bin'
        modpred_new_name = os.path.join(path_to_temp_dir,'modpred_res_{}.bin'.format(band))

        modpred_mask_name = 'modpred_mask.bin'
        modpred_mask_new_name = os.path.join(path_to_temp_dir,'modpred_res_{}_mask.bin'.format(band))

        out_name = 'outname.bin'
        new_outname = os.path.join(path_to_temp_dir,'starfm_test_{}.bin'.format(band))

        resolution= 'RESOLUTION = 30'
        new_resolution = 'RESOLUTION = 30'

        search_dist = 'MAX_SEARCH_DISTANCE = 47'
        new_search_dist = 'MAX_SEARCH_DISTANCE = {}'.format(np.min([shape_x, shape_y]))

        num_rows = 'NROWS = 47'
        new_num_rows = 'NROWS = {}'.format(shape_y)

        num_cols = 'NCOLS = 77'
        new_num_cols = 'NCOLS = {}'.format(shape_x)

        f = open(starfm_config,'r')
        filedata = f.read()
        f.close()        
        newdata = filedata.replace(ls1_name,ls1_new_name)
        newdata = newdata.replace(ls1_mask_name,ls1_mask_new_name)
        newdata = newdata.replace(mod1_name,mod1_new_name)
        newdata = newdata.replace(mod1_mask_name,mod1_mask_new_name)
        newdata = newdata.replace(modpred_name,modpred_new_name)
        newdata = newdata.replace(modpred_mask_name,modpred_mask_new_name)
        newdata = newdata.replace('NUM_IN_PAIRS = 2',numland)

        if num_in_pairs == 2: 
          numland = 'NUM_IN_PAIRS = 2'
          ls2_name = ls1_new_name
          ls2_new_name = ls1_new_name + ' ' + os.path.join(path_to_temp_dir,'ls2_res_{}.bin'.format(band))

          ls2_mask_name = ls1_mask_new_name
          ls2_mask_new_name = ls1_mask_new_name + ' ' + os.path.join(path_to_temp_dir,'ls2_res_{}_mask.bin'.format(band))
        
          mod2_name = mod1_new_name
          mod2_new_name = mod1_new_name + ' ' + os.path.join(path_to_temp_dir,'mod2_res_{}.bin'.format(band))

          mod2_mask_name = mod1_mask_new_name
          mod2_mask_new_name = mod1_mask_new_name + ' ' + os.path.join(path_to_temp_dir,'mod2_res_{}_mask.bin'.format(band))
          
          newdata = newdata.replace(ls2_name,ls2_new_name)
          newdata = newdata.replace(mod2_name,mod2_new_name)
          newdata = newdata.replace(ls2_mask_name,ls2_mask_new_name)
          newdata = newdata.replace(mod2_mask_name,mod2_mask_new_name)
          newdata = newdata.replace('NUM_IN_PAIRS = 1',numland)
        
        newdata = newdata.replace(out_name,new_outname)
        newdata = newdata.replace(resolution,new_resolution)        
        newdata = newdata.replace(num_rows,new_num_rows)
        newdata = newdata.replace(num_cols,new_num_cols)
        newdata = newdata.replace(search_dist,new_search_dist)

        
        print('running starfm')
        
        f = open(os.path.join(starfm_source, 'starfm_config.txt'),'w')
        f.write(newdata)
        f.close()
        subprocess.run([os.path.join(starfm_source, 'StarFM.exe'), os.path.join(starfm_source, 'starfm_config.txt')])
        #!./StarFM.exe input_ref_test.txt

      if os.stat(os.path.join(path_to_temp_dir,'starfm_test_{}.bin'.format(bands[0]))).st_size==0:
        continue

      #blue=rxr.open_rasterio('temp/starfm_test_{}.bin'.format(bands[0]))
      #blue.name='blue'
      red=rxr.open_rasterio(os.path.join(path_to_temp_dir,'starfm_test_{}.bin'.format(bands[0])))
      red.name='red'
      nir=rxr.open_rasterio(os.path.join(path_to_temp_dir,'starfm_test_{}.bin'.format(bands[1])))
      nir.name='nir'
      combined = xr.concat([red, nir], dim='band')
      combined.attrs['long_name'] = bands
      combined = combined.rio.write_crs('EPSG:4326')
      combined.rio.transform = common_transform
      
      combined.rio.to_raster(os.path.join(path_to_temp_dir,'out_sfm.tif'), dtype=np.uint16, driver='GTiff')
      gs_write_blob(os.path.join(path_to_temp_dir,'out_sfm.tif'), row['outname'], bucket)
      
      subprocess.run(['rm', '-r', os.path.join(path_to_temp_dir,'*')])
      
      print('done')


if __name__ == "__main__":
  parser=argparse.ArgumentParser()
  
  parser.add_argument("--in_modis_dir", help="directory with input modis imagery")
  parser.add_argument("--in_landsat_dir", help="directory with input landsat imagery")
  parser.add_argument("--sfm_out_dir", help="directory to output STARFM")
  parser.add_argument("--bucket_name", help="bucket name")

  args=parser.parse_args()
  
  runlist=get_runlist(args.in_modis_dir, args.in_landsat_dir, args.sfm_out_dir, args.bucket_name)
  run_starfm(runlist, args.bucket_name)

