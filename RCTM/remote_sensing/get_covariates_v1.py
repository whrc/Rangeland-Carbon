import ee
import utils
import time 
import numpy as np
from functools import partial
import preprocess_starfm_imagery as psi
from google.cloud import storage
from datetime import datetime as dt
import pandas as pd

utils.authorize()
utils.print_root_assets()
bucket_name = 'rangelands'



path_to_footprint_list = 'res/sites_v2_rest.txt'
#roi_asset_path = 'projects/rangelands-explo-1571664594580/assets/Ameriflux_RS/KFS/KFS'
starfm_dir='Ameriflux_sites/KFS_starfm/starfm_test_v2/'
daymet_outdir='Ameriflux_sites/KFS_starfm/daymet/'

start_date = '2002-01-01'
end_date = '2023-01-01'
#start_date = '2002-01-01'
#end_date = '20022-01-01'

daymet_collection_path = 'NASA/ORNL/DAYMET_V4'

def get_daymet(start_date, end_date, roi_asset_path, daymet_collection_path):

  roi=ee.FeatureCollection(roi_asset_path)
  DAYMET = ee.ImageCollection(daymet_collection_path)
  
  temp = DAYMET.filter(ee.Filter.date(start_date, end_date)).filterBounds(roi.geometry()).select(["tmax","tmin"])
  
  def calcmean(img):
    tmean = img.expression(
      'tmean = ((Tmax + Tmin) / 2)', {
      'Tmax': img.select("tmax"),
      'Tmin': img.select("tmin")
    });
    return(img.addBands(tmean).select('tmean'))
  
  tavg = temp.map(calcmean)
  
  return tavg
  
def get_gee_queue():

  ops = ee.data.listOperations()
  ids = [ op['name'] for op in ops ]
  state = [ op['metadata']['state'] for op in ops ]
  op_df = pd.DataFrame({'name': ids, 'state': state})
  queue_len = len(op_df[(op_df['state']=='PENDING') | (op_df['state']=='RUNNING')])
  available_slots = 3000-queue_len
  
  return available_slots

def export_daymet(daymet_collection, bucket_name, roi_asset_path, starfm_dir, out_directory_path, overwrite=False):
  storage_client = storage.Client.from_service_account_json('gee_key.json')
  bucket = storage_client.get_bucket(bucket_name)
  
  roi=ee.FeatureCollection(roi_asset_path)
  count = int(daymet_collection.size().getInfo())
  names = daymet_collection.aggregate_array('system:index').getInfo()
  
  sfm_df = psi.get_sfm_date_df(starfm_dir, bucket)
  sfm_dates = sfm_df['im_date'].dt.strftime('%Y%m%d').to_list()
  
  print(sfm_dates[0:10])
  print(names[0:10])
  
  available_slots = get_gee_queue()
  
  for i in range(0, count):
    
    exists = storage.Blob(bucket=bucket, name='{}DAYMET_{}.tif'.format(out_directory_path, names[i])).exists(storage_client)

    if exists:
      if overwrite==False:
        print('DAYMET_{} already exists! Skipping'.format(names[i]))
        continue
    
    if available_slots<=0:
      print('GEE queue full, waiting')
      time.sleep(300)
      available_slots = get_gee_queue()
      
    if available_slots>0:
    
      if names[i] in sfm_dates:
          
        image = daymet_collection.filter(ee.Filter.eq('system:index', names[i])).first()
        print()
        print('exporting {} of {}: {}'.format(i+1, count, names[i]))
        task = ee.batch.Export.image.toCloudStorage(
                  image=image,
                  scale=30,
                  #crs=image.projection().crs(),
                  #crsTransform = image.projection().getInfo()['transform'],
                  crs='EPSG:4326',
                  region=roi.geometry(),
                  bucket=bucket_name,
                  fileNamePrefix='{}DAYMET_{}'.format(out_directory_path, names[i]))
                
        task.start()
        available_slots-=1
    
  return

with open(path_to_footprint_list) as f:
  
  roi_paths = [line.rstrip('\n') for line in f]
  #print('sites: {}'.format(sites))
  geometries=[]
  
  for roi_path in roi_paths:
    print(roi_path)    
    prefix='Ameriflux_sites'
    site = roi_path.split('/')[-1]
    roi=ee.FeatureCollection(roi_path)
    geometries.append(roi.geometry().getInfo())
    
    daymet_outdir = f'{prefix}/{site}_starfm/daymet/'
    starfm_dir= f'{prefix}/{site}_starfm/starfm_test_smooth_v2/'
    daymet_collection=get_daymet(start_date, end_date, roi_path, daymet_collection_path)
    export_daymet(daymet_collection, bucket_name, roi_path, starfm_dir, daymet_outdir)