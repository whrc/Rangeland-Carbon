import ee
import time 
import numpy as np
from functools import partial
import starfm_preprocessing as psi
from google.cloud import storage
from datetime import datetime as dt
import pandas as pd
import sys
sys.path.insert(1, '../utils')
import utils


utils.authorize()
#utils.print_root_assets()
bucket_name = 'rangelands'

#roi_asset_path = 'projects/rangelands-explo-1571664594580/assets/Ameriflux_RS/Rws/Rws'
#starfm_dir='Ameriflux_sites/Rws_starfm/starfm_test/'
#daymet_outdir='Ameriflux_sites/Rws_starfm/daymet/'

covariate_mapping = {
  'meteorology': 'NASA/ORNL/DAYMET_V4', #daymet
  'soil_texture': 'projects/rangelands-explo-1571664594580/assets/Covariates/Soils/SoilGrid100',
  'soil_temperature': 'projects/rangelands-explo-1571664594580/assets/Covariates/Soils/NLDAST', #NLDAST
  'soil_moisture': 'projects/rangelands-explo-1571664594580/assets/Covariates/Soils/NLDAS', #NLDAS
  'shortwave': 'NASA/NLDAS/FORA0125_H002',
  'NLCD_landcover': 'USGS/NLCD_RELEASES/2019_REL/NLCD',
  'RAP_landcover': 'projects/rap-data-365417/assets/vegetation-cover-v3'
}

#Filter out images that are missing required bands
def band_num(img):
  return img.set("band_number", img.bandNames().size())

def vpdFunc(img):
    vpd = img.select("vpd").multiply(1000)
    img = img.select(['tmmn','tmmx','srad','pr']).addBands(vpd)
    return img.select(['tmmn','tmmx','srad','pr','vpd'])

def tminFunc(img):
    tmmn = img.select("tmmn").subtract(273.15)
    img = img.select(['vpd','tmmx','srad','pr']).addBands(tmmn)
    return img.select(['tmmn','tmmx','srad','pr','vpd'])
    
def tmaxFunc(img):
    tmmx = img.select("tmmx").subtract(273.15)
    img = img.select(['vpd','tmmn','srad','pr']).addBands(tmmx)
    return img.select(['tmmn','tmmx','srad','pr','vpd'])
    
def getVPD(img): # final unit is millibar
    ea = img.select('vp').divide(100.0).rename('ea')
    es_min = img.expression(
    'es_min = 6.1078 * exp(17.269*T / (237.3+T ))', { 
      'T': img.select('tmin')
    })
    es_max = img.expression(
    'es_max = 6.1078 * exp(17.269*T / (237.3+T ))', {
      'T': img.select('tmax')
    })
    img = img.addBands([ea, es_min, es_max])
    es_avg = img.expression(
    'es_avg = 6.1078 * exp(17.269*T / (237.3+T))',{
      'T':img.select('tavg')
    })
    img = img.addBands(es_avg)
    vpd = img.expression(
    'vpd = es - ea',{
      'es':img.select('es_avg'),
      'ea':img.select('ea')
    }).addBands(img.metadata('system:time_start').divide(1e18).rename('time'))
    img = img.addBands(vpd)
    return img.select(['time','vpd']) # unit is hpa - millibar
    
def getSrad(img):
    mj = img.select("srad").divide(1000000/86400)
    img = img.select('tmax').addBands(mj)
    return img.select('srad')

def getTavg(img):
    avg = img.expression(
    'tavg = (tmin + tmax)/2', 
    {
        'tmin': img.select('tmin'),
        'tmax': img.select('tmax')})
    return img.select(['tmin','tmax']).addBands(avg)
    
def tsoilFunc(img):
    ts = img.expression(
    'tsoil = T - 273.15', 
    {
      'T': img.select(["B0"],["TS"])
    })
    return ts.copyProperties(img, ['system:time_start'])

def getSM1(img, soilgridbd):
    img = img.addBands(soilgridbd)
    sm1 = img.expression(
    'sm1 = (surfacesm/100)/(1 - (bd1/2650))',
    {
      'surfacesm':img.select('B0'),
      'bd1':img.select('0_10')
    })
    return img.addBands(sm1).select('sm1')
    
def getSM2(img, soilgridbd):
    img = img.addBands(soilgridbd)
    img = img.addBands(img.expression(
      'bd2 = (b2 + b3)/2',
      {
        'b2':img.select('10_20'),
        'b3':img.select('20_40')
    })).select('B0','B1','bd2');
    
    sm2 = img.expression(
    'sm2 = ((rootsm+surfacesm)/400)/(1 - (bd2/2650))',
    {
      'surfacesm':img.select('B0'),
      'rootsm':img.select('B1'),
      'bd2':img.select('bd2')
    })
    return img.addBands(sm2).select("sm2")
    
def nldasFunc(img):
    swir = img.expression(
    'SW_IN_NLDAS = (SWIN/1000000*86400/12)', 
    { #//convert to average hourly during the daylight time to match with Daymet4
      'SWIN': img.select('shortwave_radiation')
    })
    return swir.copyProperties(img, ['system:time_start'])


def combInputs(date, daymet_inputs, smColl1, smColl2, tsoil, NLDAS, tavg, tmin, prcp, clay):
    date_f1 = dt.strptime(date, '%Y%m%d')
    date_f2 = date_f1.strftime('%Y-%m-%d')
    date = ee.Date(date_f2)
    date2 = date.advance(1, 'days')
    daymet = daymet_inputs.filterDate(date).mean().toFloat()
    sm1 = smColl1.filterDate(date).mean().toFloat()
    sm2 = smColl2.filterDate(date).mean().toFloat()
    st = tsoil.filterDate(date).mean().toFloat()
    swir = NLDAS.filter(ee.Filter.date(date, date2)).sum().toFloat()
    Tavg = tavg.filterDate(date).mean().toFloat()
    Tmn = tmin.filterDate(date).mean().toFloat()
    ppt = prcp.filterDate(date).mean().toFloat()
    return daymet.addBands(st).addBands(sm1).addBands(sm2).addBands(swir).addBands(Tavg).addBands(Tmn).addBands(ppt).addBands(clay.toFloat())
    
def clip(image, geometry):
  return image.clip(geometry)

def get_datasets(covariate_mapping, roi_asset_path, start_date, end_date):

  ###### acquisition #######
  roi=ee.FeatureCollection(roi_asset_path).geometry()
  
  daymet = ee.ImageCollection(covariate_mapping['meteorology'])
  daymet = daymet.filterDate(start_date, end_date).filterBounds(roi).select(['tmin','tmax','srad','vp','prcp'])
  
  SoilGrid = ee.ImageCollection(covariate_mapping['soil_texture'])
  SoilGrid = SoilGrid.filterBounds(roi)
  bd1 = SoilGrid.filterMetadata('system:index', 'equals','bd_d2').toBands().rename("bd1") 
  bd2 = SoilGrid.filterMetadata('system:index', 'equals','bd_d3').toBands().rename("bd2") 
  bd3 = SoilGrid.filterMetadata('system:index', 'equals','bd_d4').toBands().rename("bd3") 
  clay1 = SoilGrid.filterMetadata('system:index', 'equals','clay_d2').toBands().rename("clay1") 
  clay2 = SoilGrid.filterMetadata('system:index', 'equals','clay_d3').toBands().rename("clay2") 
  clay3 = SoilGrid.filterMetadata('system:index', 'equals','clay_d4').toBands().rename("clay3") 
  soilgridbd = ee.Image(bd1).addBands(bd2).addBands(bd3)\
               .select(['bd1','bd2','bd3'],['0_10','10_20','20_40'])
  soilgridclay = ee.Image(clay1).addBands(clay2).addBands(clay3)\
               .select(['clay1','clay2','clay3'],['0_10','10_20','20_40'])
  clay = (soilgridclay.select('0_10').add(soilgridclay.select('10_20'))\
              .add(soilgridclay.select('20_40'))).divide(3).rename("clay")\
              .set('system:time_start', ee.Date('2019-01-01').millis())
  
  smColl = ee.ImageCollection(covariate_mapping['soil_moisture'])
  smColl = smColl.filterDate(start_date, end_date)\
                .filterBounds(roi).select(['B0','B1'])
  
  stColl = ee.ImageCollection(covariate_mapping['soil_temperature'])
  stColl = stColl.filterDate(start_date, end_date)\
                .filterBounds(roi).select(['B0'])
                
  NLDAS = ee.ImageCollection(covariate_mapping['shortwave'])  
  NLDAS = NLDAS.filterDate(start_date, end_date)\
                .filterBounds(roi).select(['shortwave_radiation']).map(nldasFunc)
                
  ###### processing #######
  temp = daymet.map(getTavg)
  tmin = temp.select('tmin')
  tmax = temp.select('tmax')
  tavg = temp.select('tavg')
  srad = daymet.map(getSrad)
  vpdColl = daymet.select('vp').combine(temp)
  VPD = vpdColl.map(getVPD)
  prcp = daymet.select('prcp')
  daymet_inputs = srad.combine(VPD)
  tsoil = stColl.map(tsoilFunc)
  
  smColl1_fill = partial(getSM1, soilgridbd=soilgridbd)
  smColl2_fill = partial(getSM2, soilgridbd=soilgridbd)
  smColl1 = smColl.map(smColl1_fill)
  smColl2 = smColl.map(smColl2_fill)
  
  return daymet_inputs, smColl1, smColl2, tsoil, NLDAS, tavg, tmin, prcp, clay

def get_daymet(start_date, end_date, roi_asset_path, daymet_collection_path):
  """Fetch DAYMET mean temperature from Earth Engine

  Args:
    start_date (string): acquisition start date
    end_date (string): acquisition end date
    roi_asset_path (string): path to ROI stored as Earth Engine Asset
    daymet_collection_path (string): google storage bucket name

  Returns:
     None
  """ 
  roi=ee.FeatureCollection(roi_asset_path)
  DAYMET = ee.ImageCollection(daymet_collection_path)
  
  temp = DAYMET.filter(ee.Filter.date(start_date, end_date)).filterBounds(roi.geometry()).select(['tmax','tmin','srad','vp','prcp'])
  
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
  i=0
  for op in ops:
    if op['metadata']['state']=='SUCCEDED':
      print(op)
      i+=1
    if i>200:
      break
  ids = [ op['name'] for op in ops ]
  state = [ op['metadata']['state'] for op in ops ]
  op_df = pd.DataFrame({'name': ids, 'state': state})
  queue_len = len(op_df[(op_df['state']=='PENDING') | (op_df['state']=='RUNNING')])
  available_slots = 3000-queue_len
  
  return available_slots

def export_landcover(bucket_name, roi_asset_path, out_directory):

  roi=ee.FeatureCollection(roi_asset_path)
  
  NLCD = ee.ImageCollection(covariate_mapping['NLCD_landcover']).filterBounds(roi.geometry()).filterDate('2019-01-01', '2020-01-01').first().clip(roi)
  RAP = ee.ImageCollection(covariate_mapping['RAP_landcover']).filterBounds(roi.geometry()).filterDate('2019-01-01', '2020-01-01').first().clip(roi)
  print('exporting NLCD')
  task = ee.batch.Export.image.toCloudStorage(
                  image=NLCD,
                  scale=30,
                  crs='EPSG:4326',
                  region=roi.geometry(),
                  bucket=bucket_name,
                  fileNamePrefix='{}NLCD_2019'.format(out_directory))
  task.start()
    
  print('exporting RAP')              
  task = ee.batch.Export.image.toCloudStorage(
                  image=RAP,
                  scale=30,
                  crs='EPSG:4326',
                  region=roi.geometry(),
                  bucket=bucket_name,
                  fileNamePrefix='{}RAP_2019'.format(out_directory))
                  
  task.start()
  
  return

def export_daymet(covariate_tuple, bucket_name, roi_asset_path, modis_dir, out_directory_path):

  daymet_inputs, smColl1, smColl2, tsoil, NLDAS, tavg, tmin, prcp, clay = covariate_tuple
  
  storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
  bucket = storage_client.get_bucket(bucket_name)
  
  roi=ee.FeatureCollection(roi_asset_path)
  count = int(daymet_inputs.size().getInfo())
  names = daymet_inputs.aggregate_array('system:index').getInfo()
  
  modis_df = psi.get_modis_date_df(modis_dir, bucket)
  modis_dates = modis_df['im_date'].dt.strftime('%Y%m%d').to_list()
  
  available_slots = get_gee_queue()
  
  for i in range(0, count):
  
    if available_slots<=0:
      print('GEE queue full, waiting')
      time.sleep(300)
      available_slots = get_gee_queue()
      
    if available_slots>0:
    
      if names[i] in modis_dates:
       
        image = combInputs(names[i], daymet_inputs, smColl1, smColl2, tsoil, NLDAS, tavg, tmin, prcp, clay)
        image = image.clip(roi.geometry())
        print()
        print('exporting {} of {}: {}'.format(i+1, count, names[i]))
        task = ee.batch.Export.image.toCloudStorage(
                  image=image,
                  scale=30,
                  crs='EPSG:4326',
                  region=roi.geometry(),
                  bucket=bucket_name,
                  fileNamePrefix='{}covariates_{}'.format(out_directory_path, names[i]))
                
        task.start()
        available_slots-=1
    
  return
  
start_date = "2002-01-01"
end_date = "2023-01-01"
roi_asset_path = 'projects/rangelands-explo-1571664594580/assets/Ameriflux_RS/A32/A32'
#roi_asset_path = 'projects/rangelands-explo-1571664594580/assets/Ranches/HLD/G1'
#landcover_dir = '/HLD/G1/landcover/'
#export_landcover(bucket_name, roi_asset_path, landcover_dir)
#out_dir = 'Ranch_Runs/HLD/G1/covariates/'
#out_dir= 'Ameriflux_sites/A32_starfm/covariates/'
#modis_dir= 'Ameriflux_sites/A32_starfm/modis_test_v2/'
#modis_dir= 'Ranch_Runs/HLD/G1/modis/'
#covariate_tuple = get_datasets(covariate_mapping, roi_asset_path, start_date, end_date)
#export_daymet(covariate_tuple, bucket_name, roi_asset_path, modis_dir, out_dir)