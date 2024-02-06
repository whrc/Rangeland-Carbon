import ee
import time 
import numpy as np
from functools import partial
import starfm_preprocessing as psi
from google.cloud import storage
from datetime import datetime as dt
import sys
import argparse
import os
import pandas as pd
import sys
sys.path.insert(1, '../utils')
import utils

utils.authorize(high_endpoint=False)
utils.print_root_assets()

#roi_asset_path = 'projects/rangelands-explo-1571664594580/assets/Ameriflux_RS/Rws/Rws'

start_date = '2002-01-01'
end_date = '2023-01-01'

LS8_collection = 'LANDSAT/LC08/C02/T1_L2'
LS7_collection = 'LANDSAT/LE07/C02/T1_L2'
LS5_collection = 'LANDSAT/LT05/C02/T1_L2'

MODIS_collection = 'MODIS/061/MCD43A4'
COMMON_BAND_NAMES = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'QA_PIXEL', 'QA_RADSAT']



#slope and intercept citation: Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity, Remote Sensing of Environment, 185, 57-70.(http://dx.doi.org/10.1016/j.rse.2015.12.024); Table 2 - reduced major axis (RMA) regression coefficients
#Table 2 OLS regression coefficients. Band-respective coefficients are defined in the following dictionary with slope (slopes) and intercept (itcps) image constants
coefficients = {'itcps': ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254])
                            .multiply(10000),
  'slopes': ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937])
}

def etmToOli(img):
  """Applies linear correction to Landsat ETM bands to match Landsat OLI

  Args:
    img (ee.Image): Landsat ETM image.

  Returns:
     ee.Image
  """ 
  corrected_bands = img.select(['blue', 'green', 'red', 'nir', 'swir1']).multiply(coefficients['slopes']).add(coefficients['itcps']).round().toShort()
  return img.addBands(corrected_bands, None, True)

    
def topoCorrect(image):
  """Applies topographic correction to Landsat images

  Args:
    img (ee.Image): Landsat ETM image.

  Returns:
     ee.Image
  """ 
  
  elevation = ee.Image('USGS/NED')
  elevation = elevation.reduceResolution(
      reducer= ee.Reducer.mean(),
      maxPixels= 1024
    ).reproject(
      crs=  image.projection()
    )
  
  slope = ee.Terrain.slope(elevation).divide(180).multiply(np.pi);
  aspect = ee.Terrain.aspect(elevation).divide(180).multiply(np.pi);
  
  #generate slope and aspect rasters
  slope = slope.updateMask(image.select(['blue']).mask())
  aspect = aspect.updateMask(image.select(['blue']).mask())
  elevation = elevation.updateMask(image.select(['blue']).mask())
  solar_z = ee.Number(90).subtract(image.getNumber("SUN_ELEVATION")).divide(180).multiply(np.pi)
  solar_a = image.getNumber("SUN_AZIMUTH").divide(180).multiply(np.pi)
  incidence_angle = slope.expression(
    'acos(cos(slope) * cos( solar_z ) + (sin( solar_z ) * sin( slope ) * cos(solar_a - aspect)))', {
    'slope' : slope,
    'solar_z' : solar_z,
    'solar_a' : solar_a,
    'aspect' : aspect
  })
  
  #incidence_mask = incidence_angle.select('slope').lte(1.22);
  cos_incidence_arr = incidence_angle.cos()
  DN_arr = image
  combinedImage = cos_incidence_arr.addBands(DN_arr)
  
  def corr(bandName, acc):

    bandName = ee.String(bandName) #Must cast from ee.Element to actual type
    linearFit = combinedImage.select(['slope', bandName]).reduceRegion(
      reducer =  ee.Reducer.linearFit(),
      geometry = image.geometry(),
      scale = 30,
      )
   
    #store emperical parameters
    m = linearFit.getNumber('scale')
    b = linearFit.getNumber('offset')
    c = b.divide(m);
        
    #correct surface reflectance
    image_corr = image.select([bandName]).expression(
    'sr * ((cos(solar_z) + c) / (cos(incidence_angle) + c))', {
    'sr' : image.select([bandName]),
    'solar_z' : solar_z,
    'incidence_angle' : incidence_angle,
    'c' : c
    })

    return ee.Image(acc).addBands(image_corr, None, True)
  
  bandNames = ee.List(image.select(['blue', 'green', 'red', 'nir', 'swir1']).bandNames())

  correctedImage = ee.Image(bandNames.iterate(corr, image))  
  
  return image.addBands(correctedImage, None, True).toInt16()

def gapfill_landsat(image, reference_collection, kernelSize=50): #kernel size 50-100 good for clouds/etc
  """Applies gap-filling to Landsat images to fill in scanline errors and cloud gaps. Applies linear regression between reference image (without gaps) and target image (with gaps). 
  Reference image is based on a median composite image.

  Args:
    img (ee.Image): Landsat ETM image.
    reference_collection (ee.ImageCollection): Landsat image collection to base gap-filling on.
    kernel_size (int), optional: Window size for linear regression. Larger kernel means larger gaps are able to be filled, but the regression may have reduced accuracy.

  Returns:
     ee.Image
  """ 
  
  kernel = ee.Kernel.square(kernelSize * 30, 'meters', False)

  start = image.date().advance(-1, 'year')
  end = image.date().advance(1, 'year')
  im_month = image.date().get('month')
  min_month=im_month.subtract(1)
  max_month=im_month.add(1)
  fill = reference_collection.filter(ee.Filter.calendarRange(min_month, max_month, 'month')).filterMetadata('CLOUD_COVER', 'less_than', 5)
  #fill = reference_collection.filterDate(start, end).filterMetadata('CLOUD_COVER', 'less_than', 5)
  #fill = reference_collection.filterDate(start, end).filterMetadata('CloudSnowMaskedPercent', 'less_than', 5)
  fill = fill.median()#.clip(roi)
  #fill = reference_collection.median().clip(roi)

  regress = fill.addBands(image)

  #extract band pairs from fill image (e.g. 'blue') and input image (e.g. 'blue_1')
  #blue_bands = regress.select(['blue','blue_1'])
  #green_bands = regress.select(['green','green_1'])
  red_bands = regress.select(['red','red_1'])
  nir_bands = regress.select(['nir','nir_1'])
  #swir1_bands = regress.select(['swir1','swir1_1'])

  #for each band, fit linear regression for each kernel window, with fill values as x and input image values as y
  #next, apply regression coefficients to fill image

  #blue_fit = blue_bands.reduceNeighborhood(ee.Reducer.linearFit(), kernel, None, False); #linear regression coefficients are calculated 
  #fill = fill.addBands(fill.select(['blue']).multiply(blue_fit.select('scale')).add(blue_fit.select('offset')), ['blue'], True)
  #green_fit = green_bands.reduceNeighborhood(ee.Reducer.linearFit(), kernel, None, False); #linear regression coefficients are calculated 
  #fill = fill.addBands(fill.select(['green']).multiply(green_fit.select('scale')).add(green_fit.select('offset')), ['green'], True)

  red_fit = red_bands.reduceNeighborhood(ee.Reducer.linearFit(), kernel, None, False); #linear regression coefficients are calculated
  fill = fill.addBands(fill.select(['red']).multiply(red_fit.select('scale')).add(red_fit.select('offset')), ['red'], True) 

  nir_fit = nir_bands.reduceNeighborhood(ee.Reducer.linearFit(), kernel, None, False); #linear regression coefficients are calculated
  fill = fill.addBands(fill.select(['nir']).multiply(nir_fit.select('scale')).add(nir_fit.select('offset')), ['nir'], True) 

  #swir1_fit = swir1_bands.reduceNeighborhood(ee.Reducer.linearFit(), kernel, None, False); #linear regression coefficients are calculated
  #fill = fill.addBands(fill.select(['swir1']).multiply(swir1_fit.select('scale')).add(swir1_fit.select('offset')), ['swir1'], True) 

  return image.unmask(fill).toInt16()
  
def landsat_mask(img, roi_asset_path):
  """Bit masking for Landsat to remove clouds, cloud shadows, and snow

  Args:
    img (ee.Image): Landsat image.

  Returns:
     ee.Image
  """ 
  roi=ee.FeatureCollection(roi_asset_path)
  geometry = roi.geometry()
  
  DilatedCloud = 1 << 1
  Cirrus = 1 << 2
  CloudBitMask = 1 << 3
  CloudShadowBitMask = 1 << 4
  SnowMask = 1 << 5
  WaterMask = 1 << 7
  CloudShadowConfFirstBit = 1 << 10
  CloudShadowConfFirstBit = 1 << 11
  
  qa_b1 = 1 << 0
  qa_b2 = 1 << 1
  qa_b3 = 1 << 2
  qa_b4 = 1 << 3
  qa_b5 = 1 << 4

  qa = img.select('QA_PIXEL')                                      
  mask = qa.bitwiseAnd(DilatedCloud).eq(0).And(
         qa.bitwiseAnd(Cirrus).eq(0)).And(
         qa.bitwiseAnd(CloudBitMask).eq(0)).And(
         qa.bitwiseAnd(CloudShadowBitMask).eq(0)).And(  
         qa.bitwiseAnd(SnowMask).eq(0))
         
  buffered_mask = mask.focal_min(radius= 15, units= 'pixels')
  
  ndsi = (img.select('green').subtract(img.select('swir1'))).divide(img.select('green').add(img.select('swir1')))
  snow_cover_mask = ndsi.lte(0.4);
  buffered_mask = buffered_mask.updateMask(snow_cover_mask)
  
  water_mask = qa.bitwiseAnd(WaterMask).eq(0).focal_min(radius= 2, units= 'pixels')
  
  #qa_radsat = img.select('QA_RADSAT')                                      
  #mask_radsat = qa_radsat.bitwiseAnd(DilatedCloud).eq(0).And(
         #qa_radsat.bitwiseAnd(Cirrus).eq(0)).And(
         #qa_radsat.bitwiseAnd(CloudBitMask).eq(0)).And(
         #qa_radsat.bitwiseAnd(CloudShadowBitMask).eq(0))#.And(  
         #qa_radsat.bitwiseAnd(SnowMask).eq(0))
  
  maskedMask = buffered_mask.updateMask(buffered_mask)
  
  maskedCount = maskedMask.select(['QA_PIXEL']) \
        .reduceRegion(reducer=ee.Reducer.count(),
                      geometry=geometry,
                      scale=ee.Number(30),
                      maxPixels=ee.Number(4e10))

   #count the total number of pixels
  origCount = img.select(['blue']) \
        .reduceRegion(reducer=ee.Reducer.count(),
                      geometry=geometry,
                      scale=ee.Number(30),
                      maxPixels=ee.Number(4e10))

   #calculate the percent of masked pixels
  percent = ee.Number(origCount.get('blue'))\
        .subtract(maskedCount.get('QA_PIXEL'))\
        .divide(origCount.get('blue'))\
        .multiply(100)\
        .round()

   # Return the masked image with new property and time stamp
  return img.updateMask(buffered_mask).updateMask(water_mask).set('CloudSnowMaskedPercent', percent)
        
  #return img.updateMask(buffered_mask).updateMask(mask_radsat)
def landsat_water_snow_mask(img, roi_asset_path):
  """Bit masking for Landsat to remove water, and snow

  Args:
    img (ee.Image): Landsat image.

  Returns:
     ee.Image
  """ 
  roi=ee.FeatureCollection(roi_asset_path)
  geometry = roi.geometry()
  
  SnowMask = 1 << 5
  WaterMask = 1 << 7

  qa = img.select('QA_PIXEL')                                      
  mask = qa.bitwiseAnd(WaterMask).eq(0)
  
  #ndsi = (img.select('green').subtract(img.select('swir1'))).divide(img.select('green').add(img.select('swir1')))
  #snow_cover_mask = ndsi.lte(0.4);
  #mask = mask.updateMask(snow_cover_mask)
  
  #qa_radsat = img.select('QA_RADSAT')                                      
 # mask_radsat = qa_radsat.bitwiseAnd(SnowMask).eq(0)

   # Return the masked image with new property and time stamp
  return img.updateMask(mask).select(['red', 'nir'])
        
  #return img.updateMask(buffered_mask).updateMask(mask_radsat)

def prepOli(img, roi_asset_path):
  """Prepares Landsat OLI images for merging with ETM Image Collection. 
     Rescales to scale factor of 10000, maps bands to common names, and masks image

  Args:
    img (ee.Image): Landsat OLI image.

  Returns:
     ee.Image
  """ 
  #apply collection 2 scaling factors, then rescale to 10000 to match MODIS scale factor
  #opticalBands = img.select('SR_B.').multiply(0.0000275).add(-0.2)
  #opticalBands = opticalBands.multiply(10000).toInt16()
  opticalBands = img.select('SR_B.').multiply(0.275).add(-2000).toInt16()
  img = img.addBands(opticalBands, None, True)
  
  #filter out thermal bands
  img=img.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','QA_PIXEL','QA_RADSAT'], COMMON_BAND_NAMES)
  orig = img
  img = landsat_mask(img, roi_asset_path)

  #return ee.Image(img.copyProperties(orig, orig.propertyNames())).select(['blue', 'green', 'red', 'nir', 'swir1'])
  return ee.Image(img.copyProperties(orig, orig.propertyNames())).select(['green', 'red', 'nir', 'swir1', 'QA_PIXEL','QA_RADSAT'])
 
def prepEtm(img, roi_asset_path):
  """Prepares Landsat ETM images for merging with OLI Image Collection. 
     Rescales to scale factor of 10000, maps bands to common names, masks image, and harmonizes with OLI

  Args:
    img (ee.Image): Landsat ETM image.

  Returns:
     ee.Image
  """ 
  #apply collection 2 scaling factors, then rescale to 10000 to match MODIS scale factor
  #opticalBands = img.select('SR_B.').multiply(0.0000275).add(-0.2)
  #opticalBands = opticalBands.multiply(10000).toInt16()
  opticalBands = img.select('SR_B.').multiply(0.275).add(-2000).toInt16() 
  img = img.addBands(opticalBands, None, True)
  
  #filter out thermal bands
  img=img.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','QA_PIXEL', 'QA_RADSAT'], COMMON_BAND_NAMES)
  #img=img.select(['B1','B2','B3','B4','B5','B7','pixel_qa'], COMMON_BAND_NAMES)  
  
  orig = img;
  img = landsat_mask(img, roi_asset_path);
  #img = etmToOli(img);
  
  #return ee.Image(img.copyProperties(orig, orig.propertyNames())).select(['blue', 'green', 'red', 'nir', 'swir1'])
  return ee.Image(img.copyProperties(orig, orig.propertyNames())).select(['green', 'red', 'nir', 'swir1', 'QA_PIXEL','QA_RADSAT'])

def filterColl (SOURCE, geometry, startdate, enddate):
  """Retrieves Image Collection from Earth Engine Catalog

  Args:
    source (string): path to Earth Engine collection.
    geometry (table): geometry to filter images.
    startdate (string): start date for filtering in format YYYY-MM-DD
    enddate (string): end date for filtering in format YYYY-MM-DD

  Returns:
     ee.ImageImageCollection
  """ 

  coll = ee.ImageCollection(SOURCE).filter(ee.Filter.bounds(geometry)).filterDate(startdate, enddate)
  return(coll)

def clip(image, roi_asset_path):
  roi=ee.FeatureCollection(roi_asset_path)
  return image.clip(roi.geometry()).toInt16()

def get_landsat(LS8_collection, LS7_collection, LS5_collection, roi_asset_path, start_date, end_date):
  """Retrieves Landsat Image Collection from Earth Engine Catalog

  Args:
    LS8_collection (string): path to Landsat 8 Earth Engine collection.
    LS7_collection (string): path to Landsat 7 Earth Engine collection.
    LS5_collection (string): path to Landsat 5 Earth Engine collection.
    roi_asset_path (string): path to roi asset.
    startdate (string): start date for filtering in format YYYY-MM-DD
    enddate (string): end date for filtering in format YYYY-MM-DD

  Returns:
     ee.ImageImageCollection
  """ 
  
  roi=ee.FeatureCollection(roi_asset_path)
  
  oliCol = filterColl(LS8_collection, roi, start_date, end_date)
  etmCol = filterColl(LS7_collection, roi, start_date, end_date)
  tmCol = filterColl(LS5_collection, roi, start_date, end_date)
  
  clip_fill = partial(clip, roi_asset_path=roi_asset_path)
  oli_fill = partial(prepOli, roi_asset_path=roi_asset_path)
  etm_fill = partial(prepEtm, roi_asset_path=roi_asset_path)
  tm_fill = partial(prepEtm, roi_asset_path=roi_asset_path)
  
  snow_water_mask_fill = partial(landsat_water_snow_mask, roi_asset_path=roi_asset_path)
  
  #Prepare them for merging by first running the transformation functions.
  
  oliColR = oliCol.map(oli_fill)
  oliColR = oliColR.map(clip_fill)
  
  etmColR = etmCol.map(etm_fill)
  etmColR = etmColR.map(clip_fill)
  
  tmColR = tmCol.map(tm_fill)
  tmColR = tmColR.map(clip_fill)
  
  # Merge the collections.
  mercollection = oliColR.merge(etmColR).merge(tmColR).sort('system:time_start').filter(ee.Filter.lt('CloudSnowMaskedPercent',60))
  #mercollection = mercollection.map(topoCorrect)
  
  fill_function = partial(gapfill_landsat, reference_collection=mercollection, kernelSize=50)
  gap_filled_collection = mercollection.map(fill_function)
  
  gap_filled_collection = gap_filled_collection.map(snow_water_mask_fill)

  #return gap_filled_collection
  return gap_filled_collection

def MODIS_mask(img):
  """Bit masking for MODIS to remove low quality retrievals

  Args:
    img (ee.Image): MODIS image.

  Returns:
     ee.Image
  """ 
  QualityMask = 1 << 0

  b1 = img.select('BRDF_Albedo_Band_Mandatory_Quality_Band1')                                      
  b1_mask = b1.bitwiseAnd(QualityMask).eq(0)

  b2 = img.select('BRDF_Albedo_Band_Mandatory_Quality_Band2')                                      
  b2_mask = b2.bitwiseAnd(QualityMask).eq(0)

  b3 = img.select('BRDF_Albedo_Band_Mandatory_Quality_Band3')                                      
  b3_mask = b3.bitwiseAnd(QualityMask).eq(0)

  b4 = img.select('BRDF_Albedo_Band_Mandatory_Quality_Band4')                                      
  b4_mask = b4.bitwiseAnd(QualityMask).eq(0)

  b5 = img.select('BRDF_Albedo_Band_Mandatory_Quality_Band5')                                      
  b5_mask = b5.bitwiseAnd(QualityMask).eq(0)

  b6 = img.select('BRDF_Albedo_Band_Mandatory_Quality_Band6')                                      
  b6_mask = b6.bitwiseAnd(QualityMask).eq(0)

  b7 = img.select('BRDF_Albedo_Band_Mandatory_Quality_Band7')                                      
  b7_mask = b7.bitwiseAnd(QualityMask).eq(0)

  img = img.select(['Nadir_Reflectance_Band1','Nadir_Reflectance_Band2','Nadir_Reflectance_Band3',
                'Nadir_Reflectance_Band4','Nadir_Reflectance_Band5','Nadir_Reflectance_Band6',
                'Nadir_Reflectance_Band7']).updateMask(b1_mask).updateMask(b2_mask).updateMask(b3_mask).updateMask(b4_mask).updateMask(b5_mask).updateMask(b6_mask).updateMask(b7_mask)
                
  ndsi = (img.select('Nadir_Reflectance_Band4').subtract(img.select('Nadir_Reflectance_Band6'))).divide(img.select('Nadir_Reflectance_Band4').add(img.select('Nadir_Reflectance_Band6')))
  snow_cover_mask = ndsi.lte(0.4);
  img = img.updateMask(snow_cover_mask)
  
  return img
  

def get_MODIS(MODIS_collection, roi_asset_path, start_date, end_date):
  """Retrieves MODIS Image Collection from Earth Engine Catalog

  Args:
    LS8_collection (string): path to Landsat 8 Earth Engine collection.
    LS7_collection (string): path to Landsat 7 Earth Engine collection.
    LS5_collection (string): path to Landsat 5 Earth Engine collection.
    roi_asset_path (string): path to roi asset.
    startdate (string): start date for filtering in format YYYY-MM-DD
    enddate (string): end date for filtering in format YYYY-MM-DD

  Returns:
     ee.ImageImageCollection
  """ 
  roi=ee.FeatureCollection(roi_asset_path)
  MODIS = filterColl(MODIS_collection, roi, start_date, end_date)
  MODIS = MODIS.map(MODIS_mask)
  #MODIS_export = MODIS.select(['Nadir_Reflectance_Band3','Nadir_Reflectance_Band4','Nadir_Reflectance_Band1', 'Nadir_Reflectance_Band2', 'Nadir_Reflectance_Band6'],['blue', 'green', 'red', 'nir', 'swir1'])
  
  MODIS_export = MODIS.select(['Nadir_Reflectance_Band1', 'Nadir_Reflectance_Band2'],['red', 'nir'])
  
  clip_fill = partial(clip, roi_asset_path=roi_asset_path)
  
  MODIS_export = MODIS_export.map(clip_fill)
  
  names = MODIS.sort('system:time_start').aggregate_array('system:index').getInfo()

  return MODIS_export

def get_gee_queue():

  ops = ee.data.listOperations()
  ids = [ op['name'] for op in ops ]
  state = [ op['metadata']['state'] for op in ops ]
  op_df = pd.DataFrame({'name': ids, 'state': state})
  queue_len = len(op_df[(op_df['state']=='PENDING') | (op_df['state']=='RUNNING')])
  available_slots = 3000-queue_len
  
  return available_slots

def export_landsat_collection(collection, roi_asset_path, bucket_name, directory_path, overwrite=True):
  """Exports landsat image collection to google storage bucket

  Args:
    collection (ee.ImageCollection): Landsat image collection.
    roi_asset_path (string): path to roi Earth Engine Asset
    bucket_name (string): google storage bucket name
    directory_path (string): path to directory to store

  Returns:
     None
  """ 
  storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
  bucket = storage_client.get_bucket(bucket_name)
  
  roi=ee.FeatureCollection(roi_asset_path)
  
  names_gee = collection.aggregate_array('system:index')
  names=names_gee.getInfo()
  count = len(names)
  export_names = ['_'.join(name.split('_')[-3:]) for name in names]
  if len(names)==0:
    return ([], [])
    
  df_im_paths = pd.DataFrame({'name':names})
  dates = pd.to_datetime(df_im_paths['name'].str[-8:], format='%Y%m%d').dt.strftime('%Y_%m_%d').to_list()
  available_slots = get_gee_queue()
  task_ids=[]
  
  for i in range(0, count):
  
    exists = storage.Blob(bucket=bucket, name='{}{}'.format(directory_path, export_names[i]+'.tif')).exists(storage_client)

    if exists:
      if overwrite==False:
        print('{} already exists! Skipping'.format(export_names[i]))
        continue
        
    if available_slots<=0:
      print('GEE queue full, waiting')
      time.sleep(300)
      available_slots = get_gee_queue()
      
    if available_slots>0:
      image = collection.filter(ee.Filter.eq('system:index', names[i])).first()
      print()
      print('submitting landsat export task for {}'.format(export_names[i]))
      task = ee.batch.Export.image.toCloudStorage(
                image=image,
                scale=30,
                crs='EPSG:4326',
                region=roi.geometry(),
                bucket=bucket_name,
                fileNamePrefix='{}{}'.format(directory_path, export_names[i]))
              
      task.start()
      task_ids.append(task.id)
      available_slots-=1
  
  return dates, task_ids

def export_modis_collection(modis_collection, roi_asset_path, bucket_name, out_directory_path, landsat_dates=[], landsat_dir=None, overwrite=False):
  """Exports MODIS image collection to google storage bucket. Exports images matching landsat acquisitions and every five days.

  Args:
    collection (ee.ImageCollection): MODIS image collection.
    roi_asset_path (string): path to roi Earth Engine Asset
    bucket_name (string): google storage bucket name
    out_directory_path (string): path to directory to store

  Returns:
     None
  """ 
  storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
  bucket = storage_client.get_bucket(bucket_name)
  
  roi=ee.FeatureCollection(roi_asset_path)
  
  names_gee = modis_collection.aggregate_array('system:index')
  names=names_gee.getInfo()
  count = len(names)
  
  if len(names)==0:
    return []
  
  task_ids=[]
  
  if landsat_dir is not None:
    landsat_df = psi.get_landsat_date_df(landsat_dir, bucket)
    landsat_dates = landsat_df['im_date'].dt.strftime('%Y_%m_%d').to_list()
  
  available_slots = get_gee_queue()
  landsat_dates
  for i in range(0, count):
  
    exists = storage.Blob(bucket=bucket, name='{}MCD43A4_{}.tif'.format(out_directory_path, names[i])).exists(storage_client)

    if exists:
      if overwrite==False:
        print('MCD43A4_{} already exists! Skipping'.format(names[i]))
        continue
    
    if names[i] in landsat_dates or i%5==0:
      
      if available_slots<=0:
        print('GEE queue full, waiting')
        time.sleep(300)
        available_slots = get_gee_queue()
        
      if available_slots>0:
        image = modis_collection.filter(ee.Filter.eq('system:index', names[i])).first()
        print()
        print('submitting MODIS export task for {}'.format(names[i]))
        task = ee.batch.Export.image.toCloudStorage(
                  image=image,
                  scale=30,
                  #crs=image.projection().crs(),
                  #crsTransform = image.projection().getInfo()['transform'],
                  crs='EPSG:4326',
                  region=roi.geometry(),
                  bucket=bucket_name,
                  fileNamePrefix='{}MCD43A4_{}'.format(out_directory_path, names[i]))
        task.start()
        task_ids.append(task.id)
        available_slots-=1
        
  return task_ids
  


if __name__ == "__main__":
  parser=argparse.ArgumentParser()
  parser.add_argument("--which", help="whether to fetch and download landsat, modis, or both")
  parser.add_argument("--out_directory", help="path to export imagery subdirectories")
  parser.add_argument("--bucket_name", help="bucket name")
  parser.add_argument("--roi_asset_path", help="path to roi Earth Engine Asset")

  args=parser.parse_args()
  
  if args.which=='landsat':
    landsat = get_landsat(LS8_collection, LS7_collection, LS5_collection, args.roi_asset_path, start_date, end_date)
    export_landsat_collection(landsat, args.roi_asset_path, args.bucket_name, os.path.join(args.out_directory, 'landsat_v2/'))
    
  if args.which=='modis':
    modis = get_MODIS(MODIS_collection, args.roi_asset_path, start_date, end_date)
    export_modis_collection(modis, args.roi_asset_path, args.bucket_name, os.path.join(args.out_directory, 'modis/'), landsat_dir = os.path.join(args.out_directory, 'landsat_v2/'))
    
  if args.which=='both':
    
    landsat = get_landsat(LS8_collection, LS7_collection, LS5_collection, args.roi_asset_path, start_date, end_date)
    landsat_dates, landsat_task_ids = export_landsat_collection(landsat, args.roi_asset_path, args.bucket_name, os.path.join(args.out_directory, 'landsat_v2/'))
    
    modis = get_MODIS(MODIS_collection, args.roi_asset_path, start_date, end_date)
    modis_task_ids = export_modis_collection(modis, args.roi_asset_path, args.bucket_name, os.path.join(args.out_directory, 'modis/'), landsat_dates=landsat_dates)  
    
        

#python GEE_fetch_starfm_imagery.py --which='landsat' --out_directory='Ameriflux_sites/Cop_starfm/' --bucket_name='rangelands' --roi_asset_path='projects/rangelands-explo-1571664594580/assets/Ameriflux_RS/Cop/Cop'

#python GEE_fetch_starfm_imagery.py --which='modis' --out_directory='Ameriflux_sites/BRG_starfm/' --bucket_name='rangelands' --roi_asset_path='projects/rangelands-explo-1571664594580/assets/Ameriflux_RS/BRG/BRG'

#modis = get_MODIS(MODIS_collection, roi_asset_path, start_date, end_date)
#export_modis_collection(modis, roi_asset_path, 'rangelands', 'Ameriflux_sites/Rws_starfm/modis_test/', 'Ameriflux_sites/Rws_starfm/landsat_test/')
  


#landsat = get_landsat(LS8_collection, LS7_collection, LS5_collection, roi_asset_path, start_date, end_date)
#export_landsat_collection(landsat, roi_asset_path, 'rangelands', 'Ameriflux_sites/Rws_starfm/landsat_test_no_topo/')


       
