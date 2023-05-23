import ee
import utils
import time 

utils.authorize()
utils.print_root_assets()

roi_asset_path = 'projects/rangelands-explo-1571664594580/assets/Ameriflux_RS/Rws/Rws'

#start_date = '2002-01-01'
#end_date = '2023-01-01'
start_date = '2002-01-01'
end_date = '2003-01-01'

LS8_collection = 'LANDSAT/LC08/C02/T1_L2'
LS7_collection = 'LANDSAT/LE07/C02/T1_L2'
LS5_collection = 'LANDSAT/LT05/C02/T1_L2'

#LS8_collection = 'LANDSAT/LC08/C01/T1_SR'
#LS7_collection = 'LANDSAT/LE07/C01/T1_SR'
#LS5_collection = 'LANDSAT/LT05/C01/T1_SR'

MODIS_collection = 'MODIS/061/MCD43A4'
COMMON_BAND_NAMES = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'QA_PIXEL', 'QA_RADSAT']
#COMMON_BAND_NAMES = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'pixel_qa']

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
  
def landsat_mask(img):
  """Bit masking for Landsat to remove clouds, cloud shadows, and snow

  Args:
    img (ee.Image): Landsat image.

  Returns:
     ee.Image
  """ 
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
         qa.bitwiseAnd(SnowMask).eq(0)).And(  
         qa.bitwiseAnd(WaterMask).eq(0))
         
  buffered_mask = mask.focal_min(radius= 15, units= 'pixels')
  #buffered_mask = mask 
  qa_radsat = img.select('QA_RADSAT')                                      
  mask_radsat = qa_radsat.bitwiseAnd(DilatedCloud).eq(0).And(
         qa_radsat.bitwiseAnd(Cirrus).eq(0)).And(
         qa_radsat.bitwiseAnd(CloudBitMask).eq(0)).And(
         qa_radsat.bitwiseAnd(CloudShadowBitMask).eq(0)).And(  
         qa_radsat.bitwiseAnd(SnowMask).eq(0))
        
  return img.updateMask(buffered_mask).updateMask(mask_radsat)
'''  
def landsat_mask(img):
  CloudShadowBitMask = 1 << 3
  CloudBitMask = 1 << 5
  SnowConMask = 1 << 4
  CloudConf = 1 << 7
  Cirrus = 1 << 9
  qa = img.select('pixel_qa')                                      
  mask = qa.bitwiseAnd(CloudShadowBitMask).eq(0).And( 
                qa.bitwiseAnd(CloudConf).eq(0)).And(
                qa.bitwiseAnd(Cirrus).eq(0)).And(
               qa.bitwiseAnd(CloudBitMask).eq(0)).And(
                 qa.bitwiseAnd(SnowConMask).eq(0))                    
  return img.updateMask(mask)
'''  
def prepOli(img):
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
  #img=img.select(['B2','B3','B4','B5','B6','B7','pixel_qa'], COMMON_BAND_NAMES) 
  orig = img
  img = landsat_mask(img)
  
  return ee.Image(img.copyProperties(orig, orig.propertyNames())).select(['blue', 'green', 'red', 'nir', 'swir1'])
 
def prepEtm(img):
  """Prepares Landsat OLI images for merging with ETM Image Collection. 
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
  img = landsat_mask(img);
  img = etmToOli(img);
  
  return ee.Image(img.copyProperties(orig, orig.propertyNames())).select(['blue', 'green', 'red', 'nir', 'swir1'])

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
  #Prepare them for merging by first running the transformation functions.
  oliColR = oliCol.map(prepOli)
  etmColR = etmCol.map(prepEtm)
  tmColR = tmCol.map(prepEtm)
  
  # Merge the collections.
  mercollection = oliColR.merge(etmColR).merge(tmColR).sort('system:time_start')

  return mercollection  

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
  MODIS_export = MODIS.select(['Nadir_Reflectance_Band3','Nadir_Reflectance_Band4','Nadir_Reflectance_Band1', 'Nadir_Reflectance_Band2', 'Nadir_Reflectance_Band6'],['blue', 'green', 'red', 'nir', 'swir1'])
  
  names = MODIS.sort('system:time_start').aggregate_array('system:index').getInfo()

  return MODIS_export

def export_landsat_collection(collection, roi_asset_path, bucket_name, directory_path):
  """Exports landsat image collection to google storage bucket

  Args:
    collection (ee.ImageCollection): Landsat image collection.
    roi_asset_path (string): path to roi Earth Engine Asset
    bucket_name (string): google storage bucket name
    directory_path (string): path to directory to store

  Returns:
     None
  """ 
  roi=ee.FeatureCollection(roi_asset_path)
  count = int(collection.size().getInfo())
  names = collection.aggregate_array('system:index').getInfo()
  
  for i in range(0, count):
  
    image = collection.filter(ee.Filter.eq('system:index', names[i])).first()
    print()
    print('exporting {} of {}: {}'.format(i+1, count, names[i]))
    task = ee.batch.Export.image.toCloudStorage(
              image=image,
              scale=30,
              crs='EPSG:4326',
              region=roi.geometry(),
              bucket=bucket_name,
              fileNamePrefix='{}{}'.format(directory_path, names[i]))
            
    task.start()
  
  return
    
def export_modis_collection(collection, roi_asset_path, bucket_name, directory_path):
  """Exports MODIS image collection to google storage bucket

  Args:
    collection (ee.ImageCollection): MODIS image collection.
    roi_asset_path (string): path to roi Earth Engine Asset
    bucket_name (string): google storage bucket name
    directory_path (string): path to directory to store

  Returns:
     None
  """ 
  
  roi=ee.FeatureCollection(roi_asset_path)
  count = int(collection.size().getInfo())
  names = collection.aggregate_array('system:index').getInfo()
  
  for i in range(0, count):
  
    image = collection.filter(ee.Filter.eq('system:index', names[i])).first()
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
              fileNamePrefix='{}MCD43A4_{}'.format(directory_path, names[i]))
            
    task.start()
    
  return
  
  
#modis, ids = get_MODIS(MODIS_collection, roi_asset_path, start_date, end_date)
#export_modis_collection(modis, ids, roi_asset_path, 'rangelands', 'Ameriflux_sites/Rws_starfm/modis_test/')

landsat = get_landsat(LS8_collection, LS7_collection, LS5_collection, roi_asset_path, start_date, end_date)
export_landsat_collection(landsat, roi_asset_path, 'rangelands', 'Ameriflux_sites/Rws_starfm/landsat_test/')


       
