
from google.cloud import storage
import ee

def print_root_assets():
  folders = ee.data.getAssetRoots()
  for folder in folders:
    print(folder['id'] + ', type: ' + folder['type'])
            
def authorize():
  service_account = 'rctm-07dea6fa718fdd0921f124263@rangelands-explo-1571664594580.iam.gserviceaccount.com'
  credentials = ee.ServiceAccountCredentials(service_account, 'gee_key.json')
  
  # Initialize the library.
  ee.Initialize(credentials, project='rangelands-explo-1571664594580')
  print('Earth Engine Initialized')
  
def gs_listdir(directory, bucket):
  """returns list of all items in a directory within a google cloud bucket

  Args:
    directory (string): path to directory to list (should not include bucket name, should not start with forward slash)
    bucket (storage.bucket): google cloud bucket

  Returns:
     list
  """ 
  listdir = [blob.name for blob in list(bucket.list_blobs(prefix=directory))]
  return listdir

def gs_write_blob(local_filepath, gs_filepath, bucket):
  """copies file from disk to location in google cloud bucket

  Args:
    local_filepath (string): path to file to copy
    gs_filepath (string): path to location in google bucket (should not include bucket name, should not start with forward slash)
    bucket (storage.bucket): google cloud bucket

  Returns:
     None
  """ 
  blob = bucket.blob(gs_filepath)
  blob.upload_from_filename(local_filepath)



