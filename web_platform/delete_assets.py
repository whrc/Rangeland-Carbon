from google.auth.transport.requests import AuthorizedSession
import ee
from google.auth import compute_engine

project_folder = 'rangelands-explo-1571664594580'

#ee.data.deleteAsset(f'projects/{project_folder}/assets/' + asset_path)

credentials = compute_engine.Credentials()
ee.Initialize(credentials, project='rangelands-explo-1571664594580')
session = AuthorizedSession(credentials)

ameriflux_folders = ee.data.listAssets(f'projects/{project_folder}/assets/Ameriflux_RS')['assets']

for ameriflux_folder in ameriflux_folders:

  folder = ameriflux_folder['name']
  
  for item in ee.data.listAssets(folder)['assets']:
  
    if item['type'] == 'IMAGE_COLLECTION':
      
      for image in ee.data.listAssets(item['id'])['assets']:
        
        ee.data.deleteAsset(image['id'])
  

