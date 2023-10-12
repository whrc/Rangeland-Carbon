import ee
import utils
import time 
import numpy as np
from functools import partial
import preprocess_starfm_imagery as psi
from google.cloud import storage
from datetime import datetime as dt
import sys
import argparse
import os
import pandas as pd

utils.authorize()
utils.print_root_assets()

#def verify_ameriflux_rois(path_to_footprint_list):
def verify_ameriflux_rois(path_to_footprint_list = 'res/spatial_sites.txt'):
  
  with open(path_to_footprint_list) as f:
  
      sites = [line.rstrip('\n') for line in f]
      print('sites: {}'.format(sites))
  
  total_area=0
  
  for site in sites:
    if site=='Kon':
      continue
    #roi_asset_path = 'projects/rangelands-explo-1571664594580/assets/Ameriflux_RS/{}/{}'.format(site, site)
    roi_asset_path = 'projects/rangelands-explo-1571664594580/assets/Shapefiles/{}'.format(site, site)
    roi=ee.FeatureCollection(roi_asset_path)
    area = roi.geometry().area().getInfo()
    print(area)
    total_area+=area
    
  print(f'total area: {total_area}')

#verify_ameriflux_rois(path_to_footprint_list = 'res/sites.txt')    
verify_ameriflux_rois(path_to_footprint_list='res/spatial_sites.txt')
