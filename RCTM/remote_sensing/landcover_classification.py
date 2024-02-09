import rioxarray as rxr
from google.cloud import storage
import yaml
import numpy as np
from scipy import stats
import pandas as pd
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from matplotlib import pyplot as plt

bucket_name='rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')

#get param file for pft specifications
param_file = '/home/amullen/Rangeland-Carbon/modeling/RCTM_params.yaml'
params=[]
with open(param_file, 'r') as file:
  params = yaml.safe_load(file)
  
pft_key = {
  'grassland': ['AR1', 'AR2','ARbc', 'Aud', 'Bkg', 'BMM', 'BRG', 'Ctn', 'Dia', 'FPe', 'Hn2_3',  'IB2', 'KFS', 'KLS', 'KM2_3', 'KM4', 'Kon', 'LS1_2', 'RFW', 'Ro4','SCg', 'SdH', 'Seg', 'Shd', 'Var', 'Wjs', 'xAE', 'xCP', 'xDC', 'xKA', 'xKZ', 'xNG', 'xWD', 'A32', 'Tx1', 'Tx2', 'Snd', 'Snf'],
  'grass-shrub': ['Cop', 'LS1_2', 'SRM', 'Wkg', 'Jo1', 'Rws', 'Ses', 'xJR', 'xMB', 'xNQ', 'Wdn', 'Rwe_Rms', 'Rls','RIs', 'Rwf', 'xYE','Hn2_3'],
  'grass-tree': ['CZ1_xSJ', 'Fwf', 'Mpj', 'Ton', 'xCL'],
  'hay-pasture': ['A32', 'Tx1', 'Tx2', 'Snd', 'Snf']}
  
# Create a list of dictionaries
data = [{'site': value, 'true_pft': key} for key, values in pft_key.items() for value in values]

# Create a Pandas DataFrame
pft_key_df = pd.DataFrame(data)
  
NLCD_LUT = {'41': 'grass-tree',
            '42': 'grass-tree',
            '43': 'grass-tree',
            '51': 'grass-shrub',
            '52': 'grass-shrub',
            '71': 'grassland',
            '72': 'grassland',
            '81': 'grassland',
            '82': 'grassland',
            '90': 'grass-shrub',
            '95': 'grassland'}
            
RAP_thresholds = {'grassland':
                    {'FG':[30,100],
                     'SHR':[0,15],
                     'TRE':[0,10]},
                  'grass-shrub':
                    {'FG':[0,10],
                     'SHR':[30,100],
                     'TRE':[0,10]},
                  'grass-tree':
                    {'FG':[0,40],
                     'SHR':[0,10],
                     'TRE':[40,100]}
                  }

#open footprint list
path_to_footprint_list = '/home/amullen/Rangeland-Carbon/res/site_footprints/sites.txt'    
with open(path_to_footprint_list) as f:
  
  rois = [line.rstrip('\n') for line in f]
  sites = [roi.split('/')[-1] for roi in rois]
  NLCD_classes = []
  RAP_classes= []
  
  for i, roi, in enumerate(rois):
    print(sites[i])
    
    #open RAP
    landcover_RAP_path = 'Ameriflux_sites/{}_starfm/landcover/RAP_2019.tif'.format(sites[i])
    landcover_RAP = rxr.open_rasterio('gs://' + bucket_name + '/' + landcover_RAP_path, masked=True)
    #open NLCD
    NLDC_valid_covers = [41, 42, 43, 51, 52, 71, 72, 81, 82, 90, 95]
    landcover_NLCD_path = 'Ameriflux_sites/{}_starfm/landcover/NLCD_2019.tif'.format(sites[i])
    landcover_NLCD = rxr.open_rasterio('gs://' + bucket_name + '/' + landcover_NLCD_path, masked=True)
    landcover_NLCD[0] = landcover_NLCD[0].where(landcover_NLCD[0].isin(NLDC_valid_covers))
    
    #average RAP
    RAP_grass = landcover_RAP[0].mean().values + landcover_RAP[3].mean().values
    RAP_shrub = landcover_RAP[4].mean().values
    RAP_tree = landcover_RAP[5].mean().values
    
    RAP_result = np.array([RAP_grass, RAP_shrub, RAP_tree])
    RAP_result = (RAP_result / np.sum(RAP_result))*100
    
    print(RAP_result)
    
    classes = ['grassland', 'grass-shrub', 'grass-tree']
    
    RAP_class=''
    for pft in RAP_thresholds.keys():
      if (RAP_result[0] >= RAP_thresholds[pft]['FG'][0] and RAP_result[0] <= RAP_thresholds[pft]['FG'][1] and
          RAP_result[1] >= RAP_thresholds[pft]['SHR'][0] and RAP_result[1] <= RAP_thresholds[pft]['SHR'][1] and
          RAP_result[2] >= RAP_thresholds[pft]['TRE'][0] and RAP_result[2] <= RAP_thresholds[pft]['TRE'][1]):
          RAP_class = pft

      print(RAP_thresholds[pft]['FG'][0])
    
    #RAP_class = classes[np.argmax(RAP_result)]
    print(f'RAP class: {RAP_class}')
    
    RAP_classes.append(RAP_class)
    
    
    #assign class based on thresholds
    
    
    #mode NLCD class
    NLCD_class = int(stats.mode(landcover_NLCD[0].values, axis=None, keepdims=False, nan_policy = 'omit')[0])
    NLCD_class_str = NLCD_LUT[str(NLCD_class)]
    print(f'NLCD class: {NLCD_class}, {NLCD_class_str}')
    NLCD_classes.append(NLCD_class_str)
    
df_predicted_class = pd.DataFrame({'site': sites, 'nlcd_predicted_pft': NLCD_classes, 'rap_predicted_pft': RAP_classes})
  
df_merged_classes = pft_key_df.merge(df_predicted_class, on='site')

cm_nlcd = confusion_matrix(df_merged_classes['true_pft'], df_merged_classes['nlcd_predicted_pft'])
cm_rap = confusion_matrix(df_merged_classes['true_pft'], df_merged_classes['rap_predicted_pft'])

fig, ax = plt.subplots()
cm_display = ConfusionMatrixDisplay(confusion_matrix = cm_nlcd, display_labels = ['grass-shrub', 'grass-tree', 'grassland', 'hay-pasture'])

cm_display.plot()
fig.tight_layout()
plt.savefig('output/landcover_classification/nlcd_cm.jpg', dpi=300)
plt.show()

fig, ax = plt.subplots()
cm_display = ConfusionMatrixDisplay(confusion_matrix = cm_rap, display_labels = ['grass-shrub', 'grass-tree', 'grassland', 'hay-pasture'])

cm_display.plot()
fig.tight_layout()
plt.savefig('output/landcover_classification/rap_cm.jpg', dpi=300)
plt.show()
    
    