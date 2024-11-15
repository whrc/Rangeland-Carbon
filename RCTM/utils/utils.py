
from google.cloud import storage
import ee
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import geemap
import geopandas as gpd
import time
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
#import tensorflow as tf
#from tensorflow import keras
#from tensorflow.keras import layers
from sklearn.metrics import r2_score
import xarray as xr
from matplotlib import pyplot as plt
import seaborn as sns
import os
import itertools
import logging
  
def xr_dataset_to_data_array(dataset):

  band_names = list(dataset.attrs['long_name'])
  band_name_dict = {}
  for i, band_name in enumerate(band_names):
    band_name_dict[i+1] = band_name 
  dataset = dataset.to_dataset(dim='band')
  dataset = dataset.rename_vars(name_dict=band_name_dict)
  return dataset
  
def image_average_variables(ds, variable_list, time_index='time', plot_dir=None, gs=False):
  """Take spatial average of xarray Dataset for variables in variable_list. Returns pandas dataframe

  Args:
    ds (xarray.Dataset): dataset with variables to average and dimensions x, y, and time
    variable_list (list of strings): variable names to spatially average
    plot (string): full path to save timeseries plot to, will only plot if this argument is set

  Returns:
     pandas.DataFrame with a column for time, and spatial average of each variable in variable_list
  """
  if time_index=='time':
    site_dates = pd.to_datetime(pd.DatetimeIndex(ds.indexes[time_index].to_datetimeindex()))
  else:
    site_dates = ds.indexes[time_index]
    
  mean_indices=ds[variable_list].mean(dim=['y', 'x']).to_pandas()
  mean_indices = mean_indices.sort_values(by=time_index, ascending=True)

  for variable in variable_list:
    if plot_dir!=None:
      fig, ax = plt.subplots()
      sns.lineplot(x=site_dates, y=mean_indices[variable].values, label=variable)
      plt.ylabel(variable)
      print('saving plot')
      
      if gs==True:
        plt.savefig(os.path.join(path_to_temp, 'temp_plot.jpg'), dpi=300)
        utils.gs_write_blob(os.path.join(path_to_temp, 'temp_plot.jpg'), plot_dir+variable+'.jpg', bucket)
      else:
        plt.savefig(os.path.join(plot_dir, variable+'.jpg'), dpi=300)
        
      plt.show()
      plt.clf()

  return mean_indices

def cancel_running_tasks():
  tasks = ee.batch.Task.list()
  for task in tasks:
      task_id = task.status()['id']
      task_state = task.status()['state']
      if task_state == 'RUNNING' or task_state == 'READY':
          task.cancel()
          print('Task {} canceled'.format(task_id))
      else:
          print('Task {} state is {}'.format(task_id, task_state))

def get_task_status(task_id):
  return ee.data.getTaskStatus(task_id)[0]
  #tasks = ee.batch.Task.list()
  #for task in tasks:
      #print(task.status())
      #task_id = task.status()['id']
      #task_state = task.status()['state']
      #if task_state == 'FAILED':
      #    print(task.status())

          
def print_root_assets():
  folders = ee.data.getAssetRoots()
  for folder in folders:
    print(folder['id'] + ', type: ' + folder['type'])
            
def authorize(gcloud_project = None, gee_key_json = None, service_account = None, high_endpoint=False):
  credentials = ee.ServiceAccountCredentials(service_account, gee_key_json)
  
  # Initialize the library.
  if high_endpoint==True:
    ee.Initialize(credentials, project=gcloud_project, opt_url='https://earthengine-highvolume.googleapis.com')
  else:
    ee.Initialize(credentials, project=gcloud_project)
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

def write_read(bucket_name, blob_name):
    """Write and read a blob from GCS using file-like IO"""
    # The ID of your GCS bucket
    # bucket_name = "your-bucket-name"

    # The ID of your new GCS object
    # blob_name = "storage-object-name"

    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(blob_name)

    # Mode can be specified as wb/rb for bytes mode.
    # See: https://docs.python.org/3/library/io.html

    with blob.open("r") as f:
        print(f.read())
        
def asset_exists(asset_id: str) -> bool:
    """Check the existence of the asset"""
    try:
        ee.data.getAsset(asset_id)
    except ee.EEException:
        exists = False
    else:
        exists = True
    return exists
    
def append_col_overwrite(df, value_arrays, cols, id_col):
    value_arrays = list(map(list, itertools.zip_longest(*value_arrays, fillvalue=None)))
    #create new dataframe from arrays of values and list of columns
    new_data = pd.DataFrame(data = value_arrays, columns = cols)
    
    #filter items from old dataframe containing id from new dataframe
    df = df.loc[~df[id_col].isin(new_data[id_col])]
    
    #concatenate two dataframes   
    df = pd.concat([df, new_data], ignore_index=True)
    
    return df
    
def principal_component_regression(dataframe, column_names, target_col, n_components=3):

    dataframe = dataframe.dropna(subset=column_names + [target_col])
    # Extract the predictor variables (X) and the target variable (y)
    X = dataframe[column_names].values
    y = dataframe[target_col].values

    # Standardize the predictor variables (optional, but recommended for PCA)
    X_standardized = (X - X.mean(axis=0)) / X.std(axis=0)

    # Perform PCA to reduce the dimensionality of the predictor variables
    pca = PCA()
    X_pca = pca.fit_transform(X_standardized)

    # Use the first n principal components to create the transformed dataset
    n_components = X.shape[1]  # Use all components initially
    X_pca_final = X_pca[:, :n_components]

    # Fit a linear regression model on the transformed dataset
    model = LinearRegression()
    model.fit(X_pca_final, y)

    # Predict the target variable using the fitted model
    y_pred = model.predict(X_pca_final)

    # Calculate the R-squared score
    r2 = r2_score(y, y_pred)

    return r2, model, pca

def ee_export_table_to_asset(ee_feature, name, asset_id):
  
  try:
    task = ee.batch.Export.table.toAsset(ee_feature, name, asset_id)
    task.start()
    print(f"Uploaded {asset_id}")
            
  except Exception as e:
    print(f"Error uploading {asset_id}: {e}")
    
  return task.id

def upload_features_to_ee(geojson_path, asset_dir, name_col = 'Name', drop_cols = ['Shape_Leng', 'Shape_Area'], overwrite=False):
    # Load your shapefile using GeoPandas
    gdf = gpd.read_file(geojson_path)

    #print(gdf.head())
    asset_paths = []  # To store uploaded asset paths
    crs=gdf.crs
    # Loop through each row in the GeoDataFrame
    task_ids = []
    
    for index, row in gdf.iterrows():
        if name_col not in gdf.columns:
           name = f'geometry{index}'
        else:
           name = row[name_col]
        row = gpd.GeoDataFrame(row.to_dict(), index=[0])
        #row = row.drop(columns = drop_cols)
        row.crs = crs

        ee_geometry = geemap.geopandas_to_ee(row)
        ee_feature = ee.FeatureCollection(ee_geometry)
        
        # Generate an asset ID for the feature
        asset_id = os.path.join(asset_dir, name)
        task_id = ''
        #print(asset_id)
        if asset_exists(asset_id):
          
          if overwrite:
            print(f"{asset_id} exists, overwriting")
            ee.data.deleteAsset(asset_id)
            task_id = ee_export_table_to_asset(ee_feature, name, asset_id)
          
          else:
            print(f"{asset_id} exists, skipping")
            task_id = None
            
        # Upload the feature to Earth Engine as an asset
        else:
          task_id = ee_export_table_to_asset(ee_feature, name, asset_id)
            
        asset_paths.append(asset_id)
        task_ids.append(task_id)
            
        return asset_paths, task_ids

def predict_column_with_neural_network(dataframe, features_columns, target_column, epochs=50, batch_size=32):
    dataframe = dataframe.dropna(subset=features_columns + [target_column])
    # Split the data into features and target
    features = dataframe[features_columns]
    target = dataframe[target_column]

    # Split the data into training and testing sets
    features_train, features_test, target_train, target_test = train_test_split(features, target, test_size=0.4, random_state=42)

    # Standardize features
    scaler = StandardScaler()
    features_train_scaled = scaler.fit_transform(features_train)
    features_test_scaled = scaler.transform(features_test)

    # Build the neural network model
    model = keras.Sequential([
        layers.Input(shape=(len(features_columns),)),
        layers.Dense(64, activation='relu'),
        layers.Dense(32, activation='relu'),
        layers.Dense(1)  # Output layer
    ])

    model.compile(optimizer='adam', loss='mean_squared_error')

    # Train the model
    model.fit(features_train_scaled, target_train, epochs=epochs, batch_size=batch_size, verbose=1)

    # Predict the target column on the test data
    predictions = model.predict(features_test_scaled)

    # Calculate Mean Squared Error
    mse = mean_squared_error(target_test, predictions)
    r2 = r2_score(target_test, predictions)
    print(f"Mean Squared Error: {mse}")
    print(f"r2: {r2}")

    return model, scaler

def get_sfm_date_df(in_sfm_dir, bucket):
  """gets dataframe of MODIS images with date in Y_m_d format and milliseconds, image name, and full path

  Args:
    in_modis_dir (string): path to google storage directory to list (should not include bucket name, should not start with forward slash)
    bucket (storage.bucket): google cloud bucket

  Returns:
     pandas dataframe
  """ 
  im_paths = gs_listdir(in_sfm_dir, bucket)
  im_paths = [k for k in im_paths if '.tif' in k]
  
  df_im_paths = pd.DataFrame({'im_path':im_paths})
  
  df_im_paths['im_name'] = [i[-1] for i in df_im_paths['im_path'].str.split('/').to_list()]

  df_im_paths['im_date'] = pd.to_datetime(df_im_paths['im_path'].str[-14:-4], format='%Y-%m-%d')
  df_im_paths['millis'] = df_im_paths['im_date'].astype(np.int64) / int(1e6)
  df_im_paths = df_im_paths.set_index(df_im_paths['im_date'])
  df_im_paths = df_im_paths.sort_index()
  
  return df_im_paths

def get_modis_date_df(in_modis_dir, bucket):
  """gets dataframe of MODIS images with date in Y_m_d format and milliseconds, image name, and full path

  Args:
    in_modis_dir (string): path to google storage directory to list (should not include bucket name, should not start with forward slash)
    bucket (storage.bucket): google cloud bucket

  Returns:
     pandas dataframe
  """ 
  im_paths = gs_listdir(in_modis_dir, bucket)
  im_paths = [k for k in im_paths if '.tif' in k]
  
  df_im_paths = pd.DataFrame({'im_path':im_paths})
  df_im_paths['im_name'] = [i[-1] for i in df_im_paths['im_path'].str.split('/').to_list()]
  df_im_paths['im_date'] = pd.to_datetime(df_im_paths['im_path'].str[-14:-4], format='%Y_%m_%d')
  df_im_paths['millis'] = df_im_paths['im_date'].astype(np.int64) / int(1e6)
  df_im_paths = df_im_paths.set_index(df_im_paths['im_date'])
  df_im_paths = df_im_paths.sort_index()
  
  return df_im_paths
  
def get_landsat_date_df(in_landsat_dir, bucket):
  """gets dataframe of MODIS images with date in Y_m_d format and milliseconds, image name, and full path

  Args:
    in_modis_dir (string): path to google storage directory to list (should not include bucket name, should not start with forward slash)
    bucket (storage.bucket): google cloud bucket

  Returns:
     pandas dataframe
  """ 
  im_paths = gs_listdir(in_landsat_dir, bucket)
  im_paths = [k for k in im_paths if '.tif' in k]
  
  df_im_paths = pd.DataFrame({'im_path':im_paths})
  df_im_paths['im_name'] = [i[-1] for i in df_im_paths['im_path'].str.split('/').to_list()]
  df_im_paths['im_date'] = pd.to_datetime(df_im_paths['im_path'].str[-12:-4], format='%Y%m%d')
  df_im_paths['millis'] = df_im_paths['im_date'].astype(np.int64) / int(1e6)
  df_im_paths = df_im_paths.set_index(df_im_paths['im_date'])
  df_im_paths = df_im_paths.sort_index()
  
  return df_im_paths

def addLoggingLevel(levelName, levelNum, methodName=None):
    """
    Comprehensively adds a new logging level to the `logging` module and the
    currently configured logging class.

    `levelName` becomes an attribute of the `logging` module with the value
    `levelNum`. `methodName` becomes a convenience method for both `logging`
    itself and the class returned by `logging.getLoggerClass()` (usually just
    `logging.Logger`). If `methodName` is not specified, `levelName.lower()` is
    used.

    To avoid accidental clobberings of existing attributes, this method will
    raise an `AttributeError` if the level name is already an attribute of the
    `logging` module or if the method name is already present 

    Example
    -------
    >>> addLoggingLevel('TRACE', logging.DEBUG - 5)
    >>> logging.getLogger(__name__).setLevel("TRACE")
    >>> logging.getLogger(__name__).trace('that worked')
    >>> logging.trace('so did this')
    >>> logging.TRACE
    5

    """
    if not methodName:
        methodName = levelName.lower()

    if hasattr(logging, levelName):
       raise AttributeError('{} already defined in logging module'.format(levelName))
    if hasattr(logging, methodName):
       raise AttributeError('{} already defined in logging module'.format(methodName))
    if hasattr(logging.getLoggerClass(), methodName):
       raise AttributeError('{} already defined in logger class'.format(methodName))

    # This method was inspired by the answers to Stack Overflow post
    # http://stackoverflow.com/q/2183233/2988730, especially
    # http://stackoverflow.com/a/13638084/2988730
    def logForLevel(self, message, *args, **kwargs):
        if self.isEnabledFor(levelNum):
            self._log(levelNum, message, args, **kwargs)
    def logToRoot(message, *args, **kwargs):
        logging.log(levelNum, message, *args, **kwargs)

    logging.addLevelName(levelNum, levelName)
    setattr(logging, levelName, levelNum)
    setattr(logging.getLoggerClass(), methodName, logForLevel)
    setattr(logging, methodName, logToRoot)


