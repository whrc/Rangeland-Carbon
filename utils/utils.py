
from google.cloud import storage
import ee
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import utils
import geemap
import geopandas as gpd
import time
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.metrics import r2_score
import xarray as xr
from matplotlib import pyplot as plt
import seaborn as sns
import os

path_to_temp = '/home/amullen/temp/'
bucket_name = 'rangelands'
storage_client = storage.Client.from_service_account_json('/home/amullen/Rangeland-Carbon/res/gee_key.json')
bucket = storage_client.get_bucket(bucket_name)

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

def get_task_status():
  tasks = ee.batch.Task.list()
  for task in tasks:
      print(task.status())
      #task_id = task.status()['id']
      #task_state = task.status()['state']
      #if task_state == 'FAILED':
      #    print(task.status())

          
def print_root_assets():
  folders = ee.data.getAssetRoots()
  for folder in folders:
    print(folder['id'] + ', type: ' + folder['type'])
            
def authorize(high_endpoint=False):
  service_account = 'rctm-07dea6fa718fdd0921f124263@rangelands-explo-1571664594580.iam.gserviceaccount.com'
  credentials = ee.ServiceAccountCredentials(service_account, '/home/amullen/Rangeland-Carbon/res/gee_key.json')
  
  # Initialize the library.
  if high_endpoint==True:
    ee.Initialize(credentials, project='rangelands-explo-1571664594580', opt_url='https://earthengine-highvolume.googleapis.com')
  else:
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


def upload_features_to_earth_engine(shapefile_path, username, asset_prefix, output_txt):
    # Initialize Earth Engine
    authorize()
    
    # Load your shapefile using GeoPandas
    gdf = gpd.read_file(shapefile_path)
    print(gdf.head())
    asset_paths = []  # To store uploaded asset paths
    crs=gdf.crs
    # Loop through each row in the GeoDataFrame
    ids = []
    for index, row in gdf.iterrows():
        name = row['Grid']
        row = gpd.GeoDataFrame(row.to_dict(), index=[0])
        row = row.drop(columns = ['Shape_Leng', 'Shape_Area'])
        row.crs = crs

        ee_geometry = geemap.geopandas_to_ee(row)
        ee_feature = ee.FeatureCollection(ee_geometry)
        
        # Generate an asset ID for the feature
        asset_id = 'projects/{}/assets/{}/{}'.format(username, asset_prefix, name)
        print(asset_id)
        # Upload the feature to Earth Engine as an asset
        try:
            task = ee.batch.Export.table.toAsset(ee_feature, name, asset_id)
            task.start()
            print(f"Uploaded {asset_id}")
            asset_paths.append(asset_id)
            ids.append(task.id)
        except Exception as e:
            print(f"Error uploading {asset_id}: {e}")
    
    # Wait for assets to upload and export paths to a text file
    with open(output_txt, 'w') as txt_file:
        for i, asset_path in enumerate(asset_paths):
          
          while True:
            status=ee.data.getTaskStatus(ids[i])[0]
            print(status)
            
            if status['state'] == 'COMPLETED':
              txt_file.write(asset_path + '\n')
              break
              
            if status['state'] == 'FAILED':
              txt_file.write(asset_path + '\n')
              break
        
            time.sleep(5)  # Wait for 5 seconds before checking again

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





