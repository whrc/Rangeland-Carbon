import logging
import sys
import os
import pandas as pd
import omegaconf
from google.cloud import storage

from RCTM.remote_sensing.GEE_covariates import get_covariates, export_covariates
from RCTM.remote_sensing.GEE_landcover import export_landcover
from RCTM.remote_sensing.GEE_starfm import get_landsat, export_landsat_collection, get_MODIS, export_modis_collection
from RCTM.config.RCTM_config import RCTMConfig as Config
from RCTM.utils.utils import authorize, asset_exists, upload_features_to_ee, append_col_overwrite, get_task_status, cancel_running_tasks



class GEEPipeline(object):

  LANDCOVER_STATUS_COLS = ['ee_geom_path', 'source', 'task_id', 'task_status', 'image_full_path', 'task_submitted', 'task_completed', 'image_name']
  IMAGE_STATUS_COLS = ['ee_geom_path', 'date', 'task_id', 'task_status', 'image_full_path', 'task_submitted', 'task_completed', 'image_name']
  WORKFLOW_STATUS_COLS = ['ee_geom_path', 'ee_geom_task_id', 'ee_geom_upload_status', 'landsat_imagery', 'modis_imagery', 'covariates', 'landcover']
  #GEOM_UPLOAD_STATUS_COLS = ['ee_geom_path', 'task_id', 'status','time_completed', 'err']
  
  def __init__(
                self,
                config_filename: str = None,
                default_config: str = 'templates/watermask_default.yaml'
            ):
           
    # Configuration file intialization
    if config_filename is None:
      
      config_filename = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), default_config)
      logging.info(f'Loading default config: {config_filename}')
    
    self.conf = self._read_config(config_filename, Config)
    
    # Authorize EE
    authorize(gcloud_project = self.conf['gcloud_project'], gee_key_json = self.conf['gee_key_json'], service_account = self.conf['service_account'])
    
    # Authenticate gcloud session
    self.storage_client = storage.Client.from_service_account_json(self.conf['gee_key_json'])
    
    # create new status file if doesn't already exists
    self.workflow_status_path = os.path.join(self.conf['workflows_path'], 'workflow_process_status.csv')
    self.landsat_im_status_path = os.path.join(self.conf['workflows_path'], 'landsat_downloads.csv')
    self.modis_im_status_path = os.path.join(self.conf['workflows_path'], 'modis_downloads.csv')
    self.covariate_status_path = os.path.join(self.conf['workflows_path'], 'covariate_downloads.csv')
    self.landcover_status_path = os.path.join(self.conf['workflows_path'], 'landcover_downloads.csv')
    
    if not os.path.isfile(self.workflow_status_path):
      workflow_status = pd.DataFrame(columns = self.WORKFLOW_STATUS_COLS)
      workflow_status.to_csv(self.workflow_status_path, index=False)
      workflow_status = None
      
    if not os.path.isfile(self.landsat_im_status_path):
      landsat_im_status = pd.DataFrame(columns = self.IMAGE_STATUS_COLS)
      landsat_im_status.to_csv(self.landsat_im_status_path, index=False)
      landsat_im_status = None
      
    if not os.path.isfile(self.modis_im_status_path):
      modis_im_status = pd.DataFrame(columns = self.IMAGE_STATUS_COLS)
      modis_im_status.to_csv(self.modis_im_status_path, index=False)
      modis_im_status = None
      
    if not os.path.isfile(self.covariate_status_path):
      covariate_status = pd.DataFrame(columns = self.IMAGE_STATUS_COLS)
      covariate_status.to_csv(self.covariate_status_path, index=False)
      covariate_status = None
      
    if not os.path.isfile(self.landcover_status_path):
      landcover_status = pd.DataFrame(columns = self.LANDCOVER_STATUS_COLS)
      landcover_status.to_csv(self.landcover_status_path, index=False)
      landcover_status = None

  
  def _read_config(self, filename: str, config_class=Config):
        """
        Read configuration filename and initiate objects
        """
        # Configuration file initialization
        schema = omegaconf.OmegaConf.structured(config_class)
        conf = omegaconf.OmegaConf.load(filename)
        try:
            conf = omegaconf.OmegaConf.merge(schema, conf)
        except BaseException as err:
            sys.exit(f"ERROR: {err}")
        return conf
        
  def clean_status_files(self):
    ##
    # deletes all status files to start from clean slate
    ##
    
    os.remove(self.workflow_status_path)
    os.remove(self.landsat_im_status_path)
    os.remove(self.modis_im_status_path)
    os.remove(self.covariate_status_path)
    os.remove(self.landcover_status_path)
  
  # -------------------------------------------------------------------------
  # Getters
  # ------------------------------------------------------------------------- 
  def get_config(self):
    
    return self.conf
    
  def get_workflow_status(self):
    
    return pd.read_csv(self.workflow_status_path)
    
  def get_landsat_status(self):
  
    return pd.read_csv(self.landsat_im_status_path)
    
  def get_modis_status(self):
  
    return pd.read_csv(self.modis_im_status_path)
    
  def get_covariate_status(self):
  
    return pd.read_csv(self.covariate_status_path)
    
  def get_landcover_status(self):
  
    return pd.read_csv(self.landcover_status_path)
    
  def get_geom_upload_status(self):
    
    return pd.read_csv(self.geometry_uploads_path)
    
    
  # -------------------------------------------------------------------------
  # Writing
  # -------------------------------------------------------------------------  
  
  def write_workflow_status(self, workflow_status_file):
    
    workflow_status_file.to_csv(self.workflow_status_path, index=False)
    
    return
  
  def write_geom_upload_status(self, geom_status_file):
    
    geom_status_file.to_csv(self.geometry_uploads_path, index=False)
    
    return
    
  def write_landsat_status(self, landsat_status_file):
    
    landsat_status_file.to_csv(self.landsat_im_status_path, index=False)
    
    return
    
  def write_modis_status(self, modis_status_file):
    
    modis_status_file.to_csv(self.modis_im_status_path, index=False)
    
    return
    
  def write_covariate_status(self, covariate_status_file):
    
    covariate_status_file.to_csv(self.covariate_status_path, index=False)
    
    return
    
  def write_landcover_status(self, landcover_status_file):
    
    landcover_status_file.to_csv(self.landcover_status_path, index=False)
    
    return
    
  # -------------------------------------------------------------------------
  # Update
  # ------------------------------------------------------------------------- 
    
  def update_status(self):
    
    #update geometry upload statuses
    print('updating geometry status')
    workflow_status = self.get_workflow_status()
    workflow_status.loc[workflow_status[self.WORKFLOW_STATUS_COLS[2]]!='COMPLETED', self.WORKFLOW_STATUS_COLS[2]] = [get_task_status(task_id)['state'] for task_id in workflow_status.loc[workflow_status[self.WORKFLOW_STATUS_COLS[2]]!='COMPLETED', self.WORKFLOW_STATUS_COLS[1]]]
    self.write_workflow_status(workflow_status)
    workflow_status = None
    
    #update landsat download statuses
    print('updating landsat status')
    landsat_status = self.get_landsat_status()
    landsat_status.loc[landsat_status[self.IMAGE_STATUS_COLS[3]]!='COMPLETED', self.IMAGE_STATUS_COLS[3]] = [get_task_status(task_id)['state'] for task_id in landsat_status.loc[landsat_status[self.IMAGE_STATUS_COLS[3]]!='COMPLETED', self.IMAGE_STATUS_COLS[2]]]
    self.write_landsat_status(landsat_status)
    landsat_status = None
    
    #update modis download status
    print('updating modis status')
    modis_status = self.get_modis_status()
    modis_status.loc[modis_status[self.IMAGE_STATUS_COLS[3]]!='COMPLETED', self.IMAGE_STATUS_COLS[3]] = [get_task_status(task_id)['state'] for task_id in modis_status.loc[modis_status[self.IMAGE_STATUS_COLS[3]]!='COMPLETED', self.IMAGE_STATUS_COLS[2]]]
    self.write_modis_status(modis_status)
    modis_status = None
    
    #update covariate download statuses
    print('updating covariate status')
    covariate_status = self.get_covariate_status()
    covariate_status.loc[covariate_status[self.IMAGE_STATUS_COLS[3]]!='COMPLETED', self.IMAGE_STATUS_COLS[3]] = [get_task_status(task_id)['state'] for task_id in covariate_status.loc[covariate_status[self.IMAGE_STATUS_COLS[3]]!='COMPLETED', self.IMAGE_STATUS_COLS[2]]]
    self.write_covariate_status(covariate_status)
    covariate_status = None
    
    #update landcover download statuses
    print('updating landcover status')
    landcover_status = self.get_landcover_status()
    landcover_status.loc[landcover_status[self.LANDCOVER_STATUS_COLS[3]]!='COMPLETED', self.LANDCOVER_STATUS_COLS[3]] = [get_task_status(task_id)['state'] for task_id in landcover_status.loc[landcover_status[self.LANDCOVER_STATUS_COLS[3]]!='COMPLETED', self.LANDCOVER_STATUS_COLS[2]]]
    self.write_landcover_status(landcover_status)
    landcover_status = None
  
  def init_pipeline_from_geometry(self):
    ## Initialize pipeline from local .geojson file or existing Earth Engine geometry asset as specified in config file ##
    # 
    # Function populates the status_file.csv file located in the workflows_path dir specified in the config file.
    # If local .geojson file is specified, this function will upload it as an EE asset for future use.
    # This function does not wait for completion of the ROI upload.
    #
    ##
    
    if self.conf.path_to_existing_ee_geom is None:
      
      # upload geojson as ee asset
      print(f'no existing ee asset supplied, uploading asset from provided geojson: {self.conf.path_to_geometry_local}')
      asset_paths, task_ids = upload_features_to_ee(self.conf.path_to_geometry_local, self.conf.ee_geometry_save_dir, self.conf.shape_name_col, drop_cols = ['Shape_Leng', 'Shape_Area'], overwrite=self.conf.overwrite_geom_assets)
      
      # filter asset paths, task ids for new uploads (i.e. task ids that are not None)
      asset_paths = [asset_path for asset_path, task_id in zip(asset_paths, task_ids) if task_id is not None]
      task_ids = [task_id for task_id in task_ids if task_id is not None]
      
      # query initial statuses of submitted tasks
      upload_statuses = [get_task_status(task_id)['state'] for task_id in task_ids]
      
      # add column(s) to status file with geometry asset path, if it does not exist
      workflow_status = self.get_workflow_status()
      workflow_status = append_col_overwrite(df = workflow_status, value_arrays = [asset_paths, task_ids, upload_statuses], cols = self.WORKFLOW_STATUS_COLS[:3], id_col = self.WORKFLOW_STATUS_COLS[0])    
      self.write_workflow_status(workflow_status)
      workflow_status = None
      
    else:
      # asset already exists in earth engine, as specified in conf file
      print('Initializing using existing asset')
      
      if not asset_exists(self.conf.path_to_existing_ee_geom):
        print('config provided asset name: {self.conf.path_to_existing_ee_geom} does not exist')
      
      else:
        # write to status file with path of asset, and other required fields for asset management
        workflow_status = self.get_workflow_status()
        workflow_status = append_col_overwrite(df = workflow_status, value_arrays = [[self.conf.path_to_existing_ee_geom], ['asset already uploaded'], ['COMPLETED']], cols = self.WORKFLOW_STATUS_COLS[:3], id_col = self.WORKFLOW_STATUS_COLS[0])    
        self.write_workflow_status(workflow_status)
        workflow_status = None
    
    return
    
  def download_landsat_imagery(self):
    
    ## download landsat imagery ##
    
    #check for geometries present in workflow status file
    workflow_status = self.get_workflow_status()
    
    if len(workflow_status.loc[workflow_status[self.WORKFLOW_STATUS_COLS[2]]=='COMPLETED']) == 0:
      
      print('no uploaded geometries present in status file, have you run init_pipeline_from_geometry()?')
      
      return
    
    # iterate geometries  
    for path_to_geometry in workflow_status.loc[workflow_status[self.WORKFLOW_STATUS_COLS[2]]=='COMPLETED', self.WORKFLOW_STATUS_COLS[0]]:
      
      print(f'downloading landsat imagery for: {path_to_geometry}')
      print(f'date range: {self.conf.sfm_imagery_download_date_range}')
      
      # retrieve landsat imagery from ee and export
      landsat = get_landsat(self.conf.LS8_collection, self.conf.LS7_collection, self.conf.LS5_collection, path_to_geometry, self.conf.sfm_imagery_download_date_range[0], self.conf.sfm_imagery_download_date_range[1])
      dates, task_ids, filepaths = export_landsat_collection(landsat, path_to_geometry, self.conf.bucket_name, self.conf.landsat_save_dir, storage_client = self.storage_client, overwrite = self.conf.landsat_overwrite)
      landsat = None
      
      # get initial download statuses
      download_statuses = [get_task_status(task_id)['state'] for task_id in task_ids]
      geom_paths = [path_to_geometry]*len(download_statuses)
      
      # update landsat status file
      landsat_status = self.get_landsat_status()
      landsat_status = append_col_overwrite(df = landsat_status, value_arrays = [geom_paths, dates, task_ids, download_statuses, filepaths], cols = self.IMAGE_STATUS_COLS[:5], id_col = self.IMAGE_STATUS_COLS[1])
      self.write_landsat_status(landsat_status)
      landsat_status = None

    return
    
  def download_modis_imagery(self):
    
    ## download MODIS imagery ##
    
    # get landsat image dates
    landsat_status = self.get_landsat_status()
    landsat_dates = landsat_status[self.IMAGE_STATUS_COLS[1]]
    
    # check for geometries present in workflow status file
    workflow_status = self.get_workflow_status()
    
    if len(workflow_status.loc[workflow_status[self.WORKFLOW_STATUS_COLS[2]]=='COMPLETED']) == 0:
      
      print('no uploaded geometries present in status file, have you run init_pipeline_from_geometry()?')
      
      return
    
    # iterate geometries  
    for path_to_geometry in workflow_status.loc[workflow_status[self.WORKFLOW_STATUS_COLS[2]]=='COMPLETED', self.WORKFLOW_STATUS_COLS[0]]:
      
      print(f'downloading modis imagery for: {path_to_geometry}')
      print(f'date range: {self.conf.sfm_imagery_download_date_range}')
      
      # retrieve modis imagery from ee and export
      modis = get_MODIS(self.conf.MODIS_collection, path_to_geometry, self.conf.sfm_imagery_download_date_range[0], self.conf.sfm_imagery_download_date_range[1])
      dates, task_ids, filepaths = export_modis_collection(modis, path_to_geometry, self.conf.bucket_name, self.conf.modis_save_dir, self.storage_client, landsat_dates=landsat_dates, date_download_interval = self.conf.modis_date_download_interval, overwrite=self.conf.modis_overwrite)
      modis = None
      
      # get initial download statuses
      download_statuses = [get_task_status(task_id)['state'] for task_id in task_ids]
      geom_paths = [path_to_geometry]*len(download_statuses)
      
      # update modis status file
      modis_status = self.get_modis_status()
      modis_status = append_col_overwrite(df = modis_status, value_arrays = [geom_paths, dates, task_ids, download_statuses, filepaths], cols = self.IMAGE_STATUS_COLS[:5], id_col = self.IMAGE_STATUS_COLS[1])
      self.write_modis_status(modis_status)
      modis_status = None

    return
    
  def download_covariates(self):
    
    ## download covariates imagery ##
    # currently downloads imagery for modis dates
    # TODO: change to collect covariates regardless of modis acquisition
    
    # get modis dates
    modis_status = self.get_modis_status()
    modis_dates = modis_status[self.IMAGE_STATUS_COLS[1]]
    
    # check for geometries present in workflow status file
    workflow_status = self.get_workflow_status()
    
    if len(workflow_status.loc[workflow_status[self.WORKFLOW_STATUS_COLS[2]]=='COMPLETED']) == 0:
      
      print('no uploaded geometries present in status file, have you run init_pipeline_from_geometry()?')
      
      return
      
    # iterate geometries  
    for path_to_geometry in workflow_status.loc[workflow_status[self.WORKFLOW_STATUS_COLS[2]]=='COMPLETED', self.WORKFLOW_STATUS_COLS[0]]:
      
      print(f'downloading covariates for: {path_to_geometry}')
      print(f'date range: {self.conf.sfm_imagery_download_date_range}')
      
      # retrieve covariates from ee and export
      covariate_tuple = get_covariates(path_to_geometry, self.conf.sfm_imagery_download_date_range[0], self.conf.sfm_imagery_download_date_range[1])
      dates, task_ids, filepaths = export_covariates(covariate_tuple, self.conf.bucket_name, path_to_geometry, modis_dates, self.conf.covariates_save_dir, self.storage_client, overwrite=self.conf.covariates_overwrite)
      covariate_tuple = None
      
      # get initial statuses
      download_statuses = [get_task_status(task_id)['state'] for task_id in task_ids]
      geom_paths = [path_to_geometry]*len(download_statuses)
      
      # update covariate status file
      covariate_status = self.get_covariate_status()
      covariate_status = append_col_overwrite(df = covariate_status, value_arrays = [geom_paths, dates, task_ids, download_statuses, filepaths], cols = self.IMAGE_STATUS_COLS[:5], id_col = self.IMAGE_STATUS_COLS[1])
      self.write_covariate_status(covariate_status)
      covariate_status = None

    return
    
  def download_landcover(self):
    
    ## download lancover maps ##
    
    # check for geometries present in workflow status file
    workflow_status = self.get_workflow_status()
    
    if len(workflow_status.loc[workflow_status[self.WORKFLOW_STATUS_COLS[2]]=='COMPLETED']) == 0:
      
      print('no uploaded geometries present in status file, have you run init_pipeline_from_geometry()?')
      
      return
      
    # iterate geometries  
    for path_to_geometry in workflow_status.loc[workflow_status[self.WORKFLOW_STATUS_COLS[2]]=='COMPLETED', self.WORKFLOW_STATUS_COLS[0]]:
      
      print(f'downloading landcover imagery for: {path_to_geometry}')
      
      # retrieve landcover from ee and export
      im_names, task_ids, filepaths = export_landcover(self.conf.bucket_name, path_to_geometry, self.conf.landcover_save_dir, self.storage_client, overwrite=self.conf.landcover_overwrite)
      
      # get initial statuses
      download_statuses = [get_task_status(task_id)['state'] for task_id in task_ids]
      geom_paths = [path_to_geometry]*len(download_statuses)
      
      # update landcover status file
      landcover_status = self.get_landcover_status()
      landcover_status = append_col_overwrite(df = landcover_status, value_arrays = [geom_paths, im_names, task_ids, download_statuses, filepaths], cols = self.LANDCOVER_STATUS_COLS[:5], id_col = self.LANDCOVER_STATUS_COLS[1])
      self.write_landcover_status(landcover_status)
      landcover_status = None
          

pipe = GEEPipeline(config_filename='/home/amullen/Rangeland-Carbon/examples/config/test_config.yaml')
#pipe.clean_status_files()
#pipe.init_pipeline_from_geometry() #populates rows for each geometry in workflow status file to dictate pocessing
pipe.download_landsat_imagery() #downloads landsat imagery for regions in workflow status file
pipe.download_modis_imagery() #uses landsat image dates from landsat status file to download every five days of modis + landsat dates
pipe.download_covariates() #uses modis status file to download covariates for each modis date
pipe.download_landcover()#downloads landcover for regions in workflow status file
#cancel_running_tasks()
pipe.update_status()