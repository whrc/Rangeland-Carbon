import pandas as pd
import omegaconf
from google.cloud import storage
import os
import sys

from RCTM.remote_sensing.starfm_preprocessing import smooth_modis_col
from RCTM.remote_sensing.starfm import get_runlist, run_starfm
from RCTM.remote_sensing.starfm_postprocessing import gen_covariates, image_average_variables, aggregate_for_spinup, get_spatial_RCTM_params
from RCTM.config.RCTM_config import RCTMConfig, RCTM_params
from RCTM.utils import utils

class RCTMPrePipeline(object):
  
  PREPROCESS_STATUS_COLS = ['ee_geom_path', 'modis_smoothing', 'starfm', 'input_formatting', 'landcover_fusion']
  
  def __init__(
                self,
                config_filename: str = None,
                rctm_param_filename: str = None,
                default_config: str = ''
            ):
            
    # Configuration file intialization
    if config_filename is None:
        
      config_filename = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), default_config)
      logging.info(f'Loading default config: {config_filename}')
      
    self.conf = self._read_config(config_filename, RCTMConfig)
    self.RCTM_params = self._read_config(rctm_param_filename, config_class = RCTM_params)
      
    # Authenticate gcloud session
    self.storage_client = storage.Client.from_service_account_json(self.conf['gee_key_json'])
    self.bucket = self.storage_client.get_bucket(self.conf.bucket_name)
      
    # create new status file if doesn't already exists
    self.preprocess_status_path = os.path.join(self.conf['workflows_path'], 'preprocess_status.csv')
    
    if not os.path.isfile(self.preprocess_status_path):
      preprocess_status = pd.DataFrame(columns = self.PREPROCESS_STATUS_COLS)
      preprocess_status.to_csv(self.preprocess_status_path, index=False)
      preprocess_status = None
    
  def _read_config(self, filename: str, config_class=RCTMConfig):
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
        
  # -------------------------------------------------------------------------
  # Getters
  # -------------------------------------------------------------------------       
  def get_preprocess_status(self):
    
    return pd.read_csv(self.preprocess_status_path)
  
  # -------------------------------------------------------------------------
  # Writing
  # -------------------------------------------------------------------------    
  def write_preprocess_status(self, preprocess_status_file):
    
    preprocess_status_file.to_csv(self.preprocess_status_path, index=False)
    
    return
    
  # -------------------------------------------------------------------------
  # Processes
  # -------------------------------------------------------------------------    
  def smooth_modis(self):
      
    smooth_modis_col(self.conf.modis_smooth_in_dir, self.conf.modis_smooth_out_dir, self.bucket, self.conf.path_to_temp_dir, windowsize='20d', min_periods=1)
    
    return
    
  def starfm(self):
  
    runlist=get_runlist(self.conf.starfm_in_modis_dir, self.conf.starfm_in_landsat_dir, self.conf.starfm_out_dir, self.bucket)
    run_starfm(runlist, self.conf.path_to_temp_dir, self.bucket, self.conf.starfm_source)
    
    return
  
  def starfm_postprocessing(self):
  
    #gen covariates
    ds = gen_covariates(self.conf.starfm_out_dir, self.conf.covariates_save_dir,  self.conf.RCTM_input_dir, 'RCTM_inputs.nc', self.bucket, '2002-01-01', '2005-12-31', self.conf.path_to_temp_dir, gap_fill = True)
    #average indices
    df = image_average_variables(ds, ['ndvi','srad','vpd','tsoil','sm1','sm2','shortwave_radition','tavg','tmin','prcp','clay'], self.bucket, self.conf.path_to_temp_dir, plot_dir=os.path.join(self.conf.RCTM_input_dir, 'transient_figs/'))
    df.to_csv('gs://' + self.conf.bucket_name + '/' + self.conf.RCTM_input_dir + 'RCTM_inputs.csv')
    #gen spinup
    spin_ds = aggregate_for_spinup(ds, self.conf.RCTM_input_dir, 'RCTM_spin_inputs.nc', '2002-01-01', '2005-12-31', self.bucket,  self.conf.path_to_temp_dir, period = 5)
    df = image_average_variables(spin_ds, ['ndvi','srad','vpd','tsoil','sm1','sm2','shortwave_radition','tavg','tmin','prcp','clay'], self.bucket, self.conf.path_to_temp_dir, plot_dir=os.path.join(self.conf.RCTM_input_dir, 'spin_figs/'))
    df.to_csv('gs://' + self.conf.bucket_name + '/' + self.conf.RCTM_input_dir + 'RCTM_spin_inputs.csv')
    #get spatial params
  
  def gen_RCTM_params(self):
  
    get_spatial_RCTM_params(os.path.join(self.conf.landcover_save_dir, 'NLCD_2019.tif'), os.path.join(self.conf.landcover_save_dir, 'RAP_2019.tif'), self.RCTM_params, self.conf.fused_landcover_outname, self.conf.spatial_param_outname, self.bucket, self.conf.path_to_temp_dir, param_type='starfm')
    
pipe=RCTMPrePipeline(config_filename='/home/amullen/Rangeland-Carbon/examples/config/test_config.yaml', rctm_param_filename='/home/amullen/Rangeland-Carbon/RCTM/modeling/RCTM_params.yaml')
#pipe.smooth_modis()
#pipe.starfm()
#pipe.starfm_postprocessing()
#pipe.gen_RCTM_params()
