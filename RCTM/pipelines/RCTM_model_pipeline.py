import omegaconf
from google.cloud import storage
import os
import sys
import pandas as pd

from RCTM.modeling.RCTM import main
from RCTM.config.RCTM_config import RCTMConfig, RCTM_params, C_stock_inits


class RCTMPipeline(object):
  
  RCTM_STATUS_COLS = ['ee_geom_path', 'mode', 'force_pft', 'spinup_years', 'transient_start', 'transient_end']
  
  def __init__(
                self,
                config_filename: str = None,
                default_config: str = ''
            ):
            
    # Configuration file intialization
    if config_filename is None:
      logging.info('No config file supplied')
      return
      
    self.conf = self._read_config(config_filename, RCTMConfig)
    self.RCTM_params = self._read_config(self.conf.path_to_RCTM_params, config_class = RCTM_params).__dict__
    self.C_stock_init_yaml = self._read_config(self.conf.C_stock_inits_yaml, config_class = C_stock_inits).__dict__
      
    # Authenticate gcloud session
    self.storage_client = storage.Client.from_service_account_json(self.conf.gee_key_json)
    self.bucket = self.storage_client.get_bucket(self.conf.bucket_name)
      
    # create new status file if doesn't already exists
    self.RCTM_status_path = os.path.join(self.conf.workflows_path, 'RCTM_status.csv')
    
    if not os.path.isfile(self.RCTM_status_path):
      RCTM_status = pd.DataFrame(columns = self.RCTM_STATUS_COLS)
      RCTM_status.to_csv(self.RCTM_status_path, index=False)
      RCTM_status = None
    
  def _read_config(self, filename: str, config_class=RCTMConfig):
    """
    Read configuration filename and initiate objects
    """
    # Configuration file initialization
    conf = omegaconf.OmegaConf.load(filename)
    conf_dict = omegaconf.OmegaConf.to_container(conf, resolve=True)

    try:
      conf = config_class(**conf_dict)

    except BaseException as err:
      sys.exit(f"ERROR: {err}")

    return conf

  def run_RCTM(self):
    ## runs model ##
    
    main(force_pft=self.conf.force_pft, 
         point_mode=self.conf.point_mode, 
         spin_years=self.conf.spin_years, 
         run_transient=self.conf.run_transient,
         transient_date_range=self.conf.RCTM_transient_run_date_range, 
         init_C_stocks_with_image=self.conf.init_C_stocks_with_image,
         RCTM_params_point = self.RCTM_params, 
         path_to_RCTM_spatial_params = self.conf.path_to_RCTM_spatial_params, 
         path_to_spin_covariates_point = self.conf.path_to_spin_covariates_point, 
         path_to_spin_covariates_spatial = self.conf.path_to_spin_covariates_spatial, 
         C_stock_spin_out_path_point = self.conf.C_stock_spin_out_path_point,
         C_stock_spin_out_path = self.conf.C_stock_spin_out_path,
         C_stock_inits_yaml = self.C_stock_init_yaml, 
         C_stock_init_image = self.conf.C_stock_init_image, 
         spatial_spin_fig_path = self.conf.spatial_spin_fig_path, 
         transient_covariate_path = self.conf.transient_covariate_path,
         transient_C_stock_hist = self.conf.transient_C_stock_hist, 
         transient_flux_hist = self.conf.transient_flux_hist,
         bucket = self.bucket,
         path_to_temp = self.conf.path_to_temp_dir)

