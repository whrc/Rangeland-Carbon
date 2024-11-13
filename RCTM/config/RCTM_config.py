from typing import List, Optional
from dataclasses import dataclass, field
import os

@dataclass
class RCTMConfig:
    """
    RCTM config class.
    """

    ############### setup params #########################
    
    workflows_path: str
    path_to_temp_dir: str
    
    bucket_name: str
    gcloud_workflow_base_dir: str
    
    gcloud_project: str
    service_account: str
    gee_key_json: str
    
    ############### GEE pipeline params #########################
    
    geometry_polygon: Optional[str] = 'single'
    path_to_geometry_local: Optional[str] = None
    ee_geometry_save_dir: Optional[str] = None
    shape_name_col: Optional[str] = None
    path_to_existing_ee_geom: Optional[str] = None
    overwrite_geom_assets: Optional[bool] = False
    
    sfm_imagery_download_date_range: List[str] = field(default_factory=lambda: ['2002-01-01', '2003-01-01'])
    LS8_collection: Optional[str] = 'LANDSAT/LC08/C02/T1_L2'
    LS7_collection: Optional[str] = 'LANDSAT/LE07/C02/T1_L2'
    LS5_collection: Optional[str] = 'LANDSAT/LT05/C02/T1_L2'
    landsat_overwrite: Optional[bool] = False
    
    MODIS_collection: Optional[str] = 'MODIS/061/MCD43A4'
    modis_date_download_interval: Optional[int] = 5
    
    modis_overwrite: Optional[bool] = False
    covariates_overwrite: Optional[bool] = False
    landcover_overwrite: Optional[bool] = False

    #optional
    landsat_save_dir: Optional[str] = None
    modis_save_dir: Optional[str] = None
    covariates_save_dir: Optional[str] = None
    landcover_save_dir: Optional[str] = None
    
    ################ RCTM pipeline params ########################

    RCTM_spinup_date_range: List[str] = field(default_factory=lambda: ['2002-01-01', '2005-12-31'])
    RCTM_transient_input_gen_date_range: List[str] = field(default_factory=lambda: ['2002-01-01', '2023-12-31'])
    RCTM_transient_run_date_range: List[str] = None

    starfm_source: Optional[str] = os.path.join(os.getenv('RCTMPATH'), 'remote_sensing/starfm_source/')
    starfm_config: Optional[str] = os.path.join(os.getenv('RCTMPATH'), 'config/input_ref.txt')

    #optional
    modis_smooth_in_dir: Optional[str] = None
    modis_smooth_out_dir: Optional[str] = None
    
    starfm_in_modis_dir : Optional[str] = None
    starfm_in_landsat_dir : Optional[str] = None
    starfm_out_dir: Optional[str] = None
    RCTM_input_dir: Optional[str] = None
    spatial_param_outname: Optional[str] = None
    fused_landcover_outname: Optional[str] = None
    
    #model config
    point_mode: Optional[bool] = False
    spin_years: Optional[int] = 100
    run_transient: Optional[bool] = True
    init_C_stocks_with_image: Optional[bool] = False
    force_pft: Optional[str] = None
    
    #params
    path_to_RCTM_params: Optional[str] = os.path.join(os.getenv('RCTMPATH'), 'templates/RCTM_params.yaml')
    
    #optional
    path_to_RCTM_spatial_params: Optional[str] = None

    #C stock initialization files
    C_stock_inits_yaml: Optional[str] = os.path.join(os.getenv('RCTMPATH'), 'templates/C_stock_inits.yaml')
    C_stock_init_image: Optional[str] = None

    #spinup covariate files
    path_to_spin_covariates_point: Optional[str] = None
    path_to_spin_covariates_spatial: Optional[str] = None
        
    #transient covariate inputs
    transient_covariate_path: Optional[str] = None

    #spinup C stock output
    C_stock_spin_out_path_point: Optional[str] = None
    C_stock_spin_out_path: Optional[str] = None
        
    #spinup figure output
    spatial_spin_fig_path: Optional[str] = None
        
    #transient results output
    transient_C_stock_hist: Optional[str] = None
    transient_flux_hist: Optional[str] = None

    def __post_init__(self):
        
        ############### GEE pipeline params #########################
        self.landsat_save_dir=define_if_none(self.landsat_save_dir, os.path.join(self.gcloud_workflow_base_dir, 'landsat/'))
        self.modis_save_dir=define_if_none(self.modis_save_dir, os.path.join(self.gcloud_workflow_base_dir, 'modis/'))
        self.covariates_save_dir=define_if_none(self.covariates_save_dir, os.path.join(self.gcloud_workflow_base_dir, 'covariates/'))
        self.landcover_save_dir=define_if_none(self.landcover_save_dir, os.path.join(self.gcloud_workflow_base_dir, 'landcover/'))

        #RCTM pipeline params
        self.modis_smooth_in_dir=define_if_none(self.modis_smooth_in_dir, os.path.join(self.gcloud_workflow_base_dir, 'modis/'))
        self.modis_smooth_out_dir=define_if_none(self.modis_smooth_out_dir, os.path.join(self.gcloud_workflow_base_dir, 'modis_smooth/'))

        self.starfm_in_modis_dir=define_if_none(self.starfm_in_modis_dir, os.path.join(self.gcloud_workflow_base_dir, 'modis_smooth/'))
        self.starfm_in_landsat_dir=define_if_none(self.starfm_in_landsat_dir, os.path.join(self.gcloud_workflow_base_dir, 'landsat/'))
        self.starfm_out_dir=define_if_none(self.starfm_out_dir, os.path.join(self.gcloud_workflow_base_dir, 'starfm/'))
        self.RCTM_input_dir=define_if_none(self.RCTM_input_dir, os.path.join(self.gcloud_workflow_base_dir, 'RCTM_ins/'))
        self.spatial_param_outname=define_if_none(self.spatial_param_outname, os.path.join(self.gcloud_workflow_base_dir, 'params/spatial_params.tif'))
        self.fused_landcover_outname=define_if_none(self.fused_landcover_outname, os.path.join(self.gcloud_workflow_base_dir, 'landcover/fused_landcover.tif'))

        #params
        self.path_to_RCTM_spatial_params=define_if_none(self.path_to_RCTM_spatial_params, os.path.join(self.gcloud_workflow_base_dir, 'params/spatial_params.tif'))
        
        #spinup covariate files
        self.path_to_spin_covariates_point=define_if_none(self.path_to_spin_covariates_point, os.path.join(self.gcloud_workflow_base_dir, 'RCTM_ins/RCTM_spin_inputs.csv'))
        self.path_to_spin_covariates_spatial=define_if_none(self.path_to_spin_covariates_spatial, os.path.join(self.gcloud_workflow_base_dir, 'RCTM_ins/RCTM_spin_inputs.nc'))

        #transient covariate inputs
        self.transient_covariate_path=define_if_none(self.transient_covariate_path, os.path.join(self.gcloud_workflow_base_dir, 'RCTM_ins/RCTM_inputs.nc'))

        #spinup C stock output
        self.C_stock_spin_out_path_point=define_if_none(self.C_stock_spin_out_path_point, os.path.join(self.gcloud_workflow_base_dir, 'RCTM_output/spinup/RCTM_C_stocks_spin_outputs.csv'))
        self.C_stock_spin_out_path=define_if_none(self.C_stock_spin_out_path, os.path.join(self.gcloud_workflow_base_dir, 'RCTM_output/spinup/RCTM_C_stocks_spin_output.tif'))

        #spinup figure output
        self.spatial_spin_fig_path=define_if_none(self.spatial_spin_fig_path, os.path.join(self.gcloud_workflow_base_dir, 'RCTM_output/spinup/figs/spin_fig.jpg'))

        #transient results output
        self.transient_C_stock_hist=define_if_none(self.transient_C_stock_hist, os.path.join(self.gcloud_workflow_base_dir, 'RCTM_output/transient/C_stock_hist.nc'))
        self.transient_flux_hist=define_if_none(self.transient_flux_hist, os.path.join(self.gcloud_workflow_base_dir, 'RCTM_output/transient/flux_hist.nc'))

def define_if_none(var, default_val):
    
    if var is None:
        return default_val   
    return var

@dataclass
class RCTM_params:
    RCTM_params: Optional[dict] = None
    
@dataclass
class C_stock_inits: #gC/m2
    AGB: Optional[float] = 5.91 
    BGB: Optional[float] = 29.99
    AGL: Optional[float] = 50.12
    BGL: Optional[float] = 48.63
    POC: Optional[float] = 1562.44
    HOC: Optional[float] = 2907.64

    

    
    
