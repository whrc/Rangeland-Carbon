from typing import List, Optional
from dataclasses import dataclass, field

@dataclass
class RCTMConfig:
    """
    RCTM config class.
    """

    ############### setup params #########################
    
    workflows_path: Optional[str] = None
    path_to_temp_dir: Optional[str] = None
    
    bucket_name: Optional[str] = None
    
    gcloud_project: Optional[str] = None
    service_account: Optional[str] = None
    gee_key_json: Optional[str] = None
    
    ############### GEE pipeline params #########################
    
    geometry_polygon: Optional[str] = None
    path_to_geometry_local: Optional[str] = None
    ee_geometry_save_dir: Optional[str] = None
    shape_name_col: Optional[str] = None
    path_to_existing_ee_geom: Optional[str] = None
    overwrite_geom_assets: Optional[bool] = None
    
    sfm_imagery_download_date_range: List[str] = field(default_factory=lambda: ['2002-01-01', '2003-01-01'])
    LS8_collection: Optional[str] = 'LANDSAT/LC08/C02/T1_L2'
    LS7_collection: Optional[str] = 'LANDSAT/LE07/C02/T1_L2'
    LS5_collection: Optional[str] = 'LANDSAT/LT05/C02/T1_L2'
    landsat_save_dir: Optional[str] = None
    landsat_overwrite: Optional[bool] = None
    
    MODIS_collection: Optional[str] = 'MODIS/061/MCD43A4'
    modis_date_download_interval: Optional[int] = 5
    modis_save_dir: Optional[str] = None
    modis_overwrite: Optional[bool] = None
    
    covariates_save_dir: Optional[str] = None
    covariates_overwrite: Optional[bool] = None
    
    landcover_save_dir: Optional[str] = None
    landcover_overwrite: Optional[bool] = None
    
    ################ RCTM pipeline params ########################
    
    modis_smooth_in_dir: Optional[str] = None
    modis_smooth_out_dir: Optional[str] = None
    
    starfm_in_modis_dir : Optional[str] = None
    starfm_in_landsat_dir : Optional[str] = None
    starfm_out_dir: Optional[str] = None
    starfm_source: Optional[str] = None
    starfm_config: Optional[str] = None
    
    RCTM_input_dir: Optional[str] = None
    spatial_param_outname: Optional[str] = None
    fused_landcover_outname: Optional[str] = None
    
    #model config
    point_mode: Optional[bool] = None
    spin_years: Optional[int] = None
    run_transient: Optional[bool] = None
    init_C_stocks_with_image: Optional[bool] = None
    force_pft: Optional[str] = None
    
    #params
    path_to_RCTM_params: Optional[str] = None
    path_to_RCTM_spatial_params: Optional[str] = None
    
    #C stock initialization files
    C_stock_inits_yaml: Optional[str] = None
    C_stock_init_image: Optional[str] = None
    
    #spinup covariate files
    path_to_spin_covariates_point: Optional[str] = None
    path_to_spin_covariates_spatial: Optional[str] = None
    
    #spinup C stock output
    C_stock_spin_out_path_point: Optional[str] = None
    C_stock_spin_out_path: Optional[str] = None
    
    #spinup figure output
    spatial_spin_fig_path: Optional[str] = None
    
    #transient covariate inputs
    transient_covariate_path: Optional[str] = None
    
    #transient results output
    transient_C_stock_hist: Optional[str] = None
    transient_flux_hist: Optional[str] = None

@dataclass
class RCTM_params:
    RCTM_params: Optional[dict] = None
    
@dataclass
class C_stock_inits:
    AGB: Optional[float] = None
    BGB: Optional[float] = None
    AGL: Optional[float] = None
    BGL: Optional[float] = None
    POC: Optional[float] = None
    HOC: Optional[float] = None
    

    
    
