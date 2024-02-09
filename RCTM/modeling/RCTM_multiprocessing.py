import RCTM
import os
from multiprocessing import Pool

roi_file='/home/amullen/Rangeland-Carbon/res/site_footprints/HLD_tiles.txt'
args=[]

with open(roi_file) as f:
  sites = [line.rstrip('\n') for line in f]
  print(sites)
  print(sites.reverse())
  
  
  for site in sites:
    ranch='/'.join(site.split('/')[-2:])

    for pft in ['grassland', 'grass-shrub', 'grass-tree']:
      if ranch=='HLD/B5' and pft=='grassland':
      
        force_pft=pft
        point_mode=False
        spin_years=0
        run_transient=True
        init_C_stocks_with_image=True
        ranch=ranch
        
        C_stock_inits_yaml = f'Ranch_Runs/{ranch}/covariates_nc/C_stocks_{force_pft}.yaml'
        path_to_RCTM_params = '/home/amullen/Rangeland-Carbon/modeling/RCTM_params.yaml'
        path_to_spin_covariates_point = f'gs://rangelands/Ranch_Runs/{ranch}/covariates_nc/covariates_spin.csv'
        path_to_RCTM_spatial_params = f'gs://rangelands/Ranch_Runs/{ranch}/covariates_nc/params.tif'
        path_to_spin_covariates_spatial = f'gs://rangelands/Ranch_Runs/{ranch}/covariates_nc/covariates_spin.nc'
        C_stock_spin_out_path_point = f'Ranch_Runs/{ranch}/covariates_nc/C_stocks_{force_pft}.yaml'
        C_stock_spin_out_path = f'Ranch_Runs/{ranch}/covariates_nc/C_stocks_{force_pft}.tif'
        C_stock_init_image = 'gs://rangelands/' + C_stock_spin_out_path
        spatial_spin_fig_path = f'output/RCTM/{ranch}_spinup_stocks.jpg'
        transient_covariate_path = f'gs://rangelands/Ranch_Runs/{ranch}/covariates_nc/covariates.nc'
        transient_C_stock_hist = f'Ranch_Runs/{ranch}/results/C_stocks_hist_{force_pft}.nc'
        transient_flux_hist = f'Ranch_Runs/{ranch}/results/flux_hist_{force_pft}.nc'
        
        args.append((force_pft, point_mode, spin_years, run_transient, init_C_stocks_with_image, ranch,
         path_to_RCTM_params, path_to_RCTM_spatial_params, path_to_spin_covariates_point,
         path_to_spin_covariates_spatial, C_stock_spin_out_path_point, C_stock_spin_out_path,
         C_stock_inits_yaml, C_stock_init_image, spatial_spin_fig_path, transient_covariate_path,
         transient_C_stock_hist, transient_flux_hist))
        
        RCTM.main(force_pft, point_mode, spin_years, run_transient, init_C_stocks_with_image, ranch,
         path_to_RCTM_params, path_to_RCTM_spatial_params, path_to_spin_covariates_point,
         path_to_spin_covariates_spatial, C_stock_spin_out_path_point, C_stock_spin_out_path,
         C_stock_inits_yaml, C_stock_init_image, spatial_spin_fig_path, transient_covariate_path,
         transient_C_stock_hist, transient_flux_hist)
     
#pool = Pool(2, maxtasksperchild=1)                      # Create a multiprocessing Pool
#pool.starmap(RCTM.main, args)
#pool.close()
#pool.join()