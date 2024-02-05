import starfm_postprocessing
import os

roi_file='/home/amullen/Rangeland-Carbon/res/site_footprints/HLD_tiles.txt'

with open(roi_file) as f:
  sites = [line.rstrip('\n') for line in f]
  
  args=[]
  
  for site in sites:
    ranch='/'.join(site.split('/')[-2:])
                                
    args = {'starfm_in_dir': f'Ranch_Runs/{ranch}/starfm_v2/',
            'covariate_in_dir': f'Ranch_Runs/{ranch}/covariates/',
            'out_dir': f'Ranch_Runs/{ranch}/covariates_nc/',
            'NLCD_in_dir': f'Ranch_Runs/{ranch}/landcover/NLCD_2019.tif',
            'RAP_in_dir': f'Ranch_Runs/{ranch}/landcover/RAP_2019.tif',
            'param_file': '/home/amullen/Rangeland-Carbon/modeling/RCTM_params.yaml',
            'bucket_name': 'rangelands'}
    
    covariate_nc_outname = 'covariates.nc'
    covariate_csv_outname = 'covariates.csv'
    spin_nc_outname = 'covariates_spin.nc'
    spin_csv_outname = 'covariates_spin.csv'
    param_outname = args['out_dir'] + 'params.tif'
            
    ds = starfm_postprocessing.gen_covariates(args['starfm_in_dir'], args['covariate_in_dir'],  args['out_dir'], covariate_nc_outname, args['bucket_name'], '2002-01-01', '2023-01-01', gap_fill = True)
    #average indices
    df = starfm_postprocessing.image_average_variables(ds, ['ndvi','srad','vpd','tsoil','sm1','sm2','shortwave_radition','tavg','tmin','prcp','clay'], plot_dir=os.path.join(args['out_dir'], 'transient_figs/'))
    df.to_csv('gs://' + args['bucket_name'] + '/' + args['out_dir'] + covariate_csv_outname)
    #gen spinup
    spin_ds = starfm_postprocessing.aggregate_for_spinup(ds, args['out_dir'], spin_nc_outname, '2002-01-01', '2005-12-31', args['bucket_name'],  period = 5)
    df = starfm_postprocessing.image_average_variables(spin_ds, ['ndvi','srad','vpd','tsoil','sm1','sm2','shortwave_radition','tavg','tmin','prcp','clay'], plot_dir=os.path.join(args['out_dir'], 'spin_figs/'))
    df.to_csv('gs://' + args['bucket_name'] + '/' + args['out_dir'] + spin_csv_outname)
    #get spatial params
    starfm_postprocessing.get_spatial_RCTM_params(args['NLCD_in_dir'], args['RAP_in_dir'], args['param_file'], param_outname, args['bucket_name'], param_type='starfm')
    