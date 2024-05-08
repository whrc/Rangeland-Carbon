# Development Code Table

A list of open-source codes provided by this work for developing the Rangeland Carbon Tracking and Monitoring (RCTM) system. Each step is numbered to indicate the corresponding deposited codes.

| Function<sup>a</sup>                             | Step | Platform<sup>b</sup> | Description                                                                                       | Optional |
|--------------------------------------|------|----------|---------------------------------------------------------------------------------------------------|----------|
| Extract land cover information for flux tower sites | 1    | GEE      | Extract land cover information for all relevant Ameriflux and NEON sites to identify grasslands or rangelands | N        |
| Download and convert flux tower data to a daily time step | 2a   | R        | Example for processing typical Ameriflux site data and visualize results                           | N        |
|                                      | 2b   | R        | Examples for processing Ameriflux site data originally provided in Matlab or netCDF formats     | Y        |
|                                      | 3    | R        | Batch downloading NEON datasets                                                                 | N        |
|                                      | 4    | R        | Batch processing NEON datasets                                                                   | N        |
|                                      | 5    | R        | Example for converting NEON data to daily time step and visualize results                         | N        |
| Quality control of flux tower datasets | 6    | R        | Example for conducting quality control for GPP measurements                                      | N        |
|                                      | 7    | R        | Example for conducting quality control for NEE measurements                                      | N        |
| Extract environmental covariates for flux tower sites | 8a   | GEE      | Extract environmental covariates for GPP model                                                   | N        |
|                                      | 8b   | GEE      | Example for extracting environmental covariates for GPP model using alternative climate inputs due to incomplete temporal data coverage | Y        |
|                                      | 9a   | GEE      | Extract environmental covariates for the NEE model                                               | N        |
|                                      | 9b   | GEE      | Example for extracting environmental covariates for NEE model using alternative climate inputs due to incomplete temporal data coverage | Y        |
|                                      | 10   | GEE      | Extract environmental covariates model spin up                                                    | N        |
| Combine datasets for vegetation type-based model calibration and validation | 11   | Colab    | Combine quality-controlled flux tower measurements-derived GPP and covariates for all sites     | N        |
|                                      | 12   | Colab    | Combine quality-controlled flux tower measurements of NEE and covariates for all sites            | N        |
|                                      | 13   | Colab    | Combine covariate datasets for model spin up                                                      | N        |
|                                      | 14   | R        | Combine remote sensing inputs with environmental covariates and measured data                     | Y        |
| Conduct GPP model calibration and validation | 15   | R        | Fit GPP model with observations                                                                  | N        |
|                                      | 16   | R        | Model calibration and validation for each individual site                                         | N         |
|                                      | 17   | R        | Model calibration and validation using cross-validation                                            | N         |
|                                      | 18   | R        | Calculate error metrics associated with cumulative GPP                                             | N         |
| Conduct NEE model calibration and validation | 19   | R        | Initialize carbon pools to reduce computation time                                                 | Y        |
|                                      | 20   | R        | Fit NEE model with observations                                                                  | N        |
|                                      | 21   | R        | Model calibration and validation for each individual site                                         | N         |
|                                      | 22   | R        | Model calibration and validation using cross-validation                                            | N         |
|                                      | 23   | R        | Calculate error metrics associated with cumulative NEE                                             | N         |
| SOC model validation                 | 24   | R        | Validate model performance for estimating SOC stocks                                               | N         |

<sup>a</sup>STARFM: Spatial and Temporal Adaptive Reflectance Fusion Model; GPP: Gross primary Productivity; NEE: net ecosystem exchange of CO2; SOC: soil organic carbon.
<sup>b</sup>GEE: Google Earth Engine; Colab: Google Colaboratory.

