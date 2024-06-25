### Table SI

A list of open-source codes provided by this work for designing a targeted field sampling plan and generating high-resolution estimates of soil organic carbon stocks. Each step is numbered to indicate the corresponding deposited codes.

| Function<sup>a</sup> | Step | Platform<sup>b</sup> | Description | Optional |
| --- | --- | --- | --- | --- |
| Generate targeted points for initial soil sampling at case study sites | 1 | GEE | Extract pixel-based environmental covariates for preliminary sampling | Y |
|  | 2a | R | Determine initial sampling points based on environmental covariates | Y |
|  | 2b | R | Determine initial sampling points based on environmental covariates and distance to road | Y |
| Download and process legacy (KSSL) dataset | 3 | R | Download KSSL dataset for SOCC | N |
|  | 4 | GEE | Extract environmental covariates for KSSL dataset | N |
|  | 5 | R | Match KSSL measures with corresponding depth layers | N |
|  | 6 | R | Merge KSSL data with different sets of covariates | N |
|  | 7 | R | Obtain BD, SOCD, and SOCS datasets from KSSL | N |
| Generate high-resolution SOCC maps using different modeling methods | 8 | GEE | Extract environmental covariates for the sampling points of the study site | Y |
|  | 9 | R | Select point-based calibration dataset from KSSL using covariates | N |
|  | 10 | R | Select site-based calibration dataset from KSSL using covariates | Y |
|  | 11 | R | Build point and site-based models without using local samples for spiking | N |
|  | 12 | GEE | Extract pixel-based environmental covariates for the study site | N |
|  | 13 | R | Prepare pixel-based covariates for SOCC mapping | N |
|  | 14 | R | Carry out pixel-based modeling for SOCC mapping | N |
| Determine the number of local samples needed for a targeted field sampling design | 15 | R | Select different sets of local samples for model spiking and validation | N |
|  | 9 | R | Select pixel-based calibration set from KSSL using covariates | N |
|  | 10 | R | Select site-based calibration set from KSSL using covariates | Y |
|  | 16 | R | Build and compare models using different spiking and modeling schemes for determining the optimal number of local samples needed for model spiking | N |
| Determine targeted points for repeated soil sampling at case study sites | 12 | GEE | Extract pixel-based environmental covariates for the re-sampling campaign | N |
|  | 17 | R | Determine points for re-sampling campaign based on environmental covariates, site accessibility, and distance to roads | N |
| Model calibration and validation using independent samples | 15 | R | Select different sets of local samples for model spiking | Y |
|  | 18 | R | Select pixel-based calibration set from KSSL using covariates for SOCC, BD, and SOCS | N |
|  | 19 | R | Select site-based calibration set from KSSL using covariates for SOCC, BD, and SOCS | Y |
|  | 20 | R | Build and validate models for SOCC | N |
|  | 21 | R | Build and validate models for BD | N |
|  | 22 | R | Build and validate models for SOCS | N |
| Test model improvement strategies for SOCC and BD | 23 | R | Test the use of a subset of the covariates for SOCC model building | Y |
|  | 24 | R | Test the use of a subset of the covariates for BD model building | Y |
|  | 25 | R | Select point-based calibration dataset from broader-scale KSSL using covariates | Y |
|  | 26 | R | Test the use of broader-scale KSSL dataset for SOCC model building | Y |
|  | 27 | R | Test the use of broader-scale KSSL dataset for BD model building | Y |
|  | 28 | R | Select point-based calibration dataset by depth from KSSL using covariates | Y |
|  | 29 | R | Test the estimates of SOCC by building models separately for different depth layers | Y |
|  | 30 | R | Test the estimates of BD by building models separately for different depth layers | Y |
|  | 31 | R | Select point-based calibration dataset by depth from KSSL using subset of covariates | Y |
|  | 32 | R | Test the estimates of SOCC by building models separately for different depth layers and using subset of covariates for calibration set selection and modeling | Y |
|  | 33 | R | Test the estimates of BD by building models separately for different depth layers and using subset of covariates for calibration set selection and modeling | Y |
| Test alternative SOCS modeling strategies | 34 | R | Select point and site-based calibration sets from KSSL using covariates for SOCD | Y |
|  | 35 | R | Build and validate models for SOCD | Y |
|  | 36 | R | Compare different methods for modeling field and cluster-based SOCS, namely modeling SOCS directly, modeling SOCC and BD separately before combining them to estimate SOCS, and modeling SOCD first before integrating the results to estimate SOCS | Y |
|  | 37 | R | Examine performance for estimating SOCS using SoilGrid+ model outputs | Y |
| Generate high-resolution SOCS maps | 12 | GEE | Extract pixel-based environmental covariates for high-resolution SOCS mapping | N |
|  | 38 | R | Prepare site and pixel-based calibration sets for SOCS mapping | N |
|  | 39 | R | Carry out pixel-based modeling for SOCS mapping using local and KSSL samples | N |
|  | 40 | R | Carry out site-based modeling for SOCS mapping using local and KSSL samples | Y |
|  | 41 | R | Carry out SOCS mapping using only local samples | Y |
|  | 42 | R | Carry out pixel-based meta-modeling using local samples and KSSL samples | Y |

<sup>a</sup>BD = Soil bulk density; KSSL = Kellogg Soil Survey Laboratory; SOCC = soil organic carbon concentration; SOCD = soil organic carbon density; SOCS = soil organic carbon stocks.

<sup>b</sup>GEE = Google Earth Engine; R = R programming language.
