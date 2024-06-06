Description (v1.2.1)
====================

MODIS and Landsat surface reflectance products have complementary characteristics 
in terms of spatial and temporal resolution. To fully exploit these individual 
datasets, a Spatial and Temporal Adaptive Reflectance Fusion Model (STARFM) was 
developed by Gao et al. (2006). The STARFM approach combines the spatial resolution of 
Landsat with the temporal frequency of MODIS. STARFM uses comparisons of one or more 
pairs of observed Landsat/MODIS maps, collected on the same day, to predict maps at 
Landsat-scale on other MODIS observation dates. This program implements the STARFM 
approach using C. It has been tested in Linux system using gcc compiler 
(require v4.4 or newer in order to support multiprocessing OpenMP 3.0). 

changes for RCTM:
- fixed minor bug by uncommenting line 396 in StarFM_util.c.

v1.1.2: 
- original version (first release)
- accepts one or two image pairs
- one prediction per run

v1.2.1: 
- includes data fusion for a sub-window (defined in input file)
- allows combining of multiple predictions together if the pair images are the same
- implements parallel computing on the multiple cores system using OpenMP 3.0
  (remove the OpenMP option if you don't want to use multiple cores)
- disabled LUT and highest occurrence option    

Reference
  
Gao, F., Masek, J., Schwaller, M., and Hall, F. On the Blending of the 
Landsat and MODIS Surface Reflectance: Predict Daily Landsat Surface Reflectance, 
IEEE Transactions on Geoscience and Remote Sensing. 44(8):2207-2218. 2006.

Gao, F., Hilker, T., Zhu, X., Anderson, M. A., Masek, J., Wang, P. and Yang, Y. 
Fusing Landsat and MODIS data for vegetation monitoring, IEEE Geoscience and Remote 
Sensing Magazine. 3(3):47-60. doi: 10.1109/MGRS.2015.2434351. 2015.


Installation
============

extract package: 
> tar -xvzf StarFM.tar.gz

StarFM/source
StarFM.h
StarFM_main.c
StarFM_compute.c
StarFM_util.c
StarFM_alloc.c
Makefile

StarFM/simulation_test
StarFM/simulation_test/change_veg  
StarFM/simulation_test/road_veg  
StarFM/simulation_test/water_veg
(include simulated testing files)

StarFM/reflectance_test (full version only)
include MODIS and Landsat testing files)

then run 
> make 
under source directory it will create an executive program "StarFM.exe"
(OpenMP can be disabled by removing "-fopenmp" option in the Makefile)

Usage 
=====

> StarFM.exe input.txt
where input.txt contains the input files and parameters (see details below)


Inputs
======

All Landsat and MODIS data should be first pre-processed and co-registered and saved in 
  - same resolution   (Landsat resolution)
  - same image size   
  - same map projection
  - same image extent (orthorectifed and precisely registered)
  - same scale factor (10000 for MODIS reflectance products)
  - "short int" type (2 byte/pixel) for the scaled reflectance data
  - "unsigned char" type (1 byte/pixel) for optional mask and classification file 

Major inputs include:
  - Landsat and MODIS surface reflectance image pair in binary format
  - MODIS surface reflectance for the prediction date
  - optional Landsat classification map for input pair 
  - optional Landsat and MODIS cloud and poor quality data mask
  - control parameters 

Detailed input parameters and explanations:

STARFM_PARAMETER_START
# number of input pairs, maximum = 2
        NUM_IN_PAIRS = 
# input MODIS files
        IN_PAIR_MODIS_FNAME = 
# optional mask file, same order MODIS files and use NONE if not exist
        IN_PAIR_MODIS_MASK = 
# input Landsat files in the same order as MODIS
        IN_PAIR_LANDSAT_FNAME = 
# optional mask file for Landsat
        IN_PAIR_LANDSAT_MASK = 
# optional classification map in the same order as input pairs
# classification map has to be generated from input Landsat data 
# with as many separable classes as possible
        IN_PAIR_CLASSIFICATION_MAP =
# MODIS input for the prediction date 
        IN_PDAY_MODIS_FNAME = (accepts multiple dates MODIS files separated by space)
# optional mask file for MODIS
        IN_PDAY_MODIS_MASK = 
# output Landsat prediction file
        OUT_PDAY_LANDSAT_FNAME = (include same number of MODIS prediction dates) 
# number of rows 
        NROWS =
# number of columns
        NCOLS =
# define data fusion processing window (default is entire image)
  	START_IROW =
	START_ICOL =
	END_IROW =
	END_ICOL = 
# spatial resolution
        RESOLUTION =
# scale factor for input reflectance file 
        SCALE_FACTOR =
# fill value for Landsat surface reflectance
        LANDSAT_FILLV =
# Landsat data range
        LANDSAT_DATA_RANGE =
# uncertainty in the same unit as the scaled reflectance
        LANDSAT_UNCERTAINTY = 
# fill value for MODIS surface reflectance
        MODIS_FILLV =
# MODIS data range
        MODIS_DATA_RANGE =
# uncertainty for MODIS 
        MODIS_UNCERTAINTY =
# spatial information flag, "ON" is strongly suggested
        USE_SPATIAL_FLAG = ON
# maximum search distance for the spectral similar neighbor pixels 
        MAX_SEARCH_DISTANCE =
# number of slice for the spectral similar test (pure pixel)
        NUM_SLICE_PURE_TEST =
STARFM_PARAMETER_END

(use # for comments)

MODIS inputs can choose daily surface reflectance products (MOD09) or nadir BRDF-adjusted 
reflectance (NBAR) (MOD43). They can be reprojected and resampled to Landsat projection 
and resolution using MODIS Reprojection Tool (MRT). The MRT can be downloaded from 
http://edcdaac.usgs.gov/landdaac/tools/modis/index.asp 


Output
======

Landsat surface reflectance for the prediction date in binary format. 
(comes with an associated ENVI header file.)


Testing Data
============

1) simulated data fusion

- StarFM/simulation_test/water_veg (see Fig. 2 in Gao et al.)

"input_use_spatial.txt" 
use t1 image pair to predict reflectance for t2 using spatial information. 

"input_no_spatial.txt"  
use t1 image pair to predict reflectance for t2 without using spatial information. 

This example shows the difference with and without using spatial information. 
The "USE_SPATIAL_FLAG" is set to "ON" for the rest of tests. We strongly recommend
to use spatial information in the prediction.  

- StarFM/simulation_test/change_veg (see Fig. 3 in Gao et al.)

"input_t2.txt" 
uses t1 and t4 pairs to predict t2

"input_t3.txt" 
uses t1 and t4 pairs to predict t3

This example shows how StarFM algorithm handles changing objects.

- StarFM/simulation_test/road_veg (see Fig. 5 in Gao et al.)

"input_t2.txt" 
uses t1 and t4 image pairs to predict t2 

"input_t3.txt" 
uses t1 and t4 image pairs to predict t3

This example shows how StarFM handles linear objects.

2) Actual Landsat data fusion 
(full version only, compact version need to download separately)

"input.06-04-01.*.txt"
uses 05-24-01 image pair to predict 06-04-01 Landsat surface reflectance
(used default slice approach for spectral similar testing)

"input.07-11-01.*.txt"
uses 05-24-01 image pair to predict 07-11-01 Landsat surface reflectance
(used defined classification map for spectral similar testing)

"input2.07-11-01.*.txt"
uses 05-24-01 and 08-12-01 image pairs to predict 07-11-01
(removed replacement option, slight different from v1.1.2)

"inputm.pair_05-24-01.*.txt"
uses 05-24-01 image pair to predict Landsat surface reflectance for later days
on 06-04-01, 07-11-01 and 08-12-01 in one combined run 
(results are included for code verification only)


Appendix A
==========

1) MODIS data preparation
MODIS inputs can be either daily surface reflectance (MOD09) or 16-day nadir BRDF-adjusted
surface reflectance (NBAR) product (MOD43B4 or MCD43A4). They need to be pre-processed
for STARFM inputs.
a) order MODIS tiles that cover whole Landsat scene on the Landsat acquisition date and
prediction dates
b) use MODIS reprojection tool (MRT) to mosaic and reproject MODIS tiles to Landsat 
projection, resolution and extents. Save otuputs in binary format.

2) Landsat data preparation
Landsat DN values need to be calibrated and atmosphericly corrected and saved in the
same format as MODIS surface reflectance (2-byte integer with scale factor of 10000) 
with one file per band. 
Note that MODIS bands use different band number. You can relate them with 
Landsat  MODIS
   1       3
   2       4
   3       1
   4       2
   5       6
   7       7

 




