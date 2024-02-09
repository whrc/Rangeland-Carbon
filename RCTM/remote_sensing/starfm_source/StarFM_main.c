/**   
 * ! Description
 *
 *   This program predicts Landsat surface reflectance based on the STARFM algorithm
 *   by Gao et al. (2006). It blends Landsat and MODIS surface reflectance by fusing 
 *   high-frequency temporal information from MODIS and high spatial resolution 
 *   information from Landsat. Note that this approach has several reasonable assumptions.
 *   Please make sure your inputs satisfy the following conditions:  
 *     1. atmospherically corrected surface reflectances are assumed to be comparable 
 *   from time to time and location to location. 
 *     2. similar adjacent land cover areas are assumed to have similar spectral patterns 
 *   and temporal change patterns over a small area. 
 *     3. surface reflectance of a homogeneous land cover type is assumed to be identical 
 *   for both coarse and fine spatial resolution.
 * 
 * ! Input
 *
 *   - Landsat and MODIS input pairs in surface reflectance 
 *   - MODIS surface reflectance on the prediction date
 *   - prediction control parameters
 *
 * ! Ouput
 *
 *   - Landsat surface reflectance on the prediction date
 *   - prediction quality saves in basis point (0-10000) of samples used for prediction 
 *
 * ! Usage
 *
 *   StarFM.exe input.txt
 *   (where input.txt contains all necessary input files and control parameters)
 * 
 * ! Reference
 *
 *   Gao, F., J. Masek, M. Schwaller and H. Forrest, (2006). On the Blending of the 
 *   Landsat and MODIS Surface Reflectance: Predict Daily Landsat Surface Reflectance, 
 *   IEEE Transactions on Geoscience and Remote Sensing, Vol. 44, No. 8, pp. 2207-2218.
 *
 * ! Revision History
 *
 *   Revision 1.2.1   10/2014   Feng Gao
 *     - enabled OpenMP for multi-thread processing   
 *
 *   Revision 1.2.0   06/2014   Feng Gao
 *     - modified to allow multiple-date predictions 
 *     - does not support highest occurences
 *     - does not support LUT improvement
 *
 *   Revision 1.1.3   01/2009   Feng Gao
 *     - added new prediction approach based on the samples with highest occurences 
 *       The new approach can be used to handle discontinuous patches such as cropland 
 *       While weighting approach is good for large continous pathches such as forest    
 *
 *   Revision 1.1.2   10/2008   Feng Gao
 *     - improved weighting algorithm. If input contains only one pair of Landsat and MODIS, then
 *       just use spectral (S_ijk) and spatial (d_ijk) differences since temporal difference (T_ijk)
 *       doesn't carry appropriate information as it was first designed for the multiple input pairs.
 *
 *   Revision 1.1.1   03/2008   Feng Gao
 *     - smoothed best predictions in the look-up-table to be used for replacing poor prediction
 *     - allowed to deal with negative input data range
 *     - fixed a bug in creating candidates from multiple input pairs under no spatial info option
 *     - constrained candidates by excluding invalid MODIS and Landsat inputs  
 *     - constrained cases in the look-up-table for improving poor quality prediction
 *     - changed MIN_SAMPLE_BASIS_POINT to MIN_SAMPLE_PERCENT (input parameter) for easier use
 *
 *   Revision 1.1     01/2008     Feng Gao   
 *     public release version 
 *     - completely restructured code 
 *     - simplified processing based on per band 
 *     - computed and saved prediction quality (in basis point of samples used)
 *     - allowed defining cloud and poor quality data mask to exclude data from prediction 
 *     - checked and improved poor prediction ( NEW to the paper )
 *       (if MIN_SAMPLE_BASIS_POINT is defined in control parameters)
 *       this is useful for a small object to utlize temporal information from 
 *       a same type large object that is beyond the searching distance
 *     - changed to accept binary input (in 2 bytes per pixel) instead of HDF format
 * 
 *   Revision 1.0     02/2005     Feng Gao
 *     original version (research version) 
 */

#include "StarFM.h"

int main(int argc, char *argv[])
{
  int i;

  SENSOR_PAIR *psensor[MAX_NPAIRS];
  CONTROL_PARAMETER *par;
  BEST_PREDICTION **bplut;

  /* verify command line usage */
  if(argc!=2) {
    usage(argv[0]);
    exit(1);
  }

  /* allocate memory for variable */
  for(i=0; i<MAX_NPAIRS; i++)
    if(!(psensor[i] = malloc(sizeof(SENSOR_PAIR)))) exit(1);
  if(!(par = malloc(sizeof(CONTROL_PARAMETER)))) exit(1);

  /* retrieve input parameters */
  printf("parse parameters");
  if(parseParameters(argv[1], psensor, par) == FAILURE) {
    printf("Retrieve input parameters error! \n");
    usage(argv[0]);
    exit(1);
  }

  
  /* allocate memory for look-up-table that will hold the best prediction */
  alloc_2dim_contig((void ***) (&bplut),  psensor[0]->landsat.range[1]-psensor[0]->landsat.range[0], 
		    par->NUM_PAIRS, sizeof(BEST_PREDICTION));  
  
  /* initialize variables and allocate memory for working */
  printf("\ninitialize variables");
  initializeVariables(psensor, par, bplut);
  
  /* only do statistics if spatial information will be used */
  if(par->USE_SPATIAL_FLAG) {
    /* do statistics on the input Landsat (finer resolution) data */      
    printf("\ndoing statistics on Landsat ...");
    doStatistics(psensor, par);
  }

  /* do prediction based on data pair and the daily MODIS observation */
  printf("\ndoing prediction ...");
  if(doPrediction(psensor, par, bplut) == FAILURE) {
    printf("\nDoing prediction error\n");
    exit(1);
  }
 
  /* write ENVI header for output */
  printf("\nwrite ENVI header file for output");
  if(writeHeader(psensor, par) == FAILURE) {
    printf("\nwrite ENVI file header error");
    exit(1);
  }

  /* cleanup memory */
  printf("\nfree memories and close files");
  cleanup(psensor, par, bplut);

  printf("\nfinished successfully\n");
  return SUCCESS;
}


