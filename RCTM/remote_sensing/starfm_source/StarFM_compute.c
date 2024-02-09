/** 
 * computation functions for StarFM
 *   Revision 1.2.2   07/2017   Feng Gao
 *     - fixed a bug in processing MODIS fillvalue
 *
 *   Revision 1.2.1   10/2014   Feng Gao
 *     - enabled OpenMP for multi-thread processing   
 *
 *   Revision 1.2.0   06/2014   Feng Gao
 *     - modified to allow multiple-date predictions 
 *     - does not support highest occurences
 *     - does not support LUT improvement
 *
 *    Revision 1.1.x  11/2012  Feng Gao added processing for subwindow 
 *    Revision 1.1    01/2008  Feng Gao
 */

#include "StarFM.h"

/**
 * do statistics on the input Landsat data and
 * find min, max, mean and stdev and then
 * compute slice distance based on the defined
 * number of classes (NUM_CLS) for homogeneity testing
 */
void doStatistics(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par) 
{
  int ip, irow, icol;
  unsigned char tmpchar;
  short int tmpint, min, max;
  double sumx, sumx2, num;
 
  /* compute mean and stdev for each input Landsat SR file */
  for(ip=0; ip<par->NUM_PAIRS; ip++) {

#ifdef DEBUG
    printf("\n\tfor input %s ", psensor[ip]->landsat.fname);
#endif

    min = psensor[ip]->landsat.range[1];
    max = psensor[ip]->landsat.range[0];
    sumx = 0.0;
    sumx2 = 0.0;
    num = 0.0;
    for(irow=0; irow<psensor[ip]->landsat.nrows; irow++) {
      for(icol=0; icol<psensor[ip]->landsat.ncols; icol++) {
	 fread(&tmpint, sizeof(short int), 1, psensor[ip]->landsat.fp);
	 if(psensor[ip]->landsat.mfp != NULL) 
	   fread(&tmpchar, sizeof(char), 1, psensor[ip]->landsat.mfp);
	 else
	   /* assume all valid if no mask file */
	   tmpchar = VALID;

	 if(irow<psensor[ip]->start_irow||irow>psensor[ip]->end_irow||icol<psensor[ip]->start_icol||icol>psensor[ip]->end_icol) continue;

	 /* only check valid value */
	 if(tmpchar == 1 && tmpint > psensor[ip]->landsat.range[0] && tmpint < psensor[ip]->landsat.range[1]) {
	   if(tmpint > max) max = tmpint;
	   if(tmpint < min) min = tmpint;	  
	   sumx += tmpint;
	   sumx2 += tmpint * tmpint;
	   num++;
	 } /* end of valid pixels */	
      } /* end of icol */
    } /* end of irow */
    
    /* compute basic statistics */
    psensor[ip]->landsat.min = min;
    psensor[ip]->landsat.max = max;
    psensor[ip]->landsat.mean  = sumx/num;
    psensor[ip]->landsat.stdev = sqrt(sumx2/num - (sumx/num) * (sumx/num));
    /* compute slice value for each input in 2 stdev * 2 sides / num_slice / 2 */
    psensor[ip]->landsat.slice_value = (int)(psensor[ip]->landsat.stdev * 2.0 / par->NUM_CLS + 0.5);

#ifdef DEBUG
    printf("\n\tmin=%4d  max=%4d  mean=%4d  stdev=%4d  slice_dis=%4d", 
	   psensor[ip]->landsat.min, psensor[ip]->landsat.max, 
	   psensor[ip]->landsat.mean, psensor[ip]->landsat.stdev, psensor[ip]->landsat.slice_value);
#endif    

    rewind(psensor[ip]->landsat.fp);
    if(psensor[ip]->landsat.mfp != NULL)
      rewind(psensor[ip]->landsat.mfp);
  } 
}


/**
 * do prediction based on the StarFM algorithm (Gao et al., 2006)
 * write prediction and sample basis point (1/10000) to outputs
 * save best prediction to look-up-table for later improvement
 */
int doPrediction(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, BEST_PREDICTION **bplut) 
{
  int i, j, k;
  short int **row_psr;                  /* predicted Landsat SR in short int (for store to file) */

  /*FILE *fp[MAX_NPAIRS];
    char tmpname[200];*/

  /* open qa file for write */
  /* for(k=0; k<par->NUM_PREDICTIONS; k++) {
     sprintf(tmpname, "%s.psam", psensor[par->NUM_PAIRS+k]->landsat.fname);
     if((fp[k]=fopen(tmpname, "wb"))==NULL) {
     printf("\nOpen QA file %s error!", tmpname);
     return FAILURE;
     }
     }*/

  alloc_2dim_contig ((void ***) (&row_psr), psensor[0]->landsat.ncols, par->NUM_PREDICTIONS, sizeof(short int));
  
  loadFirstBlock(psensor, par);

  printf("\n\tprocessing irow ... ");
  
  for(i=0; i<psensor[0]->landsat.nrows; i++) {

    if(i<psensor[0]->start_irow||i>psensor[0]->end_irow) continue;
    printf("%5d\b\b\b\b\b", i+1);

#pragma omp parallel for shared(psensor, par)
    for(j=0; j<psensor[0]->landsat.ncols; j++) {

      /* process each pixel (i,j) */
      if(j<psensor[0]->start_icol||j>psensor[0]->end_icol) continue;
      doOnePixel(psensor, par, i, j, row_psr[j]);

    } /* end of icol loop */

    /* save results for one row after all parallel computings end */
    for(j=0; j<psensor[0]->landsat.ncols; j++)
      for(k=0; k<par->NUM_PREDICTIONS; k++)
	fwrite(&(row_psr[j][k]), sizeof(short int), 1, psensor[par->NUM_PAIRS+k]->landsat.fp);

    loadNextRow(psensor, par, i);
    
  } /* end of irow (i) */ 
  
  free_2dim_contig((void **)row_psr);
  /*for(k=0; k<par->NUM_PREDICTIONS; k++)
    fclose(fp[k]);*/

  return SUCCESS;
}


/**
 * do prediction for one pixel based on the StarFM algorithm (Gao et al., 2006)
 * separated for OpenMP processing
 */
int doOnePixel(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, int i, int j, short int *all_psr) 
{
  int ican;         /* index of selected candidates */
  int rel_n;        /* relative n (index of column within searching window) */
  int ip;           /* index of input periods (0, NPAIRS) */
  int k, m, n, modv, lndv;
  int nsam[MAX_NPAIRS];
  short int diff;
  short int psam;         /* basis point (1/10000) of samples used within searching distance */
  int total_candidates;   /* total number of potential Landsat and MODIS SR pairs */
  int usable_nsam;        /* total number of usable samples */
  int is_similar;         /* 1 means a spectral similar type */               
  float s_dis, t_dis, dis_x, dis_y, sum_reverse_dis;  /* temporary variables used to compute spatial and temporal distance */
  float fsr;                                  /* final predicted Landsat SR in float */
  short int psr;                              /* predicted Landsat SR in short int (for store to file) */
  short int min_r_msr[MAX_NPAIRS], min_r_lsr; /* minimum SR differences on r_msr and r_lsr */
  int r_msr_uncertain, r_lsr_uncertain;       /* uncertainty of SR differences (assume independent to each other) */
  int msr_uncer1, msr_uncer2, lsr_uncer;      /* temporary variable to compute uncertainty */
  int combined_uncertain;
  CANDIDATE_PIXEL candidate[MAX_NSAMPLES];

  /* decide surface reflectance uncertainties from QA and AOT */
  /* here I just use predefined values */
  msr_uncer1 = psensor[0]->modis.uncertainty;
  msr_uncer2 = psensor[0]->modis.uncertainty;
  lsr_uncer = psensor[0]->landsat.uncertainty;

  /* compute uncertainty for the combined variables (assume independent) */
  r_msr_uncertain = sqrt(msr_uncer1*msr_uncer1 + msr_uncer2*msr_uncer2);
  r_lsr_uncertain = sqrt(msr_uncer1*msr_uncer1 + lsr_uncer*lsr_uncer);
  combined_uncertain = sqrt(r_msr_uncertain*r_msr_uncertain +  r_lsr_uncertain* r_lsr_uncertain); 
      
  ican = 0;

  /* save working pixel to index 0, all other candidates start from 1 */
  candidate[ican].loc[0] = j;
  candidate[ican].loc[1] = i;
  for(k=0; k<par->NUM_PREDICTIONS; k++)
    /* save MODIS data for later replacement if no prediction (temparory storage) */ 	  
    candidate[ican].r_msr[k] = psensor[par->NUM_PAIRS+k]->modis.data[par->WIN_SIZE/2][j];
  candidate[ican].lsr = psensor[0]->landsat.data[par->WIN_SIZE/2][j];
  ican++;

  for(k=0; k<par->NUM_PREDICTIONS; k++) 
    min_r_msr[k] = psensor[0]->modis.range[1];
  min_r_lsr = psensor[0]->landsat.range[1];
      
  for(ip=0; ip<par->NUM_PAIRS; ip++) {

    /* save current location to the head of array, ip=0 to ican=1 and ip=1 to ican=2 */
    candidate[ican].loc[0] = j;
    candidate[ican].loc[1] = i;
    candidate[ican].msr = psensor[ip]->modis.data[par->WIN_SIZE/2][j];
    candidate[ican].lsr = psensor[ip]->landsat.data[par->WIN_SIZE/2][j];
    candidate[ican].mask = psensor[ip]->modis.mask[par->WIN_SIZE/2][j] 
      * psensor[ip]->landsat.mask[par->WIN_SIZE/2][j];

    for(k=0; k<par->NUM_PREDICTIONS; k++) {
      candidate[ican].r_msr[k] = psensor[ip]->modis.data[par->WIN_SIZE/2][j] -
	psensor[par->NUM_PAIRS+k]->modis.data[par->WIN_SIZE/2][j];
      if(abs(candidate[ican].r_msr[k]) < min_r_msr[k]) min_r_msr[k] = abs(candidate[ican].r_msr[k]);
    }

    candidate[ican].r_lsr = psensor[ip]->modis.data[par->WIN_SIZE/2][j] -
      psensor[ip]->landsat.data[par->WIN_SIZE/2][j];
    if(abs(candidate[ican].r_lsr) < min_r_lsr) min_r_lsr = abs(candidate[ican].r_lsr);
    
    ican++;      
  }

  if(par->USE_SPATIAL_FLAG) {
    /* if use spatial info, then search all pixels include current one, reset index to 1 */
    ican = 1;
    
    /* search small window and find potential candidates */
    for(ip=0; ip<par->NUM_PAIRS; ip++) {
      
      nsam[ip] = 0;
      for(m=0; m<par->WIN_SIZE; m++)
	for(n=0; n<par->WIN_SIZE; n++) {
	  
	  rel_n = n - par->WIN_SIZE/2 + j;   /* convert to actual icol */
	      
	  /* valid image range test */
	  if(rel_n < 0 || rel_n >= psensor[0]->landsat.ncols) 
	    continue;
	      
	  /* check MODIS data */ 
	  modv = psensor[ip]->modis.data[m][rel_n];
	  lndv = psensor[ip]->landsat.data[m][rel_n];
	  if( modv == psensor[0]->modis.fillv ||
	      modv < psensor[0]->modis.range[0] ||
	      modv > psensor[0]->modis.range[1])
	    continue;

	  /* check Landsat data */
	  if(lndv == psensor[0]->landsat.fillv ||
	     lndv < psensor[0]->landsat.range[0] ||
	     lndv > psensor[0]->landsat.range[1] )
	    continue;
		
	  is_similar = 0;
	  if(par->USE_CLUSTER_MAP) {
	    /* use defined classification map */
	    if(psensor[ip]->cls[m][rel_n] == psensor[ip]->cls[par->WIN_SIZE/2][j])
	      is_similar = 1;
	  } 
	  else {
	    /* do spectral similar test if classification map is not defined */
	    diff = abs(psensor[ip]->landsat.data[m][rel_n] - psensor[ip]->landsat.data[par->WIN_SIZE/2][j]); 		
	    /* if spectral similar - I just use a single band slice test here */
	    if(diff < psensor[ip]->landsat.slice_value) 			
	      is_similar = 1;
	  }
		
	  if(is_similar) {
	    /* put to the candidate list */
	    candidate[ican].loc[0] = rel_n;
	    candidate[ican].loc[1] = m - par->WIN_SIZE/2 + i;
	    candidate[ican].msr = psensor[ip]->modis.data[m][rel_n];
	    candidate[ican].lsr = psensor[ip]->landsat.data[m][rel_n];
	    candidate[ican].mask = psensor[ip]->modis.mask[m][rel_n] * psensor[ip]->landsat.mask[m][rel_n];

	    for(k=0; k<par->NUM_PREDICTIONS; k++) {
	      modv = psensor[par->NUM_PAIRS+k]->modis.data[m][rel_n];
	      if(modv == psensor[0]->modis.fillv || 
		 modv < psensor[0]->modis.range[0] || modv > psensor[0]->modis.range[1])
		candidate[ican].r_msr[k] = psensor[0]->modis.fillv;
	      else	  
		candidate[ican].r_msr[k] = psensor[ip]->modis.data[m][rel_n] - modv;
	    }
	    candidate[ican].r_lsr = candidate[ican].msr - candidate[ican].lsr;
	    
	    ican++;		
	    nsam[ip]++;

	  } /* endif spectral similar */ 

	} /* end of moving window searching */
      
    } /* end of foreach input pair */

  } /* end of use spatial (neighbor) information */

  total_candidates = ican;

  /* do predictions for all MODIS dates (new feature v1.2.0) */
  for(k=0; k<par->NUM_PREDICTIONS; k++) {

    /* data valid check */
    for(m=1; m<total_candidates; m++) {      
      candidate[m].valid = 1;      
      /* data range checking */
      modv =  candidate[m].r_msr[k];
      /* accepts the small differences within one stdev */
      /* added check on MODIS data (Feng 06/23/2017) */
      if((abs(candidate[m].r_lsr) > (min_r_lsr+r_lsr_uncertain) &&
	  abs(modv) > (min_r_msr[k]+r_msr_uncertain)) 
	 || modv == psensor[0]->modis.fillv || candidate[m].mask == INVALID)
	candidate[m].valid = 0;
	  
      /*if( candidate[m].mask == INVALID || modv == psensor[0]->modis.fillv 
	|| abs(modv) > (min_r_msr[k]+r_msr_uncertain)) 
	candidate[m].valid = 0;*/
    }      
	
    /* do prediction using weighting approach (default) */
    sum_reverse_dis = 0.0;	
    usable_nsam = 0;
    for(m=1; m<total_candidates; m++) {
	
      if(candidate[m].valid == 1) {
	    
	/* compute distance to the central pixel */
	dis_x = candidate[m].loc[0]-candidate[0].loc[0];
	dis_y = candidate[m].loc[1]-candidate[0].loc[1];
	s_dis = sqrt(dis_x*dis_x + dis_y*dis_y);
	
	/* consider a small uncertainty for the difference of reflectance */
	/* to avoid zero from either of them (and so give a biased large weight) */
	if(par->NUM_PAIRS == 1)	    
	  /* if input only contains one pair of Landsat and MODIS, then only use 
	     the difference between MODIS and Landsat - Feng (10/21/08) */
	  t_dis = (abs(candidate[m].r_lsr) + 1);
	else
	  t_dis = (abs(candidate[m].r_msr[k]) + 1) * (abs(candidate[m].r_lsr) + 1);

	if(t_dis < combined_uncertain)
	  /* very close (less than uncertainty), gets maximum weight */
	  candidate[m].r_dis = 1.0;	      
	else 
	  /* otherwise, compute combined weight */
	  candidate[m].r_dis = 1.0/(t_dis*(1.0+s_dis * psensor[0]->landsat.res / par->MAX_SEARCH_DIS));
	  
	sum_reverse_dis += candidate[m].r_dis;
	usable_nsam++;
	    
      }
    }
      
    if(usable_nsam == 0)
      /* use MODIS SR if no good candidate */
      /*fsr = candidate[0].r_msr[k];*/ 	  
      /* use fill value */
      fsr = psensor[0]->landsat.fillv;
    else {
      /* initialize value for sum */
      fsr = 0.0;      
      for(m=1; m<total_candidates; m++) {
	/* normalize weight */
	if(candidate[m].valid == 1)
	  candidate[m].weight = candidate[m].r_dis/sum_reverse_dis;
	else
	  candidate[m].weight = 0.0;     
	/* predict reflectance using weighting function */
	fsr += candidate[m].weight * (candidate[m].lsr - candidate[m].r_msr[k]);	  
      }
    } /* end of prediction using weighting approach */
    
    /* convert to integer */
    if(fsr > psensor[0]->landsat.range[0] && fsr < psensor[0]->landsat.range[1])
      psr = (int)(fsr + 0.5);
    else
      psr = psensor[0]->landsat.fillv;

    /* compute basis point of samples used (keep precision with basis point instead of percentage) */
    psam = (short int)((10000.0*usable_nsam)/(par->WIN_SIZE*par->WIN_SIZE*par->NUM_PAIRS) + 0.5);

    /* save prediction for k date */
    all_psr[k] = psr;

#ifdef DEBUG
    if(i==DEBUG_irow && j==DEBUG_icol) {
      printf("\nfile=%s icol=%d irow=%d",  psensor[par->NUM_PAIRS+k]->modis.fname, candidate[0].loc[0], candidate[0].loc[1]);
      printf("\npsr=%f", fsr);
      printf("\nmodis=%d min_r_msr=%d  min_r_lsr=%d prediction=%d use_nsam=%d psam=%d(basis_point)",
	     candidate[0].r_msr[k], min_r_msr[k], min_r_lsr, psr, usable_nsam, psam);
      printf("\nr_msr_uncertain=%d r_lsr_uncertain=%d combined_uncertain=%d\n", 
	     r_msr_uncertain, r_lsr_uncertain, combined_uncertain);
      for (m=1; m<total_candidates; m++)
	/*if(candidate[m].valid == 1) */
	printf("%3d x=%4d y=%4d modis=%4d diff_modis=%4d landsat=%4d weight=%6.4f mask=%d valid=%d\n", 
	       m, candidate[m].loc[0], candidate[m].loc[1], candidate[m].msr, candidate[m].r_msr[k], 
	       candidate[m].lsr, candidate[m].weight, candidate[m].mask, candidate[m].valid);
    }
#endif

  } /* end of k loop for each MODIS date */	      

  return SUCCESS;
}

