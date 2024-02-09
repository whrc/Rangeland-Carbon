/** 
 * utility functions for StarFM
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

# include "StarFM.h"

/**
 * display command line usage including parameters   
 */
void usage(char *command) {
  printf("Usage: %s <IN_PARAMETER_FILE>\n", command);
  printf("<IN_PARAMETER_FILE> is a text input file which contains:\n");
  printf("STARFM_PARAMETER_START\n");
  printf("\tNUM_IN_PAIRS = \n");
  printf("\tIN_PAIR_MODIS_FNAME = \n");
  printf("\tIN_PAIR_MODIS_MASK = (optional)\n");
  printf("\tIN_PAIR_LANDSAT_FNAME = \n");
  printf("\tIN_PAIR_LANDSAT_MASK = (optional)\n");
  printf("\tIN_PAIR_CLASSIFICATION_MAP = (optional)\n");
  printf("\tIN_PDAY_MODIS_FNAME = (accept mulitple days) \n");
  printf("\tIN_PDAY_MODIS_MASK = (optional)\n");
  printf("\tOUT_PDAY_LANDSAT_FNAME = (corresponding mulitple outputs) \n");
  printf("\tNROWS = \n");
  printf("\tNCOLS = \n");
  printf("\tSTART_ICOL = (optional)\n");
  printf("\tSTART_IROW = (optional)\n");
  printf("\tEND_ICOL = (optional)\n");
  printf("\tEND_IROW = (optional)\n");
  printf("\tRESOLUTION = \n");
  printf("\tSCALE_FACTOR = \n");
  printf("\tLANDSAT_FILLV = \n");
  printf("\tLANDSAT_DATA_RANGE = \n");
  printf("\tLANDSAT_UNCERTAINTY = \n");
  printf("\tMODIS_FILLV = \n");
  printf("\tMODIS_DATA_RANGE = \n");
  printf("\tMODIS_UNCERTAINTY = \n");
  printf("\tUSE_SPATIAL_FLAG = ON\n");
  printf("\tMAX_SEARCH_DISTANCE = \n");
  printf("\tNUM_SLICE_PURE_TEST = \n");
  printf("\tMIN_SAMPLE_PERCENT = (optional)\n");
  printf("\tPREDICTION_METHOD = 1 or 2 (1=weighting; 2=maximum occurence; default=1)\n"); 
  printf("\tUSE_LUT_IMPROVEMENT =\n"); 
  printf("STARFM_PARAMETER_END\n");
}


/** 
 * parse parameters from input file to variables
 */
int parseParameters(char *fname, SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par)
{
  int   i, k, total_pairs=0;
  char  buffer[MAX_STR_LEN] = "\0";
  char  *label = NULL;
  char  *tokenptr = NULL;
  char  *separator = "= ,";
  FILE  *in; 

  /* read from input parameter file */
  if((in=fopen(fname,"r"))==NULL) {
    printf("Can't open input %s\n", fname);
    return FAILURE;
  }

  fscanf(in, "%s", buffer);
  if(strcasecmp(buffer, "STARFM_PARAMETER_START") != 0) {
    printf("This is not a valid input file\n");
    return FAILURE;
  }

  /* set default values */
  par->USE_SPATIAL_FLAG = 1;
  par->USE_CLUSTER_MAP = 0;
  par->USE_LUT_IMPROVEMENT = 0;
  for(i=0; i<MAX_NPAIRS; i++) {
    strcpy(psensor[i]->landsat.mname, "none");
    strcpy(psensor[i]->modis.mname, "none");
    psensor[i]->start_irow = -1;
    psensor[i]->start_icol = -1;
    psensor[i]->end_irow = -1;
    psensor[i]->end_icol = -1;
  }
  /* set default minimum acceptable basis point */
  par->MIN_BASIS_POINT = 0;
  /* set default prediction method (weighting) */
  par->PREDICTION_METHOD = 1;

  par->NUM_PREDICTIONS = 0;
  /* process line by line */
  while(fgets(buffer, MAX_STR_LEN, in) != NULL) {

    if(strcasecmp(buffer, "STARFM_PARAMETER_END") == 0) break;
    
    /* get string token */
    tokenptr = strtok(buffer, separator);
    label=tokenptr;
    
    /* skip comment line */
    if(strcmp(label,"#") == 0) continue;

    while(tokenptr != NULL) {
      
      tokenptr = strtok(NULL, separator);
   
      /* extract parameters by comparing label */   
      if(strcasecmp(label, "NUM_IN_PAIRS") == 0) 
	par->NUM_PAIRS = atoi(tokenptr);
      else if(strcasecmp(label, "IN_PAIR_LANDSAT_FNAME") == 0)
	for(i=0; i<par->NUM_PAIRS; i++) {
	  sscanf(tokenptr, "%s", psensor[i]->landsat.fname);
	  tokenptr = strtok(NULL, separator);
	}
      else if(strcasecmp(label, "IN_PAIR_LANDSAT_MASK") == 0)
	for(i=0; i<par->NUM_PAIRS; i++) {
	  sscanf(tokenptr, "%s", psensor[i]->landsat.mname);
	  tokenptr = strtok(NULL, separator);
	}
      else if(strcasecmp(label, "IN_PAIR_MODIS_FNAME") == 0)
	for(i=0; i<par->NUM_PAIRS; i++) {
	  sscanf(tokenptr, "%s", psensor[i]->modis.fname);
	  tokenptr = strtok(NULL, separator);
	}      
      else if(strcasecmp(label, "IN_PAIR_MODIS_MASK") == 0)
	for(i=0; i<par->NUM_PAIRS; i++) {
	  sscanf(tokenptr, "%s", psensor[i]->modis.mname);
	  tokenptr = strtok(NULL, separator);
	}
      else if(strcasecmp(label, "IN_PAIR_CLASSIFICATION_MAP") == 0) {
	for(i=0; i<par->NUM_PAIRS; i++) {
	  sscanf(tokenptr, "%s", psensor[i]->cname);
	  tokenptr = strtok(NULL, separator);
	}     
	par->USE_CLUSTER_MAP = 1;
      }
      else if(strcasecmp(label, "IN_PDAY_MODIS_FNAME") == 0) {
	k = 0;
	do {
	  sscanf(tokenptr, "%s", psensor[par->NUM_PAIRS+k]->modis.fname);
	  tokenptr = strtok(NULL, separator);
	  k++;
	}  while(tokenptr != NULL);
	if(par->NUM_PREDICTIONS == 0)
	  par->NUM_PREDICTIONS = k;
	else if(k != par->NUM_PREDICTIONS) {
	  printf("\nnumber of IN_PDAY_MODIS_FNAME does not match OUT_PDAY_LANDSAT_FNAME\n");
	  return FAILURE;
	}
	total_pairs = par->NUM_PAIRS + par->NUM_PREDICTIONS;
      }
      else if(strcasecmp(label, "IN_PDAY_MODIS_MASK") == 0) {	
	k= 0;
	do {
	  sscanf(tokenptr, "%s", psensor[par->NUM_PAIRS+k]->modis.mname);
	  tokenptr = strtok(NULL, separator);
	  k++;
	}  while(tokenptr != NULL);
	if(par->NUM_PREDICTIONS == 0)
	  par->NUM_PREDICTIONS = k;
	else if(k != par->NUM_PREDICTIONS) {
	  printf("\nnumber of IN_PDAY_MODIS_MASK does not match IN_PDAY_MODIS_FNAME\n");
	  return FAILURE;
	}
	total_pairs = par->NUM_PAIRS+par->NUM_PREDICTIONS;
      }
      else if(strcasecmp(label, "OUT_PDAY_LANDSAT_FNAME") == 0) {
	k = 0;
	do {
	  sscanf(tokenptr, "%s", psensor[par->NUM_PAIRS+k]->landsat.fname);
	  tokenptr = strtok(NULL, separator);
	  k++;
	}  while(tokenptr != NULL);
	if(par->NUM_PREDICTIONS == 0)
	  par->NUM_PREDICTIONS = k;
	else if(k != par->NUM_PREDICTIONS) {
	  printf("\nnumber of OUT_PDAY_LANDSAT_FNAME does not match IN_PDAY_MODIS_FNAME\n");
	  return FAILURE;
	}
	total_pairs = par->NUM_PAIRS+par->NUM_PREDICTIONS;
      }
      else if(strcasecmp(label, "NROWS") == 0)
	for(i=0; i<MAX_NPAIRS; i++) { 
	  psensor[i]->landsat.nrows = atoi(tokenptr);
	  psensor[i]->modis.nrows = atoi(tokenptr);
	}
      else if(strcasecmp(label, "NCOLS") == 0)
	for(i=0; i<MAX_NPAIRS; i++) { 
	  psensor[i]->landsat.ncols = atoi(tokenptr);
	  psensor[i]->modis.ncols = atoi(tokenptr);
	}

      /* define sub-window for processing (11/2012) */
      else if(strcasecmp(label, "START_IROW") == 0)
	for(i=0; i<MAX_NPAIRS; i++)
	  psensor[i]->start_irow = atoi(tokenptr);
      else if(strcasecmp(label, "START_ICOL") == 0)
	for(i=0; i<MAX_NPAIRS; i++)
	  psensor[i]->start_icol = atoi(tokenptr);
      else if(strcasecmp(label, "END_IROW") == 0)
	for(i=0; i<MAX_NPAIRS; i++)
	  psensor[i]->end_irow = atoi(tokenptr);
      else if(strcasecmp(label, "END_ICOL") == 0)
	for(i=0; i<MAX_NPAIRS; i++)
	  psensor[i]->end_icol = atoi(tokenptr);

      else if(strcasecmp(label, "RESOLUTION") == 0)
	for(i=0; i<MAX_NPAIRS; i++) { 
	  psensor[i]->landsat.res = atof(tokenptr);
	  psensor[i]->modis.res = atof(tokenptr);
	}
      else if(strcasecmp(label, "SCALE_FACTOR") == 0)
	for(i=0; i<MAX_NPAIRS; i++) {
	  psensor[i]->landsat.scale = atof(tokenptr);    
	  psensor[i]->modis.scale = atof(tokenptr);
	}    
      
      else if(strcasecmp(label, "LANDSAT_FILLV") == 0)
	for(i=0; i<MAX_NPAIRS; i++) 
	  psensor[i]->landsat.fillv = atoi(tokenptr);    
      else if(strcasecmp(label, "LANDSAT_DATA_RANGE") == 0) {
	for(i=0; i<MAX_NPAIRS; i++)
	  psensor[i]->landsat.range[0] = atoi(tokenptr);    
	tokenptr = strtok(NULL, separator);
	for(i=0; i<MAX_NPAIRS; i++)
	  psensor[i]->landsat.range[1] = atoi(tokenptr);    	
      }
      else if(strcasecmp(label, "LANDSAT_UNCERTAINTY") == 0)
	for(i=0; i<MAX_NPAIRS; i++) 
	  psensor[i]->landsat.uncertainty = atof(tokenptr);    

      else if(strcasecmp(label, "MODIS_FILLV") == 0)
	for(i=0; i<MAX_NPAIRS; i++) 
	  psensor[i]->modis.fillv = atoi(tokenptr);    
      else if(strcasecmp(label, "MODIS_DATA_RANGE") == 0) {
	for(i=0; i<MAX_NPAIRS; i++)
	  psensor[i]->modis.range[0] = atoi(tokenptr);    
	tokenptr = strtok(NULL, separator);
	for(i=0; i<MAX_NPAIRS; i++)
	  psensor[i]->modis.range[1] = atoi(tokenptr);    	
      }
      else if(strcasecmp(label, "MODIS_UNCERTAINTY") == 0)
	for(i=0; i<MAX_NPAIRS; i++) 
	  psensor[i]->modis.uncertainty = atof(tokenptr);    

      else if(strcasecmp(label, "USE_SPATIAL_FLAG") == 0)
	par->USE_SPATIAL_FLAG = atoi(tokenptr);
      else if(strcasecmp(label, "USE_LUT_IMPROVEMENT") == 0)
	par->USE_LUT_IMPROVEMENT = atoi(tokenptr);
      else if(strcasecmp(label, "MAX_SEARCH_DISTANCE") == 0)
	par->MAX_SEARCH_DIS = atoi(tokenptr);
      else if(strcasecmp(label, "NUM_SLICE_PURE_TEST") == 0)
	par->NUM_CLS = atoi(tokenptr);
      else if(strcasecmp(label, "MIN_SAMPLE_PERCENT") == 0)
	/* convert to basis point */
	par->MIN_BASIS_POINT = atof(tokenptr) * 100.0;
      else if(strcasecmp(label, "PREDICTION_METHOD") == 0)
	par->PREDICTION_METHOD = atof(tokenptr);

      /* in case label (key words) is no the first word in a line */
      label = tokenptr;

    } /* while token */
  } /* while line */

  for(i=0; i<total_pairs; i++) {
    if(psensor[i]->start_irow == -1)
      psensor[i]->start_irow = 0;
    if(psensor[i]->start_icol == -1)
      psensor[i]->start_icol = 0;
    if(psensor[i]->end_irow == -1)
      psensor[i]->end_irow = psensor[i]->landsat.nrows-1;
    if(psensor[i]->end_icol == -1)
      psensor[i]->end_icol = psensor[i]->landsat.ncols-1;
  }

#ifdef DEBUG
  printf("\tSTARFM_PARAMETER_START\n");
  printf("\tNUM_IN_PAIRS = %d\n", par->NUM_PAIRS);
  printf("\tIN_PAIR_MODIS_FNAME = ");
  for(i=0; i<par->NUM_PAIRS; i++)
    printf("%s ", psensor[i]->modis.fname);
  printf("\n");
  printf("\tIN_PAIR_MODIS_MASK = ");
  for(i=0; i<par->NUM_PAIRS; i++)
    printf("%s ", psensor[i]->modis.mname);
  printf("\n");
  printf("\tIN_PAIR_LANDSAT_FNAME = ");
  for(i=0; i<par->NUM_PAIRS; i++)
    printf("%s ", psensor[i]->landsat.fname);
  printf("\n");
  printf("\tIN_PAIR_LANDSAT_MASK = ");
  for(i=0; i<par->NUM_PAIRS; i++)
    printf("%s ", psensor[i]->landsat.mname);
  printf("\n");
  if(par->USE_CLUSTER_MAP) {
    printf("\tIN_PAIR_CLASSIFICATION_MAP = ");
    for(i=0; i<par->NUM_PAIRS; i++)
      printf("%s ", psensor[i]->cname);
    printf("\n");
  }

  printf("\tIN_PDAY_MODIS_FNAME = ");
  for(k=par->NUM_PAIRS; k<total_pairs; k++)
    printf("%s ", psensor[k]->modis.fname);
  printf("\n\tIN_PDAY_MODIS_MASK = ");
  for(k=par->NUM_PAIRS; k<total_pairs; k++)
    printf("%s ", psensor[k]->modis.mname);
  printf("\n\tOUT_PDAY_LANDSAT_FNAME = ");
  for(k=par->NUM_PAIRS; k<total_pairs; k++)
    printf("%s ", psensor[k]->landsat.fname);

  printf("\n\tNROWS = %d\n", psensor[0]->landsat.nrows);
  printf("\tNCOLS = %d\n", psensor[0]->landsat.ncols);
  if(psensor[0]->start_irow!=0 || psensor[0]->start_icol!=0 || 
     psensor[0]->end_irow!=psensor[0]->landsat.nrows-1 || 
     psensor[0]->end_icol!=psensor[0]->landsat.ncols-1) {
    printf("\tSTART_ICOL = %d\n", psensor[0]->start_icol);
    printf("\tSTART_IROW = %d\n", psensor[0]->start_irow);
    printf("\tEND_ICOL = %d\n", psensor[0]->end_icol);
    printf("\tEND_IROW = %d\n", psensor[0]->end_irow);
  }
  printf("\tRESOLUTION = %5.1f\n", psensor[0]->landsat.res);
  printf("\tSCALE_FACTOR = %6.1f\n", psensor[0]->landsat.scale); 
  printf("\tLANDSAT_FILLV = %d\n", psensor[0]->landsat.fillv);
  printf("\tLANDSAT_DATA_RANGE = %d %d\n", psensor[0]->landsat.range[0], psensor[0]->landsat.range[1]);
  printf("\tLANDSAT_UNCERTAINTY = %5.1f\n", psensor[0]->landsat.uncertainty);
  printf("\tMODIS_FILLV = %d\n", psensor[0]->modis.fillv);
  printf("\tMODIS_DATA_RANGE = %d %d\n", psensor[0]->modis.range[0], psensor[0]->modis.range[1]);
  printf("\tMODIS_UNCERTAINTY = %5.1f\n", psensor[0]->modis.uncertainty);
  printf("\tUSE_SPATIAL_FLAG = %d\n", par->USE_SPATIAL_FLAG);
  printf("\tMAX_SEARCH_DISTANCE = %d\n", par->MAX_SEARCH_DIS);
  printf("\tNUM_SLICE_PURE_TEST = %d\n", par->NUM_CLS);
  printf("\tMIN_SAMPLE_PERCENT = %5.2f\n", par->MIN_BASIS_POINT/100.0);
  printf("\tPREDICTION_METHOD = %d\n", par->PREDICTION_METHOD);
  printf("\tUSE_LUT_IMPROVEMENT = %d\n", par->USE_LUT_IMPROVEMENT);
  printf("\tSTARFM_PARAMETER_END");

  printf("\n\nin_npairs=%d  npredicts=%d total_modis=%d\n", par->NUM_PAIRS, par->NUM_PREDICTIONS, total_pairs);
#endif

  fclose(in);
  return SUCCESS;
}



/**
 * initialize variables and allocate memory 
 */
int initializeVariables(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, BEST_PREDICTION **bplut)
{
  int i, j, r0;

  /* save low range for using look-up-table, just in case range[0] < 0 */
  r0 = psensor[0]->landsat.range[0];

  /* compute searching window size according to the search distance */
  par->WIN_SIZE = par->MAX_SEARCH_DIS * 2 / psensor[0]->landsat.res;

  /* compute maximum number of allowed candidate pixels for prediction */
  par->MAX_NUM_CAN = par->WIN_SIZE * par->WIN_SIZE * par->NUM_PAIRS + 1;

  /* initialize look-up-table to hold best predictions */
  for(i=psensor[0]->landsat.range[0]; i<psensor[0]->landsat.range[1]; i++)
    for(j=0; j<par->NUM_PAIRS; j++) {
      bplut[i-r0][j].psr = psensor[0]->landsat.fillv;
      bplut[i-r0][j].nsam = 0;
      bplut[i-r0][j].weight = 0.0;
    }

#ifdef DEBUG
  printf("\n\tWIN_SIZE = %d", par->WIN_SIZE);
  printf("\n\tMAX_NUM_CAN = %d", par->MAX_NUM_CAN);
#endif

  /* open input pairs */
  for(i=0; i<par->NUM_PAIRS+par->NUM_PREDICTIONS; i++) {

    if(i < par->NUM_PAIRS) {
      /* open Landsat reflectance file */
      if((psensor[i]->landsat.fp = fopen(psensor[i]->landsat.fname, "rb"))==NULL) {
	printf("Open Landsat input file %s error\n", psensor[i]->landsat.fname);
	return FAILURE;
      }
      /* open Landsat mask file */
      if((psensor[i]->landsat.mfp = fopen(psensor[i]->landsat.mname, "rb"))==NULL)
	printf("\nCannot open Landsat mask file %s, Landsat mask will not be used\n", psensor[i]->landsat.mname);

      /* open Landsat classification file if defined */
      if(par->USE_CLUSTER_MAP)
	if((psensor[i]->cfp = fopen(psensor[i]->cname, "rb"))==NULL) {
	  printf("Open Landsat classification file %s error\n", psensor[i]->cname);
	  return FAILURE;
	}
    }
    else {
      /* store predicted Landsat result */
      if((psensor[i]->landsat.fp = fopen(psensor[i]->landsat.fname, "wb"))==NULL) {
	printf("Open Landsat input file %s error\n", psensor[i]->landsat.fname);
	return FAILURE;
      }
    }
    
    /* allocate memories for Landsat */
    alloc_2dim_contig((void ***) (&psensor[i]->landsat.data), par->WIN_SIZE, 
		      psensor[0]->landsat.ncols, sizeof(short int));  
    alloc_2dim_contig((void ***) (&psensor[i]->landsat.mask), par->WIN_SIZE, 
		      psensor[0]->landsat.ncols, sizeof(char));  
    alloc_2dim_contig((void ***) (&psensor[i]->cls), par->WIN_SIZE, 
		      psensor[0]->landsat.ncols, sizeof(char));  

    /* open MODIS reflectance file */
    if((psensor[i]->modis.fp = fopen(psensor[i]->modis.fname, "rb"))==NULL) {
      printf("Open Landsat input file %s error\n", psensor[i]->modis.fname);
      return FAILURE;
    }
    /* open MODIS mask file */
    if((psensor[i]->modis.mfp = fopen(psensor[i]->modis.mname, "rb"))==NULL) 
      printf("Open MODIS mask file %s error, MODIS mask will not be used\n", psensor[i]->modis.mname);
    
    /* allocate memories for MODIS */
    alloc_2dim_contig((void ***) (&psensor[i]->modis.data),  par->WIN_SIZE, 
		      psensor[0]->modis.ncols, sizeof(short int));  
    alloc_2dim_contig((void ***) (&psensor[i]->modis.mask),  par->WIN_SIZE, 
		      psensor[0]->modis.ncols, sizeof(char));  
    
  }
  
  return SUCCESS;
}


/**
 * load first block of data to memory and keep current processing row in the middle of array
 * so the first half is filled fill values and the second half contains beginning rows
 */
void loadFirstBlock(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par)
{
  int ip;          /* index of input periods (0, NPAIRS) */
  int m, n;
  short int tmpint;
  unsigned char tmpchar, tmpcls=0;

  for(ip=0; ip<par->NUM_PAIRS+par->NUM_PREDICTIONS; ip++) {
  
    /* load block of rows (par->WIN_SIZE) from all inputs into memory */
    for(m=0; m<par->WIN_SIZE; m++) 
      for(n=0; n<psensor[0]->landsat.ncols; n++) { 
	/* read MODIS SR (NUM_PAIRS) */
	if(m >= par->WIN_SIZE/2) {
	  /* read second half */ 
	  fread(&tmpint, sizeof(short int), 1, psensor[ip]->modis.fp);
	  if(psensor[ip]->modis.mfp != NULL)
	    fread(&tmpchar, sizeof(char), 1, psensor[ip]->modis.mfp);
	  else
	    /* assume all valid if no mask file */
	    tmpchar = VALID;
	}
	else {
	  /* fill first half */
	  tmpint = psensor[0]->modis.fillv;
	  tmpchar = INVALID;
	}
	psensor[ip]->modis.data[m][n] = tmpint;
	psensor[ip]->modis.mask[m][n] = tmpchar;

	/* read Landsat SR (NUM_PAIRS - 1), as last one is used for prediction */
	if(ip < par->NUM_PAIRS) {
	  if(m >= par->WIN_SIZE/2) { 
	    fread(&tmpint, sizeof(short int), 1, psensor[ip]->landsat.fp);
	    /* read mask data */
	    if(psensor[ip]->landsat.mfp != NULL)
	      fread(&tmpchar, sizeof(char), 1, psensor[ip]->landsat.mfp);
	    else
	      tmpchar = VALID;
	    /* read classification data */
	    if(par->USE_CLUSTER_MAP)
	      fread(&tmpcls, sizeof(char), 1, psensor[ip]->cfp);
	  }
	  else {
	    tmpint = psensor[0]->landsat.fillv;
	    tmpchar = INVALID;
	  }
	  psensor[ip]->landsat.data[m][n] = tmpint;
	  psensor[ip]->landsat.mask[m][n] = tmpchar;
	  psensor[ip]->cls[m][n] = tmpcls;	  
	}
      }
  }
}


/**
 * load next processing line to memory
 * first move the current array one line up 
 * then load next image row to the end line of array
 */  
void loadNextRow(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, int irow)
{
  int ip;          /* index of input periods (0, NPAIRS) */
  int m, n;
  short int tmpint;
  unsigned char tmpchar, tmpcls=0;

  /* read next processing line */
  for(ip=0; ip<par->NUM_PAIRS+par->NUM_PREDICTIONS; ip++) {
    
    for(m=0; m<par->WIN_SIZE; m++)
      
      for(n=0; n<psensor[0]->landsat.ncols; n++) { 
	
	/* move one row up in data array */	      
	if(m < par->WIN_SIZE-1) {
	  psensor[ip]->modis.data[m][n] = psensor[ip]->modis.data[m+1][n];
	  psensor[ip]->modis.mask[m][n] = psensor[ip]->modis.mask[m+1][n];
	  if(ip < par->NUM_PAIRS) {
	    psensor[ip]->landsat.data[m][n] = psensor[ip]->landsat.data[m+1][n];
	    psensor[ip]->landsat.mask[m][n] = psensor[ip]->landsat.mask[m+1][n];
	    psensor[ip]->cls[m][n] = psensor[ip]->cls[m+1][n];
	  }
	}
	/* then read next line to the end of array */
	else {
	  /* read MODIS SR (NUM_PAIRS) */
	  if(irow + par->WIN_SIZE/2 < psensor[0]->modis.nrows) {
	    /* read MODIS reflectance */
	    fread(&tmpint, sizeof(short int), 1, psensor[ip]->modis.fp);
	    /* read MODIS mask */
	    if(psensor[ip]->modis.mfp != NULL)
	      fread(&tmpchar, sizeof(char), 1, psensor[ip]->modis.mfp);
	    else
	      tmpchar = VALID;
	  }
	  else {
	    tmpint = psensor[0]->modis.fillv;
	    tmpchar = INVALID;
	  }
	  psensor[ip]->modis.data[m][n] = tmpint;
	  psensor[ip]->modis.mask[m][n] = tmpchar;

	  /* read Landsat SR (NUM_PAIRS - 1) */
	  if(ip < par->NUM_PAIRS) {
	    if(irow + par->WIN_SIZE/2 < psensor[ip]->landsat.nrows) {
	      /* read Landsat reflectance */
	      fread(&tmpint, sizeof(short int), 1, psensor[ip]->landsat.fp);
	      /* read Landsat mask */
	      if(psensor[ip]->landsat.mfp != NULL)
		fread(&tmpchar, sizeof(char), 1, psensor[ip]->landsat.mfp);
	      else
		tmpchar = VALID;
	      /* read Landsat classification data */
	      if(par->USE_CLUSTER_MAP)
		fread(&tmpcls, sizeof(char), 1, psensor[ip]->cfp);
	    }
	    else {
	      tmpint = psensor[0]->landsat.fillv;
	      tmpchar = INVALID;
	    }
	    psensor[ip]->landsat.data[m][n] = tmpint;
	    psensor[ip]->landsat.mask[m][n] = tmpchar;
	    psensor[ip]->cls[m][n] = tmpcls;
	  }
	    
	} /* end if( m < par->WIN_SIZE-1) */
      } /* end of icol (n) */
  } /* end loop for all input pairs */
}


/**
 * write ENVI header for two outputs (prediction and quality) 
 */
int writeHeader(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par)
{
  int i, k;
  FILE *fp;
  char tmpname[200];

  for(k=0; k<par->NUM_PREDICTIONS; k++) {
    /* write ENVI header for predicted surface reflectance */  
    sprintf(tmpname, "%s.hdr", psensor[par->NUM_PAIRS+k]->landsat.fname);
    if((fp=fopen(tmpname, "w"))==NULL) {
      printf("\nOpen header file %s error!", tmpname);
      return FAILURE;
    }
    
    fprintf(fp, "ENVI\n");
    fprintf(fp, "description = {\n"); 
    fprintf(fp, "predicted fine resolution reflectance for %s\n", psensor[par->NUM_PAIRS+k]->modis.fname);
    fprintf(fp, "using following image pairs:\n");
    for(i=0; i<par->NUM_PAIRS; i++)
      fprintf(fp, "%s %s\n", psensor[i]->landsat.fname, psensor[i]->modis.fname);
    fprintf(fp, "}\n");

    fprintf(fp, "samples = %d\n", psensor[0]->end_icol-psensor[0]->start_icol+1);
    fprintf(fp, "lines = %d\n",   psensor[0]->end_irow-psensor[0]->start_irow+1);
    fprintf(fp, "x start = %d\n", psensor[0]->start_icol+1);
    fprintf(fp, "y start = %d\n", psensor[0]->start_irow+1);

    fprintf(fp, "bands = 1\n");
    fprintf(fp, "header offset = 0\n");
    fprintf(fp, "file type = ENVI Standard\n");
    fprintf(fp, "data type = 2\n");
    fprintf(fp, "interleave = bsq\n");
    //MPG
    //fprintf(fp, "map info = {Geographic Lat/Lon, 1, 1, -114.3904801329673404, 46.9917708275762962, 0.000269494585235856, 0.000269494585235856,WGS-84}\n");

    //BGR
    //fprintf(fp, "map info = {Geographic Lat/Lon, 1, 1, -105.239791491284, 47.4493726333068, 0.000269494585235856, 0.000269494585235856,WGS-84}\n");
    
    //HPGRS
    //fprintf(fp, "map info = {Geographic Lat/Lon, 1, 1, -104.8773212741416074, 41.2186578226537748, 0.000269494585235856, 0.000269494585235856,WGS-84}\n");
    
    //S66
    //fprintf(fp, "map info = {Geographic Lat/Lon, 1, 1, -100.9187153116121181, 46.7715937514386013, 0.000269494585235856, 0.000269494585235856,WGS-84}\n");
    
    //fprintf(fp, "coordinate system string = {GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]}\n");
    fprintf(fp, "sensor type = Unknown\n");
    fprintf(fp, "byte order = 0\n");
    fprintf(fp, "wavelength units = Unknown\n");
    fprintf(fp, "data ignore value = 0");
    fclose(fp);

    /* write ENVI header for quality (basis point used within searching distance) */  
    /*sprintf(tmpname, "%s.psam.hdr", psensor[par->NUM_PAIRS+k]->landsat.fname);
    if((fp=fopen(tmpname, "w"))==NULL) {
      printf("\nOpen header file %s error!", tmpname);
      return FAILURE;
    }
  
    fprintf(fp, "ENVI\n");
    fprintf(fp, "description = {basis point (0-10000) of samples used within searching distance\n");
    fprintf(fp, "32767 means replacement using the highest psam}\n"); 
    fprintf(fp, "samples = %d\n", psensor[0]->end_icol-psensor[0]->start_icol+1);
    fprintf(fp, "lines = %d\n",   psensor[0]->end_irow-psensor[0]->start_irow+1);
    fprintf(fp, "x start = %d\n", psensor[0]->start_icol+1);
    fprintf(fp, "y start = %d\n", psensor[0]->start_irow+1);
    fprintf(fp, "bands = 1\n");
    fprintf(fp, "header offset = 0\n");
    fprintf(fp, "file type = ENVI Standard\n");
    fprintf(fp, "data type = 2\n");
    fprintf(fp, "interleave = bsq\n");
    fprintf(fp, "sensor type = Unknown\n");
    fprintf(fp, "byte order = 0\n");
    fprintf(fp, "wavelength units = Unknown\n");

    fclose(fp);*/
  }

  return SUCCESS;
}



/**
 * free memory and close file pointer
 */
void cleanup(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, BEST_PREDICTION **bplut)
{
  int i;

  /* close input pairs and output */
  for(i=0; i<par->NUM_PAIRS+par->NUM_PREDICTIONS; i++) {
    fclose(psensor[i]->landsat.fp);
    fclose(psensor[i]->modis.fp);
    if(psensor[i]->landsat.mfp != NULL && i<par->NUM_PAIRS) 
      fclose(psensor[i]->landsat.mfp);
    if(psensor[i]->modis.mfp != NULL) 
      fclose(psensor[i]->modis.mfp);  
    if(par->USE_CLUSTER_MAP && i<par->NUM_PAIRS)
      fclose(psensor[i]->cfp);
  }

  /* free allocated memory */
  for(i=0; i<par->NUM_PAIRS+par->NUM_PREDICTIONS; i++) {
    free_2dim_contig((void **) psensor[i]->landsat.data);
    free_2dim_contig((void **) psensor[i]->modis.data);
    free_2dim_contig((void **) psensor[i]->landsat.mask);
    free_2dim_contig((void **) psensor[i]->modis.mask);
    free_2dim_contig((void **) psensor[i]->cls);
  }
  
  for(i=0; i<MAX_NPAIRS; i++)
    free(psensor[i]);
  free(par);
  free_2dim_contig((void **) bplut);

}


