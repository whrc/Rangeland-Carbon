/** 
 *   header file for StarFM
 *   Revision 1.2.1   10/2014   Feng Gao
 *     - enabled OpenMP for multi-thread processing   
 *
 *   Revision 1.2.0   06/2014   Feng Gao
 *     - modified to allow multiple-date predictions 
 *     - does not support highest occurences
 *     - does not support LUT improvement
 *
 *   Revision 1.1   01/2008   Feng Gao
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <omp.h>

#define MAX_STR_LEN 100000
#define SUCCESS     0
#define FAILURE     -1
#define VALID   1
#define INVALID 0

#define MAX_NPAIRS  365      /* maximum number of input and output pairs */
#define MIN_TWEIGHT 0.0005   /* minimum acceptable sum of weight */
#define MAX_NSAMPLES 10000   /* maximum number of neighboring pixels with similar spectral signal */

/* data struct for surface reflectance */
typedef struct {
  char fname[MAX_STR_LEN];   /* data file name */
  char mname[MAX_STR_LEN];   /* mask file name */
  int  nrows;                /* number of rows */
  int  ncols;                /* number of columns */
  int  fillv;                /* fill value used */
  int  range[2];             /* sensor data range */
  FILE *fp;                  /* file pointer for data */
  FILE *mfp;                 /* file pointer for mask */
  float res;                 /* spatial resolution */
  float scale;               /* scale factor */
  float uncertainty;         /* data accuracy */
  short int **data;          /* sensor data: i=irow; j=icol */
  unsigned char **mask;      /* mask data: 1=use(VALID) 0=not use(INVALID); i=irow; j=icol */
  short int min;             /* minimum value for ipair */ 
  short int max;             /* maximum value for ipair */
  short int mean;            /* mean value of valid data */
  short int stdev;           /* standard deviation of valid data */
  short int slice_value;     /* computed slice value based on the (max-min)/NUM_CLS */
} SENSOR;


/* struct for input image pair (Landsat and MODIS) */
typedef struct {
  SENSOR landsat;           /* for finer resolution data */
  SENSOR modis;             /* for coarse resolution data */
  char cname[MAX_STR_LEN];  /* classification file name */
  FILE *cfp;                /* file pointer to classification map */
  unsigned char **cls;      /* cluster map: i=irow; j=icol */
  int start_irow, end_irow;
  int start_icol, end_icol;
} SENSOR_PAIR;


/* struct for control parameters */
typedef struct {
  int USE_SPATIAL_FLAG;      /* 1=use spatial information; 0=not use */ 
  int USE_CLUSTER_MAP;       /* 1=use defined cluster map if defined; 0=not use if not defined */
  int MAX_SEARCH_DIS;        /* maximum searching distance in the same unit as resolution */
  int NUM_CLS;               /* number of slice for homogeneity testing */
  int NUM_PAIRS;             /* number of input data pair */ 
  int NUM_PREDICTIONS;       /* number of input prediction dates */
  int WIN_SIZE;              /* window size determined by the MAX_SEARCH_DIS */
  int MAX_NUM_CAN;           /* maximum number of candidate pixels determined by WIN_SIZE */ 
  int MIN_BASIS_POINT;       /* minimum acceptable basis point (1/10000) within searching window */
  int PREDICTION_METHOD;     /* 1=weighting; 2=maximum frequency */
  int USE_LUT_IMPROVEMENT;   /* 1=use LUT value to improve poor prediction; 0=not use (added on Nov 2011) */
} CONTROL_PARAMETER;


/* struct for the selected prediction samples */
typedef struct {
  short int loc[2];        /* sample location */
  short int msr;           /* MODIS surface reflectance */
  short int lsr;           /* Landsat surface reflectance */
  short int r_lsr;         /* relative difference between Landsat and MODIS */
  short int r_msr[MAX_NPAIRS];        /* relative difference between two MODIS data */
  short int valid;
  unsigned char mask;
  float r_dis;             /* reversed distance */
  float weight;            /* sample normalized weight */
} CANDIDATE_PIXEL;


/* struct for look-up-table to hold the best prediction */
typedef struct {
  short int psr;       /* predicted surface reflectance */
  short int nsam;      /* samples used within searching window in basis point (*10000) */
  float weight;        /* total weight (before normalize) for this prediction */
} BEST_PREDICTION;


/* functions in StarFM_util.c */
void usage(char *command);
void loadFirstBlock(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par);
void loadNextRow(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, int irow);
void cleanup(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, BEST_PREDICTION **bplut);
int  parseParameters(char *fname, SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par);
int  initializeVariables(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, BEST_PREDICTION **bplut);
int  writeHeader(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par);

/* functions in StarFM_compute.c */
void doStatistics(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par);
int  doPrediction(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, BEST_PREDICTION **bplut);
int  improvePrediction(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, BEST_PREDICTION **bplut);
void smoothLUT(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, BEST_PREDICTION **bplut);
int  doOnePixel(SENSOR_PAIR *psensor[], CONTROL_PARAMETER *par, int i, int j, short int *all_psr) ;

/* functions in StarFM_alloc.c */
void alloc_1dim_contig (void **, int, int);
void alloc_2dim_contig (void ***, int, int, int);
void alloc_3dim_contig (void ****, int, int, int, int);
void alloc_4dim_contig (void *****, int, int, int, int, int);
void free_2dim_contig  (void **);
void free_3dim_contig  (void ***);
void free_4dim_contig  (void ****);

/* use DEBUG mode to debug specific pixel */
/*#define DEBUG*/
/*#define DEBUG_LUT*/
#ifdef DEBUG
  #define DEBUG_icol 10
  #define DEBUG_irow 10
#endif
