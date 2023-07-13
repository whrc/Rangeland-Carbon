/** 
 * memory allocation functions for StarFM
 * revised from MODIS code
 * Revision 1.0   02/2005  Feng Gao
 */

#include "StarFM.h"

void alloc_1dim_contig (void **ptr, int d1, int elsize)
{
   void *p = NULL;

   p = calloc (d1, elsize);
   if (!p) {
      printf ("Memory allocation error in alloc_1dim_contig\n");
      exit(1);
   }
   *ptr = p;
   return;
}


void alloc_2dim_contig (void ***ptr, int d1, int d2, int elsize)
{
   void *p = NULL;
   void **pp = NULL;
   int i = 0;

   /* alloc array for data */
   alloc_1dim_contig ((void **) (&p), d1 * d2, elsize);

   /* alloc array for pointers */
   alloc_1dim_contig ((void **) (&pp), d1, sizeof (void *));

   /* Set the pointers to indicate the appropriate elements of the data array. */
   for (i = 0; i < d1; i++) {
      pp[i] = (char *) p + (i * d2 * elsize);
   }

   *ptr = pp;

   return;
}


void alloc_3dim_contig (void ****ptr, int d1, int d2, int d3, int elsize)
{
   void *p = NULL;
   void **pp = NULL;
   void ***ppp = NULL;
   int i = 0;

   /* allocate the data array */
   alloc_1dim_contig ((void **) (&p), d1 * d2 * d3, elsize);

   /* alloc the double pointers */
   alloc_1dim_contig ((void **) (&pp), d1 * d2, sizeof (void *));

   /* and again for the first dimensional pointers */
   alloc_1dim_contig ((void **) (&ppp), d1, sizeof (void **));

   /* first set the d1 pointers */
   for (i = 0; i < d1; i++) {
      ppp[i] = pp + (i * d2);
   }

   /* next set all of the d2 pointers */
   for (i = 0; i < d1 * d2; i++) {
      pp[i] = (char *) p + (i * d3 * elsize);
   }

   *ptr = ppp;

   return;
}


void alloc_4dim_contig (void *****ptr, int d1, int d2, int d3, int d4, int elsize)
{
   void *p = NULL;
   void **pp = NULL;
   void ***ppp = NULL;
   void ****pppp = NULL;
   int i = 0;

   /* allocate the data array */
   alloc_1dim_contig ((void **) (&p), d1 * d2 * d3 * d4, elsize);

   /* alloc the double pointers */
   alloc_1dim_contig ((void **) (&pp), d1 * d2 * d3, sizeof (void *));

   /* and again for the triple pointers */
   alloc_1dim_contig ((void **) (&ppp), d1 * d2, sizeof (void **));

   alloc_1dim_contig ((void **) (&pppp), d1, sizeof (void ***));

   for (i = 0; i < d1; i++) {
      pppp[i] = ppp + (i * d2);
   }

   for (i = 0; i < d1 * d2; i++) {
      ppp[i] = pp + (i * d3);
   }

   for (i = 0; i < d1 * d2 * d3; i++) {
      pp[i] = (char *) p + (i * d4 * elsize);
   }

   *ptr = pppp;

   return;
}


void free_2dim_contig (void **a)
{
  free (a[0]);
  free (a);
  return;
} 

void free_3dim_contig (void ***a)
{
   free (a[0][0]);
   free (a[0]);
   free (a);
   return;
}

void free_4dim_contig (void ****a)
{
   free (a[0][0][0]);
   free (a[0][0]);
   free (a[0]);
   free (a);
   return;
}



