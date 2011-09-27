static char rcsid[] = "$Id: chrsegment.c,v 1.6 2005/06/21 18:39:08 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "chrsegment.h"

#include <stdlib.h>
#include <string.h>		/* For strlen */
#include <math.h>
#include "bool.h"
#include "mem.h"
#include "intlist.h"
#include "nr-x.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define PVALUE_THRESHOLD 1e-10
#define MINEVAL 500000U
#define MININTERVAL 200000U
#define MINCOUNT 5
#define NPSEUDO 0


#define T Chrsegment_T
typedef struct T *T;
struct T {
  int start;
  int end;

  int parent_start;
  int parent_end;

  unsigned int left_startpos;
  unsigned int left_endpos;
  unsigned int middle_startpos;
  unsigned int middle_endpos;
  unsigned int right_startpos;
  unsigned int right_endpos;

  double pvalue;
  double mean_left;
  double mean_middle;
  double mean_right;

  T left;
  T middle;
  T right;
};


static void
Chrsegment_free (T *old) {
  T left, middle, right;

  if (*old != NULL) {
    right = (*old)->right;
    Chrsegment_free(&right);
    middle = (*old)->middle;
    Chrsegment_free(&middle);
    left = (*old)->left;
    Chrsegment_free(&left);
    FREE(*old);
  }
  return;
}


static T
Chrsegment_new (int start, int end, int parent_start, int parent_end,
		double pvalue, double mean_left, double mean_middle, 
		double mean_right, unsigned int *chrpositions) {
  T new = (T) MALLOC(sizeof(*new));

  new->start = start;
  new->end = end;
  new->parent_start = parent_start;
  new->parent_end = parent_end;

  new->left_startpos = chrpositions[parent_start];
  new->left_endpos = chrpositions[start-1];
  
  new->middle_startpos = chrpositions[start];
  new->middle_endpos = chrpositions[end-1];

  new->right_startpos = chrpositions[end];
  new->right_endpos = chrpositions[parent_end-1];

  new->pvalue = pvalue;
  new->mean_left = mean_left;
  new->mean_middle = mean_middle;
  new->mean_right = mean_right;

  new->left = (T) NULL;
  new->middle = (T) NULL;
  new->right = (T) NULL;

  return new;
}

static void
Chrsegment_print (FILE *fp, T this) {
  if (this != NULL) {
    if (this->left == NULL) {
      /* If this->start == this->parent_start, signifies no breakpoint */
      if (this->start != this->parent_start) {
	fprintf(fp,"%g %d..%d %u..%u %f\n",
		this->pvalue,this->parent_start,this->start-1,
		this->left_startpos,this->left_endpos,this->mean_left);
      }
    } else {
      Chrsegment_print(fp,this->left);
    }

    if (this->middle == NULL) {
      /* If this->end == this->start, signifies a single breakpoint */
      if (this->end != this->start) {
	fprintf(fp,"%g %d..%d %u..%u %f\n",
		this->pvalue,this->start,this->end-1,
		this->middle_startpos,this->middle_endpos,this->mean_middle);
      }
    } else {
      Chrsegment_print(fp,this->middle);
    }

    if (this->right == NULL) {
      /* If this->parent_end == this->end, signifies no breakpoint */
      if (this->parent_end != this->end) {
	fprintf(fp,"%g %d..%d %u..%u %f\n",
		this->pvalue,this->end,this->parent_end-1,
		this->right_startpos,this->right_endpos,this->mean_right);
      }
    } else {
      Chrsegment_print(fp,this->right);
    }
  }
  return;
}

static void
Chrsegment_dump (Uintlist_T *startpositions, Uintlist_T *endpositions, Doublelist_T *means,
		 T this) {
  if (this != NULL) {
    if (this->left == NULL) {
      /* If this->start == this->parent_start, signifies no breakpoint */
      if (this->start != this->parent_start) {
	*startpositions = Uintlist_push(*startpositions,this->left_startpos);
	*endpositions = Uintlist_push(*endpositions,this->left_endpos);
	*means = Doublelist_push(*means,this->mean_left);
      }
    } else {
      Chrsegment_dump(&(*startpositions),&(*endpositions),&(*means),this->left);
    }

    if (this->middle == NULL) {
      /* If this->end == this->start, signifies a single breakpoint */
      if (this->end != this->start) {
	*startpositions = Uintlist_push(*startpositions,this->middle_startpos);
	*endpositions = Uintlist_push(*endpositions,this->middle_endpos);
	*means = Doublelist_push(*means,this->mean_middle);
      }
    } else {
      Chrsegment_dump(&(*startpositions),&(*endpositions),&(*means),this->middle);
    }

    if (this->right == NULL) {
      /* If this->parent_end == this->end, signifies no breakpoint */
      if (this->parent_end != this->end) {
	*startpositions = Uintlist_push(*startpositions,this->right_startpos);
	*endpositions = Uintlist_push(*endpositions,this->right_endpos);
	*means = Doublelist_push(*means,this->mean_right);
      }
    } else {
      Chrsegment_dump(&(*startpositions),&(*endpositions),&(*means),this->right);
    }
  }
  return;
}

static void
Chrsegment_breakpoints (Intlist_T *breakpoints, T this) {
  if (this != NULL) {
    if (this->left == NULL) {
      /* If this->start == this->parent_start, signifies no breakpoint */
      if (this->start != this->parent_start) {
	*breakpoints = Intlist_push(*breakpoints,this->start);
      }
    } else {
      Chrsegment_breakpoints(&(*breakpoints),this->left);
    }

    if (this->middle == NULL) {
      /* If this->end == this->start, signifies a single breakpoint */
      if (this->end != this->start) {
	*breakpoints = Intlist_push(*breakpoints,this->end);
      }
    } else {
      Chrsegment_breakpoints(&(*breakpoints),this->middle);
    }

    if (this->right == NULL) {
      /* If this->parent_end == this->end, signifies no breakpoint */
      if (this->parent_end != this->end) {
	*breakpoints = Intlist_push(*breakpoints,this->parent_end);
      }
    } else {
      Chrsegment_breakpoints(&(*breakpoints),this->right);
    }
  }
  return;
}


static double
interval_detect (int *best_start, int *best_end, 
		 double *best_meanx_left, double *best_meanx_middle, double *best_meanx_right,
		 double *sumx, double *sumxx, unsigned int *chrpositions, int start, int end) {
  double pvalue = 1.0, fscore, tstat, min_ssr_sep;
  int pos, pos1, pos2;
  double sumx_pseudo, meanx, meanx_outside, meanx_inside, meanx_left, meanx_right,
    ssr_outside, ssr_inside, ssr_left, ssr_right, ssr_sep, ssr;
  int n, n_outside, n_inside, n_left, n_right;

  *best_start = start;
  *best_end = end;

  if ((n = end - start) < MINCOUNT) {
    return 1.0;
  }
  meanx = (sumx[end] - sumx[start])/(double) n;
  min_ssr_sep = ssr = (sumxx[end] - sumxx[start]) - (double) n*meanx*meanx;
  debug(printf("meanx = %f, min_ssr_sep = %f\n",meanx,min_ssr_sep));
  if (min_ssr_sep == 0.0) {
    return 1.0;
  }

  sumx_pseudo = NPSEUDO * meanx;

  debug(printf("start = %d, end = %d\n",start,end));

  /* Test for single break */
  for (pos = start+2; pos <= end-2; pos++) {
    if (chrpositions[pos-1] - chrpositions[start] >= MININTERVAL && chrpositions[end-1] - chrpositions[pos] >= MININTERVAL) {
      n_left = pos - start;
      n_right = end - pos;
      meanx_left = (sumx[pos] - sumx[start])/(double) n_left;
      meanx_right = (sumx[end] - sumx[pos])/(double) n_right;
      ssr_left = sumxx[pos] - sumxx[start] - (double) n_left * meanx_left * meanx_left;
      ssr_right = sumxx[end] - sumxx[pos] - (double) n_right * meanx_right * meanx_right;
      ssr_sep = ssr_left + ssr_right;
      
      debug(printf("  one: %d..%d..%d\t%f\t%f\t%f\t%f\t%f\n",
		   start,pos,end,meanx_left,meanx_right,ssr_left,ssr_right,ssr_sep));

      if (ssr_sep >= 0.0 && ssr_sep < min_ssr_sep) {
	min_ssr_sep = ssr_sep;
	*best_start = pos;
	*best_end = pos;
	*best_meanx_left = meanx_left;
	*best_meanx_middle = 0.0;
	*best_meanx_right = meanx_right;
      }
    }
  }

  /* Test for double break */
  for (pos1 = start+2; pos1 <= end-2; pos1++) {
    if (chrpositions[pos1-1] - chrpositions[start] >= MININTERVAL) {
      for (pos2 = pos1+2; pos2 <= end-2; pos2++) {
	if (chrpositions[pos2-1] - chrpositions[pos1] >= MININTERVAL && chrpositions[end-1] - chrpositions[pos2] >= MININTERVAL) {
	  n_outside = (pos1 - start) + (end - pos2);
	  n_inside = pos2 - pos1;
	  meanx_outside = ((sumx[pos1] - sumx[start]) + (sumx[end] - sumx[pos2]) + sumx_pseudo)/
	    ((double) (n_outside + NPSEUDO));
	  meanx_inside = (sumx[pos2] - sumx[pos1] + sumx_pseudo)/((double) (n_inside + NPSEUDO));
	  ssr_outside = (sumxx[pos1] - sumxx[start]) + (sumxx[end] - sumxx[pos2]) - 
	    (double) n_outside * meanx_outside * meanx_outside;
	  ssr_inside = (sumxx[pos2] - sumxx[pos1]) - (double) n_inside * meanx_inside * meanx_inside;
	  ssr_sep = ssr_outside + ssr_inside;
	  
	  debug(printf("  two: %d..%d..%d..%d\t%f\t%f\t%f\t%f\t%f\n",
		       start,pos1,pos2,end,meanx_outside,meanx_inside,
		       ssr_outside,ssr_inside,ssr_sep));

	  /* fscore = (n-2)*(rss - rss_sep)/rss_sep = (n-2)*(rss/rss_sep -
	     1) is maximized when rss_sep is minimized */

	  if (ssr_sep >= 0.0 && ssr_sep < min_ssr_sep) {
	    min_ssr_sep = ssr_sep;
	    *best_start = pos1;
	    *best_end = pos2;
	    *best_meanx_left = (sumx[pos1] - sumx[start])/((double) (pos1 - start));
	    *best_meanx_middle = meanx_inside;
	    *best_meanx_right = (sumx[end] - sumx[pos2])/((double) (end - pos2));
	  }
	}
      }
    }
  }

  if (min_ssr_sep == 0.0) {
    pvalue = 0.0;
  } else {
    fscore = ((double) (n - 2))*(ssr - min_ssr_sep)/min_ssr_sep;
    /* The F test has k and n-2k degrees of freedom.  Here, k = 1, so
       we have 1 and n-2.  For 1 df, t stat is sqrt of F stat. */

    tstat = sqrt(fscore);
    pvalue = NR_ptc(tstat,(double) (n - 2));
    debug(printf("fscore %f at 1,%d df = tstat %f at %d df => pvalue %g\n",
		 fscore,n-2,tstat,n-2,pvalue));
  }

  return pvalue;
}


static T
chrsegment_compute_aux (int start, int end, double *sumx, double *sumxx, unsigned int *chrpositions) {
  int new_start, new_end;
  double pvalue, mean_left, mean_middle, mean_right;
  T chrsegment = NULL, left, middle, right;

  if (chrpositions[end-1] - chrpositions[start] >= MINEVAL) {
    pvalue = interval_detect(&new_start,&new_end,&mean_left,&mean_middle,&mean_right,
			     sumx,sumxx,chrpositions,start,end);
    if (pvalue < PVALUE_THRESHOLD && 
	(fabs(mean_middle-mean_left) > 0.5 || fabs(mean_middle-mean_right) > 0.5)) {
      chrsegment = Chrsegment_new(new_start,new_end,start,end,pvalue,
				  mean_left,mean_middle,mean_right,chrpositions);
      chrsegment->left = chrsegment_compute_aux(start,new_start,sumx,sumxx,chrpositions);
      chrsegment->middle = chrsegment_compute_aux(new_start,new_end,sumx,sumxx,chrpositions);
      chrsegment->right = chrsegment_compute_aux(new_end,end,sumx,sumxx,chrpositions);
    }
  }

  return chrsegment;
}


static double
breakpoint_pvalue (int start, int pos, int end, double *sumx, double *sumxx) {
  double pvalue, ssr, ssr_sep, ssr_left, ssr_right, meanx, meanx_left, meanx_right,
    fscore, tstat;
  int n, n_left, n_right;

  if ((n = end - start) < MINCOUNT) {
    pvalue = 1.0;
  } else {
    meanx = (sumx[end] - sumx[start])/(double) n;
    ssr = (sumxx[end] - sumxx[start]) - (double) n*meanx*meanx;
    
    n_left = pos - start;
    n_right = end - pos;
    meanx_left = (sumx[pos] - sumx[start])/(double) n_left;
    meanx_right = (sumx[end] - sumx[pos])/(double) n_right;
    ssr_left = (sumxx[pos] - sumxx[start]) - (double) n_left*meanx_left*meanx_left;
    ssr_right = (sumxx[end] - sumxx[pos]) - (double) n_right*meanx_right*meanx_right;
    ssr_sep = ssr_left + ssr_right;
      
    fscore = ((double) (n - 2))*(ssr - ssr_sep)/ssr_sep;
    tstat = sqrt(fscore);
    pvalue = NR_ptc(tstat,(double) (n-2));
  }

  return pvalue;
}


void
Chrsegment_compute (Uintlist_T *startpositions, Uintlist_T *endpositions, Doublelist_T *means, 
		    double *x, unsigned int *chrpositions, int nvalues) {
  double *sumx, *sumxx, mean, pvalue, max_pvalue;
  int i;
  int best_start, best_end;
  T chrsegment;
  bool changedp;

  Intlist_T breakpoints = NULL, p, q, r, worst, prev;
  int start, pos, end;

  sumx = (double *) CALLOC(nvalues+1,sizeof(double));
  sumxx = (double *) CALLOC(nvalues+1,sizeof(double));

  sumx[0] = 0.0;
  sumxx[0] = 0.0;
  for (i = 1; i <= nvalues; i++) {
    sumx[i] = sumx[i-1] + x[i-1];
    sumxx[i] = sumxx[i-1] + x[i-1]*x[i-1];
  }

  chrsegment = chrsegment_compute_aux(0,nvalues,sumx,sumxx,chrpositions);
  if (chrsegment == NULL) {
    mean = (sumx[nvalues] - sumx[0])/(double) nvalues;
    chrsegment = Chrsegment_new(/*start*/0,/*end*/nvalues,/*parent_start*/0,/*parent_end*/nvalues,
				/*pvalue*/1.0,/*mean_left*/0.0,mean,/*mean_right*/0.0,chrpositions);
  }

  breakpoints = Intlist_push(NULL,0);
  Chrsegment_breakpoints(&breakpoints,chrsegment);
  breakpoints = Intlist_reverse(breakpoints);

  /* Re-check breakpoints */
  changedp = true;
  while (changedp == true) {
    prev = worst = (Intlist_T) NULL;
    max_pvalue = 0.0;
    p = breakpoints;
    q = Intlist_next(p);
    r = Intlist_next(q);
    while (r != NULL) {
      start = Intlist_head(p);
      pos = Intlist_head(q);
      end = Intlist_head(r);
      pvalue = breakpoint_pvalue(start,pos,end,sumx,sumxx);
      if (pvalue > max_pvalue) {
	max_pvalue = pvalue;
	prev = p;
	worst = q;
      }
      fprintf(stderr,"%d %g\n",pos,pvalue);

      p = Intlist_next(p);
      q = Intlist_next(q);
      r = Intlist_next(r);
    }
    if (worst == NULL || max_pvalue < 1e-10) {
      changedp = false;
    } else {
      fprintf(stderr,"Deleting breakpoint %d\n\n",Intlist_head(worst));
      Intlist_delete(prev,worst);
      changedp = true;
    }
  }

  /* Set up output */
  *startpositions = (Uintlist_T) NULL;
  *endpositions = (Uintlist_T) NULL;
  *means = (Doublelist_T) NULL;

  p = breakpoints;
  q = Intlist_next(p);
  while (q != NULL) {
    start = Intlist_head(p);
    end = Intlist_head(q);

    *startpositions = Uintlist_push(*startpositions,chrpositions[start]);
    *endpositions = Uintlist_push(*endpositions,chrpositions[end-1]);
    *means = Doublelist_push(*means,(sumx[end] - sumx[start])/(double) (end - start));

    p = Intlist_next(p);
    q = Intlist_next(q);
  }

  Chrsegment_free(&chrsegment);
  FREE(sumxx);
  FREE(sumx);
  return;
}



#ifdef TEST
int
main (int argc, char *argv[]) {
  double pvalue, value, *x;
  Doublelist_T valuelist = NULL;
  int n;

  while (scanf("%lf",&value) == 1) {
    valuelist = Doublelist_push(valuelist,value);
  }
  valuelist = Doublelist_reverse(valuelist);
  x = Doublelist_to_array(&n,valuelist);

  Chrsegment_compute(x,n);

  return 0;
}
#endif


