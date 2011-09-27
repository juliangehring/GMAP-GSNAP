static char rcsid[] = "$Id: oligoindex.c,v 1.106 2007/09/17 20:35:53 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "oligoindex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy and memset */
#include "bool.h"
#include "mem.h"
#include "types.h"
#include "genomicpos.h"
#include "list.h"
#include "diag.h"
#include "intlistdef.h"
#include "orderstat.h"


#define DOMINANCE_END_EQUIV 20
#define EXTRA_BOUNDS 120	/* Number of extra diaglines to make active */

#define MAX_DIAGONALS 20

#define THETADIFF1 20.0
#define THETADIFF2 20.0
#define REPOLIGOCOUNT 8

#define OVERABUNDANCE_CHECK 50
#define OVERABUNDANCE_PCT 0.97
#define OVERABUNDANCE_MIN 200
#define CHANGEPOINT 1		/* Needed for bad sequences like BM926731 */

#ifdef PMAP
#define NAMINOACIDS_STAGE2 20
#define NT_PER_MATCH 3
#else
#define NT_PER_MATCH 1
#endif

/* To ensure that we search for initial and final exons */
#if 0
#ifdef PMAP
#define ACTIVE_BUFFER 30
#else
#define ACTIVE_BUFFER 90
#endif
#endif

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* For trimming ends of query sequence */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Segments.  May want to turn on DEBUG2 in stage2.c. */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Diagonals */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


#ifdef PMAP
/* A short oligomer, representing 3 amino acids, or 9 nt, requiring 18
   bits or 2.25 bytes */
#else
/* A short oligomer, typically 8 nt (but could by 9 or 10 nt),
   requiring 16 bits or 2 bytes */
#endif

typedef UINT4 Shortoligomer_T;

#define T Oligoindex_T
struct T {

#ifdef PMAP
  int indexsize_aa;
  Shortoligomer_T msb;
#else
  int indexsize;
  Shortoligomer_T mask;
#endif

  bool query_evaluated_p;
  int oligospace;
  bool *overabundant;
  bool *inquery;
  int *counts;
  int *relevant_counts;
  Genomicpos_T **positions;
  Genomicpos_T **pointers;
};


#ifdef PMAP

#if NAMINOACIDS_STAGE2 == 20
static int aa_index_table[128] =
  { -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,

  /*     *                                  */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,

  /* A,  B,  C,  D,  E,  F,  G,  H,  I,  J, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* K,  L,  M,  N,  O,  P,  Q,  R,  S,  T  */
     8,  9, 10, 11, -1, 12, 13, 14, 15, 16,

  /* U,  V,  W,  X,  Y,  Z  */
    -1, 17, 18, -1, 19, -1,

    -1, -1, -1, -1, -1, -1,

  /* a,  b,  c,  d,  e,  f,  g,  h,  i,  j, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* k,  l,  m,  n,  o,  p,  q,  r,  s,  t  */
     8,  9, 10, 11, -1, 12, 13, 14, 15, 16,

  /* u,  v,  w,  x,  y,  z  */
    -1, 17, 18, -1, 19, -1,

    -1, -1, -1, -1, -1};

#elif NAMINOACIDS_STAGE2 == 18

static int aa_index_table[128] =
  { -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,

  /*     *                                  */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,

  /* A,  B,  C,  D,  E,  F,  G,  H,  I,  J, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* K,  L,  M,  N,  O,  P,  Q,  R,  S,  T  */
     8,  7,  9, 10, -1, 11, 12, 13, 14, 15,

  /* U,  V,  W,  X,  Y,  Z  */
    -1,  7, 16, -1, 17, -1,

    -1, -1, -1, -1, -1, -1,

  /* a,  b,  c,  d,  e,  f,  g,  h,  i,  j, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* k,  l,  m,  n,  o,  p,  q,  r,  s,  t  */
     8,  7,  9, 10, -1, 11, 12, 13, 14, 15,

  /* u,  v,  w,  x,  y,  z  */
    -1,  7, 16, -1, 17, -1,

    -1, -1, -1, -1, -1};

#elif NAMINOACIDS_STAGE2 == 16

static int aa_index_table[128] =
  { -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,

  /*     *                                  */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,

  /* A,  B,  C,  D,  E,  F,  G,  H,  I,  J, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* K,  L,  M,  N,  O,  P,  Q,  R,  S,  T  */
     8,  7,  7,  9, -1, 10, 11,  8, 12, 13,

  /* U,  V,  W,  X,  Y,  Z  */
    -1,  7, 14, -1, 15, -1,

    -1, -1, -1, -1, -1, -1,

  /* a,  b,  c,  d,  e,  f,  g,  h,  i,  j, */
     0, -1,  1,  2,  3,  4,  5,  6,  7, -1,

  /* k,  l,  m,  n,  o,  p,  q,  r,  s,  t  */
     8,  7,  7,  9, -1, 10, 11,  8, 12, 13,

  /* u,  v,  w,  x,  y,  z  */
    -1,  7, 14, -1, 15, -1,

    -1, -1, -1, -1, -1};

#endif	/* NAMINOACIDS_STAGE2 */

#endif	/* PMAP */


static int
power (int base, int exponent) {
  int result = 1, i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


T
Oligoindex_new (int indexsize) {
  T new = (T) MALLOC(sizeof(*new));

#ifdef PMAP
  new->indexsize_aa = indexsize;
  new->msb = power(NAMINOACIDS_STAGE2,indexsize-1);
  debug(printf("msb for indexsize %d is %u\n",indexsize,new->msb));
  new->oligospace = power(NAMINOACIDS_STAGE2,indexsize);
#else
  new->indexsize = indexsize;
  new->mask = ~(~0U << 2*indexsize);
  new->oligospace = power(4,indexsize);
#endif

  new->query_evaluated_p = false;
  new->overabundant = (bool *) CALLOC(new->oligospace,sizeof(bool));
  new->inquery = (bool *) CALLOC(new->oligospace,sizeof(bool));
  new->counts = (int *) CALLOC(new->oligospace,sizeof(int));
  new->relevant_counts = (int *) CALLOC(new->oligospace,sizeof(int));
  new->positions = (Genomicpos_T **) CALLOC(new->oligospace,sizeof(Genomicpos_T *));
  new->pointers = (Genomicpos_T **) CALLOC(new->oligospace,sizeof(Genomicpos_T *));

  return new;
}


/*                87654321 */
#define RIGHT_A 0x00000000
#define RIGHT_C 0x00000001
#define RIGHT_G 0x00000002
#define RIGHT_T 0x00000003

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

static char *
shortoligo_nt (Shortoligomer_T oligo, int oligosize) {
  char *nt;
  int i, j;
  Shortoligomer_T lowbits;

  nt = (char *) CALLOC(oligosize+1,sizeof(char));
  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  return nt;
}


#ifdef PMAP
#if NAMINOACIDS_STAGE2 == 20
static char aa_table[NAMINOACIDS_STAGE2] = "ACDEFGHIKLMNPQRSTVWY";
#elif NAMINOACIDS_STAGE2 == 18
static char aa_table[NAMINOACIDS_STAGE2] = "ACDEFGHIKMNPQRSTWY";
#elif NAMINOACIDS_STAGE2 == 16
static char aa_table[NAMINOACIDS_STAGE2] = "ACDEFGHIKNPQSTWY";
#endif

static char *
aaindex_aa (unsigned int aaindex, int indexsize_aa) {
  char *aa;
  int i, j;

  aa = (char *) CALLOC(indexsize_aa+1,sizeof(char));
  j = indexsize_aa-1;
  for (i = 0; i < indexsize_aa; i++) {
    aa[j] = aa_table[aaindex % NAMINOACIDS_STAGE2];
    aaindex /= NAMINOACIDS_STAGE2;
    j--;
  }

  return aa;
}

#endif


#define NPSEUDO 0.0

/* -1 means edge is on 5', +1 means edge is on 3'  0 means no edge found. */
static int
edge_detect (int *edge, int *sumx, int *sumxx, int length) {
  int side = 0;
  int pos, sumx_left, sumxx_left, sumx_right, sumxx_right, n_left, n_right;
  double theta, sumx_pseudo, x_pseudo, theta_left, theta_right, rss, rss_left, rss_right, rss_sep;
  double fscore, min_rss_sep;

  debug1(printf("\n*** Start of edge_detect\n"));

  sumx_right = sumx[length] - sumx[0];
  sumxx_right = sumxx[length] - sumxx[0];

  theta = (double) sumx_right/(double) length;
  sumx_pseudo = NPSEUDO * theta;
  min_rss_sep = rss = sumxx_right - sumx_right*theta;
  debug1(printf("theta: %d/%d = %f\n",sumx_right,length,theta));
  debug1(printf("rss: %f\n",rss)); 

  debug1(printf("%s %s %s %s %s %s %s %s %s %s %s\n",
		"pos","x","sumx.left","n.left","sumx.right","n.right",
		"theta.left","theta.right","rss.left","rss.right","fscore"));

  n_left = 1;
  n_right = length-1;
  for (pos = 1; pos < length; pos++) {
    sumx_left = sumx[pos] - sumx[0];
    sumxx_left = sumxx[pos] - sumxx[0];
    sumx_right = sumx[length] - sumx[pos];
    sumxx_right = sumxx[length] - sumxx[pos];

    theta_left = ((double) sumx_left + sumx_pseudo)/((double) n_left + NPSEUDO);
    theta_right = ((double) sumx_right + sumx_pseudo)/((double) n_right + NPSEUDO);
    rss_left = sumxx_left - sumx_left*theta_left;
    rss_right = sumxx_right - sumx_right*theta_right;
    rss_sep = rss_left + rss_right;

    if (rss_sep == 0) {
      debug1(printf("%d %d %d %d %d %d %f %f %f %f NA\n",
		    pos,sumx[pos]-sumx[pos-1],sumx_left,n_left,sumx_right,n_right,
		    theta_left,theta_right,rss_left,rss_right));
    } else {
      debug1(      
	     fscore = ((double) (length - 2))*(rss - rss_sep)/rss_sep;
	     debug1(printf("%d %d %d %d %d %d %f %f %f %f %f\n",
			   pos,sumx[pos]-sumx[pos-1],sumx_left,n_left,sumx_right,n_right,
			   theta_left,theta_right,rss_left,rss_right,fscore));
	     );
    }
    /* fscore = (n-2)*(rss - rss_sep)/rss_sep = (n-2)*(rss/rss_sep -
       1) is maximized when rss_sep is minimized */

    if (theta_left > theta_right + THETADIFF1) {
      if (rss_sep < min_rss_sep) {
	min_rss_sep = rss_sep;
	*edge = pos;
	side = -1;
	debug1(printf("Set edge to %d\n",pos));
      }
    } else if (theta_right > theta_left + THETADIFF1) {
      if (rss_sep < min_rss_sep) {
	min_rss_sep = rss_sep;
	*edge = pos;
	side = +1;
	debug1(printf("Set edge to %d\n",pos));
      }
    }

    n_left += 1;
    n_right -= 1;
  }

  debug1(printf("*** End of edge_detect.  Returning %d\n\n",side));

  return side;
}

static int
trim_start_detect (int start, int end, int *sumx, int *sumxx) {
  int edge = -1;
  int pos, sumx_left, sumxx_left, sumx_right, sumxx_right, n_left, n_right;
  double theta, sumx_pseudo, x_pseudo, theta_left, theta_right, rss, rss_left, rss_right, rss_sep;
  double fscore, min_rss_sep;

  debug1(printf("\n*** Start of trim_start_detect\n"));

  sumx_right = sumx[end] - sumx[start];
  sumxx_right = sumxx[end] - sumxx[start];

  if (end <= start) {
    return -1;
  }
  theta = (double) sumx_right/(double) (end - start);
  sumx_pseudo = NPSEUDO * theta;
  min_rss_sep = rss = sumxx_right - sumx_right*theta;
  debug1(printf("%d/%d = %f\n",sumx_right,end-start,theta));
  
  debug1(printf("%s %s %s %s %s %s %s %s %s %s %s\n",
		"pos","counts","sumx.left","n.left","sumx.right","n.right",
		"theta.left","theta.right","rss.left","rss.right","fscore"));

  n_left = 1;
  n_right = end - (start+1);
  for (pos = start+1; pos < end; pos++) {
    sumx_left = sumx[pos] - sumx[start];
    sumxx_left = sumxx[pos] - sumxx[start];
    sumx_right = sumx[end] - sumx[pos];
    sumxx_right = sumxx[end] - sumxx[pos];

    theta_left = ((double) sumx_left + sumx_pseudo)/((double) n_left + NPSEUDO);
    theta_right = ((double) sumx_right + sumx_pseudo)/((double) n_right + NPSEUDO);
    rss_left = sumxx_left - sumx_left*theta_left;
    rss_right = sumxx_right - sumx_right*theta_right;
    rss_sep = rss_left + rss_right;

    if (rss_sep == 0) {
      debug1(printf("%d %d %d %d %d %f %f %f %f NA\n",
		    pos,sumx_left,n_left,sumx_right,n_right,
		    theta_left,theta_right,rss_left,rss_right));
    } else {
      debug1(      
	     fscore = ((double) (end - start - 2))*(rss - rss_sep)/rss_sep;
	     debug1(printf("%d %d %d %d %d %f %f %f %f %f\n",
			   pos,sumx_left,n_left,sumx_right,n_right,
			   theta_left,theta_right,rss_left,rss_right,fscore));
	     );
    }
    /* fscore = (n-2)*(rss - rss_sep)/rss_sep = (n-2)*(rss/rss_sep -
       1) is maximized when rss_sep is minimized */

    if (theta_left < theta_right) {
      debug1(printf("trim_detect aborting early with edge=%d\n",edge));
      return edge;
    } else if (theta_left > theta_right + THETADIFF2) {
      if (rss_sep < min_rss_sep) {
	min_rss_sep = rss_sep;
	edge = pos;
	debug1(printf("Set trim_start to %d\n",pos));      
      }
    }

    n_left += 1;
    n_right -= 1;
  }

  debug1(printf("trim_start_detect returning %d\n",edge));
  return edge;
}

static int
trim_end_detect (int start, int end, int *sumx, int *sumxx) {
  int edge = -1;
  int pos, sumx_left, sumxx_left, sumx_right, sumxx_right, n_left, n_right;
  double theta, sumx_pseudo, x_pseudo, theta_left, theta_right, rss, rss_left, rss_right, rss_sep;
  double fscore, min_rss_sep;

  debug1(printf("\n*** Start of trim_end_detect\n"));

  sumx_right = sumx[end] - sumx[start];
  sumxx_right = sumxx[end] - sumxx[start];

  if (end <= start) {
    return -1;
  }
  theta = (double) sumx_right/(double) (end - start);
  sumx_pseudo = NPSEUDO * theta;
  min_rss_sep = rss = sumxx_right - sumx_right*theta;
  debug1(printf("%d/%d = %f\n",sumx_right,end-start,theta));
  
  debug1(printf("%s %s %s %s %s %s %s %s %s %s %s\n",
		"pos","counts","sumx.left","n.left","sumx.right","n.right",
		"theta.left","theta.right","rss.left","rss.right","fscore"));

  n_left = end - (start+1);
  n_right = 1;
  for (pos = end-1; pos > start; --pos) {
    sumx_left = sumx[pos] - sumx[start];
    sumxx_left = sumxx[pos] - sumxx[start];
    sumx_right = sumx[end] - sumx[pos];
    sumxx_right = sumxx[end] - sumxx[pos];

    theta_left = ((double) sumx_left + sumx_pseudo)/((double) n_left + NPSEUDO);
    theta_right = ((double) sumx_right + sumx_pseudo)/((double) n_right + NPSEUDO);
    rss_left = sumxx_left - sumx_left*theta_left;
    rss_right = sumxx_right - sumx_right*theta_right;
    rss_sep = rss_left + rss_right;

    if (rss_sep == 0) {
      debug1(printf("%d %d %d %d %d %f %f %f %f NA\n",
		    pos,sumx_left,n_left,sumx_right,n_right,
		    theta_left,theta_right,rss_left,rss_right));
    } else {
      debug1(      
	     fscore = ((double) (end - start - 2))*(rss - rss_sep)/rss_sep;
	     debug1(printf("%d %d %d %d %d %f %f %f %f %f\n",
			   pos,sumx_left,n_left,sumx_right,n_right,
			   theta_left,theta_right,rss_left,rss_right,fscore));
	     );
    }
    /* fscore = (n-2)*(rss - rss_sep)/rss_sep = (n-2)*(rss/rss_sep -
       1) is maximized when rss_sep is minimized */

    if (theta_right < theta_left) {
      debug1(printf("trim_detect aborting early with edge=%d\n",edge));
      return edge;
    } else if (theta_right > theta_left + THETADIFF2) {
      if (rss_sep < min_rss_sep) {
	min_rss_sep = rss_sep;
	edge = pos;
	debug1(printf("Set trim_end to %d\n",pos));      
      }
    }

    n_left -= 1;
    n_right += 1;
  }

  debug1(printf("trim_end_detect returning %d\n",edge));
  return edge;
}


/* Run query sequence through this procedure.  First, we count occurrences
 * of each oligo in queryuc (upper case version of queryseq).  This
 * allows us to scan genomicseg intelligently, because then we know
 * whether to store positions for that oligo. */

double
Oligoindex_set_inquery (int *badoligos, int *repoligos, int *trimoligos, int *trim_start, int *trim_end,
			T this, Sequence_T queryuc, bool trimp) {
  double oligodepth;
  int nunique = 0, ngoodoligos, nrepoligos, x, *sumx, *sumxx, sumx0 = 0, sumxx0 = 0, querylength;
  int i = 0, noligos = 0, edge, side;
  int in_counter = 0, querypos;
  Shortoligomer_T oligo = 0U, masked;
  char *p, *ptr;
  bool prevunique = false;
#ifdef PMAP
  int indexsize = this->indexsize_aa, index;
  unsigned int aaindex = 0U;
  char *aa;
#else
  int indexsize = this->indexsize;
  char *nt;
#endif

  this->query_evaluated_p = true; /* Set this flag so we don't redo this part */

  memset((void *) this->inquery,false,this->oligospace*sizeof(bool));
  memset((void *) this->counts,0,this->oligospace*sizeof(int));

  querylength = Sequence_fulllength(queryuc);
  if (querylength <= indexsize) {
    *badoligos = 0;
    *trim_start = 0;
    *trim_end = querylength;
    return 1.0;
  }
    
  for (i = 0, p = Sequence_fullpointer(queryuc); i < querylength; i++, p++) {
    in_counter++;

#ifdef PMAP
    if ((index = aa_index_table[*p]) < 0) {
      aaindex = 0U;
      in_counter = 0;
    } else {
      aaindex = aaindex % this->msb;
      aaindex = aaindex * NAMINOACIDS_STAGE2 + index;
    }
#else
    switch (*p) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }
#endif

    if (in_counter == indexsize) {
#ifdef PMAP
      masked = aaindex;
      debug(aa = aaindex_aa(masked,indexsize); /* Should be indexsize and not indexsize/3 */
	    printf("At querypos %d, aa %s seen (masked = %u)\n",i,aa,masked);
	    FREE(aa));
#else
      masked = oligo & this->mask;
      noligos++;
      debug(nt = shortoligo_nt(oligo,indexsize);
	    printf("At querypos %d, oligo %s seen\n",i,nt);
	    FREE(nt));
#endif
      this->counts[masked] += 1;
      if (this->inquery[masked] == false) {
	nunique += 1;
	this->inquery[masked] = true;
      }
      in_counter--;
    }
  }

  if (trimp == false) {
    *badoligos = 0;
    *trim_start = 0;
    *trim_end = querylength;
    return 1.0;
  }

#ifdef PMAP
  /* trimp should be set to false for PMAP */
  abort();

#else
  /* Determine where to trim using a changepoint analysis */
  sumx = (int *) CALLOC(querylength - indexsize + 1,sizeof(int));
  sumxx = (int *) CALLOC(querylength - indexsize + 1,sizeof(int));

  in_counter = 0;
  querypos = -indexsize;
  for (i = 0, p = Sequence_fullpointer(queryuc); i < querylength; i++, p++) {
    in_counter++;
    querypos++;

    switch (*p) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }

    if (in_counter == indexsize) {
      x = this->counts[oligo&this->mask];
      in_counter--;
    } else {
      x = 1;
    } 

    if (querypos >= 0) {
      sumx0 += x;
      sumxx0 += x*x;
      sumx[querypos] = sumx0;
      sumxx[querypos] = sumxx0;
    }
  }
  sumx[querylength-indexsize] = sumx0;
  sumxx[querylength-indexsize] = sumxx0;


  *trim_start = 0;
  *trim_end = querylength-1;
  if ((side = edge_detect(&edge,sumx,sumxx,querylength-indexsize)) == -1) {
    *trim_start = edge+1;
    if ((edge = trim_end_detect(*trim_start,querylength-indexsize,sumx,sumxx)) >= 0) {
      *trim_end = edge+1;
    }
  } else if (side == +1) {
    *trim_end = edge+1;
    if ((edge = trim_start_detect(0,*trim_end,sumx,sumxx)) >= 0) {
      *trim_start = edge;
    }
  }

  FREE(sumxx);
  FREE(sumx);

  debug1(printf("trim_start = %d, trim_end = %d\n",*trim_start,*trim_end));
#endif

#if 0
  /* If not using changepoint */

  *trim_start = -1;
  *trim_end = 0;

  in_counter = 0;
  querypos = -indexsize;
  for (i = 0, p = Sequence_fullpointer(queryuc); i < querylength; i++, p++) {
    in_counter++;
    querypos++;

    switch (*p) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }

    debug1(if (querypos < 0) {
	     /* Print nothing */
	   } else if (in_counter == indexsize) {
	     printf("%d\n",this->counts[oligo&mask]);
	   } else {
	     printf("0\n");
	   });

    if (in_counter == indexsize) {
      if (this->counts[oligo&this->mask] == 1) {
	if (prevunique == true) {
	  if (*trim_start < 0) {
	    *trim_start = querypos; /* trim 5' at first unique 8-mer */
	  }
	  *trim_end = querypos;	/* trim 3' at last unique 8-mer */
	}
	prevunique = true;
      } else {
	prevunique = false;
      }
      in_counter--;
    } else {
      prevunique = false;
    }
  }

  if (*trim_start < 0) {
    *trim_start = 0;
  } else {
    *trim_start -= 1;		/* To compensate for prevunique */
  }
  *trim_end += indexsize;
#endif

#ifdef PMAP
  *trimoligos = (*trim_end - indexsize) - (*trim_start) + 1;
  *badoligos = 0;
  return 1000000.0;
#else
  /* Count good oligos in trimmed region */
  ngoodoligos = nrepoligos = 0;
  in_counter = 0;
  ptr = Sequence_fullpointer(queryuc);
  for (querypos = (*trim_start)-indexsize, p = &(ptr[*trim_start]);
       querypos < (*trim_end)-indexsize; querypos++, p++) {
    in_counter++;

    switch (*p) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }

    if (in_counter == indexsize) {
      ngoodoligos++;
      if (this->counts[oligo&this->mask] >= REPOLIGOCOUNT) {
	nrepoligos++;
      }
      in_counter--;
    }
  }

  *trimoligos = (*trim_end - indexsize) - (*trim_start) + 1;
  *badoligos = (*trimoligos) - ngoodoligos;
  *repoligos = nrepoligos;

  if (nunique == 0) {
    return 1000000.0;
  } else {
    /* trimlength = *trim_end - *trim_start; */
    oligodepth = (double) noligos/(double) nunique;
    return oligodepth;
  }
#endif
}


/************************************************************************/

#ifdef PMAP

/*              87654321 */
#define CHAR1 0x00000030
#define CHAR2 0x0000000C
#define CHAR3 0x00000003

#define A1    0x00000000	/* 0 */
#define C1    0x00000010	/* 16 */
#define G1    0x00000020	/* 32 */
#define T1    0x00000030	/* 48 */

#define A2    0x00000000	/* 0 */
#define C2    0x00000004	/* 4 */
#define G2    0x00000008	/* 8 */
#define T2    0x0000000C	/* 12 */

#define A3    0x00000000	/* 0 */
#define C3    0x00000001	/* 1 */
#define G3    0x00000002	/* 2 */
#define T3    0x00000003	/* 3 */

/* Do not include stop codon, because oligospace will not allow it */
#if NAMINOACIDS_STAGE2 == 20
typedef enum {AA_A, AA_C, AA_D, AA_E, AA_F, 
	      AA_G, AA_H, AA_I, AA_K, AA_L,
	      AA_M, AA_N, AA_P, AA_Q, AA_R,
	      AA_S, AA_T, AA_V, AA_W, AA_Y} Aminoacid_T;
#elif NAMINOACIDS_STAGE2 == 18
typedef enum {AA_A, AA_C, AA_D, AA_E, AA_F, 
	      AA_G, AA_H, AA_ILV, AA_K,
	      AA_M, AA_N, AA_P, AA_Q, AA_R,
	      AA_S, AA_T, AA_W, AA_Y} Aminoacid_T;
#elif NAMINOACIDS_STAGE2 == 16
typedef enum {AA_A, AA_C, AA_D, AA_E, AA_F, 
	      AA_G, AA_H, AA_ILMV, AA_KR,
	      AA_N, AA_P, AA_Q,
	      AA_S, AA_T, AA_W, AA_Y} Aminoacid_T;
#else
#error The given value of NAMINOACIDS_STAGE2 is not supported by oligoindex.c
#endif

static int
get_last_codon (Shortoligomer_T oligo) {
  switch (oligo & CHAR2) {
  case T2:
    switch (oligo & CHAR1) {
    case T1:
      switch (oligo & CHAR3) {
      case T3: case C3: return AA_F;
      default: /* case A3: case G3: */ 
#if NAMINOACIDS_STAGE2 == 20
	return AA_L;
#elif NAMINOACIDS_STAGE2 == 18
	return AA_ILV;
#elif NAMINOACIDS_STAGE2 == 16
	return AA_ILMV;
#endif
      }	/* CHAR3 */
#if NAMINOACIDS_STAGE2 == 20
    case A1: 
      switch (oligo & CHAR3) {
      case G3: return AA_M;
      default: /* case T3: case A3: case C3: */ return AA_I;
      }
    case C1: return AA_L;
    default: /* case G1: */ return AA_V;
#elif NAMINOACIDS_STAGE2 == 18
    case A1: 
      switch (oligo & CHAR3) {
      case G3: return AA_M;
      default: /* case T3: case A3: case C3: */ return AA_ILV;
      }
    default: /* case C1: case G1: */ return AA_ILV;
#elif NAMINOACIDS_STAGE2 == 16
    default: /* case A1: case C1: case G1: */ return AA_ILMV;
#endif
  } /* CHAR1 */

  case C2:
    switch (oligo & CHAR1) {
    case C1: return AA_P;
    case T1: return AA_S;
    case A1: return AA_T;
    default: /* case G1: */ return AA_A;
    } /* CHAR1 */
  case A2:
    switch (oligo & CHAR1) {
    case T1:
      switch (oligo & CHAR3) {
      case T3: case C3: return AA_Y;
      default: /* case A3: case G3: */ return -1; /* STOP */
      }
    case C1:
      switch (oligo & CHAR3) {
      case T3: case C3: return AA_H;
      default: /* case A3: case G3: */ return AA_Q;
      }
    case A1:
      switch (oligo & CHAR3) {
      case T3: case C3: return AA_N;
      default: /* case A3: case G3: */
#if NAMINOACIDS_STAGE2 == 20
	return AA_K;
#elif NAMINOACIDS_STAGE2 == 18
	return AA_K;
#elif NAMINOACIDS_STAGE2 == 16
	return AA_KR;
#endif
      }
    default: /* case G1: */
      switch (oligo & CHAR3) {
      case T3: case C3: return AA_D;
      default: /* case A3: case G3: */ return AA_E;
      }
    }
  default: /* case G2: */
    switch (oligo & CHAR1) {
    case T1:
      switch (oligo & CHAR3) {
      case T3: case C3: return AA_C;
      case A3: return -1; /* STOP */
      default: /* case G3: */ return AA_W;
      }
    case C1:
#if NAMINOACIDS_STAGE2 == 20
      return AA_R;
#elif NAMINOACIDS_STAGE2 == 18
      return AA_R;
#elif NAMINOACIDS_STAGE2 == 16
      return AA_KR;
#endif
    case A1:
      switch (oligo & CHAR3) {
      case T3: case C3: return AA_S;
      default: /* case A3: case G3: */
#if NAMINOACIDS_STAGE2 == 20
	return AA_R;
#elif NAMINOACIDS_STAGE2 == 18
	return AA_R;
#elif NAMINOACIDS_STAGE2 == 16
	return AA_KR;
#endif
      }
    default: /* case G1: */ return AA_G;
    }
  }

  abort();
  return -1;
}

#endif


/************************************************************************/


/* Second, run genomicuc through this procedure, to determine the genomic positions for each oligo.  */
static int
allocate_positions (Genomicpos_T **pointers, Genomicpos_T **positions, bool *overabundant,
		    bool *inquery, int *counts, int *relevant_counts, int oligospace, int indexsize,
#ifdef PMAP
		    Shortoligomer_T msb,
#else
		    Shortoligomer_T mask,
#endif
		    char *sequence, int seqlength) {
#ifdef PMAP
  int index;
  int frame = 2;
  unsigned int aaindex0 = 0U, aaindex1 = 0U, aaindex2 = 0U;
  int in_counter_0 = 0, in_counter_1 = 0, in_counter_2 = 0;
#endif
  int i = 0, n, sequencepos;
  int in_counter = 0, count;
  Shortoligomer_T oligo = 0U, masked;
  char *p;
  int availslot = 0;
  int nqueryoligos = 0;
  Genomicpos_T *ptr;
  int totalcounts;
  int overabundance_threshold;

#ifdef DEBUG
#ifdef PMAP
  char *aa;
#else
  char *nt;
#endif
#endif

  sequencepos = -indexsize;
  for (i = 0, p = sequence; i < seqlength; i++, p++) {
#ifdef PMAP
    in_counter_0++;
    in_counter_1++;
    in_counter_2++;
#else
    in_counter++;
#endif
    sequencepos++;

    debug(printf("At genomicpos %u, char is %c\n",sequencepos,*p));

    switch (*p) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; 
#ifdef PMAP
      in_counter_0 = in_counter_1 = in_counter_2 = 0; 
#else
      in_counter = 0;
#endif
      break;
    }

#ifdef PMAP
    index = get_last_codon(oligo);
    /* debug(printf("%d %d %d => %d\n",oligo&CHAR1,oligo&CHAR2,oligo&CHAR3,index)); */
    if (frame == 0) {
      frame = 1;
      if (index < 0) {
	aaindex1 = 0;
	in_counter = in_counter_1 = 0;
      } else {
	aaindex1 = aaindex1 % msb;
	masked = aaindex1 = aaindex1 * NAMINOACIDS_STAGE2 + index;
	in_counter = in_counter_1;
      }
    } else if (frame == 1) {
      frame = 2;
      if (index < 0) {
	aaindex2 = 0;
	in_counter = in_counter_2 = 0;
      } else {
	aaindex2 = aaindex2 % msb;
	masked = aaindex2 = aaindex2 * NAMINOACIDS_STAGE2 + index;
	in_counter = in_counter_2;
      }
    } else {
      frame = 0;
      if (index < 0) {
	aaindex0 = 0;
	in_counter = in_counter_0 = 0;
      } else {
	aaindex0 = aaindex0 % msb;
	masked = aaindex0 = aaindex0 * NAMINOACIDS_STAGE2 + index;
	in_counter = in_counter_0;
      }
    }
    debug(printf("frame = %d, masked = %d, in_counters = %d/%d/%d\n",frame,masked,in_counter_0,in_counter_1,in_counter_2));
#endif

    if (in_counter == indexsize) {
#ifndef PMAP
      masked = oligo & mask;
#endif
      if (overabundant[masked] == true) {
	/* Don't bother */
#ifdef PMAP
	debug(aa = aaindex_aa(masked,indexsize/3);
	      printf("At genomicpos %u, aa %s is overabundant\n",sequencepos,aa);
	      FREE(aa));
#else
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s is overabundant\n",sequencepos,nt);
	      FREE(nt));
#endif

      } else if (inquery[masked] == false) {
	/* Don't bother, because it's not in the query sequence */
#ifdef PMAP
	debug(aa = aaindex_aa(masked,indexsize/3);
	      printf("At genomicpos %u, aa %s wasn't seen in querypos\n",
		     sequencepos,aa,counts[masked]);
	      FREE(aa));
#else
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s wasn't seen in querypos\n",sequencepos,nt,counts[masked]);
	      FREE(nt));
#endif

      } else {
	counts[masked] += 1;

#ifdef PMAP
	debug(aa = aaindex_aa(masked,indexsize/3);
	      printf("At genomicpos %u, aa %s seen, counts is now %d\n",
		     sequencepos,aa,counts[masked]);
	      FREE(aa));
#else
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s seen, counts is now %d\n",sequencepos,nt,counts[masked]);
	      FREE(nt));

#endif
      }
#ifndef PMAP
      in_counter--;
#endif
    }
#ifdef PMAP
    if (in_counter_0 == indexsize) {
      in_counter_0--;
    }
    if (in_counter_1 == indexsize) {
      in_counter_1--;
    }
    if (in_counter_2 == indexsize) {
      in_counter_2--;
    }
#endif
  }

  n = 0;
  for (i = 0; i < oligospace; i++) {
    if (counts[i] > 0) {
      relevant_counts[n++] = counts[i];
    }
  }

  totalcounts = 0;
  if (n < OVERABUNDANCE_CHECK) {
    debug(printf("only %d entries => don't use orderstat\n",n));

    for (i = 0; i < oligospace; i++) {
      totalcounts += counts[i];
    }

  } else {
    overabundance_threshold = Orderstat_int_pct_inplace(relevant_counts,n,OVERABUNDANCE_PCT);
    debug(printf("overabundance threshold is %d\n",overabundance_threshold));
    if (overabundance_threshold < OVERABUNDANCE_MIN) {
      overabundance_threshold = OVERABUNDANCE_MIN;
      debug(printf("  => resetting to %d\n",overabundance_threshold));
    }

    for (i = 0; i < oligospace; i++) {
      if (counts[i] > overabundance_threshold) {
	overabundant[i] = true;
	counts[i] = 0;
      } else {
	totalcounts += counts[i];
      }
    }
  }

  if (totalcounts == 0) {
    positions[0] = (Genomicpos_T *) NULL;
  } else {
    ptr = (Genomicpos_T *) CALLOC(totalcounts,sizeof(Genomicpos_T));
    for (i = 0; i < oligospace; i++) {
      positions[i] = ptr;
      ptr += counts[i];
    }
    memcpy((void *) pointers,positions,oligospace*sizeof(Genomicpos_T *));
  }

  /* printf("nqueryoligos = %d\n",nqueryoligos); */

  return totalcounts;
}


/* Third, run genomicuc through this procedure, to determine the genomic positions for each oligo.  */
/* Logic of this procedure should match that of allocate_positions */
static void
store_positions (Genomicpos_T **pointers, bool *overabundant, 
		 bool *inquery, int oligospace, int indexsize,
#ifdef PMAP
		 Shortoligomer_T msb,
#else
		 Shortoligomer_T mask,
#endif
		 char *sequence, int seqlength) {
#ifdef PMAP
  int index;
  int frame = 2;
  unsigned int aaindex0 = 0U, aaindex1 = 0U, aaindex2 = 0U;
  int in_counter_0 = 0, in_counter_1 = 0, in_counter_2 = 0;
#endif
  int i = 0, sequencepos;
  int in_counter = 0, count;
  Shortoligomer_T oligo = 0U, masked;
  char *p;
  int availslot = 0;
  int nqueryoligos = 0;

#ifdef DEBUG
#ifdef PMAP
  char *aa;
#else
  char *nt;
#endif
#endif

  sequencepos = -indexsize;
  for (i = 0, p = sequence; i < seqlength; i++, p++) {
#ifdef PMAP
    in_counter_0++;
    in_counter_1++;
    in_counter_2++;
#else
    in_counter++;
#endif
    sequencepos++;

    debug(printf("At genomicpos %u, char is %c\n",sequencepos,*p));

    switch (*p) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; 
#ifdef PMAP
      in_counter_0 = in_counter_1 = in_counter_2 = 0; 
#else
      in_counter = 0;
#endif
      break;
    }

#ifdef PMAP
    index = get_last_codon(oligo);
    /* debug(printf("%d %d %d => %d\n",oligo&CHAR1,oligo&CHAR2,oligo&CHAR3,index)); */
    if (frame == 0) {
      frame = 1;
      if (index < 0) {
	aaindex1 = 0;
	in_counter = in_counter_1 = 0;
      } else {
	aaindex1 = aaindex1 % msb;
	masked = aaindex1 = aaindex1 * NAMINOACIDS_STAGE2 + index;
	in_counter = in_counter_1;
      }
    } else if (frame == 1) {
      frame = 2;
      if (index < 0) {
	aaindex2 = 0;
	in_counter = in_counter_2 = 0;
      } else {
	aaindex2 = aaindex2 % msb;
	masked = aaindex2 = aaindex2 * NAMINOACIDS_STAGE2 + index;
	in_counter = in_counter_2;
      }
    } else {
      frame = 0;
      if (index < 0) {
	aaindex0 = 0;
	in_counter = in_counter_0 = 0;
      } else {
	aaindex0 = aaindex0 % msb;
	masked = aaindex0 = aaindex0 * NAMINOACIDS_STAGE2 + index;
	in_counter = in_counter_0;
      }
    }
    debug(printf("frame = %d, masked = %d, in_counters = %d/%d/%d\n",frame,masked,in_counter_0,in_counter_1,in_counter_2));
#endif

    if (in_counter == indexsize) {
#ifndef PMAP
      masked = oligo & mask;
#endif
      if (overabundant[masked] == true) {
	/* Don't bother */

      } else if (inquery[masked] == false) {
	/* Don't bother, because it's not in the query sequence */

      } else {
	pointers[masked][0] = (Genomicpos_T) sequencepos;
	pointers[masked]++;
      }
#ifndef PMAP
      in_counter--;
#endif
    }
#ifdef PMAP
    if (in_counter_0 == indexsize) {
      in_counter_0--;
    }
    if (in_counter_1 == indexsize) {
      in_counter_1--;
    }
    if (in_counter_2 == indexsize) {
      in_counter_2--;
    }
#endif
  }

  /* printf("nqueryoligos = %d\n",nqueryoligos); */

  return;
}


static void
dump_positions (Genomicpos_T **positions, int *counts, int oligospace, int indexsize) {
  int i;
  char *aa, *nt;

#ifdef PMAP
  for (i = 0; i < oligospace; i++) {
    aa = aaindex_aa(i,indexsize/3);
    if (counts[i] >= 1) {
      printf("AA %s => %d entries: %u...%u\n",
	     aa,counts[i],positions[i][0],positions[i][counts[i]-1]);
    }
    FREE(aa);
  }
#else
  for (i = 0; i < oligospace; i++) {
    nt = shortoligo_nt(i,indexsize);
    if (counts[i] >= 1) {
      printf("Oligo %s => %d entries: %u...%u\n",
	     nt,counts[i],positions[i][0],positions[i][counts[i]-1]);
    }
    FREE(nt);
  }
#endif

  return;
}


#ifndef PMAP
#define POLY_A 0x0000
#define POLY_C 0x5555
#define POLY_G 0xAAAA
#define POLY_T 0xFFFF
#endif

void
Oligoindex_tally (T this, Sequence_T genomicuc, Sequence_T queryuc) {
  int badoligos, repoligos, trimoligos, trim_start, trim_end;

  if (this->query_evaluated_p == false) {
    Oligoindex_set_inquery(&badoligos,&repoligos,&trimoligos,&trim_start,&trim_end,this,queryuc,/*trimp*/false);
  }

  memset((void *) this->counts,0,this->oligospace*sizeof(int));
  memset((void *) this->overabundant,0,this->oligospace*sizeof(bool));
#ifndef PMAP
  /* These values will prevent oligoindex from getting mappings later */
  this->overabundant[POLY_A & this->mask] = true;
  this->overabundant[POLY_C & this->mask] = true;
  this->overabundant[POLY_G & this->mask] = true;
  this->overabundant[POLY_T & this->mask] = true;
#endif

  if ((allocate_positions(this->pointers,this->positions,this->overabundant,
			  this->inquery,this->counts,this->relevant_counts,this->oligospace,
#ifdef PMAP
			  3*this->indexsize_aa,this->msb,
#else
			  this->indexsize,this->mask,
#endif
			  Sequence_trimpointer(genomicuc),
			  Sequence_trimlength(genomicuc))) > 0) {

    store_positions(this->pointers,this->overabundant,this->inquery,this->oligospace,
#ifdef PMAP
		    3*this->indexsize_aa,this->msb,
#else
		    this->indexsize,this->mask,
#endif
		    Sequence_trimpointer(genomicuc),
		    Sequence_trimlength(genomicuc));

#ifdef PMAP
    debug(dump_positions(this->positions,this->counts,this->oligospace,3*this->indexsize_aa));
#else
    debug(dump_positions(this->positions,this->counts,this->oligospace,this->indexsize));
#endif
  }

  return;
}
  
void
Oligoindex_untally (T this) {
  if (this->positions[0] != NULL) {
    FREE(this->positions[0]);
  }
  return;
}

void
Oligoindex_clear_inquery (T this) {
  this->query_evaluated_p = false;
  return;
}



void
Oligoindex_free (T *old) {
  if (*old) {
    FREE((*old)->pointers);
    FREE((*old)->positions);
    FREE((*old)->relevant_counts);
    FREE((*old)->counts);
    FREE((*old)->inquery);
    FREE((*old)->overabundant);
    FREE(*old);
  }
  return;
}

static Genomicpos_T *
lookup (int *nhits, T this, Shortoligomer_T masked) {
#ifdef PMAP
  char *aa;
#else
  char *nt;
#endif

  if ((*nhits = this->counts[masked]) >= 1) {
#ifdef PMAP
    debug(aa = aaindex_aa(masked,this->indexsize_aa);
	  printf("masked %s => %d entries: %u...%u\n",
		 aa,*nhits,this->positions[masked][0],this->positions[masked][*nhits-1]);
	  FREE(aa));
#else
    debug(nt = shortoligo_nt(masked,this->indexsize);
	  printf("masked %s => %d entries: %u...%u\n",
		 nt,*nhits,this->positions[masked][0],this->positions[masked][*nhits-1]);
	  FREE(nt));
#endif
    return this->positions[masked];
  } else {
#ifdef PMAP
    debug(aa = aaindex_aa(masked,this->indexsize_aa);
	  printf("masked %s not found\n",aa);
	  FREE(aa));
#else
    debug(nt = shortoligo_nt(masked,this->indexsize);
	  printf("masked %s not found\n",nt);
	  FREE(nt));
#endif
    /* Warning: *nhits might be -1 here, but shouldn't affect anything */
    return NULL;
  }
}


static bool
consecutivep (int prev_querypos, unsigned int *prev_mappings, int prev_nhits,
	      int cur_querypos, unsigned int *cur_mappings, int cur_nhits) {
  int genomicdist, i, j;

  if (prev_nhits > 0 && cur_nhits > 0) {
    j = i = 0;
    genomicdist = NT_PER_MATCH*(cur_querypos - prev_querypos);
    while (j < prev_nhits && i < cur_nhits) {
      /* printf("Comparing %u with %u + %d\n",cur_mappings[i],prev_mappings[j],NT_PER_MATCH); */
      if (cur_mappings[i] == prev_mappings[j] + genomicdist) {
	/* printf("true\n"); */
	return true;
      } else if (cur_mappings[i] < prev_mappings[j] + genomicdist) {
	i++;
      } else {
	j++;
      }
    }
  }
  /* printf("false\n"); */
  return false;
}


static int
abs_compare (const void *x, const void *y) {
  int absa, absb;
  int a = * (int *) x;
  int b = * (int *) y;

  absa = (a >= 0) ? a : -a;
  absb = (b >= 0) ? b : -b;

  if (absa < absb) {
    return -1;
  } else if (absa > absb) {
    return +1;
  } else if (a >= 0 && b < 0) {
    /* positives before negatives */
    return -1;
  } else if (a < 0 && b >= 0) {
    return +1;
  } else {
    return 0;
  }
}


/* Assume array sorted by querystart */
static Diag_T *
compute_dominance (int *nunique, Diag_T *array, int ndiagonals, int indexsize,
		   bool debug_graphic_p) {
  int superstart, superend, expected_nconsecutive, threshold;
  int i, j, k;
  Diag_T super, sub;

  if (debug_graphic_p == true) {
    for (i = 0; i < ndiagonals; i++) {
      sub = array[i];
#ifdef PMAP
      printf("segments(%d,%d+%d,%d,%d+%d,col=\"black\")  # nconsecutive = %d\n",
	     Diag_querystart(sub),Diag_diagonal(sub),3*Diag_querystart(sub),
	     Diag_queryend(sub),Diag_diagonal(sub),3*Diag_queryend(sub),Diag_nconsecutive(sub));
#else
      printf("segments(%d,%d+%d,%d,%d+%d,col=\"black\")  # nconsecutive = %d\n",
	     Diag_querystart(sub),Diag_diagonal(sub),Diag_querystart(sub),
	     Diag_queryend(sub),Diag_diagonal(sub),Diag_queryend(sub),Diag_nconsecutive(sub));
#endif
    }
  }

  qsort(array,ndiagonals,sizeof(Diag_T),Diag_compare_nconsecutive);

  *nunique = ndiagonals;
  i = 0;
  while (i < *nunique) {
    super = array[i];
    superstart = Diag_querystart(super);
    superend = Diag_queryend(super);

    expected_nconsecutive = superend + 1 - superstart;
    if (Diag_nconsecutive(super) > expected_nconsecutive - 10) {
      threshold = Diag_nconsecutive(super) - DOMINANCE_END_EQUIV;
      for (j = i+1; j < *nunique; j++) {
	sub = array[j];
	if (Diag_querystart(sub) >= superstart && Diag_queryend(sub) <= superend && Diag_nconsecutive(sub) < threshold) {
	  Diag_set_dominatedp(sub);
	}
      }

      /* Shift array to contain non-dominated diagonals */
      k = i+1;
      for (j = i+1; j < *nunique; j++) {
	sub = array[j];
	if (Diag_dominatedp(sub) == false) {
	  array[k++] = array[j];
	}
      }
      *nunique = k;
    }
    i++;
  }

  return array;
}
  

static int
diagonal_coverage (int *clear_coverage, Diag_T *array, int nunique) {
  int coverage = 0, regionstart, clear_regionstart, i, j = 0, status;
  int *events, nevents;
  Diag_T diag;
  
  *clear_coverage = 0;
  nevents = nunique+nunique;

  events = (int *) CALLOC(nevents,sizeof(int));
  for (i = 0; i < nunique; i++) {
    diag = array[i];
    events[j++] = +Diag_querystart(diag);
    events[j++] = -Diag_queryend(diag);
  }
  qsort(events,nevents,sizeof(int),abs_compare);

  status = 0;
  for (j = 0; j < nevents; j++) {
    if (events[j] >= 0) {
      status++;
      if (status == 1) {
	/* Start of region.  Now in state 1..MAX_DIAGONALS */
	clear_regionstart = events[j];
	regionstart = events[j];
      } else if (status == MAX_DIAGONALS + 1) {
	/* End of region.  Now in state MAX_DIAGONALS+1 .. Inf */
	*clear_coverage += (events[j]) - clear_regionstart + 1;
      }
    } else {
      status--;
      if (status == MAX_DIAGONALS) {
	/* Start of region.  Now in state 1..MAX_DIAGONALS */
	clear_regionstart = -events[j];
      } else if (status == 0) {
	/* End of region.  Now in state 0 */
	coverage += (-events[j]) - regionstart + 1;
	*clear_coverage += (-events[j]) - clear_regionstart + 1;
      }
    }
  }

  FREE(events);
  
  return coverage;
}


static int
compute_bounds (double *pct_coverage, double *pct_clear_coverage, 
		unsigned int *minactive, unsigned int *maxactive, List_T diagonals,
		int genomiclength, int querylength, int indexsize, bool debug_graphic_p) {
  int nunique, ndiagonals, diagonal, i, j;
  int querypos, position;
  Diag_T *array, diag;
  int coverage, clear_coverage;
  int activestart, activeend;

  /* Sort diagonals */
  if ((ndiagonals = List_length(diagonals)) == 0) {
    nunique = 0;
    *pct_coverage = 0.0;
    for (querypos = 0; querypos < querylength; querypos++) {
      minactive[querypos] = 0U;
      maxactive[querypos] = genomiclength;
    }
  } else {
    if (debug_graphic_p == true) {
      printf("plot(c(%d,%d),c(%d,%d),type=\"n\",xlab=\"Query\",ylab=\"Genomic\")\n",0,querylength,0,genomiclength);
    }

    array = (Diag_T *) List_to_array(diagonals,NULL);
    array = compute_dominance(&nunique,array,ndiagonals,indexsize,debug_graphic_p);

    coverage = diagonal_coverage(&clear_coverage,array,nunique);
    *pct_coverage = (double) coverage/(double) querylength;
    *pct_clear_coverage = (double) clear_coverage/(double) querylength;
    debug2(fprintf(stderr,"coverage = %d/%d = %f\n",coverage,querylength,*pct_coverage));
    debug2(fprintf(stderr,"clear coverage = %d/%d = %f\n",clear_coverage,querylength,*pct_clear_coverage));

    qsort(array,nunique,sizeof(Diag_T),Diag_compare_diagonal);
    if (debug_graphic_p == true) {
      for (i = 0; i < nunique; i++) {
	diag = array[i];
#ifdef PMAP
	printf("segments(%d,%d+%d,%d,%d+%d,col=\"red\")  # nconsecutive = %d\n",
	       Diag_querystart(diag),Diag_diagonal(diag),3*Diag_querystart(diag),
	       Diag_queryend(diag),Diag_diagonal(diag),3*Diag_queryend(diag),Diag_nconsecutive(diag));
#else
	printf("segments(%d,%d+%d,%d,%d+%d,col=\"red\")  # nconsecutive = %d\n",
	       Diag_querystart(diag),Diag_diagonal(diag),Diag_querystart(diag),
	       Diag_queryend(diag),Diag_diagonal(diag),Diag_queryend(diag),Diag_nconsecutive(diag));
#endif
      }
    }

    /* Set minactive */
#ifdef ACTIVE_BUFFER
    /* Allow buffer on 5' end to make sure we identify the best initial exon */
    if ((activestart = Diag_querystart(array[0])) < ACTIVE_BUFFER) {
      activestart = ACTIVE_BUFFER;
      if (activestart > querylength) {
	activestart = querylength;
      }
    }
#else
    activestart = Diag_querystart(array[0]);
#endif

    for (querypos = 0; querypos < activestart; querypos++) {
      minactive[querypos] = 0U;
    }

    diagonal = Diag_diagonal(array[0]);
    for ( ; querypos < Diag_queryend(array[0]); querypos++) {
#ifdef PMAP
      if ((position = diagonal + 3*querypos - EXTRA_BOUNDS) < 0) {
	minactive[querypos] = 0U;
      } else {
	minactive[querypos] = (unsigned int) position;
      }
#else
      if ((position = diagonal + querypos - EXTRA_BOUNDS) < 0) {
	minactive[querypos] = 0U;
      } else {
	minactive[querypos] = (unsigned int) position;
      }
#endif
    }

    i = 0;
    while (i < nunique) {
      j = i+1;
      while (j < nunique && Diag_queryend(array[j]) <= Diag_queryend(array[i])) {
	j++;
      }
      diagonal = Diag_diagonal(array[i]);
      if (j < nunique) {
	for ( ; querypos <= Diag_queryend(array[j]); querypos++) {
#ifdef PMAP
	  if ((position = diagonal + 3*querypos - EXTRA_BOUNDS) < 0) {
	    minactive[querypos] = 0U;
	  } else {
	    minactive[querypos] = (unsigned int) position;
	  }
#else
	  if ((position = diagonal + querypos - EXTRA_BOUNDS) < 0) {
	    minactive[querypos] = 0U;
	  } else {
	    minactive[querypos] = (unsigned int) position;
	  }
#endif
	}
      }
      i = j;
    }

    for ( ; querypos < querylength; querypos++) {
#ifdef PMAP
      if ((position = diagonal + 3*querypos - EXTRA_BOUNDS) < 0) {
	minactive[querypos] = 0U;
      } else {
	minactive[querypos] = (unsigned int) position;
      }
#else
      if ((position = diagonal + querypos - EXTRA_BOUNDS) < 0) {
	minactive[querypos] = 0U;
      } else {
	minactive[querypos] = (unsigned int) position;
      }
#endif
    }

    /* Set maxactive */
#ifdef ACTIVE_BUFFER
    /* Allow buffer on 3' end to make sure we identify the best final exon */
    if ((activeend = Diag_queryend(array[nunique-1])) > querylength-ACTIVE_BUFFER) {
      activeend = querylength-ACTIVE_BUFFER;
      if (activeend < 0) {
	activeend = 0;
      }
    }
#else
    activeend = Diag_queryend(array[nunique-1]);
#endif

    for (querypos = querylength-1; querypos > activeend; --querypos) {
      maxactive[querypos] = genomiclength;
    }

    diagonal = Diag_diagonal(array[nunique-1]);
    for ( ; querypos > Diag_querystart(array[nunique-1]); --querypos) {
#ifdef PMAP
      if ((position = diagonal + 3*querypos + EXTRA_BOUNDS) > genomiclength) {
	maxactive[querypos] = genomiclength;
      } else {
	maxactive[querypos] = (unsigned int) position;
      }
#else
      if ((position = diagonal + querypos + EXTRA_BOUNDS) > genomiclength) {
	maxactive[querypos] = genomiclength;
      } else {
	maxactive[querypos] = (unsigned int) position;
      }
#endif
    }

    i = nunique-1;
    while (i >= 0) {
      j = i-1;
      while (j >= 0 && Diag_querystart(array[j]) > Diag_querystart(array[i])) {
	--j;
      }
      diagonal = Diag_diagonal(array[i]);
      if (j >= 0) {
	for ( ; querypos >= Diag_querystart(array[j]); --querypos) {
#ifdef PMAP
	  if ((position = diagonal + 3*querypos + EXTRA_BOUNDS) > genomiclength) {
	    maxactive[querypos] = genomiclength;
	  } else {
	    maxactive[querypos] = (unsigned int) position;
	  }
#else
	  if ((position = diagonal + querypos + EXTRA_BOUNDS) > genomiclength) {
	    maxactive[querypos] = genomiclength;
	  } else {
	    maxactive[querypos] = (unsigned int) position;
	  }
#endif
	}
      }
      i = j;
    }

    for ( ; querypos >= 0; --querypos) {
#ifdef PMAP
      if ((position = diagonal + 3*querypos + EXTRA_BOUNDS) > genomiclength) {
	maxactive[querypos] = genomiclength;
      } else {
	maxactive[querypos] = (unsigned int) position;
      }
#else
      if ((position = diagonal + querypos + EXTRA_BOUNDS) > genomiclength) {
	maxactive[querypos] = genomiclength;
      } else {
	maxactive[querypos] = (unsigned int) position;
      }
#endif
    }

    FREE(array);
  }

  return nunique;
}


struct Genomicdiag_T {
  int i;
  int querypos;
  int best_nconsecutive;
  int nconsecutive;
  int best_consecutive_start;
  int consecutive_start;
  int best_consecutive_end;
};
typedef struct Genomicdiag_T *Genomicdiag_T;


/* Third, retrieves appropriate oligo information for a given querypos and
   copies it to that querypos */
/* Note: Be careful on totalpositions, because nhits may be < 0 */
unsigned int **
Oligoindex_get_mappings (int *nunique, int *maxnconsecutive, double *pct_coverage, double *pct_clear_coverage,
			 int **npositions, int *totalpositions,
			 unsigned int *minactive, unsigned int *maxactive,
			 T this, Sequence_T queryuc, Sequence_T genomicuc, Diagpool_T diagpool,
			 int local_lookback, bool debug_graphic_p) {
  unsigned int **mappings, position;
  int nhits, hit, diagi, diagi_adjustment, i;
  int querylength, genomiclength, trimlength, trim_start;
#ifdef PREV_MAXCONSECUTIVE
  int prev_querypos, prev_nhits;
  unsigned int *prev_mappings;
  int ngoodconsecutive;
#endif

  char *p;
  int in_counter = 0, querypos;
  Shortoligomer_T oligo = 0U, masked;
  Intlist_T prevdiag;

  struct Genomicdiag_T *genomicdiag;
  Genomicdiag_T ptr;
  List_T diagonals = NULL, good_genomicdiags = NULL, q;

#ifdef PMAP
  int indexsize = this->indexsize_aa, index;
  unsigned int aaindex = 0U;
#else
  int indexsize = this->indexsize;
#endif

  querylength = Sequence_fulllength(queryuc);
  genomiclength = Sequence_fulllength(genomicuc);
  trimlength = Sequence_trimlength(queryuc);
  trim_start = Sequence_trim_start(queryuc);
  mappings = (unsigned int **) CALLOC(querylength,sizeof(unsigned int *));
  *npositions = (int *) CALLOC(querylength,sizeof(int));
  *totalpositions = 0;

  *maxnconsecutive = 0;
#ifdef PMAP
  genomicdiag = (struct Genomicdiag_T *) CALLOC(3*querylength+genomiclength,sizeof(struct Genomicdiag_T));
  for (diagi = 0; diagi < 3*querylength+genomiclength; diagi++) {
    genomicdiag[diagi].i = diagi;
    genomicdiag[diagi].querypos = -local_lookback; /* guarantees first check won't be consecutive */
  }
#else
  genomicdiag = (struct Genomicdiag_T *) CALLOC(querylength+genomiclength,sizeof(struct Genomicdiag_T));
  for (diagi = 0; diagi < querylength+genomiclength; diagi++) {
    genomicdiag[diagi].i = diagi;
    genomicdiag[diagi].querypos = -local_lookback; /* guarantees first check won't be consecutive */
  }
#endif

  querypos = trim_start - indexsize;
  for (i = 0, p = Sequence_trimpointer(queryuc); i < trimlength; i++, p++) {
    in_counter++;
    querypos++;
    
#ifdef PMAP
    if ((index = aa_index_table[*p]) < 0) {
      aaindex = 0U;
      in_counter = 0;
    } else {
      aaindex = aaindex % this->msb;
      aaindex = aaindex * NAMINOACIDS_STAGE2 + index;
    }
#else
    switch (*p) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }
#endif

    if (in_counter == indexsize) {
#ifdef PMAP
      masked = aaindex;
#else
      masked = oligo & this->mask;
#endif
      mappings[querypos] = lookup(&nhits,this,masked);
      (*npositions)[querypos] = nhits;
      if (nhits > 0) {
	*totalpositions += nhits;

#ifdef PMAP
	/* diagonal is (position - 3*querypos); diagi is (position - 3*querypos) + 3*querylength */
	diagi_adjustment = 3*(querylength - querypos);
#else
	/* diagonal is (position - querypos); diagi is (position - querypos) + querylength */
	diagi_adjustment = querylength - querypos;
#endif
	for (hit = 0; hit < nhits; hit++) {
	  diagi = mappings[querypos][hit] + diagi_adjustment;
	  ptr = &(genomicdiag[diagi]);
	  /* Must use >= here, so querypos 0 - (-local_lookback) will fail */
	  if (querypos - ptr->querypos >= local_lookback) {
	    debug3(printf("At diagi %d (checking querypos %d to %d), no consecutive\n",diagi,ptr->querypos,querypos));
	    ptr->nconsecutive = 0;
	    ptr->consecutive_start = querypos;

	  } else if (++ptr->nconsecutive > ptr->best_nconsecutive) {
	    ptr->best_consecutive_start = ptr->consecutive_start;
	    ptr->best_consecutive_end = querypos;
	    ptr->best_nconsecutive = ptr->nconsecutive;
	    debug3(printf("At diagi %d (checking querypos %d to %d), best consecutive of %d from %d to %d\n",
			  diagi,ptr->querypos,querypos,ptr->best_nconsecutive,
			  ptr->best_consecutive_start,ptr->best_consecutive_end));
	    if (ptr->best_nconsecutive == SUFFNCONSECUTIVE) {
	      good_genomicdiags = List_push(good_genomicdiags,(void *) ptr);
	    }
	    if (ptr->best_nconsecutive > *maxnconsecutive) {
	      *maxnconsecutive = ptr->best_nconsecutive;
	    }
	  }
	  ptr->querypos = querypos;
	}
      }
      in_counter--;
    }
  }

  while (good_genomicdiags != NULL) {
    good_genomicdiags = List_pop(good_genomicdiags,(void **) &ptr);
#ifdef PMAP
    diagonals = Diagpool_push(diagonals,diagpool,/*diagonal*/(ptr->i - 3*querylength),
			      ptr->best_consecutive_start,ptr->best_consecutive_end,
			      ptr->best_nconsecutive+1);
#else
    diagonals = Diagpool_push(diagonals,diagpool,/*diagonal*/(ptr->i - querylength),
			      ptr->best_consecutive_start,ptr->best_consecutive_end,
			      ptr->best_nconsecutive+1);
#endif
  }

  FREE(genomicdiag);

  *nunique = compute_bounds(&(*pct_coverage),&(*pct_clear_coverage),minactive,maxactive,diagonals,
			    genomiclength,querylength,indexsize,debug_graphic_p);

  debug2(fprintf(stderr,"%d diagonals (%d not dominated), maxnconsecutive = %d\n",
		 List_length(diagonals),*nunique,*maxnconsecutive));

  Diagpool_reset(diagpool);
  /* No need to free diagonals */

  return mappings;
}


