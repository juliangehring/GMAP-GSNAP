static char rcsid[] = "$Id: oligoindex.c 79521 2012-11-19 22:11:24Z twu $";
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
#include "assert.h"
#include "mem.h"
#include "list.h"
#include "diag.h"
#include "intlistdef.h"
#include "orderstat.h"


#ifndef PMAP
#define STRAIGHT_MASK_2  0x0000000F /* 2-mer: 1111 */
#define STRAIGHT_MASK_3  0x0000003F /* 3-mer: 11 1111 */
#define STRAIGHT_MASK_4  0x000000FF /* 4-mer: 1111 1111 */
#define STRAIGHT_MASK_5  0x000003FF /* 5-mer: 11 1111 1111 */
#define STRAIGHT_MASK_6  0x00000FFF /* 6-mer: 1111 1111 1111 */
#define STRAIGHT_MASK_7  0x00003FFF /* 7-mer: 11 1111 1111 1111 */
#define STRAIGHT_MASK_8  0x0000FFFF /* 8-mer: 1111 1111 1111 1111 */
#define STRAIGHT_MASK_9  0x0003FFFF /* 9-mer: 11 1111 1111 1111 1111 */
#define STRAIGHT_MASK_10 0x000FFFFF /* 10-mer: 1111 1111 1111 1111 1111 */
#define STRAIGHT_MASK_11 0x003FFFFF /* 11-mer: 11 1111 1111 1111 1111 1111 */

#define WOBBLE_MASK_2  0x0000000F /* 2-mer: 1111 */
#define WOBBLE_MASK_3  0x0000003C /* 3-mer: 11 1100 */
#define WOBBLE_MASK_4  0x000000F3 /* 4-mer: 1111 0011 */
#define WOBBLE_MASK_5  0x000003CF /* 5-mer: 11 1100 1111 */
#define WOBBLE_MASK_6  0x00000F3C /* 6-mer: 1111 0011 1100 */
#define WOBBLE_MASK_7  0x00003CF3 /* 7-mer: 11 1100 1111 0011 */
#define WOBBLE_MASK_8  0x0000F3CF /* 8-mer: 1111 0011 1100 1111 */
#define WOBBLE_MASK_9  0x0003CF3C /* 9-mer: 11 1100 1111 0011 1100 */
#define WOBBLE_MASK_10 0x000F3CF3 /* 10-mer: 1111 0011 1100 1111 0011 */
#define WOBBLE_MASK_11 0x003CF3CF /* 11-mer: 11 1100 1111 0011 1100 1111 */
#endif


#ifdef PMAP
#define NOLIGOINDICES_MAJOR 1
static int indexsizes_major[NOLIGOINDICES_MAJOR] = {6}; /* In aa */
static int diag_lookbacks_major[NOLIGOINDICES_MAJOR] = {40};
static int suffnconsecutives_major[NOLIGOINDICES_MAJOR] = {10};

#define NOLIGOINDICES_MINOR 1
static int indexsizes_minor[NOLIGOINDICES_MINOR] = {6}; /* In aa */
static int diag_lookbacks_minor[NOLIGOINDICES_MINOR] = {40};
static int suffnconsecutives_minor[NOLIGOINDICES_MINOR] = {10};

#elif defined(GSNAP)

/* Have fewer to enable speedup.  Note: Including 7-mers causes an 8x
   increase in run-time for score_querypos, and including 6-mers causes a
   30x increase. */
#define NOLIGOINDICES_MAJOR 1
static int indexsizes_major[NOLIGOINDICES_MAJOR] = {8};
static Shortoligomer_T masks_major[NOLIGOINDICES_MAJOR] = {STRAIGHT_MASK_8};
static int diag_lookbacks_major[NOLIGOINDICES_MAJOR] = {120};
static int suffnconsecutives_major[NOLIGOINDICES_MAJOR] = {20};

#define NOLIGOINDICES_MINOR 3
static int indexsizes_minor[NOLIGOINDICES_MINOR] = {8, 7, 6};
static Shortoligomer_T masks_minor[NOLIGOINDICES_MINOR] = {STRAIGHT_MASK_8, STRAIGHT_MASK_7, STRAIGHT_MASK_6};
static int diag_lookbacks_minor[NOLIGOINDICES_MINOR] = {120, 60, 30};
static int suffnconsecutives_minor[NOLIGOINDICES_MINOR] = {20, 15, 10};

#else

#define NOLIGOINDICES_MAJOR 3
static int indexsizes_major[NOLIGOINDICES_MAJOR] = {8, 7, 6};
static Shortoligomer_T masks_major[NOLIGOINDICES_MAJOR] = {STRAIGHT_MASK_8, STRAIGHT_MASK_7, STRAIGHT_MASK_6};
static int diag_lookbacks_major[NOLIGOINDICES_MAJOR] = {120, 60, 30};
static int suffnconsecutives_major[NOLIGOINDICES_MAJOR] = {20, 15, 10};

#define NOLIGOINDICES_MINOR 3
static int indexsizes_minor[NOLIGOINDICES_MINOR] = {8, 7, 6};
static Shortoligomer_T masks_minor[NOLIGOINDICES_MINOR] = {STRAIGHT_MASK_8, STRAIGHT_MASK_7, STRAIGHT_MASK_6};
static int diag_lookbacks_minor[NOLIGOINDICES_MINOR] = {120, 60, 30};
static int suffnconsecutives_minor[NOLIGOINDICES_MINOR] = {20, 15, 10};

#endif

#if 0
/* This fails badly on NM_001142699 */
#define NOLIGOINDICES_MINOR 1
static int indexsizes_minor[NOLIGOINDICES_MINOR] = {4};
static Shortoligomer_T masks_minor[NOLIGOINDICES_MINOR] = {STRAIGHT_MASK_4};
static int diag_lookbacks_minor[NOLIGOINDICES_MINOR] = {30};
static int suffnconsecutives_minor[NOLIGOINDICES_MINOR] = {5};
#endif


#define THETADIFF1 20.0
#define THETADIFF2 20.0
#define REPOLIGOCOUNT 8

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


/* Treats 'N' as 'A', just as Oligoindex_hr_tally does */
/* #define EXTRACT_GENOMICSEG 1 */


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

/* Diagonals */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* dump_positions */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif


#ifdef PMAP
/* A short oligomer, representing 3 amino acids, or 9 nt, requiring 18
   bits or 2.25 bytes */
#else
/* A short oligomer, typically 8 nt (but could by 9 or 10 nt),
   requiring 16 bits or 2 bytes */
#endif

#define T Oligoindex_T


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

int
Oligoindex_indexsize (T this) {
#ifdef PMAP
  return this->indexsize_aa;
#else
  return this->indexsize;
#endif
}


static T
Oligoindex_new (int indexsize, int diag_lookback, int suffnconsecutive
#ifndef PMAP
		, Shortoligomer_T mask
#endif
		) {
  T new = (T) MALLOC(sizeof(*new));

#ifdef PMAP
  new->indexsize_aa = indexsize;
  new->msb = power(NAMINOACIDS_STAGE2,indexsize-1);
  debug(printf("msb for indexsize %d is %u\n",indexsize,new->msb));
  new->oligospace = power(NAMINOACIDS_STAGE2,indexsize);
#else
  new->indexsize = indexsize;
  new->mask = mask;
  new->oligospace = power(4,indexsize);
#endif

  new->diag_lookback = diag_lookback;
  new->suffnconsecutive = suffnconsecutive;

  new->query_evaluated_p = false;
  new->overabundant = (bool *) CALLOC(new->oligospace,sizeof(bool));
  new->inquery = (bool *) CALLOC(new->oligospace,sizeof(bool));
  new->counts = (int *) CALLOC(new->oligospace,sizeof(int));
  new->relevant_counts = (int *) CALLOC(new->oligospace,sizeof(int));
  new->positions = (Genomicpos_T **) CALLOC(new->oligospace,sizeof(Genomicpos_T *));
  new->pointers = (Genomicpos_T **) CALLOC(new->oligospace,sizeof(Genomicpos_T *));

  return new;
}


int *
Oligoindex_counts_copy (T this) {
  int *counts;

  counts = (int *) CALLOC(this->oligospace,sizeof(int));
  memcpy(counts,this->counts,this->oligospace*sizeof(int));
  return counts;
}


T *
Oligoindex_new_major (int *noligoindices) {
  T *oligoindices;
  int source;

  *noligoindices = NOLIGOINDICES_MAJOR;
  oligoindices = (T *) CALLOC(NOLIGOINDICES_MAJOR,sizeof(T));
  for (source = 0; source < NOLIGOINDICES_MAJOR; source++) {
#ifdef PMAP
    oligoindices[source] = Oligoindex_new(indexsizes_major[source],diag_lookbacks_major[source],suffnconsecutives_major[source]);
#else
    oligoindices[source] = Oligoindex_new(indexsizes_major[source],diag_lookbacks_major[source],suffnconsecutives_major[source],masks_major[source]);
#endif
  }
  
  return oligoindices;
}

T *
Oligoindex_new_minor (int *noligoindices) {
  T *oligoindices;
  int source;

  *noligoindices = NOLIGOINDICES_MINOR;
  oligoindices = (T *) CALLOC(NOLIGOINDICES_MINOR,sizeof(T));
  for (source = 0; source < NOLIGOINDICES_MINOR; source++) {
#ifdef PMAP
    oligoindices[source] = Oligoindex_new(indexsizes_minor[source],diag_lookbacks_minor[source],suffnconsecutives_minor[source]);
#else
    oligoindices[source] = Oligoindex_new(indexsizes_minor[source],diag_lookbacks_minor[source],suffnconsecutives_minor[source],masks_minor[source]);
#endif
  }
  
  return oligoindices;
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


#ifndef PMAP
void
Oligoindex_dump (T this) {
  Oligospace_T i;
  char *nt;

  for (i = 0; i < this->oligospace; i++) {
    if (this->counts[i] == 0) {
    } else {
      nt = shortoligo_nt(i,this->indexsize);
      printf("%s %d\n",nt,this->counts[i]);
      FREE(nt);
    }
  }
  return;
}


void
Oligoindex_counts_dump (T this, int *counts) {
  Oligospace_T i;
  char *nt;

  for (i = 0; i < this->oligospace; i++) {
    if (counts[i] == 0 && this->counts[i] == 0) {
    } else {
      nt = shortoligo_nt(i,this->indexsize);
      printf("%s %d %d\n",nt,counts[i],this->counts[i]);
      FREE(nt);
    }
  }
  return;
}


bool
Oligoindex_counts_equal (T this, int *counts) {
  Oligospace_T i;
  char *nt;

  for (i = 0; i < this->oligospace; i++) {
    if (counts[i] != this->counts[i]) {
      nt = shortoligo_nt(i,this->indexsize);
      fprintf(stderr,"At %s, counts = %d, this->counts = %d\n",
	      nt,counts[i],this->counts[i]);
      FREE(nt);
      /* Oligoindex_counts_dump(this,counts); */
      return false;
    }
  }
  return true;
}
#endif



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

#ifndef PMAP
/* -1 means edge is on 5', +1 means edge is on 3'  0 means no edge found. */
static int
edge_detect (int *edge, int *sumx, int *sumxx, int length) {
  int side = 0;
  int pos, sumx_left, sumxx_left, sumx_right, sumxx_right, n_left, n_right;
  double theta, sumx_pseudo, theta_left, theta_right, rss_left, rss_right, rss_sep;
  double min_rss_sep;
#ifdef DEBUG1
  double fscore;
#endif

  debug1(printf("\n*** Start of edge_detect\n"));

  sumx_right = sumx[length] - sumx[0];
  sumxx_right = sumxx[length] - sumxx[0];

  theta = (double) sumx_right/(double) length;
  sumx_pseudo = NPSEUDO * theta;
  min_rss_sep = sumxx_right - sumx_right*theta;
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

    debug1(
	   if (rss_sep > 0.0) {
	     fscore = ((double) (length - 2))*(rss - rss_sep)/rss_sep;
	     printf("%d %d %d %d %d %d %f %f %f %f %f\n",
		    pos,sumx[pos]-sumx[pos-1],sumx_left,n_left,sumx_right,n_right,
		    theta_left,theta_right,rss_left,rss_right,fscore);
	   } else {
	     printf("%d %d %d %d %d %d %f %f %f %f NA\n",
		    pos,sumx[pos]-sumx[pos-1],sumx_left,n_left,sumx_right,n_right,
		    theta_left,theta_right,rss_left,rss_right);
	   });
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
#endif

#ifndef PMAP
static int
trim_start_detect (int start, int end, int *sumx, int *sumxx) {
  int edge = -1;
  int pos, sumx_left, sumxx_left, sumx_right, sumxx_right, n_left, n_right;
  double theta, sumx_pseudo, theta_left, theta_right, rss_left, rss_right, rss_sep;
  double min_rss_sep;
#ifdef DEBUG1
  double fscore;
#endif

  debug1(printf("\n*** Start of trim_start_detect\n"));

  sumx_right = sumx[end] - sumx[start];
  sumxx_right = sumxx[end] - sumxx[start];

  if (end <= start) {
    return -1;
  }
  theta = (double) sumx_right/(double) (end - start);
  sumx_pseudo = NPSEUDO * theta;
  min_rss_sep = sumxx_right - sumx_right*theta;
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

    debug1(
	   if (rss_sep > 0.0) {
	     fscore = ((double) (end - start - 2))*(rss - rss_sep)/rss_sep;
	     printf("%d %d %d %d %d %f %f %f %f %f\n",
		    pos,sumx_left,n_left,sumx_right,n_right,
		    theta_left,theta_right,rss_left,rss_right,fscore);
	   } else {
	     printf("%d %d %d %d %d %f %f %f %f NA\n",
		    pos,sumx_left,n_left,sumx_right,n_right,
		    theta_left,theta_right,rss_left,rss_right);
	   });
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
#endif

#ifndef PMAP
static int
trim_end_detect (int start, int end, int *sumx, int *sumxx) {
  int edge = -1;
  int pos, sumx_left, sumxx_left, sumx_right, sumxx_right, n_left, n_right;
  double theta, sumx_pseudo, theta_left, theta_right, rss_left, rss_right, rss_sep;
  double min_rss_sep;
#ifdef DEBUG1
  double fscore;
#endif

  debug1(printf("\n*** Start of trim_end_detect\n"));

  sumx_right = sumx[end] - sumx[start];
  sumxx_right = sumxx[end] - sumxx[start];

  if (end <= start) {
    return -1;
  }
  theta = (double) sumx_right/(double) (end - start);
  sumx_pseudo = NPSEUDO * theta;
  min_rss_sep = sumxx_right - sumx_right*theta;
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

    debug1(
	   if (rss_sep == 0) {
	     printf("%d %d %d %d %d %f %f %f %f NA\n",
		    pos,sumx_left,n_left,sumx_right,n_right,
		    theta_left,theta_right,rss_left,rss_right);
	   } else {
	     fscore = ((double) (end - start - 2))*(rss - rss_sep)/rss_sep;
	     printf("%d %d %d %d %d %f %f %f %f %f\n",
		    pos,sumx_left,n_left,sumx_right,n_right,
		    theta_left,theta_right,rss_left,rss_right,fscore);
	   });
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
#endif


/* Run query sequence through this procedure.  First, we count occurrences
 * of each oligo in queryuc (upper case version of queryseq).  This
 * allows us to scan genomicseg intelligently, because then we know
 * whether to store positions for that oligo. */

double
Oligoindex_set_inquery (int *badoligos, int *repoligos, int *trimoligos, int *trim_start, int *trim_end,
			T this, char *queryuc_ptr, int querylength, bool trimp) {
#ifndef PMAP
  double oligodepth;
  int ngoodoligos, nrepoligos, x, *sumx, *sumxx, sumx0 = 0, sumxx0 = 0;
  int edge, side;
  int querypos;
  Shortoligomer_T oligo = 0U;
  char *ptr;
#endif
  int nunique = 0;
  int i, noligos = 0;
  int in_counter = 0;
  Shortoligomer_T masked;
  char *p;
#ifdef PMAP
  int indexsize = this->indexsize_aa, index;
  unsigned int aaindex = 0U;
#else
  int indexsize = this->indexsize;
#endif
#if 0
  bool prevunique = false;
#endif
#ifdef DEBUG
  char *aa, *nt;
#endif

  if (this->query_evaluated_p == true) {
    return 1.0;
  } else {
    this->query_evaluated_p = true; /* Set this flag so we don't redo this part */
  }

  memset((void *) this->inquery,false,this->oligospace*sizeof(bool));
  memset((void *) this->counts,0,this->oligospace*sizeof(int));
#if 0
  /* Test for thread safety */
  for (i = 0; i < this->oligospace; i++) {
    if (this->inquery[i] != false) {
      abort();
    }
  }
  for (i = 0; i < this->oligospace; i++) {
    if (this->counts[i] != 0) {
      abort();
    }
  }
#endif

  if (querylength <= indexsize) {
    *badoligos = 0;
    *trim_start = 0;
    *trim_end = querylength;
    return 1.0;
  }
    
  for (i = 0, p = queryuc_ptr; i < querylength; i++, p++) {
    in_counter++;

#ifdef PMAP
    if ((index = aa_index_table[(int) *p]) < 0) {
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
    *badoligos = (querylength + 1 - indexsize) - noligos;
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
  for (i = 0, p = queryuc_ptr; i < querylength; i++, p++) {
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
      x = this->counts[oligo & this->mask];
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
  for (i = 0, p = queryuc_ptr; i < querylength; i++, p++) {
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
      if (this->counts[oligo & this->mask] == 1) {
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
  ptr = queryuc_ptr;
  p = &(ptr[*trim_start]);
  for (querypos = (*trim_start)-indexsize; querypos < (*trim_end)-indexsize;
       querypos++, p++) {
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
      if (this->counts[oligo & this->mask] >= REPOLIGOCOUNT) {
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
		    char *sequence, int seqlength, int sequencepos) {
#ifdef PMAP
  int index;
  int frame = 2;
  unsigned int aaindex0 = 0U, aaindex1 = 0U, aaindex2 = 0U;
  int in_counter_0 = 0, in_counter_1 = 0, in_counter_2 = 0;
#endif
  int i = 0, n;
  int in_counter = 0;
  Shortoligomer_T oligo = 0U, masked;
  char *p;
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

  sequencepos -= indexsize;
  for (i = 0, p = sequence; i < seqlength; i++, p++) {
#ifdef PMAP
    in_counter_0++;
    in_counter_1++;
    in_counter_2++;
#else
    in_counter++;
#endif
    sequencepos++;

    switch (*p) {
    case 'A':
#ifdef EXTRACT_GENOMICSEG
    case 'N':
#endif
      oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; 
#ifdef PMAP
      in_counter_0 = in_counter_1 = in_counter_2 = 0; 
#else
      in_counter = 0;
#endif
    }

    debug(printf("At genomicpos %u, char is %c, oligo is %04X\n",
		 sequencepos,*p,oligo));

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
      debug(printf("%04X\n",masked));
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
	      printf("At genomicpos %u, oligo %s wasn't seen in querypos\n",sequencepos,nt);
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

  return totalcounts;
}


/* Third, run genomicuc through this procedure, to determine the genomic positions for each oligo.  */
/* Logic of this procedure should match that of allocate_positions */
static int
store_positions (Genomicpos_T **pointers, bool *overabundant, 
		 bool *inquery, Oligospace_T oligospace, int indexsize,
#ifdef PMAP
		 Shortoligomer_T msb,
#else
		 Shortoligomer_T mask,
#endif
		 char *sequence, int seqlength, int sequencepos) {
#ifdef PMAP
  int index;
  int frame = 2;
  unsigned int aaindex0 = 0U, aaindex1 = 0U, aaindex2 = 0U;
  int in_counter_0 = 0, in_counter_1 = 0, in_counter_2 = 0;
#endif
  int nstored = 0;
  int i = 0;
  int in_counter = 0;
  Shortoligomer_T oligo = 0U, masked;
  char *p;

#ifdef DEBUG
#ifdef PMAP
  char *aa;
#else
  char *nt;
#endif
#endif

  sequencepos -= indexsize;
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
    case 'A':
#ifdef EXTRACT_GENOMICSEG
    case 'N':
#endif
      oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; 
#ifdef PMAP
      in_counter_0 = in_counter_1 = in_counter_2 = 0; 
#else
      in_counter = 0;
#endif
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
	if (masked >= oligospace) {
	  abort();
	}
	pointers[masked][0] = (Genomicpos_T) sequencepos;
	pointers[masked]++;
	nstored++;
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

  return nstored;
}


#ifdef DEBUG9
static void
dump_positions (Genomicpos_T **positions, int *counts, int oligospace, int indexsize) {
  int i;
#ifdef PMAP
  char *aa;
#else
  char *nt;
#endif

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
#endif


#ifndef PMAP
#define POLY_A 0x0000
#define POLY_C 0x5555
#define POLY_G 0xAAAA
#define POLY_T 0xFFFF
#endif

void
Oligoindex_tally (T this, char *genomicuc_trimptr, int genomicuc_trimlength,
		  char *queryuc_ptr, int querylength, int sequencepos) {
  int badoligos, repoligos, trimoligos, trim_start, trim_end;
  int nallocated, nstored;

  Oligoindex_set_inquery(&badoligos,&repoligos,&trimoligos,&trim_start,&trim_end,this,
			 queryuc_ptr,querylength,/*trimp*/false);

  memset((void *) this->counts,0,this->oligospace*sizeof(int));
  memset((void *) this->overabundant,false,this->oligospace*sizeof(bool));
#if 0
  /* Test for thread safety */
  for (i = 0; i < this->oligospace; i++) {
    if (this->counts[i] != 0) {
      abort();
    }
  }
  for (i = 0; i < this->oligospace; i++) {
    if (this->overabundant[i] != false) {
      abort();
    }
  }
#endif

#ifndef PMAP
  /* These values will prevent oligoindex from getting mappings later */
  this->overabundant[POLY_A & this->mask] = true;
  this->overabundant[POLY_C & this->mask] = true;
  this->overabundant[POLY_G & this->mask] = true;
  this->overabundant[POLY_T & this->mask] = true;
#endif

  debug9(printf("sequencepos is %d\n",sequencepos));
  if ((nallocated = allocate_positions(this->pointers,this->positions,this->overabundant,
				       this->inquery,this->counts,this->relevant_counts,this->oligospace,
#ifdef PMAP
				       3*this->indexsize_aa,this->msb,
#else
				       this->indexsize,this->mask,
#endif
				       &(genomicuc_trimptr[sequencepos]),genomicuc_trimlength,
				       sequencepos)) > 0) {

    nstored = store_positions(this->pointers,this->overabundant,this->inquery,this->oligospace,
#ifdef PMAP
			      3*this->indexsize_aa,this->msb,
#else
			      this->indexsize,this->mask,
#endif
			      &(genomicuc_trimptr[sequencepos]),genomicuc_trimlength,
			      sequencepos);

#ifdef PMAP
    debug9(dump_positions(this->positions,this->counts,this->oligospace,3*this->indexsize_aa));
#else
    debug9(dump_positions(this->positions,this->counts,this->oligospace,this->indexsize));
#endif

    if (nstored != nallocated) {
      fprintf(stderr,"Bug in Oligoindex_tally: %d allocated, but %d stored\n",
	      nallocated,nstored);
      abort();
    }
  }

  return;
}
  
void
Oligoindex_clear_inquery (T this) {

  this->query_evaluated_p = false;
  return;
}

void
Oligoindex_untally (T this) {

  /* The following line is necessary */
  Oligoindex_clear_inquery(this);
  if (this->positions[0] != NULL) {
    FREE(this->positions[0]);
  }

  return;
}



static void
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

void
Oligoindex_free_array (T **oligoindices, int noligoindices) {
  int source;

  for (source = 0; source < noligoindices; source++) {
    Oligoindex_free(&((*oligoindices)[source]));
  }
  FREE(*oligoindices);
  return;
}


static Genomicpos_T *
lookup (int *nhits, T this, Shortoligomer_T masked) {
#ifdef DEBUG
  char *aa, *nt;
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


#if 0
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
#endif


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
List_T
Oligoindex_get_mappings (List_T diagonals,
			 bool *coveredp, Genomicpos_T **mappings, int *npositions,
			 int *totalpositions, bool *oned_matrix_p, int *maxnconsecutive, 
			 T this, char *queryuc_ptr, int querylength,
			 Genomicpos_T chrstart, Genomicpos_T chrend,
			 Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp,
			 Diagpool_T diagpool) {
  int nhits, hit, diagi_adjustment, i;
  Genomicpos_T diagi;
  int diag_lookback, suffnconsecutive;
#ifdef PREV_MAXCONSECUTIVE
  int prev_querypos, prev_nhits;
  unsigned int *prev_mappings;
  int ngoodconsecutive;
#endif
#ifndef PMAP
  Shortoligomer_T oligo = 0U;
#endif

  char *p;
  int in_counter = 0, querypos;
  Shortoligomer_T masked;
  Genomicpos_T genomiclength, chrinit;

  void *item;
  struct Genomicdiag_T *genomicdiag;
  Genomicdiag_T ptr;
  List_T good_genomicdiags = NULL;

#ifdef PMAP
  int indexsize = this->indexsize_aa, index;
  unsigned int aaindex = 0U;
#else
  int indexsize = this->indexsize;
#endif

  diag_lookback = this->diag_lookback;
  suffnconsecutive = this->suffnconsecutive;
  genomiclength = chrend - chrstart;
  if (plusp == true) {
    chrinit = chrstart;
  } else {
    chrinit = (chrhigh - chroffset) - chrend;
  }

#ifdef PMAP
  genomicdiag = (struct Genomicdiag_T *) CALLOC(3*querylength+genomiclength+1,sizeof(struct Genomicdiag_T));
#else
  genomicdiag = (struct Genomicdiag_T *) CALLOC(querylength+genomiclength+1,sizeof(struct Genomicdiag_T));
#endif

#if 0
  /* Too time consuming.  Just initialize when we see [diagi] for the first time. */
#ifdef PMAP
  for (diagi = 0; diagi < 3*querylength+genomiclength; diagi++) {
    genomicdiag[diagi].i = diagi;
    genomicdiag[diagi].querypos = -diag_lookback; /* guarantees first check won't be consecutive */
  }
#else
  for (diagi = 0; diagi < querylength+genomiclength; diagi++) {
    genomicdiag[diagi].i = diagi;
    genomicdiag[diagi].querypos = -diag_lookback; /* guarantees first check won't be consecutive */
  }
#endif
#endif


  querypos = -indexsize;
  *oned_matrix_p = true;
  for (i = 0, p = queryuc_ptr; i < querylength; i++, p++) {
    in_counter++;
    querypos++;
    
#ifdef PMAP
    if ((index = aa_index_table[(int) *p]) < 0) {
      aaindex = 0U;
      in_counter = 0;
    } else {
      aaindex = aaindex % this->msb;
      aaindex = aaindex * NAMINOACIDS_STAGE2 + index;
    }
#else
    switch (*p) {
    case 'A':
#ifdef EXTRACT_GENOMICSEG
    case 'N':
#endif
      oligo = (oligo << 2); break;
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
      if (coveredp[querypos] == false) {
	mappings[querypos] = lookup(&nhits,this,masked);
	npositions[querypos] = nhits;
	if (nhits > 0) {
	  *totalpositions += nhits;
	  if (*totalpositions < 0) {
	    /* fprintf(stderr,"totalpositions %d is negative for masked oligo %u\n",*totalpositions,masked); */
	    *oned_matrix_p = false;
	  }

#ifdef PMAP
	  /* diagonal is (position - 3*querypos); diagi is (position - 3*querypos) + 3*querylength */
	  diagi_adjustment = 3*(querylength - querypos);
#else
	  /* diagonal is (position - querypos); diagi is (position - querypos) + querylength */
	  diagi_adjustment = querylength - querypos;
#endif
	  for (hit = 0; hit < nhits; hit++) {
	    diagi = mappings[querypos][hit] + diagi_adjustment - chrinit;
	    ptr = &(genomicdiag[diagi]);

#ifdef PMAP
	    assert(diagi <= 3*querylength+genomiclength);
#else
	    assert(diagi <= querylength+genomiclength);
#endif

	    if (ptr->i == 0) {
	      /* Initialize */
	      ptr->i = diagi;
	      ptr->querypos = -diag_lookback; /* guarantees first check won't be consecutive */
	    }

	    /* Must use >= here, so querypos 0 - (-diag_lookback) will fail */
	    if (querypos - ptr->querypos >= diag_lookback) {
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
	      if (ptr->best_nconsecutive == suffnconsecutive) {
		/* Need to check for ==, not >=, because this will store the ptr once */
		good_genomicdiags = List_push(good_genomicdiags,(void *) ptr);
	      }
	      if (ptr->best_nconsecutive > *maxnconsecutive) {
		*maxnconsecutive = ptr->best_nconsecutive;
	      }
	    }
	    ptr->querypos = querypos;
	  }
	}
      }
      in_counter--;
    }
  }

  while (good_genomicdiags != NULL) {
    good_genomicdiags = List_pop(good_genomicdiags,&item);
    ptr = (Genomicdiag_T) item;
#ifdef USE_DIAGPOOL
#ifdef PMAP
    diagonals = Diagpool_push(diagonals,diagpool,/*diagonal*/(ptr->i - 3*querylength),
			      ptr->best_consecutive_start,ptr->best_consecutive_end,
			      ptr->best_nconsecutive+1);
#else
    diagonals = Diagpool_push(diagonals,diagpool,/*diagonal*/(ptr->i - querylength),
			      ptr->best_consecutive_start,ptr->best_consecutive_end,
			      ptr->best_nconsecutive+1);
#endif
#else
#ifdef PMAP
    diagonals = List_push(diagonals,(void *) Diag_new(/*diagonal*/(ptr->i - 3*querylength),
						      ptr->best_consecutive_start,ptr->best_consecutive_end,
						      ptr->best_nconsecutive+1));
#else
    diagonals = List_push(diagonals,(void *) Diag_new(/*diagonal*/(ptr->i - querylength),
						      ptr->best_consecutive_start,ptr->best_consecutive_end,
						      ptr->best_nconsecutive+1));
#endif
#endif
  }

  FREE(genomicdiag);

  return diagonals;
}


