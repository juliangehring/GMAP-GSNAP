static char rcsid[] = "$Id: oligoindex.c,v 1.61 2006/03/05 03:20:46 twu Exp $";
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

#ifdef PMAP
#define NAMINOACIDS_STAGE2 20
#define NT_PER_MATCH 3
#else
#define NT_PER_MATCH 1
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
  bool *overabundant;
  bool *inquery;
  int *counts;
  Genomicpos_T **positions;
  Genomicpos_T *data;
};

static int oligospace;
#ifdef PMAP
static int indexsize_aa;
static int indexsize_nt;
static Shortoligomer_T msb;
#else
static int indexsize;
static Shortoligomer_T mask;
#endif

#define THETADIFF1 2.0
#define THETADIFF2 2.0
#define BADOLIGOCOUNT 8
#define ALLOCSIZE 200
#define MAXHITS 200		/* Certain code below depends on whether this is greater than ALLOCSIZE */
#define MAXHITS_GT_ALLOCSIZE 0


#ifdef PMAP
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

#endif


static int
power (int base, int exponent) {
  int result = 1, i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}

void
Oligoindex_init (int indexsize0) {
#ifdef PMAP
  indexsize_nt = indexsize0;
  indexsize_aa = indexsize_nt / 3;
  msb = power(NAMINOACIDS_STAGE2,indexsize_aa-1);
  oligospace = power(NAMINOACIDS_STAGE2,indexsize_aa);
#else
  indexsize = indexsize0;
  mask = ~(~0U << 2*indexsize);
  oligospace = power(4,indexsize);
#endif
  return;
}

T
Oligoindex_new (void) {
  T new = (T) MALLOC(sizeof(*new));

  new->overabundant = (bool *) CALLOC(oligospace,sizeof(bool));
  new->inquery = (bool *) CALLOC(oligospace,sizeof(bool));
  new->counts = (int *) CALLOC(oligospace,sizeof(int));
  new->positions = (Genomicpos_T **) CALLOC(oligospace,sizeof(Genomicpos_T *));
  new->data = (Genomicpos_T *) CALLOC(oligospace*ALLOCSIZE,sizeof(Genomicpos_T));
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
static char aa_table[NAMINOACIDS_STAGE2] = "ACDEFGHIKLMNPQRSTVWY";

static char *
aaindex_aa (unsigned int aaindex) {
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

  sumx_right = sumx[length] - sumx[0];
  sumxx_right = sumxx[length] - sumxx[0];

  theta = (double) sumx_right/(double) length;
  sumx_pseudo = NPSEUDO * theta;
  min_rss_sep = rss = sumxx_right - sumx_right*theta;
  debug1(printf("theta: %d/%d = %f\n",sumx_right,length,theta));
  debug1(printf("rss: %f\n",rss)); 

  debug1(printf("%s %s %s %s %s %s %s %s %s %s\n",
		"pos","sumx.left","n.left","sumx.right","n.right",
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
      debug1(printf("%d %d %d %d %d %f %f %f %f NA\n",
		    pos,sumx_left,n_left,sumx_right,n_right,
		    theta_left,theta_right,rss_left,rss_right));
    } else {
      debug1(      
	     fscore = ((double) (length - 2))*(rss - rss_sep)/rss_sep;
	     debug1(printf("%d %d %d %d %d %f %f %f %f %f\n",
			   pos,sumx_left,n_left,sumx_right,n_right,
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

  return side;
}

static int
trim_start_detect (int start, int end, int *sumx, int *sumxx) {
  int edge = -1;
  int pos, sumx_left, sumxx_left, sumx_right, sumxx_right, n_left, n_right;
  double theta, sumx_pseudo, x_pseudo, theta_left, theta_right, rss, rss_left, rss_right, rss_sep;
  double fscore, min_rss_sep;

  sumx_right = sumx[end] - sumx[start];
  sumxx_right = sumxx[end] - sumxx[start];

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

  debug1(printf("trim_detect returning %d\n",edge));
  return edge;
}

static int
trim_end_detect (int start, int end, int *sumx, int *sumxx) {
  int edge = -1;
  int pos, sumx_left, sumxx_left, sumx_right, sumxx_right, n_left, n_right;
  double theta, sumx_pseudo, x_pseudo, theta_left, theta_right, rss, rss_left, rss_right, rss_sep;
  double fscore, min_rss_sep;

  sumx_right = sumx[end] - sumx[start];
  sumxx_right = sumxx[end] - sumxx[start];

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

  debug1(printf("trim_detect returning %d\n",edge));
  return edge;
}



/* Run query sequence through this procedure.  First, we count occurrences
 * of each oligo in queryuc (upper case version of queryseq).  This
 * allows us to scan genomicseg intelligently, because then we know
 * whether to store positions for that oligo. */

double
Oligoindex_set_inquery (int *badoligos, int *trimoligos, int *trim_start, int *trim_end,
			T this, Sequence_T queryuc, bool trimp) {
  double oligodepth;
  int nunique = 0, ngoodoligos = 0, x, *sumx, *sumxx, sumx0 = 0, sumxx0 = 0, querylength;
  int i = 0, noligos = 0, edge, side;
  int in_counter = 0, querypos;
  Shortoligomer_T oligo = 0U, masked;
  char *p, *ptr;
  bool prevunique = false;
#ifdef PMAP
  int indexsize = indexsize_aa, index;
  unsigned int aaindex = 0U;
  char *aa;
#else
  char *nt;
#endif

  memset((void *) this->inquery,false,oligospace*sizeof(bool));
  memset((void *) this->counts,0,oligospace*sizeof(int));

  querylength = Sequence_fulllength(queryuc);
  for (i = 0, p = Sequence_fullpointer(queryuc); i < querylength; i++, p++) {
    in_counter++;

#ifdef PMAP
    if ((index = aa_index_table[*p]) < 0) {
      aaindex = 0U;
      in_counter = 0;
    } else {
      aaindex = aaindex % msb;
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
      debug(aa = aaindex_aa(masked);
	    printf("At querypos %d, aa %s seen (masked = %u)\n",i,aa,masked);
	    FREE(aa));
#else
      masked = oligo & mask;
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

  /*
  printf("Saw %d oligos out of %d expected\n",noligos,seqlength-indexsize+1);
  */

#ifdef CHANGEPOINT
  /* Determine where to trim */
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
      x = this->counts[oligo&mask];
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
  *trim_end = querylength-indexsize;
  if ((side = edge_detect(&edge,sumx,sumxx,querylength-indexsize)) == -1) {
    *trim_start = edge;
    if ((edge = trim_end_detect(*trim_start,querylength-indexsize,sumx,sumxx)) >= 0) {
      *trim_end = edge;
    }
  } else if (side == +1) {
    *trim_end = edge;
    if ((edge = trim_start_detect(0,*trim_end,sumx,sumxx)) >= 0) {
      *trim_start = edge;
    }
  }

  FREE(sumxx);
  FREE(sumx);

  *trim_end += indexsize;

  debug1(printf("trim_start = %d, trim_end = %d\n",*trim_start,*trim_end));
#else

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
      if (this->counts[oligo&mask] == 1) {
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

#ifdef COUNTBADOLIGOS
  /* Count good oligos */
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
      if (this->counts[oligo&mask] < BADOLIGOCOUNT) {
	ngoodoligos++;
      }
      in_counter--;
    }
  }

  *trimoligos = (*trim_end - indexsize) - (*trim_start) + 1;
  *badoligos = (*trimoligos) - ngoodoligos;

  if (nunique == 0) {
    return 1000000.0;
  } else {
    /* trimlength = *trim_end - *trim_start; */
    oligodepth = (double) noligos/(double) nunique;
    return oligodepth;
  }
#else
  *trimoligos = (*trim_end - indexsize) - (*trim_start) + 1;
  *badoligos = 0;
  return 1000000.0;
#endif

#endif
}


/************************************************************************/

#ifdef PMAP

/*              87654321 */
#define CHAR1 0x00000030
#define CHAR2 0x0000000C
#define CHAR3 0x00000003

#define A1    0x00000000
#define C1    0x00000010
#define G1    0x00000020
#define T1    0x00000030

#define A2    0x00000000
#define C2    0x00000004
#define G2    0x00000008
#define T2    0x0000000C

#define A3    0x00000000
#define C3    0x00000001
#define G3    0x00000002
#define T3    0x00000003

/* Do not include stop codon, because oligospace will not allow it */
typedef enum {AA_A, AA_C, AA_D, AA_E, AA_F, 
	      AA_G, AA_H, AA_I, AA_K, AA_L,
	      AA_M, AA_N, AA_P, AA_Q, AA_R,
	      AA_S, AA_T, AA_V, AA_W, AA_Y} Aminoacid_T;

static int
get_last_codon (Shortoligomer_T oligo) {
  switch (oligo & CHAR2) {
  case T2:
    switch (oligo & CHAR1) {
    case T1:
      switch (oligo & CHAR3) {
      case T3: case C3: return AA_F;
      default: /* case A3: case G3: */ return AA_L;
      }
    case C1: return AA_L;
    case A1: 
      switch (oligo & CHAR3) {
      case G3: return AA_M;
      default: /* case T3: case A3: case C3: */ return AA_I;
      }
    default: /* case G1: */ return AA_V;
  }
  case C2:
    switch (oligo & CHAR1) {
    case T1: return AA_S;
    case C1: return AA_P;
    case A1: return AA_T;
    default: /* case G1: */ return AA_A;
    }
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
      default: /* case A3: case G3: */ return AA_K;
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
    case C1: return AA_R;
    case A1:
      switch (oligo & CHAR3) {
      case T3: case C3: return AA_S;
      default: /* case A3: case G3: */ return AA_R;
      }
    default: /* case G1: */ return AA_G;
    }
  }

  abort();
  return -1;
}

#endif


/************************************************************************/



/* Run genomicuc through this procedure */
static void
store_positions (Genomicpos_T **positions, Genomicpos_T *data, 
		 bool *overabundant, bool *inquery,
		 int *counts, char *sequence, int seqlength) {
#ifdef PMAP
  int indexsize = indexsize_nt, index;
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
#if MAXHITS_GT_ALLOCSIZE
  Genomicpos_T *ptr;
#endif

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
	debug(aa = aaindex_aa(masked);
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
	debug(aa = aaindex_aa(masked);
	      printf("At genomicpos %u, aa %s wasn't seen in querypos\n",
		     sequencepos,aa,counts[masked]);
	      FREE(aa));
#else
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s wasn't seen in querypos\n",sequencepos,nt,counts[masked]);
	      FREE(nt));
#endif

      } else if ((count = counts[masked]) == 0) {
	/* Initial allocation; use pre-allocated space */
	positions[masked] = &(data[availslot]);
	availslot += ALLOCSIZE;
	positions[masked][0] = (Genomicpos_T) sequencepos;
	counts[masked] = 1;
#ifdef PMAP
	debug(aa = aaindex_aa(masked);
	      printf("At genomicpos %u, aa %s seen, counts is now %d\n",
		     sequencepos,aa,counts[masked]);
	      FREE(aa));
#else
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s seen, counts is now %d\n",sequencepos,nt,counts[masked]);
	      FREE(nt));
#endif

      } else if (count >= MAXHITS) {
	/* Overabundant: Turn this oligo off */
	overabundant[masked] = true;
	counts[masked] = 0;	/* This is okay as long as check for overabundant is first */
	/* This also prevents allocated memory from being freed later */
	/* Free previously allocated new space */
	/* if (count > ALLOCSIZE) { -- Should always be true if MAXHITS > ALLOCSIZE */
#if MAXHITS_GT_ALLOCSIZE
	FREE(positions[masked]);
#endif

#if MAXHITS_GT_ALLOCSIZE
      } else if ((count % ALLOCSIZE) == 0) {
	/* Change to new space */
	/*
	  positions[masked] = 
	  (Genomicpos_T *) RESIZE(positions[masked],(count+ALLOCSIZE)*sizeof(Genomicpos_T)); 
	*/
	ptr = (Genomicpos_T *) CALLOC(count+ALLOCSIZE,sizeof(Genomicpos_T));
	memcpy((void *) ptr,positions[masked],count*sizeof(Genomicpos_T));
	if (count > ALLOCSIZE) {
	  /* Free previously allocated new space */
	  FREE(positions[masked]);
	}

	positions[masked] = ptr;
	positions[masked][count] = (Genomicpos_T) sequencepos;
	counts[masked] += 1;
#ifdef PMAP
	debug(aa = aaindex_aa(masked);
	      printf("At genomicpos %u, aa %s seen, counts is now %d\n",
		     sequencepos,aa,counts[masked]);
	      FREE(aa));
#else
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s seen, counts is now %d\n",sequencepos,nt,counts[masked]);
	      FREE(nt));

#endif
#endif	/* MAXHITS_GT_ALLOCSIZE */
      } else {
	positions[masked][count] = (Genomicpos_T) sequencepos;
	counts[masked] += 1;
#ifdef PMAP
	debug(aa = aaindex_aa(masked);
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

#ifdef PMAP
  debug(
	for (i = 0; i < oligospace; i++) {
	  aa = aaindex_aa(i);
	  if (counts[i] >= 1) {
	    printf("AA %s => %d entries: %u...%u\n",
		   aa,counts[i],positions[i][0],positions[i][counts[i]-1]);
	  }
	  FREE(aa);
	}
	);
#else
  debug(
	for (i = 0; i < oligospace; i++) {
	  nt = shortoligo_nt(i,indexsize);
	  if (counts[i] >= 1) {
	    printf("Oligo %s => %d entries: %u...%u\n",
		   nt,counts[i],positions[i][0],positions[i][counts[i]-1]);
	  }
	  FREE(nt);
	}
	);
#endif

  /* printf("nqueryoligos = %d\n",nqueryoligos); */

  return;
}

#if MAXHITS_GT_ALLOCSIZE
/* Run queryuc through this procedure */
static void
clear_positions (Genomicpos_T **positions, int *counts, char *queryuc_ptr, int seqlength) {
  int i = 0;
  int in_counter = 0;
  Shortoligomer_T oligo = 0U, masked;
  char *p;
#ifdef PMAP
  int indexsize = indexsize_aa;
  unsigned int aaindex = 0U;
#endif

  for (i = 0, p = queryuc_ptr; i < seqlength; i++, p++) {
    in_counter++;

#ifdef PMAP
    if ((index = aa_index_table[*p]) < 0) {
      aaindex = 0U;
      in_counter = 0;
    } else {
      aaindex = aaindex % msb;
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
      masked = oligo & mask;
#endif
      if (counts[masked] > ALLOCSIZE) {
	FREE(positions[masked]);
	counts[masked] = 0;	/* To avoid freeing list twice */
      }
      in_counter--;
    }
  }

  return;
}
#endif


#ifndef PMAP
#define POLY_A 0x0000
#define POLY_C 0x5555
#define POLY_G 0xAAAA
#define POLY_T 0xFFFF
#endif


/* Returns nqueryoligos */
void
Oligoindex_tally (T this, Sequence_T genomicuc) {

  memset((void *) this->counts,0,oligospace*sizeof(int));
  memset((void *) this->overabundant,0,oligospace*sizeof(bool));

#ifndef PMAP
  /* These values will prevent oligoindex from getting mappings later */
  this->overabundant[POLY_A & mask] = true;
  this->overabundant[POLY_C & mask] = true;
  this->overabundant[POLY_G & mask] = true;
  this->overabundant[POLY_T & mask] = true;
#endif

  store_positions(this->positions,this->data,this->overabundant,this->inquery,
		  this->counts,Sequence_trimpointer(genomicuc),
		  Sequence_trimlength(genomicuc));
  return;
}
  
void
Oligoindex_cleanup (T this, Sequence_T queryuc) {

#if MAXHITS_GT_ALLOCSIZE
  clear_positions(this->positions,this->counts,Sequence_trimpointer(queryuc),Sequence_trimlength(queryuc));
#endif

  return;
}


void
Oligoindex_free (T *old) {
  if (*old) {
    FREE((*old)->data);
    FREE((*old)->positions);
    FREE((*old)->counts);
    FREE((*old)->inquery);
    FREE((*old)->overabundant);
    FREE(*old);
  }
  return;
}

static Genomicpos_T *
lookup (int *nhits, T this, Shortoligomer_T masked) {
  if ((*nhits = this->counts[masked]) >= 1) {
    debug(printf("masked %04X => %d entries: %u...%u\n",
		 masked,*nhits,this->positions[masked][0],this->positions[masked][*nhits-1]));
    return this->positions[masked];
  } else {
    debug(printf("masked %04X not found\n",masked));
    /* Warning: *nhits might be -1 here, but shouldn't affect anything */
    return NULL;
  }
}


static bool
consecutivep (unsigned int *prev_mappings, int prev_nhits, unsigned int *cur_mappings, int cur_nhits) {
  int i, j;

  if (prev_nhits > 0 && cur_nhits > 0) {
    j = i = 0;
    while (j < prev_nhits && i < cur_nhits) {
      /* printf("Comparing %u with %u + %d\n",cur_mappings[i],prev_mappings[j],NT_PER_MATCH); */
      if (cur_mappings[i] == prev_mappings[j] + NT_PER_MATCH) {
	/* printf("true\n"); */
	return true;
      } else if (cur_mappings[i] < prev_mappings[j] + NT_PER_MATCH) {
	i++;
      } else {
	j++;
      }
    }
  }
  /* printf("false\n"); */
  return false;
}


/* Retrieves appropriate oligo information for a given position and
   copies it to that position */
/* Note: Be careful on totalpositions, because nhits may be < 0 */
unsigned int **
Oligoindex_get_mappings (int *nconsecutive, int **npositions, int *totalpositions, T this, 
			 Sequence_T queryuc) {
  unsigned int **mappings;
  int prev_nhits = 0, nhits, i;
  int querylength, trimlength, trim_start;

  char *p;
  int in_counter = 0, prev_querypos = 0, querypos;
  Shortoligomer_T oligo = 0U, masked;

#ifdef PMAP
  int indexsize = indexsize_aa, index;
  unsigned int aaindex = 0U;
#endif

  querylength = Sequence_fulllength(queryuc);
  trimlength = Sequence_trimlength(queryuc);
  trim_start = Sequence_trim_start(queryuc);
  mappings = (unsigned int **) CALLOC(querylength,sizeof(unsigned int *));
  *npositions = (int *) CALLOC(querylength,sizeof(int));
  *totalpositions = 0;

  *nconsecutive = 0;
  querypos = trim_start - indexsize;
  for (i = 0, p = Sequence_trimpointer(queryuc); i < trimlength; i++, p++) {
    in_counter++;
    querypos++;
    
#ifdef PMAP
    if ((index = aa_index_table[*p]) < 0) {
      aaindex = 0U;
      in_counter = 0;
    } else {
      aaindex = aaindex % msb;
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
      masked = oligo & mask;
#endif
      mappings[querypos] = lookup(&nhits,this,masked);
      (*npositions)[querypos] = nhits;
      if (nhits > 0) {
	*totalpositions += nhits;
      }
      in_counter--;

      if (consecutivep(mappings[prev_querypos],prev_nhits,mappings[querypos],nhits) == true) {
	(*nconsecutive)++;
      }
      prev_nhits = nhits;
      prev_querypos = querypos;
    }
  }

  return mappings;
}

