static char rcsid[] = "$Id: stage2.c,v 1.186 2006/11/28 00:47:13 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mem.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "listdef.h"

/* #define SQUARE 1 */

#define SUFFMAPFRACTION 0.85
#define MINMAPFRACTION 0.50

#define MIN_CONSECUTIVE 3	/* Used for trimming off ends */


/* Penalty for genomic distances */

#define INTRON_PENALTY_UNKNOWN 8
#define INTRON_PENALTY_INCONSISTENT 16

#define NINTRON_PENALTY 8
#define NINTRON_PENALTY_MISMATCH 8

#define ENOUGH_CONSECUTIVE 32

#define INFINITE 1000000

#ifdef PMAP
#define EQUAL_DISTANCE 7
#else
#define EQUAL_DISTANCE 4
#endif
#define INTRON_DEFN 9		/* Cannot exceed 9 */

#define DEAD_PENALTY 3

#ifdef PMAP
#define SAMPLE_INTERVAL 1
#define NT_PER_MATCH 3
#define CONSEC_POINTS_PER_MATCH 3 /* Possible increase to reward consecutiveness */
#else
#define SAMPLE_INTERVAL 2	/* For cases where adjacentp == false.
				   Means that we can find islands of
				   9-mers */
#define NT_PER_MATCH 1
#define CONSEC_POINTS_PER_MATCH 1 /* Possible increase to reward consecutiveness */
#endif

/* Dynamic programming */
#ifdef DEBUG
#define debug(x) x
#else 
#define debug(x)
#endif

/* Final results of stage 2 */
#ifdef DEBUG0
#define debug0(x) x
#else 
#define debug0(x)
#endif

/* Print all links */
#ifdef DEBUG1
#define debug1(x) x
#else 
#define debug1(x)
#endif

/* Shifted canonical */
#ifdef DEBUG2
#define debug2(x) x
#else 
#define debug2(x)
#endif

/* For generating a graph */
#ifdef DEBUG3
#define debug3(x) x
#else 
#define debug3(x)
#endif

/* For determining indexsize */
#ifdef DEBUG4
#define debug4(x) x
#else 
#define debug4(x)
#endif

/* Converting to nucleotides */
#ifdef DEBUG5
#define debug5(x) x
#else 
#define debug5(x)
#endif


#define T Stage2_T
struct T {
  List_T path;
  int nfwdintrons;
  int ncanonical;
  int nnoncanonical;
  int cdna_direction;

  double runtime;
  int indexsize;
  double mapfraction;
  int maxconsecutive;		/* Based on oligoindex */
  double defectrate;
};

List_T
Stage2_path (T this) {
  return this->path;
}

int
Stage2_ncanonical (T this) {
  return this->ncanonical;
}

int
Stage2_nnoncanonical (T this) {
  return this->nnoncanonical;
}

int
Stage2_cdna_direction (T this) {
  return this->cdna_direction;
}

double
Stage2_runtime (T this) {
  return this->runtime;
}

int
Stage2_indexsize (T this) {
  return this->indexsize;
}

double
Stage2_mapfraction (T this) {
  return this->mapfraction;
}

int
Stage2_maxconsecutive (T this) {
  return this->maxconsecutive;
}

double
Stage2_defectrate (T this) {
  return this->defectrate;
}


int
Stage2_pathlength (T this) {
  return List_length(this->path);
}

static T
Stage2_new (List_T path, int nfwdintrons, int nrevintrons, int nunkintrons,
	    int cdna_direction, int indexsize_nt, double mapfraction,
	    int maxconsecutive, double defectrate) {
  T new = (T) MALLOC(sizeof(*new));

  new->path = path;
  new->defectrate = defectrate;
  if (cdna_direction > 0) {
    new->ncanonical = nfwdintrons;
    new->nnoncanonical = nrevintrons + nunkintrons;
  } else if (cdna_direction < 0) {
    new->ncanonical = nrevintrons;
    new->nnoncanonical = nfwdintrons + nunkintrons;
  } else {
    abort();
  }
  new->cdna_direction = cdna_direction;
  new->indexsize = indexsize_nt;
  new->maxconsecutive = maxconsecutive;
  new->mapfraction = mapfraction;

  return new;
}
	    
void
Stage2_free (T *old) {
  if (*old) {
    /* These cells are taken from Pairpool_T and should not be freed. */
    /* List_free(&(*old)->path); */
    FREE(*old);
  }
  return;
}


/************************************************************************/

typedef struct Link_T *Link_T;
struct Link_T {
  int fwd_consecutive;
  int fwd_score;
  int fwd_pos;
  int fwd_hit;
  int fwd_intronnfwd;
  int fwd_intronnrev;
  int fwd_intronnunk;
#ifndef PMAP
  int rev_consecutive;
  int rev_score;
  int rev_pos;
  int rev_hit;
  int rev_intronnfwd;
  int rev_intronnrev;
  int rev_intronnunk;
#endif
};

/* lengths2 is has length1 entries.  Note that lengths2 may have
   negative entries */
static struct Link_T **
Linkmatrix_new (int length1, int *lengths2, int totallength) {
  struct Link_T **links;
  int i;
  
  links = (struct Link_T **) CALLOC(length1,sizeof(struct Link_T *));
  links[0] = (struct Link_T *) CALLOC(totallength,sizeof(struct Link_T));
  for (i = 1; i < length1; i++) {
    if (lengths2[i-1] < 0) {
      links[i] = links[i-1];
    } else {
      links[i] = &(links[i-1][lengths2[i-1]]);
    }
  }
  return links;
}

static void
Linkmatrix_free (struct Link_T ***links) {
  FREE((*links)[0]);
  FREE(*links);
  return;
}


/* For PMAP, indexsize is in aa */
static void
Linkmatrix_print_fwd (struct Link_T **links, unsigned int **mappings, int length1, int *npositions,
		      char *queryseq_ptr, int indexsize) {
  int i, j, lastpos;
  char *oligo;

  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
  lastpos = length1 - indexsize;

  for (i = 0; i < lastpos; i++) {
    strncpy(oligo,&(queryseq_ptr[i]),indexsize);
    printf("Querypos %d (%s, %d positions):",i,oligo,npositions[i]);
    for (j = 0; j < npositions[i]; j++) {
      printf(" %u:%d(%d,%d)",mappings[i][j],links[i][j].fwd_score,links[i][j].fwd_pos,links[i][j].fwd_hit);
    }
    printf("\n");
  }
  printf("\n");

  FREE(oligo);
  return;
}

#ifndef PMAP
static void
Linkmatrix_print_rev (struct Link_T **links, unsigned int **mappings, int length1, int *npositions,
		      char *queryseq_ptr, int indexsize) {
  int i, j;
  char *oligo;

  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
  for (i = 0; i < length1-indexsize; i++) {
    strncpy(oligo,&(queryseq_ptr[i]),indexsize);
    printf("Querypos %d (%s, %d positions):",i,oligo,npositions[i]);
    for (j = 0; j < npositions[i]; j++) {
      printf(" %u:%d(%d,%d)",mappings[i][j],links[i][j].rev_score,links[i][j].rev_pos,links[i][j].rev_hit);
    }
    printf("\n");
  }
  printf("\n");

  FREE(oligo);
  return;
}

static void
Linkmatrix_print_both (struct Link_T **links, unsigned int **mappings, int length1, int *npositions,
		      char *queryseq_ptr, int indexsize) {
  int i, j;
  char *oligo;

  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
  for (i = 0; i < length1-indexsize; i++) {
    strncpy(oligo,&(queryseq_ptr[i]),indexsize);
    printf("Querypos %d (%s, %d positions):",i,oligo,npositions[i]);
    for (j = 0; j < npositions[i]; j++) {
      printf(" %u:%d(%d,%d)-%d(%d,%d)",
	     mappings[i][j],links[i][j].fwd_score,
	     links[i][j].fwd_pos,links[i][j].fwd_hit,
	     links[i][j].rev_score,
	     links[i][j].rev_pos,links[i][j].rev_hit);
    }
    printf("\n");
  }
  printf("\n");

  FREE(oligo);
  return;
}
#endif
    


/* The value of this procedure determines the number of matches worth
   having after a gap.  For some reason, indexsize_nt - 4 finds exons of
   size indexsize_nt, but indexsize_nt - 3 does not. */
#define GAPEQUALIZER 4

#define TEN_THOUSAND 10000.0
#define HUNDRED_THOUSAND 100000.0

static int
gendist_penalty (unsigned int gendistance, bool crossspeciesp) {
  if (gendistance < TEN_THOUSAND) {
    return 0;
  } else {
    return (int) 10.0*log10((double) gendistance/TEN_THOUSAND);
  }
}


static int
gendist_penalty_dead (unsigned int gendistance, bool crossspeciesp) {
  if (gendistance < HUNDRED_THOUSAND) {
    return 0;
  } else {
    return (int) 10.0*log10((double) gendistance/TEN_THOUSAND);
  }
}


static int
querydist_mismatch_penalty (int querydistance, bool crossspeciesp) {
  if (crossspeciesp == true) {
    if (querydistance <= 24) {
      return (querydistance+7)/4;
    } else if (querydistance <= 40) {
      return (querydistance+7)/2;
    } else if (querydistance <= 56) {
      return (querydistance+7);
    } else {
      return (querydistance+7)*3/2;
    }
  } else {
    if (querydistance <= 24) {
      return (querydistance+7)/8;
    } else if (querydistance <= 40) {
      return (querydistance+7)/4;
    } else if (querydistance <= 56) {
      return (querydistance+7)/2;
    } else {
      return (querydistance+7)*3/4;
    }
  }
}

static int
diffdist_mismatch_penalty (int diffdistance, bool crossspeciesp) {
  if (crossspeciesp == true) {
    if (diffdistance <= 24) {
      return (diffdistance+7)/4;
    } else if (diffdistance <= 40) {
      return (diffdistance+7)/2;
    } else if (diffdistance <= 56) {
      return (diffdistance+7);
    } else {
      return (diffdistance+7)*3/2;
    }
  } else {
    if (diffdistance <= 24) {
      return (diffdistance+7)/8;
    } else if (diffdistance <= 40) {
      return (diffdistance+7)/4;
    } else if (diffdistance <= 56) {
      return (diffdistance+7)/2;
    } else {
      return (diffdistance+7)*3/4;
    }
  }
}


/************************************************************************
 *  Procedures for finding canonical introns quickly
 ************************************************************************/

/* A specialized Boyer-Moore algorithm */
static void
find_canonical_dinucleotides (Gbuffer_T gbuffer, char *genomicuc_ptr, int genomiclength) {
  int *lastGT, *lastAG, *lastCT, *lastAC;
  int pos, GT = -1, AG = -1, CT = -1, AC = -1;
  char c1, c2;

  lastGT = Gbuffer_lastGT(gbuffer);
  lastAG = Gbuffer_lastAG(gbuffer);
  lastCT = Gbuffer_lastCT(gbuffer);
  lastAC = Gbuffer_lastAC(gbuffer);

  lastAG[0] = lastAG[1] = lastAG[2] = -1;
  lastAC[0] = lastAC[1] = lastAC[2] = -1;

  pos = 0;
  c1 = genomicuc_ptr[pos+1];
  while (pos <= genomiclength-4) {
    if ((c2 = genomicuc_ptr[pos+2]) == 'T') {
      if (c1 == 'C') {
	CT = pos;
      } else if (c1 == 'G') {
	GT = pos;
      }
      
      lastGT[pos] = lastGT[pos+1] = GT;
      lastAG[pos+3] = lastAG[pos+4] = AG;
      lastCT[pos] = lastCT[pos+1] = CT;
      lastAC[pos+3] = lastAC[pos+4] = AC;

      pos += 2;
      c1 = genomicuc_ptr[pos+1];

    } else {
      if (c2 == 'C') {
	if (c1 == 'A') {
	  AC = pos + 3;
	}
      } else if (c2 == 'G') {
	if (c1 == 'A') {
	  AG = pos + 3;
	}
      }
	
      lastGT[pos] = GT;
      lastAG[pos+3] = AG;
      lastCT[pos] = CT;
      lastAC[pos+3] = AC;

      pos++;
      c1 = c2;
    }


  }

  /* pos is now either genomiclength-3 (one more dinucleotide) or genomiclength-2 (none) */
  if (pos == genomiclength-3) {
    if ((c2 = genomicuc_ptr[pos+2]) == 'T') {
      if (c1 == 'C') {
	CT = pos;
      } else if (c1 == 'G') {
	GT = pos;
      }
    } else {
      if (c2 == 'C') {
	if (c1 == 'A') {
	  AC = pos + 3;
	}
      } else if (c2 == 'G') {
	if (c1 == 'A') {
	  AG = pos + 3;
	}
      }
    }
    
    lastGT[pos] = GT;
    lastAG[pos+3] = AG;
    lastCT[pos] = CT;
    lastAC[pos+3] = AC;
  }

  return;
}


static bool
find_shifted_canonical (char *genomicuc_ptr, unsigned int prevposition, unsigned int position,
			int querydistance, int indexsize_nt, int *last_leftdi, int *last_rightdi) {
  unsigned int leftpos, rightpos;
  int shift, leftmiss, rightmiss;

  leftpos = prevposition + querydistance - 1;
  rightpos = position;
  shift = 0;

  debug2(printf("Looking for shifted canonical at prevposition %d to position %d\n",prevposition,position));
  while (shift <= querydistance - indexsize_nt) {
    leftmiss = leftpos - last_leftdi[leftpos];
    rightmiss = rightpos - last_rightdi[rightpos];
    debug2(printf("%d/%d/%d",shift,leftmiss,rightmiss));
    if (leftmiss == 0 && rightmiss == 0) {
      debug2(printf(" => Success\n\n"));
      return true;
    } else if (leftmiss >= rightmiss) {
      shift += leftmiss;
      leftpos -= leftmiss;
      rightpos -= leftmiss;
    } else {
      shift += rightmiss;
      leftpos -= rightmiss;
      rightpos -= rightmiss;
    }
    debug2(printf("\n"));
  }
  debug2(printf("\n"));
  
  return false;
}


/* For PMAP, indexsize is in aa. */
static void
score_querypos (Link_T currlink, int querypos, unsigned int position,
		struct Link_T **links, unsigned int **mappings, int *npositions, 
		int *fwd_restrict_hit, 
#ifndef PMAP
		int *rev_restrict_hit,
#endif
		char *genomicuc_ptr, int *lastGT, int *lastAG, int *lastCT, int *lastAC, 
		int indexsize, int sufflookback, int nsufflookback, int maxintronlen,
		bool deadp, bool crossspeciesp, double mapfraction) {
  Link_T prevlink;
  int best_fwd_consecutive = 0;
  int best_fwd_score = 0, fwd_score;
  int best_fwd_prevpos = -1, best_fwd_prevhit = -1;
  int best_fwd_intronnfwd = 0, best_fwd_intronnrev = 0, best_fwd_intronnunk = 0;
#ifndef PMAP
  int best_rev_consecutive = 0;
  int best_rev_score = 0, rev_score;
  int best_rev_prevpos = -1, best_rev_prevhit = -1;
  int best_rev_intronnfwd = 0, best_rev_intronnrev = 0, best_rev_intronnunk = 0;
#endif
  bool adjacentp = false;
  int prev_querypos, intronpos_lower, intronpos_upper, prevhit, lastprevhit;
  unsigned int *prev_mappings;
  unsigned int prevposition, gendistance;
  int querydistance, diffdistance, lookback, nlookback, sample_interval, nseen, indexsize_nt;
  int canonicalsgn, mismatch_start, fwd_gendistance_penalty;
#ifndef PMAP
  int rev_gendistance_penalty;
#endif
  int enough_consecutive;

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

  if (crossspeciesp == true) {
    enough_consecutive = 56;
  } else {
    enough_consecutive = 32;
  }

  /* A. Evaluate adjacent position (at querypos - 1) */
  for (prevhit = 0; prevhit < npositions[querypos-1]; prevhit++) {
    if (mappings[querypos-1][prevhit] == position - NT_PER_MATCH) {
      prevlink = &(links[querypos-1][prevhit]);
      best_fwd_consecutive = prevlink->fwd_consecutive + NT_PER_MATCH;
      best_fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_MATCH;
#ifndef PMAP
      best_rev_consecutive = prevlink->rev_consecutive + NT_PER_MATCH;
      best_rev_score = prevlink->rev_score + CONSEC_POINTS_PER_MATCH;
#endif

      best_fwd_prevpos = 
#ifndef PMAP
	best_rev_prevpos = 
#endif
	querypos-1;
      best_fwd_prevhit = 
#ifndef PMAP
	best_rev_prevhit = 
#endif
	prevhit;
      best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
      best_fwd_intronnrev = prevlink->fwd_intronnrev;
      best_fwd_intronnunk = prevlink->fwd_intronnunk;
#ifndef PMAP
      best_rev_intronnfwd = prevlink->rev_intronnfwd;
      best_rev_intronnrev = prevlink->rev_intronnrev;
      best_rev_intronnunk = prevlink->rev_intronnunk;
#endif
      adjacentp = true;
#ifdef PMAP
      debug(printf("\tA. Adjacent qpos %d,%d at %u (scores = %d -> %d, consec = %d, intr = %d-%d-%d)\n",
		   querypos-1,prevhit,position - NT_PER_MATCH,prevlink->fwd_score,
		   best_fwd_score,best_fwd_consecutive,best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk));
#else
      debug(printf("\tA. Adjacent qpos %d,%d at %u (scores = %d/%d -> %d/%d, consec = %d/%d, intr = %d-%d-%d/%d-%d-%d)\n",
		   querypos-1,prevhit,position - NT_PER_MATCH,prevlink->fwd_score,prevlink->rev_score,
		   best_fwd_score,best_rev_score,best_fwd_consecutive,best_rev_consecutive,
		   best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		   best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk));
#endif
    }
  }

  if (best_fwd_consecutive < enough_consecutive
#ifndef PMAP
      || best_rev_consecutive < enough_consecutive
#endif
      ) {
    /* B. Evaluate for intron (at querypos - indexsize).  
       Don't update best_consecutive, because a reasonable alternative 
       might be a mismatch, regardless of the quality of the intron. */

#ifdef PMAP
    intronpos_lower = querypos - indexsize - 1;
    intronpos_upper = querypos - indexsize;
#else
    intronpos_lower = intronpos_upper = querypos - indexsize;
#endif
    for (prev_querypos = intronpos_upper; prev_querypos >= intronpos_lower; --prev_querypos) {
      if (prev_querypos >= 0 && npositions[prev_querypos] > 0) {
#ifdef PMAP
	querydistance = (querypos - prev_querypos)*3;
#else
	querydistance = querypos - prev_querypos;
#endif
	prev_mappings = &(mappings[prev_querypos][0]);
	for (prevhit = 0; prevhit < npositions[prev_querypos]; prevhit++) {
	  prevposition = *(prev_mappings++);
	  if (position >= prevposition + indexsize_nt) {
	    prevlink = &(links[prev_querypos][prevhit]);
	    gendistance = position - prevposition;
	    diffdistance = abs(gendistance - querydistance);

	    if (diffdistance > maxintronlen) {
	      canonicalsgn = 0;
	    } else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,
					      indexsize_nt,lastGT,lastAG) == true) {
	      canonicalsgn = +1;
	    } else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,
					      indexsize_nt,lastCT,lastAC) == true) {
	      canonicalsgn = -1;
	    } else {
	      canonicalsgn = 0;
	    }

	    /* Need to both compensate for the oligomer lost across
	       the intron and penalize for each additional intron */
#ifdef PMAP
	    fwd_score = prevlink->fwd_score + querydistance - NINTRON_PENALTY;
#else
	    fwd_score = prevlink->fwd_score + querydistance - NINTRON_PENALTY;
	    rev_score = prevlink->rev_score + querydistance - NINTRON_PENALTY;
#endif

	    /* Penalty for introns. */
	    if (diffdistance > maxintronlen) {
	      fwd_gendistance_penalty = 
#ifndef PMAP
		rev_gendistance_penalty = 
#endif
		INFINITE;
	    } else if (deadp == true) {
	      fwd_gendistance_penalty = 
#ifndef PMAP
		rev_gendistance_penalty = 
#endif
		gendist_penalty_dead(gendistance,crossspeciesp) + INTRON_PENALTY_UNKNOWN + DEAD_PENALTY;
	    } else {
	      if (canonicalsgn == 0) {
	      /* Extra penalty for known non-canonical intron */
		fwd_gendistance_penalty = 
#ifndef PMAP
		  rev_gendistance_penalty = 
#endif
		  gendist_penalty(gendistance,crossspeciesp) + INTRON_PENALTY_UNKNOWN;
	      } else if (canonicalsgn == +1) {
		fwd_gendistance_penalty = gendist_penalty(gendistance,crossspeciesp) /*- INTRON_REWARD_CONSISTENT*/;
#ifndef PMAP
		rev_gendistance_penalty = gendist_penalty(gendistance,crossspeciesp) + INTRON_PENALTY_INCONSISTENT;
#endif
	      } else if (canonicalsgn == -1) {
#ifdef PMAP
		fwd_gendistance_penalty = gendist_penalty(gendistance,crossspeciesp) + INTRON_PENALTY_UNKNOWN;
#else
		fwd_gendistance_penalty = gendist_penalty(gendistance,crossspeciesp) + INTRON_PENALTY_INCONSISTENT;
		rev_gendistance_penalty = gendist_penalty(gendistance,crossspeciesp) /*- INTRON_REWARD_CONSISTENT*/;
#endif
	      }
	    }
	    fwd_score -= fwd_gendistance_penalty;
#ifndef PMAP
	    rev_score -= rev_gendistance_penalty;
#endif

#ifdef PMAP
	    debug(printf("\tB. Intron qpos %d,%d at %u (scores = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, pen %d)",
			 prev_querypos,prevhit,prevposition,
			 prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,
			 best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			 gendistance,fwd_gendistance_penalty));
#else
	    debug(printf("\tB. Intron qpos %d,%d at %u (scores = %d/%d -> %d/%d, consec = %d/%d, intr = %d-%d-%d/%d-%d-%d, gendist %u, pen %d/%d)",
			 prev_querypos,prevhit,prevposition,
			 prevlink->fwd_score,prevlink->rev_score,fwd_score,rev_score,
			 prevlink->fwd_consecutive,prevlink->rev_consecutive,
			 best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			 best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk,
			 gendistance,fwd_gendistance_penalty,rev_gendistance_penalty));
#endif
	  
	    debug(printf(" =>"));
	    if (fwd_score > best_fwd_score) {
#ifdef PMAP
	      best_fwd_consecutive = indexsize*3;
#else
	      best_fwd_consecutive = indexsize;
#endif
	      best_fwd_score = fwd_score;
	      best_fwd_prevpos = prev_querypos;
	      best_fwd_prevhit = prevhit;
	      best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	      best_fwd_intronnrev = prevlink->fwd_intronnrev;
	      best_fwd_intronnunk = prevlink->fwd_intronnunk;
	      switch (canonicalsgn) {
	      case 1: best_fwd_intronnfwd++; break;
	      case -1: best_fwd_intronnrev++; break;
	      case 0: best_fwd_intronnunk++; break;
	      }
	      debug(printf(" Best fwd at %d",fwd_score));
	    }

#ifndef PMAP
	    if (rev_score > best_rev_score) {
	      best_rev_consecutive = indexsize;
	      best_rev_score = rev_score;
	      best_rev_prevpos = prev_querypos;
	      best_rev_prevhit = prevhit;
	      best_rev_intronnfwd = prevlink->rev_intronnfwd;
	      best_rev_intronnrev = prevlink->rev_intronnrev;
	      best_rev_intronnunk = prevlink->rev_intronnunk;
	      switch (canonicalsgn) {
	      case 1: best_rev_intronnfwd++; break;
	      case -1: best_rev_intronnrev++; break;
	      case 0: best_rev_intronnunk++; break;
	      }
	      debug(printf(" Best rev at %d",rev_score));
	    }
#endif
	    debug(printf("\n"));
	  }
	}
      }
    }
    debug(printf("\n"));

    /* Use intronpos_upper instead of intronpos_lower because correct position is uncertain in PMAP */
    mismatch_start = intronpos_upper - 1;

    /* C (fwd). Evaluate for mismatches (all other previous querypos) */
    /* Set parameters */
    if (adjacentp == true) {
      /* Sample at larger interval and not so far back */
      nlookback = 1;
      lookback = sufflookback/2;
      sample_interval = indexsize;
    } else if (deadp == true) {
      /* Sample at smaller interval and farther back */
      nlookback = 2*nsufflookback;
      lookback = 2*sufflookback;
      sample_interval = SAMPLE_INTERVAL;
    } else {
      /* Sample at smaller interval and farther back */
      nlookback = nsufflookback;
      lookback = sufflookback;
      sample_interval = SAMPLE_INTERVAL;
    }

    /*printf("        mismatch_start = %d, nlookback = %d, lookback = %d, sample_interval = %d\n",
      mismatch_start,nlookback,lookback,sample_interval); */

    /* The nseen check protects against occurrences of N's in
       the sequences, by making sure we see at least some hits */
    nseen = 0;
    for (prev_querypos = mismatch_start;
	 prev_querypos >= 0 && best_fwd_consecutive < enough_consecutive && 
	   (nseen < nlookback || (querypos - prev_querypos) < lookback);
	 prev_querypos -= sample_interval) {
      if (npositions[prev_querypos] > 0) {
	nseen++;
#ifdef PMAP
	querydistance = (querypos - prev_querypos)*3;
#else
	querydistance = querypos - prev_querypos;
#endif
	if ((prevhit = fwd_restrict_hit[prev_querypos]) < 0) {
	  prevhit = 0;	/* Not restricted */
	  lastprevhit = npositions[prev_querypos] - 1;
	} else {
	  lastprevhit = prevhit; /* Restricts loop to one iteration */
	}
	prev_mappings = &(mappings[prev_querypos][prevhit]);
	for ( ; prevhit <= lastprevhit; prevhit++) {
	  prevposition = *(prev_mappings++);
	  if (position >= prevposition + indexsize_nt) {
	    prevlink = &(links[prev_querypos][prevhit]);

	    gendistance = position - prevposition;
	    diffdistance = abs(gendistance - querydistance);
	    fwd_score = prevlink->fwd_score;
	    if (diffdistance > maxintronlen) {
	      fwd_score -= INFINITE;
	    } else if (querydistance < INTRON_DEFN) {
	      fwd_score += querydistance;
	    } else if (diffdistance < EQUAL_DISTANCE) {
	      fwd_score += querydistance - diffdist_mismatch_penalty(diffdistance,crossspeciesp);
	      fwd_score -= querydist_mismatch_penalty(querydistance,crossspeciesp);
#ifdef PMAP
	      if (diffdistance % 3 != 0) {
		/* fwd_score -= non_codon_penalty; */
	      }
#endif
	    } else if (diffdistance < INTRON_DEFN) {
	      fwd_score -= querydist_mismatch_penalty(querydistance,crossspeciesp);
	    } else {
	      fwd_score -= querydist_mismatch_penalty(querydistance,crossspeciesp) + NINTRON_PENALTY_MISMATCH;
	    }

	    /* For deadp use diffdistance, not gendistance */
	    if (deadp == true) {
	      fwd_score -= gendist_penalty_dead(diffdistance,crossspeciesp) + DEAD_PENALTY;
	    } else {
	      fwd_score -= gendist_penalty(diffdistance,crossspeciesp);
	    }

	    debug(printf("\tC+. Fwd mismatch qpos %d,%d at %u (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d)",
			 prev_querypos,prevhit,prevposition,prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,
			 best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			 gendistance,querydistance));
	    
	    if (fwd_score > best_fwd_score) {
#ifdef PMAP
	      best_fwd_consecutive = indexsize*3;
#else
	      best_fwd_consecutive = indexsize;
#endif
	      best_fwd_score = fwd_score;
	      best_fwd_prevpos = prev_querypos;
	      best_fwd_prevhit = prevhit;
	      best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	      best_fwd_intronnrev = prevlink->fwd_intronnrev;
	      best_fwd_intronnunk = prevlink->fwd_intronnunk;
	      if (diffdistance > INTRON_DEFN) {
		best_fwd_intronnunk++;
	      }
#ifdef PMAP
	      /* sample_interval = indexsize; */
#else
	      if (crossspeciesp == false) {
		sample_interval = indexsize;
	      }
#endif
	      debug(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	    } else {
	      debug(printf(" => Loses to %d\n",best_fwd_score));
	    }
	  }
	}
      }
    }

#ifndef PMAP
    /* C (rev). Evaluate for mismatches (all other previous querypos) */
    /* Set parameters */
    if (adjacentp == true) {
      /* Sample at larger interval and not so far back */
      nlookback = 1;
      lookback = sufflookback/2;
      sample_interval = indexsize;
    } else if (deadp == true) {
      /* Sample at smaller interval and farther back */
      nlookback = 2*nsufflookback;
      lookback = 2*sufflookback;
      sample_interval = SAMPLE_INTERVAL;
    } else {
      /* Sample at smaller interval and farther back */
      nlookback = nsufflookback;
      lookback = sufflookback;
      sample_interval = SAMPLE_INTERVAL;
    }

    /* The nseen check protects against occurrences of N's in
       the sequences, by making sure we see at least some hits */
    nseen = 0;
    for (prev_querypos = mismatch_start;
	 prev_querypos >= 0 && best_rev_consecutive < enough_consecutive && 
	   (nseen < nlookback || (querypos - prev_querypos) < lookback);
	 prev_querypos -= sample_interval) {
      if (npositions[prev_querypos] > 0) {
	nseen++;
	querydistance = querypos - prev_querypos;
	if ((prevhit = rev_restrict_hit[prev_querypos]) < 0) {
	  prevhit = 0;	/* Not restricted */
	  lastprevhit = npositions[prev_querypos] - 1;
	} else {
	  lastprevhit = prevhit; /* Restricts loop to one iteration */
	}
	prev_mappings = &(mappings[prev_querypos][prevhit]);
	for ( ; prevhit <= lastprevhit; prevhit++) {
	  prevposition = *(prev_mappings++);
	  if (position >= prevposition + indexsize_nt) {
	    prevlink = &(links[prev_querypos][prevhit]);
	    
	    gendistance = position - prevposition;
	    diffdistance = abs(gendistance - querydistance);
	    rev_score = prevlink->rev_score;
	    if (diffdistance > maxintronlen) {
	      rev_score -= INFINITE;
	    } else if (querydistance < INTRON_DEFN) {
	      rev_score += querydistance;
	    } else if (diffdistance < EQUAL_DISTANCE) {
	      rev_score += querydistance - diffdist_mismatch_penalty(diffdistance,crossspeciesp);
	      rev_score -= querydist_mismatch_penalty(querydistance,crossspeciesp);
	    } else if (diffdistance < INTRON_DEFN) {
	      rev_score -= querydist_mismatch_penalty(querydistance,crossspeciesp);
	    } else {
	      rev_score -= querydist_mismatch_penalty(querydistance,crossspeciesp) + NINTRON_PENALTY_MISMATCH;
	    }

	    /* For deadp, use diffdistance, not gendistance */
	    if (deadp == true) {
	      rev_score -= gendist_penalty_dead(diffdistance,crossspeciesp) + DEAD_PENALTY;
	    } else {
	      rev_score -= gendist_penalty(diffdistance,crossspeciesp);
	    }
	    
	    debug(printf("\tC-. Rev mismatch qpos %d,%d at %u (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d)",
			 prev_querypos,prevhit,prevposition,prevlink->rev_score,rev_score,prevlink->rev_consecutive,
			 best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk,
			 gendistance,querydistance));
	    
	    if (rev_score > best_rev_score) {
	      best_rev_consecutive = indexsize;
	      best_rev_score = rev_score;
	      best_rev_prevpos = prev_querypos;
	      best_rev_prevhit = prevhit;
	      best_rev_intronnfwd = prevlink->rev_intronnfwd;
	      best_rev_intronnrev = prevlink->rev_intronnrev;
	      best_rev_intronnunk = prevlink->rev_intronnunk;
	      if (diffdistance > INTRON_DEFN) {
		best_rev_intronnunk++;
	      }
	      if (crossspeciesp == false) {
		sample_interval = indexsize;
	      }
	      debug(printf(" => Best rev at %d (consec = %d)\n",rev_score,best_rev_consecutive));
	    } else {
	      debug(printf(" => Loses to %d\n",best_rev_score));
	    }
	  }
	}
      }
    }
#endif
  }

  /* Best_score needs to beat something positive to prevent a
     small local extension from beating a good canonical intron.
     If querypos is too small, don't insert an intron.  */
  /* linksconsecutive already assigned above */
  currlink->fwd_consecutive = best_fwd_consecutive;
  currlink->fwd_pos = best_fwd_prevpos;
  currlink->fwd_hit = best_fwd_prevhit;
  currlink->fwd_score = best_fwd_score;
  currlink->fwd_intronnfwd = best_fwd_intronnfwd;
  currlink->fwd_intronnrev = best_fwd_intronnrev;
  currlink->fwd_intronnunk = best_fwd_intronnunk;

#ifndef PMAP
  /* linksconsecutive already assigned above */
  currlink->rev_consecutive = best_rev_consecutive;
  currlink->rev_pos = best_rev_prevpos;
  currlink->rev_hit = best_rev_prevhit;
  currlink->rev_score = best_rev_score;
  currlink->rev_intronnfwd = best_rev_intronnfwd;
  currlink->rev_intronnrev = best_rev_intronnrev;
  currlink->rev_intronnunk = best_rev_intronnunk;
#endif

#ifdef PMAP
  debug(printf("\tChose %d,%d with score %d (fwd)\n",
	       currlink->fwd_pos,currlink->fwd_hit,currlink->fwd_score));
#else
  debug(printf("\tChose %d,%d with score %d (fwd) and %d,%d with score %d (rev)\n",
	       currlink->fwd_pos,currlink->fwd_hit,currlink->fwd_score,currlink->rev_pos,currlink->rev_hit,currlink->rev_score));
#endif
  debug3(printf("%d %d  %d %d  1\n",querypos,hit,best_prevpos,best_prevhit));
  return;
}


/* For PMAP, indexsize is in aa. */
static Link_T
align_compute_scores (bool *fwdp, struct Link_T **links,
		      unsigned int **mappings, int *npositions, int querystart, int queryend,
		      int querylength, Sequence_T genomicuc, Gbuffer_T gbuffer, int indexsize,
		      int sufflookback, int nsufflookback, int maxintronlen, bool crossspeciesp,
		      double mapfraction, char *queryseq_ptr) {
  Link_T termlink = NULL, currlink;
  int fwd_score, *fwd_restrict_hit;
#ifndef PMAP
  int rev_score, *rev_restrict_hit;
#endif
  int querypos, indexsize_nt, hit;
  int best_score, best_hit, second_score, score;
  unsigned int position;
  char *genomicuc_ptr;
  int *lastGT, *lastAG, *lastCT, *lastAC;
  debug(char *oligo);

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

  debug(oligo = (char *) CALLOC(indexsize+1,sizeof(char)));
  debug0(printf("Querystart = %d, queryend = %d\n",querystart,queryend));

  genomicuc_ptr = Sequence_fullpointer(genomicuc);
  find_canonical_dinucleotides(gbuffer,genomicuc_ptr,Sequence_fulllength(genomicuc));
  lastGT = Gbuffer_lastGT(gbuffer);
  lastAG = Gbuffer_lastAG(gbuffer);
  lastCT = Gbuffer_lastCT(gbuffer);
  lastAC = Gbuffer_lastAC(gbuffer);

  fwd_restrict_hit = (int *) CALLOC(querylength,sizeof(int));
#ifndef PMAP
  rev_restrict_hit = (int *) CALLOC(querylength,sizeof(int));
#endif

  /* Treat querypos 0 as a special case */
  fwd_restrict_hit[querystart] = 
#ifndef PMAP
    rev_restrict_hit[querystart] = 
#endif
    -1;

  for (hit = 0; hit < npositions[querystart]; hit++) {
#ifdef PMAP    
    links[querystart][hit].fwd_pos = links[querystart][hit].fwd_hit = -1;
    links[querystart][hit].fwd_score = indexsize*3;
    links[querystart][hit].fwd_consecutive = indexsize*3;
#else
    links[querystart][hit].fwd_pos = links[querystart][hit].fwd_hit = -1;
    links[querystart][hit].rev_pos = links[querystart][hit].rev_hit = -1;
    links[querystart][hit].fwd_score = indexsize;
    links[querystart][hit].rev_score = indexsize;
    links[querystart][hit].fwd_consecutive = indexsize;
    links[querystart][hit].rev_consecutive = indexsize;
#endif
  }

  for (querypos = querystart+1; querypos <= queryend - indexsize; querypos++) {
    best_score = 0;
    for (hit = 0; hit < npositions[querypos]; hit++) {
      currlink = &(links[querypos][hit]);
      position = mappings[querypos][hit];

      debug(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
      debug(printf("Finding link at querypos %d,%d at %u (%s)\n",querypos,hit,position,oligo));

      score_querypos(currlink,querypos,position,links,mappings,npositions,fwd_restrict_hit,
#ifndef PMAP
		     rev_restrict_hit,
#endif
		     genomicuc_ptr,lastGT,lastAG,lastCT,lastAC,indexsize,sufflookback,nsufflookback,
		     maxintronlen,/*deadp*/false,crossspeciesp,mapfraction);
      if (currlink->fwd_score > best_score) {
	best_score = currlink->fwd_score;
      }
#ifndef PMAP
      if (currlink->rev_score > best_score) {
	best_score = currlink->rev_score;
      }
#endif
    }
      
    if (best_score <= 0) {
      for (hit = 0; hit < npositions[querypos]; hit++) {
	currlink = &(links[querypos][hit]);
	position = mappings[querypos][hit];

	debug(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
	debug(printf("Redoing link at querypos %d,%d at %u, without gendistance penalty\n",querypos,hit,position));

	score_querypos(currlink,querypos,position,links,mappings,npositions,fwd_restrict_hit,
#ifndef PMAP
		       rev_restrict_hit,
#endif
		       genomicuc_ptr,lastGT,lastAG,lastCT,lastAC,indexsize,sufflookback,nsufflookback,
		       maxintronlen,/*deadp*/true,crossspeciesp,mapfraction);
      }
    }

    /* Compute best hit over all hits at this querypos, if clearly a best */
    best_hit = -1;
    best_score = second_score = -100000;
    for (hit = 0; hit < npositions[querypos]; hit++) {
      fwd_score = links[querypos][hit].fwd_score;
      if (fwd_score >= best_score) {
	second_score = best_score;
	best_score = fwd_score;
	best_hit = hit;
      } else if (fwd_score >= second_score) {
	second_score = fwd_score;
      }
    }
    if (best_score - second_score > indexsize_nt - GAPEQUALIZER) {
      fwd_restrict_hit[querypos] = best_hit;
    } else {
      fwd_restrict_hit[querypos] = -1;
    }

#ifndef PMAP
    /* Compute best hit over all hits at this querypos, if clearly a best */
    best_hit = -1;
    best_score = second_score = -100000;
    for (hit = 0; hit < npositions[querypos]; hit++) {
      rev_score = links[querypos][hit].rev_score;
      if (rev_score >= best_score) {
	second_score = best_score;
	best_score = rev_score;
	best_hit = hit;
      } else if (rev_score >= second_score) {
	second_score = rev_score;
      }
    }
    if (best_score - second_score > indexsize_nt - GAPEQUALIZER) {
      rev_restrict_hit[querypos] = best_hit;
    } else {
      rev_restrict_hit[querypos] = -1;
    }
#endif
  }

  /* Find best match */
  best_score = 0;
  termlink = NULL;
  /* Go backwards, because best score is likely to be there, and so
     ties favor longer path */

#ifdef PMAP
  *fwdp = true;
  for (querypos = queryend - indexsize; querypos > querystart; --querypos) {
    /* Don't use restrict_hit here, because it doesn't save time */
    for (hit = 0; hit < npositions[querypos]; hit++) {
      if ((score = links[querypos][hit].fwd_score) > best_score) {
	best_score = score;
	termlink = &(links[querypos][hit]);
      }
    }
  }
  debug(printf("==> Found best score %d at fwd %d,%d\n",best_score,termlink->fwd_pos,termlink->fwd_hit));

#else

  for (querypos = queryend - indexsize; querypos > querystart; --querypos) {
    /* Don't use restrict_hit here, because it doesn't save time */
    for (hit = 0; hit < npositions[querypos]; hit++) {
      if ((score = links[querypos][hit].fwd_score) > best_score) {
	best_score = score;
	termlink = &(links[querypos][hit]);
	*fwdp = true;
      }
      if ((score = links[querypos][hit].rev_score) > best_score) {
	best_score = score;
	termlink = &(links[querypos][hit]);
	*fwdp = false;
      }
    }
  }

  debug(
	if (*fwdp == true) {
	  printf("==> Found best score %d at fwd %d,%d\n",best_score,termlink->fwd_pos,termlink->fwd_hit);
	} else {
	  printf("==> Found best score %d at rev %d,%d\n",best_score,termlink->rev_pos,termlink->rev_hit);
	}
	);
#endif

#ifndef PMAP
  FREE(rev_restrict_hit);
#endif
  FREE(fwd_restrict_hit);

  debug(FREE(oligo));

  return termlink;
}



/* Performs dynamic programming.  For PMAP, indexsize is in aa. */
static List_T
align_compute (int *pathlength, double *defectrate,
	       int *nfwdintrons, int *nrevintrons, int *nunkintrons, int *cdna_direction,
	       unsigned int **mappings, int *npositions, int totalpositions, 
	       Sequence_T queryseq, Sequence_T genomicseg, Sequence_T genomicuc, Gbuffer_T gbuffer, int indexsize,
	       int sufflookback, int nsufflookback, int maxintronlen, Pairpool_T pairpool, bool crossspeciesp,
	       double mapfraction) {
  List_T path = NULL;
  struct Link_T **links;
  Link_T termlink;
  bool fwdp;
  int termpos, termhit;
  int querypos, prev_querypos, hit, prevhit;
  unsigned int position, prevposition;
  int querystart, queryend, querylength;
  char *queryseq_ptr, *genomicseg_ptr;
#ifdef PMAP
  char *genomicuc_ptr;
#endif
  int nconsecutive = 0, nintrons = 0, ndefects = 0;
  debug0(char *oligo);

  debug0(oligo = (char *) CALLOC(indexsize+1,sizeof(char)));

  *pathlength = 0;

  queryseq_ptr = Sequence_fullpointer(queryseq);
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
#ifdef PMAP
  genomicuc_ptr = Sequence_fullpointer(genomicuc);
#endif

  querystart = Sequence_trim_start(queryseq);
  queryend = Sequence_trim_end(queryseq);
  querylength = Sequence_fulllength(queryseq);

  links = Linkmatrix_new(querylength,npositions,totalpositions);

  termlink = align_compute_scores(&fwdp,links,mappings,npositions,querystart,queryend,querylength,
				  genomicuc,gbuffer,indexsize,sufflookback,nsufflookback,maxintronlen,crossspeciesp,
				  mapfraction,queryseq_ptr);
  debug(printf("align_compute_scores returns fwdp = %d\n",fwdp));

  if (termlink != NULL) {
#ifdef PMAP
    termpos = termlink->fwd_pos;
    termhit = termlink->fwd_hit;
    *nfwdintrons = termlink->fwd_intronnfwd;
    *nrevintrons = termlink->fwd_intronnrev;
    *nunkintrons = termlink->fwd_intronnunk;
#else
    if (fwdp == true) {
      termpos = termlink->fwd_pos;
      termhit = termlink->fwd_hit;
      *nfwdintrons = termlink->fwd_intronnfwd;
      *nrevintrons = termlink->fwd_intronnrev;
      *nunkintrons = termlink->fwd_intronnunk;
    } else {
      termpos = termlink->rev_pos;
      termhit = termlink->rev_hit;
      *nfwdintrons = termlink->rev_intronnfwd;
      *nrevintrons = termlink->rev_intronnrev;
      *nunkintrons = termlink->rev_intronnunk;
    }
#endif

#ifdef PMAP
    debug1(Linkmatrix_print_fwd(links,mappings,querylength,npositions,queryseq_ptr,indexsize));
#else
    debug1(Linkmatrix_print_both(links,mappings,querylength,npositions,queryseq_ptr,indexsize));
#endif

    querypos = termpos;
    hit = termhit;

    /* Traceback */
    /* nconsecutive = nintrons = ndefects = 0; */
    if (querypos != -1) {
      prevposition = mappings[querypos][hit]+1;
      prev_querypos = querypos+1;
    }
    while (querypos >= 0) {
      position = mappings[querypos][hit];
      debug0(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
#ifdef PMAP
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d, intr = %d(+)/%d(-)/%d(?)\n",
		    querypos,hit,oligo,position,links[querypos][hit].fwd_score,links[querypos][hit].fwd_consecutive,
		    links[querypos][hit].fwd_intronnfwd,links[querypos][hit].fwd_intronnrev,
		    links[querypos][hit].fwd_intronnunk));
#else
      debug0(
	     if (fwdp) {
	       printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d, intr = %d(+)/%d(-)/%d(?)\n",
		      querypos,hit,oligo,position,links[querypos][hit].fwd_score,links[querypos][hit].fwd_consecutive,
		      links[querypos][hit].fwd_intronnfwd,links[querypos][hit].fwd_intronnrev,
		      links[querypos][hit].fwd_intronnunk);
	     } else {
	       printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d, intr = %d(+)/%d(-)/%d(?)\n",
		      querypos,hit,oligo,position,links[querypos][hit].rev_score,links[querypos][hit].rev_consecutive,
		      links[querypos][hit].rev_intronnfwd,links[querypos][hit].rev_intronnrev,
		      links[querypos][hit].rev_intronnunk);
	     });
#endif

#ifdef PMAP
      /* Change querypos positions from protein to nucleotide */
      path = Pairpool_push(path,pairpool,querypos*3+2,position+2,
 			   genomicuc_ptr[position+2],MATCH_COMP,
 			   genomicseg_ptr[position+2],/*dynprogindex*/0);
      path = Pairpool_push(path,pairpool,querypos*3+1,position+1,
 			   genomicuc_ptr[position+1],MATCH_COMP,
 			   genomicseg_ptr[position+1],/*dynprogindex*/0);
      path = Pairpool_push(path,pairpool,querypos*3,position,
 			   genomicuc_ptr[position],MATCH_COMP,
			   genomicseg_ptr[position],/*dynprogindex*/0);
#else
      path = Pairpool_push(path,pairpool,querypos,position,
			   queryseq_ptr[querypos],MATCH_COMP,
			   genomicseg_ptr[position],/*dynprogindex*/0);
#endif
      
      /* Be careful with parallel reads of querypos and hit */
      if (position == prevposition - NT_PER_MATCH) {
	nconsecutive++;
      } else if (prevposition - position > prev_querypos - querypos) {
	nintrons++;
      } else {
	ndefects++;
      }

      prevposition = position;
      prev_querypos = querypos;
      prevhit = hit;
      if (fwdp) {
	querypos = links[prev_querypos][prevhit].fwd_pos;
	hit = links[prev_querypos][prevhit].fwd_hit;
#ifndef PMAP
      } else {
	querypos = links[prev_querypos][prevhit].rev_pos;
	hit = links[prev_querypos][prevhit].rev_hit;
#endif
      }
      debug3(printf("%d %d  %d %d  3\n",prev_querypos,prevhit,querypos,hit));
    }

    if (nconsecutive == 0) {
      *defectrate = 1.0;
    } else {
      *defectrate = (double) (ndefects)/(double) (nconsecutive + ndefects);
    }
    *pathlength = nconsecutive + nintrons + ndefects;
  }
  Linkmatrix_free(&links);

  debug0(FREE(oligo));

  if (fwdp == true) {
    *cdna_direction = +1;
  } else {
    *cdna_direction = -1;
  }

  return path;			/* previously was List_reverse(path) */
}


/* queryseq_ptr is NULL for PMAP.  querypos here is in nt. */
static List_T
convert_to_nucleotides (List_T path, char *queryseq_ptr, char *genomicseg_ptr,
#ifdef PMAP
			char *genomicuc_ptr,
#endif
			Pairpool_T pairpool, int indexsize_nt) {
  List_T pairs = NULL;
  Pair_T pair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump, fill;

  pair = path->first;
  querypos = pair->querypos;
  genomepos = pair->genomepos;

#ifdef PMAP
  fill = indexsize_nt - 3;
#else
  fill = indexsize_nt - 1;
#endif

  lastquerypos = querypos + fill;
  lastgenomepos = genomepos + fill;
  while (lastquerypos > querypos) {
#ifdef PMAP
    pairs = Pairpool_push(pairs,pairpool,lastquerypos,lastgenomepos,
			  genomicuc_ptr[lastgenomepos],MATCH_COMP,
			  genomicseg_ptr[lastgenomepos],/*dynprogindex*/0);
    debug5(printf("Pushing %c | %c at %d,%d\n",genomicuc_ptr[lastgenomepos],genomicseg_ptr[lastgenomepos],
		  lastquerypos,lastgenomepos));
#else
    pairs = Pairpool_push(pairs,pairpool,lastquerypos,lastgenomepos,
			  queryseq_ptr[lastquerypos],MATCH_COMP,
			  genomicseg_ptr[lastgenomepos],/*dynprogindex*/0);
    debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
		  lastquerypos,lastgenomepos));
#endif
    --lastquerypos;
    --lastgenomepos;
  }

  /* Take care of first pair */
  pairs = Pairpool_push_existing(pairs,pairpool,pair);
  lastquerypos = querypos;
  lastgenomepos = genomepos;
  path = path->rest;

  while (path != NULL) {
    pair = path->first;
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - 1 - querypos;
    genomejump = lastgenomepos - 1 - genomepos;

    if (queryjump == 0 && genomejump == 0) {
      /* Do nothing */
    } else {
      debug5(printf("At querypos %d, saw queryjump of %d and genomejump of %d\n",querypos,queryjump,genomejump));

#ifdef PMAP
      fill = indexsize_nt - 3;
#else
      fill = indexsize_nt - 1;
#endif
      if (querypos + fill >= lastquerypos || genomepos + fill >= lastgenomepos) {
	abort();
      }

      lastquerypos = querypos + fill;
      lastgenomepos = genomepos + fill;
      while (lastquerypos > querypos) {
#ifdef PMAP
	pairs = Pairpool_push(pairs,pairpool,lastquerypos,lastgenomepos,
			      genomicuc_ptr[lastgenomepos],MATCH_COMP,
			      genomicseg_ptr[lastgenomepos],/*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",genomicuc_ptr[lastgenomepos],genomicseg_ptr[lastgenomepos],
		      lastquerypos,lastgenomepos));
#else
	pairs = Pairpool_push(pairs,pairpool,lastquerypos,lastgenomepos,
			      queryseq_ptr[lastquerypos],MATCH_COMP,
			      genomicseg_ptr[lastgenomepos],/*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
		      lastquerypos,lastgenomepos));
#endif
	--lastquerypos;
	--lastgenomepos;
      }
    }

    /* Take care of observed match */
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
    lastquerypos = querypos;
    lastgenomepos = genomepos;
    path = path->rest;
  }

  return List_reverse(pairs);
}


T
Stage2_compute (Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		Sequence_T queryntseq,
#endif
		Sequence_T genomicseg, Sequence_T genomicuc, Oligoindex_T *oligoindices,
		Gbuffer_T gbuffer, int maxoligohits, int minindexsize, int maxindexsize, Pairpool_T pairpool, 
		int sufflookback, int nsufflookback, int maxintronlen, bool crossspeciesp, Stopwatch_T stopwatch) {
  T stage2;
  int indexsize, best_indexsize = maxindexsize, high_indexsize = maxindexsize, low_indexsize = maxindexsize, indexsize_nt;
  List_T path, pairs;
  unsigned int ***mappings;
  int nfwdintrons, nrevintrons, nunkintrons;
  int cdna_direction = 0;
  int **npositions, *totalpositions, pathlength, maxconsecutive, best_maxconsecutive = 0;
  double defectrate, mapfraction, best_mapfraction;
#ifdef PMAP
  List_T p;
  Pair_T pair;
  char *queryntseq_ptr;
#endif
  int badoligos = 0, repoligos = 0, trimoligos, trim_start, trim_end;

  Stopwatch_start(stopwatch);

  mappings = (unsigned int ***) CALLOC(maxindexsize+1,sizeof(unsigned int **));
  npositions = (int **) CALLOC(maxindexsize+1,sizeof(int *));
  totalpositions = (int *) CALLOC(maxindexsize+1,sizeof(int));

  /* Find indexsize adaptively */
  indexsize = maxindexsize;
  Oligoindex_tally(oligoindices[indexsize],genomicuc,maxoligohits);
  mappings[indexsize] = Oligoindex_get_mappings(&mapfraction,&maxconsecutive,&(npositions[indexsize]),
						&(totalpositions[indexsize]),oligoindices[indexsize],queryuc);
  best_mapfraction = mapfraction;
  best_maxconsecutive = maxconsecutive;

  debug4(printf("Mapfraction at indexsize %d is %f\n",indexsize,mapfraction));
  debug4(printf("Maxconsecutive at indexsize %d is %d\n",indexsize,maxconsecutive));
  while (mapfraction < SUFFMAPFRACTION && (indexsize-1) >= minindexsize) {
    low_indexsize = --indexsize;
    Oligoindex_set_inquery(&badoligos,&repoligos,&trimoligos,&trim_start,&trim_end,oligoindices[indexsize],queryuc,/*trimp*/false);
    Oligoindex_tally(oligoindices[indexsize],genomicuc,maxoligohits);
    mappings[indexsize] = Oligoindex_get_mappings(&mapfraction,&maxconsecutive,&(npositions[indexsize]),
						  &(totalpositions[indexsize]),oligoindices[indexsize],queryuc);
    if (mapfraction > best_mapfraction) {
      best_indexsize = indexsize;
      best_mapfraction = mapfraction;
      best_maxconsecutive = maxconsecutive;
    }
    debug4(printf("Mapfraction at indexsize %d is %f\n",indexsize,mapfraction));
    debug4(printf("Maxconsecutive at indexsize %d is %d\n",indexsize,maxconsecutive));
  }

  indexsize = best_indexsize;

#ifdef PMAP
  indexsize_nt = 3*indexsize;
#else
  indexsize_nt = indexsize;
#endif

  if (best_mapfraction < MINMAPFRACTION) {
    debug4(printf("Quitting because best_mapfraction is only %f\n",best_mapfraction));
    path = NULL;
    defectrate = 1.0;
  } else if (totalpositions[best_indexsize] == 0) {
    debug4(printf("Quitting because totalpositions is zero\n"));
    path = NULL;
    defectrate = 1.0;
#ifndef PMAP
  } else if (best_maxconsecutive < 2*indexsize) {
    debug4(printf("Quitting because maxconsecutive is only %d\n",best_maxconsecutive));
    path = NULL;
    defectrate = 1.0;
#endif
  } else {
    if (best_mapfraction < 0.70) {
      sufflookback *= 4;
    } else if (best_mapfraction < 0.80) {
      sufflookback *= 3;
    } else if (best_mapfraction < 0.90) {
      sufflookback *= 2;
    }
    path = align_compute(&pathlength,&defectrate,
			 &nfwdintrons,&nrevintrons,&nunkintrons,&cdna_direction,
			 mappings[indexsize],npositions[indexsize],totalpositions[indexsize],
			 queryseq,genomicseg,genomicuc,gbuffer,indexsize,sufflookback,nsufflookback,
			 maxintronlen,pairpool,crossspeciesp,best_mapfraction);

#if 0
    /* pathlength is the number of links in the path */
    /* printf("pathlength %d vs trimlength/10 = %d\n",pathlength,Sequence_trimlength(queryseq)/10); */
    if (pathlength < Sequence_trimlength(queryseq)/10) {
      debug4(printf("Skipping because pathlength %d < trimlength %d/10\n",
		    pathlength,Sequence_trimlength(queryseq)));

      /* Don't free path, because its memory belongs to pairpool */
      path = (List_T) NULL;
    }
#endif
  }

  for (indexsize = low_indexsize; indexsize <= high_indexsize; indexsize++) {
    FREE(npositions[indexsize]);
    FREE(mappings[indexsize]);
    Oligoindex_cleanup(oligoindices[indexsize],queryuc);
  }
  FREE(totalpositions);
  FREE(npositions);
  FREE(mappings);

  if (path == NULL) {
    debug(printf("Couldn't find alignment in stage 2\n"));
    pairs = (List_T) NULL;
    stage2 = NULL;
    Stopwatch_stop(stopwatch);
  } else {
#ifdef PMAP
    pairs = convert_to_nucleotides(List_reverse(path),(char *) NULL,
				   Sequence_fullpointer(genomicseg),Sequence_fullpointer(genomicuc),
				   pairpool,indexsize_nt);
#else
    pairs = convert_to_nucleotides(List_reverse(path),Sequence_fullpointer(queryseq),
				   Sequence_fullpointer(genomicseg),pairpool,indexsize_nt);
#endif
    stage2 = Stage2_new(pairs,nfwdintrons,nrevintrons,nunkintrons,cdna_direction,
			indexsize_nt,best_mapfraction,best_maxconsecutive,defectrate);
    stage2->runtime = Stopwatch_stop(stopwatch);
  }

#ifdef PMAP
  /* Fill in queryntseq */
  queryntseq_ptr = Sequence_fullpointer(queryntseq);
  for (p = pairs; p != NULL; p = p->rest) {
    pair = p->first;
    queryntseq_ptr[pair->querypos] = pair->cdna;
  }
#endif

  /* Don't free path, because its memory belongs to pairpool */

  return stage2;
}
