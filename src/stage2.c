static char rcsid[] = "$Id: stage2.c,v 1.210 2007/04/17 00:16:42 twu Exp $";
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
#include "intlist.h"

/* #define SQUARE 1 */


#ifdef PMAP
#define LOCAL_LOOKBACK 40
#else
#define LOCAL_LOOKBACK 120
#endif

#define MAX_NACTIVE 100
#define MAX_GRAND_LOOKBACK 200

/* Penalty for genomic distances */

#define INTRON_PENALTY_UNKNOWN 8
#define INTRON_PENALTY_INCONSISTENT 16

#define NINTRON_PENALTY 8
#define NINTRON_PENALTY_MISMATCH 8

#define ENOUGH_CONSECUTIVE 32

#define INFINITE 1000000

#ifdef PMAP
#define EQUAL_DISTANCE 3
#else
#define EQUAL_DISTANCE 1	/* Making this too large can allow bad
				   global alignment to beat good local
				   alignment */
#endif
#define INTRON_DEFN 9		/* Cannot exceed 9 */

#define SCORE_FOR_RESTRICT 10


#ifdef PMAP
#define SAMPLE_INTERVAL 1
#define NT_PER_MATCH 3
#define CONSEC_POINTS_PER_MATCH 3 /* Possible increase to reward consecutiveness */
#define NONCODON_INDEL_PENALTY 15
#else
#define SAMPLE_INTERVAL 2	/* For cases where adjacentp == false.
				   Means that we can find islands of
				   9-mers */
#define NT_PER_MATCH 1
#define CONSEC_POINTS_PER_MATCH 1 /* Possible increase to reward consecutiveness */
#endif

/* Dynamic programming */
/* Can also define debug(x) as: if (querypos == XX) x */
#ifdef DEBUG
#define debug(x) if (querypos == 117533) x
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

/* Result of maxconsecutive.  May want to turn on DEBUG2 in oligoindex.c. */
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

/* revise_active */
#ifdef DEBUG6
#define debug6(x) x
#else 
#define debug6(x)
#endif

/* Shifted canonical */
#ifdef DEBUG7
#define debug7(x) x
#else 
#define debug7(x)
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
Stage2_defectrate (T this) {
  return this->defectrate;
}


int
Stage2_pathlength (T this) {
  return List_length(this->path);
}

static T
Stage2_new (List_T path, int nfwdintrons, int nrevintrons, int nunkintrons,
	    int cdna_direction, int indexsize_nt, double defectrate) {
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
      printf(" %d.%u:%d(%d,%d)",j,mappings[i][j],links[i][j].fwd_score,links[i][j].fwd_pos,links[i][j].fwd_hit);
    }
    printf("\n");
  }
  printf("\n");

  FREE(oligo);
  return;
}

#ifndef PMAP
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
      printf(" %d.%u:%d(%d,%d)-%d(%d,%d)",
	     j,mappings[i][j],links[i][j].fwd_score,
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

static void
mappings_dump_R (unsigned int **mappings, int *npositions, int length1,
		 int **active, int *firstactive, int indexsize, char *varname) {
  int querypos;
  int i, j, lastpos, hit;

  lastpos = length1 - indexsize;
  printf("%s <- matrix(c(\n",varname);
  for (querypos = 0; querypos < lastpos; querypos++) {
    if (firstactive) {
      if (mappings[querypos] != NULL) {
	hit = firstactive[querypos];
	while (hit != -1) {
	  /* Last elt is for score */
	  printf("%d,%d,\n",querypos,mappings[querypos][hit]);
	  hit = active[querypos][hit];
	}
      }
    } else {
      for (hit = 0; hit < npositions[querypos]; hit++) {
	printf("%d,%d,\n",querypos,mappings[querypos][hit]);
      }
    }
  }
  printf("),ncol=2,byrow=T)\n");

  return;
}
    

static void
best_path_dump_R (struct Link_T **links, unsigned int **mappings,
		  int querypos, int hit, bool fwdp, char *varname) {
  unsigned int position;
  int prev_querypos, prevhit, save_querypos, savehit;

  save_querypos = querypos;
  savehit = hit;

  printf("%s <- matrix(c(\n",varname);
  prev_querypos = querypos+1;
  while (querypos >= 0) {
    position = mappings[querypos][hit];

    printf("%d,%d,\n",querypos,position);

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
  }
  printf("),ncol=2,byrow=T)\n");

  querypos = save_querypos;
  hit = savehit;

  printf("%s <- matrix(c(\n","scores");
  prev_querypos = querypos+1;
  while (querypos >= 0) {
    position = mappings[querypos][hit];

    if (fwdp == true) {
      printf("%d,%d,\n",querypos,links[querypos][hit].fwd_score);
#ifndef PMAP
    } else {
      printf("%d,%d,\n",querypos,links[querypos][hit].rev_score);
#endif
    }

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
  }
  printf("),ncol=2,byrow=T)\n");

  return;
}

static void
active_bounds_dump_R (unsigned int *minactive, unsigned int *maxactive,
		      int querylength) {
  int querypos;

  printf("querypos <- 0:%d\n",querylength-1);
  printf("%s <- c(\n","minactive");
  for (querypos = 0; querypos < querylength; querypos++) {
    printf("%d,\n",minactive[querypos]);
  }
  printf(")\n");

  printf("%s <- c(\n","maxactive");
  for (querypos = 0; querypos < querylength; querypos++) {
    printf("%d,\n",maxactive[querypos]);
  }
  printf(")\n");

  return;
}


#define TEN_THOUSAND 10000.0
#define HUNDRED_THOUSAND 100000.0

static int
gendist_penalty (unsigned int diffdistance, int querydistance) {
  int penalty;

  /* Add 1 to get effect of querydistance multiplier below */
  penalty = diffdistance/HUNDRED_THOUSAND + 1;
  if (diffdistance < EQUAL_DISTANCE) {
    return penalty;
  } else if (querydistance < 8) {
    return penalty;
  } else {
    penalty *= (querydistance/8) * (querydistance/8);
    return penalty;
  }
}


static int
querydist_mismatch_penalty (int querydistance, bool crossspeciesp) {
  return 0;
  /* return (querydistance/8) * (querydistance/8); */

  if (querydistance <= 24) {
    return (querydistance+7);
  } else if (querydistance <= 40) {
    return 2*(querydistance+7);
  } else if (querydistance <= 56) {
    return 2*(querydistance+7);
  } else {
    return (querydistance+7)*3;
  }
}

#if 0
static int
diffdist_mismatch_penalty (int diffdistance, bool crossspeciesp) {
  int penalty;

  if (crossspeciesp == true) {
    if (diffdistance <= 24) {
      return (diffdistance+7)/2;
    } else if (diffdistance <= 40) {
      return (diffdistance+7)/2;
    } else if (diffdistance <= 56) {
      return (diffdistance+7);
    } else {
      return (diffdistance+7)*3/2;
    }
  } else {
#ifdef PMAP
    if (diffdistance <= 24) {
      penalty = (diffdistance+7)/4;
    } else if (diffdistance <= 40) {
      penalty = (diffdistance+7)/4;
    } else if (diffdistance <= 56) {
      penalty = (diffdistance+7)/2;
    } else {
      penalty = (diffdistance+7)*3/4;
    }
    return penalty;
#else
    if (diffdistance <= 24) {
      return (diffdistance+7)/8;
    } else if (diffdistance <= 40) {
      return (diffdistance+7)/4;
    } else if (diffdistance <= 56) {
      return (diffdistance+7)/2;
    } else {
      return (diffdistance+7)*3/4;
    }
#endif
  }
}
#endif


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
  if (leftpos > rightpos) {
    return false;
  }

  debug7(printf("Looking for shifted canonical at prevposition %d to position %d\n",prevposition,position));
  shift = 0;
  while (shift <= querydistance - indexsize_nt) {
    leftmiss = leftpos - last_leftdi[leftpos];
    rightmiss = rightpos - last_rightdi[rightpos];
    debug7(printf("%d/L%d/R%d/%d/%d",shift,leftpos,rightpos,leftmiss,rightmiss));
    if (leftmiss == 0 && rightmiss == 0) {
      debug7(printf(" => Success\n\n"));
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
  }
  debug7(printf("\n"));
  
  return false;
}


/* For PMAP, indexsize is in aa. */
static void
score_querypos (Link_T currlink, int querypos, unsigned int position, struct Link_T **links,
		unsigned int **mappings, int *npositions, int **active, int *firstactive, int *nactive,
		int grand_fwd_querypos, int grand_fwd_hit,
#ifndef PMAP
		int grand_rev_querypos, int grand_rev_hit,
#endif

		char *genomicuc_ptr, int *lastGT, int *lastAG, int *lastCT, int *lastAC, 
		int indexsize, int sufflookback, int nsufflookback, int maxintronlen,
		bool crossspeciesp) {
  Link_T prevlink;
  int best_fwd_consecutive = indexsize*NT_PER_MATCH;
  int best_fwd_score = indexsize*CONSEC_POINTS_PER_MATCH, fwd_score;
  int best_fwd_prevpos = -1, best_fwd_prevhit = -1;
  int best_fwd_intronnfwd = 0, best_fwd_intronnrev = 0, best_fwd_intronnunk = 0;
#ifndef PMAP
  int best_rev_consecutive = indexsize*NT_PER_MATCH;
  int best_rev_score = indexsize*CONSEC_POINTS_PER_MATCH, rev_score;
  int best_rev_prevpos = -1, best_rev_prevhit = -1;
  int best_rev_intronnfwd = 0, best_rev_intronnrev = 0, best_rev_intronnunk = 0;
#endif
  bool adjacentp = false;
  int prev_querypos, intronpos_lower, intronpos_upper, prevhit, finalprevhit;
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
  prevhit = firstactive[querypos-1];
  while (prevhit != -1 && (prevposition = mappings[querypos-1][prevhit]) + NT_PER_MATCH <= position) {
    if (prevposition + NT_PER_MATCH == position) {
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
      debug(printf("\tA. Adjacent qpos %d,%d at %ux%d (scores = %d -> %d, consec = %d, intr = %d-%d-%d)\n",
		   querypos-1,prevhit,position - NT_PER_MATCH,active[querypos-1][prevhit],prevlink->fwd_score,
		   best_fwd_score,best_fwd_consecutive,best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk));
#else
      debug(printf("\tA. Adjacent qpos %d,%d at %ux%d (scores = %d/%d -> %d/%d, consec = %d/%d, intr = %d-%d-%d/%d-%d-%d)\n",
		   querypos-1,prevhit,position - NT_PER_MATCH,active[querypos-1][prevhit],prevlink->fwd_score,prevlink->rev_score,
		   best_fwd_score,best_rev_score,best_fwd_consecutive,best_rev_consecutive,
		   best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		   best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk));
#endif
      prevhit = -1;		/* Exit loop */
    } else {
      prevhit = active[querypos-1][prevhit];
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
    if (intronpos_lower < 0) {
      intronpos_lower = 0;
    }
    for (prev_querypos = intronpos_upper; prev_querypos >= intronpos_lower; --prev_querypos) {
      if (nactive[prev_querypos] > MAX_NACTIVE) {
	prevhit = -1;
      } else {
#ifdef PMAP
	querydistance = (querypos - prev_querypos)*3;
#else
	querydistance = querypos - prev_querypos;
#endif
	prevhit = firstactive[prev_querypos];
      }
      while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
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
	} else {
	  if (canonicalsgn == 0) {
	    /* Extra penalty for known non-canonical intron */
	    fwd_gendistance_penalty = 
#ifndef PMAP
	      rev_gendistance_penalty = 
#endif
	      gendist_penalty(diffdistance,querydistance) + INTRON_PENALTY_UNKNOWN;
	  } else if (canonicalsgn == +1) {
	    fwd_gendistance_penalty = gendist_penalty(diffdistance,querydistance) /*- INTRON_REWARD_CONSISTENT*/;
#ifndef PMAP
	    rev_gendistance_penalty = gendist_penalty(diffdistance,querydistance) + INTRON_PENALTY_INCONSISTENT;
#endif
	  } else if (canonicalsgn == -1) {
#ifdef PMAP
	    fwd_gendistance_penalty = gendist_penalty(diffdistance,querydistance) + INTRON_PENALTY_UNKNOWN;
#else
	    fwd_gendistance_penalty = gendist_penalty(diffdistance,querydistance) + INTRON_PENALTY_INCONSISTENT;
	    rev_gendistance_penalty = gendist_penalty(diffdistance,querydistance) /*- INTRON_REWARD_CONSISTENT*/;
#endif
	  }
	}
	fwd_score -= fwd_gendistance_penalty;
#ifndef PMAP
	rev_score -= rev_gendistance_penalty;
#endif

#ifdef PMAP
	debug(printf("\tB. Intron qpos %d,%d at %ux%d (scores = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, pen %d)",
		     prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
		     prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,
		     best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		     gendistance,fwd_gendistance_penalty));
#else
	debug(printf("\tB. Intron qpos %d,%d at %ux%d (scores = %d/%d -> %d/%d, consec = %d/%d, intr = %d-%d-%d/%d-%d-%d, gendist %u, pen %d/%d)",
		     prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
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
	prevhit = active[prev_querypos][prevhit];
      }
    }
    debug(printf("\n"));

    /* C (fwd). Look back at grand querypos and hit */
    if (grand_fwd_querypos > 0 && grand_fwd_querypos + indexsize_nt <= querypos &&
	(prevposition = mappings[grand_fwd_querypos][grand_fwd_hit]) + indexsize_nt <= position) {
      prev_querypos = grand_fwd_querypos;
      prevhit = grand_fwd_hit;
#ifdef PMAP
      querydistance = (querypos - prev_querypos)*3;
#else
      querydistance = querypos - prev_querypos;
#endif
      if (querydistance < MAX_GRAND_LOOKBACK) {
	prevlink = &(links[prev_querypos][prevhit]);
	
	gendistance = position - prevposition;
	diffdistance = abs(gendistance - querydistance);

	if (diffdistance > maxintronlen) {
	  canonicalsgn = 0;
	} else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,
					  indexsize_nt,lastGT,lastAG) == true) {
	  canonicalsgn = +1;
	} else {
	  canonicalsgn = 0;
	}

	fwd_score = prevlink->fwd_score;
	if (diffdistance > maxintronlen) {
	  fwd_score -= INFINITE;
	} else if (diffdistance <= EQUAL_DISTANCE) {
	  fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_MATCH;
	} else if (diffdistance < INTRON_DEFN) {
	  /* fwd_score -= querydist_mismatch_penalty(querydistance,crossspeciesp); */
#ifdef PMAP
	  if (diffdistance % 3 != 0) {
	    fwd_score -= NONCODON_INDEL_PENALTY;
	  }
#endif
	} else {
	  /* fwd_score -= querydist_mismatch_penalty(querydistance,crossspeciesp); */
	  fwd_gendistance_penalty = gendist_penalty(diffdistance,querydistance);
	  fwd_score -= fwd_gendistance_penalty;
	  if (0 && canonicalsgn == +1) {
	    fwd_score += querydistance;
	  } else {
	    fwd_score -= NINTRON_PENALTY_MISMATCH;
	  }
	}

	debug(printf("\tC+. Fwd grand qpos %d,%d at %ux%d (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
		     prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
		     prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,
		     best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		     gendistance,querydistance,canonicalsgn));
	    
	if (fwd_score > best_fwd_score) {
	  if (diffdistance <= EQUAL_DISTANCE) {
	    best_fwd_consecutive = prevlink->fwd_consecutive + NT_PER_MATCH;
	  } else {
#ifdef PMAP
	    best_fwd_consecutive = indexsize*3;
#else
	    best_fwd_consecutive = indexsize;
#endif
	  }
	  best_fwd_score = fwd_score;
	  best_fwd_prevpos = prev_querypos;
	  best_fwd_prevhit = prevhit;
	  best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	  best_fwd_intronnrev = prevlink->fwd_intronnrev;
	  best_fwd_intronnunk = prevlink->fwd_intronnunk;
	  if (diffdistance > INTRON_DEFN) {
	    best_fwd_intronnunk++;
	  }
	  debug(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	} else {
	  debug(printf(" => Loses to %d\n",best_fwd_score));
	}
      }
    }

#ifndef PMAP
    /* C (rev). Look back at grand querypos and hit */
    if (grand_rev_querypos > 0 && grand_rev_querypos + indexsize_nt <= querypos &&
	(prevposition = mappings[grand_rev_querypos][grand_rev_hit]) + indexsize_nt <= position) {
      prev_querypos = grand_rev_querypos;
      prevhit = grand_rev_hit;
      querydistance = querypos - prev_querypos;
      if (querydistance < MAX_GRAND_LOOKBACK) {
	prevlink = &(links[prev_querypos][prevhit]);

	gendistance = position - prevposition;
	diffdistance = abs(gendistance - querydistance);

	if (diffdistance > maxintronlen) {
	  canonicalsgn = 0;
	} else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,
					  indexsize_nt,lastGT,lastAG) == true) {
	  canonicalsgn = +1;
	} else {
	  canonicalsgn = 0;
	}

	rev_score = prevlink->rev_score;
	if (diffdistance > maxintronlen) {
	  rev_score -= INFINITE;
	} else if (diffdistance <= EQUAL_DISTANCE) {
	  rev_score = prevlink->rev_score + CONSEC_POINTS_PER_MATCH;
	} else if (diffdistance < INTRON_DEFN) {
	  /* rev_score -= querydist_mismatch_penalty(querydistance,crossspeciesp); */
	} else {
	  /* rev_score -= querydist_mismatch_penalty(querydistance,crossspeciesp); */
	  rev_gendistance_penalty = gendist_penalty(diffdistance,querydistance);
	  rev_score -= rev_gendistance_penalty;
	  if (0 && canonicalsgn == +1) {
	    rev_score += querydistance;
	  } else {
	    rev_score -= NINTRON_PENALTY_MISMATCH;
	  }
	}

	debug(printf("\tC-. Rev grand qpos %d,%d at %ux%d (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
		     prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
		     prevlink->rev_score,rev_score,prevlink->rev_consecutive,
		     best_rev_intronnrev,best_rev_intronnrev,best_rev_intronnunk,
		     gendistance,querydistance,canonicalsgn));
	    
	if (rev_score > best_rev_score) {
	  if (diffdistance <= EQUAL_DISTANCE) {
	    best_rev_consecutive = prevlink->rev_consecutive + NT_PER_MATCH;
	  } else {
	    best_rev_consecutive = indexsize;
	  }
	  best_rev_score = rev_score;
	  best_rev_prevpos = prev_querypos;
	  best_rev_prevhit = prevhit;
	  best_rev_intronnfwd = prevlink->rev_intronnfwd;
	  best_rev_intronnrev = prevlink->rev_intronnrev;
	  best_rev_intronnunk = prevlink->rev_intronnunk;
	  if (diffdistance > INTRON_DEFN) {
	    best_rev_intronnunk++;
	  }
	  debug(printf(" => Best rev at %d (consec = %d)\n",rev_score,best_rev_consecutive));
	} else {
	  debug(printf(" => Loses to %d\n",best_rev_score));
	}
      }
    }
#endif

    /* Use intronpos_upper instead of intronpos_lower because correct position is uncertain in PMAP */
    mismatch_start = intronpos_upper - 1;

    /* D (fwd). Evaluate for mismatches (all other previous querypos) */
    /* Set parameters */
    if (adjacentp == true) {
      /* Sample at larger interval and not so far back */
      nlookback = 1;
      lookback = sufflookback/2;
      sample_interval = indexsize;
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
      if (nactive[prev_querypos] > MAX_NACTIVE) {
	prevhit = -1;
      } else if ((prevhit = firstactive[prev_querypos]) != -1) {
	nseen++;
#ifdef PMAP
	querydistance = (querypos - prev_querypos)*3;
#else
	querydistance = querypos - prev_querypos;
#endif
      }
      while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
	prevlink = &(links[prev_querypos][prevhit]);
	
	gendistance = position - prevposition;
	diffdistance = abs(gendistance - querydistance);

	if (diffdistance > maxintronlen) {
	  canonicalsgn = 0;
	} else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,
					  indexsize_nt,lastGT,lastAG) == true) {
	  canonicalsgn = +1;
	} else {
	  canonicalsgn = 0;
	}

	fwd_score = prevlink->fwd_score;
	if (diffdistance > maxintronlen) {
	  fwd_score -= INFINITE;
	} else if (diffdistance <= EQUAL_DISTANCE) {
	  fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_MATCH;
	} else if (diffdistance < INTRON_DEFN) {
	  /* fwd_score -= querydist_mismatch_penalty(querydistance,crossspeciesp); */
#ifdef PMAP
	  if (diffdistance % 3 != 0) {
	    fwd_score -= NONCODON_INDEL_PENALTY;
	  }
#endif
	} else {
	  /* fwd_score -= querydist_mismatch_penalty(querydistance,crossspeciesp); */
	  fwd_score -= gendist_penalty(diffdistance,querydistance);
	  if (0 && canonicalsgn == +1) {
	    fwd_score += querydistance;
	  } else {
	    fwd_score -= NINTRON_PENALTY_MISMATCH;
	  }
	}

	debug(printf("\tD+. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
		     prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
		     prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,
		     best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		     gendistance,querydistance,canonicalsgn));
	    
	if (fwd_score > best_fwd_score) {
	  if (diffdistance <= EQUAL_DISTANCE) {
	    best_fwd_consecutive = prevlink->fwd_consecutive + NT_PER_MATCH;
	  } else {
#ifdef PMAP
	    best_fwd_consecutive = indexsize*3;
#else
	    best_fwd_consecutive = indexsize;
#endif
	  }
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
	prevhit = active[prev_querypos][prevhit];
      }
    }

#ifndef PMAP
    /* D (rev). Evaluate for mismatches (all other previous querypos) */
    /* Set parameters */
    if (adjacentp == true) {
      /* Sample at larger interval and not so far back */
      nlookback = 1;
      lookback = sufflookback/2;
      sample_interval = indexsize;
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
      if (nactive[prev_querypos] > MAX_NACTIVE) {
	prevhit = -1;
      } else if ((prevhit = firstactive[prev_querypos]) != -1) {
	nseen++;
	querydistance = querypos - prev_querypos;
      }
      while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
	prevlink = &(links[prev_querypos][prevhit]);
	  
	gendistance = position - prevposition;
	diffdistance = abs(gendistance - querydistance);

	if (diffdistance > maxintronlen) {
	  canonicalsgn = 0;
	} else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,
					  indexsize_nt,lastCT,lastAC) == true) {
	  canonicalsgn = -1;
	} else {
	  canonicalsgn = 0;
	}
	  
	rev_score = prevlink->rev_score;
	if (diffdistance > maxintronlen) {
	  rev_score -= INFINITE;
	} else if (diffdistance <= EQUAL_DISTANCE) {
	  rev_score = prevlink->rev_score + CONSEC_POINTS_PER_MATCH;
	} else if (diffdistance < INTRON_DEFN) {
	  /* rev_score -= querydist_mismatch_penalty(querydistance,crossspeciesp); */
	} else {
	  /* rev_score -= querydist_mismatch_penalty(querydistance,crossspeciesp); */
	  rev_score -= gendist_penalty(diffdistance,querydistance);
	  if (0 && canonicalsgn == -1) {
	    rev_score += querydistance;
	  } else {
	    rev_score -= NINTRON_PENALTY_MISMATCH;
	  }
	}

	debug(printf("\tD-. Rev mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
		     prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
		     prevlink->rev_score,rev_score,prevlink->rev_consecutive,
		     best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk,
		     gendistance,querydistance,canonicalsgn));
	    
	if (rev_score > best_rev_score) {
	  if (diffdistance <= EQUAL_DISTANCE) {
	    best_rev_consecutive = prevlink->rev_consecutive + NT_PER_MATCH;
	  } else {
	    best_rev_consecutive = indexsize;
	  }
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
	prevhit = active[prev_querypos][prevhit];
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


static void
revise_active (int **active, int *firstactive, int *nactive, 
	       int low_hit, int high_hit, struct Link_T **links, int querypos) {
  int best_fwd_score, best_rev_score, fwd_threshold, rev_threshold, score;
  int hit, *ptr;

  debug6(printf("Revising querypos %d from low_hit %d to high_hit %d.  Scores:\n",querypos,low_hit,high_hit));
  if ((hit = low_hit) >= high_hit) {
    firstactive[querypos] = -1;
  } else {
    debug6(printf("At hit %d, fwd_score is %d\n",hit,links[querypos][hit].fwd_score));
    best_fwd_score = links[querypos][hit].fwd_score;
#ifndef PMAP
    best_rev_score = links[querypos][hit].rev_score;
#endif
    for (hit++; hit < high_hit; hit++) {
      debug6(printf("At hit %d, fwd_score is %d\n",hit,links[querypos][hit].fwd_score));
      if ((score = links[querypos][hit].fwd_score) > best_fwd_score) {
	best_fwd_score = score;
      }
#ifndef PMAP
      if ((score = links[querypos][hit].rev_score) > best_rev_score) {
	best_rev_score = score;
      }
#endif
    }

    fwd_threshold = best_fwd_score - SCORE_FOR_RESTRICT;
    if (fwd_threshold < 0) {
      fwd_threshold = 0;
    }

#ifndef PMAP
    rev_threshold = best_rev_score - SCORE_FOR_RESTRICT;
    if (rev_threshold < 0) {
      rev_threshold = 0;
    }
#endif

    nactive[querypos] = 0;
    ptr = &(firstactive[querypos]);
    hit = low_hit;
    while (hit < high_hit) {
      while (hit < high_hit && links[querypos][hit].fwd_score <= fwd_threshold
#ifndef PMAP
	     && links[querypos][hit].rev_score <= rev_threshold
#endif
	     ) {
	hit++;
      }
      *ptr = hit;
      if (hit < high_hit) {
	nactive[querypos] += 1;
	ptr = &(active[querypos][hit]);
	hit++;
      }
    }
    *ptr = -1;
  }

  debug6(
	 printf("Valid hits (%d) at querypos %d:",nactive[querypos],querypos);
	 hit = firstactive[querypos];
	 while (hit != -1) {
	   printf(" %d",hit);
	   hit = active[querypos][hit];
	 }
	 printf("\n");
	 );

  return;
}

static int **
intmatrix_new (int length1, int *lengths2, int totallength) {
  int **matrix;
  int i;
  
  matrix = (int **) CALLOC(length1,sizeof(int *));
  matrix[0] = (int *) CALLOC(totallength,sizeof(int));
  for (i = 1; i < length1; i++) {
    if (lengths2[i-1] <= 0) {
      matrix[i] = matrix[i-1];
    } else {
      matrix[i] = &(matrix[i-1][lengths2[i-1]]);
    }
  }
  return matrix;
}



/* For PMAP, indexsize is in aa. */
static int
align_compute_scores (bool *fwdp, int *grand_querypos, int *grand_hit,
		      struct Link_T **links, unsigned int **mappings, int *npositions, int totalpositions,
		      unsigned int *minactive, unsigned int *maxactive, int querystart, int queryend,
		      int querylength, Sequence_T genomicuc, Gbuffer_T gbuffer, int indexsize,
		      int sufflookback, int nsufflookback, int maxintronlen, bool crossspeciesp,
		      char *queryseq_ptr, bool debug_graphic_p) {
  int result;
  Link_T termlink = NULL, currlink;
  int querypos, indexsize_nt, hit, low_hit, high_hit, score;
  int grand_fwd_score, grand_fwd_querypos, grand_fwd_hit, best_fwd_hit, best_fwd_score;
#ifndef PMAP
  int grand_rev_score, grand_rev_querypos, grand_rev_hit, best_rev_hit, best_rev_score;
#endif
  int **active, *firstactive, *nactive;
  unsigned int position;
  char *genomicuc_ptr;
  int *lastGT, *lastAG, *lastCT, *lastAC;
#ifdef DEBUG
  char *oligo;
#endif

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

#ifdef DEBUG
  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
#endif
  debug0(printf("Querystart = %d, queryend = %d\n",querystart,queryend));

  genomicuc_ptr = Sequence_fullpointer(genomicuc);
  find_canonical_dinucleotides(gbuffer,genomicuc_ptr,Sequence_fulllength(genomicuc));
  lastGT = Gbuffer_lastGT(gbuffer);
  lastAG = Gbuffer_lastAG(gbuffer);
  lastCT = Gbuffer_lastCT(gbuffer);
  lastAC = Gbuffer_lastAC(gbuffer);

  active = intmatrix_new(querylength,npositions,totalpositions);
  firstactive = (int *) CALLOC(querylength,sizeof(int));
  nactive = (int *) CALLOC(querylength,sizeof(int));

  /* Initialize */
  for (querypos = 0; querypos < querystart; querypos++) {
    firstactive[querypos] = -1;
    nactive[querypos] = 0;
  }
  while (querypos <= queryend - indexsize && npositions[querypos] <= 0) {
    firstactive[querypos] = -1;
    querypos++;
  }
  if (querypos <= queryend - indexsize) {
    for (hit = 0; hit < npositions[querypos]; hit++) {
      currlink = &(links[querypos][hit]);
#ifdef PMAP    
      currlink->fwd_pos = currlink->fwd_hit = -1;
      currlink->fwd_score = indexsize_nt;
      currlink->fwd_consecutive = indexsize_nt;
#else
      currlink->fwd_pos = currlink->fwd_hit = -1;
      currlink->rev_pos = currlink->rev_hit = -1;
      currlink->fwd_score = indexsize_nt;
      currlink->rev_score = indexsize_nt;
      currlink->fwd_consecutive = indexsize_nt;
      currlink->rev_consecutive = indexsize_nt;
#endif
    }
    revise_active(active,firstactive,nactive,0,npositions[querypos],links,querypos);
  }

  grand_fwd_score = 0;
  grand_fwd_querypos = -1;
  grand_fwd_hit = -1;
#ifndef PMAP
  grand_rev_score = 0;
  grand_rev_querypos = -1;
  grand_rev_hit = -1;
#endif
  for (querypos++; querypos <= queryend - indexsize; querypos++) {
    best_fwd_score = 0;
    best_fwd_hit = -1;
#ifndef PMAP
    best_rev_score = 0;
    best_rev_hit = -1;
#endif
    
    hit = 0;
    while (hit < npositions[querypos] && mappings[querypos][hit] < minactive[querypos]) {
      hit++;
    }
    low_hit = hit;
    while (hit < npositions[querypos] && mappings[querypos][hit] < maxactive[querypos]) {
      hit++;
    }
    high_hit = hit;
    debug(printf("Querypos %d has hit %d..%d out of %d (minactive = %u, maxactive = %u)\n",
		 querypos,low_hit,high_hit-1,npositions[querypos],minactive[querypos],maxactive[querypos]));

    /* Can't use nactive yet, so use high_hit - low_hit */
    if (high_hit - low_hit >= MAX_NACTIVE) {
      firstactive[querypos] = -1;
    } else {
      for (hit = low_hit; hit < high_hit; hit++) {
	currlink = &(links[querypos][hit]);
	position = mappings[querypos][hit];

	debug(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
	debug(printf("Finding link at querypos %d,%d at %ux%d (%s)\n",
		     querypos,hit,position,active[querypos][hit],oligo));

	score_querypos(currlink,querypos,position,links,mappings,npositions,active,firstactive,nactive,
		       grand_fwd_querypos,grand_fwd_hit,
#ifndef PMAP
		       grand_rev_querypos,grand_rev_hit,
#endif
		       genomicuc_ptr,lastGT,lastAG,lastCT,lastAC,indexsize,sufflookback,nsufflookback,
		       maxintronlen,crossspeciesp);
	if (currlink->fwd_score > best_fwd_score) {
	  best_fwd_score = currlink->fwd_score;
	  best_fwd_hit = hit;
	}
#ifndef PMAP
	if (currlink->rev_score > best_rev_score) {
	  best_rev_score = currlink->rev_score;
	  best_rev_hit = hit;
	}
#endif
      }
      
      /* Use >= to favor longer path in case of ties */
      if (best_fwd_hit >= 0 && best_fwd_score >= grand_fwd_score) {
	grand_fwd_score = best_fwd_score;
	grand_fwd_querypos = querypos;
	grand_fwd_hit = best_fwd_hit;
	debug(termlink = &(links[querypos][best_fwd_hit]));
	debug(printf("At querypos %d, revising grand fwd to be hit %d with score of %d (pointing back to %d,%d)\n",
		     querypos,best_fwd_hit,best_fwd_score,termlink->fwd_pos,termlink->fwd_hit));
      }

#ifndef PMAP
      /* Use >= to favor longer path in case of ties */
      if (best_rev_hit >= 0 && best_rev_score >= grand_rev_score) {
	grand_rev_score = best_rev_score;
	grand_rev_querypos = querypos;
	grand_rev_hit = best_rev_hit;
      }
#endif

      revise_active(active,firstactive,nactive,low_hit,high_hit,links,querypos);
    }

    /* These are the final active oligomers, after pruning by score */
    if (debug_graphic_p == true) {
      mappings_dump_R(mappings,npositions,querylength,active,firstactive,indexsize,"active.mers");
    }
  }

  FREE(nactive);
  FREE(firstactive);
  FREE(active[0]);
  FREE(active);


  /* Grand winner */
#ifdef PMAP
  if (grand_fwd_querypos < 0) {
    result = 0;
  } else {
    result = 1;
    *fwdp = true;
    *grand_querypos = grand_fwd_querypos;
    *grand_hit = grand_fwd_hit;
    debug(printf("==> Found best score %d at %d,%d\n",grand_fwd_score,grand_fwd_querypos,grand_fwd_hit));
  }
#else
  if (grand_fwd_querypos < 0 && grand_rev_querypos < 0) {
    result = 0;
  } else {
    result = 1;
    if (best_fwd_score >= best_rev_score) {
      *fwdp = true;
      *grand_querypos = grand_fwd_querypos;
      *grand_hit = grand_fwd_hit;
      debug(printf("==> Found best score (fwd) %d at %d,%d\n",grand_fwd_score,grand_fwd_querypos,grand_fwd_hit));
    } else {
      *fwdp = false;
      *grand_querypos = grand_rev_querypos;
      *grand_hit = grand_rev_hit;
      debug(printf("==> Found best score (rev) %d at %d,%d\n",grand_rev_score,grand_rev_querypos,grand_rev_hit));
    }
  }
#endif

#ifdef PMAP
  if (grand_fwd_score > 3*querylength) {
    abort();
  }
#else
  if (grand_fwd_score > querylength) {
    abort();
  }
#endif

  debug(FREE(oligo));

  return result;
}



/* Performs dynamic programming.  For PMAP, indexsize is in aa. */
static List_T
align_compute (int *pathlength, int *start_querypos, int *end_querypos, double *defectrate,
	       int *nfwdintrons, int *nrevintrons, int *nunkintrons, int *cdna_direction,
	       unsigned int **mappings, int *npositions, int totalpositions, unsigned int *minactive, unsigned int *maxactive,
	       Sequence_T queryseq, Sequence_T genomicseg, Sequence_T genomicuc, Gbuffer_T gbuffer, int indexsize,
	       int sufflookback, int nsufflookback, int maxintronlen, Pairpool_T pairpool, bool crossspeciesp,
	       bool debug_graphic_p) {
  List_T path = NULL;
  struct Link_T **links;
  Link_T termlink;
  bool fwdp;
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

  /* These are all oligomers */
  if (debug_graphic_p == true) {
    mappings_dump_R(mappings,npositions,querylength,/*active*/NULL,/*firstactive*/NULL,indexsize,"all.mers");
  }
  
  if (align_compute_scores(&fwdp,&querypos,&hit,links,mappings,npositions,totalpositions,
			   minactive,maxactive,querystart,queryend,querylength,genomicuc,gbuffer,indexsize,
			   sufflookback,nsufflookback,maxintronlen,crossspeciesp,queryseq_ptr,
			   debug_graphic_p) == 0) {
    path = (List_T) NULL;
    *start_querypos = *end_querypos = 0;

  } else {

#ifdef PMAP
    debug1(Linkmatrix_print_fwd(links,mappings,querylength,npositions,queryseq_ptr,indexsize));
#else
    debug1(Linkmatrix_print_both(links,mappings,querylength,npositions,queryseq_ptr,indexsize));
#endif

    *start_querypos = querypos;
    termlink = &(links[querypos][hit]);
    if (fwdp == true) {
      *nfwdintrons = termlink->fwd_intronnfwd;
      *nrevintrons = termlink->fwd_intronnrev;
      *nunkintrons = termlink->fwd_intronnunk;
#ifndef PMAP
    } else {
      *nfwdintrons = termlink->rev_intronnfwd;
      *nrevintrons = termlink->rev_intronnrev;
      *nunkintrons = termlink->rev_intronnunk;
#endif
    }

    if (debug_graphic_p == true) {
      best_path_dump_R(links,mappings,querypos,hit,fwdp,"best.path");
      printf("plot(all.mers,col=\"black\")\n");
      printf("points(active.mers,col=\"red\")\n");
      printf("points(best.path,col=\"green\")\n");
      printf("lines(querypos,minactive,col=\"blue\")\n");
      printf("lines(querypos,maxactive,col=\"blue\")\n");
    }

    prevposition = mappings[querypos][hit]+1;
    prev_querypos = querypos+1;
    while (querypos >= 0) {
      position = mappings[querypos][hit];

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

#ifdef DEBUG0
      debug0(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
      if (fwdp == true) {
	debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d, intr = %d(+)/%d(-)/%d(?)\n",
		      querypos,hit,oligo,position,
		      links[querypos][hit].fwd_score,links[querypos][hit].fwd_consecutive,
		      links[querypos][hit].fwd_intronnfwd,links[querypos][hit].fwd_intronnrev,
		      links[querypos][hit].fwd_intronnunk));
#ifndef PMAP
      } else {
	debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d, intr = %d(+)/%d(-)/%d(?)\n",
		      querypos,hit,oligo,position,
		      links[querypos][hit].rev_score,links[querypos][hit].rev_consecutive,
		      links[querypos][hit].rev_intronnfwd,links[querypos][hit].rev_intronnrev,
		      links[querypos][hit].rev_intronnunk));
#endif
      }
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
      *end_querypos = querypos;
    }

    if (nconsecutive == 0) {
      *defectrate = 1.0;
    } else {
      *defectrate = (double) (ndefects)/(double) (nconsecutive + ndefects);
    }
    *pathlength = nconsecutive + nintrons + ndefects;

    if (fwdp == true) {
      *cdna_direction = +1;
    } else {
      *cdna_direction = -1;
    }
  }
  Linkmatrix_free(&links);

  debug0(FREE(oligo));

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
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump, fill, default_fill;

  pair = path->first;
  querypos = pair->querypos;
  genomepos = pair->genomepos;

#ifdef PMAP
  default_fill = indexsize_nt - 3;
#else
  default_fill = indexsize_nt - 1;
#endif

  lastquerypos = querypos + default_fill;
  lastgenomepos = genomepos + default_fill;
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

      if (querypos + default_fill >= lastquerypos || genomepos + default_fill >= lastgenomepos) {
	if (lastquerypos - querypos < lastgenomepos - genomepos) {
	  fprintf(stderr,"Partial fill from querypos %d to %d (genomepos goes from %u to %u)\n",
		  querypos,lastquerypos,genomepos,lastgenomepos);
	  abort();
	  fill = lastquerypos - querypos - 1;
	} else {
	  fprintf(stderr,"Partial fill from genomepos %u to %u (querypos goes from %d to %d)\n",
		  genomepos,lastgenomepos,querypos,lastquerypos);
	  abort();
	  fill = lastgenomepos - genomepos - 1;
	}
      } else {
	fill = default_fill;
      }

      lastquerypos = querypos + fill;
      lastgenomepos = genomepos + fill;
      debug5(printf("  Fill from querypos %d down to %d\n",lastquerypos,querypos));
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
		Sequence_T genomicseg, Sequence_T genomicuc, Oligoindex_T oligoindex,
		Gbuffer_T gbuffer, int maxoligohits, int minindexsize, int maxindexsize, Pairpool_T pairpool, 
		Intpool_T intpool, Diagpool_T diagpool, int sufflookback, int nsufflookback, int maxintronlen,
		bool crossspeciesp, bool debug_graphic_p, Stopwatch_T stopwatch) {
  T stage2;
  int indexsize, indexsize_nt;
  List_T path = NULL, pairs;
  unsigned int **mappings;
  int nfwdintrons, nrevintrons, nunkintrons;
  int cdna_direction = 0;
  int *npositions, totalpositions, pathlength, maxlength, trimlength,
    start_querypos, end_querypos;
  int ndiagonals;
  unsigned int *minactive, *maxactive;
  int querylength, genomiclength;
  double defectrate, pct_coverage, pct_clear_coverage;
  int badoligos = 0, repoligos = 0, trimoligos, trim_start, trim_end;
  int maxnconsecutive;

  Stopwatch_start(stopwatch);

  indexsize = minindexsize;
  querylength = Sequence_fulllength(queryuc);
  genomiclength = Sequence_fulllength(genomicseg);
  trimlength = Sequence_trimlength(queryseq);
  Oligoindex_tally(oligoindex,genomicuc,queryuc,maxoligohits);

  minactive = (unsigned int *) CALLOC(querylength,sizeof(unsigned int));
  maxactive = (unsigned int *) CALLOC(querylength,sizeof(unsigned int));

  if (debug_graphic_p == true) {
    printf("par(mfrow=c(1,2),cex=0.2)\n");
  }

  mappings = Oligoindex_get_mappings(&ndiagonals,&maxnconsecutive,&pct_coverage,&pct_clear_coverage,
				     &npositions,&totalpositions,minactive,maxactive,
				     oligoindex,queryuc,genomicuc,diagpool,
				     LOCAL_LOOKBACK,debug_graphic_p);

  if (totalpositions == 0) {
    debug2(fprintf(stderr,"Quitting because totalpositions is zero\n"));
    path = (List_T) NULL;
  } else if (trimlength > 20000 & (maxnconsecutive < SUFFNCONSECUTIVE || pct_clear_coverage < SUFF_PCTCOVERAGE)) {
    debug2(fprintf(stderr,"Quitting because maxnconsecutive is only %d or pct_coverage is only %f\n",
		   maxnconsecutive,pct_clear_coverage));
    path = (List_T) NULL;
  } else if (trimlength > 120 && (maxnconsecutive < SUFFNCONSECUTIVE || pct_coverage < SUFF_PCTCOVERAGE)) {
    debug2(fprintf(stderr,"Quitting because maxnconsecutive is only %d or pct_coverage is only %f\n",
		   maxnconsecutive,pct_coverage));
    path = (List_T) NULL;
  } else {
    debug2(fprintf(stderr,"Proceeding because maxnconsecutive is %d and pct_coverage is %f\n",
		   maxnconsecutive,pct_coverage));

    if (debug_graphic_p == true) {
      active_bounds_dump_R(minactive,maxactive,querylength);
      printf("lines(querypos,minactive,col=\"blue\")\n");
      printf("lines(querypos,maxactive,col=\"blue\")\n");
    }

    path = align_compute(&pathlength,&start_querypos,&end_querypos,&defectrate,
			 &nfwdintrons,&nrevintrons,&nunkintrons,&cdna_direction,
			 mappings,npositions,totalpositions,minactive,maxactive,
			 queryseq,genomicseg,genomicuc,gbuffer,indexsize,sufflookback,nsufflookback,
			 maxintronlen,pairpool,crossspeciesp,debug_graphic_p);
    debug4(fprintf(stderr,"aligned from %d to %d\n",start_querypos,end_querypos));
  }

  FREE(maxactive);
  FREE(minactive);
  FREE(npositions);
  FREE(mappings);		/* Don't need to free contents of mappings */
  Oligoindex_cleanup(oligoindex,queryuc);

#ifdef PMAP
  indexsize_nt = 3*indexsize;
#else
  indexsize_nt = indexsize;
#endif

  if (path == NULL) {
    debug4(printf("Couldn't find alignment in stage 2\n"));
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
    stage2 = Stage2_new(pairs,nfwdintrons,nrevintrons,nunkintrons,cdna_direction,indexsize_nt,defectrate);
    stage2->runtime = Stopwatch_stop(stopwatch);
  }

  /* Don't need to free path, because its memory belongs to pairpool */

  return stage2;
}


