static char rcsid[] = "$Id: stage2.c 33407 2011-01-06 21:52:30Z twu $";
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
#include "diag.h"

/* #define SQUARE 1 */

#define SUFF_PCTCOVERAGE_OLIGOINDEX 0.90

/* #define SUFF_PCTCOVERAGE_STAGE2 0.10 */
#define SUFF_NCOVERED 200
#define SUFF_MAXNCONSECUTIVE 20

#define MAX_NACTIVE 100	/* 100 previously considered too low, but may
			   be okay in conjunction with
			   diagonalization */
#define MAX_GRAND_LOOKBACK 200

/* Penalty for genomic distances */

#define INTRON_PENALTY_UNKNOWN 8
#define INTRON_PENALTY_INCONSISTENT 16

#define NINTRON_PENALTY_MISMATCH 8

#define ENOUGH_CONSECUTIVE 32

#define INFINITE 1000000

#ifdef PMAP
#define EQUAL_DISTANCE 3
#else
#define EQUAL_DISTANCE 6	/* Should be allowance for indels. */
#endif


#define INTRON_DEFN 9		/* Cannot exceed 9 */
#define EXON_DEFN 30
#define MAX_SKIPPED 3

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
#define NT_PER_CODON 3
#define CONSEC_POINTS_PER_MATCH 1 /* Possible increase to reward consecutiveness */
#define CONSEC_POINTS_PER_CODON 3 /* Possible increase to reward consecutiveness */
#endif

#define SHIFT_EXTRA 10

/* General */
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

/* For generating a graph */
#ifdef DEBUG3
#define debug3(x) x
#else 
#define debug3(x)
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

/* Dynamic programming */
/* Can also define debug9(x) as: if (querypos == XX) x */
#ifdef DEBUG9
#define debug9(x) x
#else 
#define debug9(x)
#endif

/* Grand winner */
#ifdef DEBUG10
#define debug10(x) x
#else 
#define debug10(x)
#endif



#define T Stage2_T
struct T {
  List_T path;
  int nfwdintrons;
  int ncanonical;
  int nnoncanonical;
  int cdna_direction;

  double diag_runtime;
  double align_runtime;
  int indexsize;
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
Stage2_diag_runtime (T this) {
  return this->diag_runtime;
}

double
Stage2_align_runtime (T this) {
  return this->align_runtime;
}

int
Stage2_indexsize (T this) {
  return this->indexsize;
}

int
Stage2_pathlength (T this) {
  return List_length(this->path);
}

#if 0
static T
Stage2_new (List_T path, int nfwdintrons, int nrevintrons, int nunkintrons,
	    int cdna_direction, int indexsize_nt) {
  T new = (T) MALLOC(sizeof(*new));

  new->path = path;
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
#endif	    

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
Linkmatrix_1d_new (int length1, int *lengths2, int totallength) {
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
Linkmatrix_1d_free (struct Link_T ***links) {
  FREE((*links)[0]);
  FREE(*links);
  return;
}


static struct Link_T **
Linkmatrix_2d_new (int length1, int *lengths2) {
  struct Link_T **links;
  int i;
  
  links = (struct Link_T **) CALLOC(length1,sizeof(struct Link_T *));
  for (i = 0; i < length1; i++) {
    if (lengths2[i] <= 0) {
      links[i] = (struct Link_T *) NULL;
    } else {
      links[i] = (struct Link_T *) CALLOC(lengths2[i],sizeof(struct Link_T));
    }
  }
  return links;
}

static void
Linkmatrix_2d_free (struct Link_T ***links, int length1) {
  int i;

  for (i = 0; i < length1; i++) {
    if ((*links)[i]) {
      FREE((*links)[i]);
    }
  }
  FREE(*links);
  return;
}



#ifdef DEBUG1
#ifdef PMAP
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

#else

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
#endif

static void
mappings_dump_R (unsigned int **mappings, int *npositions, int length1,
		 int **active, int *firstactive, int indexsize, char *varname) {
  int querypos;
  int lastpos, hit;
  bool printp = false;

  lastpos = length1 - indexsize;
  printf("%s <- matrix(c(\n",varname);
  for (querypos = 0; querypos < lastpos; querypos++) {
    if (firstactive) {
      if (mappings[querypos] != NULL) {
	hit = firstactive[querypos];
	while (hit != -1) {
	  /* Last elt is for score */
	  if (printp == false) {
	    printp = true;
	  } else {
	    printf(",\n");
	  }
	  printf("%d,%d,%d,%d",querypos,mappings[querypos][hit],
		 hit,active[querypos][hit]);
	  hit = active[querypos][hit];
	}
      }
    } else {
      for (hit = 0; hit < npositions[querypos]; hit++) {
	if (printp == false) {
	  printp = true;
	} else {
	  printf(",\n");
	}
	printf("%d,%d,%d",querypos,mappings[querypos][hit],hit);
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
  bool printp = false;

  save_querypos = querypos;
  savehit = hit;

  printf("%s <- matrix(c(\n",varname);
  prev_querypos = querypos+1;
  while (querypos >= 0) {
    position = mappings[querypos][hit];

    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    printf("%d,%d",querypos,position);

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

  printp = false;
  printf("%s <- matrix(c(\n","scores");
  prev_querypos = querypos+1;
  while (querypos >= 0) {
    position = mappings[querypos][hit];

    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    if (fwdp == true) {
      printf("%d,%d",querypos,links[querypos][hit].fwd_score);
#ifndef PMAP
    } else {
      printf("%d,%d",querypos,links[querypos][hit].rev_score);
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
  bool printp = false;

  printf("querypos <- 0:%d\n",querylength-1);
  printf("%s <- c(\n","minactive");
  for (querypos = 0; querypos < querylength; querypos++) {
    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    printf("%d",minactive[querypos]);
  }
  printf(")\n");

  printp = false;
  printf("%s <- c(\n","maxactive");
  for (querypos = 0; querypos < querylength; querypos++) {
    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    printf("%d",maxactive[querypos]);
  }
  printf(")\n");

  return;
}


#define TEN_THOUSAND 10000.0
#define HUNDRED_THOUSAND 100000.0
#define ONE_MILLION 1000000.0

static int
gendist_penalty (unsigned int diffdistance, int querydistance) {
  int penalty;

  /* Add 1 to get effect of querydistance multiplier below */
  penalty = diffdistance/HUNDRED_THOUSAND + 1;
  if (diffdistance <= EQUAL_DISTANCE) {
    return 0;			/* was penalty */
#if 0
  } else if (querydistance < 8) {
    return penalty;
#endif
  } else {
    penalty *= (querydistance/8) * (querydistance/8);
    return penalty;
  }
}

static int
diffdist_penalty (unsigned int diffdistance) {
  int penalty;

  penalty = diffdistance/TEN_THOUSAND + 1;
  return penalty;
}

/* querydistance already has indexsize_nt subtracted */
static int
querydist_penalty (int querydistance) {
#ifdef PMAP
  return querydistance/2;
#else
  return 0;
#endif
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
find_canonical_dinucleotides (int *lastGT, int *lastAG, 
#ifndef PMAP
			      int *lastCT, int *lastAC,
#endif
			      char *genomicuc_ptr, int genomiclength) {
  int pos, GT = -1, AG = -1;
#ifndef PMAP
  int CT = -1, AC = -1;
#endif
  char c1, c2;

  lastAG[0] = lastAG[1] = lastAG[2] = -1;
#ifndef PMAP
  lastAC[0] = lastAC[1] = lastAC[2] = -1;
#endif

  pos = 0;
  c1 = genomicuc_ptr[pos+1];
  while (pos <= genomiclength-4) {
    if ((c2 = genomicuc_ptr[pos+2]) == 'T') {
      if (c1 == 'G') {
	GT = pos;
#ifndef PMAP
      } else if (c1 == 'C') {
	CT = pos;
#endif
      }
      
      lastGT[pos] = lastGT[pos+1] = GT;
      lastAG[pos+3] = lastAG[pos+4] = AG;
#ifndef PMAP
      lastCT[pos] = lastCT[pos+1] = CT;
      lastAC[pos+3] = lastAC[pos+4] = AC;
#endif

      pos += 2;
      c1 = genomicuc_ptr[pos+1];

    } else {
      if (c2 == 'G') {
	if (c1 == 'A') {
	  AG = pos + 3;
	}
#ifndef PMAP
      } else if (c2 == 'C') {
	if (c1 == 'A') {
	  AC = pos + 3;
	}
#endif
      }
	
      lastGT[pos] = GT;
      lastAG[pos+3] = AG;
#ifndef PMAP
      lastCT[pos] = CT;
      lastAC[pos+3] = AC;
#endif

      pos++;
      c1 = c2;
    }
  }

  /* pos is now either genomiclength-3 (one more dinucleotide) or genomiclength-2 (none) */
  if (pos == genomiclength-3) {
    if ((c2 = genomicuc_ptr[pos+2]) == 'T') {
      if (c1 == 'G') {
	GT = pos;
#ifndef PMAP
      } else if (c1 == 'C') {
	CT = pos;
#endif
      }
    } else {
      if (c2 == 'G') {
	if (c1 == 'A') {
	  AG = pos + 3;
	}
#ifndef PMAP
      } else if (c2 == 'C') {
	if (c1 == 'A') {
	  AC = pos + 3;
	}
#endif
      }
    }
    
    lastGT[pos] = GT;
    lastAG[pos+3] = AG;
#ifndef PMAP
    lastCT[pos] = CT;
    lastAC[pos+3] = AC;
#endif
  }

#ifdef SHIFT_EXTRA
  /* Fill in rest */
  pos++;
  while (pos < genomiclength - 3 + SHIFT_EXTRA) {
    lastGT[pos] = GT;
    lastAG[pos+3] = AG;
#ifndef PMAP
    lastCT[pos] = CT;
    lastAC[pos+3] = AC;
#endif
    pos++;
  }
#endif

  return;
}


/* Need this procedure because we are skipping some oligomers */
static bool
find_shifted_canonical (char *genomicuc_ptr, unsigned int prevposition, unsigned int position,
			int querydistance, int genomiclength, int indexsize_nt, int *last_leftdi, int *last_rightdi,
			bool skip_repetitive_p) {
  unsigned int leftpos, rightpos;
  int shift, leftmiss, rightmiss;

  leftpos = prevposition + querydistance + indexsize_nt - 1;
  rightpos = position;

  debug7(printf("Looking for shifted canonical at prevposition %d to position %d\n",prevposition,position));
  debug7(printf("leftpos = %u, rightpos = %u\n",leftpos,rightpos));

  if (leftpos > genomiclength || rightpos > genomiclength) {
    return false;
  }

  if (skip_repetitive_p == false) {
    return (leftpos == last_leftdi[leftpos] && rightpos == last_rightdi[rightpos]);
  }

#ifdef SHIFT_EXTRA
  /* Allow canonical to be to right of match */
  leftpos += SHIFT_EXTRA;
  rightpos += SHIFT_EXTRA;
  debug7(printf("after shift, leftpos = %u, rightpos = %u\n",leftpos,rightpos));
#endif

  shift = 0;
  while (shift <= querydistance
#ifdef SHIFT_EXTRA
	 + SHIFT_EXTRA + SHIFT_EXTRA
#endif
	 ) {
    if (leftpos > rightpos) {
      return false;
    }

    if (last_leftdi[leftpos] < 0) {
      debug7(printf("\n"));
      return false;
    } else {
      leftmiss = leftpos - last_leftdi[leftpos];
    }

    if (last_rightdi[rightpos] < 0) {
      debug7(printf("\n"));
      return false;
    } else {
      rightmiss = rightpos - last_rightdi[rightpos];
    }

    debug7(printf("shift %d/left %d (miss %d)/right %d (miss %d)\n",shift,leftpos,leftmiss,rightpos,rightmiss));
    if (leftmiss == rightmiss) {  /* was leftmiss == 0 && rightmiss == 0, which doesn't allow for a shift */
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

		char *genomicuc_ptr, int *lastGT, int *lastAG,
#ifndef PMAP
		int *lastCT, int *lastAC, 
#endif
		int indexsize, Intlist_T processed, int sufflookback, int nsufflookback, int maxintronlen, 
		int genomiclength, bool localp, bool skip_repetitive_p) {
  Link_T prevlink;
  int best_fwd_consecutive = indexsize*NT_PER_MATCH;
#if 0
  int best_fwd_score = indexsize*CONSEC_POINTS_PER_MATCH, fwd_score;
#else
  int best_fwd_score = 0, fwd_score;
#endif
  int best_fwd_prevpos = -1, best_fwd_prevhit = -1;
  int best_fwd_intronnfwd = 0, best_fwd_intronnrev = 0, best_fwd_intronnunk = 0;
#ifndef PMAP
  int best_rev_consecutive = indexsize*NT_PER_MATCH;
#if 0
  int best_rev_score = indexsize*CONSEC_POINTS_PER_MATCH, rev_score;
#else
  int best_rev_score = 0, rev_score;
#endif
  int best_rev_prevpos = -1, best_rev_prevhit = -1;
  int best_rev_intronnfwd = 0, best_rev_intronnrev = 0, best_rev_intronnunk = 0;
#endif
  bool adjacentp = false, donep;
  int prev_querypos, prevhit;
  Intlist_T p;
  unsigned int prevposition, gendistance;
  int querydistance, diffdistance, lookback, nlookback, nseen, indexsize_nt;
  int canonicalsgn;
  int enough_consecutive;

#ifdef PMAP
  int intronpos_lower, intronpos_upper;
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

  enough_consecutive = 32;

  /* A. Evaluate adjacent position (at last one processed) */
  if (processed != NULL) {
    prev_querypos = Intlist_head(processed);
#ifdef PMAP
    querydistance = (querypos - prev_querypos)*3;
#else
    querydistance = querypos - prev_querypos;
#endif
    prevhit = firstactive[prev_querypos];
    while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + querydistance <= position) {
      if (prevposition + querydistance == position) {
	prevlink = &(links[prev_querypos][prevhit]);
	best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	best_fwd_score = prevlink->fwd_score + querydistance;
#ifndef PMAP
	best_rev_consecutive = prevlink->rev_consecutive + querydistance;
	best_rev_score = prevlink->rev_score + querydistance;
#endif

	best_fwd_prevpos = 
#ifndef PMAP
	  best_rev_prevpos = 
#endif
	  prev_querypos;
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
	debug9(printf("\tA. Adjacent qpos %d,%d at %ux%d (scores = %d -> %d, consec = %d, intr = %d-%d-%d)\n",
		      prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],prevlink->fwd_score,
		      best_fwd_score,best_fwd_consecutive,best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk));
#else
	debug9(printf("\tA. Adjacent qpos %d,%d at %ux%d (scores = %d/%d -> %d/%d, consec = %d/%d, intr = %d-%d-%d/%d-%d-%d)\n",
		      prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],prevlink->fwd_score,prevlink->rev_score,
		      best_fwd_score,best_rev_score,best_fwd_consecutive,best_rev_consecutive,
		      best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		      best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk));
#endif
	prevhit = -1;		/* Exit loop */
      } else {
	prevhit = active[prev_querypos][prevhit];
      }
    }
  }

#if 0
  /* AA. Evaluate querypos - 3 (for wobble) */
  if (adjacentp == false && querypos >= 3) {
    prevhit = firstactive[querypos-3];
    while (prevhit != -1 && (prevposition = mappings[querypos-3][prevhit]) + NT_PER_CODON <= position) {
      if (prevposition + NT_PER_CODON == position) {
	prevlink = &(links[querypos-3][prevhit]);
	best_fwd_consecutive = prevlink->fwd_consecutive + NT_PER_CODON;
	best_fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_CODON;
	best_rev_consecutive = prevlink->rev_consecutive + NT_PER_CODON;
	best_rev_score = prevlink->rev_score + CONSEC_POINTS_PER_CODON;

	best_fwd_prevpos = best_rev_prevpos = querypos-3;
	best_fwd_prevhit = best_rev_prevhit = prevhit;
	best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	best_fwd_intronnrev = prevlink->fwd_intronnrev;
	best_fwd_intronnunk = prevlink->fwd_intronnunk;
	best_rev_intronnfwd = prevlink->rev_intronnfwd;
	best_rev_intronnrev = prevlink->rev_intronnrev;
	best_rev_intronnunk = prevlink->rev_intronnunk;
	adjacentp = true;

	debug9(printf("\tAA. Codon-adjacent qpos %d,%d at %u (scores = %d/%d -> %d/%d, consec = %d/%d, intr = %d-%d-%d/%d-%d-%d)\n",
		     querypos-3,prevhit,position - NT_PER_CODON,prevlink->fwd_score,prevlink->rev_score,
		     best_fwd_score,best_rev_score,best_fwd_consecutive,best_rev_consecutive,
		     best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		     best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk));
	prevhit = -1; /* Exit loop */
      } else {
	prevhit = active[querypos-3][prevhit];
      }
    }
  }
#endif


  if (best_fwd_consecutive < enough_consecutive
#ifndef PMAP
      || best_rev_consecutive < enough_consecutive
#endif
      ) {

#ifdef BLUETARP
    /* It looks like case D can cover cases B and C */

    /* B. Evaluate for intron (at querypos - indexsize).  
       Don't update best_consecutive, because a reasonable alternative 
       might be a mismatch, regardless of the quality of the intron. */

#ifdef PMAP
    /* Use indexsize, not indexsize_nt */
    intronpos_lower = querypos - indexsize - 1;
    intronpos_upper = querypos - indexsize;
#else
    intronpos_lower = intronpos_upper = querypos - indexsize_nt;
#endif
    if (intronpos_lower < 0) {
      intronpos_lower = 0;
    }

    for (prev_querypos = intronpos_upper; prev_querypos >= intronpos_lower; --prev_querypos) {
      if (0 && skip_repetitive_p && nactive[prev_querypos] > MAX_NACTIVE) { /* Turn off.  If too many active, then firstactive will be -1 */
	debug9(printf("Not evaluating for intron because nactive at prev_querypos %d is %d > max %d\n",
		      prev_querypos,nactive[prev_querypos],MAX_NACTIVE));
	prevhit = -1;
      } else {
#ifdef PMAP
	querydistance = (querypos - prev_querypos)*3 - indexsize_nt
#else
	querydistance = (querypos - prev_querypos) - indexsize_nt;
#endif
	prevhit = firstactive[prev_querypos];
      }

      debug9(printf("Evaluating for intron at prev_querypos %d, prevhit %d\n",prev_querypos,prevhit));
      while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
	prevlink = &(links[prev_querypos][prevhit]);
	if (prevlink->fwd_consecutive > EXON_DEFN
#ifndef PMAP
	    || prevlink->rev_consecutive > EXON_DEFN
#endif
	    ) {
	  gendistance = position - prevposition - indexsize_nt;
	  diffdistance = abs(gendistance - querydistance);
	
	  if (diffdistance > maxintronlen) {
	    canonicalsgn = 0;
	  } else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,genomiclength,indexsize_nt,
					    lastGT,lastAG,skip_repetitive_p) == true) {
	    canonicalsgn = +1;
#ifndef PMAP
	  } else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,genomiclength,indexsize_nt,
					    lastCT,lastAC,skip_repetitive_p) == true) {
	    canonicalsgn = -1;
#endif
	  } else {
	    canonicalsgn = 0;
	  }

	  /* Need to both compensate for the oligomer lost across
	     the intron and penalize for each additional intron */
#ifdef PMAP
	  fwd_score = prevlink->fwd_score + querydistance + indexsize_nt;
#else
	  fwd_score = prevlink->fwd_score + querydistance + indexsize_nt;
	  rev_score = prevlink->rev_score + querydistance + indexsize_nt;
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
		gendist_penalty(diffdistance,querydistance/*-indexsize_nt*/) + INTRON_PENALTY_UNKNOWN;
	    } else if (canonicalsgn == +1) {
	      /* Don't penalize for intron */
	      fwd_gendistance_penalty = gendist_penalty(diffdistance,querydistance/*-indexsize_nt*/) /* - INTRON_REWARD_CONSISTENT*/;
#ifndef PMAP
	      rev_gendistance_penalty = gendist_penalty(diffdistance,querydistance/*-indexsize_nt*/) + INTRON_PENALTY_INCONSISTENT;
#endif
#ifndef PMAP
	    } else if (canonicalsgn == -1) {
	      fwd_gendistance_penalty = gendist_penalty(diffdistance,querydistance) + INTRON_PENALTY_INCONSISTENT;
	      /* Don't penalize for intron */
	      rev_gendistance_penalty = gendist_penalty(diffdistance,querydistance) /* - INTRON_REWARD_CONSISTENT*/;
#endif
	    }
	  }
	  fwd_score -= fwd_gendistance_penalty;
#ifndef PMAP
	  rev_score -= rev_gendistance_penalty;
#endif

	  if (fwd_score > best_fwd_score) {
	    best_fwd_consecutive = indexsize_nt;
	    /* best_fwd_score = fwd_score; -- Put below so debug9 works */
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
	  }

#ifndef PMAP
	  if (rev_score > best_rev_score) {
	    best_rev_consecutive = indexsize_nt;
	    /* best_rev_score = rev_score; -- Put below so debug9 works */
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
	  }
#endif

#ifdef PMAP
	  debug9(printf("\tB. Intron qpos %d,%d at %ux%d (scores = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, pen %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance,fwd_gendistance_penalty));
#else
	  debug9(printf("\tB. Intron qpos %d,%d at %ux%d (scores = %d/%d -> %d/%d, consec = %d/%d, intr = %d-%d-%d/%d-%d-%d, gendist %u, pen %d/%d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,prevlink->rev_score,fwd_score,rev_score,
			prevlink->fwd_consecutive,prevlink->rev_consecutive,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk,
			gendistance,fwd_gendistance_penalty,rev_gendistance_penalty));
#endif
	  debug9({
	    printf(" =>");
	    if (fwd_score > best_fwd_score) {
	      printf(" Best fwd at %d",fwd_score);
	    } else {
	      printf(" Fwd loses to %d",best_fwd_score);
	    }
#ifndef PMAP
	    if (rev_score > best_rev_score) {
	      printf(" Best rev at %d",rev_score);
	    } else {
	      printf(" Rev loses to %d\n",best_rev_score);
	    }
#endif
	    printf("\n");
	  });
	}

	/* Put these down here so debug9 works */
	if (fwd_score > best_fwd_score) {
	  best_fwd_score = fwd_score;
	}
#ifndef PMAP
	if (rev_score > best_rev_score) {
	  best_rev_score = rev_score;
	}
#endif
	prevhit = active[prev_querypos][prevhit];
      }
    }

#endif /* bluetarp */


#ifdef BLUETARP
    /* It looks like case D can cover cases B and C */

    /* C (fwd). Look back at grand querypos and hit */
    debug9(printf("\tGrand fwd querypos is %d,%d at position %u\n",
		  grand_fwd_querypos,grand_fwd_hit,grand_fwd_hit >= 0 ? mappings[grand_fwd_querypos][grand_fwd_hit] : 0));
    if (grand_fwd_querypos > 0 && grand_fwd_querypos + indexsize_nt <= querypos &&
	(prevposition = mappings[grand_fwd_querypos][grand_fwd_hit]) + indexsize_nt <= position) {
      prev_querypos = grand_fwd_querypos;
      prevhit = grand_fwd_hit;
#ifdef PMAP
      querydistance = (querypos - prev_querypos) - indexsize_nt;
#else
      querydistance = querypos - prev_querypos - indexsize_nt;
#endif
      if (querydistance < MAX_GRAND_LOOKBACK) {
	prevlink = &(links[prev_querypos][prevhit]);
	if (1 || prevlink->fwd_consecutive > EXON_DEFN) {
	  gendistance = position - prevposition - indexsize_nt;
	  diffdistance = abs(gendistance - querydistance);

	  if (diffdistance <= EQUAL_DISTANCE) {
	    canonicalsgn = 9;
	  } else if (diffdistance > maxintronlen) {
	    canonicalsgn = 0;
	  } else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,genomiclength,indexsize_nt,
					    lastGT,lastAG,skip_repetitive_p) == true) {
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
	    /* fwd_score -= querydist_mismatch_penalty(querydistance); */
#ifdef PMAP
	    if (diffdistance % 3 != 0) {
	      fwd_score -= NONCODON_INDEL_PENALTY;
	    }
#endif
	  } else if (canonicalsgn == +1) {
	    /* Don't penalize for intron.  Don't penalize for querydistance to grand querypos. */
	    /* fwd_score -= gendist_penalty(diffdistance,querydistance); */
	    fwd_score -= diffdist_penalty(diffdistance);

	  } else {
	    /* fwd_score -= querydist_mismatch_penalty(querydistance); */
	    /* fwd_gendistance_penalty = gendist_penalty(diffdistance,querydistance); */
	    fwd_gendistance_penalty = diffdist_penalty(diffdistance);
	    fwd_score -= fwd_gendistance_penalty;
	    fwd_score -= NINTRON_PENALTY_MISMATCH;
	  }

	  debug9(printf("\tC+. Fwd grand qpos %d,%d at %ux%d (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance,querydistance,canonicalsgn));
	    
	  if (fwd_score > best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + NT_PER_MATCH;
	    } else {
	      best_fwd_consecutive = indexsize_nt;
	    }
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
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }
	}
      }
    }

#ifndef PMAP
    /* C (rev). Look back at grand querypos and hit */
    debug9(printf("\tGrand rev querypos is %d,%d at position %u\n",
		  grand_rev_querypos,grand_rev_hit,grand_rev_hit >= 0 ? mappings[grand_rev_querypos][grand_rev_hit] : 0));
    if (grand_rev_querypos > 0 && grand_rev_querypos + indexsize_nt <= querypos &&
	(prevposition = mappings[grand_rev_querypos][grand_rev_hit]) + indexsize_nt <= position) {
      prev_querypos = grand_rev_querypos;
      prevhit = grand_rev_hit;
#ifdef PMAP
      querydistance = (querypos - prev_querypos)*3 - indexsize_nt;
#else
      querydistance = querypos - prev_querypos - indexsize_nt;
#endif
      if (querydistance < MAX_GRAND_LOOKBACK) {
	prevlink = &(links[prev_querypos][prevhit]);
	if (1 || prevlink->rev_consecutive > EXON_DEFN) {
	  gendistance = position - prevposition - indexsize_nt;
	  diffdistance = abs(gendistance - querydistance);

	  if (diffdistance <= EQUAL_DISTANCE) {
	    canonicalsgn = 9;
	  } else if (diffdistance > maxintronlen) {
	    canonicalsgn = 0;
	  } else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,genomiclength,indexsize_nt,
					    lastCT,lastAC,skip_repetitive_p) == true) {
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
	    /* rev_score -= querydist_mismatch_penalty(querydistance); */
	  } else if (canonicalsgn == -1) {
	    /* Don't penalize for intron.  Don't penalize for querydistance back to grand querypos. */
	    /* rev_score -= gendist_penalty(diffdistance,querydistance); */
	    rev_score -= diffdist_penalty(diffdistance);

	  } else {
	    /* rev_score -= querydist_mismatch_penalty(querydistance); */
	    /* rev_gendistance_penalty = gendist_penalty(diffdistance,querydistance); */
	    rev_gendistance_penalty = diffdist_penalty(diffdistance);
	    rev_score -= rev_gendistance_penalty;
	    rev_score -= NINTRON_PENALTY_MISMATCH;
	  }

	  debug9(printf("\tC-. Rev grand qpos %d,%d at %ux%d (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->rev_score,rev_score,prevlink->rev_consecutive,
			best_rev_intronnrev,best_rev_intronnrev,best_rev_intronnunk,
			gendistance,querydistance,canonicalsgn));
	    
	  if (rev_score > best_rev_score) {
	    if (diffdistance <= EQUAL_DISTANCE) {
	      best_rev_consecutive = prevlink->rev_consecutive + NT_PER_MATCH;
	    } else {
	      best_rev_consecutive = indexsize_nt;
	    }
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
	    debug9(printf(" => Best rev at %d (consec = %d)\n",rev_score,best_rev_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_rev_score));
	  }
	}
      }
    }
#endif

#endif /* bluetarp */


    /* D (fwd). Evaluate for mismatches (all other previous querypos) */
    /* Set parameters */
    if (adjacentp == true) {
      /* Look not so far back */
      nlookback = 1;
      lookback = sufflookback/2;
    } else {
      /* Look farther back */
      nlookback = nsufflookback;
      lookback = sufflookback;
    }

    donep = false;
    for (p = processed, nseen = 0;
	 p != NULL && best_fwd_consecutive < enough_consecutive && donep == false;
	 p = Intlist_next(p), nseen++) {
				  
      prev_querypos = Intlist_head(p);
#ifdef PMAP
      querydistance = (querypos - prev_querypos)*3 - indexsize_nt;
#else
      querydistance = querypos - prev_querypos - indexsize_nt;
#endif

      if (nseen > nlookback && querydistance > lookback) {
	donep = true;
      }

      prevhit = firstactive[prev_querypos];
      while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
	prevlink = &(links[prev_querypos][prevhit]);

	gendistance = position - prevposition - indexsize_nt;
	diffdistance = abs(gendistance - querydistance);

	if (diffdistance < maxintronlen) {
	  if (diffdistance <= EQUAL_DISTANCE) {
	    canonicalsgn = 9;
	    fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_MATCH;
#ifdef PMAP
	    if (diffdistance % 3 != 0) {
	      fwd_score -= NONCODON_INDEL_PENALTY;
	    }
#endif
	  } else if (prevlink->fwd_consecutive < EXON_DEFN) {
	    canonicalsgn = 0;
	    fwd_score = prevlink->fwd_score - diffdist_penalty(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;

	  } else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,genomiclength,indexsize_nt,
					    lastGT,lastAG,skip_repetitive_p) == true) {
	    canonicalsgn = +1;
	    fwd_score = prevlink->fwd_score - diffdist_penalty(diffdistance) - querydist_penalty(querydistance);
	  } else {
	    canonicalsgn = 0;
	    fwd_score = prevlink->fwd_score - diffdist_penalty(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;
	  }

	  debug9(printf("\tD+. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance,querydistance,canonicalsgn));
	    
	  if (fwd_score > best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + (querydistance + indexsize_nt);
	    } else {
	      best_fwd_consecutive = 0;
	    }
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
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }
	}

	prevhit = active[prev_querypos][prevhit];
      }
    }

#ifndef PMAP
    /* D (rev). Evaluate for mismatches (all other previous querypos) */
    /* Parameters already set in D (fwd). */

    donep = false;
    for (p = processed, nseen = 0;
	 p != NULL && best_rev_consecutive < enough_consecutive && donep == false;
	 p = Intlist_next(p), nseen++) {

      prev_querypos = Intlist_head(p);
#ifdef PMAP
      querydistance = (querypos - prev_querypos)*3 - indexsize_nt;
#else
      querydistance = querypos - prev_querypos - indexsize_nt;
#endif

      if (nseen > nlookback && querydistance > lookback) {
	donep = true;
      }

      prevhit = firstactive[prev_querypos];
      while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
	prevlink = &(links[prev_querypos][prevhit]);

	gendistance = position - prevposition - indexsize_nt;
	diffdistance = abs(gendistance - querydistance);

	if (diffdistance < maxintronlen) {
	  if (diffdistance <= EQUAL_DISTANCE) {
	    canonicalsgn = 9;
	    rev_score = prevlink->rev_score + CONSEC_POINTS_PER_MATCH;

	  } else if (prevlink->rev_consecutive < EXON_DEFN) {
	    canonicalsgn = 0;
	    rev_score = prevlink->rev_score - diffdist_penalty(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;

	  } else if (find_shifted_canonical(genomicuc_ptr,prevposition,position,querydistance,genomiclength,indexsize_nt,
					    lastCT,lastAC,skip_repetitive_p) == true) {
	    canonicalsgn = -1;
	    rev_score = prevlink->rev_score - diffdist_penalty(diffdistance) - querydist_penalty(querydistance);
	  } else {
	    canonicalsgn = 0;
	    rev_score = prevlink->rev_score - diffdist_penalty(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;
	  }
	  
	  debug9(printf("\tD-. Rev mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->rev_score,rev_score,prevlink->rev_consecutive,
			best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk,
			gendistance,querydistance,canonicalsgn));
	    
	  if (rev_score > best_rev_score) {
	    if (diffdistance <= EQUAL_DISTANCE) {
	      best_rev_consecutive = prevlink->rev_consecutive + (querydistance + indexsize_nt);
	    } else {
	      best_rev_consecutive = 0;
	    }
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
	    debug9(printf(" => Best rev at %d (consec = %d)\n",rev_score,best_rev_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_rev_score));
	  }
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
  if (localp && currlink->fwd_pos < 0) {
    currlink->fwd_score = indexsize_nt;
  } else {
    currlink->fwd_score = best_fwd_score;
  }
  currlink->fwd_intronnfwd = best_fwd_intronnfwd;
  currlink->fwd_intronnrev = best_fwd_intronnrev;
  currlink->fwd_intronnunk = best_fwd_intronnunk;

#ifndef PMAP
  /* linksconsecutive already assigned above */
  currlink->rev_consecutive = best_rev_consecutive;
  currlink->rev_pos = best_rev_prevpos;
  currlink->rev_hit = best_rev_prevhit;
  if (localp && currlink->rev_pos < 0) {
    currlink->rev_score = indexsize_nt;
  } else {
    currlink->rev_score = best_rev_score;
  }
  currlink->rev_intronnfwd = best_rev_intronnfwd;
  currlink->rev_intronnrev = best_rev_intronnrev;
  currlink->rev_intronnunk = best_rev_intronnunk;
#endif

#ifdef PMAP
  debug9(printf("\tChose %d,%d with score %d (fwd)\n",
	       currlink->fwd_pos,currlink->fwd_hit,currlink->fwd_score));
#else
  debug9(printf("\tChose %d,%d with score %d (fwd) and %d,%d with score %d (rev)\n",
	       currlink->fwd_pos,currlink->fwd_hit,currlink->fwd_score,currlink->rev_pos,currlink->rev_hit,currlink->rev_score));
#endif
  debug3(printf("%d %d  %d %d  1\n",querypos,hit,best_prevpos,best_prevhit));
  return;
}


static void
revise_active (int **active, int *firstactive, int *nactive, 
	       int low_hit, int high_hit, struct Link_T **links, int querypos) {
  int best_score, threshold, score;
  int hit, *ptr;

  debug6(printf("Revising querypos %d from low_hit %d to high_hit %d.  Scores:\n",querypos,low_hit,high_hit));
  if ((hit = low_hit) >= high_hit) {
    firstactive[querypos] = -1;
  } else {
    debug6(printf("At hit %d, fwd_score is %d",hit,links[querypos][hit].fwd_score));
    best_score = links[querypos][hit].fwd_score;
#ifndef PMAP
    debug6(printf(" and rev_score is %d",links[querypos][hit].rev_score));
    if ((score = links[querypos][hit].rev_score) > best_score) {
      best_score = score;
    }
#endif
    debug6(printf("\n"));

    for (hit++; hit < high_hit; hit++) {
      debug6(printf("At hit %d, fwd_score is %d",hit,links[querypos][hit].fwd_score));
      if ((score = links[querypos][hit].fwd_score) > best_score) {
	best_score = score;
      }
#ifndef PMAP
      debug6(printf(" and rev_score is %d",links[querypos][hit].rev_score));
      if ((score = links[querypos][hit].rev_score) > best_score) {
	best_score = score;
      }
#endif
      debug6(printf("\n"));
    }

    threshold = best_score - SCORE_FOR_RESTRICT;
    if (threshold < 0) {
      threshold = 0;
    }

    nactive[querypos] = 0;
    ptr = &(firstactive[querypos]);
    hit = low_hit;
    while (hit < high_hit) {
      while (hit < high_hit && links[querypos][hit].fwd_score <= threshold
#ifndef PMAP
	     && links[querypos][hit].rev_score <= threshold
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
intmatrix_1d_new (int length1, int *lengths2, int totallength) {
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

static void
intmatrix_1d_free (int ***matrix) {
  FREE((*matrix)[0]);
  FREE(*matrix);
  return;
}


static int **
intmatrix_2d_new (int length1, int *lengths2) {
  int **matrix;
  int i;
  
  matrix = (int **) CALLOC(length1,sizeof(int *));
  for (i = 0; i < length1; i++) {
    if (lengths2[i] <= 0) {
      matrix[i] = (int *) NULL;
    } else {
      matrix[i] = (int *) CALLOC(lengths2[i],sizeof(int));
    }
  }
  return matrix;
}

static void
intmatrix_2d_free (int ***matrix, int length1) {
  int i;

  for (i = 0; i < length1; i++) {
    if ((*matrix)[i]) {
      FREE((*matrix)[i]);
    }
  }
  FREE(*matrix);
  return;
}


static bool
zero_score_p (struct Link_T *links, int npositions) {
  int hit;
  Link_T currlink;

  for (hit = 0; hit < npositions; hit++) {
    currlink = &(links[hit]);
    if (currlink->fwd_score > 0) {
      return false;
    }
#ifndef PMAP
    if (currlink->rev_score > 0) {
      return false;
    }
#endif
  }
  return true;
}


/* For PMAP, indexsize is in aa. */
static bool
align_compute_scores (bool *fwdp, int *grand_querypos, int *grand_hit,
		      struct Link_T **links, unsigned int **mappings, int *npositions, int totalpositions,
		      bool oned_matrix_p, unsigned int *minactive, unsigned int *maxactive,
		      int querystart, int queryend, int querylength,
		      char *genomicuc_ptr, int genomiclength, int indexsize,
		      int sufflookback, int nsufflookback, int maxintronlen,
		      char *queryseq_ptr, bool localp, bool skip_repetitive_p, bool debug_graphic_p) {
  bool result;
  Link_T currlink, prevlink;
  int querypos, indexsize_nt, hit, low_hit, high_hit;
  int nskipped, min_hits, specific_querypos, specific_low_hit, specific_high_hit, next_querypos;
  Intlist_T processed = NULL;
  int grand_fwd_score, grand_fwd_querypos, grand_fwd_hit, best_fwd_hit, best_fwd_score;
  int best_score;
#ifndef PMAP
  int grand_rev_score, grand_rev_querypos, grand_rev_hit, best_rev_hit, best_rev_score;
#endif
  int **active, *firstactive, *nactive;
  unsigned int position, prevposition;
  int *lastGT, *lastAG;
#ifndef PMAP
  int *lastCT, *lastAC;
#endif
#ifdef DEBUG9
  Link_T termlink = NULL;
  char *oligo;
#endif

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

#ifdef DEBUG9
  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
#endif
  debug0(printf("Querystart = %d, queryend = %d, indexsize = %d\n",querystart,queryend,indexsize));

  if (oned_matrix_p == true) {
    active = intmatrix_1d_new(querylength,npositions,totalpositions);
  } else {
    active = intmatrix_2d_new(querylength,npositions);
  }

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

#ifdef SHIFT_EXTRA
#ifdef PMAP
  lastGT = (int *) CALLOC(genomiclength+SHIFT_EXTRA+1,sizeof(int));
  lastAG = (int *) CALLOC(genomiclength+SHIFT_EXTRA+1,sizeof(int));
  find_canonical_dinucleotides(lastGT,lastAG,genomicuc_ptr,genomiclength);
#else
  lastGT = (int *) CALLOC(genomiclength+SHIFT_EXTRA+1,sizeof(int));
  lastAG = (int *) CALLOC(genomiclength+SHIFT_EXTRA+1,sizeof(int));
  lastCT = (int *) CALLOC(genomiclength+SHIFT_EXTRA+1,sizeof(int));
  lastAC = (int *) CALLOC(genomiclength+SHIFT_EXTRA+1,sizeof(int));
  find_canonical_dinucleotides(lastGT,lastAG,lastCT,lastAC,genomicuc_ptr,genomiclength);
#endif
#else
#ifdef PMAP
  lastGT = (int *) CALLOC(genomiclength+1,sizeof(int));
  lastAG = (int *) CALLOC(genomiclength+1,sizeof(int));
  find_canonical_dinucleotides(lastGT,lastAG,genomicuc_ptr,genomiclength);
#else
  lastGT = (int *) CALLOC(genomiclength+1,sizeof(int));
  lastAG = (int *) CALLOC(genomiclength+1,sizeof(int));
  lastCT = (int *) CALLOC(genomiclength+1,sizeof(int));
  lastAC = (int *) CALLOC(genomiclength+1,sizeof(int));
  find_canonical_dinucleotides(lastGT,lastAG,lastCT,lastAC,genomicuc_ptr,genomiclength);
#endif
#endif

  grand_fwd_score = 0;
  grand_fwd_querypos = -1;
  grand_fwd_hit = -1;
#ifndef PMAP
  grand_rev_score = 0;
  grand_rev_querypos = -1;
  grand_rev_hit = -1;
#endif

  nskipped = 0;
  min_hits = 1000000;
  specific_querypos = -1;

  querypos += 1;
  while (querypos <= queryend - indexsize) {
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
    debug9(printf("Querypos %d has hit %d..%d out of %d (minactive = %u, maxactive = %u)\n",
		 querypos,low_hit,high_hit-1,npositions[querypos],minactive[querypos],maxactive[querypos]));

    /* Can't use nactive yet, so use high_hit - low_hit */
    if (skip_repetitive_p && high_hit - low_hit >= MAX_NACTIVE && nskipped <= MAX_SKIPPED) { /* Previously turned off */
      debug6(printf("Too many active (%d - %d) at querypos %d.  Setting firstactive to be -1\n",high_hit,low_hit,querypos));
      firstactive[querypos] = -1;
      nskipped++;
      debug9(printf("  %d skipped because of %d hits\n",nskipped,high_hit - low_hit + 1));

      /* Store most specific querypos in section of skipped */
      if (high_hit - low_hit < min_hits) {
	min_hits = high_hit - low_hit;
	specific_querypos = querypos;
	specific_low_hit = low_hit;
	specific_high_hit = high_hit;
      }
      querypos++;

    } else {
      if (nskipped > MAX_SKIPPED) {
	debug9(printf("Too many skipped.  Going back to specific querypos %d\n",specific_querypos));
	next_querypos = querypos;
	querypos = specific_querypos;
	low_hit = specific_low_hit;
	high_hit = specific_high_hit;
      } else {
	next_querypos = querypos + 1;
      }

      for (hit = low_hit; hit < high_hit; hit++) {
	currlink = &(links[querypos][hit]);
	position = mappings[querypos][hit];

	debug9(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
	debug9(printf("Finding link at querypos %d,%d at %ux%d (%s).  prev_querypos was %d\n",
		      querypos,hit,position,active[querypos][hit],oligo,processed ? Intlist_head(processed) : -1));

	score_querypos(currlink,querypos,position,links,mappings,npositions,active,firstactive,nactive,
		       grand_fwd_querypos,grand_fwd_hit,
#ifndef PMAP
		       grand_rev_querypos,grand_rev_hit,
#endif
		       genomicuc_ptr,lastGT,lastAG,
#ifndef PMAP
		       lastCT,lastAC,
#endif
		       indexsize,processed,sufflookback,nsufflookback,maxintronlen,genomiclength,
		       localp,skip_repetitive_p);

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
      nskipped = 0;
      min_hits = 1000000;
      specific_querypos = -1;
      
#ifdef PMAP
      debug9(printf("Overall result at querypos %d yields best_fwd_hit %d\n",
		   querypos,best_fwd_hit));
#else
      debug9(printf("Overall result at querypos %d yields best_fwd_hit %d and best_rev_hit %d\n",
		   querypos,best_fwd_hit,best_rev_hit));
#endif

      if (best_fwd_hit >= 0 && links[querypos][best_fwd_hit].fwd_hit < 0 && 
	  grand_fwd_querypos >= 0 && querypos >= grand_fwd_querypos + indexsize_nt) {
	prevlink = &(links[grand_fwd_querypos][grand_fwd_hit]);
	if ((best_fwd_score = prevlink->fwd_score - (querypos - grand_fwd_querypos)) > 0) {
	  prevposition = mappings[grand_fwd_querypos][grand_fwd_hit];
	  for (hit = low_hit; hit < high_hit; hit++) {
	    currlink = &(links[querypos][hit]);
	    if ((position = mappings[querypos][hit]) >= prevposition + indexsize_nt) {
	      currlink->fwd_consecutive = indexsize_nt;
	      currlink->fwd_pos = grand_fwd_querypos;
	      currlink->fwd_hit = grand_fwd_hit;
	      currlink->fwd_score = best_fwd_score;
	      currlink->fwd_intronnfwd = prevlink->fwd_intronnfwd;
	      currlink->fwd_intronnrev = prevlink->fwd_intronnrev;
	      currlink->fwd_intronnunk = prevlink->fwd_intronnunk + 1;
	    }
	  }
	  debug9(printf("At querypos %d, setting all fwd hits to point back to grand_fwd %d,%d with a score of %d\n",
		       querypos,grand_fwd_querypos,grand_fwd_hit,prevlink->fwd_score));
	}
      }

      /* Use >= to favor longer path in case of ties */
      if (best_fwd_hit >= 0 && best_fwd_score >= grand_fwd_score && 
	  links[querypos][best_fwd_hit].fwd_consecutive > EXON_DEFN) {
	grand_fwd_score = best_fwd_score;
	grand_fwd_querypos = querypos;
	grand_fwd_hit = best_fwd_hit;
	debug9(termlink = &(links[querypos][best_fwd_hit]));
	debug9(printf("At querypos %d, revising grand fwd to be hit %d with score of %d (pointing back to %d,%d)\n",
		     querypos,best_fwd_hit,best_fwd_score,termlink->fwd_pos,termlink->fwd_hit));
      }

#ifndef PMAP
      if (best_rev_hit >= 0 && links[querypos][best_rev_hit].rev_hit < 0 && 
	  grand_rev_querypos >= 0 && querypos >= grand_rev_querypos + indexsize_nt) {
	prevlink = &(links[grand_rev_querypos][grand_rev_hit]);
	if ((best_rev_score = prevlink->rev_score - (querypos - grand_rev_querypos)) > 0) {
	  prevposition = mappings[grand_rev_querypos][grand_rev_hit];
	  for (hit = low_hit; hit < high_hit; hit++) {
	    currlink = &(links[querypos][hit]);
	    if ((position = mappings[querypos][hit]) >= prevposition + indexsize_nt) {
	      currlink->rev_consecutive = indexsize_nt;
	      currlink->rev_pos = grand_rev_querypos;
	      currlink->rev_hit = grand_rev_hit;
	      currlink->rev_score = best_rev_score;
	      currlink->rev_intronnrev = prevlink->rev_intronnfwd;
	      currlink->rev_intronnrev = prevlink->rev_intronnrev;
	      currlink->rev_intronnunk = prevlink->rev_intronnunk + 1;
	    }
	  }
	  debug9(printf("At querypos %d, setting all rev hits to point back to grand_rev %d,%d with a score of %d\n",
		       querypos,grand_rev_querypos,grand_rev_hit,prevlink->rev_score));
	}
      }

      /* Use >= to favor longer path in case of ties */
      if (best_rev_hit >= 0 && best_rev_score >= grand_rev_score &&
	  links[querypos][best_rev_hit].rev_consecutive > EXON_DEFN) {
	grand_rev_score = best_rev_score;
	grand_rev_querypos = querypos;
	grand_rev_hit = best_rev_hit;
      }
#endif

      revise_active(active,firstactive,nactive,low_hit,high_hit,links,querypos);
      processed = Intlist_push(processed,querypos);
      querypos = next_querypos;
    }
  }

  Intlist_free(&processed);

#ifndef PMAP
  FREE(lastAC);
  FREE(lastCT);
#endif
  FREE(lastAG);
  FREE(lastGT);

  /* These are the final active oligomers, after pruning by score */
  if (debug_graphic_p == true) {
    mappings_dump_R(mappings,npositions,querylength,active,firstactive,indexsize,"active.mers");
  }

  FREE(nactive);
  FREE(firstactive);

  if (oned_matrix_p == true) {
    intmatrix_1d_free(&active);
  } else {
    intmatrix_2d_free(&active,querylength);
  }


  /* Grand winner */
  debug10(printf("Finding grand winner, using local/global method\n"));
  /* Local/global */
  querypos = queryend - indexsize;
  while (querypos >= 0 && (npositions[querypos] <= 0 || 
			   zero_score_p(links[querypos],npositions[querypos]) == true)) {
    querypos--;
  }
  if (querypos < 0) {
    result = false;
  } else {
    result = true;

    best_score = -1000000;
    for (hit = 0; hit < npositions[querypos]; hit++) {
      currlink = &(links[querypos][hit]);
      if (currlink->fwd_score > best_score) {
	*grand_querypos = querypos;
	*grand_hit = hit;
	*fwdp = true;
	best_score = currlink->fwd_score;
      }
#ifndef PMAP
      if (currlink->rev_score > best_score) {
	*grand_querypos = querypos;
	*grand_hit = hit;
	*fwdp = false;
	best_score = currlink->rev_score;
      }
#endif
    }
    
    if (grand_fwd_score > best_score) {
      best_score = grand_fwd_score;
      *grand_querypos = grand_fwd_querypos;
      *grand_hit = grand_fwd_hit;
      *fwdp = true;
    }
    
#ifndef PMAP
    if (grand_rev_score > best_score) {
      best_score = grand_rev_score;
      *grand_querypos = grand_rev_querypos;
      *grand_hit = grand_rev_hit;
      *fwdp = false;
    }
#endif
  }


#if 0
  /* Global */
  querypos = queryend - indexsize;
  while (querypos >= 0 && npositions[querypos] <= 0) {
    querypos--;
  }

  if (querypos < 0) {
    result = false;
  } else {
    result = true;
    best_score = -1000000;
    for (hit = 0; hit < npositions[querypos]; hit++) {
      currlink = &(links[querypos][hit]);
      if (currlink->fwd_score > best_score) {
	*fwdp = true;
	best_score = currlink->fwd_score;
	*grand_hit = hit;
	debug10(printf("==> Found best score (fwd) %d at %d,%d\n",best_score,querypos,hit));
      }
#ifndef PMAP
      if (currlink->rev_score > best_score) {
	*fwdp = false;
	best_score = currlink->rev_score;
	*grand_hit = hit;
	debug10(printf("==> Found best score (rev) %d at %d,%d\n",best_score,querypos,hit));
      }
#endif
    }
    *grand_querypos = querypos;
  }
#endif


#if 0
  /* Local method.  This is a bad use of localp; misses short exons at 3' end. */
  debug10(printf("Finding grand winner, using local method\n"));
#ifdef PMAP
  if (grand_fwd_querypos < 0) {
    result = false;
  } else {
    result = true;
    *fwdp = true;
    *grand_querypos = grand_fwd_querypos;
    *grand_hit = grand_fwd_hit;
    debug10(printf("==> Found best score %d at %d,%d\n",grand_fwd_score,grand_fwd_querypos,grand_fwd_hit));
  }
#else
  if (grand_fwd_querypos < 0 && grand_rev_querypos < 0) {
    result = false;
  } else {
    result = true;
    debug10(printf("==> Found best score (fwd) %d at %d,%d\n",grand_fwd_score,grand_fwd_querypos,grand_fwd_hit));
    debug10(printf("==> Found best score (rev) %d at %d,%d\n",grand_rev_score,grand_rev_querypos,grand_rev_hit));
    if (grand_fwd_score >= grand_rev_score) {
      debug10(printf("Fwd wins\n"));
      *fwdp = true;
      *grand_querypos = grand_fwd_querypos;
      *grand_hit = grand_fwd_hit;
    } else {
      debug10(printf("Rev wins\n"));
      *fwdp = false;
      *grand_querypos = grand_rev_querypos;
      *grand_hit = grand_rev_hit;
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


#endif	/* Local method */


  debug9(FREE(oligo));

  return result;
}



/* Performs dynamic programming.  For PMAP, indexsize is in aa. */
static List_T
align_compute (int *start_querypos, int *end_querypos,
	       int *nfwdintrons, int *nrevintrons, int *nunkintrons, int *cdna_direction,
	       unsigned int **mappings, int *npositions, int totalpositions, bool oned_matrix_p,
	       unsigned int *minactive, unsigned int *maxactive,
	       char *queryseq_ptr, int querylength, int queryseq_trim_start, int queryseq_trim_end,
	       char *genomicseg_ptr, char *genomicuc_ptr, int genomiclength, int indexsize,
	       int sufflookback, int nsufflookback, int maxintronlen, Pairpool_T pairpool,
	       bool localp, bool skip_repetitive_p, bool debug_graphic_p) {
  List_T path = NULL;
  struct Link_T **links;
  Link_T termlink;
  bool fwdp;
  int querypos, prev_querypos, hit, prevhit;
  unsigned int position;
  int querystart, queryend;
  debug0(char *oligo);

  debug0(oligo = (char *) CALLOC(indexsize+1,sizeof(char)));

  querystart = queryseq_trim_start;
  queryend = queryseq_trim_end;

  if (oned_matrix_p == true) {
    links = Linkmatrix_1d_new(querylength,npositions,totalpositions);
  } else {
    links = Linkmatrix_2d_new(querylength,npositions);
  }

  /* These are all oligomers */
  if (debug_graphic_p == true) {
    mappings_dump_R(mappings,npositions,querylength,/*active*/NULL,/*firstactive*/NULL,indexsize,"all.mers");
  }
  
  if (align_compute_scores(&fwdp,&querypos,&hit,links,mappings,npositions,totalpositions,oned_matrix_p,
			   minactive,maxactive,querystart,queryend,querylength,genomicuc_ptr,genomiclength,indexsize,
			   sufflookback,nsufflookback,maxintronlen,queryseq_ptr,
			   localp,skip_repetitive_p,debug_graphic_p) == false) {
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
      printf("plot(all.mers,col=\"black\",pch=\".\",xlab=\"Query\",ylab=\"Genomic\")\n");
      printf("points(active.mers,col=\"red\",pch=\".\")\n");
      printf("points(best.path,col=\"green\",pch=\".\")\n");
      printf("lines(querypos,minactive,col=\"blue\")\n");
      printf("lines(querypos,maxactive,col=\"blue\")\n");
    }

    /* prevposition = mappings[querypos][hit]+1; */
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

      /* prevposition = position; */
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

    if (fwdp == true) {
      *cdna_direction = +1;
    } else {
      *cdna_direction = -1;
    }
  }

  if (oned_matrix_p == true) {
    Linkmatrix_1d_free(&links);
  } else {
    Linkmatrix_2d_free(&links,querylength);
  }

  debug0(FREE(oligo));

  return path;			/* previously was List_reverse(path) */
}



/* queryseq_ptr is NULL for PMAP.  querypos here is in nt. */
static List_T
convert_to_nucleotides (List_T path,
#ifndef PMAP
			char *queryseq_ptr, char *queryuc_ptr, 
#endif
			char *genomicseg_ptr, char *genomicuc_ptr,
			int query_offset, int genomic_offset,
			Pairpool_T pairpool, int indexsize_nt) {
  List_T pairs = NULL, pairptr;
  Pair_T pair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump, fill, default_fill;

  pairptr = path;
  path = Pairpool_pop(path,&pair);
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
    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos+genomic_offset,
			  genomicuc_ptr[lastgenomepos],MATCH_COMP,
			  genomicseg_ptr[lastgenomepos],/*dynprogindex*/0);
    debug5(printf("Pushing %c | %c at %d,%d\n",genomicuc_ptr[lastgenomepos],genomicseg_ptr[lastgenomepos],
		  lastquerypos,lastgenomepos));
#else
    if (queryuc_ptr[lastquerypos] == genomicuc_ptr[lastgenomepos]) {
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos+genomic_offset,
			    queryseq_ptr[lastquerypos],MATCH_COMP,
			    genomicseg_ptr[lastgenomepos],/*dynprogindex*/0);
      debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
		    lastquerypos+query_offset,lastgenomepos+genomic_offset));
    } else {
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos+genomic_offset,
			    queryseq_ptr[lastquerypos],MISMATCH_COMP,
			    genomicseg_ptr[lastgenomepos],/*dynprogindex*/0);
      debug5(printf("Pushing %c   %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
		    lastquerypos+query_offset,lastgenomepos+genomic_offset));
    }
#endif
    --lastquerypos;
    --lastgenomepos;
  }

  /* Take care of first pair */
  pair->querypos += query_offset; /* Revise coordinates */
  pair->genomepos += genomic_offset; /* Revise coordinates */
#ifdef WASTE
  pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
  pairs = List_push_existing(pairs,pairptr);
#endif
  lastquerypos = querypos;
  lastgenomepos = genomepos;

  while (path != NULL) {
    pairptr = path;
    path = Pairpool_pop(path,&pair);
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
#if 0
	  /* This can occur with wobble mask */
	  fprintf(stderr,"Partial fill from querypos %d to %d (genomepos goes from %u to %u)\n",
		  querypos,lastquerypos,genomepos,lastgenomepos);
	  abort();
#endif
	  fill = lastquerypos - querypos - 1;
	} else {
#if 0
	  /* This can occur with wobble mask */
	  fprintf(stderr,"Partial fill from genomepos %u to %u (querypos goes from %d to %d)\n",
		  genomepos,lastgenomepos,querypos,lastquerypos);
	  abort();
#endif
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
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos+genomic_offset,
			      genomicuc_ptr[lastgenomepos],MATCH_COMP,
			      genomicseg_ptr[lastgenomepos],/*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",genomicuc_ptr[lastgenomepos],genomicseg_ptr[lastgenomepos],
		      lastquerypos+query_offset,lastgenomepos+genomic_offset));
#else
	if (queryuc_ptr[lastquerypos] == genomicuc_ptr[lastgenomepos]) {
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos+genomic_offset,
				queryseq_ptr[lastquerypos],MATCH_COMP,
				genomicseg_ptr[lastgenomepos],/*dynprogindex*/0);
	  debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
			lastquerypos+query_offset,lastgenomepos+genomic_offset));
	} else {
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos+genomic_offset,
				queryseq_ptr[lastquerypos],MISMATCH_COMP,
				genomicseg_ptr[lastgenomepos],/*dynprogindex*/0);
	  debug5(printf("Pushing %c   %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
			lastquerypos+query_offset,lastgenomepos+genomic_offset));
	}
#endif
	--lastquerypos;
	--lastgenomepos;
      }
    }

    /* Take care of observed match */
    pair->querypos += query_offset; /* Revise coordinates */
    pair->genomepos += genomic_offset; /* Revise coordinates */
#ifdef WASTE
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
    pairs = List_push_existing(pairs,pairptr);
#endif
    lastquerypos = querypos;
    lastgenomepos = genomepos;
  }

  return List_reverse(pairs);
}




/* Returns ncovered */
int
Stage2_scan (int *stage2_source,
	     char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,
	     char *genomicseg_ptr, char *genomicuc_ptr, int genomiclength, int genomic_offset,
	     Oligoindex_T *oligoindices, int noligoindices,
	     Diagpool_T diagpool, bool debug_graphic_p, bool diagnosticp) {
  int ncovered;
  int source;
  int indexsize;
  Oligoindex_T oligoindex;
  unsigned int **mappings;
  bool *coveredp, oned_matrix_p;
  int *npositions, totalpositions;
  double pct_coverage;
  int maxnconsecutive;
  /* double diag_runtime; */
  List_T diagonals;
#ifndef USE_DIAGPOOL
  List_p;
  Diag_T diag;
#endif

  if (debug_graphic_p == true) {
    /* printf("par(mfrow=c(1,2),cex=0.2)\n"); */
    printf("par(cex=0.3)\n");
    printf("layout(matrix(c(1,2),1,2),widths=c(0.5,0.5),heights=c(1))\n");
  }

  coveredp = (bool *) CALLOC(querylength,sizeof(bool));
  mappings = (unsigned int **) CALLOC(querylength,sizeof(unsigned int *));
  npositions = (int *) CALLOC(querylength,sizeof(int));
  totalpositions = 0;
  maxnconsecutive = 0;

  source = 0;
  pct_coverage = 0.0;
  Diagpool_reset(diagpool);
  diagonals = (List_T) NULL;
  while (source < noligoindices && pct_coverage < SUFF_PCTCOVERAGE_OLIGOINDEX) {
    oligoindex = oligoindices[source];
    indexsize = Oligoindex_indexsize(oligoindex); /* Different sources can have different indexsizes */
    Oligoindex_tally(oligoindex,genomicuc_ptr,genomiclength,queryuc_ptr,querylength);
    diagonals = Oligoindex_get_mappings(diagonals,coveredp,mappings,npositions,&totalpositions,
					&oned_matrix_p,&maxnconsecutive,oligoindex,queryuc_ptr,
					querylength,genomiclength,diagpool);
    pct_coverage = Diag_update_coverage(coveredp,&ncovered,diagonals,querylength);
    if (diagnosticp) {
      printf("source = %d, ncovered = %d, pct_coverage = %f\n",source,ncovered,pct_coverage);
    }
    source++;
  }
  *stage2_source = source;

#ifdef USE_DIAGPOOL
  /* No need to free diagonals */
#else
    for (p = diagonals; p != NULL; p = List_next(p)) {
      diag = (Diag_T) List_head(p);
      Diag_free(&diag);
    }
    List_free(&diagonals);
#endif

  FREE(npositions);
  FREE(coveredp);
  FREE(mappings);		/* Don't need to free contents of mappings */

  for (source = 0; source < noligoindices; source++) {
    oligoindex = oligoindices[source];
    Oligoindex_untally(oligoindex);
  }

  return ncovered;
}


List_T
Stage2_compute (int *stage2_source, int *stage2_indexsize,
		char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,
		char *genomicseg_ptr, char *genomicuc_ptr, int genomiclength, int genomic_offset,
		Oligoindex_T *oligoindices, int noligoindices, double proceed_pctcoverage,
		Pairpool_T pairpool, Diagpool_T diagpool, int sufflookback, int nsufflookback,
		int maxintronlen, bool localp, bool skip_repetitive_p, bool debug_graphic_p,
		bool diagnosticp, Stopwatch_T stopwatch, bool diag_debug) {
  int indexsize, indexsize_nt;
  Oligoindex_T oligoindex;
  List_T path = NULL, pairs;
  unsigned int **mappings;
  bool *coveredp, oned_matrix_p;
  int source;
  int nfwdintrons, nrevintrons, nunkintrons;
  int cdna_direction = 0;
  int *npositions, totalpositions, start_querypos, end_querypos;
  unsigned int *minactive, *maxactive;
  int ncovered;
  double pct_coverage;
  int maxnconsecutive;
  /* double diag_runtime; */
  List_T diagonals;
#ifndef USE_DIAGPOOL
  List_T p;
  Diag_T diag;
#endif
#ifdef DEBUG
  int nunique;
#endif

  Stopwatch_start(stopwatch);

  if (debug_graphic_p == true) {
    /* printf("par(mfrow=c(1,2),cex=0.2)\n"); */
    printf("par(cex=0.3)\n");
    printf("layout(matrix(c(1,2),1,2),widths=c(0.5,0.5),heights=c(1))\n");
  }

  coveredp = (bool *) CALLOC(querylength,sizeof(bool));
  mappings = (unsigned int **) CALLOC(querylength,sizeof(unsigned int *));
  npositions = (int *) CALLOC(querylength,sizeof(int));
  totalpositions = 0;
  maxnconsecutive = 0;

  source = 0;
  pct_coverage = 0.0;
#ifdef USE_DIAGPOOL
  Diagpool_reset(diagpool);
#endif
  diagonals = (List_T) NULL;
  while (source < noligoindices && pct_coverage < SUFF_PCTCOVERAGE_OLIGOINDEX) {
    oligoindex = oligoindices[source];
    indexsize = Oligoindex_indexsize(oligoindex); /* Different sources can have different indexsizes */
    Oligoindex_tally(oligoindex,genomicuc_ptr,genomiclength,queryuc_ptr,querylength);
    diagonals = Oligoindex_get_mappings(diagonals,coveredp,mappings,npositions,&totalpositions,
					&oned_matrix_p,&maxnconsecutive,oligoindex,queryuc_ptr,
					querylength,genomiclength,diagpool);
    pct_coverage = Diag_update_coverage(coveredp,&ncovered,diagonals,querylength);
    if (diagnosticp) {
      printf("source = %d, ncovered = %d, pct_coverage = %f\n",source,ncovered,pct_coverage);
    }
    source++;
  }
  *stage2_source = source;
  *stage2_indexsize = indexsize;

  /* diag_runtime = */ Stopwatch_stop(stopwatch);


  minactive = (unsigned int *) CALLOC(querylength,sizeof(unsigned int));
  maxactive = (unsigned int *) CALLOC(querylength,sizeof(unsigned int));


  Stopwatch_start(stopwatch);

  if (diag_debug == true) {
    /* Do nothing */
  } else if (totalpositions == 0) {
    debug(printf("Quitting because totalpositions is zero\n"));
    path = (List_T) NULL;
  } else if (pct_coverage < proceed_pctcoverage && ncovered < SUFF_NCOVERED) {
    debug(printf("Quitting because pct_coverage is only %f and ncovered is only %d, maxnconsecutive = %d\n",
		  pct_coverage,ncovered,maxnconsecutive));
    path = (List_T) NULL;
  } else {
    debug(printf("Proceeding because maxnconsecutive is %d and pct_coverage is %f > %f or ncovered = %d > %d\n",
		 maxnconsecutive,pct_coverage,proceed_pctcoverage,ncovered,SUFF_NCOVERED));

    debug(printf("Performing diag on genomiclength %u\n",genomiclength));
    Diag_compute_bounds(minactive,maxactive,diagonals,genomiclength,querylength,
			indexsize,debug_graphic_p,diagnosticp,queryuc_ptr,genomicuc_ptr);
    
    debug(
	  nunique = Diag_compute_bounds(minactive,maxactive,diagonals,genomiclength,querylength,
					indexsize,debug_graphic_p,diagnosticp,queryuc_ptr,genomicuc_ptr);
	  fprintf(stderr,"%d diagonals (%d not dominated), maxnconsecutive = %d\n",
		  List_length(diagonals),nunique,maxnconsecutive);
	  );

    if (debug_graphic_p == true) {
      active_bounds_dump_R(minactive,maxactive,querylength);
      printf("lines(querypos,minactive,col=\"blue\")\n");
      printf("lines(querypos,maxactive,col=\"blue\")\n");
    }

    path = align_compute(&start_querypos,&end_querypos,
			 &nfwdintrons,&nrevintrons,&nunkintrons,&cdna_direction,
			 mappings,npositions,totalpositions,oned_matrix_p,minactive,maxactive,
			 queryseq_ptr,querylength,/*query_trim_start*/0,/*query_trim_end*/querylength,
			 genomicseg_ptr,genomicuc_ptr,genomiclength,
			 indexsize,sufflookback,nsufflookback,
			 maxintronlen,pairpool,localp,skip_repetitive_p,debug_graphic_p);

    debug(printf("aligned from %d to %d\n",start_querypos,end_querypos));
  }

  FREE(maxactive);
  FREE(minactive);
  FREE(npositions);
  FREE(coveredp);
  FREE(mappings);		/* Don't need to free contents of mappings */

  for (source = 0; source < noligoindices; source++) {
    oligoindex = oligoindices[source];
    Oligoindex_untally(oligoindex);
  }

#ifdef PMAP
  indexsize_nt = 3*indexsize;
#else
  indexsize_nt = indexsize;
#endif

  Stopwatch_stop(stopwatch);

  if (diag_debug == true) {
    return diagonals;
  } else {

#ifdef USE_DIAGPOOL
  /* No need to free diagonals */
#else
    for (p = diagonals; p != NULL; p = List_next(p)) {
      diag = (Diag_T) List_head(p);
      Diag_free(&diag);
    }
    List_free(&diagonals);
#endif

    if (path == NULL) {
      debug(printf("Couldn't find alignment in stage 2\n"));
      return (List_T) NULL;
    } else {
      pairs = convert_to_nucleotides(List_reverse(path),
#ifndef PMAP
				     queryseq_ptr,queryuc_ptr,
#endif
				     genomicseg_ptr,genomicuc_ptr,
				     query_offset,genomic_offset,pairpool,indexsize_nt);
      /* Don't need to free path, because its memory belongs to pairpool */
      return pairs;
    }
  }

}


