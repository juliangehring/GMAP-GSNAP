static char rcsid[] = "$Id: stage2.c,v 1.154 2006/04/07 01:21:34 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "listdef.h"

/* Penalty for genomic distances */

#define INTRON_PENALTY_CONSISTENT 0
#ifdef PMAP
#define INTRON_PENALTY_UNKNOWN 0 /* protein boundaries may not reveal GT-AG */
#else
#define INTRON_PENALTY_UNKNOWN 8
#define INTRON_PENALTY_INCONSISTENT 16
#endif

#define NINTRON_PENALTY 8
#define NINTRON_PENALTY_MISMATCH 8

#define ENOUGH_CONSECUTIVE 24

#define INFINITY 1000000

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
#else
#define SAMPLE_INTERVAL 2	/* For cases where adjacentp == false.
				   Means that we can find islands of
				   9-mers */
#define NT_PER_MATCH 1
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

/* Print all links.  Warning: lots of data */
#ifdef DEBUG1
#define debug1(x) x
#else 
#define debug1(x)
#endif

/* Not used */
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


#define T Stage2_T
struct T {
  List_T path;
  double defect_rate;
  int nfwdintrons;
  int nrevintrons;
  int nnonintrons;
  int nunkintrons;
};

List_T
Stage2_path (T this) {
  return this->path;
}

double
Stage2_defect_rate (T this) {
  return this->defect_rate;
}

int
Stage2_nfwdintrons (T this) {
  return this->nfwdintrons;
}

int
Stage2_nrevintrons (T this) {
  return this->nrevintrons;
}

int
Stage2_nunkintrons (T this) {
  return this->nunkintrons;
}

int
Stage2_pathlength (T this) {
  return List_length(this->path);
}

static T
Stage2_new (List_T path, double defect_rate, 
	    int nfwdintrons, int nrevintrons, int nunkintrons) {
  T new = (T) MALLOC(sizeof(*new));

  new->path = path;
  new->defect_rate = defect_rate;
  new->nfwdintrons = nfwdintrons;
  new->nrevintrons = nrevintrons;
  new->nunkintrons = nunkintrons;

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
  int consecutive;
  int fwd_score;
  int fwd_pos;
  int fwd_hit;
  int fwd_intronnfwd;
  int fwd_intronnrev;
  int fwd_intronnunk;
#ifndef PMAP
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


static void
Linkmatrix_print_fwd (struct Link_T **links, int length1, int *npositions) {
  int i, j;

  for (i = 0; i < length1; i++) {
    printf("Row %d (%d positions):",i,npositions[i]);
    for (j = 0; j < npositions[i]; j++) {
      printf(" %d",links[i][j].fwd_pos);
    }
    printf("\n");
  }
  printf("\n");
  return;
}

#ifndef PMAP
static void
Linkmatrix_print_rev (struct Link_T **links, int length1, int *npositions) {
  int i, j;

  for (i = 0; i < length1; i++) {
    printf("Row %d (%d positions):",i,npositions[i]);
    for (j = 0; j < npositions[i]; j++) {
      printf(" %d",links[i][j].rev_pos);
    }
    printf("\n");
  }
  printf("\n");
  return;
}
#endif
    


/* The value of this procedure determines the number of matches worth
   having after a gap.  For some reason, indexsize - 4 finds exons of
   size indexsize, but indexsize - 3 does not. */
#define GAPEQUALIZER 4

#define distpenalty_consistent(intronlength,crossspeciesp) (crossspeciesp == true ? intronlength >> 9 : intronlength >> 10)
#define distpenalty_unknown(intronlength,crossspeciesp) (crossspeciesp == true ? intronlength >> 6 : intronlength >> 10)
#define distpenalty_inconsistent(intronlength,crossspeciesp) (crossspeciesp == true ? intronlength >> 6 : intronlength >> 10)
#define distpenalty_dead(intronlength,crossspeciesp) (crossspeciesp == true ? intronlength >> 13 : intronlength >> 15)


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


static void
score_querypos (Link_T currlink, int querypos, unsigned int position, char acceptor1, char acceptor2,
		struct Link_T **links, unsigned int **mappings, int *npositions, 
		int *fwd_restrict_hit, 
#ifndef PMAP
		int *rev_restrict_hit,
#endif
		char *genomicuc_ptr, int indexsize, int sufflookback, 
		int nsufflookback, int maxintronlen, bool deadp, bool crossspeciesp) {
  Link_T prevlink;
  int best_consecutive = 0, save_consecutive;
  int best_fwd_score = 0, fwd_score;
  int best_fwd_prevpos = -1, best_fwd_prevhit = -1;
  int best_fwd_intronnfwd = 0, best_fwd_intronnrev = 0, best_fwd_intronnunk = 0;
#ifndef PMAP
  int best_rev_score = 0, rev_score;
  int best_rev_prevpos = -1, best_rev_prevhit = -1;
  int best_rev_intronnfwd = 0, best_rev_intronnrev = 0, best_rev_intronnunk = 0;
#endif
  bool adjacentp = false;
  int prev_querypos, prevhit, lastprevhit;
  unsigned int *prev_mappings;
  unsigned int prevposition, gendistance;
  int querydistance, diffdistance, lookback, nlookback, sample_interval, nseen;
  char donor1;
  int canonicalsgn, mismatch_start, fwd_gendistance_penalty;
#ifndef PMAP
  int rev_gendistance_penalty;
#endif

  /* A. Evaluate adjacent position (at querypos - 1) */
  for (prevhit = 0; prevhit < npositions[querypos-1]; prevhit++) {
    if (mappings[querypos-1][prevhit] == position - NT_PER_MATCH) {
      prevlink = &(links[querypos-1][prevhit]);
      best_consecutive = currlink->consecutive = prevlink->consecutive + NT_PER_MATCH;
      best_fwd_score = prevlink->fwd_score + NT_PER_MATCH;
#ifndef PMAP
      best_rev_score = prevlink->rev_score + NT_PER_MATCH;
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
      debug(printf("\tAdjacent qpos %d,%d at %u (scores = %d -> %d, consec = %d, intr = %d-%d-%d)\n",
		   querypos-1,prevhit,position - NT_PER_MATCH,prevlink->fwd_score,
		   best_fwd_score,best_consecutive,best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk));
#else
      debug(printf("\tAdjacent qpos %d,%d at %u (scores = %d/%d -> %d/%d, consec = %d, intr = %d-%d-%d/%d-%d-%d)\n",
		   querypos-1,prevhit,position - NT_PER_MATCH,prevlink->fwd_score,prevlink->rev_score,
		   best_fwd_score,best_rev_score,best_consecutive,
		   best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		   best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk));
#endif
    }
  }

  if (best_consecutive < ENOUGH_CONSECUTIVE) {
    /* B. Evaluate for intron (at querypos - indexsize).  
       Don't update best_consecutive, because a reasonable alternative 
       might be a mismatch, regardless of the quality of the intron. */
#ifdef PMAP
    prev_querypos = querypos - indexsize/3;
#else
    prev_querypos = querypos - indexsize;
#endif
    if (prev_querypos >= 0 && npositions[prev_querypos] > 0) {
#ifdef PMAP
      querydistance = (querypos - prev_querypos)*3;
#else
      querydistance = querypos - prev_querypos;
#endif
      prev_mappings = &(mappings[prev_querypos][0]);
      for (prevhit = 0; prevhit < npositions[prev_querypos]; prevhit++) {
	prevposition = *(prev_mappings++);
	if (position >= prevposition + indexsize) {
	  prevlink = &(links[prev_querypos][prevhit]);
	  if ((/* donor2 = */ genomicuc_ptr[prevposition+indexsize+1]) != 'T' || acceptor2 != 'A') {
	    canonicalsgn = 0;
	  } else if ((donor1 = genomicuc_ptr[prevposition+indexsize]) == 'G' && acceptor1 == 'G') {
	    canonicalsgn = +1;
	  } else if (donor1 == 'C' && acceptor1 == 'C') {
	    canonicalsgn = -1;
	  } else {
	    canonicalsgn = 0;
	  }

	  /* Add indexsize to compensate for the 8-mer lost across
	     an intron.  Subtract NINTRON_PENALTY to penalize for each
	     additional intron */
	  fwd_score = prevlink->fwd_score + indexsize - NINTRON_PENALTY;
#ifndef PMAP
	  rev_score = prevlink->rev_score + indexsize - NINTRON_PENALTY;
#endif

	  /* Penalty for introns. */
	  gendistance = position - prevposition;
	  if (gendistance > maxintronlen) {
	    fwd_gendistance_penalty = rev_gendistance_penalty = INFINITY;

	  } else if (deadp == true) {
	    fwd_gendistance_penalty = 
#ifndef PMAP
	      rev_gendistance_penalty = 
#endif
	      distpenalty_dead(gendistance,crossspeciesp) + INTRON_PENALTY_UNKNOWN + DEAD_PENALTY;
	  } else {
	    if (canonicalsgn == 0) {
	      /* Extra penalty for known non-canonical intron */
	      fwd_gendistance_penalty = 
#ifndef PMAP
		rev_gendistance_penalty = 
#endif
		distpenalty_unknown(gendistance,crossspeciesp) + INTRON_PENALTY_UNKNOWN;
	    } else if (canonicalsgn == +1) {
	      fwd_gendistance_penalty = distpenalty_consistent(gendistance,crossspeciesp) + INTRON_PENALTY_CONSISTENT;
#ifndef PMAP
	      rev_gendistance_penalty = distpenalty_inconsistent(gendistance,crossspeciesp) + INTRON_PENALTY_INCONSISTENT;
#endif
	    } else if (canonicalsgn == -1) {
#ifdef PMAP
	      fwd_gendistance_penalty = distpenalty_unknown(gendistance,crossspeciesp) + INTRON_PENALTY_UNKNOWN;
#else
	      fwd_gendistance_penalty = distpenalty_inconsistent(gendistance,crossspeciesp) + INTRON_PENALTY_INCONSISTENT;
	      rev_gendistance_penalty = distpenalty_consistent(gendistance,crossspeciesp) + INTRON_PENALTY_CONSISTENT;
#endif
	    }
	  }
	  fwd_score -= fwd_gendistance_penalty;
#ifndef PMAP
	  rev_score -= rev_gendistance_penalty;
#endif

#ifdef PMAP
	  debug(printf("\tIntron (%d:%c%c-%c%c) qpos %d,%d at %u (scores = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, distpens = %d)",
		       canonicalsgn,genomicuc_ptr[prevposition+indexsize],genomicuc_ptr[prevposition+indexsize+1],acceptor2,acceptor1,
		       prev_querypos,prevhit,prevposition,
		       prevlink->fwd_score,fwd_score,prevlink->consecutive,
		       best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		       gendistance,fwd_gendistance_penalty));
#else
	  debug(printf("\tIntron (%d:%c%c-%c%c) qpos %d,%d at %u (scores = %d/%d -> %d/%d, consec = %d, intr = %d-%d-%d/%d-%d-%d, gendist %u, distpens = %d/%d)",
		       canonicalsgn,genomicuc_ptr[prevposition+indexsize],genomicuc_ptr[prevposition+indexsize+1],acceptor2,acceptor1,
		       prev_querypos,prevhit,prevposition,
		       prevlink->fwd_score,prevlink->rev_score,fwd_score,rev_score,prevlink->consecutive,
		       best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		       best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk,
		       gendistance,fwd_gendistance_penalty,rev_gendistance_penalty));
#endif
	  
	  debug(printf(" =>"));
	  if (fwd_score > best_fwd_score) {
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
	  debug(printf("\n\n"));
	}
      }
    }
    mismatch_start = prev_querypos - 1;
    save_consecutive = best_consecutive;

    /* C (fwd). Evaluate for mismatches (all other previous querypos) */
    /* Set parameters */
    if (adjacentp == true) {
      /* Sample at larger interval and not so far back */
      nlookback = 1;
      lookback = sufflookback/2;
#ifdef PMAP
      sample_interval = indexsize/3;
#else
      sample_interval = indexsize;
#endif
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
	 prev_querypos >= 0 && best_consecutive < ENOUGH_CONSECUTIVE && 
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
	  if (position >= prevposition + indexsize) {
	    prevlink = &(links[prev_querypos][prevhit]);

	    gendistance = position - prevposition;
	    diffdistance = abs(gendistance - querydistance);
	    fwd_score = prevlink->fwd_score;
	    if (gendistance > maxintronlen) {
	      fwd_score -= INFINITY;
	    } else if (querydistance < INTRON_DEFN) {
	      fwd_score += querydistance;
	    } else if (diffdistance < EQUAL_DISTANCE) {
	      fwd_score += querydistance - diffdist_mismatch_penalty(diffdistance,crossspeciesp);
	    } else if (diffdistance < INTRON_DEFN) {
	      fwd_score -= querydist_mismatch_penalty(querydistance,crossspeciesp);
	    } else {
	      fwd_score -= querydist_mismatch_penalty(querydistance,crossspeciesp) + NINTRON_PENALTY_MISMATCH;
	    }
	    if (deadp == true) {
	      fwd_score -= distpenalty_dead(diffdistance,crossspeciesp) + DEAD_PENALTY;
	    } else {
	      fwd_score -= distpenalty_unknown(diffdistance,crossspeciesp);
	    }

	    debug(printf("\tFwd mismatch qpos %d,%d at %u (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d)",
			 prev_querypos,prevhit,prevposition,prevlink->fwd_score,fwd_score,prevlink->consecutive,
			 best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			 gendistance,querydistance));
	    
	    if (fwd_score > best_fwd_score) {
	      best_consecutive = prevlink->consecutive;
	      best_fwd_score = fwd_score;
	      best_fwd_prevpos = prev_querypos;
	      best_fwd_prevhit = prevhit;
	      best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	      best_fwd_intronnrev = prevlink->fwd_intronnrev;
	      best_fwd_intronnunk = prevlink->fwd_intronnunk;
	      if (diffdistance < INTRON_DEFN) {
		best_fwd_intronnunk++;
	      }
#ifdef PMAP
	      sample_interval = indexsize/3;
#else
	      sample_interval = indexsize;
#endif
	      debug(printf(" => Best fwd at %d (best_consec = %d)\n",fwd_score,best_consecutive));
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
    best_consecutive = save_consecutive;
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
	 prev_querypos >= 0 && best_consecutive < ENOUGH_CONSECUTIVE && 
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
	  if (position >= prevposition + indexsize) {
	    prevlink = &(links[prev_querypos][prevhit]);
	    
	    gendistance = position - prevposition;
	    diffdistance = abs(gendistance - querydistance);
	    rev_score = prevlink->rev_score;
	    if (gendistance > maxintronlen) {
	      rev_score -= INFINITY;
	    } else if (querydistance < INTRON_DEFN) {
	      rev_score += querydistance;
	    } else if (diffdistance < EQUAL_DISTANCE) {
	      rev_score += querydistance - diffdist_mismatch_penalty(diffdistance,crossspeciesp);
	    } else if (diffdistance < INTRON_DEFN) {
	      rev_score -= querydist_mismatch_penalty(querydistance,crossspeciesp);
	    } else {
	      rev_score -= querydist_mismatch_penalty(querydistance,crossspeciesp) + NINTRON_PENALTY_MISMATCH;
	    }
	    if (deadp == true) {
	      rev_score -= distpenalty_dead(diffdistance,crossspeciesp) + DEAD_PENALTY;
	    } else {
	      rev_score -= distpenalty_unknown(diffdistance,crossspeciesp);
	    }
	    
	    debug(printf("\tRev mismatch qpos %d,%d at %u (score = %d -> %d, consec = %d, intr = %d-%d-%d, gendist %u, querydist %d)",
			 prev_querypos,prevhit,prevposition,prevlink->rev_score,rev_score,prevlink->consecutive,
			 best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk,
			 gendistance,querydistance));
	    
	    if (rev_score > best_rev_score) {
	      best_consecutive = prevlink->consecutive;
	      best_rev_score = rev_score;
	      best_rev_prevpos = prev_querypos;
	      best_rev_prevhit = prevhit;
	      best_rev_intronnfwd = prevlink->rev_intronnfwd;
	      best_rev_intronnrev = prevlink->rev_intronnrev;
	      best_rev_intronnunk = prevlink->rev_intronnunk;
	      if (diffdistance < INTRON_DEFN) {
		best_rev_intronnunk++;
	      }
	      sample_interval = indexsize;
	      debug(printf(" => Best rev at %d (best_consec = %d)\n",rev_score,best_consecutive));
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
  currlink->fwd_pos = best_fwd_prevpos;
  currlink->fwd_hit = best_fwd_prevhit;
  currlink->fwd_score = best_fwd_score;
  currlink->fwd_intronnfwd = best_fwd_intronnfwd;
  currlink->fwd_intronnrev = best_fwd_intronnrev;
  currlink->fwd_intronnunk = best_fwd_intronnunk;

#ifndef PMAP
  /* linksconsecutive already assigned above */
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
align_compute_scores (Link_T *termlink, bool *fwdp, struct Link_T **links,
		      unsigned int **mappings, int *npositions, int querystart, int queryend,
		      int querylength, Sequence_T genomicuc, int indexsize, 
		      int sufflookback, int nsufflookback, int maxintronlen, bool crossspeciesp) {
  Link_T currlink;
  int fwd_score, *fwd_restrict_hit;
#ifndef PMAP
  int rev_score, *rev_restrict_hit;
#endif
  int querypos, hit;
  int best_score, best_hit, second_score, score;
  unsigned int position;
  char *genomicuc_ptr;
  char acceptor2, acceptor1;

  genomicuc_ptr = Sequence_fullpointer(genomicuc);

  debug0(printf("Querystart = %d, queryend = %d\n",querystart,queryend));

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
    links[querystart][hit].fwd_pos = links[querystart][hit].fwd_hit = -1;
#ifndef PMAP
    links[querystart][hit].rev_pos = links[querystart][hit].rev_hit = -1;
#endif
    /* links[querystart][hit].score = 0; */
    /* links[querystart][hit].consecutive = 0; */
  }

  for (querypos = querystart+1; 
#ifdef PMAP
       querypos <= queryend - indexsize/3; 
#else
       querypos <= queryend - indexsize;
#endif
       querypos++) {

    for (hit = 0; hit < npositions[querypos]; hit++) {
      currlink = &(links[querypos][hit]);
      position = mappings[querypos][hit];
      if (position < 1) {
	acceptor2 = acceptor1 = 'X';
      } else if (position == 1) {
	acceptor2 = 'X';
	acceptor1 = genomicuc_ptr[0];
      } else {
	acceptor2 = genomicuc_ptr[position-2];
	acceptor1 = genomicuc_ptr[position-1];
      }

      debug(printf("Finding link at querypos %d,%d at %u\n",querypos,hit,position));

      score_querypos(currlink,querypos,position,acceptor1,acceptor2,
		     links,mappings,npositions,fwd_restrict_hit,
#ifndef PMAP
		     rev_restrict_hit,
#endif
		     genomicuc_ptr,indexsize,sufflookback,nsufflookback,maxintronlen,
		     /*deadp*/false,crossspeciesp);
      
      if (currlink->fwd_score <= 0
#ifndef PMAP
	  || currlink->rev_score <= 0
#endif
	  ) {
	debug(printf("Redoing link at querypos %d,%d at %u, without gendistance penalty\n",querypos,hit,position));

	score_querypos(currlink,querypos,position,acceptor1,acceptor2,
		       links,mappings,npositions,fwd_restrict_hit,
#ifndef PMAP
		       rev_restrict_hit,
#endif
		       genomicuc_ptr,indexsize,sufflookback,nsufflookback,maxintronlen,
		       /*deadp*/true,crossspeciesp);
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
    if (best_score - second_score > indexsize - GAPEQUALIZER) {
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
    if (best_score - second_score > indexsize - GAPEQUALIZER) {
      rev_restrict_hit[querypos] = best_hit;
    } else {
      rev_restrict_hit[querypos] = -1;
    }
#endif
  }

  /* Find best match */
  best_score = 0;
  *termlink = NULL;
  /* Go backwards, because best score is likely to be there, and so
     ties favor longer path */

#ifdef PMAP
  *fwdp = true;
  for (querypos = queryend - indexsize/3; querypos > querystart; --querypos) {
    /* Don't use restrict_hit here, because it doesn't save time */
    for (hit = 0; hit < npositions[querypos]; hit++) {
      if ((score = links[querypos][hit].fwd_score) > best_score) {
	best_score = score;
	*termlink = &(links[querypos][hit]);
      }
    }
  }
  debug(printf("==> Found best score %d at fwd %d,%d\n",best_score,(*termlink)->fwd_pos,(*termlink)->fwd_hit));

#else

  for (querypos = queryend - indexsize; querypos > querystart; --querypos) {
    /* Don't use restrict_hit here, because it doesn't save time */
    for (hit = 0; hit < npositions[querypos]; hit++) {
      if ((score = links[querypos][hit].fwd_score) > best_score) {
	best_score = score;
	*termlink = &(links[querypos][hit]);
	*fwdp = true;
      }
      if ((score = links[querypos][hit].rev_score) > best_score) {
	best_score = score;
	*termlink = &(links[querypos][hit]);
	*fwdp = false;
      }
    }
  }

  debug(
	if (*fwdp == true) {
	  printf("==> Found best score %d at fwd %d,%d\n",best_score,(*termlink)->fwd_pos,(*termlink)->fwd_hit);
	} else {
	  printf("==> Found best score %d at rev %d,%d\n",best_score,(*termlink)->rev_pos,(*termlink)->rev_hit);
	}
	);
#endif

#ifndef PMAP
  FREE(rev_restrict_hit);
#endif
  FREE(fwd_restrict_hit);

  return;
}



/* Performs dynamic programming. */
static List_T
align_compute (int *pathlength, double *defect_rate, int *nfwdintrons, int *nrevintrons, int *nunkintrons,
	       unsigned int **mappings, int *npositions, 
	       int totalpositions, Sequence_T queryseq, Sequence_T genomicseg, Sequence_T genomicuc,
	       int indexsize, int sufflookback, int nsufflookback, int maxintronlen, Pairpool_T pairpool,
	       bool crossspeciesp) {
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

  queryseq_ptr = Sequence_fullpointer(queryseq);
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
#ifdef PMAP
  genomicuc_ptr = Sequence_fullpointer(genomicuc);
#endif

  querystart = Sequence_trim_start(queryseq);
  queryend = Sequence_trim_end(queryseq);
  querylength = Sequence_fulllength(queryseq);

  links = Linkmatrix_new(querylength,npositions,totalpositions);

  align_compute_scores(&termlink,&fwdp,links,mappings,npositions,querystart,queryend,querylength,
		       genomicuc,indexsize,sufflookback,nsufflookback,maxintronlen,crossspeciesp);

  if (termlink != NULL) {
    if (fwdp) {
      termpos = termlink->fwd_pos;
      termhit = termlink->fwd_hit;
      *nfwdintrons = termlink->fwd_intronnfwd;
      *nrevintrons = termlink->fwd_intronnrev;
      *nunkintrons = termlink->fwd_intronnunk;
#ifndef PMAP
    } else {
      termpos = termlink->rev_pos;
      termhit = termlink->rev_hit;
      *nfwdintrons = termlink->rev_intronnfwd;
      *nrevintrons = termlink->rev_intronnrev;
      *nunkintrons = termlink->rev_intronnunk;
#endif
    }

#ifdef PMAP
    debug1(Linkmatrix_print_fwd(links,querylength,npositions));
#else
    debug1(
	   if (fwdp) {
	     Linkmatrix_print_fwd(links,querylength,npositions);
	   } else {
	     Linkmatrix_print_rev(links,querylength,npositions);
	   }
	   );
#endif

    querypos = termpos;
    hit = termhit;

    /* Traceback */
    /* nconsecutive = nintrons = ndefects = 0; */
    if (querypos != -1) {
      prevposition = mappings[querypos][hit]+1;
    }
    while (querypos != -1) {
      position = mappings[querypos][hit];
#ifdef PMAP
      debug0(printf("Pushing %d,%d at %u, score = %d, intr = %d(+)/%d(-)/%d(?)\n",
		    querypos,hit,position,links[querypos][hit].fwd_score,
		    links[querypos][hit].fwd_intronnfwd,links[querypos][hit].fwd_intronnrev,
		    links[querypos][hit].fwd_intronnunk));
#else
      debug0(
	     if (fwdp) {
	       printf("Pushing %d,%d at %u, score = %d, intr = %d(+)/%d(-)/%d(?)\n",
		      querypos,hit,position,links[querypos][hit].fwd_score,
		      links[querypos][hit].fwd_intronnfwd,links[querypos][hit].fwd_intronnrev,
		      links[querypos][hit].fwd_intronnunk);
	     } else {
	       printf("Pushing %d,%d at %u, score = %d, intr = %d(+)/%d(-)/%d(?)\n",
		      querypos,hit,position,links[querypos][hit].rev_score,
		      links[querypos][hit].rev_intronnfwd,links[querypos][hit].rev_intronnrev,
		      links[querypos][hit].rev_intronnunk);
	     });
#endif

#ifdef PMAP
      /* Change querypos positions from protein to nucleotide */
      path = Pairpool_push(path,pairpool,querypos*3+2,position+2,
 			   genomicuc_ptr[position+2],MATCH_COMP,
 			   genomicseg_ptr[position+2],/*gapp*/false);
      path = Pairpool_push(path,pairpool,querypos*3+1,position+1,
 			   genomicuc_ptr[position+1],MATCH_COMP,
 			   genomicseg_ptr[position+1],/*gapp*/false);
      path = Pairpool_push(path,pairpool,querypos*3,position,
 			   genomicuc_ptr[position],MATCH_COMP,
			   genomicseg_ptr[position],/*gapp*/false);
#else
      path = Pairpool_push(path,pairpool,querypos,position,
			   queryseq_ptr[querypos],MATCH_COMP,
			   genomicseg_ptr[position],/*gapp*/false);
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
      *defect_rate = 1.0;
    } else {
      *defect_rate = (double) (ndefects)/(double) (nconsecutive + ndefects);
    }
    *pathlength = nconsecutive + nintrons + ndefects;
  }
  Linkmatrix_free(&links);

  return path;			/* previously was List_reverse(path) */
}


static List_T
fill_oligo (List_T pairs, int querypos, int genomepos, char *queryseq_ptr, char *genomicseg_ptr,
#ifdef PMAP
	    char *genomicuc_ptr,
#endif
	    Pairpool_T pairpool, int indexsize) {
  int lastquerypos, lastgenomepos;
  
#ifdef PMAP
  lastquerypos = querypos + indexsize - 2;
  lastgenomepos = genomepos + indexsize - 2;
#else
  lastquerypos = querypos + indexsize;
  lastgenomepos = genomepos + indexsize;
#endif

  while (querypos < lastquerypos - 1 &&
	 genomepos < lastgenomepos - 1) {
#ifdef PMAP
    pairs = Pairpool_push(pairs,pairpool,lastquerypos-1,lastgenomepos-1,
			  genomicuc_ptr[lastgenomepos-1],MATCH_COMP,
			  genomicseg_ptr[lastgenomepos-1],/*gapp*/false);
#else
    pairs = Pairpool_push(pairs,pairpool,lastquerypos-1,lastgenomepos-1,
			  queryseq_ptr[lastquerypos-1],MATCH_COMP,
			  genomicseg_ptr[lastgenomepos-1],/*gapp*/false);
#endif

    /*
    printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos-1],genomicseg_ptr[lastgenomepos-1],
	   lastquerypos-1,lastgenomepos-1);
    */

    lastquerypos--;
    lastgenomepos--;
  }

  return pairs;
}

/* queryseq_ptr is NULL for PMAP */
static List_T
convert_to_nucleotides (List_T path, char *queryseq_ptr, char *genomicseg_ptr,
#ifdef PMAP
			char *genomicuc_ptr,
#endif
			Pairpool_T pairpool, int indexsize) {
  List_T pairs;
  Pair_T pair, firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump;

  firstpair = path->first;
  querypos = firstpair->querypos;
  genomepos = firstpair->genomepos;
  pairs = fill_oligo(NULL,querypos,genomepos,queryseq_ptr,genomicseg_ptr,
#ifdef PMAP
		     genomicuc_ptr,
#endif
		     pairpool,indexsize);

  /* Take care of first pair */
  pairs = Pairpool_push_existing(pairs,pairpool,firstpair);
  lastquerypos = querypos;
  lastgenomepos = genomepos;
  path = path->rest;

  while (path != NULL) {
    pair = path->first;
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - querypos - 1;
    genomejump = lastgenomepos - genomepos - 1;
    if (queryjump == 0 && genomejump == 0) {
      /* Do nothing */
    } else if (querypos+indexsize-1 < lastquerypos && genomepos+indexsize-1 < lastgenomepos) {
      pairs = fill_oligo(pairs,querypos,genomepos,queryseq_ptr,genomicseg_ptr,
#ifdef PMAP
			 genomicuc_ptr,
#endif
			 pairpool,indexsize);
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
Stage2_compute (int *maxconsecutive, Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		Sequence_T queryntseq,
#endif
		Sequence_T genomicseg, Sequence_T genomicuc, 
		Oligoindex_T oligoindex, int indexsize, Pairpool_T pairpool, 
		int sufflookback, int nsufflookback, int maxintronlen, int badoligos,
		bool crossspeciesp) {
  T stage2;
  List_T path, pairs;
  unsigned int **mappings;
  int nfwdintrons, nrevintrons, nunkintrons;
  int *npositions, totalpositions, nconsecutive, pathlength;
  double defect_rate;
#ifdef PMAP
  List_T p;
  Pair_T pair;
  char *queryntseq_ptr;
#endif

  Oligoindex_tally(oligoindex,genomicuc);
  mappings = Oligoindex_get_mappings(&nconsecutive,&npositions,&totalpositions,oligoindex,queryuc);

  /* printf("nconsecutive = %d vs maxconsecutive = %d\n",nconsecutive,*maxconsecutive); */
  if (nconsecutive > *maxconsecutive) {
    *maxconsecutive = nconsecutive;
  } else if (nconsecutive < 8*(*maxconsecutive)/10) {
    FREE(npositions);
    FREE(mappings);
    Oligoindex_cleanup(oligoindex,queryuc);
    return NULL;
  }

  if (totalpositions == 0) {
    path = NULL;
    defect_rate = 1.0;
  } else {
    if (badoligos > 16) {
      sufflookback = sufflookback*3/2;
    }
    path = align_compute(&pathlength,&defect_rate,&nfwdintrons,&nrevintrons,&nunkintrons,
			 mappings,npositions,totalpositions,
			 queryseq,genomicseg,genomicuc,
			 indexsize,sufflookback,nsufflookback,maxintronlen,pairpool,
			 crossspeciesp);
    /* printf("pathlength %d vs trimlength/10 = %d\n",pathlength,Sequence_trimlength(queryseq)/10); */
    if (pathlength < Sequence_trimlength(queryseq)/10) {
      /* Don't free path, because its memory belongs to pairpool */
      path = (List_T) NULL;
    }
  }
  Oligoindex_cleanup(oligoindex,queryuc);

  if (path == NULL) {
    debug(printf("Couldn't find alignment in stage 2\n"));
    pairs = (List_T) NULL;
    stage2 = NULL;
  } else {
#ifdef PMAP
    pairs = convert_to_nucleotides(List_reverse(path),(char *) NULL,
				   Sequence_fullpointer(genomicseg),Sequence_fullpointer(genomicuc),
				   pairpool,indexsize);
#else
    pairs = convert_to_nucleotides(List_reverse(path),Sequence_fullpointer(queryseq),
				   Sequence_fullpointer(genomicseg),pairpool,indexsize);
#endif
    stage2 = Stage2_new(pairs,defect_rate,nfwdintrons,nrevintrons,nunkintrons);
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
  FREE(npositions);
  FREE(mappings);

  return stage2;
}
