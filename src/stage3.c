static char rcsid[] = "$Id: stage3.c,v 1.329 2010/02/26 23:11:40 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage3.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include <math.h>		/* For pow() */
#include "mem.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "listdef.h"
#include "comp.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "smooth.h"
#include "scores.h"
#include "intron.h"
#include "pbinom.h"
#include "changepoint.h"
#include "translation.h"
#ifdef PMAP
#include "backtranslation.h"
#endif
#include "complement.h"
#include "iit-read.h"
#include "stage2.h"


/* The following are the same as dividing by 2048 and 1024 */
#define goodness_intronlen(x) (x >> 11)
#define goodness_nonintronlen(x) (x >> 10)


#define MAXITER 2
#define INTRON_PENALTY_INCONSISTENT 16
#define NONCANONICAL_PENALTY 12

#define SINGLESLEN 9	/* Should be same as MININTRONLEN */
#define MININTRONLEN 9		/* Determines when Dynprog_genome_gap gets called vs Dynprog_single_gap */
#define MININTRONLEN_FINAL 50	/* Determines when to perform final
				   pass to find canonical introns */
#define MINENDEXON 12

#define MIN_NONINTRON 50

#define SUFF_MATCHES_KEEP 300

#define SUFFCONSECUTIVE 5
#define MAXINCURSION 5

/* For Stage3_append */
#define MERGELENGTH 100000
#define LONG_MERGELENGTH 500000	/* For strong donor and acceptor splice sites */
#define DONOR_THRESHOLD 0.90
#define ACCEPTOR_THRESHOLD 0.90

#define MINCOVERAGE 0.10  /* Not used anymore */
#define MIN_MATCHES 20

#define DYNPROGINDEX_MAJOR -1
#define DYNPROGINDEX_MINOR +1

#define DUAL_BREAK_PROB_THRESHOLD 0.90

#define THETA_SLACK 0.10
#define TRIM_END_PVALUE 1e-4

#define MAXPEELBACK_SCORE 5	/* For determining goodness of intron */
#define MAXPEELBACK_END 1000

static const Except_T gapcheck_error = {"Gap check failed"};
static const Except_T coordinate_error = {"Coordinate error"};

#define SHORTCUT 1		/* Skips re-solving introns if already canonical */

/* In debug mode, probably want to activate debug in pairpool.c and
   dynprog.c also */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Pair dump */
#ifdef DEBUG1
#define debug1(x) x
#else 
#define debug1(x)
#endif

/* Chimeras */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Before smoothing */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Show result of singles */
#ifdef DEBUG4
#define debug4(x) x
#else 
#define debug4(x)
#endif

/* Show result of removing adjacent dynamic programming */
#ifdef DEBUG5
#define debug5(x) x
#else 
#define debug5(x)
#endif

/* Show result of alignment before trimming of bad exons. */
#ifdef DEBUG6
#define debug6(x) x
#else 
#define debug6(x)
#endif

/* assign_gap_types and fill_in_gaps */
#ifdef DEBUG7
#define debug7(x) x
#else 
#define debug7(x)
#endif

/* changepoint */
#ifdef DEBUG8
#define debug8(x) x
#else 
#define debug8(x)
#endif


/************************************************************************
 *   Stage 3 merges cDNA-genomic pairs from stage 2 (called "path")
 *   and from dynamic programming into a single list.  In this
 *   process, stage 3 may also have to pop a pair off the path, insert
 *   the dynamic programming results (called "gappairs") and then push
 *   the stored pair onto the list.  The relevant pointers are
 *   leftquerypos and leftgenomepos, which refer to the left (stored) pair;
 *   rightquerypos and rightgenomepos, which refer to the right (top) pair on
 *   the list (which may represent what the list should have for the
 *   purposes of dynamic programming); and querydp5, genomedp5,
 *   querydp3, and genomedp3, which refer to the dynamic programming
 *   indices, inclusive.
 *
 *   path has the end of the query sequence as its car.
 *   pairs has the beginning of the query sequence as its car.
 *
 *   Most procedures take the top of path and put it onto pairs:
 *
 *	 <- <- path  =====>  pairs -> ->
 *	   leftpair	     rightpair
 *
 *   For stage 1 and stage 2, the poly-A/T tails, if any, was stripped
 *   off, but for stage 3, we try to extend these tails if possible.
 *   Therefore, we use the full length and full sequence here, and add
 *   the offset to the path from stage 2.
 ************************************************************************/

#define T Stage3_T
struct T {
  struct Pair_T *pairarray;	/* The array version of pairs_fwd or pairs_rev, with the gaps substituted */
  int npairs;

  List_T pairs;			/* Winning set of pairs */

  int straintype;
  char *strain;
  Chrnum_T chrnum;
  Genomicpos_T chroffset;	/* Position on genome of start of chromosome chrnum */
  Genomicpos_T chrpos;		/* Position on chromosome of start of genomic segment  */
  Genomicpos_T genomiclength;	/* Length of genomic segment */
  Genomicpos_T genomicstart;	/* Start of alignment */
  Genomicpos_T genomicend;	/* End of alignment */
  int cdna_direction;
  bool watsonp;
  double defect_rate;
  double trimmed_coverage;
  int matches;
  int unknowns;
  int mismatches;
  int qopens;
  int qindels;
  int topens;
  int tindels;
  int noncanonical;
  int goodness;

  int translation_start;
  int translation_end;
  int translation_length;

  int relaastart;
  int relaaend;

  /* Diagnostic info */
  Genomicpos_T stage1_genomicstart;
  Genomicpos_T stage1_genomiclength;
  double stage2_diag_runtime;
  double stage2_align_runtime;
  int stage2_indexsize;
  double stage2_mapfraction;
  int stage2_maxconsecutive;
  double stage3_runtime;
};


Genomicpos_T
Stage3_genomicstart (T this) {
  return this->genomicstart;
}


Genomicpos_T
Stage3_genomicend (T this) {
  return this->genomicend;
}

Genomicpos_T
Stage3_genomiclength (T this) {
  return this->genomiclength;
}



bool
Stage3_watsonp (T this) {
  return this->watsonp;
}

int
Stage3_cdna_direction (T this) {
  return this->cdna_direction;
}


int
Stage3_straintype (T this) {
  return this->straintype;
}

int
Stage3_goodness (T this) {
  debug2(printf("Overall goodness:\n"));
  debug2(printf("  %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels => %d\n",
		this->matches,this->mismatches,this->qopens,this->qindels,this->topens,this->tindels,this->goodness));

  return this->goodness;
}

struct Pair_T *
Stage3_pairarray (int *npairs, T this) {
  *npairs = this->npairs;
  return this->pairarray;
}

int
Stage3_matches (T this) {
  return this->matches;
}

int
Stage3_mismatches (T this) {
  return this->mismatches;
}

int
Stage3_indels (T this) {
  /* This should be consistent with the output from Pair_print_pathsummary */
  return this->qindels + this->tindels;
}


int
Stage3_querystart (T this) {
  return Pair_querypos(&(this->pairarray[0]));
}

int
Stage3_queryend (T this) {
  return Pair_querypos(&(this->pairarray[this->npairs-1]));
}
  
int
Stage3_translation_start (T this) {
  return this->translation_start;
}

int
Stage3_translation_end (T this) {
  return this->translation_end;
}


int
Stage3_domain (T this) {
  int querystart, queryend;

  querystart = Pair_querypos(&(this->pairarray[0]));
  queryend = Pair_querypos(&(this->pairarray[this->npairs-1]));

  return queryend - querystart + 1;
}


int
Stage3_largemargin (int *newstart, int *newend, T this, int queryntlength) {
  int leftmargin, rightmargin;
  int querystart, queryend;

  querystart = Pair_querypos(&(this->pairarray[0]));
  queryend = Pair_querypos(&(this->pairarray[this->npairs-1]));

  if ((leftmargin = querystart) < 0) {
    leftmargin = 0;
  }
  if ((rightmargin = queryntlength - queryend) < 0) {
    rightmargin = 0;
  }

  /* Return larger margin */
  *newstart = querystart;
  *newend = queryend;
  if (leftmargin > rightmargin) {
    /* Trim left */
    return leftmargin;
  } else {
    return rightmargin;
  }
}


double
Stage3_fracidentity (T this) {
  int den;

  if ((den = this->matches + this->mismatches + this->qindels + this->tindels) == 0) {
    return 1.0;
  } else {
    return (double) this->matches/(double) den;
  }
}

Genomicpos_T
Stage3_genomicpos (T this, int querypos, bool headp) {
  Genomicpos_T genomicpos;

  genomicpos = Pair_genomicpos(this->pairarray,this->npairs,querypos,headp);
  if (this->watsonp) {
    return this->chroffset + this->chrpos + genomicpos;
  } else {
    return this->chroffset + this->chrpos + (this->genomiclength - 1) - genomicpos;
  }
}


void
Stage3_pathscores (int *pathscores, T this, int querylength, cDNAEnd_T cdnaend) {
  Pair_pathscores(pathscores,this->pairarray,this->npairs,this->cdna_direction,querylength,cdnaend);
  return;
}


int
Stage3_chimeric_goodness (int *matches1, int *matches2, T part1, T part2, int breakpoint, int querylength) {
  int goodness1, goodness2, querystart, queryend;
  int unknowns1, mismatches1, qopens1, qindels1, topens1, tindels1, 
    ncanonical1, nsemicanonical1, nnoncanonical1;
  int unknowns2, mismatches2, qopens2, qindels2, topens2, tindels2, 
    ncanonical2, nsemicanonical2, nnoncanonical2;

  querystart = Pair_querypos(&(part1->pairarray[0]));
  debug2(printf("Chimeric goodness requested for part %d..%d\n",querystart+1,breakpoint));
  Pair_fracidentity_bounded(&(*matches1),&unknowns1,&mismatches1,&qopens1,&qindels1,&topens1,&tindels1,
			    &ncanonical1,&nsemicanonical1,&nnoncanonical1,
			    part1->pairarray,part1->npairs,part1->cdna_direction,
			    querystart,breakpoint);
  goodness1 = (*matches1) + MISMATCH*mismatches1 + QOPEN*qopens1 + QINDEL*qindels1 + TOPEN*topens1 + TINDEL*tindels1;
  debug2(printf("  %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels => %d\n",
		*matches1,mismatches1,qopens1,qindels1,topens1,tindels1,goodness1));

  queryend = Pair_querypos(&(part2->pairarray[part2->npairs-1]));
  debug2(printf("Chimeric goodness requested for part %d..%d\n",breakpoint+1,queryend+1));
  Pair_fracidentity_bounded(&(*matches2),&unknowns2,&mismatches2,&qopens2,&qindels2,&topens2,&tindels2,
			    &ncanonical2,&nsemicanonical2,&nnoncanonical2,
			    part2->pairarray,part2->npairs,part2->cdna_direction,
			    breakpoint,queryend);
  goodness2 = (*matches2) + MISMATCH*mismatches2 + QOPEN*qopens2 + QINDEL*qindels2 + TOPEN*topens2 + TINDEL*tindels2;
  debug2(printf("  %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels => %d\n",
		*matches2,mismatches2,qopens2,qindels2,topens2,tindels2,goodness2));

  return goodness1 + goodness2;
}


int
Stage3_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->goodness > y->goodness) {
    return -1;
  } else if (x->goodness < y->goodness) {
    return +1;
  } else if (x->straintype < y->straintype) {
    return -1;
  } else if (x->straintype > y->straintype) {
    return +1;
  } else if (x->chrnum < y->chrnum) {
    return -1;
  } else if (x->chrnum > y->chrnum) {
    return +1;
  } else if (x->chrpos < y->chrpos) {
    return -1;
  } else if (x->chrpos > y->chrpos) {
    return +1;
  } else {
    return 0;
  }
}

bool
Stage3_overlap (T x, T y) {

  if (x->straintype != y->straintype) {
    return false;
  } else if (x->watsonp != y->watsonp) {
    return false;
  } else if (x->watsonp) {
    if (x->genomicstart >= y->genomicstart && x->genomicstart <= y->genomicend) {
      return true;
    } else if (y->genomicstart >= x->genomicstart && y->genomicstart <= x->genomicend) {
      return true;
    } else {
      return false;
    }
  } else {
    if (x->genomicstart >= y->genomicend && x->genomicstart <= y->genomicstart) {
      return true;
    } else if (y->genomicstart >= x->genomicend && y->genomicstart <= x->genomicstart) {
      return true;
    } else {
      return false;
    }
  }
}

/************************************************************************
 *   Gaps
 ************************************************************************/

static List_T
check_gaps (List_T pairs, Pairpool_T pairpool) {
  List_T path = NULL;
  Pair_T pair, leftpair, rightpair;
  int queryjump, genomejump;

  debug(printf("\nBeginning check of gaps\n"));
  debug(printf("length = %d\n",List_length(pairs)));
  debug(Pair_dump_list(pairs,true));

  pairs = Pairpool_pop(pairs,&pair);
  if (pair->gapp == true) {
    fprintf(stderr,"Gap check error: Unexpected gap at start of pairs\n");
    debug(printf("Gap check error: Unexpected gap at start of pairs\n"));
#ifndef DEBUG
    Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
  } else {
    path = Pairpool_push_existing(NULL,pairpool,pair);
  }

  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == true) {
      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Observed a gap at %d..%d with queryjump = %d, genomejump = %d\n",
		   leftpair->querypos,rightpair->querypos,pair->queryjump,pair->genomejump));

      queryjump = rightpair->querypos - leftpair->querypos - 1;
      genomejump = rightpair->genomepos - leftpair->genomepos - 1;
      if (leftpair->cdna == ' ') queryjump++;
      if (leftpair->genome == ' ') genomejump++;

      if (pair->queryjump != queryjump) {
	if (rightpair->querypos >= HALFLEN && leftpair->querypos < HALFLEN) {
	  debug(printf("Accept queryjump for gap at %d..%d as probable skiplength.  It's %d, should be %d\n",
		       leftpair->querypos,rightpair->querypos,pair->queryjump,queryjump));
	} else {
	  debug(printf("Gap check error: Wrong queryjump for gap at %d..%d.  It's %d, should be %d\n",
		       leftpair->querypos,rightpair->querypos,pair->queryjump,queryjump));
#ifndef DEBUG
	  Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
	}
      }
      if (pair->genomejump != genomejump) {
	debug(printf("Gap check error: Wrong genomejump for gap at %d..%d.  It's %d, should be %d\n",
		     leftpair->querypos,rightpair->querypos,pair->genomejump,genomejump));
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
      }
      path = Pairpool_push_existing(path,pairpool,pair);

      /* Process another pair after gap */
      if (pairs == NULL) {
	fprintf(stderr,"Gap check error: Unexpected gap at end of pairs\n");
	debug(printf("Gap check error: Unexpected gap at end of pairs\n"));
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
      }
      pairs = Pairpool_pop(pairs,&pair);
      if (pair->gapp == true) {
	fprintf(stderr,"Gap check error: Unexpected gap after gap\n");
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
      }
      path = Pairpool_push_existing(path,pairpool,pair);

    } else {
      /* Not a gap */
      leftpair = path->first;
      queryjump = pair->querypos - leftpair->querypos - 1;
      genomejump = pair->genomepos - leftpair->genomepos - 1;
      if (leftpair->cdna == ' ') queryjump++;
      if (leftpair->genome == ' ') genomejump++;
	  
      if (queryjump <= 0 && genomejump <= 0) {
	path = Pairpool_push_existing(path,pairpool,pair);
      } else if (queryjump == 0 && genomejump == 0) {
	path = Pairpool_push_existing(path,pairpool,pair);
      } else {
	fprintf(stderr,"Gap check error: Unexpected missing gap at %d..%d\n",leftpair->querypos,pair->querypos);
	debug(printf("Gap check error: Unexpected missing gap at %d..%d\n",leftpair->querypos,pair->querypos));
	debug(printf("Gap check error: Pushing a gap at %d..%d because of queryjump = %d, genomejump = %d\n",
		     leftpair->querypos,pair->querypos,queryjump,genomejump));
	/* One place we need accurate queryjump and genomejump */
	path = Pairpool_push_gapholder(path,pairpool,queryjump,genomejump);
	path = Pairpool_push_existing(path,pairpool,pair);
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
#endif
      }
    }
  }	

  debug(printf("Done with check of gaps\n\n"));

  return path;
}

#if 0
static int
count_gaps (int *length5, int *length3, List_T pairs) {
  int ngaps = 0;
  List_T p;
  Pair_T pair, leftpair;
  int queryjump, genomejump;

  /* Remove all existing gaps */
  debug(printf("Beginning counting of gaps\n"));
  *length5 = *length3 = 0;

  for (p = pairs; p != NULL; p = p->rest) {
    pair = (Pair_T) p->first;
    if (p == pairs) {
      leftpair = pair;
    } else {
      queryjump = pair->querypos - leftpair->querypos - 1;
      genomejump = pair->genomepos - leftpair->genomepos - 1;
      if (leftpair->cdna == ' ') queryjump++;
      if (leftpair->genome == ' ') genomejump++;

      if (queryjump <= 0 && genomejump <= 0) {
	/* Not a gap */
	if (ngaps == 0) {
	  *length5 += 1;
	} else {
	  *length3 += 1;
	}
      } else {
	ngaps++;
	*length3 = 0;
      }
    }
    leftpair = pair;
  }

  debug(printf("Counted %d gaps.  length5 = %d, length3 = %d\n",ngaps,*length5,*length3));

  return ngaps;
}
#endif


static List_T
insert_gapholders (List_T pairs, Pairpool_T pairpool) {
  List_T path = NULL;
  Pair_T pair, leftpair;
  int queryjump, genomejump;

  /* Remove all existing gaps */
  debug(printf("Beginning deletion/insertion of gaps\n"));

  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    while (pairs != NULL && pair->gapp == true) {
      /* Discard old gap(s) */
      debug(printf("Removing a gap with queryjump = %d, genomejump = %d\n",
		   pair->queryjump,pair->genomejump));
      pairs = Pairpool_pop(pairs,&pair);
    }

    if (path == NULL) {
      path = Pairpool_push_existing(path,pairpool,pair);
      leftpair = pair;
    } else {
      queryjump = pair->querypos - leftpair->querypos - 1;
      genomejump = pair->genomepos - leftpair->genomepos - 1;
      if (leftpair->cdna == ' ') queryjump++;
      if (leftpair->genome == ' ') genomejump++;

      if (queryjump <= 0 && genomejump <= 0) {
	path = Pairpool_push_existing(path,pairpool,pair);
      } else {
	/* Insert new gap.  Need accurate queryjump and genomejump */
	debug(printf("Inserting a gap at %d..%d because of queryjump = %d, genomejump = %d\n",
		     leftpair->querypos,pair->querypos,queryjump,genomejump));
	path = Pairpool_push_gapholder(path,pairpool,queryjump,genomejump);
	path = Pairpool_push_existing(path,pairpool,pair);
      }
    }
    leftpair = pair;
  }

  debug(printf("Ending deletion/insertion of gaps\n"));

  return path;
}

#if 0
static List_T
insert_gapholders_and_singles_old (List_T path, Pairpool_T pairpool,
				   char *queryseq_ptr, char *queryuc_ptr,
				   char *genomicseg_ptr, char *genomicuc_ptr,
				   int skiplength) {
  List_T pairs = NULL;
  Pair_T pair;
  int leftquerypos, leftgenomepos, rightquerypos, rightgenomepos, queryjump, genomejump;
  char c1, c2;

  debug(printf("\nBeginning insertion of gaps and singles\n"));

  if (path == NULL) {
    return NULL;
  }

  path = Pairpool_pop(path,&pair);
  pairs = Pairpool_push_existing(NULL,pairpool,pair);
  rightquerypos = pair->querypos;
  rightgenomepos = pair->genomepos;

  Pair_dump_list(path,true);

  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == true) {
      debug(printf("Observed a gap at %d..%d with queryjump = %d, genomejump = %d\n",
		   leftquerypos,rightquerypos,pair->queryjump,pair->genomejump));
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      if (path == NULL) {
	fprintf(stderr,"Unexpected gap at beginning of query\n");
	abort();
      }
      path = Pairpool_pop(path,&pair);
      if (pair->gapp == true) {
	printf("Unexpected gap after gap\n");
	fprintf(stderr,"Unexpected gap after gap\n");
	abort();
      }
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
    } else {
      leftquerypos = pair->querypos;
      leftgenomepos = pair->genomepos;

      queryjump = rightquerypos - leftquerypos - 1;
      genomejump = rightgenomepos - leftgenomepos - 1;
      if (pair->cdna == ' ') queryjump++;
      if (pair->genome == ' ') genomejump++;

      if (skiplength > 0 && rightquerypos >= HALFLEN && leftquerypos < HALFLEN) {
	debug(printf("Introducing skiplength %d at %d..%d\n",skiplength,leftquerypos,rightquerypos));
	queryjump += skiplength;
      }
	  
      if (queryjump < 0 || genomejump < 0) {
	fprintf(stderr,"Unexpected negative queryjump or genomejump\n");
	printf("Unexpected negative queryjump %d or genomejump %d\n",queryjump,genomejump);
	Pair_dump_one(pair,true);
	/* abort() */

      } else if (queryjump == 0 && genomejump == 0) {
	pairs = Pairpool_push_existing(pairs,pairpool,pair);

      } else if (queryjump == 1 && genomejump == 1) {
	/* Special case: Single mismatch */
	c1 = queryseq_ptr[leftquerypos+1];
	c2 = genomicseg_ptr[leftgenomepos+1];
	if (queryuc_ptr[leftquerypos+1] == genomicuc_ptr[leftgenomepos+1]) {
	  /* Possible for a stretch of the same nucleotide to have a mismatch */
	  debug(printf("Stage 3: Single match at %d %d: %c | %c\n",
		       leftquerypos+1,leftgenomepos+1,c1,c2));
	  pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
				   c1,MATCH_COMP,c2,/*dynprogindex*/0);
	} else {
	  debug(printf("Stage 3: Single mismatch at %d %d: %c   %c\n",
		       leftquerypos+1,leftgenomepos+1,c1,c2));
	  pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
				   c1,MISMATCH_COMP,c2,/*dynprogindex*/0);
	}
	pairs = Pairpool_push_existing(pairs,pairpool,pair);

      } else if (queryjump == 1 && genomejump == 0) {
	/* Special case: Single cDNA insert.  This will create an apparent genomejump of -1. */
	debug(printf("Stage 3: Single cDNA insert at %d %d: %c\n",
		     leftquerypos+1,leftgenomepos+1,queryseq_ptr[leftquerypos+1]));
	/* Advance querypos to next position */
	pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
				 queryseq_ptr[leftquerypos+1],INDEL_COMP,' ',/*dynprogindex*/0);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);

#if 0
      } else if (queryjump == 0 && genomejump == 1) {
	/* Special case: Single genomic insert */
	debug(printf("Stage 3: Single genomic insert at %d %d: %c\n",
		     leftquerypos+1,leftgenomepos+1,genomicseg_ptr[leftgenomepos+1]));
	/* Advance genomepos to next position */
	pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
				 ' ',INDEL_COMP,genomicseg_ptr[leftgenomepos+1],/*dynprogindex*/0);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#endif

      } else {
	debug(printf("Pushing a gap at %d..%d because of queryjump = %d, genomejump = %d\n",
		     leftquerypos,rightquerypos,queryjump,genomejump));
	/* One place we need accurate queryjump and genomejump */
	pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }
    }

    rightquerypos = leftquerypos;
    rightgenomepos = leftgenomepos;
  }	

  debug(printf("Done with insertion of gaps and singles\n\n"));

  return pairs;
}
#endif


#if 0
static List_T
insert_singles (List_T path, Pairpool_T pairpool,
		char *queryseq_ptr, char *queryuc_ptr,
		char *genomicseg_ptr, char *genomicuc_ptr,
		int skiplength) {
  List_T pairs = NULL;
  Pair_T pair;
  int leftquerypos, leftgenomepos, rightquerypos, rightgenomepos, queryjump, genomejump;
  char c1, c2;

  debug(printf("\nBeginning insertion of gaps and singles\n"));

  if (path == NULL) {
    return NULL;
  }

  path = Pairpool_pop(path,&pair);
  pairs = Pairpool_push_existing(NULL,pairpool,pair);
  rightquerypos = pair->querypos;
  rightgenomepos = pair->genomepos;

  Pair_dump_list(path,true);

  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == true) {
      debug(printf("Observed a gap at %d..%d with queryjump = %d, genomejump = %d\n",
		   leftquerypos,rightquerypos,pair->queryjump,pair->genomejump));
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      if (path == NULL) {
	fprintf(stderr,"Unexpected gap at beginning of query\n");
	abort();
      }
      path = Pairpool_pop(path,&pair);
      if (pair->gapp == true) {
	printf("Unexpected gap after gap\n");
	fprintf(stderr,"Unexpected gap after gap\n");
	abort();
      }
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
    } else {
      leftquerypos = pair->querypos;
      leftgenomepos = pair->genomepos;

      queryjump = rightquerypos - leftquerypos - 1;
      genomejump = rightgenomepos - leftgenomepos - 1;
      if (pair->cdna == ' ') queryjump++;
      if (pair->genome == ' ') genomejump++;

      if (skiplength > 0 && rightquerypos >= HALFLEN && leftquerypos < HALFLEN) {
	debug(printf("Introducing skiplength %d at %d..%d\n",skiplength,leftquerypos,rightquerypos));
	queryjump += skiplength;
      }
	  
      if (queryjump < 0 || genomejump < 0) {
	fprintf(stderr,"Unexpected negative queryjump or genomejump\n");
	printf("Unexpected negative queryjump %d or genomejump %d\n",queryjump,genomejump);
	Pair_dump_one(pair,true);
	/* abort() */

      } else if (queryjump == 0 && genomejump == 0) {
	pairs = Pairpool_push_existing(pairs,pairpool,pair);

      } else if (queryjump == 1 && genomejump == 1) {
	/* Special case: Single mismatch */
	c1 = queryseq_ptr[leftquerypos+1];
	c2 = genomicseg_ptr[leftgenomepos+1];
	if (queryuc_ptr[leftquerypos+1] == genomicuc_ptr[leftgenomepos+1]) {
	  /* Possible for a stretch of the same nucleotide to have a mismatch */
	  debug(printf("Stage 3: Single match at %d %d: %c | %c\n",
		       leftquerypos+1,leftgenomepos+1,c1,c2));
	  pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
				   c1,MATCH_COMP,c2,/*dynprogindex*/0);
	} else {
	  debug(printf("Stage 3: Single mismatch at %d %d: %c   %c\n",
		       leftquerypos+1,leftgenomepos+1,c1,c2));
	  pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
				   c1,MISMATCH_COMP,c2,/*dynprogindex*/0);
	}
	pairs = Pairpool_push_existing(pairs,pairpool,pair);

      } else if (queryjump == 1 && genomejump == 0) {
	/* Special case: Single cDNA insert.  This will create an apparent genomejump of -1. */
	debug(printf("Stage 3: Single cDNA insert at %d %d: %c\n",
		     leftquerypos+1,leftgenomepos+1,queryseq_ptr[leftquerypos+1]));
	/* Advance querypos to next position */
	pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
				 queryseq_ptr[leftquerypos+1],INDEL_COMP,' ',/*dynprogindex*/0);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);

#if 0
      } else if (queryjump == 0 && genomejump == 1) {
	/* Special case: Single genomic insert */
	debug(printf("Stage 3: Single genomic insert at %d %d: %c\n",
		     leftquerypos+1,leftgenomepos+1,genomicseg_ptr[leftgenomepos+1]));
	/* Advance genomepos to next position */
	pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
				 ' ',INDEL_COMP,genomicseg_ptr[leftgenomepos+1],/*dynprogindex*/0);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#endif

      } else {
	debug(printf("Pushing a gap at %d..%d because of queryjump = %d, genomejump = %d\n",
		     leftquerypos,rightquerypos,queryjump,genomejump));
	/* One place we need accurate queryjump and genomejump */
	pairs = Pairpool_push_gapholder(pairs,pairpool,queryjump,genomejump);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }
    }

    rightquerypos = leftquerypos;
    rightgenomepos = leftgenomepos;
  }	

  debug(printf("Done with insertion of gaps and singles\n\n"));

  return pairs;
}
#endif



static List_T
assign_gap_types (int *npairs, int *ngaps, int *ncanonical, 
		  List_T path, Pairpool_T pairpool, char *queryseq_ptr, char *genomicseg_ptr, 
		  char *genomicuc_ptr, int cdna_direction) {
  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;
  int queryjump, genomejump, leftquerypos, leftgenomepos, rightquerypos, rightgenomepos, curquerypos,
    introntype, intronlength, genomicpos;
  char left1, left2, right2, right1, c2;

  debug(printf("\n** Starting assign_gap_types\n"));
  *npairs = 0;
  *ngaps = 0;
  *ncanonical = 0;
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == false) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      (*npairs)++;
    } else if (pairs == NULL) {
      /* Discard initial gap */
    } else {
      queryjump = pair->queryjump;
      genomejump = pair->genomejump;

      if (queryjump == 0 && genomejump == 0) {
	debug7(printf("  Gap is a non-gap\n"));
	/* Discard the gap pair */

      } else if (genomejump == 0) {
	debug7(printf("  Gap is a cDNA insertion\n"));
	/* pair->comp = INDEL_COMP; */

	leftpair = path->first;
	rightpair = pairs->first;
	leftquerypos = leftpair->querypos;
	if (leftpair->cdna == ' ') leftquerypos--;
	rightquerypos = rightpair->querypos;
	rightgenomepos = rightpair->genomepos;

	for (curquerypos = rightquerypos - 1; curquerypos > leftquerypos; --curquerypos) {
	  pairs = Pairpool_push(pairs,pairpool,curquerypos,rightgenomepos,
				queryseq_ptr[curquerypos],INDEL_COMP,' ',/*dynprogindex*/0);
	  (*npairs)++;
	}
	/* Discard the gap pair */

      } else if (queryjump > 0) {
	debug7(printf("  Gap is a dual break\n"));
	pair->comp = DUALBREAK_COMP;
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
	(*npairs)++;
	(*ngaps)++;

      } else {
	debug7(printf("  Gap is an intron\n"));

	leftpair = path->first;
	rightpair = pairs->first;

	leftquerypos = leftpair->querypos;
	leftgenomepos = leftpair->genomepos;
	if (leftpair->cdna == ' ') leftquerypos--;
	if (leftpair->genome == ' ') leftgenomepos--;
	rightquerypos = rightpair->querypos;
	rightgenomepos = rightpair->genomepos;

	intronlength = rightgenomepos - leftgenomepos - 1;
	if (intronlength < MIN_NONINTRON) {
	  debug7(printf("  Gap is short (intronlength %d).  Adding pairs from %d downto %d\n",
			intronlength,rightgenomepos-1,leftgenomepos+1));
	  for (genomicpos = rightgenomepos - 1; genomicpos > leftgenomepos; --genomicpos) {
	    c2 = genomicseg_ptr[genomicpos];
	    pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',/*comp*/SHORTGAP_COMP,c2,
				  /*dynprogindex*/0);
	    (*npairs)++;
	  }
	  debug7(printf("  Gap is a short gap, so discarding the gap pair\n"));
	  /* Discard the gap pair */

	} else {
	  debug7(printf("  Gap is not short (intronlength %d)\n",intronlength));
	  left1 = genomicuc_ptr[leftgenomepos+1];
	  left2 = genomicuc_ptr[leftgenomepos+2];
	  right2 = genomicuc_ptr[rightgenomepos-2];
	  right1 = genomicuc_ptr[rightgenomepos-1];
	  
	  debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
	  introntype = Intron_type(left1,left2,right2,right1,cdna_direction);
	  debug7(printf("  Introntype is %d (cdna_direction %d)\n",introntype,cdna_direction));
	  if (cdna_direction >= 0) {
	    switch (introntype) {
	    case GTAG_FWD: pair->comp = FWD_CANONICAL_INTRON_COMP; (*ncanonical)++; break;
	    case GCAG_FWD: pair->comp = FWD_GCAG_INTRON_COMP; break;
	    case ATAC_FWD: pair->comp = FWD_ATAC_INTRON_COMP; break;
	    case NONINTRON: pair->comp = NONINTRON_COMP; break;
	    default: 
	      printf("Unexpected intron type %d\n",introntype);
	      fprintf(stderr,"Unexpected intron type %d\n",introntype);
	      exit(9);
	    }
	  } else {
	    switch (introntype) {
	    case ATAC_REV: pair->comp = REV_ATAC_INTRON_COMP; break;
	    case GCAG_REV: pair->comp = REV_GCAG_INTRON_COMP; break;
	    case GTAG_REV: pair->comp = REV_CANONICAL_INTRON_COMP; (*ncanonical)++; break;
	    case NONINTRON: pair->comp = NONINTRON_COMP; break;
	    default: 
	      printf("Unexpected intron type %d\n",introntype);
	      fprintf(stderr,"Unexpected intron type %d\n",introntype);
	      exit(9);
	    }
	  }
	  /* Push the gap back on */
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
	  (*npairs)++;
	  (*ngaps)++;
	}
      }
    }
  }

  return pairs;
}


#if 0
/* This procedure removes gaps and returns original pairs, so need to
   call insert_gapholders after this procedure */
static List_T
remove_adjacent_dynprog (bool *removep, List_T pairs, Pairpool_T pairpool) {
  List_T path = NULL, buffer = NULL;
  Pair_T pair = NULL;
  int initial_dynprogindex;
  bool recomputep = false, dynprogp = false;
#ifdef DEBUG
  Pair_T lastpair, rightpair;
  int leftquerypos, rightquerypos, leftgenomepos, rightgenomepos;
  int queryjump, genomejump;
#endif

  debug(printf("\nBeginning removal of adjacent dynamic programming\n"));
  *removep = false;

  while (pairs != NULL) {
    debug(lastpair = pair);
    pairs = Pairpool_pop(pairs,&pair);

    if (pair->dynprogindex > 0) {
      if (dynprogp == true) {
	/* Existing region */
	buffer = Pairpool_push_existing(buffer,pairpool,pair);
	if (pair->dynprogindex != initial_dynprogindex) {
	  debug(if (recomputep == false) {
		  printf("Adjacent dynprog #%d at pair %d\n",pair->dynprogindex,pair->querypos);
		});
	  recomputep = true;
	}
      } else {
	/* New region */
	debug(printf("Start of dynprog #%d at pair %d\n",pair->dynprogindex,pair->querypos));
	debug(leftquerypos = pair->querypos);
	debug(leftgenomepos = pair->genomepos);
	buffer = Pairpool_push_existing(NULL,pairpool,pair);
	initial_dynprogindex = pair->dynprogindex;
	recomputep = false;
	dynprogp = true;
      }

    } else if (dynprogp == false) {
      if (pair->gapp == false) {
	path = Pairpool_push_existing(path,pairpool,pair);
      }

    } else if (dynprogp == true) {
      if (recomputep == false) {
	path = Pairpool_transfer(path,List_reverse(buffer));
	debug(rightpair = lastpair;
	      rightquerypos = rightpair->querypos;
	      printf("End of dynprog at pair %d\n",rightquerypos)
	      );

      } else {
	debug(
	      rightpair = lastpair;
	      rightquerypos = rightpair->querypos;
	      rightgenomepos = rightpair->genomepos;
	      printf("End of dynprog at pair %d\n",rightquerypos);
	      );
	
	debug(
	      queryjump = rightquerypos - leftquerypos + 1;
	      genomejump = rightgenomepos - leftgenomepos + 1;
	      if (pair->cdna == ' ') queryjump++;
	      if (pair->genome == ' ') genomejump++;
	      printf("Observed adjacent dynamic programming at %d..%d with queryjump = %d, genomejump = %d\n",
		     leftquerypos,rightquerypos,queryjump,genomejump);
	      Pair_dump_list(buffer,true);
	      );
	
	/* Gaps will be inserted later */
	/* path = Pairpool_push_gapholder(path,pairpool,queryjump,genomejump); */
	*removep = true;
      }
      if (pair->gapp == false) {
	path = Pairpool_push_existing(path,pairpool,pair);
      }
      dynprogp = false;
    }
  }

  if (dynprogp == true) {
    if (recomputep == false) {
      path = Pairpool_transfer(path,List_reverse(buffer));
      debug(rightpair = lastpair;
	    rightquerypos = rightpair->querypos;
	    printf("End of dynprog at pair %d\n",rightquerypos)
	    );
      
    } else {
      debug(
	    rightpair = lastpair;
	    rightquerypos = rightpair->querypos;
	    rightgenomepos = rightpair->genomepos;
	    printf("End of dynprog at pair %d\n",rightquerypos);
	    );
	
      debug(
	    queryjump = rightquerypos - leftquerypos + 1;
	    genomejump = rightgenomepos - leftgenomepos + 1;
	    if (pair->cdna == ' ') queryjump++;
	    if (pair->genome == ' ') genomejump++;
	    printf("Observed adjacent dynamic programming at %d..%d with queryjump = %d, genomejump = %d\n",
		   leftquerypos,rightquerypos,queryjump,genomejump);
	    Pair_dump_list(buffer,true);
	    );
	
      /* Gaps will be inserted later */
      /* path = Pairpool_push_gapholder(path,pairpool,queryjump,genomejump); */
      *removep = true;
    }
  }

  debug(printf("Done with removal of adjacent dynamic programming\n\n"));
  
  return List_reverse(path);
}
#endif



#ifdef PMAP
static List_T
undefine_nucleotides (char *queryseq_ptr, int querylength, List_T path, Pairpool_T pairpool, int width) {
  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;
  int leftquerypos, leftgenomepos, rightquerypos, rightgenomepos, pos;

  debug(printf("\n** Starting undefine_nucleotides\n"));

  if (path != NULL) {
    path = Pairpool_pop(path,&pair);
    pairs = Pairpool_push_existing(NULL,pairpool,pair);
    rightquerypos = pair->querypos;
    rightgenomepos = pair->genomepos;
  }

  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == true) {
      leftpair = path->first;
      rightpair = pairs->first;

      leftquerypos = leftpair->querypos;
      leftgenomepos = leftpair->genomepos;
      if (leftpair->cdna == ' ') leftquerypos--;
      if (leftpair->genome == ' ') leftgenomepos--;

      rightquerypos = rightpair->querypos;
      rightgenomepos = rightpair->genomepos;
    
      debug(printf("Undefining around rightquerypos = %d and leftquerypos = %d\n",rightquerypos,leftquerypos));
      for (pos = rightquerypos; pos < rightquerypos + width && pos < querylength; pos++) {
	queryseq_ptr[pos] = BACKTRANSLATE_CHAR;
      }
      for (pos = leftquerypos; pos > leftquerypos - width && pos >= 0; --pos) {
	queryseq_ptr[pos] = BACKTRANSLATE_CHAR;
      }
    }
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
  }

  return pairs;
}  
#endif


static List_T
add_dualbreak (List_T pairs, char *queryseq_ptr, 
#ifdef PMAP
	       char *queryaaseq_ptr,
#endif
	       char *genomicseg_ptr, char *genomicuc_ptr, int cdna_direction,
	       Pair_T leftpair, Pair_T rightpair, Pairpool_T pairpool, int ngap,
	       int diagnosticp) {
  int genomicpos, k;
  int leftquerypos, leftgenomepos, rightquerypos, rightgenomepos, gapgenomepos;
  int introntype;
  char left1, left2, right2, right1, c1, c2, comp;

  leftquerypos = leftpair->querypos;
  leftgenomepos = leftpair->genomepos;
  if (leftpair->cdna == ' ') leftquerypos--;
  if (leftpair->genome == ' ') leftgenomepos--;
  rightquerypos = rightpair->querypos;
  rightgenomepos = rightpair->genomepos;

  if (genomicuc_ptr != NULL) {
    left1 = genomicuc_ptr[leftgenomepos+1];
    left2 = genomicuc_ptr[leftgenomepos+2];
    right2 = genomicuc_ptr[rightgenomepos-2];
    right1 = genomicuc_ptr[rightgenomepos-1];
	  
    debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
    introntype = Intron_type(left1,left2,right2,right1,cdna_direction);
    switch (introntype) {
    case GTAG_FWD: comp = FWD_CANONICAL_INTRON_COMP; break;
    case GCAG_FWD: comp = FWD_GCAG_INTRON_COMP; break;
    case ATAC_FWD: comp = FWD_ATAC_INTRON_COMP; break;
    case ATAC_REV: comp = REV_ATAC_INTRON_COMP; break;
    case GCAG_REV: comp = REV_GCAG_INTRON_COMP; break;
    case GTAG_REV: comp = REV_CANONICAL_INTRON_COMP; break;
    case NONINTRON: comp = NONINTRON_COMP; break;
    default: 
      printf("Unexpected intron type %d\n",introntype);
      fprintf(stderr,"Unexpected intron type %d\n",introntype);
      exit(9);
    }
  }

  /* queryjump = rightquerypos - leftquerypos - 1; */
  /* genomejump = rightgenomepos - leftgenomepos - 1; */

  if (diagnosticp == true) {
    for (k = 0; k < ngap; k++) {
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,rightgenomepos,' ',DUALBREAK_COMP,' ',
				     /*extraexonp*/false);
    }
    for (k = 0; k < 3; k++) {
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,rightgenomepos,' ',INTRONGAP_COMP,' ',
				     /*extraexonp*/false);
    }
    for (k = 0; k < ngap; k++) {
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,rightgenomepos,' ',DUALBREAK_COMP,' ',
				     /*extraexonp*/false);
    }

  } else {

    /* First insertion */
    for (k = 0, genomicpos = rightgenomepos - 1; k < ngap; k++, --genomicpos) {
      c2 = genomicseg_ptr[genomicpos];
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,genomicpos,' ',comp,c2,/*extraexonp*/true);
    }

    /* cDNA sequence */
    gapgenomepos = genomicpos + 1;
    for (k = rightquerypos - 1; k > leftquerypos; --k) {
#ifdef PMAP
      c1 = Dynprog_codon_char(queryaaseq_ptr[k/3],k%3);
#else
      c1 = queryseq_ptr[k];
#endif
      pairs = Pairpool_push_gapalign(pairs,pairpool,k,gapgenomepos,c1,EXTRAEXON_COMP,c1,/*extraexonp*/true); /* Transfer cDNA char to genome */
    }

    /* Second insertion */
    genomicpos = leftgenomepos + ngap;
    for (k = 0; k < ngap; k++, --genomicpos) {
      c2 = genomicseg_ptr[genomicpos];
      pairs = Pairpool_push_gapalign(pairs,pairpool,leftquerypos,genomicpos,' ',comp,c2,/*extraexonp*/true);
    }
  }

  return pairs;
}


static List_T
add_intron (List_T pairs, char *queryseq_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
	    Pair_T leftpair, Pair_T rightpair,  char *comp, int ngap, int cdna_direction,
	    Pairpool_T pairpool) {
  char c2;
  int intronlength, genomicpos;
  int leftgenomepos, rightquerypos, rightgenomepos, gapgenomepos;
  int i;

  leftgenomepos = leftpair->genomepos;
  if (leftpair->genome == ' ') leftgenomepos--;
  rightquerypos = rightpair->querypos;
  rightgenomepos = rightpair->genomepos;

  intronlength = rightgenomepos - leftgenomepos - 1;

  debug7(printf("Adding gap of type %c of length %d\n",*comp,intronlength));

#if 0
  /* Should not be necessary to fix introns at this point */
  if (cdna_direction >= 0) {
    switch (*comp) {
    case FWD_CANONICAL_INTRON_COMP: case FWD_GCAG_INTRON_COMP: case FWD_ATAC_INTRON_COMP: case NONINTRON: break;
    default: 
      debug7(printf("Unexpected intron comp %c.  Need to fix.\n",*comp));

      left1 = genomicuc_ptr[leftgenomepos+1];
      left2 = genomicuc_ptr[leftgenomepos+2];
      right2 = genomicuc_ptr[rightgenomepos-2];
      right1 = genomicuc_ptr[rightgenomepos-1];
	  
      debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
      introntype = Intron_type(left1,left2,right2,right1,cdna_direction);
      debug7(printf("  Introntype is %d (cdna_direction %d)\n",introntype,cdna_direction));
      switch (introntype) {
      case GTAG_FWD: *comp = FWD_CANONICAL_INTRON_COMP; break;
      case GCAG_FWD: *comp = FWD_GCAG_INTRON_COMP; break;
      case ATAC_FWD: *comp = FWD_ATAC_INTRON_COMP; break;
      case NONINTRON:
	intronlength = rightgenomepos - leftgenomepos - 1;
	if (intronlength < MIN_NONINTRON) {
	  *comp = SHORTGAP_COMP;	/* Will be printed as INDEL_COMP, but need to score as NONINTRON_COMP */
	} else {
	  *comp = NONINTRON_COMP;
	}
      }
    }
  } else {
    switch (*comp) {
    case REV_CANONICAL_INTRON_COMP: case REV_GCAG_INTRON_COMP: case REV_ATAC_INTRON_COMP: case NONINTRON: break;
    default: 
      debug7(printf("Unexpected intron comp %c.  Need to fix.\n",*comp));

      left1 = genomicuc_ptr[leftgenomepos+1];
      left2 = genomicuc_ptr[leftgenomepos+2];
      right2 = genomicuc_ptr[rightgenomepos-2];
      right1 = genomicuc_ptr[rightgenomepos-1];
	  
      debug7(printf("  Dinucleotides are %c%c..%c%c\n",left1,left2,right2,right1));
      introntype = Intron_type(left1,left2,right2,right1,cdna_direction);
      debug7(printf("  Introntype is %d (cdna_direction %d)\n",introntype,cdna_direction));
      switch (introntype) {
      case ATAC_REV: *comp = REV_ATAC_INTRON_COMP; break;
      case GCAG_REV: *comp = REV_GCAG_INTRON_COMP; break;
      case GTAG_REV: *comp = REV_CANONICAL_INTRON_COMP; break;
      default:
	intronlength = rightgenomepos - leftgenomepos - 1;
	if (intronlength < MIN_NONINTRON) {
	  *comp = SHORTGAP_COMP;	/* Will be printed as INDEL_COMP, but need to score as NONINTRON_COMP */
	} else {
	  *comp = NONINTRON_COMP;
	}
      }
    }
  }
#endif

  if (intronlength < ngap + ngap + 3) {
    for (i = 0, genomicpos = rightgenomepos - 1; i < intronlength; i++, --genomicpos) {
      c2 = genomicseg_ptr[genomicpos];
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,genomicpos,' ',*comp,c2,/*extraexonp*/false);
    }
  } else {
    for (i = 0, genomicpos = rightgenomepos - 1; i < ngap; i++, --genomicpos) {
      c2 = genomicseg_ptr[genomicpos];
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,genomicpos,' ',*comp,c2,/*extraexonp*/false);
      debug7(printf("Pushing %c at genomicpos %d\n",c2,genomicpos));
    }
    
    gapgenomepos = genomicpos + 1;
    pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,gapgenomepos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,
				   /*extraexonp*/false);
    pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,gapgenomepos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,
				   /*extraexonp*/false);
    pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,gapgenomepos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,
				   /*extraexonp*/false);

    genomicpos = leftgenomepos + ngap;
    for (i = ngap-1; i >= 0; --i, --genomicpos) {
      c2 = genomicseg_ptr[genomicpos];
      pairs = Pairpool_push_gapalign(pairs,pairpool,rightquerypos,genomicpos,' ',*comp,c2,/*extraexonp*/false);
      debug7(printf("Pushing %c at genomicpos %d\n",c2,genomicpos));
    }
  }

  return pairs;
}


/************************************************************************
 *   Chop (trimming within end exons)
 ************************************************************************/

/* Called only by GMAP, because nucleotide matches in PMAP have several ambiguous matches. */
static List_T
clean_path_end3 (List_T path) {
  Pair_T lastpair;

  debug(printf("clean_path_end3\n"));
  /* Remove any remaining nonmatches at 3' end, which can happen rarely */
  if (path != NULL) {
    lastpair = path->first;
    while (lastpair->comp != MATCH_COMP && lastpair->comp != DYNPROG_MATCH_COMP) {
      debug(printf("Removing nonmatch at 3' end: "));
      debug(Pair_dump_one(lastpair,/*zerobasedp*/true));
      debug(printf("\n"));
      path = Pairpool_pop(path,&lastpair);
      if (path == NULL) {
	return NULL;
      } else {
	lastpair = path->first;
      }
    }
  }

#ifdef PMAP
  while (path != NULL) {
    lastpair = path->first;
    if (lastpair->querypos % 3 == 2) {
      return path;
    } else {
      debug(printf("PMAP popping querypos %d to get to codon boundary\n",lastpair->querypos));
      path = Pairpool_pop(path,&lastpair);
    }
  }
#endif

  return path;
}


static List_T
clean_pairs_end5 (List_T pairs) {
  Pair_T firstpair;

  debug(printf("clean_pairs_end5\n"));
  /* Remove any remaining nonmatches at 5' end, which can happen rarely */
  if (pairs != NULL) {
    firstpair = pairs->first;
    while (firstpair->comp != MATCH_COMP && firstpair->comp != DYNPROG_MATCH_COMP) {
      debug(printf("Removing nonmatch at 5' end: "));
      debug(Pair_dump_one(firstpair,/*zerobasedp*/true));
      debug(printf("\n"));
      pairs = Pairpool_pop(pairs,&firstpair);
      if (pairs == NULL) {
	return NULL;
      } else {
	firstpair = pairs->first;
      }
    }
  }

#ifdef PMAP
  while (pairs != NULL) {
    firstpair = pairs->first;
    if (firstpair->querypos % 3 == 0) {
      return pairs;
    } else {
      debug(printf("PMAP popping querypos %d to get to codon boundary\n",firstpair->querypos));
      pairs = Pairpool_pop(pairs,&firstpair);
    }
  }
#endif

  return pairs;
}


static List_T
chop_ends_by_changepoint (List_T pairs, Pairpool_T pairpool) {
  List_T path;
  Pair_T pair;
  int *matchscores;
  int nmatches, ntotal, nmatches_left, ntotal_left, nmatches_right, ntotal_right;
  int left_edge, right_edge, length, i;
  int side;
  double theta;
  bool chop_left_p = false, chop_right_p = false;

  if (pairs == NULL) {
    return (List_T) NULL;
  } else {
    matchscores = Pair_matchscores_list(&nmatches,&ntotal,&length,pairs);
    debug8(printf("Overall, %d matches/%d total\n",nmatches,ntotal));
  }

  left_edge = Changepoint_left(&nmatches_left,&ntotal_left,matchscores,length);
  right_edge = Changepoint_right(&nmatches_right,&ntotal_right,matchscores,length);

  debug8(printf("At left edge %d (in 0..%d), %d matches/%d total\n",left_edge,List_length(pairs),nmatches_left,ntotal_left));
  debug8(printf("At right edge %d (in 0..%d), %d matches/%d total\n",right_edge,List_length(pairs),nmatches_right,ntotal_right));

  if (right_edge <= left_edge) {
    debug8(printf("Edges cross.  Need to select one.\n"));
    /* Need to select one side to chop. */
    if (ntotal_left == 0 || ntotal - ntotal_left <= 0) {
      side = +1;		/* chop right side */
    } else if (ntotal_right == 0 || ntotal - ntotal_right <= 0) {
      side = -1;		/* chop left side */
    } else {

#if 0
      theta = (double) (nmatches - nmatches_left)/(double) (ntotal - ntotal_left);
      /* Don't have artificially high expectations for theta, e.g., 1.00 */
      theta = theta - THETA_SLACK;
      if (theta < 0.10) {
	/* Protect against negative values */
	theta = 0.10;
      }
      debug8(printf("Testing on left: Pbinom(%d,%d,%f) = %g\n",
		    nmatches_left,ntotal_left,theta,Pbinom(nmatches_left,ntotal_left,theta)));
      pbinom_left = Pbinom(nmatches_left,ntotal_left,theta);

      theta = (double) (nmatches - nmatches_right)/(double) (ntotal - ntotal_right);
      /* Don't have artificially high expectations for theta, e.g., 1.00 */
      theta = theta - THETA_SLACK;
      if (theta < 0.10) {
	/* Protect against negative values */
	theta = 0.10;
      }
      
      debug8(printf("Testing on right: Pbinom(%d,%d,%f) = %g\n",
		    nmatches_right,ntotal_right,theta,Pbinom(nmatches_right,ntotal_right,theta)));
      pbinom_right = Pbinom(nmatches_right,ntotal_right,theta);
      
      if (pbinom_left < pbinom_right) {
	if (pbinom_left < TRIM_END_PVALUE) {
	  side = -1;		/* chop left side */
	} else {
	  side = 0;
	}
      } else if (pbinom_right < pbinom_left) {
	if (pbinom_right < TRIM_END_PVALUE) {
	  side = +1;		/* chop right side */
	} else {
	  side = 0;
	}
      } else {
	side = 0;
      }
#else
      /* Pick shortest side */
      if (ntotal_left < ntotal_right) {
	debug8(printf("left side is shorter\n"));
	side = -1;		/* chop left side */
      } else {
	debug8(printf("right side is shorter\n"));
	side = +1;		/* chop right side */
      }
#endif

    }

    if (side == -1) {
      debug8(printf("Chopping %d on left.\n",left_edge));
      for (i = 0; i < left_edge; i++) {
	pairs = Pairpool_pop(pairs,&pair);
      }
      chop_left_p = true;
    } else if (side == +1) {
      debug8(printf("Chopping %d - %d on right.\n",length,right_edge));
      path = List_reverse(pairs);
      for (i = 0; i < length - right_edge; i++) {
	path = Pairpool_pop(path,&pair);
      }

      pairs = (List_T) NULL;
      while (path) {
	path = Pairpool_pop(path,&pair);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }
      chop_right_p = true;
    }

  } else {
    if (ntotal_left == 0) {
      path = List_reverse(pairs);
    } else if (ntotal - ntotal_left <= 0) {
      path = List_reverse(pairs);
    } else {
      theta = (double) (nmatches - nmatches_left)/(double) (ntotal - ntotal_left);
      /* Don't have artificially high expectations for theta, e.g., 1.00 */
      theta = theta - THETA_SLACK;
      if (theta < 0.10) {
	/* Protect against negative values */
	theta = 0.10;
      }

      debug8(printf("Testing on left: Pbinom(%d,%d,%f) = %g\n",
		    nmatches_left,ntotal_left,theta,Pbinom(nmatches_left,ntotal_left,theta)));
      if (Pbinom(nmatches_left,ntotal_left,theta) > TRIM_END_PVALUE) {
	path = List_reverse(pairs);
      } else {
	debug8(printf("Chopping %d on left.\n",left_edge));
	for (i = 0; i < left_edge; i++) {
	  pairs = Pairpool_pop(pairs,&pair);
	}
	path = (List_T) NULL;
	while (pairs) {
	  pairs = Pairpool_pop(pairs,&pair);
	  path = Pairpool_push_existing(path,pairpool,pair);
	}
	debug8(printf("path is now length %d\n",List_length(path)));
	chop_left_p = true;
      }
    }

    if (ntotal_right == 0) {
      pairs = List_reverse(path);
    } else if (ntotal - ntotal_right <= 0) {
      pairs = List_reverse(path);
    } else {
      theta = (double) (nmatches - nmatches_right)/(double) (ntotal - ntotal_right);
      /* Don't have artificially high expectations for theta, e.g., 1.00 */
      theta = theta - THETA_SLACK;
      if (theta < 0.10) {
	/* Protect against negative values */
	theta = 0.10;
      }

      debug8(printf("Testing on right: Pbinom(%d,%d,%f) = %g\n",
		    nmatches_right,ntotal_right,theta,Pbinom(nmatches_right,ntotal_right,theta)));
      if (Pbinom(nmatches_right,ntotal_right,theta) > TRIM_END_PVALUE) {
	pairs = List_reverse(path);
      } else {
	debug8(printf("Chopping %d - %d on right.\n",length,right_edge));
	for (i = 0; i < length - right_edge; i++) {
	  path = Pairpool_pop(path,&pair);
	}
	pairs = (List_T) NULL;
	while (path) {
	  path = Pairpool_pop(path,&pair);
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
	}
	chop_right_p = true;
      }
    }
  }
  
  FREE(matchscores);

  debug8(printf("Returning alignment of length %d\n",List_length(pairs)));

  return pairs;
}



/* pairs -> pairs */
static List_T
trim_short_end_exons (bool *trim5p, bool *trim3p, List_T pairs, Pairpool_T pairpool, int minendexon) {
  List_T path, exon;
  Pair_T pair;
  int exon_nmatches;

  debug8(printf("Starting trim_short_end_exons\n"));
  debug8(Pair_dump_list(pairs,true));

  /* Handle first exon */
  if (pairs == NULL) {
    *trim5p = *trim3p = false;
    return (List_T) NULL;
  } else {
    pair = pairs->first;
  }

  exon = (List_T) NULL;
  exon_nmatches = 0;
  while (pairs != NULL && !pair->gapp) {
    pairs = Pairpool_pop(pairs,&pair);
    exon = Pairpool_push_existing(exon,pairpool,pair);
    if (pair->gapp == false && pair->comp != MISMATCH_COMP && pair->comp != INDEL_COMP) {
      exon_nmatches++;
    }
  }

  if (exon_nmatches >= minendexon) {
    debug8(printf("Keeping first exon of length %d\n",exon_nmatches));
    path = exon;		/* exon already has the gap */
    *trim5p = false;
  } else if (exon_nmatches == 0) {
    debug8(printf("Trimming first exon of length %d.  firstpair must be a gap.\n",exon_nmatches));
    pairs = Pairpool_pop(pairs,&pair); /* discard gap */
    path = (List_T) NULL;
    *trim5p = false;
  } else {
    debug8(printf("Trimming first exon of length %d\n",exon_nmatches));
    path = (List_T) NULL;
    *trim5p = true;
  }

  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    path = Pairpool_push_existing(path,pairpool,pair);
  }

  /* Handle last exon */
  if (path == NULL) {
    *trim5p = *trim3p = false;
    return (List_T) NULL;
  } else {
    pair = path->first;
  }

  exon = (List_T) NULL;
  exon_nmatches = 0;
  while (path != NULL && !pair->gapp) {
    path = Pairpool_pop(path,&pair);
    exon = Pairpool_push_existing(exon,pairpool,pair);
    if (pair->gapp == false && pair->comp != MISMATCH_COMP && pair->comp != INDEL_COMP) {
      exon_nmatches++;
    }
  }

  if (exon_nmatches >= minendexon) {
    debug8(printf("Keeping last exon of length %d\n",exon_nmatches));
    pairs = exon;		/* exon already has the gap */
    *trim3p = false;
  } else if (exon_nmatches == 0) {
    debug8(printf("Trimming last exon of length %d.  firstpair must be a gap.\n",exon_nmatches));
    path = Pairpool_pop(path,&pair); /* discard gap */
    pairs = (List_T) NULL;
    *trim3p = false;
  } else {
    debug8(printf("Trimming last exon of length %d\n",exon_nmatches));
    pairs = (List_T) NULL;
    *trim3p = true;
  }

  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
  }


  debug8(printf("End of trim_short_end_exons: length = %d\n",List_length(pairs)));
  debug8(Pair_dump_list(pairs,true));
  return pairs;
}



/************************************************************************
 *   Trim
 ************************************************************************/


#if 0
static List_T
trim_end_execute (List_T pairs, int firstgap, int lastgap, Pairpool_T pairpool) {
  List_T path = NULL;
  Pair_T pair;
  int gapi = -1;

  debug8(printf("Starting trim_end_execute with firstgap = %d, lastgap = %d\n",firstgap,lastgap));
  while (pairs != NULL && gapi < firstgap) {
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == true) {
      gapi++;
    }
  }

  debug8(printf("Got to ngaps == %d\n",gapi));
  while (pairs != NULL && gapi < lastgap) {
    pairs = Pairpool_pop(pairs,&pair);
    path = Pairpool_push_existing(path,pairpool,pair);
    if (pair->gapp == true) {
      gapi++;
    }
  }

  debug8(printf("Got to ngaps == %d\n",gapi));
  while (pairs != NULL && gapi == lastgap) {
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == false) {
      path = Pairpool_push_existing(path,pairpool,pair);
    } else {
      gapi++;
    }
  }

  debug8(printf("Done with trim_end_execute\n"));

  return path;
}
#endif


#if 0
static List_T
trim_exons_end (bool *trim5p, bool *trim3p, List_T pairs, Pairpool_T pairpool,
		int ngaps, int ncanonical) {
  List_T path = NULL;
  Pair_T pair;
  double theta;
  int total_nmatches, total_nunknowns, total_nmismatches;
  int nmatches, nmismatches;	/* In each exon */
  int gapi = -1, firstgoodgap = -1, lastgoodgap = -1;
  bool goodp = false, goodexonp = false, goodlastintronp;
  bool goodintronp = false;	/* Forces first intron to be canonical */

  Pair_fracidentity_simple(&total_nmatches,&total_nunknowns,&total_nmismatches,pairs);
  if (total_nmatches == 0 && total_nmismatches == 0) {
    return List_reverse(pairs);
  } else {
    theta = ((double) total_nmatches)/((double) (total_nmatches + total_nmismatches));
    /* Don't have artificially high expectations for theta, e.g., 1.00 */
    theta = theta - THETA_SLACK;
    if (theta < 0.10) {
      /* Protect against negative values */
      theta = 0.10;
    }
  }

  debug(printf("\n** Starting trim_exons_end with ncanonical = %d, total matches %d, total mismatches %d, theta %f\n",
	       ncanonical,total_nmatches,total_nmismatches,theta));

  nmatches = nmismatches = 0;

  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    debug8(Pair_dump_one(pair,true));
    debug8(printf("\n"));
    if (pair->gapp == false) {
      if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
	nmatches++;
      } else {
	nmismatches++;
      }
      path = Pairpool_push_existing(path,pairpool,pair);

    } else {
      goodlastintronp = goodintronp;
      goodintronp = (pair->comp == FWD_CANONICAL_INTRON_COMP || pair->comp == REV_CANONICAL_INTRON_COMP);
      debug8(printf("Intron %d  ... Exon ... Intron %d\n",goodlastintronp,goodintronp));

      if (nmatches + nmismatches == 0) {
	debug8(printf("Setting goodexonp to false because nmatches + nmismatches == 0\n"));
	goodexonp = false;
      } else if (Pbinom(nmatches,nmatches + nmismatches,theta) < TRIM_END_PVALUE) {
	debug8(printf("Setting goodexonp to false because bad pct: pbinom(%d,%d,%f) = %g\n",
		      nmatches,nmatches+nmismatches,theta,Pbinom(nmatches,nmatches+nmismatches,theta)));
	goodexonp = false;
#if 0
      } else if (nmatches < MINENDEXON) {
	debug8(printf("Setting goodexonp to false because nmatches < %d\n",MINENDEXON));
	goodexonp = false;
#endif
#if 0
      } else if (goodlastintronp == false && goodintronp == false) {
	debug8(printf("Setting goodexonp to false because not surrounded by good introns\n"));
	goodexonp = false;
#endif
      } else {
	debug8(printf("Setting goodexonp to true\n"));
	goodexonp = true;
      }
      if (goodexonp == true) {
	debug(printf("Setting lastgoodgap to %d\n",gapi));
	lastgoodgap = gapi;
      }
      if (goodp == false && goodexonp == true) {
	debug(printf("Setting firstgoodgap to %d\n",gapi));
	firstgoodgap = gapi;
	goodp = true;
      }
      debug(printf("gap %d has %d matches and %d mismatches.  goodexonp = %d, %d..%d\n",
		   gapi,nmatches,nmismatches,goodexonp,firstgoodgap,lastgoodgap));

      nmatches = nmismatches = 0;

      path = Pairpool_push_existing(path,pairpool,pair);
      gapi++;
    }
  }

  goodlastintronp = goodintronp;
  goodintronp = false;		/* Forces last intron to be canonical */

  if (nmatches + nmismatches == 0) {
    debug8(printf("Setting goodexonp to false because nmatches + nmismatches == 0\n"));
    goodexonp = false;

  } else if (Pbinom(nmatches,nmatches + nmismatches,theta) < TRIM_END_PVALUE) {
    debug8(printf("Setting goodexonp to false because bad pct: pbinom(%d,%d,%f) = %g\n",
		  nmatches,nmatches+nmismatches,theta,Pbinom(nmatches,nmatches+nmismatches,theta)));
    goodexonp = false;
#if 0
  } else if (nmatches < MINENDEXON) {
    debug8(printf("Setting goodexonp to false because nmatches < %d\n",MINENDEXON));
    goodexonp = false;
#endif
#if 0
  } else if (goodlastintronp == false && goodintronp == false) {
    debug8(printf("Setting goodexonp to false because not surrounded by good introns\n"));
    goodexonp = false;
#endif
  } else {
    debug8(printf("Setting goodexonp to true\n"));
    goodexonp = true;
  }
  if (goodexonp == true) {
    debug(printf("Setting lastgoodgap to %d\n",gapi));
    lastgoodgap = gapi;
  }
  if (goodp == false && goodexonp == true) {
    debug(printf("Setting firstgoodgap to %d\n",gapi));
    firstgoodgap = gapi;
    /* goodp = true; */
  }

  debug(printf("gap %d has %d matches and %d mismatches.  goodexonp = %d, %d..%d\n",
	       gapi,nmatches,nmismatches,goodexonp,firstgoodgap,lastgoodgap));

  *trim5p = (firstgoodgap == -1) ? false : true;
  *trim3p = (lastgoodgap == gapi) ? false : true;
  if (*trim5p == true || *trim3p == true) {
    pairs = List_reverse(path);
    path = trim_end_execute(pairs,firstgoodgap,lastgoodgap,pairpool);
  }

  return path;
}
#endif


static List_T
trim_middle_execute (List_T pairs, Intlist_T bad_exons, Pairpool_T pairpool) {
  List_T path = NULL;
  Pair_T pair;
  int exoni = 0, targeti;
  /* Pair_T leftpair, rightpair */
  /* int queryjump, genomejump; */
  
  while (bad_exons != NULL) {
    bad_exons = Intlist_pop(bad_exons,&targeti);
    debug8(printf("Target is exon %d\n",targeti));

    while (pairs != NULL && exoni < targeti) {
      pairs = Pairpool_pop(pairs,&pair);
      if (pair->gapp == false) {
	path = Pairpool_push_existing(path,pairpool,pair);
      } else if (++exoni < targeti) {
	path = Pairpool_push_existing(path,pairpool,pair);
      }
    }

    debug8(printf("Got to exoni == %d\n",exoni));
    while (pairs != NULL && exoni == targeti) {
      pairs = Pairpool_pop(pairs,&pair);
      if (pair->gapp == true) {
	exoni++;
      }
    }

    /* No need to insert gapholder, since we will call insert_gapholders() afterwards */
    /*
    if (path != NULL && pairs != NULL) {
      leftpair = path->first;
      if (leftpair->gapp == false) {
	rightpair = pairs->first;
	queryjump = rightpair->querypos - leftpair->querypos - 1;
	genomejump = rightpair->genomepos - leftpair->genomepos - 1;
	if (leftpair->cdna == ' ') queryjump++;
	if (leftpair->genome == ' ') genomejump++;
	path = Pairpool_push_gapholder(path,pairpool,queryjump,genomejump);
	pair = path->first;
	pair->comp = DUALBREAK_COMP;
      }
    }
    */

  }

  debug8(printf("No more targets\n"));
  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    path = Pairpool_push_existing(path,pairpool,pair);
  }

  return path;
}


static List_T
trim_exons_middle (bool *trimp, List_T pairs, Pairpool_T pairpool, char *queryseq_ptr, char *genomicseg_ptr,
		   char *genomicuc_ptr, int cdna_direction) {
  Intlist_T bad_exons = NULL;
  List_T path = NULL;
  Pair_T pair;
  int nmatches, nmismatches, total_nmatches, total_nunknowns, total_nmismatches;
  int gapi = -1;
  bool goodexonp = false;
  int npairs, ngaps, ncanonical;
  double pbinom, theta;
  /* bool goodintronp = true; */

  debug8(printf("\n** Starting trim_exons_middle\n"));
  Pair_fracidentity_simple(&total_nmatches,&total_nunknowns,&total_nmismatches,pairs);
  theta = (double) (total_nmatches + 1)/(double) (total_nmatches + total_nmismatches + 1);
  theta = theta - THETA_SLACK;
  if (theta < 0.10) {
    /* Protect against negative values */
    theta = 0.10;
  }
  debug8(printf("%d matches, %d mismatches => theta %f\n",total_nmatches,total_nmismatches,theta));

  nmatches = nmismatches = 0;
  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == false) {
      if (pair->comp == MATCH_COMP || pair->comp == DYNPROG_MATCH_COMP || pair->comp == AMBIGUOUS_COMP) {
	nmatches++;
      } else {
	/* Note: this includes gaps */
	nmismatches++;
      }
      path = Pairpool_push_existing(path,pairpool,pair);

    } else {
      /* goodlastintronp = goodintronp; */
      /* goodintronp = (pair->comp != NONINTRON_COMP && pair->comp != DUALBREAK_COMP); */
      debug8(printf("intron comp is %c\n",pair->comp));

      if ((double) nmatches/(double) (nmatches + nmismatches) >= theta) {
	debug8(printf("Setting goodexonp to true: %d matches, %d nmismatches vs %f\n",
		      nmatches,nmismatches,theta));
	goodexonp = true;
      } else if ((pbinom = Pbinom(nmatches,nmatches+nmismatches,theta)) > 0.0001) {
	debug8(printf("Setting goodexonp to true: %d matches, %d nmismatches vs %f, pbinom %g\n",
		      nmatches,nmismatches,theta,pbinom));
	goodexonp = true;
      } else {
	debug8(printf("Setting goodexonp to false because bad pct: %d matches, %d nmismatches vs %f, pbinom %g\n",
		      nmatches,nmismatches,theta,pbinom));
	goodexonp = false;
      }
      if (goodexonp == false) {
	bad_exons = Intlist_push(bad_exons,gapi+1);
      }
      debug8(printf("exon %d has %d matches and %d mismatches.  goodexonp = %d\n",
		   gapi+1,nmatches,nmismatches,goodexonp));

      nmatches = nmismatches = 0;

      path = Pairpool_push_existing(path,pairpool,pair);
      gapi++;
    }
  }

  /* goodlastintronp = goodintronp; */
  /* goodintronp = true; */

  if ((double) nmatches/(double) (nmatches + nmismatches) >= theta) {
    debug8(printf("Setting goodexonp to true: %d matches, %d nmismatches vs %f\n",
		  nmatches,nmismatches,theta));
    goodexonp = true;
  } else if ((pbinom = Pbinom(nmatches,nmatches+nmismatches,theta)) > 0.0001) {
    debug8(printf("Setting goodexonp to true: %d matches, %d nmismatches vs %f, pbinom %g\n",
		  nmatches,nmismatches,theta,pbinom));
    goodexonp = true;
  } else {
    debug8(printf("Setting goodexonp to false because bad pct: %d matches, %d nmismatches vs %f, pbinom %g\n",
		  nmatches,nmismatches,theta,pbinom));
    goodexonp = false;
  }
  if (goodexonp == false) {
    bad_exons = Intlist_push(bad_exons,gapi+1);
  }
  debug8(printf("exon %d has %d matches and %d mismatches.  goodexonp = %d\n",
		gapi+1,nmatches,nmismatches,goodexonp));

  if (bad_exons == NULL) {
    *trimp = false;
  } else if (gapi < 1) {
    /* Only one or two exons */
    debug8(printf("gapi is %d, so not trimming middle exons\n"));
    *trimp = false;
  } else {
    *trimp = true;
    bad_exons = Intlist_reverse(bad_exons);
    pairs = List_reverse(path);
    path = trim_middle_execute(pairs,bad_exons,pairpool);
    pairs = List_reverse(path);
    path = insert_gapholders(pairs,pairpool);
    pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
    path = List_reverse(pairs);
    /* Intlist_free(&bad_exons); -- No need to free since trim_exons cleans it. */
  }

  return path;
}


/* This procedure fills in introns, so it should be called after all
   dynamic programming procedures */
static List_T 
fill_in_gaps (List_T path, Pairpool_T pairpool, char *queryseq_ptr, 
#ifdef PMAP
	      char *queryaaseq_ptr,
#endif
	      char *genomicseg_ptr, char *genomicuc_ptr, int cdna_direction,
	      int ngap, bool diagnosticp) {
  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;

  if (path == NULL) {
    return (List_T) NULL;
  } else {
    pair = path->first;
  }

  if (pair->gapp == true) {
    /* Gap at beginning of alignment.  Can occur after smoothing. */
    debug7(printf("Gap %p at beginning of alignment\n",pair));
  }

  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
#ifdef PMAP
    if (pair->cdna == BACKTRANSLATE_CHAR) {
      pair->cdna = 'N';
    }
#endif
    if (pair->comp == INDEL_COMP || pair->comp == SHORTGAP_COMP) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->gapp == false) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (path == NULL) {
      /* Gap at end of alignment.  Can occur after smoothing. */
      debug7(printf("Gap at end of alignment\n"));

    } else {
      /* Discard gap; do not push */
      leftpair = path->first;
      rightpair = pairs->first;

      if (pair->comp == DUALBREAK_COMP) {
#ifdef PMAP
	pairs = add_dualbreak(pairs,queryseq_ptr,queryaaseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
			      leftpair,rightpair,pairpool,ngap,diagnosticp);
#else
	pairs = add_dualbreak(pairs,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
			      leftpair,rightpair,pairpool,ngap,diagnosticp);
#endif
      } else {
	debug7(printf("Adding an intron at %d..%d, currently of type %c\n",leftpair->querypos,rightpair->querypos,pair->comp));
	pairs = add_intron(pairs,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,
                           leftpair,rightpair,&pair->comp,ngap,cdna_direction,pairpool);
      }
    }
  }

  debug7(printf("Final length: %d\n",List_length(pairs)));
  return pairs;
}

static List_T 
add_queryseq_offset (List_T path, int queryseq_offset, Pairpool_T pairpool) {
  List_T pairs = NULL;
  Pair_T pair;

  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == false) {
      pair->querypos += queryseq_offset;
    }
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
  }

  return pairs;
}


static void
add_skiplength (List_T pairs, int skiplength) {
  List_T p;
  Pair_T pair;

  for (p = pairs; p != NULL; p = p->rest) {
    pair = (Pair_T) p->first;
    if (pair->querypos >= HALFLEN) {
      pair->querypos += skiplength;
    }
  }
  return;
}


static struct Pair_T *
prepare_for_printing (int *npairs, List_T *pairs, int cdna_direction, Pairpool_T pairpool,
		      char *queryseq_ptr,
#ifdef PMAP
		      char *queryaaseq_ptr,
#endif
		      char *genomicseg_ptr, char *genomicuc_ptr, int ngap,
		      int subseq_offset, int skiplength, bool diagnosticp) {
  struct Pair_T *pairarray;
  List_T path, p;
  Pair_T oldpair, newpair;

  path = List_reverse(*pairs);
#ifdef PMAP
  *pairs = fill_in_gaps(path,pairpool,queryseq_ptr,queryaaseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,ngap,
			diagnosticp);
#else
  *pairs = fill_in_gaps(path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,ngap,diagnosticp);
#endif

  if (subseq_offset != 0) {
    path = List_reverse(*pairs);
    *pairs = add_queryseq_offset(path,subseq_offset,pairpool);
  }
  if (skiplength != 0) {
    add_skiplength(*pairs,skiplength);
  }
  *npairs = List_length(*pairs);

  /* Used to be Pair_block_copy */
  newpair = pairarray = (struct Pair_T *) MALLOC(*npairs*sizeof(struct Pair_T));
  for (p = *pairs; p != NULL; p = p->rest) {
    oldpair = (Pair_T) p->first;
    memcpy(newpair++,oldpair,sizeof(struct Pair_T));
  }

  /* No need to free newpairs, since they belong to pairpool */

  return pairarray;
}

	    
void
Stage3_recompute_goodness (List_T stage3list) {
  T this;
  List_T p;

#ifdef DEPEND_ON_QUALITY
  for (p = stage3list; p != NULL && high_quality_p == false; p = List_next(p)) {
    this = (T) List_head(p);
    if (this->trimmed_coverage > 0.80 && this->defect_rate < 0.10) {
      high_quality_p = true;
    }
  }

  /* Subtracts points for non-canonical introns */
  if (high_quality_p == true) {
    for (p = stage3list; p != NULL; p = List_next(p)) {
      this = (T) List_head(p);
      this->goodness = this->matches + MISMATCH*this->mismatches
	+ QOPEN*this->qopens + QINDEL*this->qindels + TOPEN*this->topens + TINDEL*this->tindels;
      /* - *CANONICAL_POINTS*this->noncanonical; */
    }
  } else {
    for (p = stage3list; p != NULL; p = List_next(p)) {
      this = (T) List_head(p);
      this->goodness = this->matches;
    }
  }

#else
  for (p = stage3list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    this->goodness = this->matches + MISMATCH*this->mismatches
      + QOPEN*this->qopens + QINDEL*this->qindels + TOPEN*this->topens + TINDEL*this->tindels;
  }
#endif
}


static T
Stage3_new (List_T pairs_fwd, List_T pairs_rev, int totalintronscore_fwd, int totalintronscore_rev,
	    Genomicpos_T stage1_genomicstart, Genomicpos_T stage1_genomiclength, double defect_rate) {
  T new;
  int ncanonical, nsemicanonical, nnoncanonical;
  int nexons;
  int winning_cdna_direction;

  if (pairs_fwd == NULL && pairs_rev == NULL) {
    return NULL;
  } else {
    
    if (pairs_rev == NULL) {
      winning_cdna_direction = +1;
    } else if (pairs_fwd == NULL) {
      winning_cdna_direction = -1;
    } else if (totalintronscore_fwd > totalintronscore_rev) {
      winning_cdna_direction = +1;
    } else if (totalintronscore_fwd < totalintronscore_rev) {
      winning_cdna_direction = -1;
    } else {
      winning_cdna_direction = +1; /* Choose fwd over rev in case of tie */
    }

    new = (T) MALLOC(sizeof(*new));

    if (winning_cdna_direction == +1) {
      new->pairs = pairs_fwd;
    } else {
      new->pairs = pairs_rev;
    }

    Pair_fracidentity(&new->matches,&new->unknowns,&new->mismatches,
		      &new->qopens,&new->qindels,&new->topens,&new->tindels,
		      &ncanonical,&nsemicanonical,&nnoncanonical,new->pairs,
		      /*cdna_direction*/winning_cdna_direction);
    new->noncanonical = nsemicanonical + nnoncanonical;
    nexons = Pair_nexons_approx(new->pairs);

    if (new->matches < MIN_MATCHES) {
      /* Too few matches */
      new->pairarray = (struct Pair_T *) NULL;
      Stage3_free(&new);
      return NULL;
    }

    if (ncanonical > 0 || nsemicanonical > 0) {
      new->cdna_direction = winning_cdna_direction;
    } else {
      new->cdna_direction = 0;
    }

#if 0
    if (nexons > 2) {
      /* Favor spliced transcripts, but only if we're sure they're
         spliced (i.e., 3 or more exons).  A random intron shouldn't
         get credit. */
      new->goodness += nexons;
    }
#endif
    
    new->translation_start = 0;
    new->translation_end = 0;
    new->translation_length = 0;

    new->stage1_genomicstart = stage1_genomicstart;
    new->stage1_genomiclength = stage1_genomiclength;

    new->stage2_diag_runtime = 0.0;
    new->stage2_align_runtime = 0.0;
    new->stage2_indexsize = 0;

    return new;
  }
}

void
Stage3_free (T *old) {

  if (*old) {
    /* Don't free strain.  Belongs to altstrain_iit. */
    FREE((*old)->pairarray);
    FREE(*old);
  }
  return;
}

/* Needed for mutation analysis */
void
Stage3_genomicbounds (Genomicpos_T *genomicstart, Genomicpos_T *genomiclength, T this) {
  *genomicstart = this->chroffset + this->chrpos;
  *genomiclength = this->genomiclength;
  return;
}


bool
Stage3_test_bounds (T this, int minpos, int maxpos) {
  if (Pairpool_count_bounded(this->pairs,minpos,maxpos) > 25) {
    return true;
  } else {
    return false;
  }
}

T
Stage3_apply_bounds (T this, int minpos, int maxpos, bool revertp) {
  int newnpairs;
  List_T p;
  Pair_T newpair, oldpair;

  this->pairs = Pairpool_transfer_bounded(NULL,this->pairs,minpos,maxpos);
  this->pairs = List_reverse(this->pairs);
  if ((newnpairs = List_length(this->pairs)) == 0) {
    if (revertp == false) {
      return NULL;
    }
  } else {
    FREE(this->pairarray);
    this->npairs = newnpairs;
      
    newpair = this->pairarray = (struct Pair_T *) MALLOC(newnpairs*sizeof(struct Pair_T));
    for (p = this->pairs; p != NULL; p = p->rest) {
      oldpair = (Pair_T) p->first;
      memcpy(newpair++,oldpair,sizeof(struct Pair_T));
    }
  }

  return this;
}


#ifdef PMAP
Stage3_T
Stage3_translate_cdna (T this, Sequence_T queryaaseq, bool strictp) {
  Translation_via_cdna(&this->translation_start,&this->translation_end,&this->translation_length,
		       &this->relaastart,&this->relaaend,
		       this->pairarray,this->npairs,Sequence_fullpointer(queryaaseq),strictp);
  return this;
}

Stage3_T
Stage3_backtranslate_cdna (T this, bool diagnosticp) {
  Backtranslation_cdna(this->pairarray,this->npairs,this->translation_start,this->translation_end,
		       diagnosticp);
  return this;
}

#else
static Stage3_T
Stage3_truncate_fulllength (Stage3_T this, bool translatep, bool strictp) {

  if (translatep == true) {
    if (this->cdna_direction < 0) {
      Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			      &this->relaastart,&this->relaaend,
			      this->pairarray,this->npairs,/*backwardsp*/true,/*revcompp*/true,/*fulllengthp*/true,
			      strictp);
    } else {
      Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			      &this->relaastart,&this->relaaend,
			      this->pairarray,this->npairs,/*backwardsp*/false,/*revcompp*/false,/*fulllengthp*/true,
			      strictp);
    }
  }

  this = Stage3_apply_bounds(this,Stage3_translation_start(this),Stage3_translation_end(this),/*revertp*/true);
  return this;
}

Stage3_T
Stage3_translate_genomic (T this, bool fulllengthp, bool truncatep, bool strictp) {
  Stage3_T new;

  if (this->cdna_direction < 0) {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairarray,this->npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp,
			    strictp);
  } else {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairarray,this->npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp,
			    strictp);
  }
  if (truncatep == false) {
    return this;
  } else {
    new = Stage3_truncate_fulllength(this,/*translatep*/false,strictp);
    if (new->cdna_direction < 0) {
      Translation_via_genomic(&new->translation_start,&new->translation_end,&new->translation_length,
			      &new->relaastart,&new->relaaend,
			      new->pairarray,new->npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp,
			      strictp);
    } else {
      Translation_via_genomic(&new->translation_start,&new->translation_end,&new->translation_length,
			      &new->relaastart,&new->relaaend,
			      new->pairarray,new->npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp,
			      strictp);
    }
    return new;
  }
}
#endif

void
Stage3_translate_cdna_via_reference (T this, T reference, bool literalrefp) {
  bool fixshiftp = !literalrefp;

  if (this->watsonp == reference->watsonp) {
    if (reference->cdna_direction < 0) {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairarray,this->npairs,this->watsonp,/*backwardsp*/true,/*revcompp*/true,
				reference->pairarray,reference->npairs,reference->watsonp,reference->genomiclength,
				fixshiftp);
    } else {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairarray,this->npairs,this->watsonp,/*backwardsp*/false,/*revcompp*/false,
				reference->pairarray,reference->npairs,reference->watsonp,reference->genomiclength,
				fixshiftp);
    }
  } else {
    if (reference->cdna_direction < 0) {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairarray,this->npairs,this->watsonp,/*backwardsp*/false,/*revcompp*/false,
				reference->pairarray,reference->npairs,reference->watsonp,reference->genomiclength,
				fixshiftp);
    } else {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairarray,this->npairs,this->watsonp,/*backwardsp*/true,/*revcompp*/true,
				reference->pairarray,reference->npairs,reference->watsonp,reference->genomiclength,
				fixshiftp);
    }
  }

  return;
}

void
Stage3_fix_cdna_direction (T this, T reference) {
  if (this->cdna_direction == 0) {
    if (reference->cdna_direction > 0) {
      if (this->watsonp == reference->watsonp) {
	this->cdna_direction = +1;
      } else {
	this->cdna_direction = -1;
      }
    } else if (reference->cdna_direction < 0) {
      if (this->watsonp == reference->watsonp) {
	this->cdna_direction = -1;
      } else {
	this->cdna_direction = +1;
      }
    }
  }
  return;
}


void
Stage3_print_pathsummary (T this, int pathnum, IIT_T chromosome_iit, IIT_T contig_iit, 
			  IIT_T altstrain_iit, Sequence_T queryseq, bool fulllengthp,
			  bool truncatep, bool strictp, char *dbversion, int maxmutations, 
			  bool zerobasedp, bool diagnosticp, bool maponlyp) {
  Pair_T start, end;
  bool referencealignp;

  if (maponlyp == true) {
    this->translation_start = 0;
    this->translation_end = 0;
    this->translation_length = 0;
  } else {
#ifdef PMAP
    Stage3_translate_cdna(this,queryseq,strictp);
    Stage3_backtranslate_cdna(this,diagnosticp);
#else
    Stage3_translate_genomic(this,fulllengthp,truncatep,strictp);
#endif
  }

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);
  referencealignp = this->straintype == 0 ? true : false;
  Pair_print_pathsummary(pathnum,start,end,this->chrnum,this->chrpos,this->chroffset,
			 chromosome_iit,referencealignp,altstrain_iit,this->strain,contig_iit,
			 dbversion,Sequence_fulllength_given(queryseq),Sequence_skiplength(queryseq),
			 Sequence_trim_start(queryseq),Sequence_trim_end(queryseq),this->genomiclength,
			 Pair_nexons(this->pairarray,this->npairs),this->matches,this->unknowns,this->mismatches,
			 this->qopens,this->qindels,this->topens,this->tindels,this->goodness,
			 this->watsonp,this->cdna_direction,this->defect_rate,
			 this->translation_start,this->translation_end,this->translation_length,
			 0,0,zerobasedp,maponlyp,diagnosticp,this->stage1_genomicstart,
			 this->stage1_genomiclength,this->stage2_diag_runtime,this->stage2_align_runtime,
			 this->stage2_indexsize,this->stage3_runtime,this->defect_rate);
  if (maponlyp == false) {
    Translation_compare(this->pairarray,this->npairs,NULL,0,this->cdna_direction,
			this->relaastart,this->relaaend,maxmutations);
  }
  printf("\n");

  return;
}

void
Stage3_print_pslformat_nt (T this, int pathnum, IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryaaseq) {
  Pair_T start, end;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_pslformat_nt(this->pairarray,this->npairs,start,end,queryaaseq,this->chrnum,this->chrpos,
			  chromosome_iit,usersegment,this->genomiclength,Pair_nexons(this->pairarray,this->npairs),
			  this->matches,this->unknowns,this->mismatches, 
			  this->qopens,this->qindels,this->topens,this->tindels,
			  this->watsonp,this->cdna_direction);
  return;
}

#ifdef PMAP
void
Stage3_print_pslformat_pro (T this, int pathnum, IIT_T chromosome_iit, Sequence_T usersegment, Sequence_T queryaaseq, bool strictp) {
  Pair_T start, end;

  Stage3_translate_cdna(this,queryaaseq,strictp);
  Stage3_backtranslate_cdna(this,/*diagnosticp*/false);

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_pslformat_pro(this->pairarray,this->npairs,start,end,queryaaseq,this->chrnum,this->chrpos,
			   chromosome_iit,usersegment,this->genomiclength,Pair_nexons(this->pairarray,this->npairs),
			   this->qopens,this->qindels,this->topens,this->tindels,
			   this->watsonp,this->cdna_direction);
  return;
}
#endif


void
Stage3_print_gff3 (T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq, char *dbversion,
		   bool diagnosticp, bool fulllengthp, bool truncatep, bool strictp,
		   bool gff_gene_format_p, char *user_genomicseg) {
  Pair_T start, end;

  if (gff_gene_format_p == true) {
#ifdef PMAP
    Stage3_translate_cdna(this,queryseq,strictp);
    Stage3_backtranslate_cdna(this,diagnosticp);
#else
    Stage3_translate_genomic(this,fulllengthp,truncatep,strictp);
#endif
  }

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_gff3(this->pairarray,this->npairs,pathnum,
		  dbversion,Sequence_accession(queryseq),start,end,
		  this->chrnum,this->chrpos,chromosome_iit,this->genomiclength,
		  this->translation_start,this->translation_end,
		  Sequence_fulllength_given(queryseq),Sequence_skiplength(queryseq),
		  this->matches,this->unknowns,this->mismatches,
		  this->qopens,this->qindels,this->topens,this->tindels,
		  this->watsonp,this->cdna_direction,gff_gene_format_p,user_genomicseg);
  return;
}


void
Stage3_print_iit_map (T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq) {
  Pair_T start, end;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_iit_map(queryseq,pathnum,Sequence_accession(queryseq),start,end,
		     this->chrnum,this->chrpos,chromosome_iit,this->genomiclength,
		     this->watsonp);
  return;
}

void
Stage3_print_iit_exon_map (T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq) {
  Pair_T start, end;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_iit_exon_map(this->pairarray,this->npairs,queryseq,pathnum,Sequence_accession(queryseq),
			  start,end,this->chrnum,this->chrpos,chromosome_iit,this->genomiclength,
			  this->watsonp);
  return;
}

void
Stage3_print_splicesites (T this, IIT_T chromosome_iit, Sequence_T queryseq) {
  Pair_print_splicesites(this->pairarray,this->npairs,Sequence_accession(queryseq),
			 Pair_nexons(this->pairarray,this->npairs),this->chrnum,this->chrpos,
			 chromosome_iit,this->genomiclength,this->watsonp);
  return;
}



void
Stage3_print_mutations (T this, T reference, IIT_T chromosome_iit, Sequence_T queryseq,
			char *dbversion, bool showalignp, bool zerobasedp, 
			bool continuousp, bool diagnosticp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength,
			int maxmutations) {
  Pair_T start, end;
  bool referencealignp;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  /*  Pair_dump_array(this->pairarray,this->npairs,false); */

  referencealignp = this->straintype == 0 ? true : false;
  Pair_print_pathsummary(/*pathnum*/1,start,end,reference->chrnum,reference->chrpos,reference->chroffset,
			 chromosome_iit,referencealignp,/*altstrain_iit*/NULL,this->strain,/*contig_iit*/NULL,
			 dbversion,Sequence_fulllength_given(queryseq),Sequence_skiplength(queryseq),
			 Sequence_trim_start(queryseq),Sequence_trim_end(queryseq),reference->genomiclength,
			 Pair_nexons(this->pairarray,this->npairs),this->matches,this->unknowns,this->mismatches,
			 this->qopens,this->qindels,this->topens,this->tindels,this->goodness,
			 this->watsonp,this->cdna_direction,this->defect_rate,
			 0,0,0,this->relaastart,this->relaaend,zerobasedp,/*maponlyp*/false,
			 diagnosticp,this->stage1_genomicstart,this->stage1_genomiclength,
			 this->stage2_diag_runtime,this->stage2_align_runtime,this->stage2_indexsize,
			 this->stage3_runtime,this->defect_rate);
  Translation_compare(this->pairarray,this->npairs,reference->pairarray,reference->npairs,
		      this->cdna_direction,this->relaastart,this->relaaend,maxmutations);
  printf("\n");

  if (showalignp == true) {
    Pair_print_alignment(this->pairarray,this->npairs,reference->chrnum,reference->chrpos,reference->chroffset,
			 chromosome_iit,reference->genomiclength,this->watsonp,this->cdna_direction,/*universalp*/false,zerobasedp,
			 diagnosticp,/*genomicprop*/false,invertmode,nointronlenp,wraplength);
  }
  debug1(Pair_dump_array(this->pairarray,this->npairs,/*zerobasedp*/true));
  debug1(Pair_check_array(this->pairarray,this->npairs));

  return;
}


static int
determine_typeint (T this, IIT_T map_iit, bool map_iit_universal_p, int map_iit_forward_type,
		   int map_iit_reverse_type, IIT_T chromosome_iit, bool map_bothstrands_p) {
  int typeint;
  char *typestring, *stranded;
  bool allocp;

  if (this->watsonp) {
    if (this->cdna_direction >= 0) {
      if (map_iit_universal_p == true) {
	typeint = map_iit_forward_type;
      } else if (map_bothstrands_p == true) {
	typestring = IIT_label(chromosome_iit,this->chrnum,&allocp);
	typeint = IIT_typeint(map_iit,typestring);
	if (allocp == true) {
	  FREE(typestring);
	}
      } else {
	typestring = IIT_label(chromosome_iit,this->chrnum,&allocp);
	stranded = (char *) CALLOC(strlen(typestring)+strlen("+")+1,sizeof(char));
	sprintf(stranded,"%s+",typestring);
	if ((typeint = IIT_typeint(map_iit,stranded)) < 0) {
	  typeint = IIT_typeint(map_iit,typestring);
	}
	if (allocp == true) {
	  FREE(typestring);
	}
	FREE(stranded);
      }
    } else {
      if (map_iit_universal_p == true) {
	typeint = map_iit_reverse_type;
      } else if (map_bothstrands_p == true) {
	typestring = IIT_label(chromosome_iit,this->chrnum,&allocp);
	typeint = IIT_typeint(map_iit,typestring);
	if (allocp == true) {
	  FREE(typestring);
	}
      } else {
	typestring = IIT_label(chromosome_iit,this->chrnum,&allocp);
	stranded = (char *) CALLOC(strlen(typestring)+strlen("-")+1,sizeof(char));
	sprintf(stranded,"%s-",typestring);
	if ((typeint = IIT_typeint(map_iit,stranded)) < 0) {
	  typeint = IIT_typeint(map_iit,typestring);
	}
	if (allocp == true) {
	  FREE(typestring);
	}
	FREE(stranded);
      }
    }

  } else {
    if (this->cdna_direction >= 0) {
      if (map_iit_universal_p == true) {
	typeint = map_iit_reverse_type;
      } else if (map_bothstrands_p == true) {
	typestring = IIT_label(chromosome_iit,this->chrnum,&allocp);
	typeint = IIT_typeint(map_iit,typestring);
	if (allocp == true) {
	  FREE(typestring);
	}
      } else {
	typestring = IIT_label(chromosome_iit,this->chrnum,&allocp);
	stranded = (char *) CALLOC(strlen(typestring)+strlen("-")+1,sizeof(char));
	sprintf(stranded,"%s-",typestring);
	if ((typeint = IIT_typeint(map_iit,stranded)) < 0) {
	  typeint = IIT_typeint(map_iit,typestring);
	}
	if (allocp == true) {
	  FREE(typestring);
	}
	FREE(stranded);
      }
    } else {
      if (map_iit_universal_p == true) {
	typeint = map_iit_forward_type;
      } else if (map_bothstrands_p == true) {
	typestring = IIT_label(chromosome_iit,this->chrnum,&allocp);
	typeint = IIT_typeint(map_iit,typestring);
	if (allocp == true) {
	  FREE(typestring);
	}
      } else {
	typestring = IIT_label(chromosome_iit,this->chrnum,&allocp);
	stranded = (char *) CALLOC(strlen(typestring)+strlen("+")+1,sizeof(char));
	sprintf(stranded,"%s+",typestring);
	if ((typeint = IIT_typeint(map_iit,stranded)) < 0) {
	  typeint = IIT_typeint(map_iit,typestring);
	}
	if (allocp == true) {
	  FREE(typestring);
	}
	FREE(stranded);
      }
    }
  }

  return typeint;
}



static void
print_map (T this, IIT_T map_iit, int *map_divint_crosstable,
	   bool map_iit_universal_p, int map_iit_forward_type,
	   int map_iit_reverse_type, IIT_T chromosome_iit, int pathnum, bool map_bothstrands_p,
	   int nflanking) {
  Genomicpos_T chrlow, chrhigh;
  Pair_T start, end;
  int chrpos1, chrpos2;
  int *iit_matches = NULL, nmatches, *leftflanks, nleftflanks, *rightflanks, nrightflanks;
  int divno, sign;

  if ((divno = map_divint_crosstable[this->chrnum]) <= 0) {
    printf("  *Map hits for path %d (0):\n\n",pathnum);
    return;
  }

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  if (this->watsonp) {
    chrlow = chrpos1 = this->chrpos + Pair_genomepos(start);
    chrhigh = chrpos2 = this->chrpos + Pair_genomepos(end);
    sign = +1;

  } else {
    chrhigh = chrpos1 = this->chrpos + (this->genomiclength - 1) - Pair_genomepos(start);
    chrlow = chrpos2 = this->chrpos + (this->genomiclength - 1) - Pair_genomepos(end);
    sign = -1;
  }

  if (map_bothstrands_p == true) {
    iit_matches = IIT_get_with_divno(&nmatches,map_iit,divno,chrlow,chrhigh,/*sortp*/false);
    if (nflanking > 0) {
      IIT_get_flanking_with_divno(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,map_iit,
				  divno,chrlow,chrhigh,nflanking,/*sign*/0);
    }
    if (nflanking > 0) {
      printf("  Map hits for path %d (%d|%d|%d):\n",pathnum,nleftflanks,nmatches,nrightflanks);
    } else {
      printf("  Map hits for path %d (%d):\n",pathnum,nmatches);
    }
    if (nflanking > 0) {
      IIT_print_header(map_iit,leftflanks,nleftflanks,/*bothstrandsp*/true,chromosome_iit,
		       /*reversep*/true,/*relativep*/false,/*left*/0U);
      printf("    ====================\n");
    }
    IIT_print_header(map_iit,iit_matches,nmatches,/*bothstrandsp*/true,chromosome_iit,
		     /*reversep*/false,/*relativep*/false,/*left*/0U);
    if (nflanking > 0) {
      printf("    ====================\n");
      IIT_print_header(map_iit,rightflanks,nrightflanks,/*bothstrandsp*/true,chromosome_iit,
		       /*reversep*/false,/*relativep*/false,/*left*/0U);
    }

  } else {
    iit_matches = IIT_get_signed_with_divno(&nmatches,map_iit,divno,chrlow,chrhigh,/*sortp*/true,sign);
    if (nflanking > 0) {
      IIT_get_flanking_with_divno(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,map_iit,
				  divno,chrlow,chrhigh,nflanking,sign);
    }
    if (nflanking > 0) {
      printf("  Map hits for path %d (%d|%d|%d):\n",pathnum,nleftflanks,nmatches,nrightflanks);
    } else {
      printf("  Map hits for path %d (%d):\n",pathnum,nmatches);
    }
    if (nflanking > 0) {
      IIT_print_header(map_iit,leftflanks,nleftflanks,/*bothstrandsp*/false,chromosome_iit,
		       /*reversep*/true,/*relativep*/false,/*left*/0U);
      printf("    ====================\n");
    }
    IIT_print_header(map_iit,iit_matches,nmatches,/*bothstrandsp*/false,chromosome_iit,
		     /*reversep*/false,/*relativep*/false,/*left*/0U);
    if (nflanking > 0) {
      printf("    ====================\n");
      IIT_print_header(map_iit,rightflanks,nrightflanks,/*bothstrandsp*/false,chromosome_iit,
		       /*reversep*/false,/*relativep*/false,/*left*/0U);
    }
  }
  printf("\n");

  if (nflanking > 0) {
    FREE(rightflanks);
    FREE(leftflanks);
  }

  FREE(iit_matches);
  return;
}


/* Doesn't handle nflanking */
static void
print_exon_map (T this, IIT_T map_iit, bool map_iit_universal_p, int map_iit_forward_type,
		int map_iit_reverse_type, IIT_T chromosome_iit, int pathnum, bool map_bothstrands_p) {
  Uintlist_T exonbounds;
  Genomicpos_T position1, position2;
  int *iit_matches = NULL, nmatches;
  int typeint, exonno = 0;
  char *divstring;

  if ((typeint = determine_typeint(this,map_iit,map_iit_universal_p,map_iit_forward_type,
				   map_iit_reverse_type,chromosome_iit,map_bothstrands_p)) < 0) {
    printf("  Map hits for path %d (0):\n\n",pathnum);
    return;
  }

  if (map_iit_universal_p == true) {
    exonbounds = Pair_exonbounds(this->pairarray,this->npairs,this->chrpos,this->chroffset,this->genomiclength,
				 this->watsonp);
    divstring = (char *) NULL;
  } else {
    exonbounds = Pair_exonbounds(this->pairarray,this->npairs,this->chrpos,/*chroffset*/0U,this->genomiclength,
				 this->watsonp);
    divstring = IIT_typestring(chromosome_iit,typeint);
  }

  while (exonbounds != NULL) {
    exonbounds = Uintlist_pop(exonbounds,&position1);
    exonbounds = Uintlist_pop(exonbounds,&position2);

    if (map_bothstrands_p == true) {
      if (position1 < position2) {
	iit_matches = IIT_get(&nmatches,map_iit,divstring,position1,position2,/*sortp*/true);
      } else {
	iit_matches = IIT_get(&nmatches,map_iit,divstring,position2,position1,/*sortp*/true);
      }
      printf("  Map hits for path %d, exon %d (%d):\n",pathnum,++exonno,nmatches);
      IIT_print(map_iit,iit_matches,nmatches,/*bothstrandsp*/true,chromosome_iit,/*levels*/NULL,
		/*reversep*/false,/*relativep*/false,/*left*/0U);
    } else {
      if (position1 < position2) {
	iit_matches = IIT_get_typed(&nmatches,map_iit,divstring,position1,position2,typeint,/*sortp*/true);
      } else {
	iit_matches = IIT_get_typed(&nmatches,map_iit,divstring,position2,position1,typeint,/*sortp*/true);
      }
      printf("  Map hits for path %d, exon %d (%d):\n",pathnum,++exonno,nmatches);
      IIT_print(map_iit,iit_matches,nmatches,/*bothstrandsp*/false,chromosome_iit,/*levels*/NULL,
		/*reversep*/false,/*relativep*/false,/*left*/0U);
    }
    printf("\n");
    FREE(iit_matches);
  }

  return;
}

void
Stage3_print_map (T this, IIT_T map_iit, int *map_divint_crosstable,
		  bool map_iit_universal_p, int map_iit_forward_type,
		  int map_iit_reverse_type, IIT_T chromosome_iit,
		  int pathnum, bool map_exons_p, bool map_bothstrands_p, int nflanking) {
  if (map_exons_p == true) {
    print_exon_map(this,map_iit,map_iit_universal_p,map_iit_forward_type,map_iit_reverse_type,
		   chromosome_iit,pathnum,map_bothstrands_p);
  } else {
    print_map(this,map_iit,map_divint_crosstable,
	      map_iit_universal_p,map_iit_forward_type,map_iit_reverse_type,
	      chromosome_iit,pathnum,map_bothstrands_p,nflanking);
  }
  return;
}



/* queryaaseq is used only by PMAP */
void
Stage3_print_alignment (T this, Sequence_T queryaaseq, Genome_T genome,
			IIT_T chromosome_iit, bool alignsummaryonlyp, bool universalp, bool zerobasedp,
			bool continuousp, bool continuous_by_exon_p, bool diagnosticp, bool strictp, bool genomefirstp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength, bool maponlyp) {
  if (continuous_by_exon_p == true) {
    Pair_print_exonsummary(this->pairarray,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			   genome,chromosome_iit,this->genomiclength,this->watsonp,universalp,zerobasedp,
			   genomefirstp,invertmode);
    Pair_print_continuous_byexon(this->pairarray,this->npairs,this->watsonp,diagnosticp,invertmode);

  } else if (continuousp == true) {
    if (maponlyp == false) {
#ifdef PMAP
      Stage3_translate_cdna(this,queryaaseq,strictp);
      Stage3_backtranslate_cdna(this,diagnosticp);
#endif
    }
    Pair_print_continuous(this->pairarray,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			  this->genomiclength,this->watsonp,this->cdna_direction,universalp,zerobasedp,
			  diagnosticp,genomefirstp,invertmode,nointronlenp);
  } else {
    /* Assumes Stage3_print_pathsummary already called */
    Pair_print_exonsummary(this->pairarray,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			   genome,chromosome_iit,this->genomiclength,this->watsonp,universalp,zerobasedp,
			   genomefirstp,invertmode);
    if (alignsummaryonlyp == false) {
      Pair_print_alignment(this->pairarray,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			   chromosome_iit,this->genomiclength,this->watsonp,
			   this->cdna_direction,universalp,zerobasedp,
			   diagnosticp,/*genomicprop*/true,invertmode,nointronlenp,wraplength);
    }
  }
  debug1(Pair_dump_array(this->pairarray,this->npairs,/*zerobasedp*/true));
  debug1(Pair_check_array(this->pairarray,this->npairs));
  return;
}


/* queryaaseq is used only by PMAP */
void
Stage3_print_coordinates (T this, Sequence_T queryaaseq, IIT_T chromosome_iit, bool zerobasedp,
			  int invertmode, bool fulllengthp, bool truncatep, bool strictp, bool maponlyp) {
  if (maponlyp == false) {
#ifdef PMAP
    Stage3_translate_cdna(this,queryaaseq,strictp);
    Stage3_backtranslate_cdna(this,/*diagnosticp*/false);
#else
    Stage3_translate_genomic(this,fulllengthp,truncatep,strictp);
#endif
  }
  Pair_print_coordinates(this->pairarray,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			 chromosome_iit,this->genomiclength,this->watsonp,
			 zerobasedp,invertmode);
  return;
}


void
Stage3_print_cdna_exons (T this, int wraplength, int ngap) {
  Pair_print_exons(this->pairarray,this->npairs,wraplength,ngap,/*cdnap*/true);
  return;
}

void
Stage3_print_genomic_exons (T this, int wraplength, int ngap) {
  Pair_print_exons(this->pairarray,this->npairs,wraplength,ngap,/*cdnap*/false);
  return;
}


void
Stage3_print_cdna (T this, Sequence_T queryaaseq, bool fulllengthp, bool truncatep, bool strictp, int wraplength) {
#ifdef PMAP
  Stage3_translate_cdna(this,queryaaseq,strictp);
  Stage3_backtranslate_cdna(this,/*diagnosticp*/false);
  Pair_print_nucleotide_cdna(this->pairarray,this->npairs,wraplength);
#else
  Stage3_translate_genomic(this,fulllengthp,truncatep,strictp);
  if (this->cdna_direction >= 0) {
    Pair_print_protein_cdna(this->pairarray,this->npairs,wraplength,/*forwardp*/true);
  } else {
    Pair_print_protein_cdna(&(this->pairarray[this->npairs-1]),this->npairs,wraplength,/*forwardp*/false);
  }
#endif
  return;
}

void
Stage3_print_protein_genomic (T this, Sequence_T queryaaseq, bool fulllengthp, bool truncatep, bool strictp,
			      int wraplength) {
#ifdef PMAP
  Stage3_translate_cdna(this,queryaaseq,strictp);
  Stage3_backtranslate_cdna(this,/*diagnosticp*/false);
#else
  Stage3_translate_genomic(this,fulllengthp,truncatep,strictp);
#endif
  if (this->cdna_direction >= 0) {
    Pair_print_protein_genomic(this->pairarray,this->npairs,wraplength,/*forwardp*/true);
  } else {
    Pair_print_protein_genomic(&(this->pairarray[this->npairs-1]),this->npairs,wraplength,/*forwardp*/false);
  }
  return;
}


void
Stage3_print_compressed (T this, Sequence_T queryseq, IIT_T chromosome_iit,
			 char *dbversion, int pathnum, int npaths,
			 bool checksump, int chimerapos, int chimeraequivpos,
			 double donor_prob, double acceptor_prob, 
			 int chimera_cdna_direction, bool zerobasedp, bool truncatep, bool strictp,
			 int worker_id) {
  Pair_T start, end;

#ifdef PMAP
  Stage3_translate_cdna(this,queryseq,strictp);
  Stage3_backtranslate_cdna(this,/*diagnosticp*/false);
#else
  if (truncatep == true) {
    Stage3_truncate_fulllength(this,/*translatep*/true,strictp);
  }
#endif

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);
  Pair_print_compressed(pathnum,npaths,start,end,queryseq,dbversion,Pair_nexons(this->pairarray,this->npairs),
			Stage3_fracidentity(this),this->pairarray,this->npairs,
			this->chrnum,this->chrpos,this->chroffset,chromosome_iit,
			Sequence_fulllength_given(queryseq),Sequence_skiplength(queryseq),
			Sequence_trim_start(queryseq),Sequence_trim_end(queryseq),
			this->genomiclength,checksump,
			chimerapos,chimeraequivpos,donor_prob,acceptor_prob,
			chimera_cdna_direction,this->strain,this->watsonp,
			this->cdna_direction,zerobasedp,worker_id);
  return;
}



#if 0
static int
compute_introntype (char left1, char left2, char right2, char right1) {
  int leftdi, rightdi;

  if (left1 == 'G' && left2 == 'T') {
    leftdi = LEFT_GT;
  } else if (left1 == 'G' && left2 == 'C') {
    leftdi = LEFT_GC;
  } else if (left1 == 'A' && left2 == 'T') {
    leftdi = LEFT_AT;
#ifndef PMAP
  } else if (left1 == 'C' && left2 == 'T') {
    leftdi = LEFT_CT;
#endif
  } else {
    leftdi = 0x00;
  }

  if (right2 == 'A' && right1 == 'G') {
    rightdi = RIGHT_AG;
  } else if (right2 == 'A' && right1 == 'C') {
    rightdi = RIGHT_AC;
#ifndef PMAP
  } else if (right2 == 'G' && right1 == 'C') {
    rightdi = RIGHT_GC;
  } else if (right2 == 'A' && right1 == 'T') {
    rightdi = RIGHT_AT;
#endif
  } else {
    rightdi = 0x00;
  }

  return leftdi & rightdi;
}
#endif

static char uppercaseCode[128] = UPPERCASE_U2T;

static List_T
peel_leftward (bool *mismatchp, List_T *peeled_path, List_T path, int *querydp5, int *genomedp5, 
	       Pairpool_T pairpool, int maxpeelback, bool throughmismatchp,
	       List_T *endgappairs, char *gapcomp, int *querydp5_medialgap, int *genomedp5_medialgap) {
  List_T peeled = NULL, rest = NULL;
  Pair_T pair, nextpair, rightpair;
  int npeelback = 0, nconsecutive = 0, init_dynprogindex = DYNPROGINDEX_MINOR;
  bool stopp;
#if 0
  int nincursion = 0;
#endif

  *mismatchp = false;
  debug(printf("Peeling leftward:"));
  if (path == NULL) {
    debug(printf(" path is empty\n"));
  } else {
    pair = path->first;
    if (pair->gapp == true) {
      /* Throw away known gap */
      debug(printf(" Known_gap"));
      path = Pairpool_pop(path,&pair);
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
    }
    rest = path->rest;
    
    stopp = false;
    while (rest != NULL && stopp == false) {
      nextpair = rest->first;
      if (nextpair->gapp == true) {
	stopp = true;
      }

      path = Pairpool_pop(path,&pair);
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
      debug(printf(" Peel [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));

      if (uppercaseCode[(int) pair->cdna] != uppercaseCode[(int) pair->genome]) {
	*mismatchp = true;
      }

      if (++npeelback >= maxpeelback) {
	stopp = true;
      }

      if (init_dynprogindex > 0 && pair->dynprogindex <= 0) {
	init_dynprogindex = pair->dynprogindex;
      }

      rest = path->rest;
    }

    if (throughmismatchp && rest != NULL && nextpair->gapp == false) {
      /* Continue to peelback through little skips and mismatches */
      debug(printf("\n||"));

      stopp = false;
      while (rest != NULL && stopp == false) {
	nextpair = rest->first;
	if (nextpair->gapp == true) {
	  stopp = true;
	}

	path = Pairpool_pop(path,&pair);
	peeled = Pairpool_push_existing(peeled,pairpool,pair);
	debug(printf(" Extrapeel [");
	      Pair_dump_one(pair,/*zerobasedp*/true);
	      printf("]"));

	if (uppercaseCode[(int) pair->cdna] != uppercaseCode[(int) pair->genome]) {
	  *mismatchp = true;
	}

	if (pair->comp == INDEL_COMP || pair->comp == MISMATCH_COMP) {
	  nconsecutive = 0;
	} else if (++nconsecutive >= SUFFCONSECUTIVE) {
	  stopp = true;
	}
      
#if 0
	if (pair->dynprogindex != init_dynprogindex) {
	  if (++nincursion >= MAXINCURSION) {
	    stopp = true;
	  }
	}
#endif

	rest = path->rest;
      }
    }
  }

#ifdef PMAP
  /* Reverse process to codon boundary.  Cases:

  X X X | X -
  0 1 2   3 4

  X X X | - X
  0 1 2   3 3

  X X - X | X
  0 1 2 2   3

  Rule: nextpair->querypos % 3 == 0 */

  debug(printf("\n<<"));
  if (peeled != NULL) {
    rest = peeled->rest;
    stopp = false;
    while (rest != NULL && stopp == false) {
      peeled = Pairpool_pop(peeled,&pair);
      path = Pairpool_push_existing(path,pairpool,pair);
      debug(printf(" Mod3putback [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));
      nextpair = rest->first;
      if (nextpair->querypos % 3 == 0) {
	stopp = true;
      }
      rest = peeled->rest;
    }
  }
#endif

  if (peeled == NULL) {
    /* Do not alter querydp5 or genomedp5 */
  } else {
    rightpair = peeled->first;
    if (rightpair->gapp == true) {
      debug(printf("Ran into gap; undoing peel\n"));
      /* was abort() */
      path = Pairpool_transfer(path,peeled);
      *peeled_path = (List_T) NULL;
      return path;
    } else {
      *querydp5 = rightpair->querypos;
      *genomedp5 = rightpair->genomepos;
    }
  }

  if (endgappairs != NULL) {
    if (path == NULL || (pair = path->first) == NULL || pair->gapp == false) {
      *endgappairs = NULL;
      *querydp5_medialgap = *querydp5;
      *genomedp5_medialgap = *genomedp5;
    } else {
      path = Pairpool_pop(path,&pair);
      *endgappairs = Pairpool_push_existing(NULL,pairpool,pair);
      debug(printf(" Peeling gap [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));
      *gapcomp = pair->comp;

      path = Pairpool_pop(path,&pair);
      *endgappairs = Pairpool_push_existing(*endgappairs,pairpool,pair);
      debug(printf(" Peeling after gap [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));

      rightpair = (*endgappairs)->first;
      if (rightpair->gapp == true) {
	debug(printf("Ran into gap; undoing peel\n"));
	/* was abort() */
	path = Pairpool_transfer(path,peeled);
	*peeled_path = (List_T) NULL;
	return path;
      } else {
	*querydp5_medialgap = rightpair->querypos;
	*genomedp5_medialgap = rightpair->genomepos;
      }
    }
  }

  debug(
	if (path == NULL) {
	  printf(" => Top of path is NULL.");
	} else {
	  pair = path->first;
	  printf(" => Top of path is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp5 = %d, genomedp5 = %d\n",*querydp5,*genomedp5);
	);

  *peeled_path = peeled;
  return path;
}


static List_T
peel_rightward (bool *mismatchp, List_T *peeled_pairs, List_T pairs, int *querydp3, int *genomedp3, 
		Pairpool_T pairpool, int maxpeelback, bool throughmismatchp,
		List_T *endgappairs, char *gapcomp, int *querydp3_medialgap, int *genomedp3_medialgap) {
  List_T peeled = NULL, rest = NULL;
  Pair_T pair, nextpair, leftpair;
  int npeelback = 0, nconsecutive = 0, init_dynprogindex = DYNPROGINDEX_MINOR;
  bool stopp;
#if 0
  int incursion = 0;
#endif

  *mismatchp = false;
  debug(printf("Peeling rightward:"));
  if (pairs == NULL) {
    debug(printf(" pairs is empty\n"));
  } else {
    pair = pairs->first;
    if (pair->gapp == true) {
      /* Throw away known gap */
      debug(printf(" Known_gap"));
      pairs = Pairpool_pop(pairs,&pair);
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
    }
    rest = pairs->rest;

    stopp = false;
    while (rest != NULL && stopp == false) {
      nextpair = rest->first;
      if (nextpair->gapp == true) {
	stopp = true;
      }
    
      pairs = Pairpool_pop(pairs,&pair);
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
      debug(printf(" Peel [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));
      
      if (uppercaseCode[(int) pair->cdna] != uppercaseCode[(int) pair->genome]) {
	*mismatchp = true;
      }

      if (++npeelback >= maxpeelback) {
	stopp = true;
      }

      if (init_dynprogindex > 0 && pair->dynprogindex <= 0) {
	init_dynprogindex = pair->dynprogindex;
      }

      rest = pairs->rest;
    }

    if (throughmismatchp && rest != NULL && nextpair->gapp == false) {
      /* Continue to peelback through little skips and mismatches */
      debug(printf("\n||"));

      stopp = false;
      while (rest != NULL && stopp == false) {
	nextpair = rest->first;
	if (nextpair->gapp == true) {
	  stopp = true;
	}
	
	pairs = Pairpool_pop(pairs,&pair);
	peeled = Pairpool_push_existing(peeled,pairpool,pair);
	debug(printf(" Extrapeel [");
	      Pair_dump_one(pair,/*zerobasedp*/true);
	      printf("]"));
	
	if (uppercaseCode[(int) pair->cdna] != uppercaseCode[(int) pair->genome]) {
	  *mismatchp = true;
	}

	if (pair->comp == INDEL_COMP || pair->comp == MISMATCH_COMP) {
	  nconsecutive = 0;
	} else if (++nconsecutive >= SUFFCONSECUTIVE) {
	  stopp = true;
	}
	
#if 0
	if (pair->dynprogindex != init_dynprogindex) {
	  if (++nincursion >= MAXINCURSION) {
	    stopp = true;
	  }
	}
#endif

	rest = pairs->rest;
      }
    }
  }

#ifdef PMAP
  /* Reverse process to codon boundary.  Cases:

  - X | X X X
  5 5   6 7 8

  X - | X X X
  5 6   6 7 8

  X | X - X X
  5   6 7 7 8

  Rule: pair->querypos % 3 == 0 */

  debug(printf("\n<<"));
  stopp = false;
  while (peeled != NULL && stopp == false) {
    peeled = Pairpool_pop(peeled,&pair);
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
    debug(printf(" Mod3putback [");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	  printf("]"));
    if (pair->querypos % 3 == 0) {
      stopp = true;
    }
  }
#endif

  if (peeled == NULL) {
    /* Do not alter querydp3 or genomedp3 */
  } else {
    leftpair = peeled->first;
    if (leftpair->gapp == true) {
      debug(printf("Ran into gap; undoing peel\n"));
      /* was abort() */
      pairs = Pairpool_transfer(pairs,peeled);
      *peeled_pairs = (List_T) NULL;
      return pairs;

    } else {
      if (leftpair->cdna == ' ') {
	*querydp3 = leftpair->querypos - 1;
      } else {
	*querydp3 = leftpair->querypos;
      }
      if (leftpair->genome == ' ') {
	*genomedp3 = leftpair->genomepos - 1;
      } else {
	*genomedp3 = leftpair->genomepos;
      }
    }
  }

  if (endgappairs != NULL) {
    if (pairs == NULL || (pair = pairs->first) == NULL || pair->gapp == false) {
      *endgappairs = NULL;
      *querydp3_medialgap = *querydp3;
      *genomedp3_medialgap = *genomedp3;
    } else {
      pairs = Pairpool_pop(pairs,&pair);
      *endgappairs = Pairpool_push_existing(NULL,pairpool,pair);
      debug(printf(" Peeling gap [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));
      *gapcomp = pair->comp;

      pairs = Pairpool_pop(pairs,&pair);
      *endgappairs = Pairpool_push_existing(*endgappairs,pairpool,pair);
      debug(printf(" Peeling after gap [");
	    Pair_dump_one(pair,/*zerobasedp*/true);
	    printf("]"));

      leftpair = (*endgappairs)->first;
      if (leftpair->gapp == true) {
	debug(printf("Ran into gap; undoing peel\n"));
	/* was abort() */
	pairs = Pairpool_transfer(pairs,peeled);
	*peeled_pairs = (List_T) NULL;
	return pairs;
      } else {
	if (leftpair->cdna == ' ') {
	  *querydp3_medialgap = leftpair->querypos - 1;
	} else {
	  *querydp3_medialgap = leftpair->querypos;
	}
	if (leftpair->genome == ' ') {
	  *genomedp3_medialgap = leftpair->genomepos - 1;
	} else {
	  *genomedp3_medialgap = leftpair->genomepos;
	}
      }
    }
  }

  debug(
	if (pairs == NULL) {
	  printf(" => Top of pairs is NULL.");
	} else {
	  pair = pairs->first;
	  printf(" => Top of pairs is ");
	  Pair_dump_one(pair,/*zerobasedp*/true);
	}
	printf("\n => querydp3 = %d, genomedp3 = %d\n",*querydp3,*genomedp3);
	);

  *peeled_pairs = peeled;
  return pairs;
}


/************************************************************************
 *  Traversal functions
 ************************************************************************/

static List_T
traverse_single_gap (bool *filledp, int *dynprogindex, List_T pairs, List_T *path, 
		     Pair_T leftpair, Pair_T rightpair,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprog,
		     int maxpeelback, int extraband_single, double defect_rate, bool forcep,
		     bool widebandp) {
  List_T gappairs, peeled_pairs, peeled_path;
  int queryjump, genomejump;
  int querydp5, genomedp5, querydp3, genomedp3;
  int nmatches, nmismatches, nopens, nindels;
  int finalscore;
  bool mismatchp = false;
  char gapcomp;
  /* int origqueryjump, origgenomejump; */

  debug(printf("\nTRAVERSE_SINGLE_GAP\n"));
  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  if (leftpair->cdna == ' ') querydp5--;
  if (leftpair->genome == ' ') genomedp5--;
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;

  /* origqueryjump = querydp3 - querydp5 + 1; */
  /* origgenomejump = genomedp3 - genomedp5 + 1; */

  /* Used to peelback only half as much as for a paired gap, to save
     on dynamic programming, but not any more. */
  pairs = peel_rightward(&mismatchp,&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
			 /*throughmismatchp*/false,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);
  *path = peel_leftward(&mismatchp,&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
			/*throughmismatchp*/false,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);

  queryjump = querydp3 - querydp5 + 1;
  genomejump = genomedp3 - genomedp5 + 1;
  
  if (queryjump <= 0 || genomejump <= 0) {
    debug(printf("Unable to perform dynamic programming\n"));
    *filledp = false;

    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);

    return pairs;

  } else {
    gappairs = Dynprog_single_gap(&(*dynprogindex),&finalscore,
				  &nmatches,&nmismatches,&nopens,&nindels,dynprog,
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				  queryjump,genomejump,querydp5,genomedp5,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  cdna_direction,pairpool,extraband_single,defect_rate,
				  /*widebandp*/true);
    debug(Pair_dump_list(gappairs,true));
  }
  debug(printf("  Score: %d\n",finalscore));

  if (!forcep && nmismatches + nopens > nmatches) {
    /* Put back peeled pairs */
    debug(printf("Bad alignment, so undoing this solution\n"));
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);

#if 0
    /* Calling procedure should put the gap onto pairs */
    leftpair = (*path)->first;
    if (leftpair->gapp == true) {
      debug(printf("Moving the gapholder over\n"));
      *path = Pairpool_pop(*path,&leftpair);
      pairs = Pairpool_push_existing(pairs,pairpool,leftpair);
    }
#endif

    *filledp = false;
  } else {
    pairs = Pairpool_transfer(pairs,gappairs);
    *filledp = true;
  }

  return pairs;
}

static List_T
traverse_cdna_gap (bool *filledp, bool *incompletep, int *dynprogindex_minor, int *dynprogindex_major,
		   List_T pairs, List_T *path, Pair_T leftpair, Pair_T rightpair,
#ifdef PMAP
		   char *queryaaseq_ptr,
#endif
		   char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		   int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, 
		   int maxpeelback, int extramaterial_paired, int extraband_paired, int extraband_single,
		   double defect_rate, bool forcep) {
  List_T gappairs, peeled_pairs = NULL, peeled_path = NULL;
  int queryjump, genomejump;
  int querydp5, genomedp5, querydp3, genomedp3;
  int finalscore;
  int nmatches, nmismatches, nopens, nindels;
  bool mismatchp = false, throughmismatchp;
  char gapcomp;

  debug(printf("\nTRAVERSE_CDNA_GAP\n"));
  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  if (leftpair->cdna == ' ') querydp5--;
  if (leftpair->genome == ' ') genomedp5--;
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;

  if (leftpair->dynprogindex < 0 && leftpair->dynprogindex == rightpair->dynprogindex) {
    debug(printf("Re-peeling prior solution\n"));
    throughmismatchp = false;
  } else {
    debug(printf("No prior solution\n"));
    throughmismatchp = true;
  }

  pairs = peel_rightward(&mismatchp,&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
			 throughmismatchp,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);
  *path = peel_leftward(&mismatchp,&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
			throughmismatchp,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);

#if 0
  if (peeled_pairs == NULL || peeled_path == NULL) {
    debug(printf("Skipping this because unable to peel\n"));
    *filledp = false;
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    return pairs;
  }
#endif

  queryjump = querydp3 - querydp5 + 1;
  genomejump = genomedp3 - genomedp5 + 1;

  if (queryjump <= genomejump + MININTRONLEN) {
    debug(printf("Really a single gap, not a cDNA gap\n"));
    gappairs = Dynprog_single_gap(&(*dynprogindex_minor),&finalscore,
				  &nmatches,&nmismatches,&nopens,&nindels,dynprogM,
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				  queryjump,genomejump,querydp5,genomedp5,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  cdna_direction,pairpool,extraband_single,defect_rate,
				  /*widebandp*/true);
    debug(Pair_dump_list(gappairs,true));
    debug(printf("  Score: %d\n",finalscore));
    pairs = Pairpool_transfer(pairs,gappairs);
    *filledp = true;

  } else {
    /* Set queryjump approximately equal to genomejump to have square
       dynamic programming matrices */
    queryjump = genomejump + extramaterial_paired;
    gappairs = Dynprog_cdna_gap(&(*dynprogindex_major),&finalscore,&(*incompletep),dynprogL,dynprogR,
				&(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				&(queryseq_ptr[querydp3]),&(queryuc_ptr[querydp3]),
				&(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				/*length1L*/queryjump,/*length1R*/queryjump,/*length2*/genomejump,
				/*offset1L*/querydp5,/*revoffset1R*/querydp3,/*offset2*/genomedp5,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				cdna_direction,pairpool,extraband_paired,defect_rate,forcep);
    debug(Pair_dump_list(gappairs,true));
    *filledp = true;
    if (gappairs == NULL) {
      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);
      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);
    } else {
      pairs = Pairpool_transfer(pairs,gappairs);
    }
  }

  return pairs;
}


/* genome_gap is usually an intron */
/* Do not set shiftp to false */
static List_T
traverse_genome_gap (bool *filledp, bool *shiftp, int *dynprogindex_minor, int *dynprogindex_major,
		     int *nintrons, int *nnonintrons, int *intronlen, int *nonintronlen, 
		     List_T pairs, List_T *path, Pair_T leftpair, Pair_T rightpair,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, Pairpool_T pairpool,
		     Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, int maxpeelback,
		     int extramaterial_paired, int extraband_paired, int extraband_single,
		     double defect_rate, bool forcep, bool finalp) {
  List_T gappairs, peeled_pairs = NULL, peeled_path = NULL, micropairs;
  int queryjump, genomejump;
  int querydp5, genomedp5, querydp3, genomedp3;
  int new_leftgenomepos, new_rightgenomepos;
  int finalscore, nmatches, nmismatches, nopens, nindels, exonhead, introntype, microintrontype;
  int acceptable_nmismatches;
  bool mismatch_rightward_p = false, mismatch_leftward_p = false, throughmismatchp;

#ifdef SHORTCUT
  char left1, left2, right2, right1;
#endif
  char gapcomp;

  debug(printf("\nTRAVERSE_GENOME_GAP\n"));

  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  if (leftpair->cdna == ' ') querydp5--;
  if (leftpair->genome == ' ') genomedp5--;
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;

  if (leftpair->dynprogindex < 0 && leftpair->dynprogindex == rightpair->dynprogindex) {
    debug(printf("Re-peeling prior solution\n"));
    throughmismatchp = false;
  } else {
    debug(printf("No prior solution\n"));
    throughmismatchp = true;
  }

#ifdef SHORTCUT
  queryjump = querydp3 - querydp5 + 1;
  genomejump = genomedp3 - genomedp5 + 1;

  if (querydp5 != querydp3 + 1) {
    pairs = peel_rightward(&mismatch_rightward_p,&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
			   throughmismatchp,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);
    *path = peel_leftward(&mismatch_leftward_p,&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
			  throughmismatchp,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);

  } else {
    left1 = genomicuc_ptr[genomedp5];
    left2 = genomicuc_ptr[genomedp5+1];
    right2 = genomicuc_ptr[genomedp3-1];
    right1 = genomicuc_ptr[genomedp3];

    introntype = Intron_type(left1,left2,right2,right1,cdna_direction);
    debug(printf("Introntype is %s\n",Intron_type_string(introntype)));

    pairs = peel_rightward(&mismatch_rightward_p,&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
			   throughmismatchp,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);
    *path = peel_leftward(&mismatch_leftward_p,&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
			  throughmismatchp,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);


    if (mismatch_rightward_p == false && mismatch_leftward_p == false) {
      debug(printf("No mismatches seen\n"));
      if ((cdna_direction > 0 && introntype == GTAG_FWD) ||
	  (cdna_direction < 0 && introntype == GTAG_REV)) {
	debug(printf("Skipping because intron is already canonical\n"));
	*filledp = false;		/* Calling procedure will replace the gap */
	pairs = Pairpool_transfer(pairs,peeled_pairs);
	*path = Pairpool_transfer(*path,peeled_path);
	return pairs;
      }
    }
  }

#else
  pairs = peel_rightward(&mismatch_rightward_p,&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
			 throughmismatchp,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);
  *path = peel_leftward(&mismatch_leftward_p,&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
			throughmismatchp,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);
#endif


  queryjump = querydp3 - querydp5 + 1;
  genomejump = genomedp3 - genomedp5 + 1;

  /* genomedp5 + genomejump - 1 >= genomedp3 - genomejump + 1) ?  but doesn't work on AA669154, chr1*/
  if (genomejump <= queryjump + MININTRONLEN) {
    debug(printf("Really a single gap, not an intron\n"));
    gappairs = Dynprog_single_gap(&(*dynprogindex_minor),&finalscore,
				  &nmatches,&nmismatches,&nopens,&nindels,dynprogM,
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				  queryjump,genomejump,querydp5,genomedp5,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  cdna_direction,pairpool,extraband_single,defect_rate,
				  /*widebandp*/true);
    debug(Pair_dump_list(gappairs,true));
    debug(printf("  Score: %d\n",finalscore));

    pairs = Pairpool_transfer(pairs,gappairs);
    *filledp = true;

  } else {
    /* Set genomejump approximately equal to queryjump to have square
       dynamic programming matrices */
    genomejump = queryjump + extramaterial_paired;
    gappairs = Dynprog_genome_gap(&(*dynprogindex_major),&finalscore,&new_leftgenomepos,&new_rightgenomepos,
				  &nmatches,&nmismatches,&nopens,&nindels,
				  &exonhead,&introntype,dynprogL,dynprogR,
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				  &(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
				  queryjump,genomejump,genomejump,
				  querydp5,genomedp5,genomedp3,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  cdna_direction,pairpool,extraband_paired,
				  defect_rate,maxpeelback,/*halfp*/false,finalp);



    if (new_leftgenomepos != leftpair->genomepos || new_rightgenomepos != rightpair->genomepos) {
      *shiftp = true;
      debug(printf("Shift in intron location from %d..%d to %d..%d\n",
		   leftpair->genomepos,rightpair->genomepos,new_leftgenomepos,new_rightgenomepos));
    } else {
      /* *shiftp = false; */
      debug(printf("No shift in intron location\n"));
    }
    debug(Pair_dump_list(gappairs,true));
    debug(printf("  Score: %d\n",finalscore));

    if (defect_rate < DEFECT_HIGHQ) {
      acceptable_nmismatches = 0;
    } else if (defect_rate < DEFECT_MEDQ) {
      acceptable_nmismatches = 2;
    } else {
      acceptable_nmismatches = 3;
    }

    debug(printf("nmismatches = %d, nopens = %d, nindels = %d.  acceptable nmismatches = %d\n",
		 nmismatches,nopens,nindels,acceptable_nmismatches));
    
    if (!forcep && finalscore < 0) {
      *filledp = false;
      /* Put back peeled pairs */
      debug(printf("Not forced and finalscore is negative\n"));
      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);
      introntype = NONINTRON;
    } else if (false && defect_rate > DEFECT_MEDQ) {
      /* Don't look for microexons in low-quality sequences */
      debug(printf("Don't look for microexon in low-quality sequence\n"));
      *filledp = true;
      pairs = Pairpool_transfer(pairs,gappairs);
    } else if (introntype != NONINTRON && nmismatches <= acceptable_nmismatches && nopens <= 1 && nindels <= 3) {
      debug(printf("introntype != NONINTRON and nmismatches, nopens,  nindels low\n"));
      *filledp = true;
      pairs = Pairpool_transfer(pairs,gappairs);
    } else {
#ifdef PMAP
      *filledp = true;
      pairs = Pairpool_transfer(pairs,gappairs);
#else      
      *filledp = true;
      debug(printf("Calling microexon because introntype == %d or nmismatches %d > acceptable %d or nopens %d > 1 or nindels %d > 3\n",
		   introntype,nmismatches,acceptable_nmismatches,nopens,nindels));

      micropairs = Dynprog_microexon_int(&(*dynprogindex_major),&microintrontype,
					 /*sequence1*/&(queryseq_ptr[querydp5]),
					 /*sequenceuc1*/&(queryuc_ptr[querydp5]),
					 /*sequence2L*/&(genomicseg_ptr[genomedp5]),
					 /*sequenceuc2L*/&(genomicuc_ptr[genomedp5]),
					 /*revsequence2R*/&(genomicseg_ptr[genomedp3]),
					 /*revsequenceuc2R*/&(genomicuc_ptr[genomedp3]),
					 /*length1*/queryjump,/*length2L*/genomejump,/*length2R*/genomejump,
					 /*offset1*/querydp5,/*offset2L*/genomedp5,/*revoffset2R*/genomedp3,
					 cdna_direction,
#ifdef PMAP
					 queryaaseq_ptr,
#endif
					 queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
					 pairpool,defect_rate);

      if (micropairs == NULL) {
	debug(printf("No microexon found\n"));
	pairs = Pairpool_transfer(pairs,gappairs);
	/* *shiftp = false; */
      } else {
	debug(printf("Microexon found\n"));
	pairs = Pairpool_transfer(pairs,micropairs);
	introntype = microintrontype;
	*shiftp = true;
      }
#endif
    }

    if (introntype == NONINTRON) {
      *nnonintrons += 1;
      *nonintronlen += new_rightgenomepos - new_leftgenomepos - 1;
    } else {
      *nintrons += 1;
      *intronlen += new_rightgenomepos - new_leftgenomepos - 1;
    }
  }

  return pairs;
}


/* Only set singlep to true, not to false */
static List_T
traverse_dual_genome_gap (bool *singlep, int *dynprogindex, List_T pairs, List_T *path, 
			  Pair_T leftpair, Pair_T rightpair, int midquerypos, int midgenomepos, 
#ifdef PMAP
			  char *queryaaseq_ptr,
#endif
			  char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
			  int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, 
			  int maxpeelback, int nullgap, int extramaterial_paired, int extraband_paired,
			  double defect_rate) {
  List_T single_gappairs, dual_gappairs_1 = NULL, dual_gappairs_2 = NULL, peeled_pairs, peeled_path;
  int queryjump, genomejump;
  int querydp5, genomedp5, querydp3, genomedp3;
  int new_leftgenomepos, new_rightgenomepos;
  int querydp5_dual, querydp3_dual, genomedp5_dual, genomedp3_dual;
  int single_nmatches = 0, dual_nmatches_1 = 0, dual_nmatches_2 = 0;
  int single_score, dual_score_1, dual_score_2, single_goodness, dual_goodness, 
    nmismatches, nopens, nindels, exonhead, right_exonhead, left_exonhead;
  int middle_exonlength, interexon_region;
  int single_introntype, dual_introntype_1, dual_introntype_2, 
    canonical_introntype, semicanonical_introntype_1, semicanonical_introntype_2;
  double middle_exonprob;
  bool mismatchp = false, single_canonical_p, dual_canonical_p;
  char gapcomp;

  debug(printf("\nTRAVERSE_DUAL_GENOME_GAP\n"));
  if (cdna_direction > 0) {
    canonical_introntype = GTAG_FWD;
    semicanonical_introntype_1 = ATAC_FWD;
    semicanonical_introntype_2 = GCAG_FWD;
  } else {
    canonical_introntype = GTAG_REV;
    semicanonical_introntype_1 = ATAC_REV;
    semicanonical_introntype_2 = GCAG_REV;
  }

  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  if (leftpair->cdna == ' ') querydp5--;
  if (leftpair->genome == ' ') genomedp5--;
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;

  pairs = peel_rightward(&mismatchp,&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
			 /*throughmismatchp*/true,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);
  *path = peel_leftward(&mismatchp,&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
			/*throughmismatchp*/true,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);

  queryjump = querydp3 - querydp5 + 1;
  genomejump = queryjump + extramaterial_paired;

  if (queryjump > nullgap) {
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);

    /* No need to compute queryjump and genomejump */
    /*
    leftpair = (*path)->first;
    rightpair = pairs->first;
    queryjump = rightpair->querypos - leftpair->querypos - 1;
    genomejump = rightpair->genomepos - leftpair->genomepos - 1;
    if (leftpair->cdna == ' ') queryjump++;
    if (leftpair->genome == ' ') genomejump++;
    */

    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);

    return pairs;
  }

  if (genomedp5 + genomejump - 1 >= genomedp3 - genomejump + 1) {
    debug(printf("Bounds don't make sense for dual intron gap: %d + %d - 1 >= %d - %d + 1\n\n",
		 genomedp5,genomejump,genomedp3,genomejump));
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);

    /* No need to compute queryjump and genomejump */
    /*
    leftpair = (*path)->first;
    rightpair = pairs->first;
    queryjump = rightpair->querypos - leftpair->querypos - 1;
    genomejump = rightpair->genomepos - leftpair->genomepos - 1;
    if (leftpair->cdna == ' ') queryjump++;
    if (leftpair->genome == ' ') genomejump++;
    */

    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);

    return pairs;
  }
  
  single_gappairs = Dynprog_genome_gap(&(*dynprogindex),&single_score,&new_leftgenomepos,&new_rightgenomepos,
				       &single_nmatches,&nmismatches,&nopens,&nindels,
				       &exonhead,&single_introntype,dynprogL,dynprogR,
				       &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				       &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				       &(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
				       queryjump,genomejump,genomejump,
				       querydp5,genomedp5,genomedp3,
#ifdef PMAP
				       queryaaseq_ptr,
#endif
				       cdna_direction,pairpool,extraband_paired,
				       defect_rate,maxpeelback,/*halfp*/false,/*finalp*/false);

  debug(Pair_check_list(single_gappairs));

  /* Okay to have one indel, because may need to shift an island */
  if (nopens <= 1) {
    single_goodness = (single_nmatches + nindels) + MISMATCH*nmismatches;
  } else {
    single_goodness = single_nmatches + MISMATCH*nmismatches + QOPEN*(nopens-1) + QINDEL*nindels;
  }

  if (single_introntype == canonical_introntype) {
    single_canonical_p = true;
  } else if (single_introntype == semicanonical_introntype_1 ||
	     single_introntype == semicanonical_introntype_2) {
    single_canonical_p = false;
  } else {
    single_canonical_p = false;
  }


  /* Right of short exon */
  querydp5_dual = midquerypos;
  genomedp5_dual = midgenomepos;
  querydp3_dual = querydp3;	/* From peel_rightward */
  genomedp3_dual = genomedp3;	/* From peel_rightward */

  queryjump = querydp3_dual - querydp5_dual + 1;
  genomejump = queryjump + extramaterial_paired;

  if (genomedp5_dual + genomejump - 1 >= genomedp3_dual) {
    /* Bounds don't make sense */
    debug(printf("Bounds don't make sense on right of dual intron gap: %d + %d - 1 >= %d\n\n",
		 genomedp5_dual,genomejump,genomedp3_dual));
    dual_gappairs_2 = NULL;

  } else {
    dual_gappairs_2 = Dynprog_genome_gap(&(*dynprogindex),&dual_score_2,&new_leftgenomepos,&new_rightgenomepos,
					 &dual_nmatches_2,&nmismatches,&nopens,&nindels,
					 &right_exonhead,&dual_introntype_2,dynprogL,dynprogR,
					 &(queryseq_ptr[querydp5_dual]),&(queryuc_ptr[querydp5_dual]),
					 &(genomicseg_ptr[genomedp5_dual]),&(genomicuc_ptr[genomedp5_dual]),
					 &(genomicseg_ptr[genomedp3_dual]),&(genomicuc_ptr[genomedp3_dual]),
					 queryjump,genomejump,genomejump,
					 querydp5_dual,genomedp5_dual,genomedp3_dual,
#ifdef PMAP
					 queryaaseq_ptr,
#endif
					 cdna_direction,pairpool,extraband_paired,
					 defect_rate,maxpeelback,/*halfp*/true,/*finalp*/false);

    dual_goodness = dual_nmatches_2 + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

    /* Left of short exon */
    querydp5_dual = querydp5;	/* From peel_leftward */
    genomedp5_dual = genomedp5;	/* From peel_leftward */
    querydp3_dual = midquerypos-1;
    genomedp3_dual = midgenomepos-1;

    queryjump = querydp3_dual - querydp5_dual + 1;
    genomejump = queryjump + extramaterial_paired;

    if (genomedp5_dual + genomejump - 1 >= genomedp3_dual) {
      /* Bounds don't make sense */
      debug(printf("Bounds don't make sense on left of dual intron gap: %d + %d - 1 >= %d\n\n",
		   genomedp5_dual,genomejump,genomedp3_dual));
      dual_gappairs_1 = NULL;

    } else {
      dual_gappairs_1 = Dynprog_genome_gap(&(*dynprogindex),&dual_score_1,&new_leftgenomepos,&new_rightgenomepos,
					   &dual_nmatches_1,&nmismatches,&nopens,&nindels,
					   &left_exonhead,&dual_introntype_1,dynprogL,dynprogR,
					   &(queryseq_ptr[querydp5_dual]),&(queryuc_ptr[querydp5_dual]),
					   &(genomicseg_ptr[genomedp5_dual]),&(genomicuc_ptr[genomedp5_dual]),
					   &(genomicseg_ptr[genomedp3_dual]),&(genomicuc_ptr[genomedp3_dual]),
					   queryjump,genomejump,genomejump,
					   querydp5_dual,genomedp5_dual,genomedp3_dual,
#ifdef PMAP
					   queryaaseq_ptr,
#endif
					   cdna_direction,pairpool,extraband_paired,
					   defect_rate,maxpeelback,/*halfp*/true,/*finalp*/false);
      
      dual_goodness += dual_nmatches_1 + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

      if (dual_introntype_1 == canonical_introntype && dual_introntype_2 == canonical_introntype) {
	dual_canonical_p = true;
      } else {
	dual_canonical_p = false;
      }
    }
  }

  if (dual_gappairs_2 == NULL || dual_gappairs_1 == NULL) {
    debug(printf("Single score wins\n"));
    debug(printf("Loser: dual_gappairs\n"));
    debug(Pair_dump_list(dual_gappairs_2,true));
    debug(Pair_dump_list(dual_gappairs_1,true));
    debug(printf("Winner: Transferring single_gappairs onto pairs\n"));
    debug(Pair_dump_list(single_gappairs,true));
    pairs = Pairpool_transfer(pairs,single_gappairs);
    *singlep = true;

  } else {
    middle_exonlength = right_exonhead-left_exonhead;
    debug(printf("Middle exon is %d - %d = %d bp in interexon region of %d bp\n",
		 right_exonhead,left_exonhead,right_exonhead-left_exonhead,new_rightgenomepos-new_leftgenomepos));
    if (middle_exonlength <= 0) {
      middle_exonprob = 0.0;
    } else {
      interexon_region = new_rightgenomepos - new_leftgenomepos;
      
      if (dual_introntype_2 == canonical_introntype) {
	middle_exonlength += DUAL_HALFCANONICAL_POINTS;
	debug(printf("Add canonical credit of %d for right intron\n",DUAL_HALFCANONICAL_POINTS));
      }
      if (dual_introntype_1 == canonical_introntype) {
	middle_exonlength += DUAL_HALFCANONICAL_POINTS;
	debug(printf("Add canonical credit of %d for left intron\n",DUAL_HALFCANONICAL_POINTS));
      }

      middle_exonprob = 1.0-pow(1.0-pow(4.0,-(double) middle_exonlength),(double) interexon_region);
      
      debug(printf("Single score = %d (%d matches).  Single canonical: %d.  Dual score = %d & %d (%d & %d matches).  Dual canonical: %d/%d.  ",
		   single_score,single_nmatches,single_canonical_p,
		   dual_score_1,dual_score_2,dual_nmatches_1,dual_nmatches_2,
		   dual_introntype_1 == canonical_introntype,
		   dual_introntype_2 == canonical_introntype));
      debug(printf("Single goodness = %d.  Dual goodness = %d.  ",
		   single_goodness,dual_goodness));
      debug(printf("Probability is %g.  ",middle_exonprob));
    }

    /* Want high threshold for accepting dual intron */
    if (dual_canonical_p == true && middle_exonprob < 0.001 &&
	(single_canonical_p == false || single_goodness <= dual_goodness)) {
      debug(printf("Dual scores win\n"));
      debug(printf("Loser: single_gappairs\n"));
      debug(Pair_dump_list(single_gappairs,true));
      debug(printf("Winner: Transferring dual_gappairs_2 onto pairs\n"));
      debug(Pair_dump_list(dual_gappairs_2,true));
      pairs = Pairpool_transfer(pairs,dual_gappairs_2);
      debug(printf("Winner: Transferring dual_gappairs_1 onto pairs\n"));
      debug(Pair_dump_list(dual_gappairs_1,true));
      pairs = Pairpool_transfer(pairs,dual_gappairs_1);
      /* *singlep = false; */
    } else {
      debug(printf("Single score wins\n"));
      debug(printf("Loser: dual_gappairs\n"));
      debug(Pair_dump_list(dual_gappairs_2,true));
      debug(Pair_dump_list(dual_gappairs_1,true));
      debug(printf("Winner: Transferring single_gappairs onto pairs\n"));
      debug(Pair_dump_list(single_gappairs,true));
      pairs = Pairpool_transfer(pairs,single_gappairs);
      *singlep = true;
    }
  }

  return pairs;
}


static List_T
distalmedial_ending5 (bool *chop_exon_p, int *dynprogindex_minor, int *dynprogindex_major, int *finalscore, List_T *pairs, 
		      int leftquerypos, int leftgenomepos, Pair_T rightpair,
#ifdef PMAP
		      char *queryaaseq_ptr,
#endif
		      char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		      int cdna_direction, int introndir, Pairpool_T pairpool,
		      Dynprog_T dynprog, int maxpeelback, int extramaterial_end,
		      int extraband_end, double defect_rate, bool end_microexons_p) {
  List_T peeled_pairs, endgappairs, continuous_gappairs_medialgap = NULL;
  int queryjump, genomejump;
  int querydp5, genomedp5, querydp3_distalgap, genomedp3_distalgap, querydp3_medialgap, genomedp3_medialgap;
  int continuous_goodness_distalgap = 0, continuous_goodness_medialgap = 0,
    nmatches, nmismatches, nopens, nindels;
  int changepoint;
  bool mismatchp = false;
  char gapcomp;

  debug(printf("\nDISTALMEDIAL_ENDING5\n"));

  querydp5 = leftquerypos + 1;
  genomedp5 = leftgenomepos + 1; /* 0 */
  querydp3_distalgap = querydp3_medialgap = rightpair->querypos - 1;
  genomedp3_distalgap = genomedp3_medialgap = rightpair->genomepos - 1;

  /* Used to peelback only half as much as for a paired gap, to save
     on dynamic programming, but not any more. */
  *pairs = peel_rightward(&mismatchp,&peeled_pairs,*pairs,&querydp3_distalgap,&genomedp3_distalgap,pairpool,maxpeelback,
			  /*throughmismatchp*/true,&endgappairs,&gapcomp,&querydp3_medialgap,&genomedp3_medialgap);
  if (endgappairs == NULL) {
    *chop_exon_p = false;
    return peeled_pairs;

  } else {
    continuous_goodness_distalgap = Pair_fracidentity_max(&changepoint,peeled_pairs,cdna_direction);
    debug(printf("continuous_goodness_distalgap at changepoint %d is %d\n",changepoint,continuous_goodness_distalgap));

    if (gapcomp == NONINTRON_COMP) {
      debug(printf("Subtracting points from continuous distal because noncanonical\n"));
      continuous_goodness_distalgap -= (CANONICAL_POINTS + CANONICAL_POINTS);
    }

    /* Solve if gap were not present */
    queryjump = querydp3_medialgap - querydp5 + 1;
    genomejump = queryjump + extramaterial_end; /* proposed */
    /* Previously, we limited genomejump = min(2*queryjump,queryjump+extramaterial_end) */

    genomedp5 = genomedp3_medialgap - genomejump + 1;
    /* Make sure we don't go past the beginning */
    if (genomedp5 < 0) {
      genomedp5 = 0;
      genomejump = genomedp3_medialgap - genomedp5 + 1;
    }

    debug(printf("Stage 3 (dir %d): traverse_ending5: Dynamic programming at 5' end (medial to gap): querydp5 = %d, querydp3 = %d, genomedp5 = %d, genomedp3 = %d\n",
		 cdna_direction,querydp5,querydp3_medialgap,genomedp5,genomedp3_medialgap));
    continuous_gappairs_medialgap = Dynprog_end5_gap(&(*dynprogindex_minor),&(*finalscore),
						     &nmatches,&nmismatches,&nopens,&nindels,dynprog,
						     &(queryseq_ptr[querydp3_medialgap]),&(queryuc_ptr[querydp3_medialgap]),
						     &(genomicseg_ptr[genomedp3_medialgap]),&(genomicuc_ptr[genomedp3_medialgap]),
						     queryjump,genomejump,querydp3_medialgap,genomedp3_medialgap,
#ifdef PMAP
						     queryaaseq_ptr,
#endif
						     cdna_direction,pairpool,extraband_end,defect_rate);
    continuous_goodness_medialgap = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

    if (continuous_goodness_distalgap > continuous_goodness_medialgap) {
      debug(printf("Continuous distal wins: %d > %d\n",continuous_goodness_distalgap,continuous_goodness_medialgap));
      *pairs = Pairpool_transfer(*pairs,endgappairs);
      *chop_exon_p = false;
      return List_truncate(peeled_pairs,changepoint);
    } else {
      debug(printf("Continuous medial wins: %d > %d\n",
		   continuous_goodness_medialgap,continuous_goodness_distalgap));
      *chop_exon_p = true;
      return continuous_gappairs_medialgap;
    }
  }

}


static List_T
extend_ending5 (int *dynprogindex_minor, int *dynprogindex_major, int *finalscore, List_T *pairs, 
		int leftquerypos, int leftgenomepos, Pair_T rightpair,
#ifdef PMAP
		char *queryaaseq_ptr,
#endif
		char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		int cdna_direction, int introndir, Pairpool_T pairpool,
		Dynprog_T dynprog, int maxpeelback, int extramaterial_end,
		int extraband_end, double defect_rate, bool end_microexons_p) {
  List_T continuous_gappairs_distalgap = NULL, peeled_pairs;
  int queryjump, genomejump;
  int querydp5, genomedp5, querydp3_distalgap, genomedp3_distalgap;
  int nmatches, nmismatches, nopens, nindels;
  bool mismatchp = false;
  char gapcomp;


  debug(printf("\nEXTEND_ENDING5\n"));

  querydp5 = leftquerypos + 1;
  genomedp5 = leftgenomepos + 1; /* 0 */
  querydp3_distalgap = rightpair->querypos - 1;
  genomedp3_distalgap = rightpair->genomepos - 1;

  /* Used to peelback only half as much as for a paired gap, to save
     on dynamic programming, but not any more. */
#if 1
  /* Don't peelback on extension */
  *pairs = peel_rightward(&mismatchp,&peeled_pairs,*pairs,&querydp3_distalgap,&genomedp3_distalgap,pairpool,maxpeelback,
			  /*throughmismatchp*/true,/*endgappairs*/NULL,&gapcomp,/*querydp3_medialgap*/NULL,/*genomedp3_medialgap*/NULL);
#endif
  
  queryjump = querydp3_distalgap - querydp5 + 1;
  genomejump = queryjump + extramaterial_end; /* proposed */
  /* Previously, we limited genomejump = min(2*queryjump,queryjump+extramaterial_end) */

  genomedp5 = genomedp3_distalgap - genomejump + 1;
  /* Make sure we don't go past the beginning */
  if (genomedp5 < 0) {
    genomedp5 = 0;
    genomejump = genomedp3_distalgap - genomedp5 + 1;
  }
  debug(printf("Stage 3 (dir %d), extend_ending5: Dynamic programming at 5' end (distal to gap): querydp5 = %d, querydp3 = %d, genomedp5 = %d, genomedp3 = %d\n",
	       cdna_direction,querydp5,querydp3_distalgap,genomedp5,genomedp3_distalgap));

  continuous_gappairs_distalgap = Dynprog_end5_gap(&(*dynprogindex_minor),&(*finalscore),
						   &nmatches,&nmismatches,&nopens,&nindels,dynprog,
						   &(queryseq_ptr[querydp3_distalgap]),&(queryuc_ptr[querydp3_distalgap]),
						   &(genomicseg_ptr[genomedp3_distalgap]),&(genomicuc_ptr[genomedp3_distalgap]),
						   queryjump,genomejump,querydp3_distalgap,genomedp3_distalgap,
#ifdef PMAP
						   queryaaseq_ptr,
#endif
						   cdna_direction,pairpool,extraband_end,defect_rate);
  debug(printf("  finalscore: %d\n",*finalscore));
  if (*finalscore < 0) {
    return (List_T) NULL;
  } else {
    return continuous_gappairs_distalgap;
  }
}


static List_T
distalmedial_ending3 (bool *chop_exon_p, int *dynprogindex_minor, int *dynprogindex_major, int *finalscore, List_T *path, 
		      Pair_T leftpair, int rightquerypos, int rightgenomepos, int querylength, Genomicpos_T genomiclength,
#ifdef PMAP
		      char *queryaaseq_ptr,
#endif
		      char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		      int cdna_direction, int introndir,
		      Pairpool_T pairpool, Dynprog_T dynprog, int maxpeelback, int extramaterial_end,
		      int extraband_end, double defect_rate, bool end_microexons_p) {
  List_T peeled_path, endgappairs, continuous_gappairs_medialgap = NULL;
  int queryjump, genomejump;
  int querydp5_distalgap, genomedp5_distalgap, querydp3, genomedp3, querydp5_medialgap, genomedp5_medialgap;
  int continuous_goodness_distalgap = 0, continuous_goodness_medialgap = 0,
    nmatches, nmismatches, nopens, nindels;
  int changepoint;
  bool mismatchp = false;
  char gapcomp;

  debug(printf("\nTRAVERSE_ENDING3\n"));
  
  querydp5_distalgap = leftpair->querypos + 1;
  genomedp5_distalgap = leftpair->genomepos + 1;
  if (leftpair->cdna == ' ') querydp5_distalgap--;
  if (leftpair->genome == ' ') genomedp5_distalgap--;
  querydp5_medialgap = querydp5_distalgap;
  genomedp5_medialgap = genomedp5_distalgap;
  querydp3 = rightquerypos - 1;
  genomedp3 = rightgenomepos - 1;

  /* Used to peelback only half as much as for a paired gap, to save
     on dynamic programming, but not any more. */
  *path = peel_leftward(&mismatchp,&peeled_path,*path,&querydp5_distalgap,&genomedp5_distalgap,pairpool,maxpeelback,
			/*throughmismatchp*/true,&endgappairs,&gapcomp,&querydp5_medialgap,&genomedp5_medialgap);
  
  if (endgappairs == NULL) {
    *chop_exon_p = false;
    return peeled_path;

  } else {
    continuous_goodness_distalgap = Pair_fracidentity_max(&changepoint,peeled_path,cdna_direction);
    debug(printf("continuous_goodness_distalgap at changepoint %d is %d\n",changepoint,continuous_goodness_distalgap));

    if (gapcomp == NONINTRON_COMP) {
      debug(printf("Subtracting points from continuous distal because noncanonical\n"));
      continuous_goodness_distalgap -= (CANONICAL_POINTS + CANONICAL_POINTS);
    }

    /* Solve if gap were not present */
    queryjump = querydp3 - querydp5_medialgap + 1;
    genomejump = queryjump + extramaterial_end; /* proposed */
    /* Previously, we limited genomejump = min(2*queryjump,queryjump+extramaterial_end) */

    genomedp3 = genomedp5_medialgap + genomejump - 1;
    /* Make sure we don't go past the end */
    if (genomedp3 > genomiclength - 1) {
      genomedp3 = genomiclength - 1;
      genomejump = genomedp3 - genomedp5_medialgap + 1;
    }
    
    debug(printf("Stage 3 (dir %d): distalmedial_ending3: Dynamic programming at 3' end (medial to gap): querydp5 = %d, querydp3 = %d, genomedp5 = %d, genomedp3 = %d\n",
		 cdna_direction,querydp5_medialgap,querydp3,genomedp5_medialgap,genomedp3));
    continuous_gappairs_medialgap = Dynprog_end3_gap(&(*dynprogindex_minor),&(*finalscore),
						     &nmatches,&nmismatches,&nopens,&nindels,dynprog,
						     &(queryseq_ptr[querydp5_medialgap]),&(queryuc_ptr[querydp5_medialgap]),
						     &(genomicseg_ptr[genomedp5_medialgap]),&(genomicuc_ptr[genomedp5_medialgap]),
						     queryjump,genomejump,querydp5_medialgap,genomedp5_medialgap,
#ifdef PMAP
						     queryaaseq_ptr,
#endif
						     cdna_direction,pairpool,extraband_end,defect_rate);
    continuous_goodness_medialgap = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

    if (continuous_goodness_distalgap > continuous_goodness_medialgap) {
      debug(printf("Continuous distal wins: %d > %d\n",continuous_goodness_distalgap,continuous_goodness_medialgap));
      *path = Pairpool_transfer(*path,endgappairs);
      *chop_exon_p = false;
      return List_truncate(peeled_path,changepoint);
    } else {
      debug(printf("Continuous medial wins: %d > %d\n",continuous_goodness_medialgap,continuous_goodness_distalgap));
      *chop_exon_p = true;
      return List_reverse(continuous_gappairs_medialgap);
    }
  }

}


static List_T
extend_ending3 (int *dynprogindex_minor, int *dynprogindex_major, int *finalscore, List_T *path, 
		Pair_T leftpair, int rightquerypos, int rightgenomepos, int querylength, Genomicpos_T genomiclength,
#ifdef PMAP
		char *queryaaseq_ptr,
#endif
		char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		int cdna_direction, int introndir,
		Pairpool_T pairpool, Dynprog_T dynprog, int maxpeelback, int extramaterial_end,
		int extraband_end, double defect_rate, bool end_microexons_p) {
  List_T continuous_gappairs_distalgap = NULL, peeled_path;
  int queryjump, genomejump;
  int querydp5_distalgap, genomedp5_distalgap, querydp3, genomedp3;
  int nmatches, nmismatches, nopens, nindels;
  bool mismatchp = false;
  char gapcomp;

  debug(printf("\nEXTEND_ENDING3\n"));
  
  querydp5_distalgap = leftpair->querypos + 1;
  genomedp5_distalgap = leftpair->genomepos + 1;
  if (leftpair->cdna == ' ') querydp5_distalgap--;
  if (leftpair->genome == ' ') genomedp5_distalgap--;
  querydp3 = rightquerypos - 1;
  genomedp3 = rightgenomepos - 1;
  debug(printf("Set dynprog 3 end to be querydp3 = %d, genomedp3 = %d\n",querydp3,genomedp3));

  /* Used to peelback only half as much as for a paired gap, to save
     on dynamic programming, but not any more. */
#if 1
  /* Don't peelback on extension */
  *path = peel_leftward(&mismatchp,&peeled_path,*path,&querydp5_distalgap,&genomedp5_distalgap,pairpool,maxpeelback,
			/*throughmismatchp*/true,/*endgappairs*/NULL,&gapcomp,/*querydp5_medialgap*/NULL,/*genomedp5_medialgap*/NULL);
#endif

  queryjump = querydp3 - querydp5_distalgap + 1;
  genomejump = queryjump + extramaterial_end; /* proposed */
  /* Previously, we limited genomejump = min(2*queryjump,queryjump+extramaterial_end) */

  genomedp3 = genomedp5_distalgap + genomejump - 1;
  /* Make sure we don't go past the end */
  if (genomedp3 > genomiclength - 1) {
    genomedp3 = genomiclength - 1;
    genomejump = genomedp3 - genomedp5_distalgap + 1;
  }
  debug(printf("Stage 3 (dir %d), extend_ending3: Dynamic programming at 3' end (distal to gap): querydp5 = %d, querydp3 = %d, genomedp5 = %d, genomedp3 = %d\n",
	       cdna_direction,querydp5_distalgap,querydp3,genomedp5_distalgap,genomedp3));
  
  continuous_gappairs_distalgap = Dynprog_end3_gap(&(*dynprogindex_minor),&(*finalscore),
						   &nmatches,&nmismatches,&nopens,&nindels,dynprog,
						   &(queryseq_ptr[querydp5_distalgap]),&(queryuc_ptr[querydp5_distalgap]),
						   &(genomicseg_ptr[genomedp5_distalgap]),&(genomicuc_ptr[genomedp5_distalgap]),
						   queryjump,genomejump,querydp5_distalgap,genomedp5_distalgap,
#ifdef PMAP
						   queryaaseq_ptr,
#endif
						   cdna_direction,pairpool,extraband_end,defect_rate);
  debug(printf("  finalscore: %d\n",*finalscore));
  if (*finalscore < 0) {
    return (List_T) NULL;
  } else {
    return List_reverse(continuous_gappairs_distalgap);
  }
}


#if 0
static List_T
traverse_dual_break_old (bool *new_exon_p, int *dynprogindex_minor, int *dynprogindex_major,
		     List_T pairs, List_T *path, Pair_T leftpair, Pair_T rightpair,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogR, 
		     int maxpeelback, int extraband_end, double defect_rate) {
  List_T gappairs, peeled_pairs = NULL, peeled_path = NULL, p;
  Pair_T pair;
  int g, gg, c, beststart_g, beststart_c, bestend_g, bestend_c, localstart_g, localstart_c, localend_g, localend_c,
    start_g, start_c;
  int maxscore = 0, localscore, score;
  int queryjump, genomejump;
  int querydp5, genomedp5, querydp3, genomedp3;
  int finalscore;
  bool mismatchp = false;
  char gapcomp;

  debug(printf("\nTRAVERSE_DUAL_BREAK\n"));
  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  if (leftpair->cdna == ' ') querydp5--;
  if (leftpair->genome == ' ') genomedp5--;
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;

  pairs = peel_rightward(&mismatchp,&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
			 /*throughmismatchp*/true,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);
  *path = peel_leftward(&mismatchp,&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
			/*throughmismatchp*/true,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);

  /* Try scanning for match */
  for (g = genomedp5; g < genomedp3 - (querydp3 - querydp5); g++) {
    score = 0;
    start_g = g;
    start_c = querydp5;
    localscore = 0;
    for (c = querydp5, gg = g; c <= querydp3; c++, gg++) {
      if (queryuc_ptr[c] == genomicuc_ptr[gg]) {
	score += 1;
      } else {
	score -= 1;
      }
      if (score <= 0) {
	score = 0;
	start_g = gg;
	start_c = c;
      } else if (score > localscore) {
	localscore = score;
	localstart_g = start_g;
	localstart_c = start_c;
	localend_g = gg;
	localend_c = c;
      }
    }
    if (localscore > maxscore) {
      maxscore = localscore;
      beststart_g = localstart_g;
      beststart_c = localstart_c;
      bestend_g = localend_g;
      bestend_c = localend_c;
    }
  }

  if (maxscore > 0) {
    debug(printf("traverse_dual_break found score of %d from %d..%d, %d..%d\n",
		 maxscore,beststart_c,bestend_c,beststart_g,bestend_g));
    for (p = peeled_pairs; p != NULL; p = p->rest) {
      pair = (Pair_T) p->first;
      if (pair->querypos > bestend_c) {
	debug(printf("Putting back originally right-peeled at %d,%d\n",pair->querypos,pair->genomepos));
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      } else {
	debug(printf("Not putting back originally right-peeled at %d,%d\n",pair->querypos,pair->genomepos));
      }
    }

    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);
    for (g = bestend_g, c = bestend_c; g >= beststart_g; g--, c--) {
      if (queryuc_ptr[c] == genomicuc_ptr[g]) {
	pairs = Pairpool_push(pairs,pairpool,c,g,queryuc_ptr[c],MATCH_COMP,genomicuc_ptr[g],*dynprogindex_major);
      } else {
	pairs = Pairpool_push(pairs,pairpool,c,g,queryuc_ptr[c],MISMATCH_COMP,genomicuc_ptr[g],*dynprogindex_major);
      }
    }
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);

    for (p = peeled_path; p != NULL; p = p->rest) {
      pair = (Pair_T) p->first;
      if (pair->querypos < beststart_c) {
	debug(printf("Putting back originally left-peeled at %d,%d\n",pair->querypos,pair->genomepos));
	*path = Pairpool_push_existing(*path,pairpool,pair);
      } else {
	debug(printf("Not putting back originally left-peeled at %d,%d\n",pair->querypos,pair->genomepos));
      }
    }

    (*dynprogindex_major) -= 1;
    *new_exon_p = true;

  } else {
    queryjump = querydp3 - querydp5 + 1;
    genomejump = genomedp3 - genomedp5 + 1; /* don't use extramaterial here */

    gappairs = Dynprog_dual_break(&(*dynprogindex_minor),&finalscore,
				  dynprogL,dynprogR,
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				  &(queryseq_ptr[querydp3]),&(queryuc_ptr[querydp3]),
				  &(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
				  queryjump,genomejump,querydp5,genomedp5,querydp3,genomedp3,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  cdna_direction,pairpool,extraband_end,defect_rate);

    debug(printf("gappairs on ends of dual break:\n"));
    debug(Pair_dump_list(gappairs,true));
    debug(printf("end\n"));
    
    if (gappairs == NULL) {
      pairs = Pairpool_transfer(pairs,peeled_pairs);
      *path = Pairpool_transfer(*path,peeled_path);
      pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);
    } else {
      pairs = Pairpool_transfer(pairs,gappairs);
    }
  }

  return pairs;
}
#endif


static List_T
traverse_dual_break (int *dynprogindex_minor, int *dynprogindex_major,
		     List_T pairs, List_T *path, Pair_T leftpair, Pair_T rightpair,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogR, 
		     int maxpeelback, Oligoindex_T *oligoindices, int noligoindices, Diagpool_T diagpool,
		     int sufflookback, int nsufflookback, int maxintronlen_bound) {
  List_T gappairs, peeled_pairs = NULL, peeled_path = NULL;
  int querydp5, genomedp5, querydp3, genomedp3;
  char gapcomp;
  bool mismatchp;

  debug(printf("\nTRAVERSE_DUAL_BREAK\n"));
  querydp5 = leftpair->querypos + 1;
  genomedp5 = leftpair->genomepos + 1;
  if (leftpair->cdna == ' ') querydp5--;
  if (leftpair->genome == ' ') genomedp5--;
  querydp3 = rightpair->querypos - 1;
  genomedp3 = rightpair->genomepos - 1;

  pairs = peel_rightward(&mismatchp,&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
			 /*throughmismatchp*/true,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);
  *path = peel_leftward(&mismatchp,&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
			/*throughmismatchp*/true,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,/*genomedp_medialgap*/NULL);

  debug(printf("genome %.*s\n",genomedp3-genomedp5+1,&(genomicseg_ptr[genomedp5])));
  debug(printf("query %.*s\n",querydp3-querydp5+1,&(queryseq_ptr[querydp5])));

  gappairs = Stage2_compute(&(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
			    querydp3-querydp5+1,/*query_offset*/querydp5,
			    &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
			    genomedp3-genomedp5+1,/*genomic_offset*/genomedp5,
			    oligoindices,noligoindices,pairpool,diagpool,
			    sufflookback,nsufflookback,maxintronlen_bound,
			    /*localp should be false*/true,/*skip_repetitive_p*/false,/*debug_graphic_p*/false,
			    /*diagnosticp*/false,/*stopwatch*/NULL);
  
  debug(printf("Internal stage2 result:\n"));
  debug(Pair_dump_list(gappairs,true));

  if (gappairs == NULL) {
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    pairs = Pairpool_push_gapholder(pairs,pairpool,/*queryjump*/UNKNOWNJUMP,/*genomejump*/UNKNOWNJUMP);
  } else {
    pairs = Pairpool_transfer(pairs,gappairs);
  }

  return pairs;
}


/************************************************************************
 *   End of traversal functions
 ************************************************************************/

static List_T
build_dual_breaks (List_T path, int *dynprogindex_minor, int *dynprogindex_major,
#ifdef PMAP
		   char *queryaaseq_ptr,
#endif
		   char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		   int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogR, 
		   int maxpeelback, Oligoindex_T *oligoindices, int noligoindices, Diagpool_T diagpool,
		   int sufflookback, int nsufflookback, int maxintronlen_bound) {

  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;

  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == true && pair->comp == DUALBREAK_COMP) {
      if (path == NULL || pairs == NULL) {
	debug(printf("Observed a dual break at the end of the alignment\n"));
      } else {
	leftpair = path->first;
	rightpair = pairs->first;
	debug(printf("Observed a dual break at %d..%d with queryjump = %d, genomejump = %d\n",
		     leftpair->querypos,rightpair->querypos,pair->queryjump,pair->genomejump));

	pairs = traverse_dual_break(&(*dynprogindex_minor),&(*dynprogindex_major),
				    pairs,&path,leftpair,rightpair,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				    cdna_direction,pairpool,dynprogL,dynprogR,maxpeelback,
				    oligoindices,noligoindices,diagpool,sufflookback,nsufflookback,
				    maxintronlen_bound);
      }
    } else {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
    }
  }

  debug(printf("After build_dual_breaks:\n"));
  debug(Pair_dump_list(pairs,true));
  return pairs;
}


  
/* Note: querypos is actually indexsize nt to the left of the last nt match.
   
	||||||||********   X  X	 XX X	X
	       ^	 <- queryjump->	 ^
	    leftquerypos		 rightquerypos

		<-     querydpspan     ->
*/
static List_T
build_path_end3 (bool *chop_exon_p, int *dynprogindex_minor, int *dynprogindex_major,
		 List_T path, int querylength, Genomicpos_T genomiclength,
#ifdef PMAP
		 char *queryaaseq_ptr,
#endif
		 char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		 int cdna_direction, int introndir, int maxpeelback, int nullgap,
		 int extramaterial_end, int extramaterial_paired, int extraband_end,
		 double defect_rate, Pairpool_T pairpool, Dynprog_T dynprogL, bool end_microexons_p,
		 bool extendp) {
  List_T gappairs;
  Pair_T leftpair;
  int queryjump, genomejump, rightquerypos;
  int finalscore;

  if (path == NULL) {
    return (List_T) NULL;
  } else {
    leftpair = path->first;
  }
  debug(printf("Stage 3 (dir %d): 3' end: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d\n",
	       cdna_direction,leftpair->querypos,querylength,leftpair->genomepos,genomiclength));
  if (leftpair->querypos < 0) {
    return (List_T) NULL;
    /* abort(); */
  }

  queryjump = querylength - leftpair->querypos - 1;
  genomejump = genomiclength - leftpair->genomepos - 1;
  if (leftpair->cdna == ' ') queryjump++;
  if (leftpair->genome == ' ') genomejump++;

  /* Note difference with 5' case.  We use queryjump+1 here instead of queryjump and genomejump */
  /* Do use nullgap here.  Truncating back to entire exon can slow algorithm down significantly. */
  if (/* 0 && */ queryjump+1 > nullgap) {
    rightquerypos = leftpair->querypos + nullgap + 1;
    debug(printf("Since queryjump+1 %d > nullgap %d, setting rightquerypos %d = %d + %d + 1\n",
		 queryjump+1,nullgap,rightquerypos,leftpair->querypos,nullgap));
  } else {
    rightquerypos = querylength;
  }

  if (extendp == true) {
    *chop_exon_p = false;
    gappairs = extend_ending3(&(*dynprogindex_minor),&(*dynprogindex_major),&finalscore,&path,
			      leftpair,rightquerypos,genomiclength,querylength,genomiclength,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			      cdna_direction,introndir,pairpool,dynprogL,maxpeelback,
			      extramaterial_end,extraband_end,defect_rate,end_microexons_p);
  } else {
    gappairs = distalmedial_ending3(&(*chop_exon_p),&(*dynprogindex_minor),&(*dynprogindex_major),&finalscore,&path,
				    leftpair,rightquerypos,genomiclength,querylength,genomiclength,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				    cdna_direction,introndir,pairpool,dynprogL,maxpeelback,
				    extramaterial_end,extraband_end,defect_rate,end_microexons_p);
  }

  debug(Pair_dump_list(gappairs,true));

  path = Pairpool_transfer(path,gappairs);

  return path;
}


/* Schematic:
   
       <- queryjump ->
       X  X  XX X   X ********||||||||
      ^		      ^
   leftquerypos	      rightquerypos

       <-    querydpspan    ->

*/
static List_T
build_pairs_end5 (bool *chop_exon_p, int *dynprogindex_minor, int *dynprogindex_major, List_T pairs,
#ifdef PMAP
		  char *queryaaseq_ptr,
#endif
		  char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		  int cdna_direction, int introndir, int maxpeelback, int nullgap,
		  int extramaterial_end, int extramaterial_paired, int extraband_end,
		  double defect_rate, Pairpool_T pairpool, Dynprog_T dynprogR, bool end_microexons_p,
		  bool extendp) {
  List_T gappairs;
  Pair_T rightpair;
  int queryjump, leftquerypos;
  int finalscore;
  /* int genomejump */

  if (pairs == NULL) {
    return (List_T) NULL;
  } else {
    rightpair = pairs->first;
  }
  debug(printf("Stage 3 (dir %d): 5' end: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d\n",
	       cdna_direction,-1,rightpair->querypos,-1,rightpair->genomepos));
  if (rightpair->querypos < 0) {
    return (List_T) NULL;
    /* abort(); */
  }

  queryjump = rightpair->querypos; /* - leftquerypos (-1) - 1 */
  /* genomejump = rightpair->genomepos; */  /* - leftgenomepos (-1) - 1 */

  /* Note difference with 3' case.  We use queryjump here instead of queryjump+1 */
  /* Do use nullgap here.  Truncating back to entire exon can slow algorithm significantly. */
  if (/*0 && */ queryjump > nullgap) {
    leftquerypos = rightpair->querypos - nullgap - 1;
  } else {
    leftquerypos = -1;
  }

  if (extendp == true) {
    *chop_exon_p = false;
    gappairs = extend_ending5(&(*dynprogindex_minor),&(*dynprogindex_major),&finalscore,
			      &pairs,leftquerypos,/*leftgenomepos*/-1,rightpair,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			      cdna_direction,introndir,pairpool,dynprogR,maxpeelback,
			      extramaterial_end,extraband_end,defect_rate,end_microexons_p);
  } else {
    gappairs = distalmedial_ending5(&(*chop_exon_p),&(*dynprogindex_minor),&(*dynprogindex_major),&finalscore,
				    &pairs,leftquerypos,/*leftgenomepos*/-1,rightpair,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				    cdna_direction,introndir,pairpool,dynprogR,maxpeelback,
				    extramaterial_end,extraband_end,defect_rate,end_microexons_p);
  }


  debug(Pair_dump_list(gappairs,true));

  pairs = Pairpool_transfer(pairs,gappairs);

  return pairs;
}


static List_T
build_pairs_singles (int *dynprogindex, List_T path,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, int maxpeelback, int nullgap,
		     int extramaterial_paired, int extraband_single, double defect_rate, bool widebandp,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;
  bool filledp;

  debug(printf("\n** Starting build_pairs_singles\n"));
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == false) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

#if 0
    } else if (pair->queryjump == 1 && pair->genomejump == 1) {
      /* Special case: Single mismatch */
      leftquerypos = pair->querypos;
      leftgenomepos = pair->genomepos;

      c1 = queryseq_ptr[leftquerypos+1];
      c2 = genomicseg_ptr[leftgenomepos+1];
      if (queryuc_ptr[leftquerypos+1] == genomicuc_ptr[leftgenomepos+1]) {
	/* Possible for a stretch of the same nucleotide to have a mismatch */
	debug(printf("Stage 3: Single match at %d %d: %c | %c\n",
		     leftquerypos+1,leftgenomepos+1,c1,c2));
	pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
			      c1,MATCH_COMP,c2,/*dynprogindex*/0);
      } else {
	debug(printf("Stage 3: Single mismatch at %d %d: %c   %c\n",
		     leftquerypos+1,leftgenomepos+1,c1,c2));
	pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
			      c1,MISMATCH_COMP,c2,/*dynprogindex*/0);
      }
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      
    } else if (pair->queryjump == 1 && pair->genomejump == 0) {
      /* Special case: Single cDNA insert.  This will create an apparent genomejump of -1. */
      leftquerypos = pair->querypos;
      leftgenomepos = pair->genomepos;

      debug(printf("Stage 3: Single cDNA insert at %d %d: %c\n",
		   leftquerypos+1,leftgenomepos+1,queryseq_ptr[leftquerypos+1]));
      /* Advance querypos to next position */
      pairs = Pairpool_push(pairs,pairpool,leftquerypos+1,leftgenomepos+1,
			    queryseq_ptr[leftquerypos+1],INDEL_COMP,' ',/*dynprogindex*/0);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#endif

    } else if (pair->queryjump > nullgap) {
      /* Large gap.  Do nothing */
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->queryjump > pair->genomejump + EXTRAQUERYGAP) {
      /* cDNA insertion.  Do nothing */
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->genomejump > pair->queryjump + SINGLESLEN) {
      /* Intron.  Do nothing */
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (path == NULL || pairs == NULL) {
      fprintf(stderr,"Single gap at beginning or end of alignment\n");
      abort();

    } else {
      /* Guarantees: queryjump <= nullgap && genomejump < queryjump - EXTRAQUERYGAP &&
	 genomejump <= queryjump + MININTRONLEN, meaning that score matrix is nearly square */
      leftpair = path->first;
      rightpair = pairs->first;
	
      debug(printf("Stage 3 (dir %d): Traversing single gap: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d\n",
		   cdna_direction,leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos));
      pairs = traverse_single_gap(&filledp,&(*dynprogindex),pairs,&path,leftpair,rightpair,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				  cdna_direction,pairpool,dynprogM,maxpeelback,extraband_single,defect_rate,
				  /*forcep*/true,widebandp);
      /* forcep needs to be true here to avoid subsequent anomalies in building dualintrons, e.g., XM_376610.2_mRNA on 7:127885572..127888991 */
      if (filledp == true) {
	/* Discard the gap */
	debug(printf("Discarding gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
      } else {
	/* Replace the gap */
	debug(printf("Replacing gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }
    }
  }

  return pairs;
}


static List_T
build_pairs_dualintrons (bool *singlep, int *dynprogindex, List_T path,
#ifdef PMAP
			 char *queryaaseq_ptr,
#endif
			 char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
			 int cdna_direction, int maxpeelback, int nullgap,
			 int extramaterial_paired, int extraband_paired, double defect_rate,
			 Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T pairs = NULL, midexon_pairs;
  Pair_T pair, leftpair, midleftpair, midpair, midrightpair, rightpair;
  int midquerypos, midgenomepos;
  bool exonp;

  debug(printf("\n** Starting build_pairs_dualintrons\n"));
  *singlep = false;
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    
    if (pair->gapp == false) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->queryjump > nullgap) {
      /* Large gap.  Do nothing */
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->queryjump > pair->genomejump + EXTRAQUERYGAP) {
      /* cDNA insertion.  Do nothing */
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->genomejump <= pair->queryjump + MININTRONLEN) {
      /* Single gap.  Do nothing */
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else {
      midpair = path->first;
      if (midpair->shortexonp == false) {
	/* Long exon; do nothing */
	pairs = Pairpool_push_existing(pairs,pairpool,pair);

      } else {
	/* Short exon */
	debug(printf("I see a short exon at %d...crossing\n",midpair->querypos));
	path = Pairpool_pop(path,&midpair);
	midexon_pairs = Pairpool_push_existing(NULL,pairpool,midpair);
	midrightpair = midpair;

	exonp = true;
	while (path != NULL && exonp) {
	  path = Pairpool_pop(path,&midpair);
	  if (midpair->gapp == true) {
	    exonp = false;
	  } else {
	    midexon_pairs = Pairpool_push_existing(midexon_pairs,pairpool,midpair);
	  }
	}
	debug(printf("Finished crossing a short exon\n"));

	if (path == NULL) {
	  /* Short exon is the first one.  Process it. */
	  pairs = Pairpool_push_existing(pairs,pairpool,pair); /* initial gap */
	  pairs = Pairpool_transfer(pairs,List_reverse(midexon_pairs));

	} else {
	  /* Perform dual intron gap */
	  midleftpair = midexon_pairs->first;
	  midgenomepos = (midleftpair->genomepos + midrightpair->genomepos)/2;
	  midquerypos = midrightpair->querypos - (midrightpair->genomepos - midgenomepos);
	  leftpair = path->first;
	  rightpair = pairs->first;
	  
	  debug(printf("Stage 3 (dir %d): Traversing dual intron gap: leftquerypos = %d, midquerypos = %d, rightquerypos = %d, leftgenomepos = %d, midgenomepos = %d, rightgenomepos = %d\n",
		       cdna_direction,leftpair->querypos,midquerypos,rightpair->querypos,
		       leftpair->genomepos,midgenomepos,rightpair->genomepos));
	  
	  pairs = traverse_dual_genome_gap(&(*singlep),&(*dynprogindex),pairs,&path,leftpair,rightpair,
					   midquerypos,midgenomepos,
#ifdef PMAP
					   queryaaseq_ptr,
#endif
					   queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
					   cdna_direction,pairpool,dynprogL,dynprogM,dynprogR,
					   maxpeelback,nullgap,extramaterial_paired,extraband_paired,
					   defect_rate);
	}
      }
    }
  }

  return pairs;
}  


static List_T
build_pairs_introns (bool *shiftp, bool *incompletep, 
		     int *nintrons, int *nnonintrons, int *intronlen, int *nonintronlen, 
		     int *dynprogindex_minor, int *dynprogindex_major, List_T path, 
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, int maxpeelback, int nullgap, int extramaterial_paired, 
		     int extraband_single, int extraband_paired, double defect_rate,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		     bool finalp, bool forcep) {
  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;
  bool filledp;
  int minintronlen;

  debug(printf("\n** Starting build_pairs_introns\n"));

  if (finalp == true) {
    minintronlen = MININTRONLEN_FINAL;
  } else {
    minintronlen = MININTRONLEN;
  }

  *shiftp = *incompletep = false;
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == false) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->queryjump > nullgap) {
      debug(leftpair = path->first;
	    rightpair = pairs->first;
	    printf("Stage 3 (dir %d): Adding large gap: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d\n",
		   cdna_direction,leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos)); 
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->queryjump > pair->genomejump + EXTRAQUERYGAP) {
      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Stage 3 (dir %d): Traversing cDNA gap: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d\n",
		   cdna_direction,leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos));
      pairs = traverse_cdna_gap(&filledp,&(*incompletep),&(*dynprogindex_minor),&(*dynprogindex_major),
				pairs,&path,leftpair,rightpair,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				pairpool,dynprogL,dynprogM,dynprogR,maxpeelback,extramaterial_paired,
				extraband_paired,extraband_single,defect_rate,forcep);

      if (filledp == true) {
	/* Discard gap */
	debug(printf("Discarding gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
      } else {
	/* Replace the gap */
	debug(printf("Replacing gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }

    } else if (pair->genomejump > pair->queryjump + minintronlen) {
      /* Previously was 2*MININTRONLEN, and comment said needed space for two introns */
      /* We will make the score matrices nearly square */
      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Stage 3 (dir %d): Traversing paired gap: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d\n",
		   cdna_direction,leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos));
      pairs = traverse_genome_gap(&filledp,&(*shiftp),&(*dynprogindex_minor),&(*dynprogindex_major),
				  &(*nintrons),&(*nnonintrons),&(*intronlen),&(*nonintronlen),pairs,&path,leftpair,rightpair,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				  pairpool,dynprogL,dynprogM,dynprogR,maxpeelback,extramaterial_paired,
				  extraband_paired,extraband_single,defect_rate,/*forcep*/true,finalp);
      /* Do forcep, because adding large gap is not a good solution */

      if (filledp == true) {
	/* Discard the gap */
	debug(printf("Discarding gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
      } else {
	/* Replace the gap */
	debug(printf("Replacing gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }

    } else if (pair->genomejump > pair->queryjump + SINGLESLEN) {
      /* Intron length shorter than MININTRONLEN_FINAL.  Just replace the gap */
      debug(printf("Short intron; not candidate for final calculation.  Replacing gap ");
	    Pair_dump_one(pair,true);
	    printf("\n"));
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else {
      /* Single gap; force fill */
      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Stage 3 (dir %d): Traversing single gap: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d.  queryjump = %d, genomejump = %d\n",
		   cdna_direction,leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos,
		   pair->queryjump,pair->genomejump));
      pairs = traverse_single_gap(&filledp,&(*dynprogindex_minor),pairs,&path,leftpair,rightpair,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				  cdna_direction,pairpool,dynprogM,maxpeelback,extraband_single,defect_rate,
				  /*forcep*/false,/*widebandp*/true);

      if (filledp == true) {
	/* Discard the gap */
	debug(printf("Discarding gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
      } else {
	/* Replace the gap */
	debug(printf("Replacing gap ");
	      Pair_dump_one(pair,true);
	      printf("\n"));
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }
    }
  }

  return pairs;
}  


static List_T
score_introns (int *totalintronscore, List_T path, int cdna_direction, int maxpeelback,
	       int nullgap, Pairpool_T pairpool) {
  List_T pairs = NULL, peeled_pairs, peeled_path;
  Pair_T pair, leftpair, rightpair;
  int minintronlen;
  int querydp5, genomedp5, querydp3, genomedp3;
  int nmatches, nmismatches;
  int nunknowns, qopens, qindels, topens, tindels, ncanonical, nsemicanonical, nnoncanonical;
  bool mismatch_leftward_p, mismatch_rightward_p;
  char gapcomp;

  debug(printf("\n** Starting score_introns\n"));

  minintronlen = MININTRONLEN_FINAL;

  *totalintronscore = 0;
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == false) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->queryjump > nullgap) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->queryjump > pair->genomejump + EXTRAQUERYGAP) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->genomejump > pair->queryjump + minintronlen) {
      leftpair = path->first;
      rightpair = pairs->first;

      querydp5 = leftpair->querypos + 1;
      genomedp5 = leftpair->genomepos + 1;
      if (leftpair->cdna == ' ') querydp5--;
      if (leftpair->genome == ' ') genomedp5--;
      querydp3 = rightpair->querypos - 1;
      genomedp3 = rightpair->genomepos - 1;

      pairs = peel_rightward(&mismatch_rightward_p,&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
			     /*throughmismatchp*/false,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,
			     /*genomedp_medialgap*/NULL);
      path = peel_leftward(&mismatch_leftward_p,&peeled_path,path,&querydp5,&genomedp5,pairpool,maxpeelback,
			   /*throughmismatchp*/false,/*endgappairs*/NULL,&gapcomp,/*querydp_medialgap*/NULL,
			   /*genomedp_medialgap*/NULL);

      debug(printf("pair->comp = %c",pair->comp));
      if (cdna_direction == +1) {
	if (pair->comp == FWD_CANONICAL_INTRON_COMP) {
	  debug(printf(" => add %d points",CANONICAL_POINTS));
	  *totalintronscore += CANONICAL_POINTS;
	}
      } else if (cdna_direction == -1) {
	if (pair->comp == REV_CANONICAL_INTRON_COMP) {
	  debug(printf(" => add %d points\n",CANONICAL_POINTS));
	  *totalintronscore += CANONICAL_POINTS;
	}
      }
      debug(printf("\n"));

      Pair_fracidentity(&nmatches,&nunknowns,&nmismatches,&qopens,&qindels,&topens,&tindels,
			&ncanonical,&nsemicanonical,&nnoncanonical,peeled_pairs,cdna_direction);
      *totalintronscore += MISMATCH*nmismatches + QOPEN*(qopens + qindels) + TOPEN*(topens + tindels);
      debug(printf("%d matches, %d nmismatches, %d+%d qgaps, %d+%d tgaps => intronscore is %d\n",
		   nmatches,nmismatches,qopens,qindels,topens,tindels,
		   MISMATCH*nmismatches + QOPEN*(qopens + qindels) + TOPEN*(topens + tindels)));

      Pair_fracidentity(&nmatches,&nunknowns,&nmismatches,&qopens,&qindels,&topens,&tindels,
			&ncanonical,&nsemicanonical,&nnoncanonical,peeled_path,cdna_direction);
      *totalintronscore += MISMATCH*nmismatches + QOPEN*(qopens + qindels) + TOPEN*(topens + tindels);
      debug(printf("%d matches, %d nmismatches, %d+%d qgaps, %d+%d tgaps => intronscore is %d\n",
		   nmatches,nmismatches,qopens,qindels,topens,tindels,
		   MISMATCH*nmismatches + QOPEN*(qopens + qindels) + TOPEN*(topens + tindels)));

      pairs = Pairpool_transfer(pairs,peeled_pairs);
      path = Pairpool_transfer(path,peeled_path);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->genomejump > pair->queryjump + SINGLESLEN) {
      /* Intron length shorter than MININTRONLEN_FINAL.  Just replace the gap */
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else {
      /* Single gap; force fill */
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
    }
  }

  return pairs;
}  


#if 0
static List_T
trim_pairs (List_T pairs, int trim5pair, int trim3pair, int matchespergap, Pairpool_T pairpool) {
  List_T path = NULL;
  Pair_T firstpair, pair;
  int i = 0, j;
  bool nongapp = false;		/* To remove initial gap */

  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    if (!pair->gapp) {
      if (i >= trim5pair && i <= trim3pair) {
	path = Pairpool_push_existing(path,pairpool,pair);
	nongapp = true;
      }
      i++;
    } else {
      for (j = 0; j < matchespergap-1; j++) {
	i++;
      }
      if (i >= trim5pair && i <= trim3pair) {
	if (nongapp == false) {
	  /* Don't push an initial gap */
	} else {
	  path = Pairpool_push_existing(path,pairpool,pair);
	}
      }
      i++;
    }
  }

  path = Pairpool_pop(path,&firstpair);
  while (path != NULL && firstpair->gapp == true) {
    path = Pairpool_pop(path,&firstpair);
  }
  pairs = Pairpool_push_existing(NULL,pairpool,firstpair);
  pairs = Pairpool_transfer(pairs,path);

  return pairs;
}
#endif

static List_T
path_compute (int *totalintronscore, double *defect_rate, int *intronlen, int *nonintronlen, 
	      List_T path, int cdna_direction, int introndir, int querylength, int skiplength, 
	      Genomicpos_T genomicstart, Genomicpos_T genomiclength,
#ifdef PMAP
	      char *queryaaseq_ptr,
#endif
	      char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
	      Genome_T genome, Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrpos, Genomicpos_T chrlength, 
	      int maxpeelback, int nullgap, int extramaterial_end, int extramaterial_paired,
	      int extraband_single, int extraband_end, int extraband_paired, int minendexon,
	      Pairpool_T pairpool, Gbuffer_T gbuffer, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
	      bool end_microexons_p, Stage3debug_T stage3debug, bool diagnosticp, bool do_final_p,
	      Oligoindex_T *oligoindices, int noligoindices, Diagpool_T diagpool, int sufflookback, int nsufflookback,
	      int maxintronlen_bound) {
  List_T pairs = NULL;
  int iter0, iter1, iter2;
  int nintrons = 0, nnonintrons = 0;
  int npairs, ngaps, ncanonical;
  bool shiftp, incompletep, singlep, shortp, deletep, trim5p, trim3p;
  int dynprogindex_minor = DYNPROGINDEX_MINOR, dynprogindex_major = DYNPROGINDEX_MAJOR;
  int matches, unknowns, mismatches, qopens, qindels, topens, tindels,
    nsemicanonical, nnoncanonical;
  int maxpeelback5, maxpeelback3;
  bool chop_exon_p, trim_middle_p = true;
#ifdef DEBUG
  Pair_T firstpair, lastpair;
#endif

  *totalintronscore = 0;
  if (path == NULL) {
    return NULL;
  }

  if (stage3debug == POST_STAGE2) {
    pairs = List_reverse(path);
    path = insert_gapholders(pairs,pairpool);
    pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
    return pairs;
  }

  iter0 = 0;
  while (trim_middle_p == true && iter0 < MAXITER) {
    /* Pass 0: Insert gaps.  pairs --> pairs */
    debug(printf("\n*** Pass 0 (dir %d): Solve single nucleotide gaps\n",cdna_direction));

    pairs = List_reverse(path);
    path = insert_gapholders(pairs,pairpool);
    pairs = List_reverse(path);

    /* Pass 1: Initial smoothing.  pairs --> path */
    debug(printf("\n*** Pass 1 (dir %d): Initial smoothing by net gap\n",cdna_direction));
    pairs = Smooth_pairs_by_netgap(&deletep,pairs,pairpool);
    if (deletep == true) {
      path = insert_gapholders(pairs,pairpool);
    } else {
      path = List_reverse(pairs);
    }

    if (stage3debug == POST_SMOOTHING) {
      pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
      return pairs;
    }

#ifdef PMAP
    /* Pass 1b: undefine nucleotides around gaps.  path --> path */
    pairs = undefine_nucleotides(queryseq_ptr,querylength,path,pairpool,/*width*/6);
    path = List_reverse(pairs);
#endif

    /* Pass 2: solve straight gaps.  path --> pairs (for defect rate) */
    debug(printf("\n*** Pass 2 (dir %d): Solve straight gaps\n",cdna_direction));
    pairs = build_pairs_singles(&dynprogindex_minor,path,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				maxpeelback,nullgap,extramaterial_paired,extraband_single,
				/*defect_rate*/0.0,/*widebandp*/false,pairpool,dynprogL,dynprogM,dynprogR);

    if (stage3debug == POST_SINGLES) {
      path = insert_gapholders(pairs,pairpool);
      pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
      return pairs;
    }

    /* Compute defect rate here */
    Pair_fracidentity(&matches,&unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
		      &ncanonical,&nsemicanonical,&nnoncanonical,pairs,/*cdna_direction*/0);
    *defect_rate = (double) mismatches/(double) (matches + mismatches);
    debug(printf("defect_rate = %f (%d matches, %d mismatches)\n",*defect_rate,matches,mismatches));
    

    /* Pass 3: introns.  pairs --> pairs */
    debug(printf("\n*** Pass 3: Smooth and solve dual introns iteratively\n"));
    iter1 = 0;
    shortp = true;
    while (shortp == true && iter1++ < MAXITER) {
      /* Pass 3a: smoothing.  pairs --> pairs */
      debug(printf("*** Pass 3a: Smoothing by size.  Iteration %d\n",iter1));
      pairs = Smooth_pairs_by_size(&shortp,&deletep,pairs,pairpool,/*stage2_indexsize*/6);
      
      /* Pass 3b: dual introns.  pairs --> pairs */
      if (deletep == false && shortp == false) {
	/* Do nothing */
      } else if (deletep == true && shortp == true) {
	path = insert_gapholders(pairs,pairpool);
	debug(printf("*** Pass 3b: Solve dual introns.  Iteration %d\n",iter1));
	pairs = build_pairs_dualintrons(&singlep,&dynprogindex_major,path,
#ifdef PMAP
					queryaaseq_ptr,
#endif
					queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
					maxpeelback,nullgap,extramaterial_paired,extraband_paired,
					*defect_rate,pairpool,dynprogL,dynprogM,dynprogR);
      } else if (deletep == true) {
	path = insert_gapholders(pairs,pairpool);
	pairs = List_reverse(path);
      } else {
	path = List_reverse(pairs);
	debug(printf("*** Pass 3b: Solve dual introns.  Iteration %d\n",iter1));
	pairs = build_pairs_dualintrons(&singlep,&dynprogindex_major,path,
#ifdef PMAP
					queryaaseq_ptr,
#endif
					queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
					maxpeelback,nullgap,extramaterial_paired,extraband_paired,
					*defect_rate,pairpool,dynprogL,dynprogM,dynprogR);
      }

      /* Pass 3c: single introns.  pairs --> pairs */
      iter2 = 0;
      shiftp = incompletep = true;
      debug(printf("\n*** Pass 3c: Solve introns iteratively\n"));
      while ((shiftp == true || incompletep == true) && iter2++ < MAXITER) {
	debug(printf("*** Pass 3c: Solve introns.  Iteration %d\n",iter2));
	path = insert_gapholders(pairs,pairpool);
	pairs = build_pairs_introns(&shiftp,&incompletep,
				    &nintrons,&nnonintrons,&(*intronlen),&(*nonintronlen),
				    &dynprogindex_minor,&dynprogindex_major,path,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				    maxpeelback,nullgap,extramaterial_paired,extraband_single,extraband_paired,
				    *defect_rate,pairpool,dynprogL,dynprogM,dynprogR,/*finalp*/false,
				    /*forcep*/iter2 == MAXITER ? true : false);
      }
    }

    if (stage3debug == POST_INTRONS) {
      path = insert_gapholders(pairs,pairpool);
      pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
      return pairs;
    }

    /* Pass 4: Final pass to solve for introns with higher rewards for canonical introns: pairs --> pairs */
    debug(printf("\n*** Pass 4: Final pass to find canonical introns\n"));
    if (do_final_p == true) {
      path = insert_gapholders(pairs,pairpool);
      pairs = build_pairs_introns(&shiftp,&incompletep,
				  &nintrons,&nnonintrons,&(*intronlen),&(*nonintronlen),
				  &dynprogindex_minor,&dynprogindex_major,path,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				  maxpeelback,nullgap,extramaterial_paired,extraband_single,extraband_paired,
				  *defect_rate,pairpool,dynprogL,dynprogM,dynprogR,/*finalp*/true,
				  /*forcep*/true);
    }

    if (stage3debug == POST_CANONICAL) {
      path = insert_gapholders(pairs,pairpool);
      pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
      return pairs;
    }


#if 1
    /* Could potentially be replaced by distal/medial computation */
    /* Pass 5: Trim distal ends. */
    path = insert_gapholders(pairs,pairpool);
    pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);

    /* Pass 5a: trim end exons: pairs -> pairs */
    debug(printf("\n*** Iteration %d, Pass 5a: Trim end exons\n",++iter1));
    pairs = chop_ends_by_changepoint(pairs,pairpool);
      
    if (stage3debug == POST_CHANGEPOINT) {
      /* Gapholders should be present */
      return pairs;
    }

    /* Pass 5b: fix ends: pairs -> pairs */
    /* Pass 5b: solve 5' end */
    debug(printf("\n*** Iteration %d, Pass 5b: Solve and clean 5' end\n",iter1));
    /* Necessary to remove gap at end, e.g., AA012859 */
    pairs = clean_pairs_end5(pairs);
    
    pairs = build_pairs_end5(&chop_exon_p,&dynprogindex_minor,&dynprogindex_major,pairs,
#ifdef PMAP
			     queryaaseq_ptr,
#endif
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			     cdna_direction,/*newintrondir*/cdna_direction,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_end,
			     *defect_rate,pairpool,dynprogR,end_microexons_p,/*extendp*/true);
    /* Necessary, e.g., AA026523 */
    pairs = clean_pairs_end5(pairs);
    path = List_reverse(pairs);
    
    debug(
	  if (pairs) {
	    firstpair = pairs->first;
	    printf("5' starts at %d\n",firstpair->querypos);
	  });
  
    /* Pass 5b: solve 3' end */
    debug(printf("\n*** Iteration %d, Pass 5b: Solve and clean 3' end\n",iter1));
    /* Necessary to remove gaps at end */
    path = clean_path_end3(path);
  
    path = build_path_end3(&chop_exon_p,&dynprogindex_minor,&dynprogindex_major,path,
			   querylength,genomiclength,
#ifdef PMAP
			   queryaaseq_ptr,
#endif
			   queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			   cdna_direction,/*newintrondir*/cdna_direction,maxpeelback,nullgap,
			   extramaterial_end,extramaterial_paired,extraband_end,
			   *defect_rate,pairpool,dynprogL,end_microexons_p,/*extendp*/true);
    /* Necessary, e.g., AA026523 */
    path = clean_path_end3(path);
    pairs = List_reverse(path);
  
    debug(
	  if (path) {
	    lastpair = path->first;
	    printf("3' ends at %d\n",lastpair->querypos);
	  });
#endif
  

    /* Pass 6: Compare distal exons with medial extensions */
    /* Need to assign gap types because traversal of ends depends on knowing if end introns are canonical */
    path = insert_gapholders(pairs,pairpool);
    pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);

#if 0
    /* If don't limit peelback, then algorithm can be very slow */
    maxpeelback5 = maxpeelback3 = npairs;
#else
    maxpeelback5 = maxpeelback3 = maxpeelback;
#endif

    /* Pass 6a: solve 5' end trying both distal and medial: pairs --> pairs */
    chop_exon_p = true;
    while (pairs && chop_exon_p == true) {
      debug(printf("\n*** Pass 6a: Solve 5' end\n"));
      pairs = build_pairs_end5(&chop_exon_p,&dynprogindex_minor,&dynprogindex_major,pairs,
#ifdef PMAP
			       queryaaseq_ptr,
#endif
			       queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			       cdna_direction,/*newintrondir*/cdna_direction,maxpeelback5,nullgap,
			       extramaterial_end,extramaterial_paired,extraband_end,
			       *defect_rate,pairpool,dynprogR,end_microexons_p,
			       /*extendp*/false);
      debug(
	    if (pairs) {
	      firstpair = pairs->first;
	      printf("5' starts at %d\n",firstpair->querypos);
	    });
    }

    /* Pass 6b: solve 3' end trying both distal and medial: pairs --> pairs */
    path = List_reverse(pairs);
    chop_exon_p = true;
    while (path && chop_exon_p == true) {
      debug(printf("\n*** Pass 6b: Solve 3' end\n"));
      path = build_path_end3(&chop_exon_p,&dynprogindex_minor,&dynprogindex_major,path,
			     querylength,genomiclength,
#ifdef PMAP
			     queryaaseq_ptr,
#endif
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			     cdna_direction,/*newintrondir*/cdna_direction,maxpeelback3,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_end,
			     *defect_rate,pairpool,dynprogL,end_microexons_p,
			     /*extendp*/false);
      debug(
	    if (path) {
	      lastpair = path->first;
	      printf("3' ends at %d\n",lastpair->querypos);
	    });
    }

    if (stage3debug == POST_DISTAL_MEDIAL) {
      pairs = List_reverse(path);
      path = insert_gapholders(pairs,pairpool);
      pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
      return pairs;
    }

    pairs = List_reverse(path);

#if 1
    /* Pass 6c: trim short exons from ends: pairs -> pairs */
    debug(printf("\n*** Pass 6c: Trim short end exons\n"));
    pairs = trim_short_end_exons(&trim5p,&trim3p,pairs,pairpool,minendexon);
#endif
    

    /* Pass 7: Trim bad middle exons */
    debug(printf("\n*** Pass 7: Trim bad middle exons\n"));

#if 0
    trim_middle_p = false;
#else
    path = insert_gapholders(pairs,pairpool);
    pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
    
    path = trim_exons_middle(&trim_middle_p,pairs,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
    if (trim_middle_p == false) {
      pairs = List_reverse(path);

    } else {
      pairs = build_dual_breaks(path,&dynprogindex_minor,&dynprogindex_major,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				cdna_direction,pairpool,dynprogL,dynprogR,maxpeelback,
				oligoindices,noligoindices,diagpool,sufflookback,nsufflookback,
				maxintronlen_bound);
    }
#endif

    if (stage3debug == POST_TRIM_MIDDLE) {
      path = insert_gapholders(pairs,pairpool);
      pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
      return pairs;
    }

    iter0++;
  }


  /* End */
  path = insert_gapholders(pairs,pairpool);
  pairs = assign_gap_types(&npairs,&ngaps,&ncanonical,path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);
  path = List_reverse(pairs);
  pairs = score_introns(&(*totalintronscore),path,cdna_direction,/*maxpeelback*/MAXPEELBACK_SCORE,
			nullgap,pairpool);

  debug(Pair_dump_list(pairs,true));
  debug(return pairs);

  if (diagnosticp == true) {
    if (pairs != NULL) {
      path = check_gaps(pairs,pairpool);
      pairs = List_reverse(path);
    }
  }

  return pairs;
}


/* Introndir shows whether we are certain about the intron direction.
   However, we try both directions anyway, using cdna_direction */
T
Stage3_compute (List_T path, Genomicpos_T genomicstart, Genomicpos_T genomiclength,
#ifdef PMAP
		Sequence_T queryaaseq,
#endif
		Sequence_T queryseq, Sequence_T queryuc,
		Sequence_T genomicseg, Sequence_T genomicuc,
		int straintype, char *strain,
		Genome_T genome, Genome_T genomealt, Chrnum_T chrnum,
		Genomicpos_T chroffset, Genomicpos_T chrpos, Genomicpos_T chrlength,
		bool watsonp, int maxpeelback, int nullgap, 
		int extramaterial_end, int extramaterial_paired,
		int extraband_single, int extraband_end, int extraband_paired,
		bool end_microexons_p, int minendexon,
		Pairpool_T pairpool, Gbuffer_T gbuffer, Dynprog_T dynprogL,
		Dynprog_T dynprogM, Dynprog_T dynprogR,
		IIT_T altstrain_iit, int ngap, Stage3debug_T stage3debug,
		bool diagnosticp, bool checkp, Stopwatch_T stopwatch, bool do_final_p, int sense,
		Oligoindex_T *oligoindices, int noligoindices, Diagpool_T diagpool,
		int sufflookback, int nsufflookback, int maxintronlen) {
  T this;
  List_T pairs_fwd, pairs_rev, path_fwd, path_rev;
  double defect_rate;
  int querylength;
  int introndir;
  int fwd_intronlen = 0, rev_intronlen = 0;
  int fwd_nonintronlen = 0, rev_nonintronlen = 0;
  int totalintronscore_fwd = 0, totalintronscore_rev = 0;
  Pair_T start, end;
#ifdef PMAP
  char *queryaaseq_ptr;
#endif
  char *queryseq_ptr, *queryuc_ptr, *genomicseg_ptr, *genomicuc_ptr;
  int *typematches, nmatches;

  Stopwatch_start(stopwatch);

  debug(printf("Stage 3: *** Starting stage 3 for genomiclength %u at %u (chrnum %d:%d)\n",
	       Sequence_fulllength(genomicseg),genomicstart,chrnum,chrpos));

  querylength = Sequence_fulllength(queryseq);
#ifdef PMAP
  queryaaseq_ptr = Sequence_fullpointer(queryaaseq);
#endif
  queryseq_ptr = Sequence_fullpointer(queryseq);
  queryuc_ptr = Sequence_fullpointer(queryuc);
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);

#ifdef PMAP
  path_fwd = path;
  path_rev = (List_T) NULL;
  pairs_rev = (List_T) NULL;
  introndir = +1;
  /* do_final_p = true; */
#else
  if (sense == 0) {
    /* Should try both even if no introns (cf, AA011563) */
    path_fwd = path;
    path_rev = Pairpool_copy(path,pairpool);
    introndir = 0;
  } else if (sense > 0) {
    path_fwd = path;
    path_rev = (List_T) NULL;
    introndir = +1;
  } else if (sense < 0) {
    path_fwd = (List_T) NULL;
    path_rev = path;
    introndir = -1;
  }
#endif

  if (path_fwd != NULL) {
    pairs_fwd = path_compute(&totalintronscore_fwd,&defect_rate,&fwd_intronlen,&fwd_nonintronlen,path_fwd,
			     /*cdna_direction*/+1,introndir,
			     querylength,Sequence_skiplength(queryseq),
			     genomicstart,genomiclength,
#ifdef PMAP
			     queryaaseq_ptr,
#endif
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			     genome,chrnum,chroffset,chrpos,chrlength,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_single,extraband_end,extraband_paired,
			     minendexon,pairpool,gbuffer,dynprogL,dynprogM,dynprogR,
			     end_microexons_p,stage3debug,diagnosticp,do_final_p,
			     oligoindices,noligoindices,diagpool,sufflookback,nsufflookback,maxintronlen);
  } else {
    pairs_fwd = NULL;
  }

#ifndef PMAP
  if (path_rev != NULL) {
    pairs_rev = path_compute(&totalintronscore_rev,&defect_rate,&rev_intronlen,&rev_nonintronlen,path_rev,
			     /*cdna_direction*/-1,introndir,
			     querylength,Sequence_skiplength(queryseq),
			     genomicstart,genomiclength,
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			     genome,chrnum,chroffset,chrpos,chrlength,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_single,extraband_end,extraband_paired,
			     minendexon,pairpool,gbuffer,dynprogL,dynprogM,dynprogR,
			     end_microexons_p,stage3debug,diagnosticp,do_final_p,
			     oligoindices,noligoindices,diagpool,sufflookback,nsufflookback,maxintronlen);
  } else {
    pairs_rev = NULL;
  }
#endif

  debug(printf("Totalintronscores: %d fwd, %d rev\n",totalintronscore_fwd,totalintronscore_rev));
  if ((this = Stage3_new(pairs_fwd,pairs_rev,totalintronscore_fwd,totalintronscore_rev,
			 genomicstart,genomiclength,defect_rate)) == NULL) {
    Stopwatch_stop(stopwatch);
    return NULL;
  } else {
    debug7(printf("Preparing for print at chrnum %d\n",chrnum));
    this->pairarray = prepare_for_printing(&this->npairs,&this->pairs,this->cdna_direction,pairpool,queryseq_ptr,
#ifdef PMAP
					   queryaaseq_ptr,
#endif
					   genomicseg_ptr,genomicuc_ptr,ngap,
					   Sequence_subseq_offset(queryseq),Sequence_skiplength(queryseq),
					   diagnosticp);
  }

  if (checkp == true && stage3debug == NO_STAGE3DEBUG && 
      Pair_check_array(this->pairarray,this->npairs) == true) {
    Pair_dump_array(this->pairarray,this->npairs,/*zerobasedp*/true);
#ifndef DEBUG
    Except_raise(&coordinate_error,__FILE__,__LINE__);
#endif
  }
#ifdef DEBUG
  Pair_dump_array(this->pairarray,this->npairs,/*zerobasedp*/true);
#endif

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  this->straintype = straintype;
  this->strain = strain;
  this->chrnum = chrnum;
  this->chrpos = chrpos;
  this->chroffset = chroffset;
  this->genomiclength = genomiclength;
  this->watsonp = watsonp;
  this->defect_rate = defect_rate;

  if (watsonp) {
    this->genomicstart = chroffset + chrpos + Pair_genomepos(start);
    this->genomicend = chroffset + chrpos + Pair_genomepos(end);
  } else {
    this->genomicstart = chroffset + chrpos + (genomiclength - 1) - Pair_genomepos(start);
    debug(printf("genomicstart: %u + %u + (%u - 1) - %u = %u\n",
		 chroffset,chrpos,genomiclength,Pair_genomepos(start),this->genomicstart));
    this->genomicend = chroffset + chrpos + (genomiclength - 1) - Pair_genomepos(end);
  }

  this->stage3_runtime = Stopwatch_stop(stopwatch);

  this->trimmed_coverage = (double) (end->querypos - start->querypos + 1)/(double) (Sequence_trimlength(queryseq) + Sequence_skiplength(queryseq));

  if (straintype == 0) {
    return this;
  } else {
    if (watsonp) {
      typematches = IIT_get_typed(&nmatches,altstrain_iit,/*divstring*/NULL,
				  this->genomicstart,this->genomicend,straintype,/*sortp*/false);
    } else {
      typematches = IIT_get_typed(&nmatches,altstrain_iit,/*divstring*/NULL,
				  this->genomicend,this->genomicstart,straintype,/*sortp*/false);
    }
    if (typematches == NULL) {
      Stage3_free(&this);
      return NULL;
    } else {
      FREE(typematches);
      return this;
    }
  }
}


static void
adjust_genomepos_pairs (List_T pairs, int delta, bool watsonp) {
  List_T p;
  Pair_T pair;

  for (p = pairs; p != NULL; p = p->rest) {
    pair = (Pair_T) List_head(p);
    if (watsonp == true) {
      pair->genomepos = delta + pair->genomepos;
    } else {
      pair->genomepos = delta - pair->genomepos;
    }
  }
  return;
}


static List_T
make_path_direct (Gregion_T gregion, char *queryseq_ptr, Pairpool_T pairpool, Genome_T genome, 
		  Gbuffer_T gbuffer, Genomicpos_T chroffset, Genomicpos_T chrpos, bool plusp) {
  List_T path = NULL;
  Pair_T pair;
  char *genomicseg_ptr;
  Sequence_T genomicseg;
  int querypos1, querypos2, queryjump, genomejump, i, j;
  Genomicpos_T position1, position2, genomepos;
  int matchsize = Gregion_matchsize(gregion);

  querypos1 = Gregion_querystart(gregion);
  querypos2 = Gregion_queryend(gregion);
  
  position1 = Gregion_genomicstart(gregion);
  position2 = position1 + Gregion_genomiclength(gregion);

  /* Handle the 5' end */
  if (plusp == true) {
    genomicseg = Genome_get_segment(genome,position1,matchsize,/*chromosome_iit*/NULL,/*revcomp*/false,
				    Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
  } else {
    genomicseg = Genome_get_segment(genome,position1-(matchsize-1U),matchsize,/*chromosome_iit*/NULL,/*revcomp*/true,
				    Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);

  genomepos = position1 - chroffset - chrpos;
  for (i = querypos1, j = 0; i < querypos1 + matchsize; i++, j++) {
    /* printf("%c %c\n",queryseq_ptr[i],genomicseg_ptr[j]); */
    path = Pairpool_push(path,pairpool,i,genomepos,queryseq_ptr[i],MATCH_COMP,genomicseg_ptr[j],/*dynprogindex*/0);
    if (plusp == true) {
      genomepos += 1U;
    } else {
      genomepos -= 1U;
    }
  }
  Sequence_free(&genomicseg);

  /* Handle the gap */
  queryjump = querypos2 - querypos1;
  genomejump = position2 - position1;
  path = Pairpool_push_gapholder(path,pairpool,queryjump,genomejump);
  pair = (Pair_T) path->first;
  pair->comp = DUALBREAK_COMP;
  /* printf("\n"); */

  /* Handle the 3' end */
  if (plusp == true) {
    genomicseg = Genome_get_segment(genome,position2,matchsize,/*chromosome_iit*/NULL,/*revcomp*/false,
				    Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
  } else {
    genomicseg = Genome_get_segment(genome,position2-(matchsize-1U),matchsize,/*chromosome_iit*/NULL,/*revcomp*/true,
				    Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);

  genomepos = position2 - chroffset - chrpos;
  for (i = querypos2, j = 0; i < querypos2 + matchsize; i++, j++) {
    /* printf("%c %c\n",queryseq_ptr[i],genomicseg_ptr[j]); */
    path = Pairpool_push(path,pairpool,i,genomepos,queryseq_ptr[i],MATCH_COMP,genomicseg_ptr[j],/*dynprogindex*/0);
    if (plusp == true) {
      genomepos += 1U;
    } else {
      genomepos -= 1U;
    }
  }
  Sequence_free(&genomicseg);

  return path;
}


T
Stage3_direct (Gregion_T gregion,
#ifdef PMAP
	       Sequence_T queryaaseq,
#endif
	       Sequence_T queryseq, Sequence_T queryuc, Pairpool_T pairpool, Genome_T genome,
	       Chrnum_T chrnum,  Genomicpos_T chroffset, Genomicpos_T chrpos, bool watsonp,
	       int ngap, Gbuffer_T gbuffer, Dynprog_T dynprogL, Dynprog_T dynprogR,
	       int extramaterial_end, int extraband_end) {
  T new;
  List_T path, pairs_fwd, gappairs;
  Pair_T start, end, leftpair, rightpair, pair;
#ifdef PMAP
  char *queryaaseq_ptr;
#endif
  Sequence_T genomicseg, genomicuc;
  Genomicpos_T genomicpos;
  char *queryseq_ptr, *queryuc_ptr, *genomicseg_ptr, *genomicuc_ptr;
  int i, querylength, queryjump, genomejump, querydp5, querydp3, genomedp5, genomedp3, endadj;
  int finalscore, nmatches, nmismatches, nopens, nindels;
  int dynprogindex_minor = +1;

  querylength = Sequence_fulllength(queryseq);
  queryseq_ptr = Sequence_fullpointer(queryseq);
#ifdef PMAP
  queryaaseq_ptr = Sequence_fullpointer(queryaaseq);
#endif
  queryuc_ptr = Sequence_fullpointer(queryuc);
  path = make_path_direct(gregion,queryseq_ptr,pairpool,genome,gbuffer,chroffset,chrpos,watsonp);

  leftpair = (Pair_T) List_head(path);
  querydp5 = leftpair->querypos + 1;
  genomicpos = chroffset + chrpos + leftpair->genomepos;
  queryjump = querylength - leftpair->querypos - 1;
  if (leftpair->cdna == ' ') queryjump++;

  genomejump = queryjump + extramaterial_end; /* proposed */
  debug(printf("For 3' end, genomicpos is %u, queryjump = %d, genomejump = %d\n",genomicpos,queryjump,genomejump));

  if (watsonp == true) {
    genomicseg = Genome_get_segment(genome,genomicpos+1U,genomejump,/*chromosome_iit*/NULL,/*revcomp*/false,
				    Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
    genomedp5 = leftpair->genomepos + 1U;
  } else {
    genomicseg = Genome_get_segment(genome,genomicpos-genomejump,genomejump,/*chromosome_iit*/NULL,/*revcomp*/true,
				    Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
    genomedp5 = leftpair->genomepos - 1U;
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc = Sequence_uppercase(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);

  gappairs = Dynprog_end3_gap(&dynprogindex_minor,&finalscore,&nmatches,&nmismatches,&nopens,&nindels,dynprogL,
			      /*sequence1*/&(queryseq_ptr[querydp5]),/*sequenceuc1*/&(queryuc_ptr[querydp5]),
			      /*sequence2*/genomicseg_ptr,/*sequenceuc2*/genomicuc_ptr,
			      /*length1*/queryjump,/*length2*/genomejump,/*offset1*/querydp5,/*offset2*/0,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      /*cdna_direction*/+1,pairpool,extraband_end,/*defect_rate*/0.0);
  Sequence_free(&genomicuc);
  Sequence_free(&genomicseg);

  debug(printf("Adjusting 3' pairs by %d\n",genomedp5));
  adjust_genomepos_pairs(gappairs,genomedp5,watsonp);

  path = Pairpool_transfer(path,List_reverse(gappairs));
  pair = path->first;
  endadj = pair->genomepos;
  debug(printf("Top of path has genomepos %d\n",pair->genomepos));

  pairs_fwd = List_reverse(path);

  rightpair = (Pair_T) List_head(pairs_fwd);
  querydp3 = rightpair->querypos - 1;
  genomicpos = chroffset + chrpos + rightpair->genomepos;
  queryjump = rightpair->querypos;
  genomejump = queryjump + extramaterial_end; /* proposed */
  debug(printf("For 5' end, genomicpos is %u, queryjump = %d, genomejump = %d\n",genomicpos,queryjump,genomejump));

  if (watsonp == true) {
    genomicseg = Genome_get_segment(genome,genomicpos-genomejump,genomejump,/*chromosome_iit*/NULL,/*revcomp*/false,
				    Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
    genomedp3 = rightpair->genomepos - 1U;
  } else {
    genomicseg = Genome_get_segment(genome,genomicpos+1U,genomejump,/*chromosome_iit*/NULL,/*revcomp*/true,
				    Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
    genomedp3 = rightpair->genomepos + 1U;
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc = Sequence_uppercase(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);

  gappairs = Dynprog_end5_gap(&dynprogindex_minor,&finalscore,&nmatches,&nmismatches,&nopens,&nindels,dynprogR,
			      &(queryseq_ptr[querydp3]),&(queryuc_ptr[querydp3]),
			      &(genomicseg_ptr[genomejump-1]),&(genomicuc_ptr[genomejump-1]),
			      queryjump,genomejump,/*offset1*/querydp3,/*offset2*/0,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      /*cdna_direction*/+1,pairpool,extraband_end,/*defect_rate*/0.0);
  Sequence_free(&genomicuc);
  Sequence_free(&genomicseg);

  debug(printf("Adjusting 5' pairs by %d\n",genomedp3));
  adjust_genomepos_pairs(gappairs,genomedp3,watsonp);
  pairs_fwd = Pairpool_transfer(pairs_fwd,gappairs);

  if (watsonp == false) {
    start = pairs_fwd->first;
    debug(printf("Adjusting all pairs by 0 (not by %d + %d)\n",start->genomepos,endadj));
    adjust_genomepos_pairs(pairs_fwd,/*start->genomepos*/0,/*watsonp*/false);
  }

  new = Stage3_new(pairs_fwd,/*pairs_rev*/NULL,/*totalintronscore_fwd*/0,/*totalintronscore_rev*/-100,
		   /*stage1_genomicstart*/0U,/*stage1_genomiclength*/0U,/*defect_rate*/0.0);
  new->cdna_direction = 0;	/* Override the +1 assigned by Stage3_new */

  new->pairarray = prepare_for_printing(&new->npairs,&new->pairs,new->cdna_direction,pairpool,queryseq_ptr,
#ifdef PMAP
					queryaaseq_ptr,
#endif
					/*genomicseg_ptr*/NULL,/*genomicuc_ptr*/NULL,ngap,
					/*subseq_offset*/0,Sequence_skiplength(queryseq),/*diagnosticp*/true);
  for (i = 0; i < new->npairs; i++) {
    pair = &(new->pairarray[i]);
    pair->aapos = 0;
    pair->aa_g = ' ';
    pair->aa_e = ' ';
    pair->aaphase_g = -1;
    pair->aaphase_e = -1;
  }

  start = &(new->pairarray[0]);
  end = &(new->pairarray[new->npairs-1]);

  new->straintype = 0;
  new->strain = (char *) NULL;
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrpos = chrpos;
  new->genomiclength = Pair_genomepos(end)+endadj+1;
  new->watsonp = watsonp;
  new->defect_rate = 0.0;

  debug(printf("chroffset = %u, chrpos = %u, genomiclength = %u\n",chroffset,chrpos,new->genomiclength));
  if (watsonp) {
    new->genomicstart = chroffset + chrpos + Pair_genomepos(start);
    new->genomicend = chroffset + chrpos + Pair_genomepos(end);
  } else {
    new->genomicstart = chroffset + chrpos + (new->genomiclength - 1) - Pair_genomepos(start);
    new->genomicend = chroffset + chrpos + (new->genomiclength - 1) - Pair_genomepos(end);
  }

  return new;
}


/************************************************************************
 *  Merging
 ************************************************************************/

bool
Stage3_mergeable (char *comp, Stage3_T firstpart, Stage3_T secondpart, int chimera_direction, 
		  double donor_prob, double acceptor_prob) {
  Pair_T end1, start2;
  bool watsonp, connectablep = false;
  Genomicpos_T endchrpos1, startchrpos2;
  int querygap;

  if (firstpart->chrnum != secondpart->chrnum) {
    return false;
  } else if (firstpart->watsonp != secondpart->watsonp) {
    return false;
  } else {
    watsonp = firstpart->watsonp;
    end1 = &(firstpart->pairarray[firstpart->npairs-1]);
    start2 = &(secondpart->pairarray[0]);
    if (watsonp == true) {
      endchrpos1 = firstpart->chrpos+end1->genomepos;
      startchrpos2 = secondpart->chrpos+start2->genomepos;
      if (endchrpos1 < startchrpos2) {
	if (startchrpos2 < endchrpos1 + MERGELENGTH) {

	  connectablep = true;
	} else if (donor_prob > DONOR_THRESHOLD && acceptor_prob >= ACCEPTOR_THRESHOLD &&
		   startchrpos2 < endchrpos1 + LONG_MERGELENGTH) {
	  connectablep = true;
	}
      }
    } else {
      endchrpos1 = firstpart->chrpos + (firstpart->genomiclength - 1) - end1->genomepos;
      startchrpos2 = secondpart->chrpos + (secondpart->genomiclength - 1) - start2->genomepos;
      if (startchrpos2 < endchrpos1) {
	if (endchrpos1 < startchrpos2 + MERGELENGTH) {
	  connectablep = true;
	} else if (donor_prob > DONOR_THRESHOLD && acceptor_prob >= ACCEPTOR_THRESHOLD &&
		   endchrpos1 < startchrpos2 + LONG_MERGELENGTH) {
	  connectablep = true;
	}
      }
    }

    if (connectablep == false) {
      return false;
    } else {
      querygap = start2->querypos - end1->querypos - 1;
      if (querygap > 0) {
	*comp = DUALBREAK_COMP;
      } else if (chimera_direction > 0) {
	*comp = FWD_CANONICAL_INTRON_COMP;
      } else if (chimera_direction < 0) {
	*comp = REV_CANONICAL_INTRON_COMP;
      } else {
	*comp = NONINTRON_COMP;
      }

      return true;
    }
  }
}


static void
adjust_genomepos (Stage3_T stage3, int delta) {
  int i;
  Pair_T pair;

  for (i = 0; i < stage3->npairs; i++) {
    pair = &(stage3->pairarray[i]);
    pair->genomepos += delta;
  }
  return;
}


void
Stage3_merge (T firstpart, T secondpart, char comp, Pairpool_T pairpool, Genome_T genome, int ngap) {
  struct Pair_T *oldpairs;
  Pair_T end1, start2, pair;
  List_T pairs = NULL;
  bool watsonp;
  int genomicpos, querypos, ptr, i;
  char *gbuffer1, *gbuffer2, *genomicseg_ptr;
  Genomicpos_T left;
  Sequence_T genomicseg;

  watsonp = firstpart->watsonp;
  end1 = &(firstpart->pairarray[firstpart->npairs-1]);
  start2 = &(secondpart->pairarray[0]);

  if (watsonp == true) {
    adjust_genomepos(secondpart,secondpart->chrpos - firstpart->chrpos);
  } else {
    adjust_genomepos(secondpart,firstpart->chrpos + firstpart->genomiclength 
		     - secondpart->chrpos - secondpart->genomiclength);
  }

  oldpairs = firstpart->pairarray;
  firstpart->pairarray = (struct Pair_T *) CALLOC(firstpart->npairs+ngap+3+ngap+secondpart->npairs,sizeof(struct Pair_T));
  memcpy(firstpart->pairarray,oldpairs,firstpart->npairs*sizeof(struct Pair_T));
  ptr = firstpart->npairs;

  genomicpos = start2->genomepos;
  querypos = start2->querypos;
  if (comp == DUALBREAK_COMP) {
    pairs = Pairpool_push_gapalign(pairs,pairpool,querypos,genomicpos,' ',comp,' ',/*extraexonp*/false);
    pair = (Pair_T) pairs->first;
    for (i = 0; i < ngap; i++) {
      memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    }
    pairs = Pairpool_push_gapalign(pairs,pairpool,querypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,/*extraexonp*/false);
    pair = (Pair_T) pairs->first;
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    pairs = Pairpool_push_gapalign(pairs,pairpool,querypos,genomicpos,' ',comp,' ',/*extraexonp*/false);
    pair = (Pair_T) pairs->first;
    for (i = 0; i < ngap; i++) {
      memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    }

  } else {
    if (watsonp == true) {
      left = firstpart->chroffset + firstpart->chrpos + end1->genomepos + 1;
    } else {
      left = firstpart->chroffset + firstpart->chrpos + (firstpart->genomiclength - 1) - end1->genomepos - ngap;
    }
    
    gbuffer1 = (char *) CALLOC(ngap+1,sizeof(char));
    gbuffer2 = (char *) CALLOC(ngap+1,sizeof(char));
    genomicseg = Genome_get_segment(genome,left,ngap,/*chromosome_iit*/NULL,!watsonp,gbuffer1,gbuffer2,ngap);
    genomicseg_ptr = Sequence_fullpointer(genomicseg);

    for (i = 0; i < ngap; i++) {
      pairs = Pairpool_push_gapalign(pairs,pairpool,querypos,genomicpos,' ',comp,genomicseg_ptr[i],/*extraexonp*/false);
      pair = (Pair_T) pairs->first;
      memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    }
    Sequence_free(&genomicseg);

    pairs = Pairpool_push_gapalign(pairs,pairpool,querypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,/*extraexonp*/false);
    pair = (Pair_T) pairs->first;
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));

    if (watsonp == true) {
      left = firstpart->chroffset + firstpart->chrpos + start2->genomepos - ngap;
    } else {
      left = firstpart->chroffset + firstpart->chrpos + (firstpart->genomiclength - 1) - start2->genomepos + 1;
    }

    genomicseg = Genome_get_segment(genome,left,ngap,/*chromosome_iit*/NULL,!watsonp,gbuffer1,gbuffer2,ngap);
    genomicseg_ptr = Sequence_fullpointer(genomicseg);

    for (i = 0; i < ngap; i++) {
      pairs = Pairpool_push_gapalign(pairs,pairpool,querypos,genomicpos,' ',comp,genomicseg_ptr[i],/*extraexonp*/false);
      pair = (Pair_T) pairs->first;
      memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    }
    Sequence_free(&genomicseg);
	
    FREE(gbuffer2);
    FREE(gbuffer1);
  }
  FREE(oldpairs);		/* Must free after last access to end1 */

  memcpy(&(firstpart->pairarray[ptr]),secondpart->pairarray,secondpart->npairs*sizeof(struct Pair_T));
  firstpart->npairs = firstpart->npairs+ngap+3+ngap+secondpart->npairs;

  firstpart->matches += secondpart->matches;
  firstpart->unknowns += secondpart->unknowns;
  firstpart->mismatches += secondpart->mismatches;
  firstpart->qopens += secondpart->qopens;
  firstpart->qindels += secondpart->qindels;
  firstpart->topens += secondpart->topens;
  firstpart->tindels += secondpart->tindels;

  return;
}

