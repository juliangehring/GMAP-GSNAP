static char rcsid[] = "$Id: stage3.c,v 1.218 2006/04/07 17:06:02 twu Exp $";
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
#include "chrnum.h"
#include "genomicpos.h"
#include "matchpair.h"
#include "smooth.h"
#include "scores.h"
#include "intron.h"
#include "translation.h"
#ifdef PMAP
#include "backtranslation.h"
#endif


/* The following are the same as dividing by 2048 and 1024 */
#define goodness_intronlen(x) (x >> 11)
#define goodness_nonintronlen(x) (x >> 10)

#define MAXITER 4
#define INTRON_PENALTY_INCONSISTENT 16
#define MININTRONLEN 6

#define MIN_NONINTRON 50

#define SUFFCONSECUTIVE 2

/* For Stage3_append */
#define MERGELENGTH 100000
#define LONG_MERGELENGTH 500000	/* For strong donor and acceptor splice sites */
#define DONOR_THRESHOLD 0.90
#define ACCEPTOR_THRESHOLD 0.90

/* Turn off CHECK for speed */
/* #define CHECK 1 */
static const Except_T gapcheck_error = {"Gap check failed"};


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

/* Show stage 2 */
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

/* Show result of smoothing.  May want to turn on DEBUG in smooth.c. */
#ifdef DEBUG5
#define debug5(x) x
#else 
#define debug5(x)
#endif

/* substitute for gaps */
#ifdef DEBUG7
#define debug7(x) x
#else 
#define debug7(x)
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
 *   For stage 1 and stage 2, the poly-A/T tails, if any, was stripped
 *   off, but for stage 3, we try to extend these tails if possible.
 *   Therefore, we use the full length and full sequence here, and add
 *   the offset to the path from stage 2.
 ************************************************************************/

#define T Stage3_T
struct T {
  struct Pair_T *pairarray;	/* The array version of pairs_fwd or pairs_rev, with the gaps substituted */
  int npairs;

  List_T pairs_fwd;		/* Pointer to list version for making a revised copy */
  List_T pairs_rev;		/* Pointer to list version for making a revised copy */

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
  double coverage;
  int matches;
  int unknowns;
  int mismatches;
  int qopens;
  int qindels;
  int topens;
  int tindels;
  int goodness;
  int nexons;
  /* int coverage_correction; */

  int translation_start;
  int translation_end;
  int translation_length;

  int relaastart;
  int relaaend;

  Matchpairend_T matchpairend;
};



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


double
Stage3_coverage (T this) {
  return this->coverage;
}

Matchpairend_T
Stage3_matchpairend (T this) {
  return this->matchpairend;
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


int *
Stage3_matchscores (T this, int querylength) {
  return Pair_matchscores(this->pairarray,this->npairs,this->cdna_direction,querylength);
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
 *  Stage3_new
 ************************************************************************/

static List_T
add_dualbreak (List_T pairs, int rightquerypos, int rightgenomepos, Pairpool_T pairpool) {
  int ngap = 3, genomicpos, querypos, k;

  for (k = 0; k < ngap; k++) {
    pairs = Pairpool_push(pairs,pairpool,rightquerypos,rightgenomepos,' ',DUALBREAK_COMP,' ',/*gapp*/true);
  }
  for (k = 0; k < 3; k++) {
    pairs = Pairpool_push(pairs,pairpool,rightquerypos,rightgenomepos,' ',INTRONGAP_COMP,' ',/*gapp*/true);
  }
  for (k = 0; k < ngap; k++) {
    pairs = Pairpool_push(pairs,pairpool,rightquerypos,rightgenomepos,' ',DUALBREAK_COMP,' ',/*gapp*/true);
  }
  return pairs;
}

static List_T
add_intron_old (List_T pairs, char *queryseq_ptr, char *genomicseg_ptr,
		int querypos, int rightquerypos, int genomepos, int rightgenomepos,
		int introntype, int ngap, Pairpool_T pairpool) {
  char comp, c2;
  int intronlength, genomicpos;
  int i;

  intronlength = rightgenomepos - genomepos;
  switch (introntype) {
  case GTAG_FWD: comp = FWD_CANONICAL_INTRON_COMP; break;
  case GCAG_FWD: comp = FWD_GCAG_INTRON_COMP; break;
  case ATAC_FWD: comp = FWD_ATAC_INTRON_COMP; break;
  case NONINTRON:
    if (intronlength < MIN_NONINTRON) {
      comp = SHORTGAP_COMP;	/* Will be printed as INDEL_COMP, but need to score as NONINTRON_COMP */
    } else {
      comp = NONINTRON_COMP;
    }
    break;
  case ATAC_REV: comp = REV_ATAC_INTRON_COMP; break;
  case GCAG_REV: comp = REV_GCAG_INTRON_COMP; break;
  case GTAG_REV: comp = REV_CANONICAL_INTRON_COMP; break;
  default: fprintf(stderr,"Unexpected intron type %d\n",introntype);
    exit(9);
  }
  debug7(printf("Adding gap of type %c of length %d\n",comp,intronlength));

  if (comp == SHORTGAP_COMP) {
    /* Treat as an insertion rather than an intron */
    for (i = 0, genomicpos = rightgenomepos - 1; i < intronlength; i++, genomicpos--) {
      c2 = genomicseg_ptr[genomicpos];
      pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',comp,c2,/*gapp*/false);
    }
  } else {
    if (ngap + ngap + 3 > intronlength) {
      for (i = 0, genomicpos = rightgenomepos - 1; i < intronlength; i++, genomicpos--) {
	c2 = genomicseg_ptr[genomicpos];
	pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',comp,c2,/*gapp*/true);
      }
    } else {
      for (i = 0, genomicpos = rightgenomepos - 1; i < ngap; i++, genomicpos--) {
	c2 = genomicseg_ptr[genomicpos];
	pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',comp,c2,/*gapp*/true);
	debug7(printf("Pushing %c at genomicpos %d\n",c2,genomicpos));
      }
    
      pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,/*gapp*/true);
      pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,/*gapp*/true);
      pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,/*gapp*/true);
    
      for (i = ngap-1, genomicpos = genomepos + ngap; i >= 0; i--, genomicpos--) {
	c2 = genomicseg_ptr[genomicpos];
	pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',comp,c2,/*gapp*/true);
	debug7(printf("Pushing %c at genomicpos %d\n",c2,genomicpos));
      }
    }
  }

  return pairs;
}

static List_T
add_intron (List_T pairs, char *queryseq_ptr, char *genomicseg_ptr,
	    int leftquerypos, int rightquerypos, int leftgenomepos, int rightgenomepos,
	    char comp, int ngap, Pairpool_T pairpool) {
  char c2;
  int intronlength, genomicpos;
  int i;

  intronlength = rightgenomepos - leftgenomepos;
  debug7(printf("Adding gap of type %c of length %d\n",comp,intronlength));

  if (ngap + ngap + 3 > intronlength) {
    for (i = 0, genomicpos = rightgenomepos - 1; i < intronlength; i++, --genomicpos) {
      c2 = genomicseg_ptr[genomicpos];
      pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',comp,c2,/*gapp*/true);
    }
  } else {
    for (i = 0, genomicpos = rightgenomepos - 1; i < ngap; i++, --genomicpos) {
      c2 = genomicseg_ptr[genomicpos];
      pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',comp,c2,/*gapp*/true);
      debug7(printf("Pushing %c at genomicpos %d\n",c2,genomicpos));
    }
    
    pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,/*gapp*/true);
    pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,/*gapp*/true);
    pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,/*gapp*/true);
    
    for (i = ngap-1, genomicpos = leftgenomepos + ngap; i >= 0; --i, --genomicpos) {
      c2 = genomicseg_ptr[genomicpos];
      pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',comp,c2,/*gapp*/true);
      debug7(printf("Pushing %c at genomicpos %d\n",c2,genomicpos));
    }
  }

  return pairs;
}

static List_T
assign_gap_types (List_T path, Pairpool_T pairpool, char *queryseq_ptr, char *genomicseg_ptr, 
		  char *genomicuc_ptr, int cdna_direction) {
  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;
  int queryjump, genomejump, leftquerypos, leftgenomepos, rightquerypos, rightgenomepos, curquerypos,
    introntype, intronlength, genomicpos;
  char left1, left2, right2, right1, c2;

  debug(printf("\n** Starting assign_gap_types\n"));
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == false) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
    } else {
      queryjump = pair->queryjump;
      genomejump = pair->genomejump;

      if (queryjump == 0 && genomejump == 0) {
	debug7(printf("  Gap is a non-gap\n"));
	/* Discard the gap pair */

      } else if (genomejump == 0) {
	debug7(printf("  Gap is a cDNA insertion\n"));
	pair->comp = INDEL_COMP;

	leftpair = path->first;
	rightpair = pairs->first;
	leftquerypos = leftpair->querypos;
	rightquerypos = rightpair->querypos;

	for (curquerypos = rightquerypos - 1; curquerypos > leftquerypos; --curquerypos) {
	  path = Pairpool_push(path,pairpool,curquerypos,rightgenomepos,
			       queryseq_ptr[curquerypos],INDEL_COMP,' ',/*gapp*/false);
	}
	/* Discard the gap pair */

      } else if (queryjump > 0) {
	debug7(printf("  Gap is a dual break\n"));
	pair->comp = DUALBREAK_COMP;
	pairs = Pairpool_push_existing(pairs,pairpool,pair);

      } else {
	debug7(printf("  Gap is an intron\n"));

	leftpair = path->first;
	rightpair = pairs->first;
	leftquerypos = leftpair->querypos;
	leftgenomepos = leftpair->genomepos;
	rightquerypos = rightpair->querypos;
	rightgenomepos = rightpair->genomepos;

	left1 = genomicuc_ptr[leftgenomepos+1];
	left2 = genomicuc_ptr[leftgenomepos+2];
	right2 = genomicuc_ptr[rightgenomepos-2];
	right1 = genomicuc_ptr[rightgenomepos-1];
	  
	introntype = Intron_type(left1,left2,right2,right1,cdna_direction);
	switch (introntype) {
	case GTAG_FWD: pair->comp = FWD_CANONICAL_INTRON_COMP; break;
	case GCAG_FWD: pair->comp = FWD_GCAG_INTRON_COMP; break;
	case ATAC_FWD: pair->comp = FWD_ATAC_INTRON_COMP; break;
	case ATAC_REV: pair->comp = REV_ATAC_INTRON_COMP; break;
	case GCAG_REV: pair->comp = REV_GCAG_INTRON_COMP; break;
	case GTAG_REV: pair->comp = REV_CANONICAL_INTRON_COMP; break;
	case NONINTRON:
	  intronlength = rightgenomepos - leftgenomepos;
	  if (intronlength < MIN_NONINTRON) {
	    pair->comp = SHORTGAP_COMP;	/* Will be printed as INDEL_COMP, but need to score as NONINTRON_COMP */
	  } else {
	    pair->comp = NONINTRON_COMP;
	  }
	  break;
	default: fprintf(stderr,"Unexpected intron type %d\n",introntype);
	  exit(9);
	}
	if (pair->comp == SHORTGAP_COMP) {
	  for (genomicpos = rightgenomepos - 1; genomicpos > leftgenomepos; --genomicpos) {
	    c2 = genomicseg_ptr[genomicpos];
	    pairs = Pairpool_push(pairs,pairpool,rightquerypos,genomicpos,' ',/*comp*/SHORTGAP_COMP,c2,/*gapp*/false);
	  }
	  /* Discard the gap pair */
	} else {
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
	}
      }
    }
  }

  return pairs;
}


static List_T 
substitute_for_gaps (List_T path, Pairpool_T pairpool, char *queryseq_ptr, char *genomicseg_ptr,
		     char *genomicuc_ptr, int cdna_direction, int ngap) {
  List_T pairs = NULL;
  int leftquerypos, leftgenomepos, rightquerypos, rightgenomepos;
  Pair_T pair, leftpair, rightpair;

  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == false) {
      debug7(printf("At %d %d, gapp is false\n",pair->querypos,pair->genomepos));
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else {
      /* Discard gap; do not push */

      if (pair->comp == INDEL_COMP) {
	/* This is already handled by assign_gap_types() */
	abort();
      } else if (pair->comp == SHORTGAP_COMP) {
	/* This is already handled by assign_gap_types() */
	abort();
      } else if (pair->comp == DUALBREAK_COMP) {
	rightpair = pairs->first;
	rightquerypos = rightpair->querypos;
	rightgenomepos = rightpair->genomepos;

	pairs = add_dualbreak(pairs,rightquerypos,rightgenomepos,pairpool);
      } else {
	leftpair = path->first;
	rightpair = pairs->first;
	leftquerypos = leftpair->querypos;
	leftgenomepos = leftpair->genomepos;
	rightquerypos = rightpair->querypos;
	rightgenomepos = rightpair->genomepos;

	pairs = add_intron(pairs,queryseq_ptr,genomicseg_ptr,leftquerypos,rightquerypos,
			  leftgenomepos,rightgenomepos,pair->comp,ngap,pairpool);
      }
    }
  }

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
		      char *queryseq_ptr, char *genomicseg_ptr, char *genomicuc_ptr, int ngap,
		      int subseq_offset, int skiplength) {
  struct Pair_T *pairarray;
  List_T path, p;
  Pair_T oldpair, newpair;

  path = List_reverse(*pairs);
  *pairs = substitute_for_gaps(path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,ngap);

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

	    
static T
Stage3_new (List_T pairs_fwd, List_T pairs_rev) {
  T new;
  int matches_fwd, matches_rev, mismatches_fwd, mismatches_rev, 
    unknowns_fwd, unknowns_rev, qopens_fwd, qindels_fwd, qopens_rev, qindels_rev, 
    topens_fwd, tindels_fwd, topens_rev, tindels_rev,
    ncanonical_fwd, ncanonical_rev, nsemicanonical_fwd, nsemicanonical_rev,
    nnoncanonical_fwd, nnoncanonical_rev;
  int goodness_fwd, goodness_rev;

  if (pairs_fwd == NULL && pairs_rev == NULL) {
    return NULL;
  } else {
    Pair_fracidentity(&matches_fwd,&unknowns_fwd,&mismatches_fwd,
		      &qopens_fwd,&qindels_fwd,&topens_fwd,&tindels_fwd,
		      &ncanonical_fwd,&nsemicanonical_fwd,&nnoncanonical_fwd,pairs_fwd,+1);
    Pair_fracidentity(&matches_rev,&unknowns_rev,&mismatches_rev,
		      &qopens_rev,&qindels_rev,&topens_rev,&tindels_rev,
		      &ncanonical_rev,&nsemicanonical_rev,&nnoncanonical_rev,pairs_rev,-1);
    /* For determining direction, a canonical intron is worth three
       mismatch.  But for comparing among paths, do not consider
       number of canonical introns, but do consider nonintronlen.
    */
    goodness_fwd = matches_fwd + MISMATCH*mismatches_fwd + QOPEN*qopens_fwd + QINDEL*qindels_fwd + TOPEN*topens_fwd + TINDEL*tindels_fwd + CANONICAL_POINTS*ncanonical_fwd + SEMICANONICAL_POINTS*nsemicanonical_fwd - CANONICAL_POINTS*nnoncanonical_fwd;
    goodness_rev = matches_rev + MISMATCH*mismatches_rev + QOPEN*qopens_rev + QINDEL*qindels_rev + TOPEN*topens_rev + TINDEL*tindels_rev + CANONICAL_POINTS*ncanonical_rev + SEMICANONICAL_POINTS*nsemicanonical_rev - CANONICAL_POINTS*nnoncanonical_rev;

    debug(printf("Stage 3: goodness_fwd %d = %d - 3*%d - 5*%d - 2*%d - 5*%d - 2*%d + 12*%d + 6*%d - 12*%d\n",
		 goodness_fwd,matches_fwd,mismatches_fwd,qopens_fwd,qindels_fwd,topens_fwd,tindels_fwd,ncanonical_fwd,nsemicanonical_fwd,nnoncanonical_fwd));
    debug(printf("Stage 3: goodness_rev %d = %d - 3*%d - 5*%d - 2*%d - 5*%d - 2*%d + 12*%d + 6*%d - 12*%d\n",
		 goodness_rev,matches_rev,mismatches_rev,qopens_rev,qindels_rev,topens_rev,tindels_rev,ncanonical_rev,nsemicanonical_rev,nnoncanonical_rev));
    debug(printf("Stage 3: Deciding between goodness_fwd (%d canonical, score %d) and goodness_rev (%d canonical, score %d)\n",
		 ncanonical_fwd,goodness_fwd,ncanonical_rev,goodness_rev));

    new = (T) MALLOC(sizeof(*new));
    new->pairs_fwd = pairs_fwd;
    new->pairs_rev = pairs_rev;

    if (pairs_rev == NULL) {
      new->cdna_direction = +1;
    } else if (pairs_fwd == NULL) {
      new->cdna_direction = -1;
    } else if (goodness_fwd > goodness_rev) {
      new->cdna_direction = +1;
    } else if (goodness_fwd < goodness_rev) {
      new->cdna_direction = -1;
    } else {
      new->cdna_direction = 0;	/* indeterminate */
    }

    if (new->cdna_direction >= 0) {
      new->matches = matches_fwd;
      new->unknowns = unknowns_fwd;
      new->mismatches = mismatches_fwd;
      new->qopens = qopens_fwd;
      new->qindels = qindels_fwd;
      new->topens = topens_fwd;
      new->tindels = tindels_fwd;
      new->goodness = goodness_fwd - CANONICAL_POINTS*ncanonical_fwd - CANONICAL_POINTS*nnoncanonical_fwd;
      new->nexons = Pair_nexons(new->pairs_fwd);
    } else {
      new->matches = matches_rev;
      new->unknowns = unknowns_rev;
      new->mismatches = mismatches_rev;
      new->qopens = qopens_rev;
      new->qindels = qindels_rev; 
      new->topens = topens_rev;
      new->tindels = tindels_rev;
      new->goodness = goodness_rev - CANONICAL_POINTS*ncanonical_rev - CANONICAL_POINTS*nnoncanonical_rev;
      new->nexons = Pair_nexons(new->pairs_rev);
    }

    if (new->nexons > 2) {
      /* Favor spliced transcripts, but only if we're sure they're
         spliced (i.e., 3 or more exons).  A random intron shouldn't
         get credit. */
      new->goodness += new->nexons;
    }
    
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
  struct Pair_T *newpairarray = NULL;
  int newnpairs;
  List_T p;
  Pair_T newpair, oldpair;

  if (this->cdna_direction >= 0) {
    if (Pairpool_count_bounded(this->pairs_fwd,minpos,maxpos) > 25) {
      return true;
    } else {
      return false;
    }
  } else {
    if (Pairpool_count_bounded(this->pairs_rev,minpos,maxpos) > 25) {
      return true;
    } else {
      return false;
    }
  }
}

T
Stage3_apply_bounds (T this, int minpos, int maxpos, bool revertp) {
  struct Pair_T *newpairarray = NULL;
  int newnpairs;
  List_T p;
  Pair_T newpair, oldpair;

  if (this->cdna_direction >= 0) {
    this->pairs_fwd = Pairpool_transfer_bounded(NULL,this->pairs_fwd,minpos,maxpos);
    this->pairs_fwd = List_reverse(this->pairs_fwd);
    if ((newnpairs = List_length(this->pairs_fwd)) == 0) {
      if (revertp == false) {
	return NULL;
      }
    } else {
      FREE(this->pairarray);
      this->npairs = newnpairs;
      
      newpair = this->pairarray = (struct Pair_T *) MALLOC(newnpairs*sizeof(struct Pair_T));
      for (p = this->pairs_fwd; p != NULL; p = p->rest) {
	oldpair = (Pair_T) p->first;
	memcpy(newpair++,oldpair,sizeof(struct Pair_T));
      }
    }
  } else {
    this->pairs_rev = Pairpool_transfer_bounded(NULL,this->pairs_rev,minpos,maxpos);
    this->pairs_rev = List_reverse(this->pairs_rev);
    if ((newnpairs = List_length(this->pairs_rev)) == 0) {
      if (revertp == false) {
	return NULL;
      }
    } else {
      FREE(this->pairarray);
      this->npairs = newnpairs;
      
      newpair = this->pairarray = (struct Pair_T *) MALLOC(newnpairs*sizeof(struct Pair_T));
      for (p = this->pairs_rev; p != NULL; p = p->rest) {
	oldpair = (Pair_T) p->first;
	memcpy(newpair++,oldpair,sizeof(struct Pair_T));
      }
    }
  }

  return this;
}


#ifdef PMAP
Stage3_T
Stage3_translate_cdna (T this, Sequence_T queryaaseq) {
  Translation_via_cdna(&this->translation_start,&this->translation_end,&this->translation_length,
		       &this->relaastart,&this->relaaend,
		       this->pairarray,this->npairs,Sequence_fullpointer(queryaaseq));
  return this;
}

Stage3_T
Stage3_backtranslate_cdna (T this, bool diagnosticp) {
  Backtranslation_cdna(this->pairarray,this->npairs,this->translation_start,this->translation_end,
		       diagnosticp);
  return this;
}

#else
Stage3_T
Stage3_truncate_fulllength (Stage3_T this, bool translatep) {
  struct Pair_T *newpairarray;
  int newnpairs;

  if (translatep == true) {
    if (this->cdna_direction < 0) {
      Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			      &this->relaastart,&this->relaaend,
			      this->pairarray,this->npairs,/*backwardsp*/true,/*revcompp*/true,/*fulllengthp*/true);
    } else {
      Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			      &this->relaastart,&this->relaaend,
			      this->pairarray,this->npairs,/*backwardsp*/false,/*revcompp*/false,/*fulllengthp*/true);
    }
  }

  this = Stage3_apply_bounds(this,Stage3_translation_start(this),Stage3_translation_end(this),/*revertp*/true);
  return this;
}

Stage3_T
Stage3_translate_genomic (T this, bool fulllengthp, bool truncatep) {
  Stage3_T new;

  if (this->cdna_direction < 0) {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairarray,this->npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp);
  } else {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairarray,this->npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp);
  }
  if (truncatep == false) {
    return this;
  } else {
    new = Stage3_truncate_fulllength(this,/*translatep*/false);
    if (new->cdna_direction < 0) {
      Translation_via_genomic(&new->translation_start,&new->translation_end,&new->translation_length,
			      &new->relaastart,&new->relaaend,
			      new->pairarray,new->npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp);
    } else {
      Translation_via_genomic(&new->translation_start,&new->translation_end,&new->translation_length,
			      &new->relaastart,&new->relaaend,
			      new->pairarray,new->npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp);
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


/* queryaaseq is used only by PMAP */
void
Stage3_print_pathsummary (T this, int pathnum, IIT_T chromosome_iit, IIT_T contig_iit, 
			  IIT_T altstrain_iit, Sequence_T queryaaseq, bool fulllengthp, bool truncatep,
			  char *dbversion, int maxmutations, bool zerobasedp, bool diagnosticp,
			  bool maponlyp) {
  Pair_T start, end;
  bool referencealignp;

  if (maponlyp == true) {
    this->translation_start = 0;
    this->translation_end = 0;
    this->translation_length = 0;
  } else {
#ifdef PMAP
    Stage3_translate_cdna(this,queryaaseq);
    Stage3_backtranslate_cdna(this,diagnosticp);
#else
    Stage3_translate_genomic(this,fulllengthp,truncatep);
#endif
  }

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);
  referencealignp = this->straintype == 0 ? true : false;
  Pair_print_pathsummary(pathnum,start,end,this->chrnum,this->chrpos,this->chroffset,
			 chromosome_iit,referencealignp,altstrain_iit,this->strain,contig_iit,
			 dbversion,this->genomiclength,
			 this->nexons,this->coverage,
			 this->matches,this->unknowns,this->mismatches,
			 this->qopens,this->qindels,this->topens,this->tindels,this->goodness,
			 this->watsonp,this->cdna_direction,this->defect_rate,
			 this->translation_start,this->translation_end,this->translation_length,
			 0,0,zerobasedp,maponlyp);
  if (maponlyp == false) {
    Translation_compare(this->pairarray,this->npairs,NULL,0,this->cdna_direction,
			this->relaastart,this->relaaend,maxmutations);
  }
  printf("\n");

  return;
}

void
Stage3_print_pslformat_nt (T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryaaseq) {
  Pair_T start, end;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_pslformat_nt(this->pairarray,this->npairs,start,end,queryaaseq,this->chrnum,this->chrpos,
			  chromosome_iit,this->genomiclength,this->nexons,this->matches,this->unknowns,this->mismatches, 
			  this->qopens,this->qindels,this->topens,this->tindels,
			  this->watsonp,this->cdna_direction);
  return;
}

#ifdef PMAP
void
Stage3_print_pslformat_pro (T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryaaseq) {
  Pair_T start, end;

  Stage3_translate_cdna(this,queryaaseq);
  Stage3_backtranslate_cdna(this,/*diagnosticp*/false);

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  Pair_print_pslformat_pro(this->pairarray,this->npairs,start,end,queryaaseq,this->chrnum,this->chrpos,
			   chromosome_iit,this->genomiclength,this->nexons,
			   this->qopens,this->qindels,this->topens,this->tindels,
			   this->watsonp,this->cdna_direction);
  return;
}
#endif

void
Stage3_print_mutations (T this, T reference, IIT_T chromosome_iit, char *dbversion,
			bool showalignp, bool zerobasedp, 
			bool continuousp, bool diagnosticp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength, int ngap,
			int maxmutations) {
  Pair_T start, end;
  bool referencealignp;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  /*  Pair_dump_array(this->pairarray,this->npairs,false); */

  referencealignp = this->straintype == 0 ? true : false;
  Pair_print_pathsummary(/*pathnum*/1,start,end,reference->chrnum,reference->chrpos,reference->chroffset,
			 chromosome_iit,referencealignp,/*altstrain_iit*/NULL,this->strain,/*contig_iit*/NULL,
			 dbversion,reference->genomiclength,
			 this->nexons,this->coverage,
			 this->matches,this->unknowns,this->mismatches,
			 this->qopens,this->qindels,this->topens,this->tindels,this->goodness,
			 this->watsonp,this->cdna_direction,this->defect_rate,
			 0,0,0,this->relaastart,this->relaaend,zerobasedp,/*maponlyp*/false);
  Translation_compare(this->pairarray,this->npairs,reference->pairarray,reference->npairs,
		      this->cdna_direction,this->relaastart,this->relaaend,maxmutations);
  printf("\n");

  if (showalignp == true) {
    Pair_print_alignment(this->pairarray,this->npairs,reference->chrnum,reference->chrpos,reference->chroffset,
			 chromosome_iit,reference->genomiclength,this->watsonp,this->cdna_direction,/*universalp*/false,zerobasedp,
			 diagnosticp,/*genomicprop*/false,invertmode,nointronlenp,wraplength,ngap);
  }
  debug1(Pair_dump_array(this->pairarray,this->npairs,/*zerobasedp*/true));
  debug1(Pair_check_array(this->pairarray,this->npairs));

  return;
}


static void
print_map (T this, IIT_T map_iit, IIT_T chromosome_iit,
	   int pathnum, bool map_bothstrands_p) {
  Genomicpos_T position1, position2;
  Pair_T start, end;
  int chrpos1, chrpos2;
  int *iit_matches = NULL, nmatches;
  int typeint;
  char *typestring;

  start = &(this->pairarray[0]);
  end = &(this->pairarray[this->npairs-1]);

  if (this->watsonp) {
    chrpos1 = this->chrpos + Pair_genomepos(start);
    chrpos2 = this->chrpos + Pair_genomepos(end);
    if (this->cdna_direction >= 0) {
      typestring = "FWD";
    } else {
      typestring = "REV";
    }
  } else {
    chrpos1 = this->chrpos + (this->genomiclength - 1) - Pair_genomepos(start);
    chrpos2 = this->chrpos + (this->genomiclength - 1) - Pair_genomepos(end);
    if (this->cdna_direction >= 0) {
      typestring = "REV";
    } else {
      typestring = "FWD";
    }
  }

  position1 = this->chroffset + chrpos1;
  position2 = this->chroffset + chrpos2;

  if (map_bothstrands_p == true) {
    if (position1 < position2) {
      iit_matches = IIT_get(&nmatches,map_iit,position1,position2);
    } else {
      iit_matches = IIT_get(&nmatches,map_iit,position2,position1);
    }
    printf("  Map hits for path %d (%d):\n",pathnum,nmatches);
    IIT_print(map_iit,iit_matches,nmatches,/*bothstrandsp*/true,chromosome_iit,/*levels*/NULL);
  } else if ((typeint = IIT_typeint(map_iit,typestring)) < 0) {
    fprintf(stderr,"Warning: no type %s in %s.  Ignoring type.\n",typestring,IIT_name(map_iit));
    if (position1 < position2) {
      iit_matches = IIT_get(&nmatches,map_iit,position1,position2);
    } else {
      iit_matches = IIT_get(&nmatches,map_iit,position2,position1);
    }
    printf("  Map hits for path %d (%d):\n",pathnum,nmatches);
    IIT_print(map_iit,iit_matches,nmatches,/*bothstrandsp*/false,chromosome_iit,/*levels*/NULL);
  } else {
    if (position1 < position2) {
      iit_matches = IIT_get_typed(&nmatches,map_iit,position1,position2,typeint);
    } else {
      iit_matches = IIT_get_typed(&nmatches,map_iit,position2,position1,typeint);
    }
    printf("  Map hits for path %d (%d):\n",pathnum,nmatches);
    IIT_print(map_iit,iit_matches,nmatches,/*bothstrandsp*/false,chromosome_iit,/*levels*/NULL);
  }
  printf("\n");

  FREE(iit_matches);
  return;
}

static void
print_exon_map (T this, IIT_T map_iit, IIT_T chromosome_iit,
		int pathnum, bool map_bothstrands_p) {
  Uintlist_T exonbounds;
  Genomicpos_T position1, position2;
  int *iit_matches = NULL, nmatches;
  int typeint, exonno = 0;
  char *typestring;

  if (this->watsonp) {
    if (this->cdna_direction >= 0) {
      typestring = "FWD";
    } else {
      typestring = "REV";
    }
  } else {
    if (this->cdna_direction >= 0) {
      typestring = "REV";
    } else {
      typestring = "FWD";
    }
  }

  exonbounds = Pair_exonbounds(this->pairarray,this->npairs,this->chrpos,this->chroffset,this->genomiclength,
			       this->watsonp);

  while (exonbounds != NULL) {
    exonbounds = Uintlist_pop(exonbounds,&position1);
    exonbounds = Uintlist_pop(exonbounds,&position2);

    if (map_bothstrands_p == true) {
      if (position1 < position2) {
	iit_matches = IIT_get(&nmatches,map_iit,position1,position2);
      } else {
	iit_matches = IIT_get(&nmatches,map_iit,position2,position1);
      }
      printf("  Map hits for path %d, exon %d (%d):\n",pathnum,++exonno,nmatches);
      IIT_print(map_iit,iit_matches,nmatches,/*bothstrandsp*/true,chromosome_iit,/*levels*/NULL);
    } else if ((typeint = IIT_typeint(map_iit,typestring)) < 0) {
      fprintf(stderr,"Warning: no type %s in %s.  Ignoring type.\n",typestring,IIT_name(map_iit));
      if (position1 < position2) {
	iit_matches = IIT_get(&nmatches,map_iit,position1,position2);
      } else {
	iit_matches = IIT_get(&nmatches,map_iit,position2,position1);
      }
      printf("  Map hits for path %d, exon %d (%d):\n",pathnum,++exonno,nmatches);
      IIT_print(map_iit,iit_matches,nmatches,/*bothstrandsp*/false,chromosome_iit,/*levels*/NULL);
    } else {
      if (position1 < position2) {
	iit_matches = IIT_get_typed(&nmatches,map_iit,position1,position2,typeint);
      } else {
	iit_matches = IIT_get_typed(&nmatches,map_iit,position2,position1,typeint);
      }
      printf("  Map hits for path %d, exon %d (%d):\n",pathnum,++exonno,nmatches);
      IIT_print(map_iit,iit_matches,nmatches,/*bothstrandsp*/false,chromosome_iit,/*levels*/NULL);
    }
    printf("\n");
    FREE(iit_matches);
  }

  return;
}

void
Stage3_print_map (T this, IIT_T map_iit, IIT_T chromosome_iit,
		  int pathnum, bool map_exons_p, bool map_bothstrands_p) {
  if (map_exons_p == true) {
    print_exon_map(this,map_iit,chromosome_iit,pathnum,map_bothstrands_p);
  } else {
    print_map(this,map_iit,chromosome_iit,pathnum,map_bothstrands_p);
  }
  return;
}



/* queryaaseq is used only by PMAP */
void
Stage3_print_alignment (T this, Sequence_T queryaaseq,
			IIT_T chromosome_iit, bool alignsummaryonlyp, bool universalp, bool zerobasedp,
			bool continuousp, bool diagnosticp, bool genomefirstp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength, int ngap, bool maponlyp) {
  if (continuousp == true) {
    if (maponlyp == false) {
#ifdef PMAP
      Stage3_translate_cdna(this,queryaaseq);
      Stage3_backtranslate_cdna(this,diagnosticp);
#endif
    }
    Pair_print_continuous(this->pairarray,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			  this->genomiclength,this->watsonp,this->cdna_direction,universalp,zerobasedp,
			  diagnosticp,genomefirstp,invertmode,nointronlenp,ngap);
  } else {
    /* Assumes Stage3_print_pathsummary already called */
    Pair_print_exonsummary(this->pairarray,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			   chromosome_iit,this->genomiclength,this->watsonp,universalp,zerobasedp,
			   genomefirstp,invertmode);
    if (alignsummaryonlyp == false) {
      Pair_print_alignment(this->pairarray,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			   chromosome_iit,this->genomiclength,this->watsonp,
			   this->cdna_direction,universalp,zerobasedp,
			   diagnosticp,/*genomicprop*/true,invertmode,nointronlenp,wraplength,ngap);
    }
  }
  debug1(Pair_dump_array(this->pairarray,this->npairs,/*zerobasedp*/true));
  debug1(Pair_check_array(this->pairarray,this->npairs));
  return;
}


/* queryaaseq is used only by PMAP */
void
Stage3_print_coordinates (T this, Sequence_T queryaaseq, IIT_T chromosome_iit, bool zerobasedp,
			  int invertmode, bool fulllengthp, bool truncatep, bool maponlyp) {
  if (maponlyp == false) {
#ifdef PMAP
    Stage3_translate_cdna(this,queryaaseq);
    Stage3_backtranslate_cdna(this,/*diagnosticp*/false);
#else
    Stage3_translate_genomic(this,fulllengthp,truncatep);
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
Stage3_print_cdna (T this, Sequence_T queryaaseq, bool fulllengthp, bool truncatep, int wraplength) {
#ifdef PMAP
  Stage3_translate_cdna(this,queryaaseq);
  Stage3_backtranslate_cdna(this,/*diagnosticp*/false);
  Pair_print_nucleotide_cdna(this->pairarray,this->npairs,wraplength);
#else
  Stage3_translate_genomic(this,fulllengthp,truncatep);
  Pair_print_protein_cdna(this->pairarray,this->npairs,wraplength);
#endif
  return;
}

void
Stage3_print_protein_genomic (T this, Sequence_T queryaaseq, bool fulllengthp, bool truncatep,
			      int wraplength) {
#ifdef PMAP
  Stage3_translate_cdna(this,queryaaseq);
  Stage3_backtranslate_cdna(this,/*diagnosticp*/false);
#else
  Stage3_translate_genomic(this,fulllengthp,truncatep);
#endif
  Pair_print_protein_genomic(this->pairarray,this->npairs,wraplength);
  return;
}


void
Stage3_print_compressed (T this, Sequence_T queryseq, IIT_T chromosome_iit,
			 char *version, int pathnum, int npaths,
			 bool checksump, int chimerapos, int chimeraequivpos,
			 double donor_prob, double acceptor_prob, 
			 int chimera_cdna_direction, bool zerobasedp, bool truncatep, int worker_id) {
#ifdef PMAP
  Stage3_translate_cdna(this,queryseq);
  Stage3_backtranslate_cdna(this,/*diagnosticp*/false);
#else
  if (truncatep == true) {
    Stage3_truncate_fulllength(this,/*translatep*/true);
  }
#endif

  Pair_print_compressed(queryseq,version,pathnum,npaths,
			this->nexons,this->coverage,Stage3_fracidentity(this),
			this->pairarray,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			chromosome_iit,this->genomiclength,checksump,
			chimerapos,chimeraequivpos,donor_prob,acceptor_prob,
			chimera_cdna_direction,this->strain,this->watsonp,
			this->cdna_direction,zerobasedp,worker_id);
  return;
}



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

static List_T
peel_forward (List_T *peeled_path, List_T path, int *querydp5, int *genomedp5, 
	      Pairpool_T pairpool, int maxpeelback, bool intronp) {
  Pair_T pair, firstpair, leftpair;
  int npeelback = 0, nconsecutive = 0;
#ifdef PMAP
  bool midcodonp = false;
#endif
  List_T peeled = NULL;

  /* Peelback */
  debug(printf("Peeling forward:"));
  if (path == NULL) {
    debug(printf(" path is empty\n"));
  } else {
    firstpair = path->first;
    if (firstpair->gapp == true) {
      /* Known gap */
      debug(printf(" Known_gap"));
      path = Pairpool_pop(path,&pair);
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
      if (path != NULL) {
	firstpair = path->first;
      }
    }
    
    while (path != NULL && firstpair->gapp == false && npeelback < maxpeelback) {
      path = Pairpool_pop(path,&pair);
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
      debug(printf(" Peel [%d %d %c %c %c]",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      npeelback++;
      if (path != NULL) {
	firstpair = path->first;
      }
    }

    debug(printf(" ||"));
    if (firstpair->gapp == true) {
      if (peeled != NULL) {
	peeled = Pairpool_pop(peeled,&pair);
	path = Pairpool_push_existing(path,pairpool,pair);
	debug(printf(" Putback [%d %d %c %c %c]",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      }
    } else if (intronp) {
      /* Continue to peelback through little skips and mismatches */
      while (path != NULL && firstpair->gapp == false && nconsecutive < SUFFCONSECUTIVE) {
	path = Pairpool_pop(path,&pair);
	peeled = Pairpool_push_existing(peeled,pairpool,pair);
	debug(printf(" Extrapeel [%d %d %c %c %c]",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	if (pair->comp == INDEL_COMP || pair->comp == MISMATCH_COMP) {
	  nconsecutive = 0;
	} else {
	  nconsecutive++;
	}
	if (path != NULL) {
	  firstpair = path->first;
	}
      }
      if (firstpair->gapp == true) {
	if (peeled != NULL) {
	  peeled = Pairpool_pop(peeled,&pair);
	  path = Pairpool_push_existing(path,pairpool,pair);
	  debug(printf(" Putback [%d %d %c %c %c]",
		       pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	}
      }
    }

    debug(printf(" ||"));
#ifdef PMAP
    /* Peel to codon boundary */
    if (peeled != NULL) {
      firstpair = peeled->first;
      while (peeled != NULL && firstpair->querypos % 3 != 0) {
	peeled = Pairpool_pop(peeled,&pair);
	path = Pairpool_push_existing(path,pairpool,pair);
	debug(printf(" Mod3putback [%d %d %c %c %c]",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	if (peeled != NULL) {
	  firstpair = peeled->first;
	}
      }
    }
#endif
  }

  debug(
	if (path == NULL) {
	  printf("Top of path is NULL\n");
	} else {
	  pair = path->first;
	  printf("Top of path is %d %d %c %c %c\n",
		 pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	}
	);

  if (path != NULL) {
    leftpair = path->first;
    *querydp5 = leftpair->querypos + 1;
    *genomedp5 = leftpair->genomepos + 1;
  } else if (peeled != NULL) {
    leftpair = peeled->first;
    *querydp5 = leftpair->querypos;
    *genomedp5 = leftpair->genomepos;
  }

  *peeled_path = peeled;
  return path;
}


static List_T
peel_back (List_T *peeled_pairs, List_T pairs, int *querydp3, int *genomedp3, 
	   Pairpool_T pairpool, int maxpeelback, bool intronp) {
  Pair_T pair, firstpair, rightpair;
  int npeelback = 0, nconsecutive = 0;
#ifdef PMAP
  bool midcodonp = false;
#endif
  List_T peeled = NULL;

  /* Peelback */
  debug(printf("Peeling back:"));
  if (pairs == NULL) {
    debug(printf(" pairs is empty\n"));
  } else {
    firstpair = pairs->first;
    if (firstpair->gapp == true) {
      /* Throw away known gap */
      debug(printf(" Known_gap"));
      pairs = Pairpool_pop(pairs,&pair);
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
      if (pairs != NULL) {
	firstpair = pairs->first;
      }
    }

    while (pairs != NULL && firstpair->gapp == false && npeelback < maxpeelback) {
      pairs = Pairpool_pop(pairs,&pair);
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
      debug(printf(" Peel [%d %d %c %c %c]",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      npeelback++;
      if (pairs != NULL) {
	firstpair = pairs->first;
      }
    }

    debug(printf(" ||"));
    if (firstpair->gapp == true) {
      if (peeled != NULL) {
	peeled = Pairpool_pop(peeled,&pair);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
	debug(printf(" Putback [%d %d %c %c %c]",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      }
    } else if (intronp) {
      /* Continue to peelback through little skips and mismatches */
      while (pairs != NULL && firstpair->gapp == false && nconsecutive < SUFFCONSECUTIVE) {
	pairs = Pairpool_pop(pairs,&pair);
	peeled = Pairpool_push_existing(peeled,pairpool,pair);
	debug(printf(" Extrapeel [%d %d %c %c %c]",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	if (pair->comp == INDEL_COMP || pair->comp == MISMATCH_COMP) {
	  nconsecutive = 0;
	} else {
	  nconsecutive++;
	}
	if (pairs != NULL) {
	  firstpair = pairs->first;
	}
      }
      if (firstpair->gapp == true) {
	if (peeled != NULL) {
	  peeled = Pairpool_pop(peeled,&pair);
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
	  debug(printf(" Putback [%d %d %c %c %c]",
		       pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	}
      }
    }

    debug(printf(" ||"));
#ifdef PMAP
    /* Peel to codon boundary */
    if (peeled != NULL) {
      firstpair = peeled->first;
      while (peeled != NULL && firstpair->querypos % 3 != 2) {
	peeled = Pairpool_pop(peeled,&pair);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
	debug(printf(" Mod3putback [%d %d %c %c %c]",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	if (peeled != NULL) {
	  firstpair = peeled->first;
	}
      }
    }
#endif
  }

  debug(
	if (pairs == NULL) {
	  printf("Top of pairs is NULL\n");
	} else {
	  pair = pairs->first;
	  printf("Top of pairs is %d %d %c %c %c\n",
		 pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	}
	);

  if (pairs != NULL) {
    rightpair = pairs->first;
    *querydp3 = rightpair->querypos - 1;
    *genomedp3 = rightpair->genomepos - 1;
  } else if (peeled != NULL) {
    rightpair = peeled->first;
    *querydp3 = rightpair->querypos;
    *genomedp3 = rightpair->genomepos;
  }

  *peeled_pairs = peeled;
  return pairs;
}


/************************************************************************
 *  Traversal functions
 ************************************************************************/

static List_T
traverse_single_gap (bool *filledp, List_T pairs, List_T *path, 
		     int leftquerypos, int leftgenomepos, int rightquerypos, int rightgenomepos,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprog,
		     int maxpeelback, int extraband_single, double defect_rate, bool forcep) {
  List_T gappairs, peeled_pairs, peeled_path;
  Pair_T leftpair, rightpair;
  int queryjump, genomejump;
  int querydp5 = leftquerypos + 1;
  int genomedp5 = leftgenomepos + 1;
  int querydp3 = rightquerypos - 1;
  int genomedp3 = rightgenomepos - 1;
  int nmatches, nmismatches, nopens, nindels;
  int finalscore;

  /* Note that we peelback only half as much as for a paired gap, to
     save on dynamic programming */
  pairs = peel_back(&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback/2,
		    /*intronp*/false);
  *path = peel_forward(&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback/2,
		       /*intronp*/false);

  queryjump = querydp3 - querydp5 + 1;
  genomejump = genomedp3 - genomedp5 + 1;
  
  if (queryjump == 0 || genomejump == 0) {
    debug(printf("Unable to perform dynamic programming\n"));
    *filledp = false;
    return NULL;
  } else {
    gappairs = Dynprog_single_gap(&finalscore,&nmatches,&nmismatches,&nopens,&nindels,dynprog,
				  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				  &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				  queryjump,genomejump,querydp5,genomedp5,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  cdna_direction,pairpool,extraband_single,defect_rate);
  }
  debug(Pair_dump_list(gappairs,true));
  debug(printf("  Score: %d\n",finalscore));

  if (!forcep && finalscore < 0) {
    *filledp = false;
    /* Put back peeled pairs */
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
  } else {
    *filledp = true;
    pairs = Pairpool_transfer(pairs,gappairs);
  }

  return pairs;
}

static List_T
traverse_cdna_gap (bool *filledp, List_T pairs, List_T *path,
		   int leftquerypos, int leftgenomepos, int rightquerypos, int rightgenomepos,
#ifdef PMAP
		   char *queryaaseq_ptr,
#endif
		   char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		   int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogR, 
		   int maxpeelback, int extramaterial_paired, int extraband_paired, double defect_rate) {
  List_T gappairs, peeled_pairs, peeled_path;
  Pair_T leftpair, rightpair;
  int queryjump, genomejump;
  int querydp5 = leftquerypos + 1;
  int genomedp5 = leftgenomepos + 1;
  int querydp3 = rightquerypos - 1;
  int genomedp3 = rightgenomepos - 1;
  int finalscore;

  pairs = peel_back(&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
		    /*intronp*/true);
  *path = peel_forward(&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
		       /*intronp*/true);

  if (peeled_pairs == NULL || peeled_path == NULL) {
    debug(printf("Skipping this because unable to peel\n"));
    *filledp = false;
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    return pairs;
  }

  genomejump = genomedp3 - genomedp5 + 1;
  /* Set queryjump approximately equal to genomejump to have square
     dynamic programming matrices */
  queryjump = genomejump + extramaterial_paired;

  /* Bounds don't make sense */
  /*
  if (querydp5 + queryjump - 1 >= querydp3 - queryjump + 1) {
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    return pairs;
  }
  */

  gappairs = Dynprog_cdna_gap(&finalscore,dynprogL,dynprogR,
			      &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
			      &(queryseq_ptr[querydp3]),&(queryuc_ptr[querydp3]),
			      &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
			      queryjump,queryjump,genomejump,
			      querydp5,querydp3,genomedp5,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      cdna_direction,pairpool,extraband_paired,defect_rate);
  debug(Pair_dump_list(gappairs,true));
  *filledp = true;
  pairs = Pairpool_transfer(pairs,gappairs);
  return pairs;
}


/* genome_gap is usually an intron */
static List_T
traverse_genome_gap (bool *filledp, int *nintrons, int *nnonintrons, int *intronlen, int *nonintronlen, 
		     List_T pairs, List_T *path,
		     int leftquerypos, int leftgenomepos, int rightquerypos, int rightgenomepos,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, Pairpool_T pairpool,
		     Dynprog_T dynprogL, Dynprog_T dynprogR, int maxpeelback,
		     int extramaterial_paired, int extraband_paired, double defect_rate, bool forcep) {
  List_T gappairs, peeled_pairs, peeled_path, micropairs;
  Pair_T leftpair, rightpair;
  int queryjump, genomejump;
  int querydp5 = leftquerypos + 1;
  int genomedp5 = leftgenomepos + 1;
  int querydp3 = rightquerypos - 1;
  int genomedp3 = rightgenomepos - 1;
  int finalscore, nmatches, nmismatches, nopens, nindels, exonhead, introntype, microintrontype;
  int acceptable_nmismatches;

  pairs = peel_back(&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
		    /*intronp*/true);
  *path = peel_forward(&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
		       /*intronp*/true);

  if (peeled_pairs == NULL || peeled_path == NULL) {
    debug(printf("Skipping this because unable to peel\n"));
    *filledp = false;
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    return pairs;
  }

  queryjump = querydp3 - querydp5 + 1;
  /* Set genomejump approximately equal to queryjump to have square
     dynamic programming matrices */
  genomejump = queryjump + extramaterial_paired;

  gappairs = Dynprog_genome_gap(&finalscore,&nmatches,&nmismatches,&nopens,&nindels,
				&exonhead,&introntype,dynprogL,dynprogR,
				&(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				&(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				&(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
				queryjump,genomejump,genomejump,
				querydp5,genomedp5,genomedp3,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				cdna_direction,pairpool,extraband_paired,false,
				defect_rate,/*returnpairsp*/true,/*addgapp*/true,maxpeelback);
  debug(Pair_dump_list(gappairs,true));
  debug(printf("  Score: %d\n",finalscore));

#if 0
  if (defect_rate < DEFECT_HIGHQ) {
    acceptable_nmismatches = 0;
  } else if (defect_rate < DEFECT_MEDQ) {
    acceptable_nmismatches = 2;
  } else {
    acceptable_nmismatches = 3;
  }
#else
  acceptable_nmismatches = 2;
#endif

  debug(printf("nmismatches = %d, nopens = %d, nindels = %d.  acceptable nmismatches = %d\n",
	       nmismatches,nopens,nindels,acceptable_nmismatches));
    
  if (!forcep && finalscore < 0) {
    *filledp = false;
    /* Put back peeled pairs */
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    introntype = NONINTRON;
  } else if (false && defect_rate > DEFECT_MEDQ) {
    /* Don't look for microexons in low-quality sequences */
    debug(printf("Don't look for microexon in low-quality sequence\n"));
    *filledp = true;
    pairs = Pairpool_transfer(pairs,gappairs);
  } else if (introntype != NONINTRON && nmismatches + nindels <= acceptable_nmismatches) {
    *filledp = true;
    pairs = Pairpool_transfer(pairs,gappairs);
  } else {
    *filledp = true;
    debug(printf("Calling microexon because introntype == %d or nmismatches %d + nindels %d > acceptable %d\n",
		 introntype,nmismatches,nindels,acceptable_nmismatches));
    micropairs = Dynprog_microexon_int(&microintrontype,
				       &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				       &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				       &(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
				       queryjump,genomejump,genomejump,
				       querydp5,genomedp5,genomedp3,cdna_direction,
#ifdef PMAP
				       queryaaseq_ptr,
#endif
				       queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,pairpool);
    if (micropairs != NULL) {
      debug(printf("Microexon found\n"));
      pairs = Pairpool_transfer(pairs,micropairs);
      introntype = microintrontype;
    } else {
      debug(printf("No microexon found\n"));
      pairs = Pairpool_transfer(pairs,gappairs);
    }
  }

  if (introntype == NONINTRON) {
    *nnonintrons += 1;
    *nonintronlen += rightgenomepos - leftgenomepos;
  } else {
    *nintrons += 1;
    *intronlen += rightgenomepos - leftgenomepos;
  }

  return pairs;
}


static List_T
traverse_dual_genome_gap (bool *singlep, List_T pairs, List_T *path, 
			  int leftquerypos, int leftgenomepos, int rightquerypos, int rightgenomepos, 
			  int midquerypos, int midgenomepos, 
#ifdef PMAP
			  char *queryaaseq_ptr,
#endif
			  char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
			  int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, 
			  int maxpeelback, int nullgap, int extramaterial_paired, int extraband_paired,
			  double defect_rate) {
  List_T single_gappairs, dual_gappairs_1, dual_gappairs_2, peeled_pairs, peeled_path;
  Pair_T leftpair, rightpair;
  int queryjump, genomejump;
  int querydp5 = leftquerypos + 1;
  int genomedp5 = leftgenomepos + 1;
  int querydp3 = rightquerypos - 1;
  int genomedp3 = rightgenomepos - 1;
  int single_score, dual_score_1, dual_score_2, single_goodness, dual_goodness, 
    nmatches, nmismatches, nopens, nindels, exonhead, right_exonhead, left_exonhead;
  int middle_exonlength, interexon_region;
  int single_introntype, dual_introntype_1, dual_introntype_2, 
    canonical_introntype, semicanonical_introntype_1, semicanonical_introntype_2;
  double middle_exonprob;

  if (cdna_direction > 0) {
    canonical_introntype = GTAG_FWD;
    semicanonical_introntype_1 = ATAC_FWD;
    semicanonical_introntype_2 = GCAG_FWD;
  } else {
    canonical_introntype = GTAG_REV;
    semicanonical_introntype_1 = ATAC_REV;
    semicanonical_introntype_2 = GCAG_REV;
  }

  pairs = peel_back(&peeled_pairs,pairs,&querydp3,&genomedp3,pairpool,maxpeelback,
		    /*intronp*/false);
  *path = peel_forward(&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback,
		       /*intronp*/false);

  queryjump = querydp3 - querydp5 + 1;
  genomejump = queryjump + extramaterial_paired;

  if (queryjump > nullgap) {
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    *singlep = false;
    return pairs;
  }

  /* Bounds don't make sense */
  if (genomedp5 + genomejump - 1 >= genomedp3 - genomejump + 1) {
    debug(printf("Bounds don't make sense for dual intron gap\n\n"));
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    *path = Pairpool_transfer(*path,peeled_path);
    *singlep = false;
    return pairs;
  }
  
  single_gappairs = Dynprog_genome_gap(&single_score,&nmatches,&nmismatches,&nopens,&nindels,
				       &exonhead,&single_introntype,dynprogL,dynprogR,
				       &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				       &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				       &(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
				       queryjump,genomejump,genomejump,
				       querydp5,genomedp5,genomedp3,
#ifdef PMAP
				       queryaaseq_ptr,
#endif
				       cdna_direction,pairpool,extraband_paired,false,
				       defect_rate,/*returnpairsp*/true,/*addgapp*/false,maxpeelback);

  debug(Pair_check_list(single_gappairs));

  single_goodness = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;
  if (single_introntype == canonical_introntype) {
    single_goodness += CANONICAL_POINTS;
  } else if (single_introntype == semicanonical_introntype_1 ||
	     single_introntype == semicanonical_introntype_2) {
    single_goodness += SEMICANONICAL_POINTS;
  }

  /* Right of short exon */
  querydp5 = midquerypos;
  genomedp5 = midgenomepos;
  querydp3 = rightquerypos;
  genomedp3 = rightgenomepos;

  queryjump = querydp3 - querydp5 + 1;
  genomejump = queryjump + extramaterial_paired;

  dual_gappairs_2 = Dynprog_genome_gap(&dual_score_2,&nmatches,&nmismatches,&nopens,&nindels,
				       &right_exonhead,&dual_introntype_2,dynprogL,dynprogR,
				       &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				       &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				       &(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
				       queryjump,genomejump,genomejump,
				       querydp5,genomedp5,genomedp3,
#ifdef PMAP
				       queryaaseq_ptr,
#endif
				       cdna_direction,pairpool,extraband_paired,false,
				       defect_rate,/*returnpairsp*/false,/*addgapp*/false,maxpeelback);

  dual_goodness = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

  /* Left of short exon */
  querydp5 = leftquerypos;
  genomedp5 = leftgenomepos;
  querydp3 = midquerypos;
  genomedp3 = midgenomepos;

  queryjump = querydp3 - querydp5 + 1;
  genomejump = queryjump + extramaterial_paired;

  dual_gappairs_1 = Dynprog_genome_gap(&dual_score_1,&nmatches,&nmismatches,&nopens,&nindels,
				       &left_exonhead,&dual_introntype_1,dynprogL,dynprogR,
				       &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
				       &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
				       &(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
				       queryjump,genomejump,genomejump,
				       querydp5,genomedp5,genomedp3,
#ifdef PMAP
				       queryaaseq_ptr,
#endif
				       cdna_direction,pairpool,extraband_paired,false,
				       defect_rate,/*returnpairsp*/false,/*addgapp*/false,maxpeelback);

  dual_goodness += nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;
  if (dual_introntype_1 == canonical_introntype && dual_introntype_1 == canonical_introntype) {
    dual_goodness += CANONICAL_POINTS;
  }
  /* Cost of extra intron */
  if (defect_rate < DEFECT_HIGHQ) {
    /* dual_goodness -= 0; */
  } else if (defect_rate < DEFECT_MEDQ) {
    dual_goodness -= 8;
  } else {
    dual_goodness -= 16;		
  }

  middle_exonlength = right_exonhead-left_exonhead;
  debug(printf("Middle exon is %d - %d = %d bp in interexon region of %d bp\n",
	       right_exonhead,left_exonhead,right_exonhead-left_exonhead,rightgenomepos-leftgenomepos));
  if (middle_exonlength <= 0) {
    middle_exonprob = 0.0;
  } else {
    interexon_region = rightgenomepos-leftgenomepos;

    if (dual_introntype_2 == canonical_introntype) {
      middle_exonlength += DUAL_HALFCANONICAL_POINTS;
      debug(printf("Add canonical credit of %d for right intron\n",DUAL_HALFCANONICAL_POINTS));
    }
    if (dual_introntype_1 == canonical_introntype) {
      middle_exonlength += DUAL_HALFCANONICAL_POINTS;
      debug(printf("Add canonical credit of %d for left intron\n",DUAL_HALFCANONICAL_POINTS));
    }

    middle_exonprob = 1.0-pow(1.0-pow(4.0,-(double) middle_exonlength),(double) interexon_region);

    debug(printf("Single score = %d.  Dual score = %d & %d.  ",single_score,dual_score_1,dual_score_2));
    debug(printf("Single goodness = %d.  Dual goodness = %d.  ",
		 single_goodness,dual_goodness));
    debug(printf("Probability is %g.  ",middle_exonprob));
  }

  if (middle_exonprob > 0.001 || single_goodness >= dual_goodness) {
    debug(printf(  "Single score wins\n"));
    pairs = Pairpool_transfer(pairs,single_gappairs);
    *singlep = true;

  } else {
    debug(printf(  "Dual scores win\n"));

    /* Put back peeled pairs */
    debug(printf("Transferring back onto pairs:\n"));
    pairs = Pairpool_transfer(pairs,peeled_pairs);
    debug(printf("Transferring back onto path:\n"));
    *path = Pairpool_transfer(*path,peeled_path);
    debug(printf("Done with transfers.\n"));
    *singlep = false;
  }

  return pairs;
}


static List_T
try_ending5_extend (int *finalscore, List_T *pairs, 
		    int leftquerypos, int leftgenomepos, int rightquerypos, int rightgenomepos,
#ifdef PMAP
		    char *queryaaseq_ptr,
#endif
		    char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		    int cdna_direction, int introndir, Pairpool_T pairpool,
		    Dynprog_T dynprog, int maxpeelback, int extramaterial_end,
		    int extraband_end, double defect_rate, bool extend_mismatch_p,
		    bool end_microexons_p) {
  Pair_T lastpair, rightpair;
  List_T intron_gappairs = NULL, continuous_gappairs = NULL, peeled_pairs, p;
  int queryjump, genomejump;
  int querydp5 = leftquerypos + 1;
  int genomedp5 = leftgenomepos + 1;
  int querydp3 = rightquerypos - 1;
  int genomedp3 = rightgenomepos - 1;
  int microintrontype, microexonlength;
  int intron_goodness, continuous_goodness, nmissed, nmatches, nmismatches, nopens, nindels, acceptable_nmismatches;

  /* Note that we peelback only half as much as for a paired gap, to
     save on dynamic programming */
  *pairs = peel_back(&peeled_pairs,*pairs,&querydp3,&genomedp3,pairpool,maxpeelback/2,
		     /*intronp*/false);
  
  queryjump = querydp3 - querydp5 + 1;
  genomejump = queryjump + extramaterial_end; /* proposed */
  /* Previously, we limited genomejump = min(2*queryjump,queryjump+extramaterial_end) */

  genomedp5 = genomedp3 - genomejump + 1;
  /* Make sure we don't go past the beginning */
  if (genomedp5 < 0) {
    genomedp5 = 0;
    genomejump = genomedp3 - genomedp5 + 1;
  }
  debug(printf("Stage 3: Dynamic programming at 5' end: querydp5 = %d, querydp3 = %d, genomedp5 = %d, genomedp3 = %d\n",
	       querydp5,querydp3,genomedp5,genomedp3));

  debug(printf("Trying to extend 5' end:\n"));
  continuous_gappairs = Dynprog_end5_gap(&(*finalscore),&nmatches,&nmismatches,&nopens,&nindels,dynprog,
					 &(queryseq_ptr[querydp3]),&(queryuc_ptr[querydp3]),
					 &(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
					 queryjump,genomejump,querydp3,genomedp3,
#ifdef PMAP
					 queryaaseq_ptr,
#endif
					 cdna_direction,pairpool,extraband_end,defect_rate,extend_mismatch_p);
  continuous_goodness = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

  /* Try to find microexons on 5' end */
  if (continuous_gappairs == NULL) {
    nmissed = querydp3;
    debug(printf("  => %d missed\n",nmissed));
  } else {
    for (p = continuous_gappairs; p != NULL; p = p->rest) {
      lastpair = p->first;
    }
    nmissed = lastpair->querypos;
    debug(printf("  => %d missed and %d mismatches\n",nmissed,nmismatches));
  }

  acceptable_nmismatches = 2;

  if (defect_rate > DEFECT_HIGHQ || introndir != cdna_direction || nmissed + nmismatches <= acceptable_nmismatches) {
    debug(printf("Continuous is acceptable\n"));
    return continuous_gappairs;
  } else {
    /* If we are sure about the intron direction, try to find microexon */
    debug(printf("Trying to find microexon at 5' end:\n"));
    intron_gappairs = Dynprog_microexon_5(&microintrontype,&microexonlength,
					  &(queryseq_ptr[querydp3]),&(queryuc_ptr[querydp3]),
					  &(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
					  queryjump,genomejump,querydp3,genomedp3,cdna_direction,
#ifdef PMAP
					  queryaaseq_ptr,
#endif
					  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
					  pairpool,end_microexons_p);
    if (intron_gappairs == NULL) {
      debug(printf("Continuous wins by default: %d\n",continuous_goodness));
      return continuous_gappairs;
    } else {
      intron_goodness = microexonlength;
      if (intron_goodness > continuous_goodness + acceptable_nmismatches) {
	debug(printf("Microexon wins: intron %d, continuous %d\n",intron_goodness,continuous_goodness));
	return intron_gappairs;
      } else {
	debug(printf("Continuous wins: intron %d, continuous %d\n",intron_goodness,continuous_goodness));
	return continuous_gappairs;
      }
    }
  }
}


static List_T
try_ending5_extend_simple (int rightquerypos, int rightgenomepos, int genomiclength,
#ifdef PMAP
			   char *queryaaseq_ptr,
#endif
			   char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
			   int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprog, int extramaterial_end,
			   int extraband_end, bool extend_mismatch_p) {
  int queryjump, genomejump;
  int querydp5 = 0;
  int genomedp5 = 0;
  int querydp3 = rightquerypos - 1;
  int genomedp3 = genomiclength - 1;
  int nmatches, nmismatches, nopens, nindels;
  int finalscore;

  queryjump = querydp3 - querydp5 + 1;
  genomejump = queryjump + extramaterial_end; /* proposed */

  genomedp5 = genomedp3 - genomejump + 1;
  /* Make sure we don't go past the beginning */
  if (genomedp5 < 0) {
    genomedp5 = 0;
    genomejump = genomedp3 - genomedp5 + 1;
  }
  debug(printf("Stage 3: Dynamic programming at 5' end: querydp5 = %d, querydp3 = %d, genomedp5 = %d, genomedp3 = %d\n",
	       querydp5,querydp3,genomedp5,genomedp3));

  debug(printf("Trying to extend 5' end:\n"));
  return Dynprog_end5_gap(&finalscore,&nmatches,&nmismatches,&nopens,&nindels,dynprog,
			  &(queryseq_ptr[querydp3]),&(queryuc_ptr[querydp3]),
			  &(genomicseg_ptr[genomedp3]),&(genomicuc_ptr[genomedp3]),
			  queryjump,genomejump,querydp3,rightgenomepos-1,
#ifdef PMAP
			  queryaaseq_ptr,
#endif
			  cdna_direction,pairpool,extraband_end,/*defect_rate*/0.0,extend_mismatch_p);
}

static List_T
try_ending3_extend (int *finalscore, List_T *path, 
		    int leftquerypos, int leftgenomepos, int rightquerypos, int rightgenomepos, 
		    int querylength, Genomicpos_T genomiclength,
#ifdef PMAP
		    char *queryaaseq_ptr,
#endif
		    char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		    int cdna_direction, int introndir,
		    Pairpool_T pairpool, Dynprog_T dynprog, int maxpeelback, int extramaterial_end,
		    int extraband_end, double defect_rate, bool extend_mismatch_p,
		    bool end_microexons_p) {
  Pair_T firstpair, leftpair;
  List_T intron_gappairs = NULL, continuous_gappairs = NULL, peeled_path;
  int queryjump, genomejump;
  int querydp5 = leftquerypos + 1;
  int genomedp5 = leftgenomepos + 1;
  int querydp3 = rightquerypos - 1;
  int genomedp3 = rightgenomepos - 1;
  int microintrontype, microexonlength;
  int intron_goodness, continuous_goodness, nmissed, nmatches, nmismatches, nopens, nindels, acceptable_nmismatches;

  /* Note that we peelback only half as much as for a paired gap, to
     save on dynamic programming */
  *path = peel_forward(&peeled_path,*path,&querydp5,&genomedp5,pairpool,maxpeelback/2,
		       /*intronp*/false);

  queryjump = querydp3 - querydp5 + 1;
  genomejump = queryjump + extramaterial_end; /* proposed */
  /* Previously, we limited genomejump = min(2*queryjump,queryjump+extramaterial_end) */

  genomedp3 = genomedp5 + genomejump - 1;
  /* Make sure we don't go past the end */
  if (genomedp3 > genomiclength - 1) {
    genomedp3 = genomiclength - 1;
    genomejump = genomedp3 - genomedp5 + 1;
  }
  debug(printf("Stage 3: Dynamic programming at 3' end: querydp5 = %d, querydp3 = %d, genomedp5 = %d, genomedp3 = %d\n",
	       querydp5,querydp3,genomedp5,genomedp3));

  debug(printf("Trying to extend 3' end:\n"));
  continuous_gappairs = Dynprog_end3_gap(&(*finalscore),&nmatches,&nmismatches,&nopens,&nindels,dynprog,
					 &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
					 &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
					 queryjump,genomejump,querydp5,genomedp5,
#ifdef PMAP
					 queryaaseq_ptr,
#endif
					 cdna_direction,pairpool,extraband_end,defect_rate,extend_mismatch_p);
  continuous_goodness = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

  /* Try to find microexons on 3' end */
  if (continuous_gappairs == NULL) {
    nmissed = querylength - 1 - querydp5;
    debug(printf("  => %d missed\n",nmissed));
  } else {
    firstpair = continuous_gappairs->first;
    nmissed = querylength - 1 - firstpair->querypos;
    debug(printf("  => %d missed and %d mismatches\n",nmissed,nmismatches));
  }

  acceptable_nmismatches = 2;
  
  if (defect_rate > DEFECT_HIGHQ || introndir != cdna_direction || nmissed + nmismatches <= acceptable_nmismatches) {
    debug(printf("Continuous is acceptable\n"));
    return List_reverse(continuous_gappairs);
  } else {
    /* If we are sure about the intron direction, try to find microexon. */
    debug(printf("Trying to find microexon at 3' end:\n"));
    intron_gappairs = Dynprog_microexon_3(&microintrontype,&microexonlength,
					  &(queryseq_ptr[querydp5]),&(queryuc_ptr[querydp5]),
					  &(genomicseg_ptr[genomedp5]),&(genomicuc_ptr[genomedp5]),
					  queryjump,genomejump,querydp5,genomedp5,cdna_direction,
#ifdef PMAP
					  queryaaseq_ptr,
#endif					  
					  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
					  genomiclength,pairpool,end_microexons_p);
    if (intron_gappairs == NULL) {
      debug(printf("Continuous wins by default: %d\n",continuous_goodness));
      return List_reverse(continuous_gappairs);
    } else {
      intron_goodness = microexonlength;
      if (intron_goodness > continuous_goodness + acceptable_nmismatches) {
	debug(printf("Microexon wins: intron %d, continuous %d\n",intron_goodness,continuous_goodness));
	return List_reverse(intron_gappairs);
      } else {
	debug(printf("Continuous wins: intron %d, continuous %d\n",intron_goodness,continuous_goodness));
	return List_reverse(continuous_gappairs);
      }
    }
  }
}

static List_T
try_ending3_extend_simple (int leftquerypos, int querylength, int leftgenomepos, Genomicpos_T genomiclength,
#ifdef PMAP
			   char *queryaaseq_ptr,
#endif
			   char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
			   int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprog, int extramaterial_end,
			   int extraband_end, bool extend_mismatch_p) {
  int queryjump, genomejump;
  int querydp5 = leftquerypos + 1;
  int genomedp5 = 0;
  int querydp3 = querylength - 1;
  int genomedp3 = genomiclength - 1;
  int nmatches, nmismatches, nopens, nindels, finalscore;

  queryjump = querydp3 + 1;
  genomejump = queryjump + extramaterial_end; /* proposed */

  genomedp3 = genomejump - 1;
  /* Make sure we don't go past the end */
  if (genomedp3 > genomiclength - 1) {
    genomedp3 = genomiclength - 1;
    genomejump = genomedp3 - genomedp5 + 1;
  }
  debug(printf("Stage 3: Dynamic programming at 3' end: querydp5 = %d, querydp3 = %d, genomedp5 = %d, genomedp3 = %d\n",
	       querydp5,querydp3,genomedp5,genomedp3));

  debug(printf("Trying to extend 3' end:\n"));
  return Dynprog_end3_gap(&finalscore,&nmatches,&nmismatches,&nopens,&nindels,dynprog,
			  /*sequence1*/&(queryseq_ptr[querydp5]),/*sequenceuc1*/&(queryuc_ptr[querydp5]),
			  /*sequence2*/genomicseg_ptr,/*sequenceuc2*/genomicuc_ptr,
			  /*length1*/queryjump,/*length2*/genomejump,/*offset1*/querydp5,/*offset2*/leftgenomepos+1,
#ifdef PMAP
			  queryaaseq_ptr,
#endif
			  cdna_direction,pairpool,extraband_end,/*defect_rate*/0.0,extend_mismatch_p);
}

  
/* Note: querypos is actually indexsize nt to the left of the last nt match.
   
        ||||||||********   X  X  XX X   X
               ^         <- queryjump->  ^
           querypos                      lastquerypos         

                <-     querydpspan     ->
*/
static List_T
build_path_end3 (List_T path, int querylength, Genomicpos_T genomiclength,
#ifdef PMAP
		 char *queryaaseq_ptr,
#endif
		 char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		 int cdna_direction, int introndir, int maxpeelback, int nullgap,
		 int extramaterial_end, int extramaterial_paired, int extraband_end,
		  double defect_rate, Pairpool_T pairpool, Dynprog_T dynprogL, bool extend_mismatch_p,
		  bool end_microexons_p) {
  List_T gappairs;
  Pair_T leftpair;
  int queryjump, genomejump, rightquerypos;
  int finalscore;

  leftpair = path->first;
  debug(printf("Stage 3: 3' end: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d\n",
	       leftpair->querypos,querylength,leftpair->genomepos,genomiclength));

  queryjump = querylength - leftpair->querypos - 1;
  genomejump = genomiclength - leftpair->genomepos - 1;

  /* Note difference with 5' case.  We use queryjump+1 here instead of queryjump and genomejump */
  if (queryjump+1 > nullgap) {
    rightquerypos = leftpair->querypos + nullgap + 1;
  } else {
    rightquerypos = querylength;
  }
  gappairs = try_ending3_extend(&finalscore,&path,
				leftpair->querypos,leftpair->genomepos,rightquerypos,genomiclength,
				querylength,genomiclength,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				cdna_direction,introndir,pairpool,dynprogL,maxpeelback,
				extramaterial_end,extraband_end,defect_rate,extend_mismatch_p,
				end_microexons_p);
  debug(Pair_dump_list(gappairs,true));
  path = Pairpool_transfer(path,gappairs);

  return path;
}


/* Schematic:
   
       <- queryjump ->
       X  X  XX X   X ********||||||||
      ^               ^
   querypos        lastquerypos

       <-    querydpspan    ->

*/
static List_T
build_pairs_end5 (List_T pairs, 
#ifdef PMAP
		  char *queryaaseq_ptr,
#endif
		  char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		  int cdna_direction, int introndir, int maxpeelback, int nullgap,
		  int extramaterial_end, int extramaterial_paired, int extraband_end,
		  double defect_rate, Pairpool_T pairpool, Dynprog_T dynprogR, bool extend_mismatch_p,
		  bool end_microexons_p) {
  List_T gappairs;
  Pair_T rightpair;
  int queryjump, genomejump, leftquerypos;
  int finalscore;

  if (pairs == NULL) {
    return NULL;
  } else {
    rightpair = pairs->first;
  }
  debug(printf("Stage 3: 5' end: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d\n",
	       -1,rightpair->querypos,-1,rightpair->genomepos));

  queryjump = rightpair->querypos; /* - leftquerypos (-1) - 1 */
  genomejump = rightpair->genomepos; /* - leftgenomepos (-1) - 1 */

  /* Note difference with 3' case.  We use queryjump here instead of queryjump+1 */
  if (queryjump > nullgap) {
    leftquerypos = rightpair->querypos - nullgap - 1;
  } else {
    leftquerypos = -1;
  }

  gappairs = try_ending5_extend(&finalscore,&pairs,leftquerypos,/*leftgenomepos*/-1,
				rightpair->querypos,rightpair->genomepos,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				cdna_direction,introndir,pairpool,dynprogR,maxpeelback,
				extramaterial_end,extraband_end,defect_rate,extend_mismatch_p,
				end_microexons_p);
  debug(Pair_dump_list(gappairs,true));
  pairs = Pairpool_transfer(pairs,gappairs);
  return pairs;
}


/************************************************************************
 *   Gaps
 ************************************************************************/

static List_T
check_gaps (List_T pairs, Pairpool_T pairpool) {
  List_T path = NULL;
  Pair_T pair, firstpair, leftpair, rightpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump;

  debug(printf("\nBeginning check of gaps\n"));

  pairs = Pairpool_pop(pairs,&firstpair);
  if (firstpair->gapp == true) {
    fprintf(stderr,"Unexpected gap at start of pairs\n");
    debug(printf("Unexpected gap at start of pairs\n"));
#ifndef DEBUG
    Except_raise(&gapcheck_error,__FILE__,__LINE__);
    abort();
#endif
  } else {
    path = Pairpool_push_existing(NULL,pairpool,firstpair);
    lastquerypos = firstpair->querypos;
    lastgenomepos = firstpair->genomepos;
  }

  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == true) {
      debug(printf("Observed a gap at %d..%d with queryjump = %d, genomejump = %d\n",
		   lastquerypos,querypos,pair->queryjump,pair->genomejump));
      leftpair = path->first;
      rightpair = pairs->first;
      queryjump = rightpair->querypos - leftpair->querypos - 1;
      genomejump = rightpair->genomepos - leftpair->genomepos - 1;

      if (pair->queryjump != queryjump) {
	if (rightpair->querypos >= HALFLEN && leftpair->querypos < HALFLEN) {
	  debug(printf("Accept queryjump for gap at %d..%d as probable skiplength.  It's %d, should be %d\n",
		       lastquerypos,querypos,pair->queryjump,queryjump));
	} else {
	  debug(printf("Wrong queryjump for gap at %d..%d.  It's %d, should be %d\n",
		       lastquerypos,querypos,pair->queryjump,queryjump));
#ifndef DEBUG
	  Except_raise(&gapcheck_error,__FILE__,__LINE__);
	  abort();
#endif
	}
      }
      if (pair->genomejump != genomejump) {
	debug(printf("Wrong genomejump for gap at %d..%d.  It's %d, should be %d\n",
		     lastquerypos,querypos,pair->genomejump,genomejump));
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
	abort();
#endif
      }
      path = Pairpool_push_existing(path,pairpool,pair);

      /* Process another pair */
      if (pairs == NULL) {
	fprintf(stderr,"Unexpected gap at end of pairs\n");
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
	abort();
#endif
      }
      pairs = Pairpool_pop(pairs,&pair);
      if (pair->gapp == true) {
	fprintf(stderr,"Unexpected gap after gap\n");
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
	abort();
#endif
      }
      path = Pairpool_push_existing(path,pairpool,pair);

    } else {
      querypos = pair->querypos;
      genomepos = pair->genomepos;

      queryjump = querypos - lastquerypos - 1;
      genomejump = genomepos - lastgenomepos - 1;
	  
      if (queryjump <= 0 && genomejump <= 0) {
	path = Pairpool_push_existing(path,pairpool,pair);
      } else if (queryjump == 0 && genomejump == 0) {
	path = Pairpool_push_existing(path,pairpool,pair);
      } else {
	fprintf(stderr,"Unexpected missing gap at %d..%d\n",lastquerypos,querypos);
	debug(printf("Unexpected missing gap at %d..%d\n",lastquerypos,querypos));
	debug(printf("Pushing a gap at %d..%d because of queryjump = %d, genomejump = %d\n",
		     lastquerypos,querypos,queryjump,genomejump));
	path = Pairpool_push_gap(path,pairpool,queryjump,genomejump);
	path = Pairpool_push_existing(path,pairpool,pair);
#ifndef DEBUG
	Except_raise(&gapcheck_error,__FILE__,__LINE__);
	abort();
#endif
      }
    }

    lastquerypos = pair->querypos;
    lastgenomepos = pair->genomepos;
  }	

  debug(printf("Done with check of gaps\n\n"));

  return path;
}

static List_T
insert_gaps (List_T pairs, Pairpool_T pairpool) {
  List_T path = NULL;
  Pair_T pair, firstpair, leftpair, rightpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump;

  debug(printf("\nBeginning insertion of gaps\n"));

  pairs = Pairpool_pop(pairs,&firstpair);
  if (firstpair->gapp == true) {
    /* Discard gap */
  } else {
    path = Pairpool_push_existing(NULL,pairpool,firstpair);
    lastquerypos = firstpair->querypos;
    lastgenomepos = firstpair->genomepos;
  }

  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == true) {
      debug(printf("Observed a gap at %d..%d with queryjump = %d, genomejump = %d\n",
		   lastquerypos,querypos,pair->queryjump,pair->genomejump));
      leftpair = path->first;
      rightpair = pairs->first;
      queryjump = rightpair->querypos - leftpair->querypos - 1;
      genomejump = rightpair->genomepos - leftpair->genomepos - 1;

      if (pair->queryjump != queryjump) {
	debug(printf("Fixing queryjump for gap at %d..%d from %d to %d\n",
		     lastquerypos,querypos,pair->queryjump,queryjump));
	pair->queryjump = queryjump;
      }
      if (pair->genomejump != genomejump) {
	debug(printf("Fixing genomejump for gap at %d..%d from %d to %d\n",
		     lastquerypos,querypos,pair->genomejump,genomejump));
	pair->genomejump = genomejump;
      }
      path = Pairpool_push_existing(path,pairpool,pair);

      /* Process another pair */
      if (pairs == NULL) {
	fprintf(stderr,"Unexpected gap at end of pairs\n");
	abort();
      }
      pairs = Pairpool_pop(pairs,&pair);
      if (pair->gapp == true) {
	fprintf(stderr,"Unexpected gap after gap\n");
	abort();
      }
      path = Pairpool_push_existing(path,pairpool,pair);

    } else {
      querypos = pair->querypos;
      genomepos = pair->genomepos;

      queryjump = querypos - lastquerypos - 1;
      genomejump = genomepos - lastgenomepos - 1;
	  
      if (queryjump <= 0 && genomejump <= 0) {
	path = Pairpool_push_existing(path,pairpool,pair);
      } else {
	debug(printf("Pushing a gap at %d..%d because of queryjump = %d, genomejump = %d\n",
		     lastquerypos,querypos,queryjump,genomejump));
	path = Pairpool_push_gap(path,pairpool,queryjump,genomejump);
	path = Pairpool_push_existing(path,pairpool,pair);
      }
    }

    lastquerypos = pair->querypos;
    lastgenomepos = pair->genomepos;
  }	

  debug(printf("Done with insertion of gaps\n\n"));

  return path;
}

static List_T
insert_gaps_and_singles (List_T pairs, Pairpool_T pairpool,
			 char *queryseq_ptr, char *queryuc_ptr,
			 char *genomicseg_ptr, char *genomicuc_ptr,
			 int skiplength) {
  List_T path = NULL;
  Pair_T pair, firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump;
  char c1, c2;

  debug(printf("\nBeginning insertion of gaps and singles\n"));

  pairs = Pairpool_pop(pairs,&firstpair);
  path = Pairpool_push_existing(NULL,pairpool,firstpair);
  lastquerypos = firstpair->querypos;
  lastgenomepos = firstpair->genomepos;

  while (pairs != NULL) {
    pairs = Pairpool_pop(pairs,&pair);
    if (pair->gapp == true) {
      debug(printf("Observed a gap at %d..%d with queryjump = %d, genomejump = %d\n",
		   lastquerypos,querypos,queryjump,genomejump));
      path = Pairpool_push_existing(path,pairpool,pair);
      if (pairs == NULL) {
	fprintf(stderr,"Unexpected gap at end of pairs\n");
	abort();
      }
      pairs = Pairpool_pop(pairs,&pair);
      if (pair->gapp == true) {
	fprintf(stderr,"Unexpected gap after gap\n");
	abort();
      }
      path = Pairpool_push_existing(path,pairpool,pair);
    } else {
      querypos = pair->querypos;
      genomepos = pair->genomepos;

      queryjump = querypos - lastquerypos - 1;
      genomejump = genomepos - lastgenomepos - 1;
      if (skiplength > 0 && querypos >= HALFLEN && lastquerypos < HALFLEN) {
	debug(printf("Introducing skiplength %d at %d..%d\n",skiplength,lastquerypos,querypos));
	queryjump += skiplength;
      }
	  
      if (queryjump < 0 || genomejump < 0) {
	fprintf(stderr,"Unexpected negative queryjump or genomejump\n");
	abort();
      } else if (queryjump == 0 && genomejump == 0) {
	path = Pairpool_push_existing(path,pairpool,pair);
      } else if (queryjump == 1 && genomejump == 1) {
	/* Special case: Single mismatch */
	c1 = queryseq_ptr[lastquerypos+1];
	c2 = genomicseg_ptr[lastgenomepos+1];
	if (queryuc_ptr[lastquerypos+1] == genomicuc_ptr[lastgenomepos+1]) {
	  /* Possible for a stretch of the same nucleotide to have a mismatch */
	  debug(printf("Stage 3: Single mismatch at %d %d: %c | %c\n",
		       lastquerypos+1,lastgenomepos+1,c1,c2));
	  path = Pairpool_push(path,pairpool,lastquerypos+1,lastgenomepos+1,
			       c1,MATCH_COMP,c2,/*gapp*/false);
	  path = Pairpool_push_existing(path,pairpool,pair);
	} else {
	  debug(printf("Stage 3: Single mismatch at %d %d: %c   %c\n",
		       lastquerypos+1,lastgenomepos+1,c1,c2));
	  path = Pairpool_push(path,pairpool,lastquerypos+1,lastgenomepos+1,
			       c1,MISMATCH_COMP,c2,/*gapp*/false);
	  path = Pairpool_push_existing(path,pairpool,pair);
	}
      } else if (queryjump == 1 && genomejump == 0) {
	/* Special case: Single cDNA insert.  This will create an apparent genomejump of -1. */
	debug(printf("Stage 3: Single cDNA insert at %d %d: %c\n",
		     lastquerypos+1,lastgenomepos,queryseq_ptr[lastquerypos+1]));
	path = Pairpool_push(path,pairpool,lastquerypos+1,lastgenomepos,
			     queryseq_ptr[lastquerypos+1],INDEL_COMP,' ',/*gapp*/false);
	path = Pairpool_push_existing(path,pairpool,pair);

#if 0
      } else if (queryjump == 0 && genomejump == 1) {
	/* Special case: Single genomic insert */
	debug(printf("Stage 3: Single genomic insert at %d %d: %c\n",
		     lastquerypos,lastgenomepos+1,genomicseg_ptr[lastgenomepos+1]));
	path = Pairpool_push(path,pairpool,lastquerypos,lastgenomepos+1,
			     ' ',INDEL_COMP,genomicseg_ptr[lastgenomepos+1],/*gapp*/false);
	path = Pairpool_push_existing(path,pairpool,pair);
#endif

      } else {
	debug(printf("Pushing a gap at %d..%d because of queryjump = %d, genomejump = %d\n",
		     lastquerypos,querypos,queryjump,genomejump));
	path = Pairpool_push_gap(path,pairpool,queryjump,genomejump);
	path = Pairpool_push_existing(path,pairpool,pair);
      }
    }

    lastquerypos = querypos;
    lastgenomepos = genomepos;
  }	

  debug(printf("Done with insertion of gaps and singles\n\n"));

  return path;
}



#ifdef PMAP
static List_T
undefine_nucleotides (char *queryseq_ptr, int querylength, List_T path, Pairpool_T pairpool, int width) {
  List_T pairs;
  Pair_T pair, firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump, pos;

  path = Pairpool_pop(path,&firstpair);
  pairs = Pairpool_push_existing(NULL,pairpool,firstpair);
  lastquerypos = firstpair->querypos;
  lastgenomepos = firstpair->genomepos;

  debug(printf("\n** Starting undefine_nucleotides\n"));
  while (path != NULL) {
    pair = path->first;
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - querypos - 1;
    genomejump = lastgenomepos - genomepos - 1;
    if (queryjump == 0 && genomejump == 0) {
      /* Do nothing */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    } else {
      debug(printf("Undefining around lastquerypos = %d and querypos = %d\n",lastquerypos,querypos));
      for (pos = lastquerypos; pos < lastquerypos + width && pos < querylength; pos++) {
	queryseq_ptr[pos] = BACKTRANSLATE_CHAR;
      }
      for (pos = querypos; pos > querypos - width && pos >= 0; --pos) {
	queryseq_ptr[pos] = BACKTRANSLATE_CHAR;
      }

      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;
    }
  }

  return pairs;
}  
#endif



static List_T
build_pairs_singles (List_T path, 
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, int maxpeelback, int nullgap,
		     int extramaterial_paired, int extraband_single, double defect_rate,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;
  bool filledp;

  debug(printf("\n** Starting build_pairs_singles\n"));
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

    } else if (pair->genomejump > pair->queryjump + MININTRONLEN) {
      /* Intron.  Do nothing */
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else {
      /* Guarantees: queryjump <= nullgap && genomejump < queryjump - EXTRAQUERYGAP &&
	 genomejump <= queryjump + MININTRONLEN, meaning that score matrix is nearly square */
      leftpair = path->first;
      rightpair = pairs->first;
	
      debug(printf("Stage 3: Traversing single gap: leftquerypos = %d, rightquerypos = %d, leftgenomepos = %d, rightgenomepos = %d\n",
		   leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos));
      pairs = traverse_single_gap(&filledp,pairs,&path,leftpair->querypos,leftpair->genomepos,rightpair->querypos,rightpair->genomepos,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				  cdna_direction,pairpool,dynprogM,maxpeelback,extraband_single,defect_rate,
				  /*forcep*/true);
      /* forcep needs to be true her to avoid subsequent anomalies in building dualintrons, e.g., XM_376610.2_mRNA on 7:127885572..127888991 */
      if (filledp == true) {
	/* Discard the gap */
      } else {
	/* Replace the gap */
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }
    }
  }

  return pairs;
}


static List_T
build_pairs_dualintrons (bool *singlep, List_T path,
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
	  
	  debug(printf("Stage 3: Traversing dual intron gap: leftquerypos = %d, midquerypos = %d, rightquerypos = %d, leftgenomepos = %d, midgenomepos = %d, rightgenomepos = %d\n",
		       leftpair->querypos,midquerypos,rightpair->querypos,
		       leftpair->genomepos,midgenomepos,rightpair->genomepos));
	  
	  pairs = traverse_dual_genome_gap(&(*singlep),pairs,&path,
					   leftpair->querypos,leftpair->genomepos,
					   rightpair->querypos,rightpair->genomepos,
					   midquerypos,midgenomepos,
#ifdef PMAP
					   queryaaseq_ptr,
#endif
					   queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
					   cdna_direction,pairpool,dynprogL,dynprogM,dynprogR,
					   maxpeelback,nullgap,extramaterial_paired,extraband_paired,
					   defect_rate);
	  if (*singlep == true) {
	    /* Single intron solved and put onto pairs.  Discard
	       initial gap (pair) and second gap (midpair) */

	  } else {
	    /* Dual intron kept; continue after first intron */
	    debug(printf("Beginning transfer of middle exon to path:\n"));
	    path = Pairpool_push_existing(path,pairpool,midpair); /* second gap */
	    path = Pairpool_transfer(path,midexon_pairs);
	    debug(printf("Done with transfer of middle exon to path\n"));
	    pairs = Pairpool_push_existing(pairs,pairpool,pair); /* initial gap */
	  }
	}
      }
    }
  }

  return pairs;
}  


static List_T
build_pairs_introns (int *nintrons, int *nnonintrons, int *intronlen, int *nonintronlen, List_T path, 
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, int maxpeelback, int nullgap, int extramaterial_paired, 
		     int extraband_single, int extraband_paired, double defect_rate,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T pairs = NULL;
  Pair_T pair, leftpair, rightpair;
  bool filledp;

  debug(printf("\n** Starting build_pairs_introns\n"));
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == false) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->queryjump > nullgap) {
      debug(leftpair = path->first;
	    rightpair = pairs->first;
	    printf("Stage 3: Adding large gap: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		   leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos)); 
     pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else if (pair->queryjump > pair->genomejump + EXTRAQUERYGAP) {
      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Stage 3: Traversing cDNA gap: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		   leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos));
      pairs = traverse_cdna_gap(&filledp,pairs,&path,leftpair->querypos,leftpair->genomepos,
				rightpair->querypos,rightpair->genomepos,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				pairpool,dynprogL,dynprogR,maxpeelback,extramaterial_paired,
				extraband_paired,defect_rate);

      if (filledp == true) {
	/* Discard gap */
      } else {
	/* Replace the gap */
	debug(leftpair = path->first;
	      rightpair = pairs->first;
	      printf("Stage 3: Adding large gap instead: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		     leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos));
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }

    } else if (pair->genomejump > pair->queryjump + 2*MININTRONLEN) {
      /* Need space for two introns */
      /* We will make the score matrices nearly square */
      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Stage 3: Traversing paired gap: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		   leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos));
      pairs = traverse_genome_gap(&filledp,&(*nintrons),&(*nnonintrons),&(*intronlen),&(*nonintronlen),
				  pairs,&path,leftpair->querypos,leftpair->genomepos,
				  rightpair->querypos,rightpair->genomepos,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				  pairpool,dynprogL,dynprogR,maxpeelback,extramaterial_paired,
				  extraband_paired,defect_rate,	/*forcep*/true);
      /* Do forcep, because adding large gap is not a good solution */

      if (filledp == true) {
	/* Discard the gap */
      } else {
	/* Replace the gap */
	debug(leftpair = path->first;
	      rightpair = pairs->first;
	      printf("Stage 3: Adding large gap instead: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		     leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos));
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      }

    } else {
      /* Single gap; force fill */
      leftpair = path->first;
      rightpair = pairs->first;
      debug(printf("Stage 3: Traversing single gap: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d.  queryjump = %d, genomejump = %d\n",
		   leftpair->querypos,rightpair->querypos,leftpair->genomepos,rightpair->genomepos,
		   pair->queryjump,pair->genomejump));
      pairs = traverse_single_gap(&filledp,pairs,&path,leftpair->querypos,leftpair->genomepos,
				  rightpair->querypos,rightpair->genomepos,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				  cdna_direction,pairpool,dynprogM,maxpeelback,extraband_single,defect_rate,
				  /*forcep*/true);
    }
  }

  return pairs;
}  


#if 0
/* This procedure is needed to fix the following problem

    GTTAATATA...ATATAATA...TATATATATA
    ||||||===...=== |===...===|| ||||
    GTTAAT   161   GA    1    ATCTATA

   which happens when a later gap intrudes upon an earlier gap */

static List_T
fix_short_gaps (List_T pairs, int lastquerypos, int lastgenomepos, List_T path,
#ifdef PMAP
		char *queryaaseq_ptr,
#endif
		char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		int cdna_direction, int maxpeelback, int nullgap, int extramaterial_paired, 
		int extraband_single, int extraband_paired, double defect_rate,
		Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T oldgappairs;
  Pair_T pair;
  int querypos, genomepos, queryjump, genomejump;
  char comp;
  bool filledp;

  debug(printf("\n** Starting fix_short_gaps\n"));
  while (path != NULL) {
    path = Pairpool_pop(path,&pair);
    if (pair->gapp == false) {
      pairs = Pairpool_push_existing(pairs,pairpool,pair);

    } else {
      comp = pair->comp;
      oldgappairs = (List_T) NULL;

      while (path != NULL && pair->gapp == true) {
	path = Pairpool_pop(path,&pair);
	if (pair->gapp == true) {
	  oldgappairs = Pairpool_push_existing(oldgappairs,pairpool,pair);
	} else {
	  /* Put back */
	  path = Pairpool_push_existing(path,pairpool,pair);
	}
      }
      querypos = pair->querypos;
      genomepos = pair->genomepos;

      queryjump = lastquerypos - querypos - 1;
      genomejump = lastgenomepos - genomepos - 1;

      if (genomejump > MININTRONLEN) {
	oldgappairs = List_reverse(oldgappairs);
	pairs = Pairpool_transfer(pairs,oldgappairs);
      } else if (comp == DUALBREAK_COMP) {
	oldgappairs = List_reverse(oldgappairs);
	pairs = Pairpool_transfer(pairs,oldgappairs);
      } else if (comp == NONINTRON_COMP) {
	pairs = traverse_single_gap(&filledp,pairs,&path,querypos,genomepos,lastquerypos,lastgenomepos,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				    cdna_direction,pairpool,dynprogM,maxpeelback,extraband_single,defect_rate,
				    /*forcep*/true);
      } else {
	/* Not clear what to do */
	oldgappairs = List_reverse(oldgappairs);
	pairs = Pairpool_transfer(pairs,oldgappairs);
      }
    }

    pair = pairs->first;
    lastquerypos = pair->querypos;
    lastgenomepos = pair->genomepos;
  }	

  return pairs;
}
#endif



static List_T
remove_end5_gap (List_T pairs, Pairpool_T pairpool) {
  Pair_T firstpair;

  if (pairs != NULL) {
    firstpair = pairs->first;
    pairs = Pairpool_pop(pairs,&firstpair);
    while ((firstpair->gapp || firstpair->comp == INDEL_COMP) && pairs != NULL) {
      pairs = Pairpool_pop(pairs,&firstpair);
    }
    if (firstpair->gapp == false && firstpair->comp != INDEL_COMP) {
      pairs = Pairpool_push_existing(pairs,pairpool,firstpair);
    }
  }
  return pairs;
}


static List_T
path_compute (int *intronlen, int *nonintronlen, 
	      List_T path, int cdna_direction, int introndir,
	      int querylength, int skiplength, Genomicpos_T genomiclength,
#ifdef PMAP
	      char *queryaaseq_ptr,
#endif
	      char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
	      int maxpeelback, int nullgap,
	      int extramaterial_end, int extramaterial_paired,
	      int extraband_single, int extraband_end, int extraband_paired, double defect_rate, 
	      Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
	      int indexsize, bool extend_mismatch_p, bool end_microexons_p) {
  List_T pairs = NULL;
  Pair_T firstpair, leftpair;
  int nshortexons, ndeleteexons, iter = 0;
  int newintrondir, nintrons = 0, nnonintrons = 0;
  bool singlep = true;

  if (path == NULL) {
    return NULL;
  }

#ifdef PMAP
  /* Pass 0: undefine nucleotides around gaps */
  pairs = undefine_nucleotides(queryseq_ptr,querylength,path,pairpool,/*width*/6);
  path = List_reverse(pairs);
#endif

  pairs = List_reverse(path);
  path = insert_gaps_and_singles(pairs,pairpool,queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,skiplength);
  debug3(pairs = assign_gap_types(path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction));
  debug3(return pairs);

  /* Pass 1: solve single gaps */
  debug(printf("*** Pass 1: Solve single gaps\n"));
  pairs = build_pairs_singles(path,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
			      maxpeelback,nullgap,extramaterial_paired,extraband_single,
			      defect_rate,pairpool,dynprogL,dynprogM,dynprogR);

  debug4(path = List_reverse(pairs));
  debug4(pairs = assign_gap_types(path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction));
  debug4(return pairs);

  /* Pass 2: smoothing and solving dual introns */
  while (singlep == true && iter < MAXITER) {
    if (iter > 0) {
      Smooth_reset(pairs);
    }
    pairs = Smooth_pairs(&nshortexons,&ndeleteexons,pairs,pairpool,indexsize);
    if (ndeleteexons > 0) {
      path = insert_gaps(pairs,pairpool);
      pairs = build_pairs_singles(path,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				  maxpeelback,nullgap,extramaterial_paired,extraband_single,
				  defect_rate,pairpool,dynprogL,dynprogM,dynprogR);
    }
    if (nshortexons == 0) {
      singlep = false;
    } else {
      path = List_reverse(pairs);
      debug(printf("*** Pass 2: Solve dual introns.  Iteration %d\n",iter));
      pairs = build_pairs_dualintrons(&singlep,path,
#ifdef PMAP
				      queryaaseq_ptr,
#endif
				      queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				      maxpeelback,nullgap,extramaterial_paired,extraband_paired,
				      defect_rate,pairpool,dynprogL,dynprogM,dynprogR);
    }
    iter++;
  }

  debug5(path = List_reverse(pairs));
  debug5(pairs = assign_gap_types(path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction));
  debug5(return pairs);

  if (pairs != NULL) {
    /* Pass 3: solve introns */
#ifdef CHECK
    path = check_gaps(pairs,pairpool);
#else
    path = List_reverse(pairs);
#endif

    debug(printf("*** Pass 3: Solve introns.\n"));
    pairs = build_pairs_introns(&nintrons,&nnonintrons,&(*intronlen),&(*nonintronlen),path,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				maxpeelback,nullgap,extramaterial_paired,extraband_single,extraband_paired,
				defect_rate,pairpool,dynprogL,dynprogM,dynprogR);

    debug(printf("nintrons = %d, nnonintrons = %d\n",nintrons,nnonintrons));
    debug(printf("Intronlen = %d, Nonintronlen = %d\n",*intronlen,*nonintronlen));

    if (nintrons > nnonintrons) {
      newintrondir = cdna_direction;
    } else {
      newintrondir = introndir;
    }

#if 0
    /* Pass 4: fix short gaps */
    path = List_reverse(pairs);

    firstpair = path->first;
    lastquerypos = firstpair->querypos;
    lastgenomepos = firstpair->genomepos;

    pairs = fix_short_gaps(NULL,lastquerypos,lastgenomepos,path,
#ifdef PMAP
			   queryaaseq_ptr,
#endif
			   queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
			   maxpeelback,nullgap,extramaterial_paired,extraband_single,extraband_paired,
			   defect_rate,pairpool,dynprogL,dynprogM,dynprogR);
#endif


    /* Pass 5: solve ends */
    pairs = build_pairs_end5(pairs,
#ifdef PMAP
			     queryaaseq_ptr,
#endif
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			     cdna_direction,newintrondir,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_end,
			     defect_rate,pairpool,dynprogR,extend_mismatch_p,
			     end_microexons_p);

    path = List_reverse(pairs);
    path = build_path_end3(path,querylength,genomiclength,
#ifdef PMAP
			   queryaaseq_ptr,
#endif
			   queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			   cdna_direction,newintrondir,maxpeelback,nullgap,
			   extramaterial_end,extramaterial_paired,extraband_end,
			   defect_rate,pairpool,dynprogL,extend_mismatch_p,
			   end_microexons_p);

    pairs = assign_gap_types(path,pairpool,queryseq_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction);

#ifdef CHECK
    path = check_gaps(pairs,pairpool);
    pairs = List_reverse(path);
#endif

  }

  return pairs;
}


/* Introndir shows whether we are certain about the intron direction.
   However, we try both directions anyway, using cdna_direction */
T
Stage3_compute (Stage2_T stage2, 
#ifdef PMAP
		Sequence_T queryaaseq,
#endif
		Sequence_T queryseq, Sequence_T queryuc,
		Sequence_T genomicseg, Sequence_T genomicuc,
		Matchpairend_T matchpairend, int straintype, char *strain,
		int chrnum, Genomicpos_T chroffset, Genomicpos_T chrpos, 
		bool watsonp, int maxpeelback, int nullgap, 
		int extramaterial_end, int extramaterial_paired,
		int extraband_single, int extraband_end, int extraband_paired,
		bool extend_mismatch_p, bool end_microexons_p, Pairpool_T pairpool, 
		Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		IIT_T altstrain_iit, int indexsize, int ngap) {
  T this;
  List_T pairs_fwd, pairs_rev, path_fwd, path_rev;
  double defect_rate;
  int querylength;
  int nfwdintrons, nrevintrons, nunkintrons, introndir;
  int fwd_intronlen = 0, rev_intronlen = 0;
  int fwd_nonintronlen = 0, rev_nonintronlen = 0;
  Genomicpos_T genomiclength;
  Pair_T start, end;
#ifdef PMAP
  char *queryaaseq_ptr;
#endif
  char *queryseq_ptr, *queryuc_ptr, *genomicseg_ptr, *genomicuc_ptr;
  int *typematches, nmatches;

  debug(printf("Stage 3: *** Starting stage 3 for genomiclength %u at chrnum %d:%d\n",
	       Sequence_fulllength(genomicseg),chrnum,chrpos));

  defect_rate = Stage2_defect_rate(stage2);
  nfwdintrons = Stage2_nfwdintrons(stage2);
  nrevintrons = Stage2_nrevintrons(stage2);
  nunkintrons = Stage2_nunkintrons(stage2);

  querylength = Sequence_fulllength(queryseq);
  genomiclength = (Genomicpos_T) Sequence_fulllength(genomicseg);;
#ifdef PMAP
  queryaaseq_ptr = Sequence_fullpointer(queryaaseq);
#endif
  queryseq_ptr = Sequence_fullpointer(queryseq);
  queryuc_ptr = Sequence_fullpointer(queryuc);
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);

#ifdef PMAP
  path_fwd = Stage2_path(stage2);
  path_rev = (List_T) NULL;
  introndir = +1;
#else
  if (nfwdintrons == 0 && nrevintrons == 0) {
    /* Should try both even if no introns (cf, AA011563) */
    path_fwd = Stage2_path(stage2);
    path_rev = Pairpool_copy(Stage2_path(stage2),pairpool);
    introndir = 0;
  } else if (nfwdintrons > nrevintrons + nunkintrons) {
    /* Forward wins */
    path_fwd = Stage2_path(stage2);
    path_rev = (List_T) NULL;
    introndir = +1;
  } else if (nrevintrons > nfwdintrons + nunkintrons) {
    /* Reverse wins */
    path_fwd = (List_T) NULL;
    path_rev = Stage2_path(stage2);
    introndir = -1;
  } else {
    /* Tie */
    path_fwd = Stage2_path(stage2);
    path_rev = Pairpool_copy(Stage2_path(stage2),pairpool);
    introndir = 0;
  }
#endif
  debug(printf("Introns: %d fwd, %d rev, %d unk => fwdp = %p, revp = %p\n",
	       nfwdintrons, nrevintrons, nunkintrons, path_fwd, path_rev));

  if (path_fwd != NULL) {
    pairs_fwd = path_compute(&fwd_intronlen,&fwd_nonintronlen,path_fwd,+1,introndir,
			     querylength,Sequence_skiplength(queryseq),genomiclength,
#ifdef PMAP
			     queryaaseq_ptr,
#endif
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_single,extraband_end,extraband_paired,
			     defect_rate,pairpool,dynprogL,dynprogM,dynprogR,indexsize,extend_mismatch_p,
			     end_microexons_p);
  } else {
    pairs_fwd = NULL;
  }
  if (path_rev != NULL) {
    pairs_rev = path_compute(&rev_intronlen,&rev_nonintronlen,path_rev,-1,introndir,
			     querylength,Sequence_skiplength(queryseq),genomiclength,
#ifdef PMAP
			     queryaaseq_ptr,
#endif
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_single,extraband_end,extraband_paired,
			     defect_rate,pairpool,dynprogL,dynprogM,dynprogR,indexsize,extend_mismatch_p,
			     end_microexons_p);
  } else {
    pairs_rev = NULL;
  }

  if ((this = Stage3_new(pairs_fwd,pairs_rev)) == NULL) {
    return NULL;
  } else if (this->cdna_direction >= 0) {
    this->pairarray = prepare_for_printing(&this->npairs,&this->pairs_fwd,this->cdna_direction,pairpool,
					   queryseq_ptr,genomicseg_ptr,genomicuc_ptr,ngap,
					   Sequence_subseq_offset(queryseq),Sequence_skiplength(queryseq));
  } else {
    this->pairarray = prepare_for_printing(&this->npairs,&this->pairs_rev,/*cdna_direction*/-1,pairpool,
					   queryseq_ptr,genomicseg_ptr,genomicuc_ptr,ngap,
					   Sequence_subseq_offset(queryseq),Sequence_skiplength(queryseq));
  }
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
  this->matchpairend = matchpairend;

  if (watsonp) {
    this->genomicstart = chroffset + chrpos + Pair_genomepos(start);
    this->genomicend = chroffset + chrpos + Pair_genomepos(end);
  } else {
    this->genomicstart = chroffset + chrpos + (genomiclength - 1) - Pair_genomepos(start);
    this->genomicend = chroffset + chrpos + (genomiclength - 1) - Pair_genomepos(end);
  }

  /*
    new->coverage_correction = 
    Sequence_count_bad(genomicseg,start->genomepos,start->querypos,-1) +
    Sequence_count_bad(genomicseg,end->genomepos,(querylength - 1) - end->querypos,+1);
  */

  /* new->coverage = (double) (end->querypos - start->querypos + 1 + new->coverage_correction)/(double) querylength; */
  this->coverage = (double) (end->querypos - start->querypos + 1)/(double) (querylength + Sequence_skiplength(queryseq));

  if (straintype == 0) {
    return this;
  } else {
    if (watsonp) {
      typematches = IIT_get_typed(&nmatches,altstrain_iit,this->genomicstart,this->genomicend,straintype);
    } else {
      typematches = IIT_get_typed(&nmatches,altstrain_iit,this->genomicend,this->genomicstart,straintype);
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


T
Stage3_direct (Matchpair_T matchpair,
#ifdef PMAP
	       Sequence_T queryaaseq,
#endif
	       Sequence_T queryseq, Sequence_T queryuc, Pairpool_T pairpool, Genome_T genome,
	       Matchpairend_T matchpairend, Chrnum_T chrnum,  Genomicpos_T chroffset, Genomicpos_T chrpos, bool watsonp,
	       int ngap, char *gbuffer1, char *gbuffer2, int gbufferlen, Dynprog_T dynprogL, Dynprog_T dynprogR,
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

  querylength = Sequence_fulllength(queryseq);
  queryseq_ptr = Sequence_fullpointer(queryseq);
#ifdef PMAP
  queryaaseq_ptr = Sequence_fullpointer(queryaaseq);
#endif
  queryuc_ptr = Sequence_fullpointer(queryuc);
  path = Matchpair_make_path(matchpair,queryseq_ptr,pairpool,genome,chroffset,chrpos,gbuffer1,gbuffer2,gbufferlen);

  leftpair = (Pair_T) List_head(path);
  querydp5 = leftpair->querypos + 1;
  genomicpos = chroffset + chrpos + leftpair->genomepos;
  queryjump = querylength - leftpair->querypos - 1;
  genomejump = queryjump + extramaterial_end; /* proposed */
  debug(printf("For 3' end, genomicpos is %u, queryjump = %d, genomejump = %d\n",genomicpos,queryjump,genomejump));

  if (watsonp == true) {
    genomicseg = Genome_get_segment(genome,genomicpos+1U,genomejump,/*revcomp*/false,gbuffer1,gbuffer2,gbufferlen);
    genomedp5 = leftpair->genomepos + 1U;
  } else {
    genomicseg = Genome_get_segment(genome,genomicpos-genomejump,genomejump,/*revcomp*/true,gbuffer1,gbuffer2,gbufferlen);
    genomedp5 = leftpair->genomepos - 1U;
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc = Sequence_uppercase(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);

  gappairs = Dynprog_end3_gap(&finalscore,&nmatches,&nmismatches,&nopens,&nindels,dynprogL,
			      /*sequence1*/&(queryseq_ptr[querydp5]),/*sequenceuc1*/&(queryuc_ptr[querydp5]),
			      /*sequence2*/genomicseg_ptr,/*sequenceuc2*/genomicuc_ptr,
			      /*length1*/queryjump,/*length2*/genomejump,/*offset1*/querydp5,/*offset2*/0,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      /*cdna_direction*/+1,pairpool,extraband_end,/*defect_rate*/0.0,/*extend_mismatch_p*/true);
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
    genomicseg = Genome_get_segment(genome,genomicpos-genomejump,genomejump,/*revcomp*/false,gbuffer1,gbuffer2,gbufferlen);
    genomedp3 = rightpair->genomepos - 1U;
  } else {
    genomicseg = Genome_get_segment(genome,genomicpos+1U,genomejump,/*revcomp*/true,gbuffer1,gbuffer2,gbufferlen);
    genomedp3 = rightpair->genomepos + 1U;
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc = Sequence_uppercase(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);

  gappairs = Dynprog_end5_gap(&finalscore,&nmatches,&nmismatches,&nopens,&nindels,dynprogR,
			      &(queryseq_ptr[querydp3]),&(queryuc_ptr[querydp3]),
			      &(genomicseg_ptr[genomejump-1]),&(genomicuc_ptr[genomejump-1]),
			      queryjump,genomejump,/*offset1*/querydp3,/*offset2*/0,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      /*cdna_direction*/+1,pairpool,extraband_end,/*defect_rate*/0.0,/*extend_mismatch_p*/true);
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

  new = Stage3_new(pairs_fwd,/*pairs_rev*/NULL);
  new->cdna_direction = 0;	/* Override the +1 assigned by Stage3_new */

  new->pairarray = prepare_for_printing(&new->npairs,&new->pairs_fwd,new->cdna_direction,pairpool,
					queryseq_ptr,/*genomicseg_ptr*/NULL,/*genomicuc_ptr*/NULL,ngap,
					/*subseq_offset*/0,Sequence_skiplength(queryseq));
  for (i = 0; i < new->npairs; i++) {
    pair = &(new->pairarray[i]);
    pair->aapos = 0;
    pair->aa_g = ' ';
    pair->aa_e = ' ';
    pair->aamarker_g = false;
    pair->aamarker_e = false;
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
  new->matchpairend = matchpairend;

  debug(printf("chroffset = %u, chrpos = %u, genomiclength = %u\n",chroffset,chrpos,new->genomiclength));
  if (watsonp) {
    new->genomicstart = chroffset + chrpos + Pair_genomepos(start);
    new->genomicend = chroffset + chrpos + Pair_genomepos(end);
  } else {
    new->genomicstart = chroffset + chrpos + (new->genomiclength - 1) - Pair_genomepos(start);
    new->genomicend = chroffset + chrpos + (new->genomiclength - 1) - Pair_genomepos(end);
  }

  /*
    new->coverage_correction = 
    Sequence_count_bad(genomicseg,start->genomepos,start->querypos,-1) +
    Sequence_count_bad(genomicseg,end->genomepos,(new->querylength - 1) - end->querypos,+1);
  */

  /* new->coverage = (double) (end->querypos - start->querypos + 1 + new->coverage_correction)/(double) new->querylength; */
  new->coverage = (double) (end->querypos - start->querypos + 1)/(double) (Sequence_fulllength(queryseq)+Sequence_skiplength(queryseq));

  return new;
}



/************************************************************************
 *  Merging
 ************************************************************************/

bool
Stage3_mergeable (char *comp, Stage3_T firstpart, Stage3_T secondpart, int chimera_direction, 
		  double donor_prob, double acceptor_prob) {
  Pair_T end1, start2, pair;
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
  List_T pairs = NULL, p;
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
    pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',comp,' ',/*gapp*/true);
    pair = (Pair_T) pairs->first;
    for (i = 0; i < ngap; i++) {
      memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    }
    pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,/*gapp*/true);
    pair = (Pair_T) pairs->first;
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',comp,' ',/*gapp*/true);
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
    genomicseg = Genome_get_segment(genome,left,ngap,!watsonp,gbuffer1,gbuffer2,ngap);
    genomicseg_ptr = Sequence_fullpointer(genomicseg);

    for (i = 0; i < ngap; i++) {
      pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',comp,genomicseg_ptr[i],/*gapp*/true);
      pair = (Pair_T) pairs->first;
      memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    }
    Sequence_free(&genomicseg);

    pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR,/*gapp*/true);
    pair = (Pair_T) pairs->first;
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));
    memcpy(&(firstpart->pairarray[ptr++]),pair,sizeof(struct Pair_T));

    if (watsonp == true) {
      left = firstpart->chroffset + firstpart->chrpos + start2->genomepos - ngap;
    } else {
      left = firstpart->chroffset + firstpart->chrpos + (firstpart->genomiclength - 1) - start2->genomepos + 1;
    }

    genomicseg = Genome_get_segment(genome,left,ngap,!watsonp,gbuffer1,gbuffer2,ngap);
    genomicseg_ptr = Sequence_fullpointer(genomicseg);

    for (i = 0; i < ngap; i++) {
      pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',comp,genomicseg_ptr[i],/*gapp*/true);
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

  firstpart->nexons += secondpart->nexons;
  firstpart->matches += secondpart->matches;
  firstpart->unknowns += secondpart->unknowns;
  firstpart->mismatches += secondpart->mismatches;
  firstpart->qopens += secondpart->qopens;
  firstpart->qindels += secondpart->qindels;
  firstpart->topens += secondpart->topens;
  firstpart->tindels += secondpart->tindels;

  return;
}

