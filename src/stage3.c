static char rcsid[] = "$Id: stage3.c,v 1.196 2005/10/25 16:50:15 twu Exp $";
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
#include "smooth.h"
#include "scores.h"
#include "translation.h"

/* The following are the same as dividing by 2048 and 1024 */
#define goodness_intronlen(x) (x >> 11)
#define goodness_nonintronlen(x) (x >> 10)

#define MAXITER 4
#define INTRON_PENALTY_INCONSISTENT 16
#define MININTRONLEN 6
#define CROSS_ONE_SHORT

#define MERGELENGTH 100000
#define LONG_MERGELENGTH 500000	/* For strong donor and acceptor splice sites */
#define DONOR_THRESHOLD 0.90
#define ACCEPTOR_THRESHOLD 0.90

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

/* Show result of dual introns.  Alters computation. */
#ifdef DEBUG6
#define debug6(x) x
#else 
#define debug6(x)
#endif


/************************************************************************
 *   Stage 3 merges cDNA-genomic pairs from stage 2 (called "path")
 *   and from dynamic programming into a single list.  In this
 *   process, stage 3 may also have to pop a pair off the path, insert
 *   the dynamic programming results (called "gappairs") and then push
 *   the stored pair onto the list.  The relevant pointers are
 *   querypos and genomepos, which refer to the stored pair;
 *   lastquerypos and lastgenomepos, which refer to the top pair on
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
  struct Pair_T *pairs;		/* The array version */
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
  int querylength;
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
  return Pair_querypos(&(this->pairs[0]));
}

int
Stage3_queryend (T this) {
  return Pair_querypos(&(this->pairs[this->npairs-1]));
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

  querystart = Pair_querypos(&(this->pairs[0]));
  queryend = Pair_querypos(&(this->pairs[this->npairs-1]));

  return queryend - querystart + 1;
}


int
Stage3_largemargin (int *newstart, int *newend, T this, int queryntlength) {
  int leftmargin, rightmargin;
  int querystart, queryend;

  querystart = Pair_querypos(&(this->pairs[0]));
  queryend = Pair_querypos(&(this->pairs[this->npairs-1]));

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

  genomicpos = Pair_genomicpos(this->pairs,this->npairs,querypos,headp);
  if (this->watsonp) {
    return this->chroffset + this->chrpos + genomicpos;
  } else {
    return this->chroffset + this->chrpos + (this->genomiclength - 1) - genomicpos;
  }
}


int *
Stage3_matchscores (T this, int querylength) {
  return Pair_matchscores(this->pairs,this->npairs,this->cdna_direction,querylength);
}

void
Stage3_pathscores (int *pathscores, T this, int querylength, cDNAEnd_T cdnaend) {
  Pair_pathscores(pathscores,this->pairs,this->npairs,this->cdna_direction,querylength,cdnaend);
  return;
}


int
Stage3_chimeric_goodness (int *matches1, int *matches2, T part1, T part2, int breakpoint, int querylength) {
  int goodness1, goodness2, querystart, queryend;
  int unknowns1, mismatches1, qopens1, qindels1, topens1, tindels1, 
    ncanonical1, nsemicanonical1, nnoncanonical1;
  int unknowns2, mismatches2, qopens2, qindels2, topens2, tindels2, 
    ncanonical2, nsemicanonical2, nnoncanonical2;

  querystart = Pair_querypos(&(part1->pairs[0]));
  debug2(printf("Chimeric goodness requested for part %d..%d\n",querystart+1,breakpoint));
  Pair_fracidentity_bounded(&(*matches1),&unknowns1,&mismatches1,&qopens1,&qindels1,&topens1,&tindels1,
			    &ncanonical1,&nsemicanonical1,&nnoncanonical1,
			    part1->pairs,part1->npairs,part1->cdna_direction,
			    querystart,breakpoint);
  goodness1 = (*matches1) + MISMATCH*mismatches1 + QOPEN*qopens1 + QINDEL*qindels1 + TOPEN*topens1 + TINDEL*tindels1;
  debug2(printf("  %d matches, %d mismatches, %d qopens, %d qindels, %d topens, %d tindels => %d\n",
		*matches1,mismatches1,qopens1,qindels1,topens1,tindels1,goodness1));

  queryend = Pair_querypos(&(part2->pairs[part2->npairs-1]));
  debug2(printf("Chimeric goodness requested for part %d..%d\n",breakpoint+1,queryend+1));
  Pair_fracidentity_bounded(&(*matches2),&unknowns2,&mismatches2,&qopens2,&qindels2,&topens2,&tindels2,
			    &ncanonical2,&nsemicanonical2,&nnoncanonical2,
			    part2->pairs,part2->npairs,part2->cdna_direction,
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

static T
Stage3_new (List_T pairs_fwd, List_T pairs_rev, Matchpairend_T matchpairend,
	    int straintype, char *strain, Chrnum_T chrnum, Genomicpos_T chrpos,
	    Genomicpos_T chroffset, int querylength, Genomicpos_T genomiclength, 
	    bool watsonp, double defect_rate, int ngap) {
  T new;
  int matches_fwd, matches_rev, mismatches_fwd, mismatches_rev, 
    unknowns_fwd, unknowns_rev, qopens_fwd, qindels_fwd, qopens_rev, qindels_rev, 
    topens_fwd, tindels_fwd, topens_rev, tindels_rev,
    ncanonical_fwd, ncanonical_rev, nsemicanonical_fwd, nsemicanonical_rev,
    nnoncanonical_fwd, nnoncanonical_rev;
  int goodness_fwd, goodness_rev;
  Pair_T start, end;

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
      new->npairs = List_length(pairs_fwd);
      new->pairs = Pair_block_copy(pairs_fwd,&new->npairs,ngap);
      new->cdna_direction = +1;
      new->matches = matches_fwd;
      new->unknowns = unknowns_fwd;
      new->mismatches = mismatches_fwd;
      new->qopens = qopens_fwd;
      new->qindels = qindels_fwd;
      new->topens = topens_fwd;
      new->tindels = tindels_fwd;
      new->goodness = goodness_fwd - CANONICAL_POINTS*ncanonical_fwd - CANONICAL_POINTS*nnoncanonical_fwd;
    } else if (pairs_fwd == NULL) {
      new->npairs = List_length(pairs_rev);
      new->pairs = Pair_block_copy(pairs_rev,&new->npairs,ngap);
      new->cdna_direction = -1;
      new->matches = matches_rev;
      new->unknowns = unknowns_rev;
      new->mismatches = mismatches_rev;
      new->qopens = qopens_rev;
      new->qindels = qindels_rev; 
      new->topens = topens_rev;
      new->tindels = tindels_rev;
      new->goodness = goodness_rev - CANONICAL_POINTS*ncanonical_rev - CANONICAL_POINTS*nnoncanonical_rev;
    } else if (goodness_fwd > goodness_rev) {
      new->npairs = List_length(pairs_fwd);
      new->pairs = Pair_block_copy(pairs_fwd,&new->npairs,ngap);
      new->cdna_direction = +1;
      new->matches = matches_fwd;
      new->unknowns = unknowns_fwd;
      new->mismatches = mismatches_fwd;
      new->qopens = qopens_fwd;
      new->qindels = qindels_fwd;
      new->tindels = tindels_fwd;
      new->topens = topens_fwd;
      new->goodness = goodness_fwd - CANONICAL_POINTS*ncanonical_fwd - CANONICAL_POINTS*nnoncanonical_fwd;
    } else if (goodness_fwd < goodness_rev) {
      new->npairs = List_length(pairs_rev);
      new->pairs = Pair_block_copy(pairs_rev,&new->npairs,ngap);
      new->cdna_direction = -1;
      new->matches = matches_rev;
      new->unknowns = unknowns_rev;
      new->mismatches = mismatches_rev;
      new->qopens = qopens_rev;
      new->qindels = qindels_rev;
      new->topens = topens_rev;
      new->tindels = tindels_rev;
      new->goodness = goodness_rev - CANONICAL_POINTS*ncanonical_rev - CANONICAL_POINTS*nnoncanonical_rev;
    } else {
      new->npairs = List_length(pairs_fwd);
      new->pairs = Pair_block_copy(pairs_fwd,&new->npairs,ngap);
      new->cdna_direction = 0;	/* indeterminate */
      new->matches = matches_fwd;
      new->unknowns = unknowns_fwd;
      new->mismatches = mismatches_fwd;
      new->qopens = qopens_fwd;
      new->qindels = qindels_fwd;
      new->tindels = tindels_fwd;
      new->topens = topens_fwd;
      new->goodness = goodness_fwd - CANONICAL_POINTS*ncanonical_fwd - CANONICAL_POINTS*nnoncanonical_fwd;
    }

    new->nexons = Pair_nexons(new->pairs,new->npairs);
    if (new->nexons > 2) {
      /* Favor spliced transcripts, but only if we're sure they're
         spliced (i.e., 3 or more exons).  A random intron shouldn't
         get credit. */
      new->goodness += new->nexons;
    }
    
    new->straintype = straintype;
    new->strain = strain;
    new->chrnum = chrnum;
    new->chrpos = chrpos;
    new->chroffset = chroffset;
    new->genomiclength = genomiclength;
    new->watsonp = watsonp;
    new->defect_rate = defect_rate;
    new->matchpairend = matchpairend;

    start = &(new->pairs[0]);
    end = &(new->pairs[new->npairs-1]);

    if (watsonp) {
      new->genomicstart = chroffset + chrpos + Pair_genomepos(start);
      new->genomicend = chroffset + chrpos + Pair_genomepos(end);
    } else {
      new->genomicstart = chroffset + chrpos + (genomiclength - 1) - Pair_genomepos(start);
      new->genomicend = chroffset + chrpos + (genomiclength - 1) - Pair_genomepos(end);
    }

    /*
    new->coverage_correction = 
      Sequence_count_bad(genomicseg,start->genomepos,start->querypos,-1) +
      Sequence_count_bad(genomicseg,end->genomepos,(querylength - 1) - end->querypos,+1);
    */

    /* new->coverage = (double) (end->querypos - start->querypos + 1 + new->coverage_correction)/(double) querylength; */
    new->querylength = querylength;
    new->coverage = (double) (end->querypos - start->querypos + 1)/(double) querylength;

    return new;
  }
}
	    
void
Stage3_free (T *old) {

  if (*old) {
    /* Don't free strain.  Belongs to altstrain_iit. */
    FREE((*old)->pairs);
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

T
Stage3_copy (T old, Pairpool_T pairpool) {
  List_T pairs_fwd, pairs_rev;

  pairs_fwd = Pairpool_transfer_copy(NULL,old->pairs_fwd,pairpool);
  pairs_rev = Pairpool_transfer_copy(NULL,old->pairs_rev,pairpool);

  pairs_fwd = List_reverse(pairs_fwd);
  pairs_rev = List_reverse(pairs_rev);

  /* ngap of -1 turns off merge_gaps() in pair.c */
  return Stage3_new(pairs_fwd,pairs_rev,old->matchpairend,old->straintype,old->strain,
		    old->chrnum,old->chrpos,old->chroffset,
		    old->querylength,old->genomiclength,
		    old->watsonp,old->defect_rate,/*ngap*/-1);
}

T
Stage3_copy_bounded (T old, int minpos, int maxpos) {
  List_T pairs_fwd, pairs_rev;

  pairs_fwd = Pairpool_transfer_bounded(NULL,old->pairs_fwd,minpos,maxpos);
  pairs_rev = Pairpool_transfer_bounded(NULL,old->pairs_rev,minpos,maxpos);

  pairs_fwd = List_reverse(pairs_fwd);
  pairs_rev = List_reverse(pairs_rev);

  /* ngap of -1 turns off merge_gaps() in pair.c */
  return Stage3_new(pairs_fwd,pairs_rev,old->matchpairend,old->straintype,old->strain,
		    old->chrnum,old->chrpos,old->chroffset,
		    old->querylength,old->genomiclength,
		    old->watsonp,old->defect_rate,/*ngap*/-1);
}


#ifdef PMAP
Stage3_T
Stage3_translate_cdna (T this, Sequence_T queryaaseq) {
  Translation_via_cdna(&this->translation_start,&this->translation_end,&this->translation_length,
		       &this->relaastart,&this->relaaend,
		       this->pairs,this->npairs,Sequence_fullpointer(queryaaseq));
  return this;
}
#else
Stage3_T
Stage3_truncate_fulllength (Stage3_T old, bool translatep) {
  Stage3_T new;

  if (translatep == true) {
    if (old->cdna_direction < 0) {
      Translation_via_genomic(&old->translation_start,&old->translation_end,&old->translation_length,
			      &old->relaastart,&old->relaaend,
			      old->pairs,old->npairs,/*backwardsp*/true,/*revcompp*/true,/*fulllengthp*/true);
    } else {
      Translation_via_genomic(&old->translation_start,&old->translation_end,&old->translation_length,
			      &old->relaastart,&old->relaaend,
			      old->pairs,old->npairs,/*backwardsp*/false,/*revcompp*/false,/*fulllengthp*/true);
    }
  }

  if ((new = Stage3_copy_bounded(old,Stage3_translation_start(old),
				 Stage3_translation_end(old))) == NULL) {
    return old;
  } else {
    Stage3_free(&old);
    return new;
  }
}

Stage3_T
Stage3_translate_genomic (T this, bool fulllengthp, bool truncatep) {
  Stage3_T new;

  if (this->cdna_direction < 0) {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairs,this->npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp);
  } else {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairs,this->npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp);
  }
  if (truncatep == false) {
    return this;
  } else {
    new = Stage3_truncate_fulllength(this,/*translatep*/false);
    if (new->cdna_direction < 0) {
      Translation_via_genomic(&new->translation_start,&new->translation_end,&new->translation_length,
			      &new->relaastart,&new->relaaend,
			      new->pairs,new->npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp);
    } else {
      Translation_via_genomic(&new->translation_start,&new->translation_end,&new->translation_length,
			      &new->relaastart,&new->relaaend,
			      new->pairs,new->npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp);
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
				this->pairs,this->npairs,this->watsonp,/*backwardsp*/true,/*revcompp*/true,
				reference->pairs,reference->npairs,reference->watsonp,reference->genomiclength,
				fixshiftp);
    } else {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairs,this->npairs,this->watsonp,/*backwardsp*/false,/*revcompp*/false,
				reference->pairs,reference->npairs,reference->watsonp,reference->genomiclength,
				fixshiftp);
    }
  } else {
    if (reference->cdna_direction < 0) {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairs,this->npairs,this->watsonp,/*backwardsp*/false,/*revcompp*/false,
				reference->pairs,reference->npairs,reference->watsonp,reference->genomiclength,
				fixshiftp);
    } else {
      Translation_via_reference(&this->relaastart,&this->relaaend,
				this->pairs,this->npairs,this->watsonp,/*backwardsp*/true,/*revcompp*/true,
				reference->pairs,reference->npairs,reference->watsonp,reference->genomiclength,
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


static void
adjust_genomepos (T this, int delta) {
  int i;
  Pair_T pair;

  for (i = 0; i < this->npairs; i++) {
    pair = &(this->pairs[i]);
    pair->genomepos += delta;
  }
  return;
}


bool
Stage3_append (T firstpart, T secondpart, int chimera_direction, 
	       double donor_prob, double acceptor_prob,
	       Pairpool_T pairpool, Genome_T genome, int ngap) {
  Pair_T end1, start2, pair;
  struct Pair_T *oldpairs;
  List_T pairs = NULL;
  bool watsonp, connectablep = false;
  Genomicpos_T endchrpos1, startchrpos2, left;
  int querygap, genomicpos, querypos, i, ptr;
  char comp, *gbuffer1, *gbuffer2, *genomicseg_ptr;
  Sequence_T genomicseg;

  if (firstpart->chrnum != secondpart->chrnum) {
    return false;
  } else if (firstpart->watsonp != secondpart->watsonp) {
    return false;
  } else {
    watsonp = firstpart->watsonp;
    end1 = &(firstpart->pairs[firstpart->npairs-1]);
    start2 = &(secondpart->pairs[0]);
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
	comp = DUALBREAK_COMP;
      } else if (chimera_direction > 0) {
	comp = FWD_CANONICAL_INTRON_COMP;
      } else if (chimera_direction < 0) {
	comp = REV_CANONICAL_INTRON_COMP;
      } else {
	comp = NONINTRON_COMP;
      }

      if (watsonp == true) {
	adjust_genomepos(secondpart,secondpart->chrpos - firstpart->chrpos);
      } else {
	adjust_genomepos(secondpart,firstpart->chrpos + firstpart->genomiclength 
			 - secondpart->chrpos - secondpart->genomiclength);
      }

      oldpairs = firstpart->pairs;
      firstpart->pairs = (struct Pair_T *) CALLOC(firstpart->npairs+ngap+3+ngap+secondpart->npairs,sizeof(struct Pair_T));
      memcpy(firstpart->pairs,oldpairs,firstpart->npairs*sizeof(struct Pair_T));
      ptr = firstpart->npairs;

      genomicpos = start2->genomepos;
      querypos = start2->querypos;
      if (comp == DUALBREAK_COMP) {
	pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',comp,' ');
	pair = (Pair_T) pairs->first;
	for (i = 0; i < ngap; i++) {
	  memcpy(&(firstpart->pairs[ptr++]),pair,sizeof(struct Pair_T));
	}
	pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
	pair = (Pair_T) pairs->first;
	memcpy(&(firstpart->pairs[ptr++]),pair,sizeof(struct Pair_T));
	memcpy(&(firstpart->pairs[ptr++]),pair,sizeof(struct Pair_T));
	memcpy(&(firstpart->pairs[ptr++]),pair,sizeof(struct Pair_T));
	pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',comp,' ');
	pair = (Pair_T) pairs->first;
	for (i = 0; i < ngap; i++) {
	  memcpy(&(firstpart->pairs[ptr++]),pair,sizeof(struct Pair_T));
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
	  pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',comp,genomicseg_ptr[i]);
	  pair = (Pair_T) pairs->first;
	  memcpy(&(firstpart->pairs[ptr++]),pair,sizeof(struct Pair_T));
	}
	Sequence_free(&genomicseg);

	pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',INTRONGAP_COMP,INTRONGAP_CHAR);
	pair = (Pair_T) pairs->first;
	memcpy(&(firstpart->pairs[ptr++]),pair,sizeof(struct Pair_T));
	memcpy(&(firstpart->pairs[ptr++]),pair,sizeof(struct Pair_T));
	memcpy(&(firstpart->pairs[ptr++]),pair,sizeof(struct Pair_T));

	if (watsonp == true) {
	  left = firstpart->chroffset + firstpart->chrpos + start2->genomepos - ngap;
	} else {
	  left = firstpart->chroffset + firstpart->chrpos + (firstpart->genomiclength - 1) - start2->genomepos + 1;
	}

	genomicseg = Genome_get_segment(genome,left,ngap,!watsonp,gbuffer1,gbuffer2,ngap);
	genomicseg_ptr = Sequence_fullpointer(genomicseg);

	for (i = 0; i < ngap; i++) {
	  pairs = Pairpool_push(pairs,pairpool,querypos,genomicpos,' ',comp,genomicseg_ptr[i]);
	  pair = (Pair_T) pairs->first;
	  memcpy(&(firstpart->pairs[ptr++]),pair,sizeof(struct Pair_T));
	}
	Sequence_free(&genomicseg);
	
	FREE(gbuffer2);
	FREE(gbuffer1);
      }
      FREE(oldpairs);		/* Must free after last access to end1 */

      memcpy(&(firstpart->pairs[ptr]),secondpart->pairs,secondpart->npairs*sizeof(struct Pair_T));
      firstpart->npairs = firstpart->npairs+ngap+3+ngap+secondpart->npairs;

      firstpart->nexons += secondpart->nexons;
      firstpart->matches += secondpart->matches;
      firstpart->unknowns += secondpart->unknowns;
      firstpart->mismatches += secondpart->mismatches;
      firstpart->qopens += secondpart->qopens;
      firstpart->qindels += secondpart->qindels;
      firstpart->topens += secondpart->topens;
      firstpart->tindels += secondpart->tindels;

      return true;
    }
  }
}


void
Stage3_print_pathsummary (T this, int pathnum, IIT_T chromosome_iit, IIT_T contig_iit, 
			  IIT_T altstrain_iit, char *dbversion, bool zerobasedp) {
  Pair_T start, end;
  bool referencealignp;

  start = &(this->pairs[0]);
  end = &(this->pairs[this->npairs-1]);
  referencealignp = this->straintype == 0 ? true : false;
  Pair_print_pathsummary(pathnum,start,end,this->chrnum,this->chrpos,this->chroffset,
			 chromosome_iit,referencealignp,altstrain_iit,this->strain,contig_iit,
			 dbversion,this->genomiclength,
			 this->nexons,this->coverage,
			 this->matches,this->unknowns,this->mismatches,
			 this->qopens,this->qindels,this->topens,this->tindels,this->goodness,
			 this->watsonp,this->cdna_direction,this->defect_rate,
			 this->translation_start,this->translation_end,this->translation_length,
			 0,0,zerobasedp);
  Translation_compare(this->pairs,this->npairs,NULL,0,this->cdna_direction,
		      this->relaastart,this->relaaend);
  printf("\n");

  return;
}

void
Stage3_print_pslformat_nt (T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq) {
  Pair_T start, end;

  start = &(this->pairs[0]);
  end = &(this->pairs[this->npairs-1]);

  Pair_print_pslformat_nt(this->pairs,this->npairs,start,end,queryseq,this->chrnum,this->chrpos,
			  chromosome_iit,this->genomiclength,this->nexons,this->matches,this->unknowns,this->mismatches, 
			  this->qopens,this->qindels,this->topens,this->tindels,
			  this->watsonp,this->cdna_direction);
  return;
}

void
Stage3_print_pslformat_pro (T this, int pathnum, IIT_T chromosome_iit, Sequence_T queryseq) {
  Pair_T start, end;

  start = &(this->pairs[0]);
  end = &(this->pairs[this->npairs-1]);

  Pair_print_pslformat_pro(this->pairs,this->npairs,start,end,queryseq,this->chrnum,this->chrpos,
			   chromosome_iit,this->genomiclength,this->nexons,this->matches,this->unknowns,this->mismatches, 
			   this->qopens,this->qindels,this->topens,this->tindels,
			   this->watsonp,this->cdna_direction);
  return;
}


void
Stage3_print_mutations (T this, T reference, IIT_T chromosome_iit, char *dbversion,
			bool showalignp, bool zerobasedp, 
			bool continuousp, bool diagnosticp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength, int ngap) {
  Pair_T start, end;
  bool referencealignp;

  start = &(this->pairs[0]);
  end = &(this->pairs[this->npairs-1]);

  /*  Pair_dump_array(this->pairs,this->npairs,false); */

  referencealignp = this->straintype == 0 ? true : false;
  Pair_print_pathsummary(/*pathnum*/1,start,end,reference->chrnum,reference->chrpos,reference->chroffset,
			 chromosome_iit,referencealignp,/*altstrain_iit*/NULL,this->strain,/*contig_iit*/NULL,
			 dbversion,reference->genomiclength,
			 this->nexons,this->coverage,
			 this->matches,this->unknowns,this->mismatches,
			 this->qopens,this->qindels,this->topens,this->tindels,this->goodness,
			 this->watsonp,this->cdna_direction,this->defect_rate,
			 0,0,0,this->relaastart,this->relaaend,zerobasedp);
  Translation_compare(this->pairs,this->npairs,reference->pairs,reference->npairs,
		      this->cdna_direction,this->relaastart,this->relaaend);
  printf("\n");

  if (showalignp == true) {
    Pair_print_alignment(this->pairs,this->npairs,reference->chrnum,reference->chrpos,reference->chroffset,
			 chromosome_iit,reference->genomiclength,this->watsonp,this->cdna_direction,/*universalp*/false,zerobasedp,
			 diagnosticp,/*genomicprop*/false,invertmode,nointronlenp,wraplength,ngap);
  }
  debug1(Pair_dump_array(this->pairs,this->npairs,/*zerobasedp*/true));
  debug1(Pair_check_array(this->pairs,this->npairs));

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

  start = &(this->pairs[0]);
  end = &(this->pairs[this->npairs-1]);

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

  exonbounds = Pair_exonbounds(this->pairs,this->npairs,this->chrpos,this->chroffset,this->genomiclength,
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
			int invertmode, bool nointronlenp, int wraplength, int ngap) {
  if (continuousp == true) {
#ifdef PMAP
    this = Stage3_translate_cdna(this,queryaaseq);
#endif
    Pair_print_continuous(this->pairs,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			  this->genomiclength,this->watsonp,this->cdna_direction,universalp,zerobasedp,
			  diagnosticp,genomefirstp,invertmode,nointronlenp,ngap);
  } else {
    Pair_print_exonsummary(this->pairs,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			   chromosome_iit,this->genomiclength,this->watsonp,universalp,zerobasedp,
			   genomefirstp,invertmode);
    if (alignsummaryonlyp == false) {
#ifdef PMAP
      this = Stage3_translate_cdna(this,queryaaseq);
#endif
      Pair_print_alignment(this->pairs,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			   chromosome_iit,this->genomiclength,this->watsonp,
			   this->cdna_direction,universalp,zerobasedp,
			   diagnosticp,/*genomicprop*/true,invertmode,nointronlenp,wraplength,ngap);
    }
  }
  debug1(Pair_dump_array(this->pairs,this->npairs,/*zerobasedp*/true));
  debug1(Pair_check_array(this->pairs,this->npairs));
  return;
}


/* queryaaseq is used only by PMAP */
void
Stage3_print_coordinates (T this, Sequence_T queryaaseq, IIT_T chromosome_iit, bool zerobasedp,
			  int invertmode) {
#ifdef PMAP
  this = Stage3_translate_cdna(this,queryaaseq);
#endif
  Pair_print_coordinates(this->pairs,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			 chromosome_iit,this->genomiclength,this->watsonp,
			 zerobasedp,invertmode);
  return;
}


void
Stage3_print_cdna_exons (T this, int wraplength) {
  Pair_print_cdna_exons(this->pairs,this->npairs,wraplength);
  return;
}


#ifdef PMAP
void
Stage3_print_nucleotide_cdna (T this, int wraplength) {
  Pair_print_nucleotide_cdna(this->pairs,this->npairs,wraplength);
  return;
}
#else
void
Stage3_print_protein_cdna (T this, int wraplength) {
  Pair_print_protein_cdna(this->pairs,this->npairs,wraplength);
  return;
}
#endif

void
Stage3_print_protein_genomic (T this, int wraplength) {
  Pair_print_protein_genomic(this->pairs,this->npairs,wraplength);
  return;
}


void
Stage3_print_compressed (T this, Sequence_T queryseq, IIT_T chromosome_iit,
			 char *version, int pathnum, int npaths,
			 bool checksump, int chimerapos, int chimeraequivpos,
			 double donor_prob, double acceptor_prob, 
			 int chimera_cdna_direction, bool zerobasedp) {
  Pair_print_compressed(queryseq,version,pathnum,npaths,
			this->nexons,this->coverage,Stage3_fracidentity(this),
			this->pairs,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			chromosome_iit,this->genomiclength,checksump,
			chimerapos,chimeraequivpos,donor_prob,acceptor_prob,
			chimera_cdna_direction,this->strain,this->watsonp,
			this->cdna_direction,zerobasedp);
  return;
}

static List_T 
add_null_gap (List_T pairs, int lastquerypos, int lastgenomepos, Pairpool_T pairpool,
	      int ngap) {
  int k;

  for (k = 0; k < ngap; k++) {
    pairs = Pairpool_push(pairs,pairpool,lastquerypos,lastgenomepos,' ',DUALBREAK_COMP,' ');
  }
  for (k = 0; k < 3; k++) {
    pairs = Pairpool_push(pairs,pairpool,lastquerypos,lastgenomepos,' ',INTRONGAP_COMP,' ');
  }
  for (k = 0; k < ngap; k++) {
    pairs = Pairpool_push(pairs,pairpool,lastquerypos,lastgenomepos,' ',DUALBREAK_COMP,' ');
  }

  return pairs;
}

static List_T
remove_null_gaps (List_T pairs, Pairpool_T pairpool) {
  List_T newpairs = NULL;
  Pair_T pair;

  while (pairs != NULL) {
    pair = pairs->first;
    if (pair->comp != DUALBREAK_COMP && pair->comp != INTRONGAP_COMP) {
      newpairs = Pairpool_push_existing(newpairs,pairpool,pair);
    }
    pairs = pairs->rest;
  }
  return List_reverse(newpairs);
}

static List_T
peel_forward (List_T *peeled_path, List_T path, int *querypos, int *genomepos, 
	      Pairpool_T pairpool, int maxpeelback, bool intronp) {
  Pair_T pair, firstpair;
  int npeelback = 0, orig_querypos, orig_genomepos;
  bool exonp = true, skipp = true, gapp = false;
#ifdef PMAP
  bool midcodonp = false;
#endif
  List_T peeled = NULL;

  orig_querypos = *querypos;
  orig_genomepos = *genomepos;

  /* Peelback */
  debug(printf("Peeling forward:"));
  if (path != NULL) {
    firstpair = path->first;
    *querypos = firstpair->querypos + 1;
    *genomepos = firstpair->genomepos + 1;
  }

  while (npeelback < maxpeelback && path != NULL && exonp) {
    path = Pairpool_pop(path,&pair);

    if (pair->gapp) {
      debug(printf(" [%d %d %c %c! %c]",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      exonp = false;
      path = Pairpool_push_existing(path,pairpool,pair);
    } else if (*querypos - pair->querypos > 1 || *genomepos - pair->genomepos > 1) {
      debug(printf(" [%d! %d! %c %c %c]",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      exonp = false;
      path = Pairpool_push_existing(path,pairpool,pair);
    } else {
      *querypos = pair->querypos;
      *genomepos = pair->genomepos;
      npeelback++;
      debug(printf(" (%d %d %c %c %c)",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
    }
  }
  
  debug(printf(" ||"));
  if (intronp) {
    /* Continue to peelback through little skips and mismatches */
    while (path != NULL && skipp) {
      path = Pairpool_pop(path,&pair);

      if (pair->comp != INDEL_COMP && pair->comp != MISMATCH_COMP) {
	debug(printf(" [%d %d %c %c! %c]",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	skipp = false;
	path = Pairpool_push_existing(path,pairpool,pair);
      } else if (*querypos - pair->querypos > 1 || *genomepos - pair->genomepos > 1) {
	debug(printf(" [%d! %d! %c %c %c]",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	skipp = false;
	path = Pairpool_push_existing(path,pairpool,pair);
      } else {
	*querypos = pair->querypos;
	*genomepos = pair->genomepos;
	debug(printf(" (%d %d %c %c %c)",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	peeled = Pairpool_push_existing(peeled,pairpool,pair);
      }
    }
  }

  debug(printf(" ||"));

#ifdef PMAP
  /* Last pair popped must be beginning of codon */
  if (peeled != NULL) {
    firstpair = peeled->first;
    if (firstpair->querypos % 3 != 0) {
      midcodonp = true;
    }
  }

  while (peeled != NULL && midcodonp) {
    peeled = Pairpool_pop(peeled,&pair);

    if (pair->querypos % 3 == 0) {
      debug(printf(" (%d %d %c %c+ %c)",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      midcodonp = false;
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
    } else {
      debug(printf(" [%d %d %c %c! %c]",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      path = Pairpool_push_existing(path,pairpool,pair);
    }
  }      
#else
  /* Last pair popped cannot be a gap */
  if (peeled != NULL) {
    firstpair = peeled->first;
    if (firstpair->comp == INDEL_COMP) {
      gapp = true;
    }
  }

  while (peeled != NULL && gapp) {
    peeled = Pairpool_pop(peeled,&pair);

    if (pair->comp != INDEL_COMP) {
      debug(printf(" (%d %d %c %c+ %c)",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      gapp = false;
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
    } else {
      debug(printf(" [%d %d %c %c! %c]",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      path = Pairpool_push_existing(path,pairpool,pair);
    }
  }      
#endif

  if (peeled == NULL) {
    *querypos = orig_querypos + 1;
    *genomepos = orig_genomepos + 1;
    debug(printf(" (start at %d %d)",*querypos,*genomepos));
  } else {
    firstpair = peeled->first;
    *querypos = firstpair->querypos;
    *genomepos = firstpair->genomepos;
  }
  debug(printf("\n"));

  debug(
	if (path == NULL) {
	  printf("Top of path is NULL\n");
	} else {
	  pair = path->first;
	  printf("Top of path is %d %d %c %c %c\n",
		 pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	}
	);

  *peeled_path = peeled;
  return path;
}

static List_T
peel_back (List_T *peeled_pairs, List_T pairs, int *lastquerypos, int *lastgenomepos, 
	   Pairpool_T pairpool, int maxpeelback, bool intronp) {
  Pair_T pair, firstpair;
  int npeelback = 0, orig_lastquerypos, orig_lastgenomepos;
  bool exonp = true, skipp = true, gapp = false;
#ifdef PMAP
  bool midcodonp = false;
#endif
  List_T peeled = NULL;

  orig_lastquerypos = *lastquerypos;
  orig_lastgenomepos = *lastgenomepos;

  /* Peelback */
  debug(printf("Peeling back:"));
  if (pairs != NULL) {
    firstpair = pairs->first;
    *lastquerypos = firstpair->querypos - 1;
    *lastgenomepos = firstpair->genomepos - 1;
  }

  while (npeelback < maxpeelback && pairs != NULL && exonp) {
    pairs = Pairpool_pop(pairs,&pair);

    if (pair->gapp) {
      debug(printf(" [%d %d %c %c! %c]",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      exonp = false;
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
    } else if (pair->querypos - *lastquerypos > 1 || pair->genomepos - *lastgenomepos > 1) {
      debug(printf(" [%d! %d! %c %c! %c]",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      exonp = false;
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
    } else {
      *lastquerypos = pair->querypos;
      *lastgenomepos = pair->genomepos;
      npeelback++;
      debug(printf(" (%d %d %c %c %c)",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
    }
  }
  
  debug(printf(" ||"));
  if (intronp) {
    /* Continue to peelback through little skips and mismatches */
    while (pairs != NULL && skipp) {
      pairs = Pairpool_pop(pairs,&pair);

      if (pair->comp != INDEL_COMP && pair->comp != MISMATCH_COMP) {
	debug(printf(" [%d %d %c %cok %c]",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	skipp = false;
	if (pair->gapp == false) {
	  /* Allow this one */
	  peeled = Pairpool_push_existing(peeled,pairpool,pair);
	} else {
	  pairs = Pairpool_push_existing(pairs,pairpool,pair);
	}
      } else if (pair->querypos - *lastquerypos > 1 || pair->genomepos - *lastgenomepos > 1) {
	debug(printf(" [%d! %d! %c %c %c]",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	skipp = false;
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
      } else {
	*lastquerypos = pair->querypos;
	*lastgenomepos = pair->genomepos;
	debug(printf(" (%d %d %c %c %c)",
		     pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
	peeled = Pairpool_push_existing(peeled,pairpool,pair);
      }
    }
  }

  debug(printf(" ||"));
#ifdef PMAP
  /* Last pair popped cannot be a gap */
  if (peeled != NULL) {
    firstpair = peeled->first;
    if (firstpair->querypos % 3 != 2) {
      midcodonp = true;
    }
  }

  while (peeled != NULL && midcodonp) {
    peeled = Pairpool_pop(peeled,&pair);

    if (pair->querypos % 3 == 2) {
      debug(printf(" (%d %d %c %c+ %c)",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      midcodonp = false;
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
    } else {
      debug(printf(" [%d %d %c %c! %c]",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
    }
  }      
#else
  /* Last pair popped cannot be a gap */
  if (peeled != NULL) {
    firstpair = peeled->first;
    if (firstpair->comp == INDEL_COMP) {
      gapp = true;
    }
  }

  while (peeled != NULL && gapp) {
    peeled = Pairpool_pop(peeled,&pair);

    if (pair->comp != INDEL_COMP) {
      debug(printf(" (%d %d %c %c+ %c)",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      gapp = false;
      peeled = Pairpool_push_existing(peeled,pairpool,pair);
    } else {
      debug(printf(" [%d %d %c %c! %c]",
		   pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome));
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
    }
  }      
#endif

  if (peeled == NULL) {
    *lastquerypos = orig_lastquerypos - 1;
    *lastgenomepos = orig_lastgenomepos - 1;
    debug(printf(" (start at %d %d)",*lastquerypos,*lastgenomepos));
  } else {
    firstpair = peeled->first;
    *lastquerypos = firstpair->querypos;
    *lastgenomepos = firstpair->genomepos;
  }
  debug(printf("\n"));

  debug(
	if (pairs == NULL) {
	  printf("Top of pairs is NULL\n");
	} else {
	  pair = pairs->first;
	  printf("Top of pairs is %d %d %c %c %c\n",
		 pair->querypos,pair->genomepos,pair->cdna,pair->comp,pair->genome);
	}
	);

  *peeled_pairs = peeled;
  return pairs;
}


static List_T
traverse_single_gap (bool *filledp, List_T pairs, List_T *path, 
		     int *querypos, int *genomepos, int lastquerypos, int lastgenomepos,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     Pairpool_T pairpool, Dynprog_T dynprog,
		     int maxpeelback, int extraband_single, double defect_rate, bool forcep) {
  List_T gappairs, peeled_pairs, peeled_path;
  int queryjump, genomejump;
  int querydp5, querydp3, genomedp5, genomedp3;
  int nmatches, nmismatches, nopens, nindels;
  int finalscore;

  /* Note that we peelback only half as much as for a paired gap, to
     save on dynamic programming */
  pairs = peel_back(&peeled_pairs,pairs,&lastquerypos,&lastgenomepos,pairpool,maxpeelback/2,
		    /*intronp*/false);
  *path = peel_forward(&peeled_path,*path,&(*querypos),&(*genomepos),pairpool,maxpeelback/2,
		       /*intronp*/false);

  /* Get dynamic programming indices, which are inclusive. */
  querydp5 = *querypos;		/* querypos is last pair popped */
  genomedp5 = *genomepos;
  querydp3 = lastquerypos;
  genomedp3 = lastgenomepos;

  queryjump = querydp3 - querydp5 + 1;
  genomejump = genomedp3 - genomedp5 + 1;
  
  if (lastquerypos == *querypos || lastgenomepos == *genomepos) {
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
				  pairpool,extraband_single,defect_rate);
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
traverse_cdna_gap (List_T pairs, List_T *path,
		   int *querypos, int *genomepos, int lastquerypos, int lastgenomepos,
#ifdef PMAP
		   char *queryaaseq_ptr,
#endif
		   char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		   int cdna_direction, Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogR, 
		   int maxpeelback, int extramaterial_paired, int extraband_paired, double defect_rate) {
  List_T gappairs, peeled_pairs, peeled_path;
  int querydp5, querydp3, genomedp5, genomedp3, queryjump, genomejump;
  int finalscore;

  pairs = peel_back(&peeled_pairs,pairs,&lastquerypos,&lastgenomepos,pairpool,maxpeelback,
		    /*intronp*/true);
  *path = peel_forward(&peeled_path,*path,&(*querypos),&(*genomepos),pairpool,maxpeelback,
		       /*intronp*/true);

  /* Get dynamic programming indices, which are inclusive. */
  querydp5 = *querypos;
  genomedp5 = *genomepos;
  querydp3 = lastquerypos;
  genomedp3 = lastgenomepos;

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
			      pairpool,extraband_paired,defect_rate);
  debug(Pair_dump_list(gappairs,true));
  pairs = Pairpool_transfer(pairs,gappairs);
  return pairs;
}


static List_T
traverse_genome_gap (bool *filledp, int *nintrons, int *nnonintrons, int *intronlen, int *nonintronlen, 
		     List_T pairs, List_T *path,
		     int *querypos, int *genomepos, int lastquerypos, int lastgenomepos,
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, int ngap, Pairpool_T pairpool,
		     Dynprog_T dynprogL, Dynprog_T dynprogR, int maxpeelback,
		     int extramaterial_paired, int extraband_paired, double defect_rate, bool forcep) {
  List_T gappairs, peeled_pairs, peeled_path, micropairs;
  int querydp5, querydp3, genomedp5, genomedp3, queryjump, genomejump;
  int finalscore, nmatches, nmismatches, nopens, nindels, exonhead, introntype, microintrontype;
  int acceptable_nmismatches;

  pairs = peel_back(&peeled_pairs,pairs,&lastquerypos,&lastgenomepos,pairpool,maxpeelback,
		    /*intronp*/true);
  *path = peel_forward(&peeled_path,*path,&(*querypos),&(*genomepos),pairpool,maxpeelback,
		       /*intronp*/true);

  /* Get dynamic programming indices, which are inclusive. */
  querydp5 = *querypos;
  genomedp5 = *genomepos;
  querydp3 = lastquerypos;
  genomedp3 = lastgenomepos;

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
				cdna_direction,ngap,pairpool,extraband_paired,false,
				defect_rate,/*returnpairsp*/true,/*addgapp*/true);
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
				       queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,pairpool,ngap);
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
    *nonintronlen += lastgenomepos - *genomepos;
  } else {
    *nintrons += 1;
    *intronlen += lastgenomepos - *genomepos;
  }

  return pairs;
}


static List_T
traverse_dual_genome_gap (bool *singlep, List_T pairs, List_T *path, 
			  int *querypos, int *genomepos, int lastquerypos, int lastgenomepos, 
			  int midquerypos, int midgenomepos, 
#ifdef PMAP
			  char *queryaaseq_ptr,
#endif
			  char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
			  int cdna_direction, int ngap,
			  Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, 
			  int maxpeelback, int nullgap, int extramaterial_paired, int extraband_paired,
			  double defect_rate) {
  List_T single_gappairs, dual_gappairs_1, dual_gappairs_2, peeled_pairs, peeled_path;
  int querydp5, querydp3, genomedp5, genomedp3, queryjump, genomejump;
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

  pairs = peel_back(&peeled_pairs,pairs,&lastquerypos,&lastgenomepos,pairpool,maxpeelback,
		    /*intronp*/false);
  *path = peel_forward(&peeled_path,*path,&(*querypos),&(*genomepos),pairpool,maxpeelback,
		       /*intronp*/false);

  /* No short exon */
  querydp5 = *querypos;		/* querypos is last pair popped */
  genomedp5 = *genomepos;
  querydp3 = lastquerypos;
  genomedp3 = lastgenomepos;

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
				       cdna_direction,ngap,pairpool,extraband_paired,false,
				       defect_rate,/*returnpairsp*/true,/*addgapp*/false);

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
  querydp3 = lastquerypos;
  genomedp3 = lastgenomepos;

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
				       cdna_direction,ngap,pairpool,extraband_paired,false,
				       defect_rate,/*returnpairsp*/false,/*addgapp*/false);

  dual_goodness = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

  /* Left of short exon */
  querydp5 = *querypos;
  genomedp5 = *genomepos;
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
				       cdna_direction,ngap,pairpool,extraband_paired,false,
				       defect_rate,/*returnpairsp*/false,/*addgapp*/false);

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
	       right_exonhead,left_exonhead,right_exonhead-left_exonhead,lastgenomepos-*genomepos));
  if (middle_exonlength <= 0) {
    middle_exonprob = 0.0;
  } else {
    interexon_region = lastgenomepos-*genomepos;

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
		    int querypos, int genomepos, int *lastquerypos, int *lastgenomepos,
#ifdef PMAP
		    char *queryaaseq_ptr,
#endif
		    char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		    int cdna_direction, int introndir, Pairpool_T pairpool,
		    Dynprog_T dynprog, int maxpeelback, int extramaterial_end,
		    int extraband_end, double defect_rate, int ngap, bool extend_mismatch_p,
		    bool end_microexons_p) {
  Pair_T lastpair;
  List_T intron_gappairs = NULL, continuous_gappairs = NULL, peeled_pairs, p;
  int queryjump, genomejump;
  int querydp5, querydp3, genomedp5, genomedp3;
  int microintrontype, microexonlength;
  int intron_goodness, continuous_goodness, nmissed, nmatches, nmismatches, nopens, nindels, acceptable_nmismatches;

  /* Note that we peelback only half as much as for a paired gap, to
     save on dynamic programming */
  *pairs = peel_back(&peeled_pairs,*pairs,&(*lastquerypos),&(*lastgenomepos),pairpool,maxpeelback/2,
		     /*intronp*/false);

  /* Get dynamic programming indices, which are inclusive. */
  querydp5 = querypos + 1;	/* Add 1 because peel_forward was not called */
  genomedp5 = genomepos + 1;
  querydp3 = *lastquerypos;
  genomedp3 = *lastgenomepos;
  
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
					 pairpool,extraband_end,defect_rate,cdna_direction,ngap,extend_mismatch_p);
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
					  pairpool,ngap,end_microexons_p);
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
try_ending3_extend (int *finalscore, List_T *pairs, 
		    int *querypos, int *genomepos, int lastquerypos, int lastgenomepos, 
		    int querylength, Genomicpos_T genomiclength,
#ifdef PMAP
		    char *queryaaseq_ptr,
#endif
		    char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		    int cdna_direction, int introndir,
		    Pairpool_T pairpool, Dynprog_T dynprog, int maxpeelback, int extramaterial_end,
		    int extraband_end, double defect_rate, int ngap, bool extend_mismatch_p,
		    bool end_microexons_p) {
  Pair_T firstpair;
  List_T intron_gappairs = NULL, continuous_gappairs = NULL, peeled_pairs;
  int queryjump, genomejump;
  int querydp5, querydp3, genomedp5, genomedp3;
  int microintrontype, microexonlength;
  int intron_goodness, continuous_goodness, nmissed, nmatches, nmismatches, nopens, nindels, acceptable_nmismatches;

  /* Note that we peelback only half as much as for a paired gap, to
     save on dynamic programming */
  *pairs = peel_forward(&peeled_pairs,*pairs,&(*querypos),&(*genomepos),pairpool,maxpeelback/2,
		       /*intronp*/false);

  /* Get dynamic programming indices, which are inclusive. */
  querydp5 = *querypos;
  genomedp5 = *genomepos;
  querydp3 = lastquerypos - 1;	/* Subtract 1 because peel_back was not called */
  genomedp3 = lastgenomepos - 1;
  
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
					 pairpool,extraband_end,defect_rate,cdna_direction,ngap,extend_mismatch_p);
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
					  genomiclength,pairpool,ngap,
					  end_microexons_p);
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

  
/* Note: querypos is actually indexsize nt to the left of the last nt match.
   
        ||||||||********   X  X  XX X   X
               ^         <- queryjump->  ^
           querypos                      lastquerypos         

                <-     querydpspan     ->
*/
static List_T
build_pairs_end3 (List_T pairs, int *querypos, int *genomepos, int querylength, Genomicpos_T genomiclength,
#ifdef PMAP
		  char *queryaaseq_ptr,
#endif
		  char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		  int cdna_direction, int introndir, int maxpeelback, int nullgap,
		  int extramaterial_end, int extramaterial_paired, int extraband_end,
		  double defect_rate, int ngap, Pairpool_T pairpool, Dynprog_T dynprogL, bool extend_mismatch_p,
		  bool end_microexons_p) {
  List_T gappairs;
  int lastquerypos, lastgenomepos;
  int queryjump, genomejump;
  int finalscore;

  lastquerypos = querylength;
  lastgenomepos = genomiclength;
  debug(printf("Stage 3: 3' end: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
	       *querypos,lastquerypos,*genomepos,lastgenomepos));

  queryjump = lastquerypos - *querypos - 1;
  genomejump = lastgenomepos - *genomepos - 1;

  /* Note difference with 5' case.  We use queryjump+1 and genomejump+1 here instead of queryjump and genomejump */
  if (genomejump+1 < 0) {
    /* Past end of genomicseg; do nothing */
  } else {
    if (queryjump+1 > nullgap) {
      lastquerypos = *querypos + nullgap + 1;
    }
    gappairs = try_ending3_extend(&finalscore,&pairs,
				  &(*querypos),&(*genomepos),lastquerypos,lastgenomepos,
				  querylength,genomiclength,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				  cdna_direction,introndir,pairpool,dynprogL,maxpeelback,
				  extramaterial_end,extraband_end,defect_rate,ngap,extend_mismatch_p,
				  end_microexons_p);
    debug(Pair_dump_list(gappairs,true));
    pairs = Pairpool_transfer(pairs,gappairs);
  }

  return pairs;
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
		  double defect_rate, int ngap, Pairpool_T pairpool, Dynprog_T dynprogR, bool extend_mismatch_p,
		  bool end_microexons_p) {
  List_T gappairs;
  Pair_T firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos;
  int queryjump, genomejump;
  int finalscore;

  if (pairs == NULL) {
    return NULL;
  } else {
    firstpair = pairs->first;
    lastquerypos = firstpair->querypos;
    lastgenomepos = firstpair->genomepos;
    querypos = -1;
    genomepos = -1;
  }
  debug(printf("Stage 3: 5' end: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
	       querypos,lastquerypos,genomepos,lastgenomepos));

  queryjump = lastquerypos - querypos - 1;
  genomejump = lastgenomepos - genomepos - 1;

  /* Note difference with 3' case.  We use queryjump here instead of queryjump+1 */
  if (genomejump < 0) {
    /* Past end of genomicseg; do nothing */
  } else {
    if (queryjump > nullgap) {
      querypos = lastquerypos - nullgap - 1;
    }

    gappairs = try_ending5_extend(&finalscore,&pairs,querypos,genomepos,&lastquerypos,&lastgenomepos,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				  cdna_direction,introndir,pairpool,dynprogR,maxpeelback,
				  extramaterial_end,extraband_end,defect_rate,ngap,extend_mismatch_p,
				  end_microexons_p);
    debug(Pair_dump_list(gappairs,true));
    pairs = Pairpool_transfer(pairs,gappairs);

    /* Remove gaps at 5' end.  Shouldn't be needed, because dynamic
       programming should eliminate these. */
    /*
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
    */

  }

  return pairs;
}


static List_T
build_pairs_null (List_T path, Pairpool_T pairpool, int ngap) {
  List_T pairs;
  Pair_T pair, firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump;

  path = Pairpool_pop(path,&firstpair);
  pairs = Pairpool_push_existing(NULL,pairpool,firstpair);
  lastquerypos = firstpair->querypos;
  lastgenomepos = firstpair->genomepos;

  while (path != NULL) {
    pair = path->first;
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - querypos - 1;
    genomejump = lastgenomepos - genomepos - 1;
    if (queryjump <= 0 && genomejump <= 0) {
      /* Do nothing */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    } else {
      pairs = add_null_gap(pairs,lastquerypos,lastgenomepos,pairpool,ngap);
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;
    }
  }

  return pairs;
}  

static List_T
build_pairs_null_complete (List_T path, Pairpool_T pairpool, int ngap) {
  List_T pairs;
  Pair_T pair, firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump;

  path = Pairpool_pop(path,&firstpair);
  pairs = Pairpool_push_existing(NULL,pairpool,firstpair);
  lastquerypos = firstpair->querypos;
  lastgenomepos = firstpair->genomepos;

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
      pairs = add_null_gap(pairs,lastquerypos,lastgenomepos,pairpool,ngap);
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;
    }
  }

  return pairs;
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
		     int extramaterial_paired, int extraband_single, double defect_rate, int ngap,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T pairs;
  Pair_T pair, firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump;
  char c1, c2;
  bool filledp;

  path = Pairpool_pop(path,&firstpair);
  pairs = Pairpool_push_existing(NULL,pairpool,firstpair);
  lastquerypos = firstpair->querypos;
  lastgenomepos = firstpair->genomepos;

  while (path != NULL) {
    pair = path->first;
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - querypos - 1;
    genomejump = lastgenomepos - genomepos - 1;
    if (queryjump == 0 && genomejump == 0) {
      /* Adjacent.  Do nothing.  Note that the condition here on
         initial pass is == 0, rather than <= 0. */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    } else if (queryjump > nullgap) {
      /* Large gap.  Do nothing */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    } else if (queryjump > genomejump + EXTRAQUERYGAP) {
      /* cDNA insertion.  Do nothing */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    } else if (genomejump > queryjump + MININTRONLEN) {
      /* Intron.  Do nothing */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    } else if (queryjump == 1 && genomejump == 1) {
      /* Special case: Single mismatch */
      c1 = queryseq_ptr[lastquerypos-1];
      c2 = genomicseg_ptr[lastgenomepos-1];
      if (queryuc_ptr[lastquerypos-1] == genomicuc_ptr[lastgenomepos-1]) {
	/* Possible for a stretch of the same nucleotide to have a mismatch */
	debug(printf("Stage 3: Single mismatch at %d %d: %c | %c\n",
		     lastquerypos-1,lastgenomepos-1,c1,c2));
	pairs = Pairpool_push(pairs,pairpool,lastquerypos-1,lastgenomepos-1,
			      c1,MATCH_COMP,c2);
      } else {
	debug(printf("Stage 3: Single mismatch at %d %d: %c   %c\n",
		     lastquerypos-1,lastgenomepos-1,c1,c2));
	pairs = Pairpool_push(pairs,pairpool,lastquerypos-1,lastgenomepos-1,
			      c1,MISMATCH_COMP,c2);
      }
      pair = pairs->first;
      lastquerypos = pair->querypos;
      lastgenomepos = pair->genomepos;

    } else if (queryjump == 1 && genomejump == 0) {
      /* Special case: Single cDNA insert */
      debug(printf("Stage 3: Single cDNA insert at %d %d: %c\n",
		   lastquerypos-1,lastgenomepos,queryseq_ptr[lastquerypos-1]));
      pairs = Pairpool_push(pairs,pairpool,lastquerypos-1,lastgenomepos,
			    queryseq_ptr[lastquerypos-1],INDEL_COMP,' ');
      pair = pairs->first;
      lastquerypos = pair->querypos;
      lastgenomepos = pair->genomepos;

    } else {
      /* Guarantees: queryjump <= nullgap && genomejump < queryjump - EXTRAQUERYGAP &&
	 genomejump <= queryjump + MININTRONLEN, meaning that score matrix is nearly square */
      debug(printf("Stage 3: Traversing single gap: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		   querypos,lastquerypos,genomepos,lastgenomepos));
      pairs = traverse_single_gap(&filledp,pairs,&path,&querypos,&genomepos,lastquerypos,lastgenomepos,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				  pairpool,dynprogM,maxpeelback,extraband_single,defect_rate,
				  /*forcep*/false);
      if (filledp == true && pairs != NULL) {
	pair = pairs->first;
	lastquerypos = pair->querypos;
	lastgenomepos = pair->genomepos;
      } else {
	path = Pairpool_pop(path,&pair);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
	lastquerypos = querypos;
	lastgenomepos = genomepos;
      }
    }
  }

  return pairs;
}  

static void
find_middle_exon_bounds (int *midquerypos_start, int *midgenomepos_start,
			 int *midquerypos_end, int *midgenomepos_end, List_T midexon_pairs) {
  List_T p;
  Pair_T pair;
  int queryjump2, genomejump2, prevquerypos, prevgenomepos;
  int longesti = -1, i, longest = 0, currlength;
  bool newexonp = true;

  debug(Pair_dump_list(midexon_pairs,true));

  i = 0;
  currlength = 0;
  p = midexon_pairs;
  while (p != NULL) {
    pair = p->first;

    if (newexonp == true) {
      /* printf("Starting new exon\n"); */
      currlength++;
      prevquerypos = pair->querypos;
      prevgenomepos = pair->genomepos;
      newexonp = false;
    } else {
      queryjump2 = prevquerypos - pair->querypos + 1;
      genomejump2 = prevgenomepos - pair->genomepos + 1;
      /* printf("%d %d: ",queryjump2,genomejump2); */
      if (queryjump2 >= 0 && genomejump2 >= 0) {
	/* printf("Continuing in exon\n"); */
	currlength++;
	prevquerypos = pair->querypos;
	prevgenomepos = pair->genomepos;

      } else {
	/* printf("Ending exon\n"); */
	if (currlength > longest) {
	  longesti = i;
	  longest = currlength;
	}
	i++;
	currlength = 0;
	newexonp = true;
      }
    }
    p = List_next(p);
  }
  if (currlength > longest) {
    longesti = i;
    longest = currlength;
  }
  
  i = 0;
  p = midexon_pairs;
  while (p != NULL) {
    pair = p->first;

    if (newexonp == true) {
      /* printf("Starting new exon\n"); */
      if (i == longesti) {
	*midquerypos_start = pair->querypos - 1;
	*midgenomepos_start = pair->genomepos - 1;
      }
      prevquerypos = pair->querypos;
      prevgenomepos = pair->genomepos;
      newexonp = false;
    } else {
      queryjump2 = prevquerypos - pair->querypos + 1;
      genomejump2 = prevgenomepos - pair->genomepos + 1;
      /* printf("%d %d: ",queryjump2,genomejump2); */
      if (queryjump2 >= 0 && genomejump2 >= 0) {
	/* printf("Continuing in exon\n"); */
	prevquerypos = pair->querypos;
	prevgenomepos = pair->genomepos;

      } else {
	/* printf("Ending exon\n"); */
	if (i == longesti) {
	  *midquerypos_end = pair->querypos;
	  *midgenomepos_end = pair->genomepos;
	}
	i++;
	newexonp = true;
      }
    }
    p = List_next(p);
  }
  if (i == longesti) {
    *midquerypos_end = pair->querypos;
    *midgenomepos_end = pair->genomepos;
  }
  
  debug(printf("Longest one is %d with length %d\n",longesti,longest));
  debug(printf("Bounds are %d..%d and %d..%d\n\n",
	       *midquerypos_start,*midquerypos_end,*midgenomepos_start,*midgenomepos_end));
  return;
}


static List_T
build_pairs_dualintrons (bool *singlep, List_T path, int nshortexons,
#ifdef PMAP
			 char *queryaaseq_ptr,
#endif
			 char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
			 int cdna_direction, int maxpeelback, int nullgap,
			 int extramaterial_paired, int extraband_paired, double defect_rate, int ngap,
			 Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T pairs, midexon_pairs;
  Pair_T pair, firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump;
  int queryjump2, genomejump2;
  int midquerypos, midgenomepos, midgenomepos_start, midgenomepos_end, prevquerypos, prevgenomepos,
    midquerypos_start;
  bool single_one = false, exonp;

  *singlep = false;

  path = Pairpool_pop(path,&firstpair);
  pairs = Pairpool_push_existing(NULL,pairpool,firstpair);
  lastquerypos = firstpair->querypos;
  lastgenomepos = firstpair->genomepos;

  while (path != NULL) {
    pair = path->first;
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - querypos - 1;
    genomejump = lastgenomepos - genomepos - 1;
    if (pair->gapp == true) {
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      /* Don't reset pointers */
      /*
      lastquerypos = querypos + 1;
      lastgenomepos = genomepos + 1;
      */

    } else if (queryjump <= 0 && genomejump <= 0) {
      /* Do nothing */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    } else if (queryjump > nullgap) {
      /* Large gap; do nothing */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    } else if (queryjump > genomejump + EXTRAQUERYGAP) {
      /* cDNA gap; do nothing */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    } else if (genomejump > queryjump + MININTRONLEN) {
      if (pair->shortexonp == true) {
	debug(printf("I see a short exon...crossing\n"));
	/* Short exon */
	path = Pairpool_pop(path,&pair);
	midexon_pairs = Pairpool_push_existing(NULL,pairpool,pair);

	prevquerypos = querypos;
	prevgenomepos = genomepos;

#ifdef CROSS_ONE_SHORT
	midquerypos_start = querypos;
	midgenomepos_start = genomepos;

	/* Cross one short exon */
	exonp = true;
	while (path != NULL && exonp) {
	  pair = path->first;

	  queryjump2 = prevquerypos - pair->querypos - 1;
	  genomejump2 = prevgenomepos - pair->genomepos - 1;
	  if (pair->gapp == true) {
	    exonp = false;
	  } else if (queryjump2 <= 0 && genomejump2 <= 0) {
	    path = Pairpool_pop(path,&pair);
	    midexon_pairs = Pairpool_push_existing(midexon_pairs,pairpool,pair);
	    prevquerypos = pair->querypos;
	    prevgenomepos = pair->genomepos;
	  } else {
	    exonp = false;
	  }
	}
	querypos = pair->querypos;
	genomepos = pair->genomepos;

	debug(printf("Finished crossing a short exon\n"));
#else
	/* Cross all short exons */
	shortp = true;
	while (path != NULL && shortp) {
	  pair = path->first;
	  if ((shortp = pair->shortexonp) == true) {
	    path = Pairpool_pop(path,&pair);
	    midexon_pairs = Pairpool_push_existing(midexon_pairs,pairpool,pair);
	  }
	}
	querypos = pair->querypos;
	genomepos = pair->genomepos;

	/* Find bounds of longest short exon */
	find_middle_exon_bounds(&midquerypos_start,&midgenomepos_start,
				&midquerypos_end,&midgenomepos_end,midexon_pairs);
#endif

	if (path == NULL) {
	  /* At beginning of query sequence; do nothing */

	} else {
	  /* Perform dual intron gap */

#ifdef CROSS_ONE_SHORT
	  midgenomepos_end = prevgenomepos;
	  midgenomepos = (midgenomepos_start + midgenomepos_end)/2;
	  midquerypos = midquerypos_start - (midgenomepos_start - midgenomepos);
#else
	  midgenomepos = (midgenomepos_start + midgenomepos_end)/2;
	  midquerypos = (midquerypos_start + midquerypos_end)/2;
#endif
	  
	  debug(printf("Stage 3: Traversing dual intron gap: querypos = %d, midquerypos = %d, lastquerypos = %d, genomepos = %d, midgenomepos = %d, lastgenomepos = %d\n",
		       querypos,midquerypos,lastquerypos,genomepos,midgenomepos,lastgenomepos));
	  
	  pairs = traverse_dual_genome_gap(&single_one,pairs,&path,&querypos,&genomepos,lastquerypos,lastgenomepos,
					   midquerypos,midgenomepos,
#ifdef PMAP
					   queryaaseq_ptr,
#endif
					   queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
					   cdna_direction,ngap,pairpool,dynprogL,dynprogM,dynprogR,
					   maxpeelback,nullgap,extramaterial_paired,extraband_paired,
					   defect_rate);
	  *singlep |= single_one;
	  if (single_one == true) {
	    pair = pairs->first;
	    lastquerypos = pair->querypos;
	    lastgenomepos = pair->genomepos;

	  } else {
	    /* Dual intron kept; continue after first intron */
	    debug(printf("Beginning transfer of middle exon to path:\n"));
	    path = Pairpool_transfer(path,midexon_pairs);
	    debug(printf("Done with transfer of middle exon to path\n"));

	    path = Pairpool_pop(path,&pair);
	    pairs = Pairpool_push_existing(pairs,pairpool,pair);
	    lastquerypos = pair->querypos;
	    lastgenomepos = pair->genomepos;
	  }
	}
      } else {
	/* Long exon; do nothing */
	path = Pairpool_pop(path,&pair);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
	lastquerypos = querypos;
	lastgenomepos = genomepos;

      }
    } else {
      /* Single gap; Do nothing */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    }
  }

  return pairs;
}  


static List_T
build_pairs_introns (int *nintrons, int *nnonintrons, int *intronlen, int *nonintronlen, 
		     List_T pairs, int lastquerypos, int lastgenomepos, List_T path, 
#ifdef PMAP
		     char *queryaaseq_ptr,
#endif
		     char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
		     int cdna_direction, int maxpeelback, int nullgap, int extramaterial_paired, 
		     int extraband_single, int extraband_paired, double defect_rate, int ngap,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  Pair_T pair;
  int querypos, genomepos, queryjump, genomejump;
  bool filledp;

  while (path != NULL) {
    pair = path->first;
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - querypos - 1;
    genomejump = lastgenomepos - genomepos - 1;

    if (pair->gapp == true) {
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      /* Don't reset pointers.  If you do, the +1 is critical here. */
      /*
      lastquerypos = querypos + 1;
      lastgenomepos = genomepos + 1;
      */

    } else if (queryjump <= 0 && genomejump <= 0) {
      /* Do nothing */
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;
      
    } else if (queryjump > nullgap) {
      debug(printf("Stage 3: Adding large gap: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		   querypos,lastquerypos,genomepos,lastgenomepos));
      pairs = add_null_gap(pairs,lastquerypos,lastgenomepos,pairpool,ngap);
      path = Pairpool_pop(path,&pair);
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
      lastquerypos = querypos;
      lastgenomepos = genomepos;

    } else if (queryjump > genomejump + EXTRAQUERYGAP) {
      debug(printf("Stage 3: Traversing cDNA gap: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		   querypos,lastquerypos,genomepos,lastgenomepos));
      pairs = traverse_cdna_gap(pairs,&path,&querypos,&genomepos,lastquerypos,lastgenomepos,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				pairpool,dynprogL,dynprogR,maxpeelback,extramaterial_paired,
				extraband_paired,defect_rate);
      pair = pairs->first;
      lastquerypos = pair->querypos;
      lastgenomepos = pair->genomepos;
      
    } else if (genomejump > queryjump + 2*MININTRONLEN) {
      /* Need space for two introns */
      /* We will make the score matrices nearly square */
      debug(printf("Stage 3: Traversing paired gap: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		   querypos,lastquerypos,genomepos,lastgenomepos));
      pairs = traverse_genome_gap(&filledp,&(*nintrons),&(*nnonintrons),&(*intronlen),&(*nonintronlen),
				  pairs,&path,&querypos,&genomepos,lastquerypos,lastgenomepos,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,ngap,
				  pairpool,dynprogL,dynprogR,maxpeelback,extramaterial_paired,
				  extraband_paired,defect_rate,	/*forcep*/true);
      /* Do forcep, because adding large gap is not a good solution */

      if (filledp == true && pairs != NULL) {
	pair = pairs->first;
	lastquerypos = pair->querypos;
	lastgenomepos = pair->genomepos;
      } else {
	querypos = pair->querypos;
	genomepos = pair->genomepos;
	
	debug(printf("Stage 3: Adding large gap instead: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		     querypos,lastquerypos,genomepos,lastgenomepos));
	pairs = add_null_gap(pairs,lastquerypos,lastgenomepos,pairpool,ngap);
	path = Pairpool_pop(path,&pair);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
	lastquerypos = querypos;
	lastgenomepos = genomepos;
      }

    } else {
      /* Single gap; force fill */
      debug(printf("Stage 3: Traversing single gap: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d.  queryjump = %d, genomejump = %d\n",
		   querypos,lastquerypos,genomepos,lastgenomepos,queryjump,genomejump));
      pairs = traverse_single_gap(&filledp,pairs,&path,&querypos,&genomepos,lastquerypos,lastgenomepos,
#ifdef PMAP
				  queryaaseq_ptr,
#endif
				  queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
				  pairpool,dynprogM,maxpeelback,extraband_single,defect_rate,
				  /*forcep*/true);
      if (pairs == NULL) {
	path = Pairpool_pop(path,&pair);
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
	lastquerypos = querypos;
	lastgenomepos = genomepos;

      } else {
	pair = pairs->first;
	lastquerypos = pair->querypos;
	lastgenomepos = pair->genomepos;
      }
    }
  }

  return pairs;
}  


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
	      int querylength, Genomicpos_T genomiclength,
#ifdef PMAP
	      char *queryaaseq_ptr,
#endif
	      char *queryseq_ptr, char *queryuc_ptr, char *genomicseg_ptr, char *genomicuc_ptr,
	      int maxpeelback, int nullgap,
	      int extramaterial_end, int extramaterial_paired,
	      int extraband_single, int extraband_end, int extraband_paired, double defect_rate, 
	      int ngap, Pairpool_T pairpool, 
	      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, bool extend_mismatch_p,
	      bool end_microexons_p) {
  List_T pairs = NULL;
  Pair_T firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, nshortexons, iter;
  int newintrondir, nintrons = 0, nnonintrons = 0;
  bool singlep;

  if (path == NULL) {
    return NULL;
  }

  debug3(pairs = build_pairs_null_complete(path,pairpool,ngap));
  debug3(return pairs);

#ifdef PMAP
  /* Pass 0: undefine nucleotides around gaps */
  pairs = undefine_nucleotides(queryseq_ptr,querylength,path,pairpool,/*width*/6);
  path = List_reverse(pairs);
#endif

  /* Pass 1: solve single gaps */
  pairs = build_pairs_singles(path,
#ifdef PMAP
			      queryaaseq_ptr,
#endif
			      queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
			      maxpeelback,nullgap,extramaterial_paired,extraband_single,
			      defect_rate,ngap,pairpool,dynprogL,dynprogM,dynprogR);
  debug4(path = List_reverse(pairs));
  debug4(pairs = build_pairs_null(path,pairpool,ngap));
  debug4(return pairs);

  /* Pass 2: solve dual introns */
  debug6(printf("Iteration 0:\n"));
  debug6(path = List_reverse(pairs));
  debug6(pairs = build_pairs_null(path,pairpool,ngap));
  debug6(Pair_debug_alignment(pairs,ngap));
  debug6(pairs = remove_null_gaps(pairs,pairpool));

  pairs = Smooth_pairs(&nshortexons,pairs,pairpool);
  /* Not necessary because no introns have been added yet */
  /* pairs = remove_end5_gap(pairs,pairpool); */

  debug5(path = List_reverse(pairs));
  debug5(pairs = build_pairs_null(path,pairpool,ngap));
  debug5(return pairs);

  if (nshortexons > 0) {
    path = List_reverse(pairs);
    pairs = build_pairs_dualintrons(&singlep,path,nshortexons,
#ifdef PMAP
				    queryaaseq_ptr,
#endif
				    queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				    maxpeelback,nullgap,extramaterial_paired,extraband_paired,
				    defect_rate,ngap,pairpool,dynprogL,dynprogM,dynprogR);
  } else {
    singlep = false;
  }

  iter = 0;
  while (singlep == true && ++iter < MAXITER) {
    debug6(printf("Iteration %d:\n",iter));
    debug6(path = List_reverse(pairs));
    debug6(pairs = build_pairs_null(path,pairpool,ngap));
    debug6(Pair_debug_alignment(pairs,ngap));
    debug6(pairs = remove_null_gaps(pairs,pairpool));

    Smooth_reset(pairs);
    pairs = Smooth_pairs(&nshortexons,pairs,pairpool);
    pairs = remove_end5_gap(pairs,pairpool);

    if (nshortexons > 0) {
      path = List_reverse(pairs);
      pairs = build_pairs_dualintrons(&singlep,path,nshortexons,
#ifdef PMAP
				      queryaaseq_ptr,
#endif
				      queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				      maxpeelback,nullgap,extramaterial_paired,extraband_paired,
				      defect_rate,ngap,pairpool,dynprogL,dynprogM,dynprogR);

    } else {
      singlep = false;
    }
  }
  debug6(printf("Done with iterations\n"));
  debug6(path = List_reverse(pairs));
  debug6(pairs = build_pairs_null(path,pairpool,ngap));
  debug6(Pair_debug_alignment(pairs,ngap));
  debug6(return pairs);

  /* Pass 3: solve introns and ends */
  if (pairs != NULL) {
    path = List_reverse(pairs);

    firstpair = path->first;
    lastquerypos = firstpair->querypos;
    lastgenomepos = firstpair->genomepos;

    pairs = build_pairs_introns(&nintrons,&nnonintrons,&(*intronlen),&(*nonintronlen),
				NULL,lastquerypos,lastgenomepos,path,
#ifdef PMAP
				queryaaseq_ptr,
#endif
				queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,cdna_direction,
				maxpeelback,nullgap,extramaterial_paired,extraband_single,extraband_paired,
				defect_rate,ngap,pairpool,dynprogL,dynprogM,dynprogR);

    debug(printf("nintrons = %d, nnonintrons = %d\n",nintrons,nnonintrons));
    debug(printf("Intronlen = %d, Nonintronlen = %d\n",*intronlen,*nonintronlen));

    if (nintrons > nnonintrons) {
      newintrondir = cdna_direction;
    } else {
      newintrondir = introndir;
    }

    pairs = build_pairs_end5(pairs,
#ifdef PMAP
			     queryaaseq_ptr,
#endif
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			     cdna_direction,newintrondir,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_end,
			     defect_rate,ngap,pairpool,dynprogR,extend_mismatch_p,
			     end_microexons_p);

    path = List_reverse(pairs);
    pairs = path;

    firstpair = pairs->first;
    querypos = firstpair->querypos;
    genomepos = firstpair->genomepos;

    pairs = build_pairs_end3(pairs,&querypos,&genomepos,querylength,genomiclength,
#ifdef PMAP
			     queryaaseq_ptr,
#endif
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,
			     cdna_direction,newintrondir,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_end,
			     defect_rate,ngap,pairpool,dynprogL,extend_mismatch_p,
			     end_microexons_p);
  }

  return List_reverse(pairs);
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
		int extraband_single, int extraband_end, int extraband_paired, int ngap,
		bool extend_mismatch_p, bool end_microexons_p, Pairpool_T pairpool, 
		Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		IIT_T altstrain_iit) {
  T this;
  List_T pairs_fwd, pairs_rev;
  double defect_rate;
  int querylength;
  int nfwdintrons, nrevintrons, nunkintrons, introndir;
  int fwd_intronlen = 0, rev_intronlen = 0;
  int fwd_nonintronlen = 0, rev_nonintronlen = 0;
  bool fwdp, revp;
  Genomicpos_T genomiclength;
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
  fwdp = true;
  revp = false;
  introndir = +1;
#else
  if (nfwdintrons == 0 && nrevintrons == 0) {
    /* Should try both even if no introns (cf, AA011563) */
    fwdp = revp = true;
    introndir = 0;
  } else if (nfwdintrons > nrevintrons + nunkintrons) {
    /* Forward wins */
    fwdp = true;
    revp = false;
    introndir = +1;
  } else if (nrevintrons > nfwdintrons + nunkintrons) {
    /* Reverse wins */
    fwdp = false;
    revp = true;
    introndir = -1;
  } else {
    /* Tie */
    fwdp = revp = true;
    introndir = 0;
  }
#endif
  debug(printf("Introns: %d fwd, %d rev, %d unk => fwdp = %d, revp = %d\n",
	       nfwdintrons, nrevintrons, nunkintrons, fwdp, revp));

  if (fwdp == true) {
    pairs_fwd = path_compute(&fwd_intronlen,&fwd_nonintronlen,Stage2_path(stage2),+1,introndir,
			     querylength,genomiclength,
#ifdef PMAP
			     queryaaseq_ptr,
#endif
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_single,extraband_end,extraband_paired,
			     defect_rate,ngap,pairpool,dynprogL,dynprogM,dynprogR,extend_mismatch_p,
			     end_microexons_p);
  } else {
    pairs_fwd = NULL;
  }
  if (revp == true) {
    pairs_rev = path_compute(&rev_intronlen,&rev_nonintronlen,Stage2_path(stage2),-1,introndir,
			     querylength,genomiclength,
#ifdef PMAP
			     queryaaseq_ptr,
#endif
			     queryseq_ptr,queryuc_ptr,genomicseg_ptr,genomicuc_ptr,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_single,extraband_end,extraband_paired,
			     defect_rate,ngap,pairpool,dynprogL,dynprogM,dynprogR,extend_mismatch_p,
			     end_microexons_p);
  } else {
    pairs_rev = NULL;
  }

  this = Stage3_new(pairs_fwd,pairs_rev,matchpairend,straintype,strain,
		    chrnum,chrpos,chroffset,querylength,genomiclength,watsonp,defect_rate,ngap);

  if (this == NULL) {
    return NULL;
  } else if (straintype == 0) {
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

