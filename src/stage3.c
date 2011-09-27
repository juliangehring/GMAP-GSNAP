static char rcsid[] = "$Id: stage3.c,v 1.161 2005/03/09 19:23:05 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For pow() */
#include "mem.h"
#include "pair.h"
#include "pairdef.h"
#include "listdef.h"
#include "smooth.h"
#include "scores.h"

/* The following are the same as dividing by 2048 and 1024 */
#define goodness_intronlen(x) (x >> 11)
#define goodness_nonintronlen(x) (x >> 10)

#define MAXITER 4
#define INTRON_PENALTY_INCONSISTENT 16
#define MININTRONLEN 6
#define CROSS_ONE_SHORT

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

/* Chimeras.  Not used currently */
#define DEBUG2 1
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
  int matches;
  int unknowns;
  int mismatches;
  int qopens;
  int qindels;
  int topens;
  int tindels;
  int goodness;
  int nexons;
  int coverage_correction;

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
Stage3_straintype (T this) {
  return this->straintype;
}

int
Stage3_goodness (T this) {
  return this->goodness;
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
  


double
Stage3_fracidentity (T this) {
  int den;

  if ((den = this->matches + this->mismatches + this->qindels + this->tindels) == 0) {
    return 1.0;
  } else {
    return (double) this->matches/(double) den;
  }
}

void
Stage3_pathscores (int *pathscores, T this, int querylength) {
  Pair_pathscores(pathscores,this->pairs,this->npairs,this->cdna_direction,querylength);
  return;
}


bool
Stage3_five_end (int *substart, int *subend, T this, Sequence_T queryseq) {
  int five_leftover, three_leftover;

  five_leftover = Pair_querypos(&(this->pairs[0]));
  three_leftover = Sequence_length(queryseq) - Pair_querypos(&(this->pairs[this->npairs-1]));
  if (five_leftover < three_leftover) {
    *substart = Pair_querypos(&(this->pairs[this->npairs-1]));
    *subend = Sequence_length(queryseq);
    return true;
  } else {
    *substart = 0;
    *subend = Pair_querypos(&(this->pairs[0]));
    return false;
  }
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
  Pair_T xstart, xend, ystart, yend;

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
Stage3_new (List_T pairs_fwd, List_T pairs_rev, 
	    Sequence_T genomicseg, Matchpairend_T matchpairend,
	    int straintype, char *strain, Chrnum_T chrnum, Genomicpos_T chrpos,
	    Genomicpos_T chroffset, int querylength, Genomicpos_T genomiclength, 
	    bool watsonp, double defect_rate, int ngap) {
  T new;
  int matches_fwd, matches_rev, mismatches_fwd, mismatches_rev, 
    unknowns_fwd, unknowns_rev, qopens_fwd, qindels_fwd, qopens_rev, qindels_rev, 
    topens_fwd, tindels_fwd, topens_rev, tindels_rev,
    ncanonical_fwd, ncanonical_rev, nsemicanonical_fwd, nsemicanonical_rev,
    nnoncanonical_fwd, nnoncanonical_rev, npairs_fwd, npairs_rev, coverage_correction;
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

    new->coverage_correction = 
      Sequence_count_bad(genomicseg,start->genomepos,start->querypos,-1) +
      Sequence_count_bad(genomicseg,end->genomepos,(querylength - 1) - end->querypos,+1);

    new->coverage = (double) (end->querypos - start->querypos + 1 + new->coverage_correction)/(double) querylength;

    return new;
  }
}
	    
void
Stage3_free (T *old) {
  List_T q;

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

void
Stage3_translate_genomic (T this, bool fulllengthp) {
  if (this->cdna_direction < 0) {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairs,this->npairs,/*backwardsp*/true,/*revcompp*/true,fulllengthp);
  } else {
    Translation_via_genomic(&this->translation_start,&this->translation_end,&this->translation_length,
			    &this->relaastart,&this->relaaend,
			    this->pairs,this->npairs,/*backwardsp*/false,/*revcompp*/false,fulllengthp);
  }
  return;
}

void
Stage3_translate_cdna (T this, T reference, bool literalrefp) {
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


void
Stage3_print_pathsummary (T this, int pathnum, IIT_T chromosome_iit, IIT_T contig_iit, 
			  char *dbversion, bool zerobasedp, int ntrimmed, bool fulllengthp) {
  Pair_T start, end;
  bool referencealignp;

  start = &(this->pairs[0]);
  end = &(this->pairs[this->npairs-1]);
  Stage3_translate_genomic(this,fulllengthp);
  referencealignp = this->straintype == 0 ? true : false;
  Pair_print_pathsummary(pathnum,start,end,this->chrnum,this->chrpos,this->chroffset,
			 chromosome_iit,referencealignp,this->strain,contig_iit,
			 dbversion,this->genomiclength,
			 this->nexons,this->coverage,this->coverage_correction,ntrimmed,
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
Stage3_print_mutations (T this, T reference, IIT_T chromosome_iit, char *dbversion, int ntrimmed, 
			bool showalignp, bool zerobasedp, 
			bool continuousp, bool diagnosticp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength) {
  Pair_T start, end;
  bool referencealignp;

  start = &(this->pairs[0]);
  end = &(this->pairs[this->npairs-1]);

  /*  Pair_dump_array(this->pairs,this->npairs,false); */

  referencealignp = this->straintype == 0 ? true : false;
  Pair_print_pathsummary(/*pathnum*/1,start,end,reference->chrnum,reference->chrpos,reference->chroffset,
			 chromosome_iit,referencealignp,this->strain,/*contig_iit*/NULL,
			 dbversion,reference->genomiclength,
			 this->nexons,this->coverage,this->coverage_correction,ntrimmed,
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
			 diagnosticp,/*genomicprop*/false,invertmode,nointronlenp,wraplength);
  }
  debug1(Pair_dump_array(this->pairs,this->npairs,/*zerobasedp*/true));
  debug1(Pair_check_array(this->pairs,this->npairs));

  return;
}



void
Stage3_print_map (T this, IIT_T map_iit, int pathnum, bool map_bothstrands_p) {
  Genomicpos_T position1, position2;
  Pair_T start, end;
  int chrpos1, chrpos2;
  Intlist_T iit_matches = NULL;
  int typeint;
  char *typestring;

  start = &(this->pairs[0]);
  end = &(this->pairs[this->npairs-1]);

  if (this->watsonp) {
    chrpos1 = this->chrpos + Pair_genomepos(start);
    chrpos2 = this->chrpos + Pair_genomepos(end);
    typestring = "FWD";
  } else {
    chrpos1 = this->chrpos + (this->genomiclength - 1) - Pair_genomepos(start);
    chrpos2 = this->chrpos + (this->genomiclength - 1) - Pair_genomepos(end);
    typestring = "REV";
  }

  position1 = this->chroffset + chrpos1;
  position2 = this->chroffset + chrpos2;

  if (map_bothstrands_p == true) {
    if (position1 < position2) {
      iit_matches = IIT_get(map_iit,position1,position2);
    } else {
      iit_matches = IIT_get(map_iit,position2,position1);
    }
    printf("  Map hits for path %d (%d):\n",pathnum,Intlist_length(iit_matches));
    IIT_print(map_iit,iit_matches,true);
  } else if ((typeint = IIT_typeint(map_iit,typestring)) < 0) {
    fprintf(stderr,"Warning: no type %s in %s.  Ignoring type.\n",typestring,IIT_name(map_iit));
    if (position1 < position2) {
      iit_matches = IIT_get(map_iit,position1,position2);
    } else {
      iit_matches = IIT_get(map_iit,position2,position1);
    }
    printf("  Map hits for path %d (%d):\n",pathnum,Intlist_length(iit_matches));
    IIT_print(map_iit,iit_matches,false);
  } else {
    if (position1 < position2) {
      iit_matches = IIT_get_typed(map_iit,position1,position2,typeint);
    } else {
      iit_matches = IIT_get_typed(map_iit,position2,position1,typeint);
    }
    printf("  Map hits for path %d (%d):\n",pathnum,Intlist_length(iit_matches));
    IIT_print(map_iit,iit_matches,false);
  }
  printf("\n");

  Intlist_free(&iit_matches);
  return;
}

void
Stage3_print_alignment (T this, IIT_T chromosome_iit,
			bool alignsummaryonlyp, bool universalp, bool zerobasedp,
			bool continuousp, bool diagnosticp, bool genomefirstp, int proteinmode,
			int invertmode, bool nointronlenp, int wraplength) {
  if (continuousp == true) {
    Pair_print_continuous(this->pairs,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			  this->genomiclength,this->watsonp,this->cdna_direction,universalp,zerobasedp,
			  diagnosticp,genomefirstp,invertmode,nointronlenp);
  } else {
    Pair_print_exonsummary(this->pairs,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			   chromosome_iit,this->genomiclength,this->watsonp,universalp,zerobasedp,
			   genomefirstp,invertmode);
    if (alignsummaryonlyp == false) {
      Pair_print_alignment(this->pairs,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			   chromosome_iit,this->genomiclength,this->watsonp,
			   this->cdna_direction,universalp,zerobasedp,
			   diagnosticp,/*genomicprop*/true,invertmode,nointronlenp,wraplength);
    }
  }
  debug1(Pair_dump_array(this->pairs,this->npairs,/*zerobasedp*/true));
  debug1(Pair_check_array(this->pairs,this->npairs));
  return;
}


void
Stage3_print_cdna_exons (T this, int wraplength) {
  Pair_print_cdna_exons(this->pairs,this->npairs,wraplength);
  return;
}


void
Stage3_print_protein_cdna (T this, int wraplength, bool fulllengthp) {
  Stage3_translate_genomic(this,fulllengthp);
  Pair_print_protein_cdna(this->pairs,this->npairs,wraplength);
  return;
}

void
Stage3_print_protein_genomic (T this, int wraplength, bool fulllengthp) {
  Stage3_translate_genomic(this,fulllengthp);
  Pair_print_protein_genomic(this->pairs,this->npairs,wraplength);
  return;
}


void
Stage3_print_compressed (T this, Sequence_T queryseq, IIT_T chromosome_iit,
			 char *version, int pathnum, int npaths,
			 bool checksump, bool chimerap, bool zerobasedp) {
  Pair_print_compressed(queryseq,version,pathnum,npaths,
			this->nexons,this->coverage,Stage3_fracidentity(this),
			this->pairs,this->npairs,this->chrnum,this->chrpos,this->chroffset,
			chromosome_iit,this->genomiclength,checksump,
			chimerap,this->strain,this->watsonp,zerobasedp);
  return;
}

static List_T 
add_null_gap (List_T pairs, int lastquerypos, int lastgenomepos, Pairpool_T pairpool,
	      int ngap) {
  int k, p;

  for (k = 0; k < ngap; k++) {
    pairs = Pairpool_push(pairs,pairpool,lastquerypos,lastgenomepos,' ','#',' ');
  }
  for (k = 0; k < 3; k++) {
    pairs = Pairpool_push(pairs,pairpool,lastquerypos,lastgenomepos,' ','.',' ');
  }
  for (k = 0; k < ngap; k++) {
    pairs = Pairpool_push(pairs,pairpool,lastquerypos,lastgenomepos,' ','#',' ');
  }

  return pairs;
}

static List_T
remove_null_gaps (List_T pairs, Pairpool_T pairpool) {
  List_T newpairs = NULL;
  Pair_T pair;

  while (pairs != NULL) {
    pair = pairs->first;
    if (pair->comp != '#' && pair->comp != '.') {
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
  int npeelback = 0, queryjump, genomejump, orig_querypos, orig_genomepos;
  bool exonp = true, skipp = true, gapp = false;
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

      if (pair->comp != '-' && pair->comp != ' ') {
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
  /* Last pair popped cannot be a gap */
  if (peeled != NULL) {
    firstpair = peeled->first;
    if (firstpair->comp == '-') {
      gapp = true;
    }
  }

  while (peeled != NULL && gapp) {
    peeled = Pairpool_pop(peeled,&pair);

    if (pair->comp != '-') {
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

  if (peeled == NULL) {
    *querypos = orig_querypos;
    *genomepos = orig_genomepos;
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
  int npeelback = 0, queryjump, genomejump, orig_lastquerypos, orig_lastgenomepos;
  bool exonp = true, skipp = true, gapp = false;
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

      if (pair->comp != '-' && pair->comp != ' ') {
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
  /* Last pair popped cannot be a gap */
  if (peeled != NULL) {
    firstpair = peeled->first;
    if (firstpair->comp == '-') {
      gapp = true;
    }
  }

  while (peeled != NULL && gapp) {
    peeled = Pairpool_pop(peeled,&pair);

    if (pair->comp != '-') {
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

  if (peeled == NULL) {
    *lastquerypos = orig_lastquerypos;
    *lastgenomepos = orig_lastgenomepos;
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
		     char *queryseq_ptr, char *genomicseg_ptr, Pairpool_T pairpool, Dynprog_T dynprog,
		     int maxpeelback, int extraband_single, double defect_rate, bool forcep) {
  List_T gappairs, peeled_pairs, peeled_path, p;
  Pair_T pair;
  int start, end, queryjump, genomejump, orig_querypos, orig_genomepos, j, i;
  int querydp5, querydp3, genomedp5, genomedp3;
  int nmatches, nmismatches, nopens, nindels;
  bool fillp = false;
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
    return NULL;
  } else {
    gappairs = Dynprog_single_gap(&finalscore,&nmatches,&nmismatches,&nopens,&nindels,dynprog,
				  &(queryseq_ptr[querydp5]),&(genomicseg_ptr[genomedp5]),
				  queryjump,genomejump,querydp5,genomedp5,pairpool,
				  extraband_single,defect_rate);
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
		   char *queryseq_ptr, char *genomicseg_ptr, int cdna_direction,
		   Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogR, int maxpeelback,
		   int extramaterial_paired, int extraband_paired, double defect_rate) {
  List_T gappairs, peeled_pairs, peeled_path, p;
  Pair_T pair;
  int querydp5, querydp3, genomedp5, genomedp3, queryjump, genomejump, j, i;
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
			      &(queryseq_ptr[querydp5]),&(queryseq_ptr[querydp3]),
			      &(genomicseg_ptr[genomedp5]),queryjump,queryjump,genomejump,
			      querydp5,querydp3,genomedp5,
			      pairpool,extraband_paired,defect_rate);
  debug(Pair_dump_list(gappairs,true));
  pairs = Pairpool_transfer(pairs,gappairs);
  return pairs;
}


static List_T
traverse_genome_gap (bool *filledp, int *nintrons, int *nnonintrons, int *intronlen, int *nonintronlen, 
		     List_T pairs, List_T *path,
		     int *querypos, int *genomepos, int lastquerypos, int lastgenomepos,
		     char *queryseq_ptr, char *genomicseg_ptr, int cdna_direction, int ngap,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogR, int maxpeelback,
		     int extramaterial_paired, int extraband_paired, double defect_rate, bool forcep) {
  List_T gappairs, peeled_pairs, peeled_path, micropairs, p;
  Pair_T pair;
  int querydp5, querydp3, genomedp5, genomedp3, queryjump, genomejump, j, i;
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
				&(queryseq_ptr[querydp5]),&(genomicseg_ptr[genomedp5]),
				&(genomicseg_ptr[genomedp3]),queryjump,genomejump,genomejump,
				querydp5,genomedp5,genomedp3,cdna_direction,ngap,
				pairpool,extraband_paired,false,
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
    micropairs = Dynprog_microexon_int(&microintrontype,&(queryseq_ptr[querydp5]),&(genomicseg_ptr[genomedp5]),
				       &(genomicseg_ptr[genomedp3]),queryjump,genomejump,genomejump,
				       querydp5,genomedp5,genomedp3,cdna_direction,
				       queryseq_ptr,genomicseg_ptr,pairpool,ngap);
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
			  int midquerypos, int midgenomepos, char *queryseq_ptr, 
			  char *genomicseg_ptr, int cdna_direction, int ngap,
			  Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, 
			  int maxpeelback, int nullgap, int extramaterial_paired, int extraband_paired,
			  double defect_rate) {
  List_T single_gappairs, dual_gappairs_1, dual_gappairs_2, peeled_pairs, peeled_path;
  int querydp5, querydp3, genomedp5, genomedp3, queryjump, genomejump, j, i;
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
				       &(queryseq_ptr[querydp5]),&(genomicseg_ptr[genomedp5]),
				       &(genomicseg_ptr[genomedp3]),queryjump,genomejump,genomejump,
				       querydp5,genomedp5,genomedp3,cdna_direction,ngap,
				       pairpool,extraband_paired,false,
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
				       &(queryseq_ptr[querydp5]),&(genomicseg_ptr[genomedp5]),
				       &(genomicseg_ptr[genomedp3]),queryjump,genomejump,genomejump,
				       querydp5,genomedp5,genomedp3,cdna_direction,ngap,
				       pairpool,extraband_paired,false,
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
				       &(queryseq_ptr[querydp5]),&(genomicseg_ptr[genomedp5]),
				       &(genomicseg_ptr[genomedp3]),queryjump,genomejump,genomejump,
				       querydp5,genomedp5,genomedp3,cdna_direction,ngap,
				       pairpool,extraband_paired,false,
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
		    char *queryseq_ptr, char *genomicseg_ptr, int cdna_direction, int introndir, Pairpool_T pairpool,
		    Dynprog_T dynprog, int maxpeelback, int extramaterial_end,
		    int extraband_end, double defect_rate, int ngap, bool extend_mismatch_p) {
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
  continuous_gappairs = Dynprog_end5_gap(&(*finalscore),&nmatches,&nmismatches,&nopens,&nindels,
					 dynprog,&(queryseq_ptr[querydp3]),&(genomicseg_ptr[genomedp3]),
					 queryjump,genomejump,querydp3,genomedp3,pairpool,
					 extraband_end,defect_rate,cdna_direction,ngap,extend_mismatch_p);
  continuous_goodness = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

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
					  &(queryseq_ptr[querydp3]),&(genomicseg_ptr[genomedp3]),
					  queryjump,genomejump,querydp3,genomedp3,cdna_direction,
					  queryseq_ptr,genomicseg_ptr,pairpool,ngap);
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
		    char *queryseq_ptr, char *genomicseg_ptr, int cdna_direction, int introndir,
		    Pairpool_T pairpool, Dynprog_T dynprog, int maxpeelback, int extramaterial_end,
		    int extraband_end, double defect_rate, int ngap, bool extend_mismatch_p) {
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
  debug(printf("Stage 3: Dynamic programming at end: querydp5 = %d, querydp3 = %d, genomedp5 = %d, genomedp3 = %d\n",
	       querydp5,querydp3,genomedp5,genomedp3));

  debug(printf("Trying to extend 3' end:\n"));
  continuous_gappairs = Dynprog_end3_gap(&(*finalscore),&nmatches,&nmismatches,&nopens,&nindels,
					 dynprog,&(queryseq_ptr[querydp5]),&(genomicseg_ptr[genomedp5]),
					 queryjump,genomejump,querydp5,genomedp5,pairpool,
					 extraband_end,defect_rate,cdna_direction,ngap,extend_mismatch_p);
  continuous_goodness = nmatches + MISMATCH*nmismatches + QOPEN*nopens + QINDEL*nindels;

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
					  &(queryseq_ptr[querydp5]),&(genomicseg_ptr[genomedp5]),
					  queryjump,genomejump,querydp5,genomedp5,cdna_direction,
					  queryseq_ptr,genomicseg_ptr,genomiclength,pairpool,ngap);
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
		  char *queryseq_ptr, char *genomicseg_ptr, int cdna_direction, int introndir,
		  int maxpeelback, int nullgap,
		  int extramaterial_end, int extramaterial_paired, int extraband_end,
		  double defect_rate, int ngap, Pairpool_T pairpool, Dynprog_T dynprogL, bool extend_mismatch_p) {
  List_T gappairs;
  int lastquerypos, lastgenomepos;
  int queryjump, genomejump, try_lastgenomepos, try;
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
				  queryseq_ptr,genomicseg_ptr,cdna_direction,introndir,
				  pairpool,dynprogL,maxpeelback,
				  extramaterial_end,extraband_end,defect_rate,ngap,extend_mismatch_p);
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
build_pairs_end5 (List_T pairs, char *queryseq_ptr, char *genomicseg_ptr, int cdna_direction, int introndir,
		  int maxpeelback, int nullgap,
		  int extramaterial_end, int extramaterial_paired, int extraband_end,
		  double defect_rate, int ngap, Pairpool_T pairpool, Dynprog_T dynprogR, bool extend_mismatch_p) {
  List_T best_gappairs, gappairs, endhits, p;
  Pair_T firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos;
  int queryjump, genomejump, try_genomepos, try;
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
				  queryseq_ptr,genomicseg_ptr,cdna_direction,introndir,
				  pairpool,dynprogR,maxpeelback,
				  extramaterial_end,extraband_end,defect_rate,ngap,extend_mismatch_p);
    debug(Pair_dump_list(gappairs,true));
    pairs = Pairpool_transfer(pairs,gappairs);

    /* Remove gaps at 5' end.  Shouldn't be needed, because dynamic
       programming should eliminate these. */
    /*
    if (pairs != NULL) {
      firstpair = pairs->first;
      pairs = Pairpool_pop(pairs,&firstpair);
      while ((firstpair->gapp || firstpair->comp == '-') && pairs != NULL) {
	pairs = Pairpool_pop(pairs,&firstpair);
      }
      if (firstpair->gapp == false && firstpair->comp != '-') {
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


static List_T
build_pairs_singles (List_T path, char *queryseq_ptr, char *genomicseg_ptr, int cdna_direction,
		     int maxpeelback, int nullgap,
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
      if (c1 == c2) {
	/* Possible for a stretch of the same nucleotide to have a mismatch */
	debug(printf("Stage 3: Single mismatch at %d %d: %c | %c\n",
		     lastquerypos-1,lastgenomepos-1,c1,c2));
	pairs = Pairpool_push(pairs,pairpool,lastquerypos-1,lastgenomepos-1,
			      c1,'|',c2);
      } else {
	debug(printf("Stage 3: Single mismatch at %d %d: %c   %c\n",
		     lastquerypos-1,lastgenomepos-1,c1,c2));
	pairs = Pairpool_push(pairs,pairpool,lastquerypos-1,lastgenomepos-1,
			      c1,' ',c2);
      }
      pair = pairs->first;
      lastquerypos = pair->querypos;
      lastgenomepos = pair->genomepos;

    } else if (queryjump == 1 && genomejump == 0) {
      /* Special case: Single cDNA insert */
      debug(printf("Stage 3: Single cDNA insert at %d %d: %c\n",
		   lastquerypos-1,lastgenomepos,queryseq_ptr[lastquerypos-1]));
      pairs = Pairpool_push(pairs,pairpool,lastquerypos-1,lastgenomepos,
			    queryseq_ptr[lastquerypos-1],'-',' ');
      pair = pairs->first;
      lastquerypos = pair->querypos;
      lastgenomepos = pair->genomepos;

    } else {
      /* Guarantees: queryjump <= nullgap && genomejump < queryjump - EXTRAQUERYGAP &&
	 genomejump <= queryjump + MININTRONLEN, meaning that score matrix is nearly square */
      debug(printf("Stage 3: Traversing single gap: querypos = %d, lastquerypos = %d, genomepos = %d, lastgenomepos = %d\n",
		   querypos,lastquerypos,genomepos,lastgenomepos));
      pairs = traverse_single_gap(&filledp,pairs,&path,&querypos,&genomepos,lastquerypos,lastgenomepos,
				  queryseq_ptr,genomicseg_ptr,pairpool,dynprogM,
				  maxpeelback,extraband_single,defect_rate,/*forcep*/false);
      if (filledp) {
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
			 char *queryseq_ptr, char *genomicseg_ptr, int cdna_direction,
			 int maxpeelback, int nullgap,
			 int extramaterial_paired, int extraband_paired, double defect_rate, int ngap,
			 Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T pairs, midexon_pairs;
  Pair_T pair, firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump, delta;
  int queryjump2, genomejump2;
  int midquerypos, midgenomepos, midgenomepos_start, midgenomepos_end, prevquerypos, prevgenomepos,
    midquerypos_start, midquerypos_end;
  bool single_one = false, exonp, shortp;

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
    if (queryjump <= 0 && genomejump <= 0) {
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
	  if (queryjump2 <= 0 && genomejump2 <= 0) {
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
					   midquerypos,midgenomepos,queryseq_ptr,genomicseg_ptr,cdna_direction,ngap,
					   pairpool,dynprogL,dynprogM,dynprogR,maxpeelback,nullgap,
					   extramaterial_paired,extraband_paired,defect_rate);
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
		     List_T pairs, int lastquerypos, int lastgenomepos, 
		     List_T path, char *queryseq_ptr, char *genomicseg_ptr, int cdna_direction,
		     int maxpeelback, int nullgap, int extramaterial_paired, 
		     int extraband_single, int extraband_paired, double defect_rate, int ngap,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  Pair_T pair;
  int querypos, genomepos, queryjump, genomejump;
  int midquerypos, midgenomepos, midgenomepos_start, midgenomepos_end, prevquerypos, prevgenomepos,
    midquerypos_start;
  bool filledp;

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
				queryseq_ptr,genomicseg_ptr,cdna_direction,
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
				  queryseq_ptr,genomicseg_ptr,cdna_direction,ngap,
				  pairpool,dynprogL,dynprogR,maxpeelback,extramaterial_paired,
				  extraband_paired,defect_rate,	/*forcep*/true);
      /* Do forcep, because adding large gap is not a good solution */

      if (filledp) {
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
				  queryseq_ptr,genomicseg_ptr,pairpool,dynprogM,
				  maxpeelback,extraband_single,defect_rate,
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
path_compute (int *intronlen, int *nonintronlen, 
	      List_T path, int cdna_direction, int introndir,
	      int querylength, Genomicpos_T genomiclength,
	      char *queryseq_ptr, char *genomicseg_ptr,  int maxpeelback, int nullgap,
	      int extramaterial_end, int extramaterial_paired,
	      int extraband_single, int extraband_end, int extraband_paired, double defect_rate, 
	      int ngap, Pairpool_T pairpool, 
	      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, bool extend_mismatch_p) {
  List_T pairs = NULL, smoothed_path;
  Pair_T firstpair;
  int querypos, genomepos, lastquerypos, lastgenomepos, nshortexons, iter;
  int newintrondir, nintrons = 0, nnonintrons = 0;
  bool singlep;

  if (path == NULL) {
    return NULL;
  }

  debug3(pairs = build_pairs_null_complete(path,pairpool,ngap));
  debug3(return pairs);

  /* Pass 1: solve single gaps */
  pairs = build_pairs_singles(path,queryseq_ptr,genomicseg_ptr,cdna_direction,
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
  debug5(path = List_reverse(pairs));
  debug5(pairs = build_pairs_null(path,pairpool,ngap));
  debug5(return pairs);

  if (nshortexons > 0) {
    path = List_reverse(pairs);
    pairs = build_pairs_dualintrons(&singlep,path,nshortexons,queryseq_ptr,genomicseg_ptr,cdna_direction,
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
    if (nshortexons > 0) {
      path = List_reverse(pairs);
      pairs = build_pairs_dualintrons(&singlep,path,nshortexons,queryseq_ptr,genomicseg_ptr,cdna_direction,
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
				queryseq_ptr,genomicseg_ptr,cdna_direction,
				maxpeelback,nullgap,extramaterial_paired,extraband_single,extraband_paired,
				defect_rate,ngap,pairpool,dynprogL,dynprogM,dynprogR);

    debug(printf("nintrons = %d, nnonintrons = %d\n",nintrons,nnonintrons));
    debug(printf("Intronlen = %d, Nonintronlen = %d\n",*intronlen,*nonintronlen));

    if (nintrons > nnonintrons) {
      newintrondir = cdna_direction;
    } else {
      newintrondir = introndir;
    }

    pairs = build_pairs_end5(pairs,queryseq_ptr,genomicseg_ptr,cdna_direction,
			     newintrondir,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_end,
			     defect_rate,ngap,pairpool,dynprogR,extend_mismatch_p);

    path = List_reverse(pairs);
    pairs = path;

    firstpair = pairs->first;
    querypos = firstpair->querypos;
    genomepos = firstpair->genomepos;

    pairs = build_pairs_end3(pairs,&querypos,&genomepos,querylength,genomiclength,
			     queryseq_ptr,genomicseg_ptr,cdna_direction,
			     newintrondir,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_end,
			     defect_rate,ngap,pairpool,dynprogL,extend_mismatch_p);

  }

  return List_reverse(pairs);
}


T
Stage3_copy (T old, Sequence_T queryseq, Genome_T genome, 
	     char *gbuffer1, char *gbuffer2, int gbufferlen, int ngap,
	     Pairpool_T pairpool) {
  T new;
  List_T pairs_fwd, pairs_rev;
  Sequence_T genomicseg;

  pairs_fwd = Pairpool_transfer_copy(NULL,old->pairs_fwd,pairpool);
  pairs_rev = Pairpool_transfer_copy(NULL,old->pairs_rev,pairpool);

  pairs_fwd = List_reverse(pairs_fwd);
  pairs_rev = List_reverse(pairs_rev);

  genomicseg = Genome_get_segment(genome,old->genomicstart,old->genomiclength,!old->watsonp,
				  gbuffer1,gbuffer2,gbufferlen);

  new = Stage3_new(pairs_fwd,pairs_rev,genomicseg,old->matchpairend,old->straintype,old->strain,
		   old->chrnum,old->chrpos,old->chroffset,
		   Sequence_length_full(queryseq),(Genomicpos_T) Sequence_length(genomicseg),
		   old->watsonp,old->defect_rate,ngap);

  Sequence_free(&genomicseg);

  return new;
}

T
Stage3_copy_bounded (T old, Sequence_T queryseq, Genome_T genome, 
		     char *gbuffer1, char *gbuffer2, int gbufferlen, int ngap,
		     int minpos, int maxpos, int maxpeelback, int nullgap, 
		     int extramaterial_end, int extramaterial_paired,
		     int extraband_single, int extraband_end, int extraband_paired,
		     Pairpool_T pairpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {

  T new;
  List_T pairs_fwd, pairs_rev;
  Sequence_T genomicseg;
  int introndir;

  if (old->pairs_fwd == NULL) {
    introndir = -1;
  } else if (old->pairs_rev == NULL) {
    introndir = +1;
  } else {
    introndir = 0;
  }

  pairs_fwd = Pairpool_transfer_copy_bounded(NULL,old->pairs_fwd,pairpool,minpos,maxpos);
  pairs_rev = Pairpool_transfer_copy_bounded(NULL,old->pairs_rev,pairpool,minpos,maxpos);

  /*
  pairs_fwd = build_pairs_end5(pairs_fwd,queryseq_ptr,genomicseg_ptr,+1,introndir,maxpeelback,nullgap,
			       extramaterial_end,extramaterial_paired,extraband_end,
			       old->defect_rate,ngap,pairpool,dynprogR);
  */

  pairs_fwd = List_reverse(pairs_fwd);
  pairs_rev = List_reverse(pairs_rev);

  genomicseg = Genome_get_segment(genome,old->genomicstart,old->genomiclength,!old->watsonp,
				  gbuffer1,gbuffer2,gbufferlen);

  new = Stage3_new(pairs_fwd,pairs_rev,genomicseg,old->matchpairend,old->straintype,old->strain,
		   old->chrnum,old->chrpos,old->chroffset,
		   Sequence_length_full(queryseq),(Genomicpos_T) Sequence_length(genomicseg),
		   old->watsonp,old->defect_rate,ngap);

  Sequence_free(&genomicseg);

  return new;
}


/* Introndir shows whether we are certain about the intron direction.
   However, we try both directions anyway, using cdna_direction */
T
Stage3_compute (Stage2_T stage2, Sequence_T queryseq, Sequence_T genomicseg, 
		Matchpairend_T matchpairend, int straintype, char *strain,
		int chrnum, Genomicpos_T chroffset, Genomicpos_T chrpos, 
		bool watsonp, int maxpeelback, int nullgap, 
		int extramaterial_end, int extramaterial_paired,
		int extraband_single, int extraband_end, int extraband_paired, int ngap,
		bool extend_mismatch_p, Pairpool_T pairpool, 
		Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		IIT_T altstrain_iit) {
  T this;
  List_T pairs_fwd, pairs_rev;
  double defect_rate;
  int querylength, querypos, genomepos;
  int nfwdintrons, nrevintrons, nunkintrons, introndir;
  int fwd_intronlen = 0, rev_intronlen = 0;
  int fwd_nonintronlen = 0, rev_nonintronlen = 0;
  bool fwdp, revp;
  Genomicpos_T genomiclength;
  char *queryseq_ptr, *genomicseg_ptr;
  Intlist_T typematches;

  debug(printf("Stage 3: *** Starting stage 3 for genomiclength %u at chrnum %d:%d\n",
	       Sequence_length(genomicseg),chrnum,chrpos);
	FREE(chrstring));

  defect_rate = Stage2_defect_rate(stage2);
  nfwdintrons = Stage2_nfwdintrons(stage2);
  nrevintrons = Stage2_nrevintrons(stage2);
  nunkintrons = Stage2_nunkintrons(stage2);

  querylength = Sequence_length_full(queryseq);
  genomiclength = (Genomicpos_T) Sequence_length(genomicseg);;
  queryseq_ptr = Sequence_pointer_full(queryseq);
  genomicseg_ptr = Sequence_pointer(genomicseg);

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
  debug(printf("Introns: %d fwd, %d rev, %d unk => fwdp = %d, revp = %d\n",
	       nfwdintrons, nrevintrons, nunkintrons, fwdp, revp));

  if (fwdp == true) {
    pairs_fwd = path_compute(&fwd_intronlen,&fwd_nonintronlen,Stage2_path(stage2),+1,introndir,
			     querylength,genomiclength,
			     queryseq_ptr,genomicseg_ptr,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_single,extraband_end,extraband_paired,
			     defect_rate,ngap,pairpool,dynprogL,dynprogM,dynprogR,extend_mismatch_p);
  } else {
    pairs_fwd = NULL;
  }
  if (revp == true) {
    pairs_rev = path_compute(&rev_intronlen,&rev_nonintronlen,Stage2_path(stage2),-1,introndir,
			     querylength,genomiclength,
			     queryseq_ptr,genomicseg_ptr,maxpeelback,nullgap,
			     extramaterial_end,extramaterial_paired,extraband_single,extraband_end,extraband_paired,
			     defect_rate,ngap,pairpool,dynprogL,dynprogM,dynprogR,extend_mismatch_p);
  } else {
    pairs_rev = NULL;
  }

  this = Stage3_new(pairs_fwd,pairs_rev,genomicseg,matchpairend,straintype,strain,
		    chrnum,chrpos,chroffset,querylength,genomiclength,watsonp,defect_rate,ngap);

  if (this == NULL) {
    return NULL;
  } else if (straintype == 0) {
    return this;
  } else {
    if (watsonp) {
      typematches = IIT_get_typed(altstrain_iit,this->genomicstart,this->genomicend,straintype);
    } else {
      typematches = IIT_get_typed(altstrain_iit,this->genomicend,this->genomicstart,straintype);
    }
    if (typematches == NULL) {
      Stage3_free(&this);
      return NULL;
    } else {
      Intlist_free(&typematches);
      return this;
    }
  }
}
