static char rcsid[] = "$Id: stage3hr.c 113087 2013-10-29 17:37:19Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage3hr.h"

#include <stdlib.h>		/* For qsort */
#include <string.h>
#include <strings.h>
#include <ctype.h>		/* For islower */
#include <math.h>		/* For exp() and log10() */
#include "assert.h"
#include "mem.h"
#include "chrnum.h"
/* #include "complement.h" */
#include "interval.h"
#include "listdef.h"
#include "substring.h"
#include "genome_hr.h"
#include "mapq.h"
#include "pair.h"		/* For Pair_print_gsnap and Pair_compute_mapq */
#include "maxent_hr.h"
#include "fastlog.h"


#define MAX_HITS 100000


#define CONCORDANT_TEXT "concordant"
#define PAIRED_TEXT "paired"
#define UNPAIRED_TEXT "unpaired"

#ifdef USE_TALLY_RATIO
#define TALLY_RATIO 2.0
#endif

#define SUBSUMPTION_SLOP 10	/* Should allow for short insert lengths */
/* #define TERMINAL_SECOND_CLASS 1 -- enabling this leads to poor results */

#define TERMINAL_COMPUTE_MINLENGTH 40


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Stage3end_new */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* printing */
#ifdef DEBUG1 
#define debug1(x) x
#else
#define debug1(x)
#endif


/* new_insertion */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* new_deletion */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Stage3end_optimal_score */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


/* Stage3_pair_up_concordant */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Stage3_pair_up_concordant, details */
#ifdef DEBUG5A
#define debug5a(x) x
#else
#define debug5a(x)
#endif

/* Stage3pair_optimal_score */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif


/* Stage3end_remove_duplicates and Stage3end_remove_overlaps */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Stage3pair_remove_overlaps */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Resolving ambiguous splice sites in Stage3pair_new */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* insert length calculation */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* genomicbound */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif

/* circular chromosomes */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif

/* Stage3pair_overlap */
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif

/* Stage3_determine_pairtype */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif



#define MAPQ_MAXIMUM_SCORE 40


static bool invert_first_p;
static bool invert_second_p;
static IIT_T genes_iit;
static int *genes_divint_crosstable;
static IIT_T tally_iit;
static int *tally_divint_crosstable;
static IIT_T runlength_iit;
static int *runlength_divint_crosstable;

static int terminal_output_minlength;
static int pairmax;

#ifdef USE_BINGO
static int expected_pairlength;
static int pairlength_deviation;
#endif

static int expected_pairlength_low;
static int expected_pairlength_high;
static int expected_pairlength_very_high;

static int localsplicing_penalty;
static int indel_penalty_middle;
static int antistranded_penalty;
static bool favor_multiexon_p;
static int gmap_min_nconsecutive;

static int ambig_end_interval;	/* For penalizing large ambiguous ends
				   in GMAP alignments, since such ends
				   should have been found */
static bool novelsplicingp;
static bool merge_samechr_p;
static bool *circularp;

/* Probably not good to use in certain genomic regions, unless we also
   use known splicesites with distance information. */
/* But sometimes need to use to get correct mapping */
static bool favor_ambiguous_p;


void
Stage3hr_setup (bool invert_first_p_in, bool invert_second_p_in,
		IIT_T genes_iit_in, int *genes_divint_crosstable_in,
		IIT_T tally_iit_in, int *tally_divint_crosstable_in,
		IIT_T runlength_iit_in, int *runlength_divint_crosstable_in,
		int terminal_output_minlength_in, bool distances_observed_p, int pairmax_in,
		Chrpos_T expected_pairlength, Chrpos_T pairlength_deviation,
		int localsplicing_penalty_in, int indel_penalty_middle_in,
		int antistranded_penalty_in, bool favor_multiexon_p_in,
		int gmap_min_nconsecutive_in, int index1part,
		int index1interval, bool novelsplicingp_in, bool merge_samechr_p_in,
		bool *circularp_in) {
  invert_first_p = invert_first_p_in;
  invert_second_p = invert_second_p_in;
  genes_iit = genes_iit_in;
  genes_divint_crosstable = genes_divint_crosstable_in;
  tally_iit = tally_iit_in;
  tally_divint_crosstable = tally_divint_crosstable_in;
  runlength_iit = runlength_iit_in;
  runlength_divint_crosstable = runlength_divint_crosstable_in;
  localsplicing_penalty = localsplicing_penalty_in;
  indel_penalty_middle = indel_penalty_middle_in;
  antistranded_penalty = antistranded_penalty_in;
  favor_multiexon_p = favor_multiexon_p_in;
  gmap_min_nconsecutive = gmap_min_nconsecutive_in;

  terminal_output_minlength = terminal_output_minlength_in;
  pairmax = pairmax_in;
  if (pairlength_deviation > expected_pairlength) {
    expected_pairlength_low = 0;
  } else {
    expected_pairlength_low = expected_pairlength - pairlength_deviation;
  }
  expected_pairlength_high = expected_pairlength + pairlength_deviation;
  expected_pairlength_very_high = expected_pairlength + 10*pairlength_deviation;

  if (distances_observed_p == true) {
    favor_ambiguous_p = false;
  } else {
    favor_ambiguous_p = true;
  }

#if 0
  ambig_end_interval = index1part + (index1interval - 1);
#else
  ambig_end_interval = 8;	/* Since GMAP uses 8-mers */
#endif

  novelsplicingp = novelsplicingp_in;
  merge_samechr_p = merge_samechr_p_in;
  circularp = circularp_in;

  return;
}



#ifdef DEBUG5
static char *
print_sense (int sense) {
  if (sense == SENSE_NULL) {
    return "sense:null";
  } else if (sense == SENSE_ANTI) {
    return "sense:anti";
  } else if (sense == SENSE_FORWARD) {
    return "sense:fwd";
  } else {
    abort();
  }
}
#endif



/* Note: Substring_T has genomiclength, but not Stage3end_T */
#define T Stage3end_T
struct T {
  Hittype_T hittype;
  int genestrand;
  bool sarrayp;			/* true if alignment found by suffix array */
  bool improved_by_gmap_p;	/* true if GMAP alignment based on this hit is better */

  Chrnum_T chrnum; /* Needed for printing paired-end results.  A chrnum of 0 indicates a distant splice. */
  Chrnum_T effective_chrnum;	/* For determining concordance */
  Chrnum_T other_chrnum;	/* 0 for non-translocations, and other chrnum besides effective_chrnum for translocations */
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  Chrpos_T chrlength;

  int querylength;		/* Needed for overlap calculations */
  int querylength_adj;		/* Adjusted for insertions */

  Univcoord_T genomicstart;
  Univcoord_T genomicend;
  bool plusp;

  Univcoord_T low;
  Univcoord_T high;
  Chrpos_T genomiclength;
  Chrpos_T guided_insertlength; /* Used only by Stage3end_eval_and_sort_guided */

  float mapq_loglik;
  int mapq_score;
  int absmq_score;		/* Absolute MAPQ, for XQ and X2 flags */

  int score;			/* Includes colordiffs and penalties */
  int ntscore;			/* Includes penalties */
  int nmatches;
  int nmatches_posttrim;
  int gmap_max_match_length;		/* Used only by GMAP */
  double gmap_min_splice_prob;		/* Used only by GMAP */

  int trim_left; /* Used by Stage3end_optimal_score for comparing terminals and non-terminals */
  int trim_right;
  bool trim_left_splicep;
  bool trim_right_splicep;

  int penalties;		/* Indel penalties */
  int score_eventrim;		/* Temporary storage used by Stage3end_optimal_score */

  Overlap_T gene_overlap;
  long int tally;

  int nmismatches_whole;
  int nmismatches_bothdiff;
  int nmismatches_refdiff;	/* Set only for display */

  int nindels;			/* for indels */
  int indel_pos;		/* for indels.  Relative to querypos 0 */
  int indel_low;		/* for indels.  Relative to chromosomal low end of read, but still 0 if no indel. */
  char *deletion;		/* for deletions */

  Chrpos_T distance;	/* for splicing or shortexon (sum of two distances) */
  Chrpos_T acceptor_distance; /* for shortexon */
  Chrpos_T donor_distance;	  /* for shortexon */

  int gmap_nindelbreaks;
  int gmap_cdna_direction;
  int gmap_nintrons;
  int sensedir;			/* for splicing */
  int sensedir_nonamb;			/* for splicing */

  bool start_ambiguous_p;
  bool end_ambiguous_p;
  int nchimera_known;
  int nchimera_novel;

  int amb_nmatches_start;	/* For splice, shortexon, and GMAP */
  int amb_nmatches_end;		/* For splice, shortexon, and GMAP */
  int amb_nmatches_donor;	/* For shortexon only */
  int amb_nmatches_acceptor;	/* For shortexon only */
  Endtype_T gmap_start_endtype;	/* For GMAP, which has no substrings */
  Endtype_T gmap_end_endtype;	/* For GMAP, which has no substrings */

  int *start_ambi;
  int *end_ambi;
  int start_nambi;
  int end_nambi;
  int *start_amb_nmismatches;
  int *end_amb_nmismatches;

  /* Single: substring1 */
  /* Indel: substring1 + substring2 */
  /* Halfsplice: substring1 */
  /* Splice: substring1 + substring2 */
  /* Shortexon: substring1 (shortexon) + substringD + substringA */

  /* Substrings should be in query order */
  Substring_T substring0;
  Substring_T substring1;	/* Main substring */
  Substring_T substring2;

  Substring_T substring_donor;	/* Just pointer to either substring1 or substring2 */
  Substring_T substring_acceptor; /* Just a pointer to either substring1 or substring2 */
  Substring_T substringD;	  /* Just a pointer to donor part of shortexon (substring0 or substring2) */
  Substring_T substringA;	  /* Just a pointer to acceptor part of shortexon (substring0 or substring2) */


  /* For GMAP alignment */
  struct Pair_T *pairarray;
  int npairs;
  int nsegments;

  Substring_T substring_low;	/* For SAM output */
  Substring_T substring_high;	/* For SAM output */

  bool paired_usedp;
  bool paired_seenp;   /* for paired-end.  set to true by Stage3_pair_up(). */
  bool concordantp;    /* for paired-end.  set to true by Stage3_pair_up(). */

  int alias;			/* -1 if below chrlength, 0 if straddles or NA (e.g., transloc), and +1 if above */
  int circularpos;		/* if alias == 0, then amount of queryseq below chrlength */
};


struct Stage3pair_T {
  Pairtype_T pairtype;
  int genestrand;

  T hit5;
  T hit3;
  bool private5p;			/* A private copy separate from hits5 and hits3, and not a pointer */
  bool private3p;

  Univcoord_T low;
  Univcoord_T high;
  int insertlength;
  int insertlength_expected_sign;	/* 1 if in (expected_pairlength_low, expected_pairlength_high),
					   0 if in (expected_pairlength_low, expected_pairlength_very_high), and
					   -1 if < expected_pairlength_low or > expected_pairlength_very_high */
				   
  Chrpos_T outerlength;

  float mapq_loglik;
  int mapq_score;
  int absmq_score;

  int score;
  int nmatches;
  int nmatches_posttrim;
  int indel_low; /* For ranking identical indel alignments, so we pick lowest coord */

  int score_eventrim;

  Overlap_T gene_overlap;
  long int tally;

#ifdef USE_ABSDIFFLENGTH
  Chrpos_T absdifflength;
#endif
#ifdef USE_BINGO
  bool absdifflength_bingo_p;
#endif
  int dir;			/* -1, 0, or +1 */
  bool sense_consistent_p;

  int nchimera_known;
  int nchimera_novel;

  bool circularp;		/* If either hit5 or hit3 are circular */
};



Hittype_T
Stage3end_hittype (T this) {
  return this->hittype;
}

static char *
hittype_string (Hittype_T hittype) {
  switch (hittype) {
  case EXACT: return "exact";
  case SUB: return "sub";
  case INSERTION: return "insertion";
  case DELETION: return "deletion";
  case HALFSPLICE_DONOR: return "donor";
  case HALFSPLICE_ACCEPTOR: return "acceptor";
  case SPLICE: return "splice";
  case SAMECHR_SPLICE: return "samechr_splice";
  case TRANSLOC_SPLICE: return "transloc_splice";
  case ONE_THIRD_SHORTEXON: return "one-third-shortexon";
  case TWO_THIRDS_SHORTEXON: return "two-thirds-shortexon";
  case SHORTEXON: return "shortexon";
  case GMAP: return "gmap";
  case TERMINAL: return "terminal";
  default: abort();
  }
}

char *
Stage3end_hittype_string (T this) {
  return hittype_string(this->hittype);
}

int
Stage3end_genestrand (T this) {
  return this->genestrand;
}

bool
Stage3end_sarrayp (T this) {
  return this->sarrayp;
}

bool
Stage3end_improved_by_gmap_p (T this) {
  return this->improved_by_gmap_p;
}

void
Stage3end_set_improved_by_gmap (T this) {
  this->improved_by_gmap_p = true;
  return;
}

bool
Stage3end_anomalous_splice_p (T this) {
  if (this->hittype == SAMECHR_SPLICE) {
    return true;
  } else {
    return false;
  }
}


Chrnum_T
Stage3end_chrnum (T this) {
  if (this == NULL) {
    return 0;
  } else {
    return this->chrnum;
  }
}

Chrnum_T
Stage3end_effective_chrnum (T this) {
  if (this == NULL) {
    return 0;
  } else {
    return this->effective_chrnum;
  }
}

Univcoord_T
Stage3end_chroffset (T this) {
  return this->chroffset;
}

Univcoord_T
Stage3end_chrhigh (T this) {
  return this->chrhigh;
}

Chrpos_T
Stage3end_chrlength (T this) {
  return this->chrlength;
}

Univcoord_T
Stage3end_genomicstart (T this) {
  return this->genomicstart;
}

Univcoord_T
Stage3end_genomicend (T this) {
  return this->genomicend;
}

/* For Goby */
int
Stage3end_query_alignment_length (T this) {
  int length;

  length = Substring_match_length(this->substring1);
  length += Substring_match_length(this->substring2);
  length += Substring_match_length(this->substring0);
  if (this->hittype == INSERTION) {
    length += this->nindels;
  }
  return length;
}

/* For Goby */
Chrpos_T
Stage3end_genomic_alignment_length (T this) {
  Chrpos_T length;

  length = Substring_genomic_alignment_length(this->substring1);
  length += Substring_genomic_alignment_length(this->substring2);
  length += Substring_genomic_alignment_length(this->substring0);
  if (this->hittype == DELETION) {
    length += (Chrpos_T) this->nindels;
  }
  return length;
}


Chrpos_T
Stage3end_chrpos_low_trim (T this) {
  if (this->plusp == true) {
    return Substring_alignstart_trim(this->substring_low) - Substring_chroffset(this->substring_low);
  } else {
    return Substring_alignend_trim(this->substring_low) - Substring_chroffset(this->substring_low);
  }
}

int
Stage3end_mapq_score (T this) {
  return this->mapq_score;
}

int
Stage3end_absmq_score (T this) {
  return this->absmq_score;
}

int
Stage3end_score (T this) {
  return this->score;
}

int
Stage3end_gmap_max_match_length (T this) {
  return this->gmap_max_match_length;
}

double
Stage3end_gmap_min_splice_prob (T this) {
  return this->gmap_min_splice_prob;
}


int
Stage3end_best_score (List_T hits) {
  List_T p;
  T hit;
  int best_score = 1000000;

  for (p = hits; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->score < best_score) {
      best_score = hit->score;
    }
  }

  return best_score;
}

bool
Stage3end_equiv_score_unpaired_p (List_T hits, int best_score) {
  List_T p;
  T hit;

  for (p = hits; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->paired_usedp == false && hit->score <= best_score) {
      return true;
    }
  }

  return false;
}

int
Stage3end_best_score_paired (List_T hits) {
  List_T p;
  T hit;
  int best_score = 1000000;

  for (p = hits; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->paired_usedp == true) {
      if (hit->score < best_score) {
	best_score = hit->score;
      }
    }
  }

  return best_score;
}


int
Stage3end_nmatches (T this) {
  return this->nmatches;
}

int
Stage3end_nmismatches_whole (T this) {
  return this->nmismatches_whole;
}

int
Stage3end_nmismatches_bothdiff (T this) {
  return this->nmismatches_bothdiff;
}

int
Stage3end_nmismatches_refdiff (T this) {
  return this->nmismatches_refdiff;
}

/* Called only for terminals */
Endtype_T
Stage3end_start_endtype (T this) {
  return Substring_start_endtype(this->substring1);
}

/* Called only for terminals */
Endtype_T
Stage3end_end_endtype (T this) {
  return Substring_end_endtype(this->substring1);
}

Endtype_T
Stage3end_gmap_start_endtype (T this) {
  return this->gmap_start_endtype;
}

Endtype_T
Stage3end_gmap_end_endtype (T this) {
  return this->gmap_end_endtype;
}

int
Stage3end_nindels (T this) {
  return this->nindels;
}

int
Stage3end_indel_pos (T this) {
  return this->indel_pos;
}


bool
Stage3end_plusp (T this) {
  return this->plusp;
}

bool
Stage3end_paired_usedp (T this) {
  return this->paired_usedp;
}

/* Accounts for ambig ends */
int
Stage3end_trim_left (T this) {
  return this->trim_left;
}

/* Accounts for ambig ends */
int
Stage3end_trim_right (T this) {
  return this->trim_right;
}

int
Stage3end_trim_left_raw (T this) {
  return this->trim_left + this->amb_nmatches_start;
}

int
Stage3end_trim_right_raw (T this) {
  return this->trim_right + this->amb_nmatches_end;
}

int
Stage3end_circularpos (T this) {
  return this->circularpos;
}


Substring_T
Stage3end_substring1 (T this) {
  return this->substring1;
}

Substring_T
Stage3end_substring2 (T this) {
  return this->substring2;
}

char *
Stage3end_deletion_string (T this) {
  return this->deletion;
}

Substring_T
Stage3end_substring_donor (T this) {
  assert(this->hittype == SPLICE || this->hittype == SAMECHR_SPLICE || this->hittype == TRANSLOC_SPLICE ||
	 this->hittype == HALFSPLICE_DONOR || this->hittype == HALFSPLICE_ACCEPTOR);
  return this->substring_donor;
}

Substring_T
Stage3end_substring_acceptor (T this) {
  assert(this->hittype == SPLICE || this->hittype == SAMECHR_SPLICE || this->hittype == TRANSLOC_SPLICE ||
	 this->hittype == HALFSPLICE_DONOR || this->hittype == HALFSPLICE_ACCEPTOR);
  return this->substring_acceptor;
}

Substring_T
Stage3end_substringD (T this) {
  assert(this->hittype == SHORTEXON || this->hittype == ONE_THIRD_SHORTEXON || this->hittype == TWO_THIRDS_SHORTEXON);
  return this->substringD;
}

Substring_T
Stage3end_substringA (T this) {
  assert(this->hittype == SHORTEXON || this->hittype == ONE_THIRD_SHORTEXON || this->hittype == TWO_THIRDS_SHORTEXON);
  return this->substringA;
}

Substring_T
Stage3end_substring_low (T this) {
  if (this == NULL) {
    return (Substring_T) NULL;
  } else {
    return this->substring_low;
  }
}

Substring_T
Stage3end_substring_high (T this) {
  if (this == NULL) {
    return (Substring_T) NULL;
  } else {
    return this->substring_high;
  }
}

Substring_T
Stage3end_substring_containing (T this, int querypos) {
  if (Substring_contains_p(this->substring1,querypos) == true) {
    return this->substring1;
  }
  if (this->substring2 != NULL && Substring_contains_p(this->substring2,querypos) == true) {
    return this->substring2;
  }
  if (this->substring0 != NULL && Substring_contains_p(this->substring0,querypos) == true) {
    return this->substring0;
  }
  return (Substring_T) NULL;
}



struct Pair_T *
Stage3end_pairarray (T this) {
  return this->pairarray;
}

int
Stage3end_npairs (T this) {
  return this->npairs;
}


Chrpos_T
Stage3end_distance (T this) {
  return this->distance;
}

Chrpos_T
Stage3end_shortexon_acceptor_distance (T this) {
  return this->acceptor_distance;
}

Chrpos_T
Stage3end_shortexon_donor_distance (T this) {
  return this->donor_distance;
}

int
Stage3end_sensedir (T this) {
  return this->sensedir;
}

int
Stage3end_sensedir_nonamb (T this) {
  return this->sensedir_nonamb;
}

int
Stage3end_cdna_direction (T this) {
  if (this == NULL) {
    return SENSE_NULL;
#if 0
  } else if (this->pairarray != NULL) {
    return this->gmap_cdna_direction;
#endif
  } else if (this->sensedir == SENSE_FORWARD) {
    return +1;
  } else if (this->sensedir == SENSE_ANTI) {
    return -1;
  } else {
#if 0
    /* Leads to non-canonical XS:A:? output in SAM format */
    return 0;
#else
    return this->gmap_cdna_direction;
#endif
  }
}

int
Stage3end_nintrons (T this) {
  return this->gmap_nintrons;
}

int
Stage3end_gmap_querystart (T this) {
  return this->pairarray[0].querypos;
}

int
Stage3end_gmap_queryend (T this) {
  return this->pairarray[this->npairs - 1].querypos;
}


int
Stage3end_terminal_trim (T this) {
  if (this->hittype != TERMINAL) {
    return 0;
  } else {
    return Substring_trim_left(this->substring1) + Substring_trim_right(this->substring1);
  }
}

int
Stage3end_terminal_length (T this) {
  assert(this->hittype == TERMINAL);
  return Substring_queryend(this->substring1) - Substring_querystart(this->substring1) + 1;
}



static Overlap_T
Stage3end_gene_overlap (T this) {
  Overlap_T overlap;
  bool foundp = false;

  if (this->hittype == GMAP) {
    return Pair_gene_overlap(this->pairarray,this->npairs,genes_iit,
			     genes_divint_crosstable[this->chrnum],favor_multiexon_p);
  } else {
    if ((overlap = Substring_gene_overlap(this->substring1,favor_multiexon_p)) == KNOWN_GENE_MULTIEXON) {
      return KNOWN_GENE_MULTIEXON;
    } else if (overlap == KNOWN_GENE) {
      if (favor_multiexon_p == false) {
	return KNOWN_GENE;
      } else {
	foundp = true;
      }
    }

    if (this->substring2 != NULL) {
      if ((overlap = Substring_gene_overlap(this->substring2,favor_multiexon_p)) == KNOWN_GENE_MULTIEXON) {
	return KNOWN_GENE_MULTIEXON;
      } else if (overlap == KNOWN_GENE) {
	if (favor_multiexon_p == false) {
	  return KNOWN_GENE;
	} else {
	  foundp = true;
	}
      }
    }

    if (this->substring0 != NULL) {
      if ((overlap = Substring_gene_overlap(this->substring0,favor_multiexon_p)) == KNOWN_GENE_MULTIEXON) {
	return KNOWN_GENE_MULTIEXON;
      } else if (overlap == KNOWN_GENE) {
	if (favor_multiexon_p == false) {
	  return KNOWN_GENE;
	} else {
	  foundp = true;
	}
      }
    }

    if (foundp == true) {
      return KNOWN_GENE;
    } else {
      return NO_KNOWN_GENE;
    }
  }
}


bool
Stage3end_contains_known_splicesite (T this) {
  /* assert(this->hittype != GMAP); */

  /* indel + splice => requires gmap
     doublesplice + splice => requires gmap
     other cases should have already been covered by gsnap
  */

  if (this->hittype == GMAP) {
    /* Possible now because performing redo of GMAP for sense inconsistency */
    return false;
  } else if (this->hittype != INSERTION && this->hittype != DELETION && this->hittype != SHORTEXON) {
    return false;
  } else if (Substring_contains_known_splicesite(this->substring1) == true) {
    return true;
  } else if (this->substring2 != NULL && Substring_contains_known_splicesite(this->substring2) == true) {
    return true;
  } else if (this->substring0 != NULL && Substring_contains_known_splicesite(this->substring0) == true) {
    return true;
  } else {
    return false;
  }
}

bool
Stage3end_indel_contains_known_splicesite (bool *leftp, bool *rightp, T this) {
  /* indel + splice => requires gmap */

  *leftp = Substring_contains_known_splicesite(this->substring1);
  *rightp = Substring_contains_known_splicesite(this->substring2);
  if (*leftp == true || *rightp == true) {
    return true;
  } else {
    return false;
  }
}


#if 0
bool
Stage3end_bad_stretch_p (T this, Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  if (this->hittype == GMAP) {
#if 0
    if (this->gmap_cdna_direction != 0 && this->sensedir == SENSE_NULL) {
      /* Doesn't work for alignments without introns */
      debug0(printf("Bad GMAP: cdna_direction %d and sense null\n",this->gmap_cdna_direction));
      return true;
    }
#endif

    if (this->gmap_nindelbreaks > 3) {
      debug0(printf("Bad GMAP: nindel breaks %d > 3\n",this->gmap_nindelbreaks));
      return true;
#if 0
    } else if (this->gmap_min_splice_prob < 0.5) {
      /* Calculation is buggy */
      debug0(printf("Bad GMAP: min splice prob %f < 0.5\n",this->gmap_min_splice_prob));
      return true;
    } else {
      return Stage3_bad_stretch_p(this->pairarray,this->npairs,/*pos5*/this->trim_left,
				  /*pos3*/this->querylength_adj - this->trim_right);
      ngoodpart = Stage3_good_part(this->pairarray,this->npairs,/*pos5*/this->trim_left,
				   /*pos3*/this->querylength_adj - this->trim_right);
      if (ngoodpart < this->querylength_adj/2) {
	return true;
      } else {
	return false;
      }
#endif
    }
  } else if (Substring_bad_stretch_p(this->substring1,query_compress_fwd,query_compress_rev) == true) {
    return true;
  } else if (this->substring2 != NULL && Substring_bad_stretch_p(this->substring2,query_compress_fwd,query_compress_rev) == true) {
    return true;
  } else if (this->substring0 != NULL && Substring_bad_stretch_p(this->substring0,query_compress_fwd,query_compress_rev) == true) {
    return true;
  } else {
    return false;
  }
}
#endif



static long int
Stage3end_compute_tally (T this) {
  long int tally = 0L;

  tally = 0L;
  if (this->substring1 != NULL) {
    tally += Substring_tally(this->substring1,tally_iit,tally_divint_crosstable);
  }
  if (this->substring2 != NULL) {
    tally += Substring_tally(this->substring2,tally_iit,tally_divint_crosstable);
  }
  if (this->substring0 != NULL) {
    tally += Substring_tally(this->substring0,tally_iit,tally_divint_crosstable);
  }

  return tally;
}

static bool
Stage3end_runlength_p (T this) {
  if (this->substring1 != NULL && Substring_runlength_p(this->substring1,runlength_iit,runlength_divint_crosstable) == true) {
    return true;
  }
  if (this->substring2 != NULL && Substring_runlength_p(this->substring2,runlength_iit,runlength_divint_crosstable) == true) {
    return true;
  }
  if (this->substring0 != NULL && Substring_runlength_p(this->substring0,runlength_iit,runlength_divint_crosstable) == true) {
    return true;
  }

  return false;
}


#if 0
/* Tries to use tally information.  Now obsolete */
static long int
Stage3end_tally (T this) {

  if (tally_iit == NULL) {
    return 0L;
  } else if (this->tally >= 0) {
    return this->tally;
  } else {
    this->tally = 0L;
    if (this->substring1 != NULL) {
      this->tally += Substring_tally(this->substring1,tally_iit,tally_divint_crosstable);
    }
    if (this->substring2 != NULL) {
      this->tally += Substring_tally(this->substring2,tally_iit,tally_divint_crosstable);
    }
    if (this->substring0 != NULL) {
      this->tally += Substring_tally(this->substring0,tally_iit,tally_divint_crosstable);
    }

    return this->tally;
  }
}
#endif


bool
Stage3end_genomicbound_from_start (Univcoord_T *genomicbound, T this, int overlap, Univcoord_T chroffset) {
  int substring_length;

  debug11(printf("Stage3end_genomicbound_from_start with overlap %d\n",overlap));
  if (this->hittype == GMAP) {
    debug11(printf("  Computing on GMAP\n"));
    *genomicbound = chroffset + Pairarray_genomicbound_from_start(this->pairarray,this->npairs,overlap);
    return true;
  } else {
    debug11(printf("  Computing on substrings\n"));
    if (this->substring0 != NULL) {
      debug11(printf("  Substring 0 has length %d\n",Substring_match_length_orig(this->substring0)));
      if ((substring_length = Substring_match_length_orig(this->substring0)) >= overlap) {
	if (this->plusp == true) {
	  *genomicbound = Substring_alignstart(this->substring0) + overlap;
	} else {
	  *genomicbound = Substring_alignstart(this->substring0) - overlap;
	}
	return true;
      } else {
	overlap -= substring_length;
      }
    }

    if ((substring_length = Substring_match_length_orig(this->substring1)) >= overlap) {
      debug11(printf("  Substring 1 has length %d\n",Substring_match_length_orig(this->substring1)));
      if (this->plusp == true) {
	*genomicbound = Substring_alignstart(this->substring1) + overlap;
      } else {
	*genomicbound = Substring_alignstart(this->substring1) - overlap;
      }
      return true;
    } else {
      overlap -= substring_length;
    }

    if (this->substring2 != NULL) {
      debug11(printf("  Substring 2 has length %d\n",Substring_match_length_orig(this->substring2)));
      if ((substring_length = Substring_match_length_orig(this->substring2)) >= overlap) {
	if (this->plusp == true) {
	  *genomicbound = Substring_alignstart(this->substring2) + overlap;
	} else {
	  *genomicbound = Substring_alignstart(this->substring2) - overlap;
	}
	return true;
      } else {
	overlap -= substring_length;
      }
    }

    debug11(printf("Still have %d of overlap\n",overlap));
    return false;
  }
}

bool
Stage3end_genomicbound_from_end (Univcoord_T *genomicbound, T this, int overlap, Univcoord_T chroffset) {
  int substring_length;

  debug11(printf("Stage3end_genomicbound_from_end with overlap %d\n",overlap));
  if (this->hittype == GMAP) {
    debug11(printf("  Computing on GMAP\n"));
    *genomicbound = chroffset + Pairarray_genomicbound_from_end(this->pairarray,this->npairs,overlap);
    return true;
  } else {
    debug11(printf("  Computing on substrings\n"));
    if (this->substring2 != NULL) {
      debug11(printf("  Substring 2 has length %d\n",Substring_match_length_orig(this->substring2)));
      if ((substring_length = Substring_match_length_orig(this->substring2)) >= overlap) {
	if (this->plusp == true) {
	  *genomicbound = Substring_alignend(this->substring2) - overlap;
	} else {
	  *genomicbound = Substring_alignend(this->substring2) + overlap;
	}
	return true;
      } else {
	overlap -= substring_length;
      }
    }

    if ((substring_length = Substring_match_length_orig(this->substring1)) >= overlap) {
      debug11(printf("  Substring 1 has length %d\n",Substring_match_length_orig(this->substring1)));
      if (this->plusp == true) {
	*genomicbound = Substring_alignend(this->substring1) - overlap;
      } else {
	*genomicbound = Substring_alignend(this->substring1) + overlap;
      }
      return true;
    } else {
      overlap -= substring_length;
    }

    if (this->substring0 != NULL) {
      debug11(printf("  Substring 0 has length %d\n",Substring_match_length_orig(this->substring0)));
      if ((substring_length = Substring_match_length_orig(this->substring0)) >= overlap) {
	if (this->plusp == true) {
	  *genomicbound = Substring_alignend(this->substring0) - overlap;
	} else {
	  *genomicbound = Substring_alignend(this->substring0) + overlap;
	}
	return true;
      } else {
	overlap -= substring_length;
      }
    }

    debug11(printf("Still have %d of overlap\n",overlap));
    return false;
  }
}


void
Stage3end_free (T *old) {
  debug0(printf("Freeing Stage3end %p of type %s\n",*old,hittype_string((*old)->hittype)));

  FREE_OUT((*old)->end_ambi);
  FREE_OUT((*old)->start_ambi);
  FREE_OUT((*old)->end_amb_nmismatches);
  FREE_OUT((*old)->start_amb_nmismatches);

  if ((*old)->deletion != NULL) {
    FREE_OUT((*old)->deletion);
  }

  if ((*old)->pairarray != NULL) {
    FREE_OUT((*old)->pairarray);
  }

  if ((*old)->substring1 != NULL) {
    Substring_free(&(*old)->substring1);
  }
  if ((*old)->substring2 != NULL) {
    Substring_free(&(*old)->substring2);
  }
  if ((*old)->substring0 != NULL) {
    Substring_free(&(*old)->substring0);
  }

  
  FREE_OUT(*old);
  return;
}


/* Used for freeing gmap_history_values in stage1hr.c */
void
Stage3end_list_free (List_T *values) {
  List_T p;
  T hit;

  for (p = *values; p != NULL; p = p->rest) {
    if ((hit = (T) p->first) != NULL) {
      Stage3end_free(&hit);
    }
  }
  List_free(&(*values));
  return;
}



bool
Stage3pair_anomalous_splice_p (Stage3pair_T this) {
  if (this->hit5 != NULL && this->hit5->hittype == SAMECHR_SPLICE) {
    return true;
  } else if (this->hit3 != NULL && this->hit3->hittype == SAMECHR_SPLICE) {
    return true;
  } else {
    return false;
  }
}


int
Stage3pair_genestrand (Stage3pair_T this) {
  return this->genestrand;
}

Stage3end_T
Stage3pair_hit5 (Stage3pair_T this) {
  return this->hit5;
}

Stage3end_T
Stage3pair_hit3 (Stage3pair_T this) {
  return this->hit3;
}

int
Stage3pair_mapq_score (Stage3pair_T this) {
  return this->mapq_score;
}

int
Stage3pair_absmq_score (Stage3pair_T this) {
  return this->absmq_score;
}

Chrpos_T
Stage3pair_pairlength (Stage3pair_T this) {
  return this->insertlength;
}

int
Stage3pair_nmatches (Stage3pair_T this) {
  return this->nmatches;
}


bool
Stage3pair_concordantp (List_T hitpairs) {
  List_T p;
  Stage3pair_T hitpair;

  for (p = hitpairs; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    if (hitpair->pairtype == CONCORDANT) {
      return true;
    }
  }
  return false;
}

List_T
Stage3pair_filter_nonconcordant (List_T hitpairs) {
  List_T filtered = NULL, p;
  Stage3pair_T hitpair;

  for (p = hitpairs; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    if (hitpair->pairtype != CONCORDANT) {
      Stage3pair_free(&hitpair);
    } else {
      filtered = List_push(filtered,(void *) hitpair);
    }
  }
  List_free(&hitpairs);
  return filtered;
}


static int
gmap5_substring3_overlap (Stage3end_T hit5, Stage3end_T hit3, Substring_T substring) {
  Univcoord_T chroffset;
  int i;

  i = 0;
  chroffset = hit3->chroffset;
  if (1 || hit5->plusp == true) {
    debug13(printf("plus goal: %u or %u\n",Substring_alignstart_trim(substring) - chroffset,Substring_alignend_trim(substring) - chroffset));
    while (i < hit5->npairs) {
      if (hit5->pairarray[i].genomepos == Substring_alignstart_trim(substring) - chroffset) {
	debug13(printf("case 1\n"));
	return (hit5->pairarray[i].querypos - hit5->pairarray[0].querypos + 1) +
	  hit3->querylength - Substring_querystart(substring) - hit3->trim_right - hit3->amb_nmatches_end;
      } else if (hit5->pairarray[i].genomepos == Substring_alignend_trim(substring) - chroffset) {
	debug13(printf("case 2\n"));
	return (hit5->pairarray[i].querypos - hit5->pairarray[0].querypos + 1) +
	  hit3->querylength - Substring_queryend(substring) - hit3->trim_right - hit3->amb_nmatches_end;
      } else {
	i++;
      }
    }
    return 0;

  } else {
    debug13(printf("minus goal: %u or %u\n",Substring_alignstart_trim(substring) - chroffset,Substring_alignend_trim(substring) - chroffset));
    while (i < hit5->npairs) {
      if (hit5->pairarray[i].genomepos == Substring_alignstart_trim(substring) - chroffset) {
	debug13(printf("case 3: genomepos %u, at GMAP %d.  substring trim %d, **%d..%d trim %d\n",
		       hit5->pairarray[i].genomepos,hit5->pairarray[i].querypos,
		       hit3->trim_left,Substring_querystart(substring),Substring_queryend(substring),hit3->trim_right));
	return (hit5->pairarray[i].querypos - hit5->pairarray[0].querypos + 1) +
	  (hit3->querylength - Substring_querystart(substring) - hit3->trim_right - hit3->amb_nmatches_end); /* trim2_1&trim2_2 */
      } else if (hit5->pairarray[i].genomepos == Substring_alignend_trim(substring) - chroffset) {
	debug13(printf("case 4: genomepos %u, at GMAP %d.  substring trim %d, %d..%d** trim %d\n",
		       hit5->pairarray[i].genomepos,hit5->pairarray[i].querypos,
		       hit3->trim_left,Substring_querystart(substring),Substring_queryend(substring),hit3->trim_right));
	return (hit5->pairarray[i].querypos - hit5->pairarray[0].querypos + 1) +
	  hit3->querylength - Substring_queryend(substring) - hit3->trim_right - hit3->amb_nmatches_end;
      } else {
	i++;
      }
    }
    return 0;
  }
}

static int
substring5_gmap3_overlap (Stage3end_T hit5, Stage3end_T hit3, Substring_T substring) {
  Univcoord_T chroffset;
  int j;

  j = 0;
  chroffset = hit5->chroffset;
  if (1 || hit5->plusp == true) {
    debug13(printf("plus goal: %u or %u\n",Substring_alignstart_trim(substring) - chroffset,Substring_alignend_trim(substring) - chroffset));
    while (j < hit3->npairs) {
      if (hit3->pairarray[j].genomepos == Substring_alignstart_trim(substring) - chroffset) {
	debug13(printf("case 4: genomepos %u, at GMAP %d.  substring trim %d, **%d..%d trim %d\n",
		       hit3->pairarray[j].genomepos,hit3->pairarray[j].querypos,
		       hit5->trim_left,Substring_querystart(substring),Substring_queryend(substring),hit5->trim_right));
	return (hit3->pairarray[hit3->npairs-1].querypos - hit3->pairarray[j].querypos + 1) + 
	  Substring_querystart(substring) - hit5->trim_left - hit5->amb_nmatches_start; /* okay */
      } else if (hit3->pairarray[j].genomepos == Substring_alignend_trim(substring) - chroffset) {
	debug13(printf("case 6: genomepos %u, at GMAP %d.  substring trim %d, %d..%d** trim %d\n",
		       hit3->pairarray[j].genomepos,hit3->pairarray[j].querypos,
		       hit5->trim_left,Substring_querystart(substring),Substring_queryend(substring),hit5->trim_right));
	return (hit3->pairarray[hit3->npairs-1].querypos - hit3->pairarray[j].querypos + 1) +
	  Substring_queryend(substring) - hit5->trim_left - hit5->amb_nmatches_start; /* okay: trim2_2&trim2_1 */
      } else {
	j++;
      }
    }
    return 0;
    
  } else {
    debug13(printf("minus goal: %u or %u\n",Substring_alignstart_trim(substring) - chroffset,Substring_alignend_trim(substring) - chroffset));
    while (j < hit3->npairs) {
      if (hit3->pairarray[j].genomepos == Substring_alignstart_trim(substring) - chroffset) {
	debug13(printf("case 7\n"));
	return (hit3->pairarray[hit3->npairs-1].querypos - hit3->pairarray[j].querypos + 1) +
	  Substring_querystart(substring) - hit5->trim_left - hit5->amb_nmatches_start;
      } else if (hit3->pairarray[j].genomepos == Substring_alignend_trim(substring) - chroffset) {
	debug13(printf("case 8\n"));
	return (hit3->pairarray[hit3->npairs-1].querypos - hit3->pairarray[j].querypos + 1) +
	  Substring_queryend(substring) - hit5->trim_left - hit5->amb_nmatches_start; /* okay: trim2_2&trim2_1 */
      } else {
	j++;
      }
    }
    return 0;
  }
}


/* Needed to compute overlap properly.  Based on pair_insert_length below, plus code for handling GMAP. */
static Chrpos_T
pair_insert_length_trimmed (Stage3end_T hit5, Stage3end_T hit3, int hit5_trimmed_length, int hit3_trimmed_length) {
  int i, j;
  int overlap;
  Univcoord_T start5, end5, start3, end3;

  assert(hit5->plusp == hit3->plusp);

  if (hit5->hittype == GMAP && hit3->hittype == GMAP) {
    i = j = 0;
    if (hit5->plusp == true) {
      while (i < hit5->npairs && j < hit3->npairs) {
	if (hit5->pairarray[i].genomepos < hit3->pairarray[j].genomepos) {
	  i++;
	} else if (hit5->pairarray[i].genomepos > hit3->pairarray[j].genomepos) {
	  j++;
	} else {
	  debug13(printf("GMAP and GMAP show overlap at position %d, querypos %d and %d\n",
			 hit5->pairarray[i].genomepos,hit5->pairarray[i].querypos,hit3->pairarray[j].querypos));
	  return (hit5->pairarray[i].querypos - hit5->pairarray[0].querypos) + (hit3->pairarray[hit3->npairs-1].querypos - hit3->pairarray[j].querypos);
	}
      }
      debug13(printf("GMAP and GMAP show no overlap\n"));
      return 0;

    } else {
      while (i < hit5->npairs && j < hit3->npairs) {
	if (hit5->pairarray[i].genomepos > hit3->pairarray[j].genomepos) {
	  i++;
	} else if (hit5->pairarray[i].genomepos < hit3->pairarray[j].genomepos) {
	  j++;
	} else {
	  debug13(printf("GMAP and GMAP show overlap at position %d, querypos %d and %d\n",
			 hit5->pairarray[i].genomepos,hit5->pairarray[i].querypos,hit3->pairarray[j].querypos));
	  return (hit5->pairarray[i].querypos - hit5->pairarray[0].querypos) + (hit3->pairarray[hit3->npairs-1].querypos - hit3->pairarray[j].querypos);
	}
      }
      debug13(printf("GMAP and GMAP show no overlap\n"));
      return 0;
    }
    
  } else if (hit5->hittype == GMAP) {
    if ((overlap = gmap5_substring3_overlap(hit5,hit3,hit3->substring1)) > 0) {
      return overlap;
    } else if (hit3->substring2 != NULL && (overlap = gmap5_substring3_overlap(hit5,hit3,hit3->substring2)) > 0) {
      return overlap;
    } else if (hit3->substring0 != NULL && (overlap = gmap5_substring3_overlap(hit5,hit3,hit3->substring0)) > 0) {
      return overlap;
    } else {
      return 0;
    }

  } else if (hit3->hittype == GMAP) {
    if ((overlap = substring5_gmap3_overlap(hit5,hit3,hit5->substring1)) > 0) {
      return overlap;
    } else if (hit5->substring2 != NULL && (overlap = substring5_gmap3_overlap(hit5,hit3,hit5->substring2)) > 0) {
      return overlap;
    } else if (hit5->substring0 != NULL && (overlap = substring5_gmap3_overlap(hit5,hit3,hit5->substring0)) > 0) {
      return overlap;
    } else {
      return 0;
    }

  } else if (hit5->plusp == true) {
    start5 = hit5->genomicstart + hit5->trim_left + hit5->amb_nmatches_start;
    end5 = hit5->genomicend - hit5->trim_right - hit5->amb_nmatches_end;
    start3 = hit3->genomicstart + hit3->trim_left + hit3->amb_nmatches_start;
    end3 = hit3->genomicend - hit3->trim_right - hit3->amb_nmatches_end;

    if (end3 < start5) {
      /* Case 1 */
      return 0;
    } else if (end5 < start3) {
      /* Case 6 */
      return 0;
    } else if (start3 < start5) {
      if (end3 < end5) {
	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug13(printf("plus case 2a: start5 %u\n",start5 - hit5->chroffset));
	if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return (Substring_alignend_trim(hit3->substring0) - start5)
	    + hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return (Substring_alignend_trim(hit3->substring1) - start5)
	    + hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return (Substring_alignend_trim(hit3->substring2) - start5)
	    /* + hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end */;
	}

	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug13(printf("plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start +
	    (end3 - Substring_alignstart_trim(hit5->substring2));
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start +
	    (end3 - Substring_alignstart_trim(hit5->substring1));
	} else if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return /* Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start + */
	    (end3 - Substring_alignstart_trim(hit5->substring0));
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug13(printf("plus case 3\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return hit5_trimmed_length + (Substring_alignend_trim(hit3->substring2) - end5)
	    /* + hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end */;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return hit5_trimmed_length + (Substring_alignend_trim(hit3->substring1) - end5)
	    + hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return hit5_trimmed_length + (Substring_alignend_trim(hit3->substring0) - end5)
	    + hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end;
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (end3 < end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug13(printf("plus case 4\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return /* (Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start) + */
	    (start3 - Substring_alignstart_trim(hit5->substring0)) + hit3_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return (Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start) +
	    (start3 - Substring_alignstart_trim(hit5->substring1)) + hit3_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return (Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start) +
	    (start3 - Substring_alignstart_trim(hit5->substring2)) + hit3_trimmed_length;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug13(printf("plus case 5a\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return /* (Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start) + */
	    (start3 - Substring_alignstart_trim(hit5->substring0)) + hit3_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return (Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start) +
	    (start3 - Substring_alignstart_trim(hit5->substring1)) + hit3_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return (Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start) +
	    (start3 - Substring_alignstart_trim(hit5->substring2)) + hit3_trimmed_length;
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug13(printf("plus case 5b\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return /* hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end + */
	    (Substring_alignend_trim(hit3->substring2) - end5) + hit5_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end +
	    (Substring_alignend_trim(hit3->substring1) - end5) + hit5_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end +
	    (Substring_alignend_trim(hit3->substring0) - end5) + hit5_trimmed_length;
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug13(printf("plus general: hit3->substring1\n"));
    if (Substring_overlap_segment_trimmed_p(hit5->substring1,hit3->substring1)) {
      return (Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start) +
	(Substring_alignend_trim(hit3->substring1) - Substring_alignstart_trim(hit5->substring1))
	+ hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end;
    } else if (hit5->substring2 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring2,hit3->substring1)) {
      return (Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start) +
	(Substring_alignend_trim(hit3->substring1) - Substring_alignstart_trim(hit5->substring2))
	+ hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end;
    } else if (hit5->substring0 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring0,hit3->substring1)) {
      return /* (Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start) + */
	(Substring_alignend_trim(hit3->substring1) - Substring_alignstart_trim(hit5->substring0))
	+ hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end;
    }

    if (hit3->substring2 != NULL) {
      debug13(printf("plus general: hit3->substring2\n"));
      if (Substring_overlap_segment_trimmed_p(hit5->substring1,hit3->substring2)) {
	return (Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start) +
	  (Substring_alignend_trim(hit3->substring2) - Substring_alignstart_trim(hit5->substring1))
	  /* + hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end */;
      } else if (hit5->substring2 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring2,hit3->substring2)) {
	return (Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start) +
	  (Substring_alignend_trim(hit3->substring2) - Substring_alignstart_trim(hit5->substring2))
	  /* + hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end */;
      } else if (hit5->substring0 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring0,hit3->substring2)) {
	return /* (Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start) + */
	  (Substring_alignend_trim(hit3->substring2) - Substring_alignstart_trim(hit5->substring0))
	  /* + hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end */;
      }
    }

    if (hit3->substring0 != NULL) {
      debug13(printf("plus general: hit3->substring0\n"));
      if (Substring_overlap_segment_trimmed_p(hit5->substring1,hit3->substring0)) {
	return (Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start) +
	  (Substring_alignend_trim(hit3->substring0) - Substring_alignstart_trim(hit5->substring1))
	  + hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end;
      } else if (hit5->substring2 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring2,hit3->substring0)) {
	return (Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start) +
	  (Substring_alignend_trim(hit3->substring0) - Substring_alignstart_trim(hit5->substring2))
	  + hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end;
      } else if (hit5->substring0 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring0,hit3->substring0)) {
	return /* (Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start) + */
	  (Substring_alignend_trim(hit3->substring0) - Substring_alignstart_trim(hit5->substring0))
	  + hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end;
      }
    }
    
    return 0;

  } else {
    /* minus */
    start5 = hit5->genomicstart - hit5->trim_left - hit5->amb_nmatches_start;
    end5 = hit5->genomicend + hit5->trim_right + hit5->amb_nmatches_end;
    start3 = hit3->genomicstart - hit3->trim_left - hit3->amb_nmatches_start;
    end3 = hit3->genomicend + hit3->trim_right + hit3->amb_nmatches_end;

    if (end3 > start5) {
      /* Case 1 */
      return 0;
    } else if (end5 > start3) {
      /* Case 6 */
      return 0;
    } else if (start3 > start5) {
      if (end3 > end5) {
	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug13(printf("plus case 2a: start5 %u\n",start5 - hit5->chroffset));
	if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return (start5 - Substring_alignend_trim(hit3->substring0))
	    + hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return (start5 - Substring_alignend_trim(hit3->substring1))
	    + hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return (start5 - Substring_alignend_trim(hit3->substring2))
	    /* + hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end */;
	}

	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug13(printf("plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start +
	    (Substring_alignstart_trim(hit5->substring2) - end3);
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start +
	    (Substring_alignstart_trim(hit5->substring1) - end3);
	} else if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return /* Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start + */
	    (Substring_alignstart_trim(hit5->substring0) - end3);
	}
	/* Fall through to general algorithm */
	  
      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug13(printf("minus case 3: end5 %u\n",end5 - hit5->chroffset));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return hit5_trimmed_length + (end5 - Substring_alignend_trim(hit3->substring2))
	    /* + hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end */;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return hit5_trimmed_length + (end5 - Substring_alignend_trim(hit3->substring1))
	    + hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return hit5_trimmed_length + (end5 - Substring_alignend_trim(hit3->substring0))
	    + hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end;
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (end3 > end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug13(printf("minus case 4\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return /* (Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start) + */
	    (Substring_alignstart_trim(hit5->substring0) - start3) + hit3_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return (Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start) +
	    (Substring_alignstart_trim(hit5->substring1) - start3) + hit3_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return (Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start) +
	    (Substring_alignstart_trim(hit5->substring2) - start3) + hit3_trimmed_length;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug13(printf("minus case 5a: start3 %u\n",start3 - hit3->chroffset));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return /* (Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start) + */
	    (Substring_alignstart_trim(hit5->substring0) - start3) + hit3_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return (Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start) +
	    (Substring_alignstart_trim(hit5->substring1) - start3) + hit3_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return (Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start) +
	    (Substring_alignstart_trim(hit5->substring2) - start3) + hit3_trimmed_length;
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug13(printf("minus case 5b: end5 %u\n",end5 - hit5->chroffset));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return /* hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end + */
	    (end5 - Substring_alignend_trim(hit3->substring2)) + hit5_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end +
	    (end5 - Substring_alignend_trim(hit3->substring1)) + hit5_trimmed_length;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end +
	    (end5 - Substring_alignend_trim(hit3->substring0)) + hit5_trimmed_length;
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug13(printf("minus general: hit3->substring1\n"));
    if (Substring_overlap_segment_trimmed_p(hit5->substring1,hit3->substring1)) {
      return (Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start) +
	(Substring_alignstart_trim(hit5->substring1) - Substring_alignend_trim(hit3->substring1))
	+ hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end;
    } else if (hit5->substring2 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring2,hit3->substring1)) {
      return (Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start) +
	(Substring_alignstart_trim(hit5->substring2) - Substring_alignend_trim(hit3->substring1))
	+ hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end;
    } else if (hit5->substring0 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring0,hit3->substring1)) {
      return /* (Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start) + */
	(Substring_alignstart_trim(hit5->substring0) - Substring_alignend_trim(hit3->substring1))
	+ hit3->querylength - Substring_queryend(hit3->substring1) - hit3->trim_right - hit3->amb_nmatches_end;
    }

    if (hit3->substring2 != NULL) {
      debug13(printf("minus general: hit3->substring2\n"));
      if (Substring_overlap_segment_trimmed_p(hit5->substring1,hit3->substring2)) {
	return (Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start) +
	  (Substring_alignstart_trim(hit5->substring1) - Substring_alignend_trim(hit3->substring2))
	  /* + hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end */;
      } else if (hit5->substring2 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring2,hit3->substring2)) {
	return (Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start) +
	  (Substring_alignstart_trim(hit5->substring2) - Substring_alignend_trim(hit3->substring2))
	  /* + hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end */;
      } else if (hit5->substring0 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring0,hit3->substring2)) {
	return /* (Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start) + */
	  (Substring_alignstart_trim(hit5->substring0) - Substring_alignend_trim(hit3->substring2))
	  /* + hit3->querylength - Substring_queryend(hit3->substring2) - hit3->trim_right - hit3->amb_nmatches_end */;
      }
    }

    if (hit3->substring0 != NULL) {
      debug13(printf("minus general: hit3->substring0\n"));
      if (Substring_overlap_segment_trimmed_p(hit5->substring1,hit3->substring0)) {
	return (Substring_querystart(hit5->substring1) - hit5->trim_left - hit5->amb_nmatches_start) +
	  (Substring_alignstart_trim(hit5->substring1) - Substring_alignend_trim(hit3->substring0))
	  + hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end;
      } else if (hit5->substring2 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring2,hit3->substring0)) {
	return (Substring_querystart(hit5->substring2) - hit5->trim_left - hit5->amb_nmatches_start) +
	  (Substring_alignstart_trim(hit5->substring2) - Substring_alignend_trim(hit3->substring0))
	  + hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end;
      } else if (hit5->substring0 != NULL && Substring_overlap_segment_trimmed_p(hit5->substring0,hit3->substring0)) {
	return /* (Substring_querystart(hit5->substring0) - hit5->trim_left - hit5->amb_nmatches_start) + */
	  (Substring_alignstart_trim(hit5->substring0) - Substring_alignend_trim(hit3->substring0))
	  + hit3->querylength - Substring_queryend(hit3->substring0) - hit3->trim_right - hit3->amb_nmatches_end;
      }
    }
	
    return 0;
  }
}



int
Stage3pair_overlap (int *hardclip5, int *hardclip3, Stage3pair_T this) {
  Stage3end_T hit5, hit3;
  int totallength, insertlength, minlength;
  int overlap;
  int clipdir;
  int hit5_trimmed_length, hit3_trimmed_length;

  hit5 = this->hit5;
  hit3 = this->hit3;

  debug13(printf("Entered Stage3pair_overlap with hittype %s and %s\n",
		 hittype_string(hit5->hittype),hittype_string(hit3->hittype)));
  if (hit5->hittype == SAMECHR_SPLICE || hit5->hittype == TRANSLOC_SPLICE) {
    return 0;
  } else if (hit3->hittype == SAMECHR_SPLICE || hit3->hittype == TRANSLOC_SPLICE) {
    return 0;
  } else {
    debug13(printf("hit5 trim_left %d + amb_start %d, trim_right %d + amb_end %d, hit3 trim_left %d + amb_start %d, trim_right %d + amb_end %d\n",
		   hit5->trim_left,hit5->amb_nmatches_start,hit5->trim_right,hit5->amb_nmatches_end,
		   hit3->trim_left,hit3->amb_nmatches_start,hit3->trim_right,hit3->amb_nmatches_end));
    if (hit5->plusp == true) {
      hit5_trimmed_length = hit5->querylength - hit5->trim_left - hit5->trim_right - hit5->amb_nmatches_start - hit5->amb_nmatches_end;
      hit3_trimmed_length = hit3->querylength - hit3->trim_left - hit3->trim_right - hit3->amb_nmatches_start - hit3->amb_nmatches_end;
      totallength = hit5_trimmed_length + hit3_trimmed_length;
      debug13(printf("totallength = %d, hit5 trimmed length = %d, hit3 trimmed length = %d\n",
		     totallength,hit5_trimmed_length,hit3_trimmed_length));

#if 0
      hit5_inside = (hit5->genomicend - 1) - hit5->trim_right;
      hit3_inside = hit3->genomicstart + hit3->trim_left; /* because hit3 is inverted */
      insertlength = (hit3_inside + hit3_trimmed_length - 1) - (hit5_inside - hit5_trimmed_length + 1) + 1;
#endif
      debug13(printf("original insertlength: %d, trim5: %d..%d, trim3: %d..%d\n",
		     this->insertlength,hit5->trim_left,hit5->trim_right,hit3->trim_left,hit3->trim_right));
      if (hit5->hittype == GMAP || hit3->hittype == GMAP) {
	insertlength = this->insertlength = pair_insert_length_trimmed(hit5,hit3,hit5_trimmed_length,hit3_trimmed_length);
      } else {
	insertlength = pair_insert_length_trimmed(hit5,hit3,hit5_trimmed_length,hit3_trimmed_length);
      }
      debug13(printf(" => trimmed insertlength: %d\n",insertlength));

      if (insertlength == 0) {
	debug13(printf("Cannot determine an insert length, so returning 0\n"));
	return 0;
      } else if ((overlap = totallength - insertlength) < 0) {
	debug13(printf("Overlap %d is negative, so returning 0\n",overlap));
	return 0;
      } else {
	debug13(printf("overlap = %d\n",overlap));
      }

      if (hit5_trimmed_length <= hit3_trimmed_length) {
	debug13(printf("Overlap is %d and hit5 is shortest\n",overlap));
	minlength = hit5_trimmed_length;
	if (insertlength < minlength) {
	  /* Clip tails equally */
	  *hardclip5 = insertlength/2;
	  *hardclip3 = insertlength - (*hardclip5);
	  *hardclip3 += hit3->trim_right;
	  *hardclip5 += hit5->trim_left;
	  debug13(printf("Clip tails with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = -1;
	} else if (overlap > minlength) {
	  /* Clip heads, but limit clipping on hit5 */
	  *hardclip5 = hit5_trimmed_length/2;
	  *hardclip3 = overlap - (*hardclip5);
	  *hardclip3 += hit3->trim_left;
	  *hardclip5 += hit5->trim_right;
	  debug13(printf("Clip heads limited with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = +1;
	} else {
	  /* Clip heads equally */
	  *hardclip5 = overlap/2;
	  *hardclip3 = overlap - (*hardclip5);
	  *hardclip3 += hit3->trim_left;
	  *hardclip5 += hit5->trim_right;
	  debug13(printf("Clip heads with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = +1;
	}

      } else {
	debug13(printf("Overlap is %d and hit3 is shortest\n",overlap));
	minlength = hit3_trimmed_length;
	if (insertlength < minlength) {
	  /* Clip tails equally */
	  *hardclip3 = insertlength/2;
	  *hardclip5 = insertlength - (*hardclip3);
	  *hardclip3 += hit3->trim_right;
	  *hardclip5 += hit5->trim_left;
	  debug13(printf("Clip tails with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = -1;
	} else if (overlap > minlength) {
	  /* Clip heads, but limit clipping on hit3 */
	  *hardclip3 = hit3_trimmed_length/2;
	  *hardclip5 = overlap - (*hardclip3);
	  *hardclip3 += hit3->trim_left;
	  *hardclip5 += hit5->trim_right;
	  debug13(printf("Clip heads limited with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = +1;
	} else {
	  /* Clip heads equally */
	  *hardclip3 = overlap/2;
	  *hardclip5 = overlap - (*hardclip3);
	  *hardclip3 += hit3->trim_left;
	  *hardclip5 += hit5->trim_right;
	  debug13(printf("Clip heads with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = +1;
	}
      }

      debug13(printf("plus, so returning clipdir %d\n",clipdir));
      return clipdir;

    } else {
      hit5_trimmed_length = hit5->querylength - hit5->trim_left - hit5->trim_right - hit5->amb_nmatches_start - hit5->amb_nmatches_end;
      hit3_trimmed_length = hit3->querylength - hit3->trim_left - hit3->trim_right - hit3->amb_nmatches_start - hit3->amb_nmatches_end;
      totallength = hit5_trimmed_length + hit3_trimmed_length;
      debug13(printf("totallength = %d, hit5 trimmed length = %d, hit3 trimmed length = %d\n",
		     totallength,hit5_trimmed_length,hit3_trimmed_length));
    
#if 0
      hit5_inside = hit5->genomicend + hit5->trim_right;
      hit3_inside = (hit3->genomicstart - 1) - hit3->trim_left; /* because hit3 is inverted */
      insertlength = (hit5_inside + hit5_trimmed_length - 1) - (hit3_inside - hit3_trimmed_length + 1) + 1;
#endif
      debug13(printf("original insertlength: %d, trim5: %d..%d, trim3: %d..%d\n",
		     this->insertlength,hit5->trim_left,hit5->trim_right,hit3->trim_left,hit3->trim_right));
      if (hit5->hittype == GMAP || hit3->hittype == GMAP) {
	insertlength = this->insertlength = pair_insert_length_trimmed(hit5,hit3,hit5_trimmed_length,hit3_trimmed_length);
      } else {
	insertlength = pair_insert_length_trimmed(hit5,hit3,hit5_trimmed_length,hit3_trimmed_length);
      }
      debug13(printf(" => trimmed insertlength: %d\n",insertlength));

      if (insertlength == 0) {
	debug13(printf("Cannot determine an insert length, so returning 0\n"));
	return 0;
      } else if ((overlap = totallength - insertlength) < 0) {
	debug13(printf("Overlap %d is negative, so returning 0\n",overlap));
	return 0;
      } else {
	debug13(printf("overlap = %d\n",overlap));
      }

      if (hit5_trimmed_length <= hit3_trimmed_length) {
	debug13(printf("Overlap is %d and hit5 is shortest\n",overlap));
	minlength = hit5_trimmed_length;
	if (insertlength < minlength) {
	  /* Clip tails equally */
	  *hardclip5 = insertlength/2;
	  *hardclip3 = insertlength - (*hardclip5);
	  *hardclip3 += hit3->trim_right;
	  *hardclip5 += hit5->trim_left;
	  debug13(printf("Clip tails with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = -1;
	} else if (overlap > minlength) {
	  /* Clip heads, but limit clipping on hit5 */
	  *hardclip5 = hit5_trimmed_length/2;
	  *hardclip3 = overlap - (*hardclip5);
	  *hardclip3 += hit3->trim_left;
	  *hardclip5 += hit5->trim_right;
	  debug13(printf("Clip heads limited with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = +1;
	} else {
	  /* Clip heads equally */
	  *hardclip5 = overlap/2;
	  *hardclip3 = overlap - (*hardclip5);
	  *hardclip3 += hit3->trim_left;
	  *hardclip5 += hit5->trim_right;
	  debug13(printf("Clip heads with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = +1;
	}

      } else {
	debug13(printf("Overlap is %d and hit3 is shortest\n",overlap));
	minlength = hit3_trimmed_length;
	if (insertlength < minlength) {
	  /* Clip tails equally */
	  *hardclip3 = insertlength/2;
	  *hardclip5 = insertlength - (*hardclip3);
	  *hardclip3 += hit3->trim_right;
	  *hardclip5 += hit5->trim_left;
	  debug13(printf("Clip tails with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = -1;
	} else if (overlap > minlength) {
	  /* Clip heads, but limit clipping on hit3 */
	  *hardclip3 = hit3_trimmed_length/2;
	  *hardclip5 = overlap - (*hardclip3);
	  *hardclip3 += hit3->trim_left;
	  *hardclip5 += hit5->trim_right;
	  debug13(printf("Clip heads limited with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = +1;
	} else {
	  /* Clip heads equally */
	  *hardclip3 = overlap/2;
	  *hardclip5 = overlap - (*hardclip3);
	  *hardclip3 += hit3->trim_left;
	  *hardclip5 += hit5->trim_right;
	  debug13(printf("Clip heads with hardclip5 %d and hardclip3 %d\n",*hardclip5,*hardclip3));
	  clipdir = +1;
	}
      }
      debug13(printf("minus, so returning -clipdir %d\n",-clipdir));
      return -clipdir;
    }
  }
}


void
Stage3pair_set_private5p (Stage3pair_T this) {
  this->private5p = true;
  return;
}

void
Stage3pair_clear_private5p (Stage3pair_T this) {
  this->private5p = false;
  return;
}

void
Stage3pair_set_private3p (Stage3pair_T this) {
  this->private3p = true;
  return;
}

void
Stage3pair_clear_private3p (Stage3pair_T this) {
  this->private3p = false;
  return;
}


void
Stage3pair_free (Stage3pair_T *old) {
  debug0(printf("Freeing pair %p with hits %p (privatep %d) and %p (privatep %d)\n",*old,(*old)->hit5,(*old)->private5p,(*old)->hit3,(*old)->private3p));
  if ((*old)->private3p == true) {
    assert((*old)->hit3 != NULL);
    debug0(printf("Freeing end3 at %p\n",(*old)->hit3));
    Stage3end_free(&(*old)->hit3);
  }
  if ((*old)->private5p == true) {
    assert((*old)->hit5 != NULL);
    debug0(printf("Freeing end5 at %p\n",(*old)->hit5));
    Stage3end_free(&(*old)->hit5);
  }
  FREE_OUT(*old);
  return;
}
  

static Overlap_T
Stage3pair_gene_overlap (Stage3pair_T this) {
  Overlap_T overlap;
  bool foundp = false;

  if (genes_iit == NULL) {
    return NO_KNOWN_GENE;
  } else {
    if ((overlap = Stage3end_gene_overlap(this->hit5)) == KNOWN_GENE_MULTIEXON) {
      return KNOWN_GENE_MULTIEXON;
    } else if (overlap == KNOWN_GENE) {
      if (favor_multiexon_p == false) {
	return KNOWN_GENE;
      } else {
	foundp = true;
      }
    }

    if ((overlap = Stage3end_gene_overlap(this->hit3)) == KNOWN_GENE_MULTIEXON) {
      return KNOWN_GENE_MULTIEXON;
    } else if (overlap == KNOWN_GENE) {
      if (favor_multiexon_p == false) {
	return KNOWN_GENE;
      } else {
	foundp = true;
      }
    }

    if (foundp == true) {
      return KNOWN_GENE;
    } else {
      return NO_KNOWN_GENE;
    }
  }
}

#if 0
static long int
Stage3pair_tally (Stage3pair_T this) {

  if (tally_iit == NULL) {
    return 0L;
  } else if (this->tally >= 0) {
    return this->tally;
  } else {
    this->tally = Stage3end_compute_tally(this->hit5) + Stage3end_compute_tally(this->hit3);
    return this->tally;
  }
}
#endif


#if 0
static char complCode[128] = COMPLEMENT_LC;

static char *
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC_OUT(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return complement;
}
#endif


const Except_T Copy_Substring = { "Substring invalid during copy" };

T
Stage3end_copy (T old) {
  T new = (T) MALLOC_OUT(sizeof(*new));

  debug0(printf("Copying Stage3end %p -> %p of type %s\n",
		old,new,hittype_string(old->hittype)));

  new->hittype = old->hittype;
  new->genestrand = old->genestrand;
  new->sarrayp = old->sarrayp;
  new->improved_by_gmap_p = old->improved_by_gmap_p;

  new->chrnum = old->chrnum;
  new->effective_chrnum = old->effective_chrnum;
  new->other_chrnum = old->other_chrnum;
  new->chroffset = old->chroffset;
  new->chrhigh = old->chrhigh;
  new->chrlength = old->chrlength;

  new->querylength = old->querylength;
  new->querylength_adj = old->querylength_adj;

  new->genomicstart = old->genomicstart;
  new->genomicend = old->genomicend;
  new->plusp = old->plusp;

  new->low = old->low;
  new->high = old->high;
  new->genomiclength = old->genomiclength;
  new->guided_insertlength = old->guided_insertlength;

  new->mapq_loglik = old->mapq_loglik;
  new->mapq_score = old->mapq_score;
  new->absmq_score = old->absmq_score;

  new->score = old->score;
  new->ntscore = old->ntscore;
  new->nmatches = old->nmatches;
  new->nmatches_posttrim = old->nmatches_posttrim;
  new->gmap_max_match_length = old->gmap_max_match_length;
  new->gmap_min_splice_prob = old->gmap_min_splice_prob;

  new->trim_left = old->trim_left;
  new->trim_right = old->trim_right;
  new->trim_left_splicep = old->trim_left_splicep;
  new->trim_right_splicep = old->trim_right_splicep;

  new->penalties = old->penalties;
  new->score_eventrim = old->score_eventrim;

  new->gene_overlap = old->gene_overlap;
  new->tally = old->tally;

  new->nmismatches_whole = old->nmismatches_whole;
  new->nmismatches_bothdiff = old->nmismatches_bothdiff;
  new->nmismatches_refdiff = old->nmismatches_refdiff;

  new->nindels = old->nindels;
  new->indel_pos = old->indel_pos;
  new->indel_low = old->indel_low;
  if (old->deletion == NULL) {
    new->deletion = (char *) NULL;
  } else {
    new->deletion = (char *) CALLOC_OUT(strlen(old->deletion)+1,sizeof(char));
    strcpy(new->deletion,old->deletion);
  }
  
  new->distance = old->distance;
  new->acceptor_distance = old->acceptor_distance;
  new->donor_distance = old->donor_distance;

  new->gmap_nindelbreaks = old->gmap_nindelbreaks;
  new->gmap_cdna_direction = old->gmap_cdna_direction;
  new->gmap_nintrons = old->gmap_nintrons;
  new->sensedir = old->sensedir;
  new->sensedir_nonamb = old->sensedir_nonamb;

  new->start_ambiguous_p = old->start_ambiguous_p;
  new->end_ambiguous_p = old->end_ambiguous_p;
  new->nchimera_known = old->nchimera_known;
  new->nchimera_novel = old->nchimera_novel;

  new->amb_nmatches_start = old->amb_nmatches_start;
  new->amb_nmatches_end = old->amb_nmatches_end;
  new->amb_nmatches_donor = old->amb_nmatches_donor;
  new->amb_nmatches_acceptor = old->amb_nmatches_acceptor;


  new->gmap_start_endtype = old->gmap_start_endtype;
  new->gmap_end_endtype = old->gmap_end_endtype;

  if ((new->start_nambi = old->start_nambi) == 0) {
    new->start_ambi = (int *) NULL;
    new->start_amb_nmismatches = (int *) NULL;
  } else {
    new->start_ambi = (int *) CALLOC_OUT(old->start_nambi,sizeof(int));
    memcpy(new->start_ambi,old->start_ambi,old->start_nambi*sizeof(int));
    new->start_amb_nmismatches = (int *) CALLOC_OUT(old->start_nambi,sizeof(int));
    memcpy(new->start_amb_nmismatches,old->start_amb_nmismatches,old->start_nambi*sizeof(int));
  }

  if ((new->end_nambi = old->end_nambi) == 0) {
    new->end_ambi = (int *) NULL;
    new->end_amb_nmismatches = (int *) NULL;
  } else {
    new->end_ambi = (int *) CALLOC_OUT(old->end_nambi,sizeof(int));
    memcpy(new->end_ambi,old->end_ambi,old->end_nambi*sizeof(int));
    new->end_amb_nmismatches = (int *) CALLOC_OUT(old->end_nambi,sizeof(int));
    memcpy(new->end_amb_nmismatches,old->end_amb_nmismatches,old->end_nambi*sizeof(int));
  }

  new->substring1 = Substring_copy(old->substring1);
  new->substring2 = Substring_copy(old->substring2);
  new->substring0 = Substring_copy(old->substring0);

  if (old->hittype == GMAP) {
    new->pairarray = Pairpool_copy_array(old->pairarray,old->npairs);
    new->npairs = old->npairs;
    new->nsegments = old->nsegments;
  } else {
    new->pairarray = (struct Pair_T *) NULL;
    new->npairs = 0;
    new->nsegments = 0;
  }

  if (old->substring_donor == NULL) {
    new->substring_donor = NULL;
  } else if (old->substring_donor == old->substring1) {
    new->substring_donor = new->substring1;
  } else if (old->substring_donor == old->substring2) {
    new->substring_donor = new->substring2;
  } else {
    fprintf(stderr,"substring_donor for type %s is not NULL, substring1, or substring2\n",
	    hittype_string(old->hittype));
    fprintf(stderr,"substring_donor %p\n",old->substring_donor);
    fprintf(stderr,"substring1 %p\n",old->substring1);
    fprintf(stderr,"substring2 %p\n",old->substring2);
    Except_raise(&Copy_Substring, __FILE__, __LINE__);
  }

  if (old->substring_acceptor == NULL) {
    new->substring_acceptor = NULL;
  } else if (old->substring_acceptor == old->substring1) {
    new->substring_acceptor = new->substring1;
  } else if (old->substring_acceptor == old->substring2) {
    new->substring_acceptor = new->substring2;
  } else {
    fprintf(stderr,"substring_acceptor for type %s is not NULL, substring1, or substring2\n",
	    hittype_string(old->hittype));
    fprintf(stderr,"substring_acceptor %p\n",old->substring_acceptor);
    fprintf(stderr,"substring1 %p\n",old->substring1);
    fprintf(stderr,"substring2 %p\n",old->substring2);
    Except_raise(&Copy_Substring, __FILE__, __LINE__);
  }

  if (old->substringD == NULL) {
    new->substringD = NULL;
  } else if (old->substringD == old->substring0) {
    new->substringD = new->substring0;
  } else if (old->substringD == old->substring2) {
    new->substringD = new->substring2;
  } else {
    fprintf(stderr,"substringD for type %s is not NULL, substring0, or substring2\n",
	    hittype_string(old->hittype));
    fprintf(stderr,"substringD %p\n",old->substringD);
    fprintf(stderr,"substring0 %p\n",old->substring0);
    fprintf(stderr,"substring2 %p\n",old->substring2);
    Except_raise(&Copy_Substring, __FILE__, __LINE__);
  }
  
  if (old->substringA == NULL) {
    new->substringA = NULL;
  } else if (old->substringA == old->substring0) {
    new->substringA = new->substring0;
  } else if (old->substringA == old->substring2) {
    new->substringA = new->substring2;
  } else {
    fprintf(stderr,"substringA for type %s is not NULL, substring0, or substring2\n",
	    hittype_string(old->hittype));
    fprintf(stderr,"substringA %p\n",old->substringA);
    fprintf(stderr,"substring0 %p\n",old->substring0);
    fprintf(stderr,"substring2 %p\n",old->substring2);
    Except_raise(&Copy_Substring, __FILE__, __LINE__);
  }


  if (old->substring_low == NULL) {
    new->substring_low = NULL;
  } else if (old->substring_low == old->substring1) {
    new->substring_low = new->substring1;
  } else if (old->substring_low == old->substring2) {
    new->substring_low = new->substring2;
  } else if (old->substring_low == old->substring0) {
    new->substring_low = new->substring0;
  } else {
    fprintf(stderr,"substring_low for type %s is not NULL, substring1, or substring2, or substring0\n",
	    hittype_string(old->hittype));
    fprintf(stderr,"substring_low %p\n",old->substring_low);
    fprintf(stderr,"substring1 %p\n",old->substring1);
    fprintf(stderr,"substring2 %p\n",old->substring2);
    fprintf(stderr,"substring0 %p\n",old->substring0);
    Except_raise(&Copy_Substring, __FILE__, __LINE__);
  }

  if (old->substring_high == NULL) {
    new->substring_high = NULL;
  } else if (old->substring_high == old->substring1) {
    new->substring_high = new->substring1;
  } else if (old->substring_high == old->substring2) {
    new->substring_high = new->substring2;
  } else if (old->substring_high == old->substring0) {
    new->substring_high = new->substring0;
  } else {
    fprintf(stderr,"substring_high for type %s is not NULL, substring1, or substring2, or substring0\n",
	    hittype_string(old->hittype));
    fprintf(stderr,"substring_high %p\n",old->substring_high);
    fprintf(stderr,"substring1 %p\n",old->substring1);
    fprintf(stderr,"substring2 %p\n",old->substring2);
    fprintf(stderr,"substring0 %p\n",old->substring0);
    Except_raise(&Copy_Substring, __FILE__, __LINE__);
  }

  new->paired_usedp = old->paired_usedp;
  new->paired_seenp = old->paired_seenp;
  new->concordantp = old->concordantp;

  new->alias = old->alias;
  new->circularpos = old->circularpos;

  return new;
}


static int
compute_circularpos (int *alias, T hit) {
  int circularpos;

  if (circularp[hit->chrnum] == false) {
    /* This also handles hit->chrnum == 0, where translocation cannot be circular */
    *alias = 0;
    return -1;

  } else if (hit->low >= hit->chroffset + hit->chrlength) {
    /* All in circular alias */
    *alias = +1;
    return -1;

  } else if (hit->high <= hit->chroffset + hit->chrlength) {
    /* All in circular proper */
    *alias = -1;
    return -1;

  } else if (hit->hittype == GMAP) {
    *alias = 0;
    return Pair_circularpos(hit->pairarray,hit->npairs,hit->chrlength,
			    hit->plusp,hit->querylength_adj);

  } else if (hit->plusp == true) {
    *alias = 0;
    if ((circularpos = Substring_circularpos(hit->substring0)) > 0) {
      return circularpos;
    } else if ((circularpos = Substring_circularpos(hit->substring1)) > 0) {
      return circularpos;
    } else if ((circularpos = Substring_circularpos(hit->substring2)) > 0) {
      return circularpos;
    } else {
      return -1;
    }

  } else {
    *alias = 0;
    if ((circularpos = Substring_circularpos(hit->substring2)) > 0) {
      return circularpos;
    } else if ((circularpos = Substring_circularpos(hit->substring1)) > 0) {
      return circularpos;
    } else if ((circularpos = Substring_circularpos(hit->substring0)) > 0) {
      return circularpos;
    } else {
      return -1;
    }
  }
}


T
Stage3end_new_exact (int *found_score, Univcoord_T left, int genomiclength, Compress_T query_compress,
		     bool plusp, int genestrand, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     Chrpos_T chrlength, bool sarrayp) {
  T new;
  Substring_T substring;
  Univcoord_T genomicstart, genomicend;

  if (plusp == true) {
    if ((genomicend = left + genomiclength) > chrhigh) {
      return (T) NULL;
    } else {
      genomicstart = left;
    }
  } else {
    if ((genomicstart = left + genomiclength) > chrhigh) {
      return (T) NULL;
    } else {
      genomicend = left;
    }
  }

  if ((substring = Substring_new(/*nmismatches*/0,chrnum,chroffset,chrhigh,chrlength,
				 left,genomicstart,genomicend,query_compress,
				 /*start_endtype*/END,/*end_endtype*/END,
				 /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				 /*alignstart*/genomicstart,/*alignend*/genomicend,
				 genomiclength,/*extraleft*/0,/*extraright*/0,/*exactp*/true,plusp,genestrand,
				 /*trim_left_p*/false,/*trim_right_p*/false,/*minlength*/0)) == NULL) {
    return (T) NULL;

  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
    debug0(printf("Stage3end_new_exact %p: left %lu, chrnum %d\n",new,left,chrnum));

    new->substring1 = substring;
    new->substring2 = (Substring_T) NULL;
    new->substring0 = (Substring_T) NULL;
    new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
    new->substringD = new->substringA = (Substring_T) NULL;
    new->substring_low = new->substring_high = new->substring1;

    new->pairarray = (struct Pair_T *) NULL;

    new->deletion = (char *) NULL;
    new->querylength_adj = new->querylength = genomiclength;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    if (genomicstart < genomicend) {
      new->low = genomicstart;
      new->high = genomicend;
    } else {
      new->low = genomicend;
      new->high = genomicstart;
    }
    new->genomiclength = new->high - new->low;
    new->guided_insertlength = 0U;
    debug0(printf("Assigned %lu to low and %lu to high\n",new->low,new->high));


    new->hittype = EXACT;
    new->genestrand = genestrand;
    new->sarrayp = sarrayp;
    new->improved_by_gmap_p = false;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
    new->chrlength = chrlength;
    new->plusp = plusp;
    new->sensedir = new->sensedir_nonamb = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring);
    new->mapq_score = 0;
    new->absmq_score = 0;
#endif

    new->nindels = 0;
    new->indel_pos = 0;
    new->indel_low = 0;
    new->nmismatches_whole = 0;
    new->nmismatches_bothdiff = 0;
    /* new->nmismatches_refdiff = 0; */
    new->ntscore = 0;
    new->score = 0;
    new->nmatches = genomiclength;
    new->nmatches_posttrim = genomiclength;

    new->trim_left = 0;
    new->trim_right = 0;
    new->trim_left_splicep = false;
    new->trim_right_splicep = false;

    new->penalties = 0;

    /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
    new->tally = -1L;
    *found_score = 0;

    new->amb_nmatches_start = new->amb_nmatches_end = 0;

    new->start_ambiguous_p = new->end_ambiguous_p = false;
    new->start_ambi = new->end_ambi = (int *) NULL;
    new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
    new->start_nambi = new->end_nambi = 0;
    new->nchimera_known = 0;
    new->nchimera_novel = 0;

    new->distance = 0U;
    new->acceptor_distance = new->donor_distance = 0U;

    new->paired_usedp = false;
    new->paired_seenp = false;
    new->concordantp = false;

    new->circularpos = compute_circularpos(&new->alias,new);

    return new;
  }
}


T
Stage3end_new_substitution (int *found_score, int nmismatches_whole, Univcoord_T left,
			    int genomiclength, Compress_T query_compress,
			    bool plusp, int genestrand, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			    Chrpos_T chrlength, bool sarrayp) {
  T new;
  Substring_T substring;
  Univcoord_T genomicstart, genomicend;

  if (plusp == true) {
    if ((genomicend = left + genomiclength) > chrhigh) {
      return (T) NULL;
    } else {
      genomicstart = left;
    }
  } else {
    if ((genomicstart = left + genomiclength) > chrhigh) {
      return (T) NULL;
    } else {
      genomicend = left;
    }
  }

  if ((substring = Substring_new(nmismatches_whole,chrnum,chroffset,chrhigh,chrlength,
				 left,genomicstart,genomicend,query_compress,
				 /*start_endtype*/END,/*end_endtype*/END,
				 /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				 /*alignstart*/genomicstart,/*alignend*/genomicend,
				 genomiclength,/*extraleft*/0,/*extraright*/0,/*exactp*/false,plusp,genestrand,
				 /*trim_left_p*/true,/*trim_right_p*/true,
				 /*minlength*/genomiclength/2)) == NULL) {
    return (T) NULL;

  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
    debug0(printf("Stage3end_new_substitution %p: left %lu, chrnum %d, nmismatches %d\n",
		  new,left,chrnum,nmismatches_whole));

    new->substring1 = substring;
    new->substring2 = (Substring_T) NULL;
    new->substring0 = (Substring_T) NULL;
    new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
    new->substringD = new->substringA = (Substring_T) NULL;
    new->substring_low = new->substring_high = new->substring1;

    new->pairarray = (struct Pair_T *) NULL;

    new->deletion = (char *) NULL;
    new->querylength_adj = new->querylength = genomiclength;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    if (genomicstart < genomicend) {
      new->low = genomicstart;
      new->high = genomicend;
    } else {
      new->low = genomicend;
      new->high = genomicstart;
    }
    new->genomiclength = new->high - new->low;
    new->guided_insertlength = 0U;

    if (nmismatches_whole == 0) {
      /* Proper hittype needed so we can eliminate identical hits */
      new->hittype = EXACT;
    } else {
      new->hittype = SUB;
    }
    new->genestrand = genestrand;
    new->sarrayp = sarrayp;
    new->improved_by_gmap_p = false;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
    new->chrlength = chrlength;
    new->plusp = plusp;
    new->sensedir = new->sensedir_nonamb = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring);
    new->mapq_score = 0;
    new->absmq_score = 0;
#endif

    new->nindels = 0;
    new->indel_pos = 0;
    new->indel_low = 0;
    new->nmismatches_whole = nmismatches_whole;
    new->ntscore = nmismatches_whole;
    new->score = nmismatches_whole;

    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(new->substring1);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1); */

#if 0
    /* This method was previously the only one correct for SNP-tolerant alignment */
    new->nmatches = Substring_match_length(new->substring1) - new->total_nmismatches;
#else
    /* This method is now correct for SNP-tolerant alignment */
    new->nmatches = Substring_nmatches(new->substring1);
    new->nmatches_posttrim = Substring_nmatches_posttrim(new->substring1);
#endif

    new->trim_left = Substring_trim_left(substring);
    new->trim_right = Substring_trim_right(substring);
    new->trim_left_splicep = Substring_trim_left_splicep(substring);
    new->trim_right_splicep = Substring_trim_right_splicep(substring);

    new->penalties = 0;

    /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
    new->tally = -1L;

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->amb_nmatches_start = new->amb_nmatches_end = 0;

    new->start_ambiguous_p = new->end_ambiguous_p = false;
    new->start_ambi = new->end_ambi = (int *) NULL;
    new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
    new->start_nambi = new->end_nambi = 0;
    new->nchimera_known = 0;
    new->nchimera_novel = 0;

    new->distance = 0U;
    new->acceptor_distance = new->donor_distance = 0U;

    new->paired_usedp = false;
    new->paired_seenp = false;
    new->concordantp = false;

    new->circularpos = compute_circularpos(&new->alias,new);

    return new;
  }
}



T
Stage3end_new_insertion (int *found_score, int nindels, int indel_pos, int nmismatches1_whole, int nmismatches2_whole,
			 Univcoord_T left, int genomiclength, Compress_T query_compress,
			 int querylength, bool plusp, int genestrand, Chrnum_T chrnum, Univcoord_T chroffset,
			 Univcoord_T chrhigh, Chrpos_T chrlength, int indel_penalty, bool sarrayp) {
  T new;
  Substring_T substring1, substring2;
  int querystart1, queryend1, querystart2, queryend2;
  Univcoord_T genomicstart, genomicend;
  Univcoord_T alignstart1, alignend1, alignstart2, alignend2;

  debug2(printf("Entered with left %lu, querylength %d, genomiclength %d, indel_pos %d\n",
		left,querylength,genomiclength,indel_pos));
  debug2(printf("q: %s\n",query));
#if 0
  debug2(printf("g: %s\n",genomicseg));
#endif

  assert(nindels > 0);

  querystart1 = 0;
  queryend1 = indel_pos;
  querystart2 = indel_pos + nindels;
  queryend2 = querylength;

  if (plusp == true) {
    if ((genomicend = left + genomiclength) > chrhigh) {
      return (T) NULL;
    } else {
      genomicstart = left;

      alignstart1 = genomicstart;
      alignend1 = alignstart2 = genomicstart + indel_pos;
      alignend2 = genomicend/* - nindels*/;
    }

  } else {
    if ((genomicstart = left + genomiclength) > chrhigh) {
      return (T) NULL;
    } else {
      genomicend = left;

      alignstart1 = genomicstart;
      alignend1 = alignstart2 = genomicstart - indel_pos;
      alignend2 = genomicend/* + nindels*/;
    }
  }

  if ((substring1 = Substring_new(nmismatches1_whole,chrnum,chroffset,chrhigh,chrlength,
				  left,genomicstart,genomicend,query_compress,
				  /*start_endtype*/END,/*end_endtype*/INS,
				  querystart1,queryend1,querylength,alignstart1,alignend1,genomiclength,
				  /*extraleft*/0,/*extraright*/0,/*exactp*/false,plusp,genestrand,
				  /*trim_left_p (previously was end1_indel_p ? false : true)*/true,
				  /*trim_right_p*/false,/*minlength*/0)) == NULL) {
    return (T) NULL;

  } else if ((substring2 = Substring_new(nmismatches2_whole,chrnum,chroffset,chrhigh,chrlength,
					 left,genomicstart,genomicend,query_compress,
					 /*start_endtype*/INS,/*end_endtype*/END,
					 querystart2,queryend2,querylength,alignstart2,alignend2,genomiclength,
					 /*extraleft*/0,/*extraright*/0,/*exactp*/false,plusp,genestrand,
					 /*trim_left_p*/false,
					 /*trim_right_p (previously was end2_indel_p ? false : true)*/true,
					 /*minlength*/0)) == NULL) {
    Substring_free(&substring1);
    return (T) NULL;

  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
    debug0(printf("Stage3end_new_insertion %p: left %lu, chrnum %d, nmismatches %d+%d, indel_pos %d, nindels %d\n",
		  new,left,chrnum,nmismatches1_whole,nmismatches2_whole,indel_pos,nindels));

    new->substring1 = substring1;
    new->substring2 = substring2;
    new->substring0 = (Substring_T) NULL;
    new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
    new->substringD = new->substringA = (Substring_T) NULL;

    new->indel_pos = indel_pos;
    if (plusp == true) {
      new->substring_low = new->substring1;
      new->substring_high = new->substring2;
      new->indel_low = indel_pos;
    } else {
      new->substring_low = new->substring2;
      new->substring_high = new->substring1;
      new->indel_low = querylength - indel_pos;
    }

    new->pairarray = (struct Pair_T *) NULL;

    new->deletion = (char *) NULL;
    new->querylength_adj = new->querylength = querylength /* - nindels */;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    if (genomicstart < genomicend) {
      new->low = genomicstart;
      new->high = genomicend;
    } else {
      new->low = genomicend;
      new->high = genomicstart;
    }
    new->genomiclength = new->high - new->low;
    new->guided_insertlength = 0U;

    new->hittype = INSERTION;
    new->genestrand = genestrand;
    new->sarrayp = sarrayp;
    new->improved_by_gmap_p = false;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
    new->chrlength = chrlength;
    new->plusp = plusp;
    new->sensedir = new->sensedir_nonamb = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring1) + Substring_mapq_loglik(substring2) + 
      MAPQ_loglik_exact(quality_string,queryend1,querystart2);
    new->mapq_score = 0;
    new->absmq_score = 0;
#endif

    new->nindels = nindels;
    new->nmismatches_whole = nmismatches1_whole + nmismatches2_whole;
    new->ntscore = indel_penalty + nmismatches1_whole + nmismatches2_whole;
    new->score = new->ntscore;

    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(new->substring1) + Substring_nmismatches_bothdiff(new->substring2);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1) + Substring_nmismatches_refdiff(new->substring2); */

#if 0
    /* This method is correct for SNP-tolerant alignment */
    new->nmatches = Substring_match_length(new->substring1) + Substring_match_length(new->substring2) - new->total_nmismatches;
#else
    /* This method is now correct for SNP-tolerant alignment */
    new->nmatches = Substring_nmatches(new->substring1) + Substring_nmatches(new->substring2);
    new->nmatches_posttrim = Substring_nmatches_posttrim(new->substring1) + Substring_nmatches_posttrim(new->substring2);
    new->nmatches_posttrim += nindels; /* for use in goodness_cmp procedures */
    new->nmatches_posttrim -= indel_penalty; /* for use in goodness_cmp procedures */
#endif

    new->trim_left = Substring_trim_left(substring1);
    new->trim_right = Substring_trim_right(substring2);
    new->trim_left_splicep = Substring_trim_left_splicep(substring1);
    new->trim_right_splicep = Substring_trim_right_splicep(substring2);

    new->penalties = indel_penalty;

    /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
    new->tally = -1L;

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->amb_nmatches_start = new->amb_nmatches_end = 0;

    new->start_ambiguous_p = new->end_ambiguous_p = false;
    new->start_ambi = new->end_ambi = (int *) NULL;
    new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
    new->start_nambi = new->end_nambi = 0;
    new->nchimera_known = 0;
    new->nchimera_novel = 0;

    new->distance = 0U;
    new->acceptor_distance = new->donor_distance = 0U;

    new->paired_usedp = false;
    new->paired_seenp = false;
    new->concordantp = false;

    new->circularpos = compute_circularpos(&new->alias,new);

    return new;
  }
}


T
Stage3end_new_deletion (int *found_score, int nindels, int indel_pos, int nmismatches1_whole, int nmismatches2_whole,
			Univcoord_T left, int genomiclength, Compress_T query_compress,
			int querylength, bool plusp, int genestrand, Chrnum_T chrnum, Univcoord_T chroffset,
			Univcoord_T chrhigh, Chrpos_T chrlength, int indel_penalty, bool sarrayp) {
  T new;
  Substring_T substring1, substring2;
  int querystart1, queryend1, querystart2, queryend2;
  Univcoord_T genomicstart, genomicend;
  Univcoord_T alignstart1, alignend1, alignstart2, alignend2;
  Univcoord_T left2;

  debug3(printf("Entered with left %lu, querylength %d, genomiclength %d, indel_pos %d\n",
		left,querylength,genomiclength,indel_pos));
#if 0
  debug3(printf("q: %s\n",query));
  debug3(printf("g: %s\n",genomicseg));
#endif

  assert(nindels > 0);

  querystart1 = 0;
  queryend1 = indel_pos;
  querystart2 = indel_pos;	/* Do not add nindels */
  queryend2 = querylength;

  if (plusp == true) {
    if ((genomicend = left + genomiclength) > chrhigh) {
      return (T) NULL;
    } else {
      genomicstart = left;

      alignstart1 = genomicstart;
      alignend1 = genomicstart + indel_pos;
      alignstart2 = alignend1 + nindels;
      alignend2 = genomicend/* + nindels*/;

      /* left1 = left; */
      left2 = left + nindels;

      debug3(printf("plusp is true.  genomicstart %lu, genomicend %lu, alignstart1 %lu, alignend1 %lu, alignstart2 %lu, alignend2 %lu, left1 %lu, left2 %lu\n",
		    genomicstart,genomicend,alignstart1,alignend1,alignstart2,alignend2,left,left2));
    }
      
  } else {
    if ((genomicstart = left + genomiclength) > chrhigh) {
      return (T) NULL;
    } else {
      genomicend = left;

      alignstart1 = genomicstart;
      alignend1 = genomicstart - indel_pos;
      alignstart2 = alignend1 - nindels;
      alignend2 = genomicend/* - nindels*/;

      /* left1 = left; */
      left2 = left + nindels;

      debug3(printf("plusp is false.  genomicstart %lu, genomicend %lu, alignstart1 %lu, alignend1 %lu, alignstart2 %lu, alignend2 %lu, left1 %lu, left2 %lu\n",
		    genomicstart,genomicend,alignstart1,alignend1,alignstart2,alignend2,left,left2));
    }
  }


  if ((substring1 = Substring_new(nmismatches1_whole,chrnum,chroffset,chrhigh,chrlength,
				  left,genomicstart,genomicend,query_compress,
				  /*start_endtype*/END,/*end_endtype*/DEL,
				  querystart1,queryend1,querylength,alignstart1,alignend1,genomiclength,
				  /*extraleft*/0,/*extraright*/0,/*exactp*/false,plusp,genestrand,
				  /*trim_left_p (previously was end1_indel_p ? false : true)*/true,
				  /*trim_right_p*/false,/*minlength*/0)) == NULL) {
    return (T) NULL;

  } else if ((substring2 = Substring_new(nmismatches2_whole,chrnum,chroffset,chrhigh,chrlength,
					 left,genomicstart,genomicend,query_compress,
					 /*start_endtype*/DEL,/*end_endtype*/END,
					 querystart2,queryend2,querylength,alignstart2,alignend2,genomiclength,
					 /*extraleft*/0,/*extraright*/0,/*exactp*/false,plusp,genestrand,
					 /*trim_left_p*/false,
					 /*trim_right_p (previously was end2_indel_p ? false : true) */true,
					 /*minlength*/0)) == NULL) {
    Substring_free(&substring1);
    return (T) NULL;
    
  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
    debug0(printf("Stage3end_new_deletion %p: left %lu, chrnum %d, nmismatches %d+%d, indel_pos %d, nindels %d\n",
		  new,left,chrnum,nmismatches1_whole,nmismatches2_whole,indel_pos,nindels));

    new->substring1 = substring1;
    new->substring2 = substring2;
    new->substring0 = (Substring_T) NULL;
    new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
    new->substringD = new->substringA = (Substring_T) NULL;

    new->pairarray = (struct Pair_T *) NULL;

#if 0
    new->deletion = (char *) CALLOC_OUT(nindels+1,sizeof(char));
    if (plusp == true) {
      strncpy(new->deletion,&(genomicseg[indel_pos]),nindels);
      new->substring_low = new->substring1;
      new->substring_high = new->substring2;
    } else {
      make_complement_buffered(new->deletion,&(genomicseg[querylength-indel_pos]),nindels);
      new->substring_low = new->substring2;
      new->substring_high = new->substring1;
    }
#else
    /* Initialize so Substring_free will not try to free */
    new->deletion = (char *) NULL;
    new->indel_pos = indel_pos;
    if (plusp == true) {
      new->substring_low = new->substring1;
      new->substring_high = new->substring2;
      new->indel_low = indel_pos;
    } else {
      new->substring_low = new->substring2;
      new->substring_high = new->substring1;
      new->indel_low = querylength - indel_pos;
    }
#endif

    new->querylength = querylength;
    new->querylength_adj = querylength + nindels;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    if (genomicstart < genomicend) {
      new->low = genomicstart;
      new->high = genomicend;
    } else {
      new->low = genomicend;
      new->high = genomicstart;
    }
    new->genomiclength = new->high - new->low;
    new->guided_insertlength = 0U;

    new->hittype = DELETION;
    new->genestrand = genestrand;
    new->sarrayp = sarrayp;
    new->improved_by_gmap_p = false;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
    new->chrlength = chrlength;
    new->plusp = plusp;
    new->sensedir = new->sensedir_nonamb = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring1) + Substring_mapq_loglik(substring2);
    new->mapq_score = 0;
    new->absmq_score = 0;
#endif

    new->nindels = nindels;
    new->nmismatches_whole = nmismatches1_whole + nmismatches2_whole;
    new->ntscore = indel_penalty + nmismatches1_whole + nmismatches2_whole;
    new->score = new->ntscore;

    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(new->substring1) + Substring_nmismatches_bothdiff(new->substring2);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1) + Substring_nmismatches_refdiff(new->substring2); */

#if 0
    /* This method is correct for SNP-tolerant alignment */
    new->nmatches = Substring_match_length(new->substring1) + Substring_match_length(new->substring2) - new->total_nmismatches;
#else
    /* This method is now correct for SNP-tolerant alignment */
    new->nmatches = Substring_nmatches(new->substring1) + Substring_nmatches(new->substring2);
    new->nmatches_posttrim = Substring_nmatches_posttrim(new->substring1) + Substring_nmatches_posttrim(new->substring2);
    new->nmatches_posttrim -= indel_penalty; /* for use in goodness_cmp procedures */
#endif

    new->trim_left = Substring_trim_left(substring1);
    new->trim_right = Substring_trim_right(substring2);
    new->trim_left_splicep = Substring_trim_left_splicep(substring1);
    new->trim_right_splicep = Substring_trim_right_splicep(substring2);

    new->penalties = indel_penalty;

    /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
    new->tally = -1L;

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->amb_nmatches_start = new->amb_nmatches_end = 0;

    new->start_ambiguous_p = new->end_ambiguous_p = false;
    new->start_ambi = new->end_ambi = (int *) NULL;
    new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
    new->start_nambi = new->end_nambi = 0;
    new->nchimera_known = 0;
    new->nchimera_novel = 0;

    new->distance = 0U;
    new->acceptor_distance = new->donor_distance = 0U;

    new->paired_usedp = false;
    new->paired_seenp = false;
    new->concordantp = false;

    new->circularpos = compute_circularpos(&new->alias,new);

    return new;
  }
}


/* Never returns NULL */
T
Stage3end_new_splice (int *found_score, int nmismatches_donor, int nmismatches_acceptor,
		      Substring_T donor, Substring_T acceptor, Chrpos_T distance,
		      bool shortdistancep, int splicing_penalty, int querylength,
		      int amb_nmatches, Intlist_T ambi_left, Intlist_T ambi_right,
		      Intlist_T amb_nmismatches_left, Intlist_T amb_nmismatches_right,
		      bool copy_donor_p, bool copy_acceptor_p, bool first_read_p, int sensedir,
		      bool sarrayp) {
  T new;
  int ignore;
  Substring_T substring_for_concordance; /* always the inner substring */
  
  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_splice %p with sensedir %d, donor substring %p and acceptor substring %p\n",
		new,sensedir,donor,acceptor));

  new->deletion = (char *) NULL;
  new->querylength_adj = new->querylength = querylength;

  if (donor == NULL) {
#ifdef DONORIS1
    new->substring1 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
    new->substring2 = (Substring_T) NULL;
    new->substring_donor = (Substring_T) NULL;
    new->substring_acceptor = new->substring1;
#else
    new->substring1 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
    new->substring2 = (Substring_T) NULL;
    new->substring_donor = (Substring_T) NULL;
    new->substring_acceptor = new->substring1;
#endif
    
  } else if (acceptor == NULL) {
#ifdef DONORIS1
    new->substring1 = copy_donor_p ? Substring_copy(donor) : donor;
    new->substring2 = (Substring_T) NULL;
    new->substring_donor = new->substring1;
    new->substring_acceptor = (Substring_T) NULL;
#else
    new->substring1 = copy_donor_p ? Substring_copy(donor) : donor;
    new->substring2 = (Substring_T) NULL;
    new->substring_donor = new->substring1;
    new->substring_acceptor = (Substring_T) NULL;
#endif

  } else {
#ifdef DONORIS1
    new->substring1 = copy_donor_p ? Substring_copy(donor) : donor;
    new->substring2 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
    new->substring_donor = new->substring1;
    new->substring_acceptor = new->substring2;
#else
    if (sensedir == SENSE_FORWARD) {
      new->substring1 = copy_donor_p ? Substring_copy(donor) : donor;
      new->substring2 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
      new->substring_donor = new->substring1;
      new->substring_acceptor = new->substring2;
    } else if (sensedir == SENSE_ANTI) {
      new->substring1 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
      new->substring2 = copy_donor_p ? Substring_copy(donor) : donor;
      new->substring_donor = new->substring2;
      new->substring_acceptor = new->substring1;
    } else {
      abort();
    }
#endif
  }
  new->substring0 = (Substring_T) NULL;
  new->substringD = new->substringA = (Substring_T) NULL;
  new->nindels = 0;
  new->indel_pos = 0;
  new->indel_low = 0;

  new->pairarray = (struct Pair_T *) NULL;

  new->sarrayp = sarrayp;
  new->improved_by_gmap_p = false;
  if (donor == NULL) {
    new->hittype = HALFSPLICE_ACCEPTOR;
    new->genestrand = Substring_genestrand(acceptor);
    new->chrnum = Substring_chrnum(acceptor);
    new->chroffset = Substring_chroffset(acceptor);
    new->chrhigh = Substring_chrhigh(acceptor);
    new->chrlength = Substring_chrlength(acceptor);
    new->plusp = Substring_plusp(acceptor);

  } else if (acceptor == NULL) {
    new->hittype = HALFSPLICE_DONOR;
    new->genestrand = Substring_genestrand(donor);
    new->chrnum = Substring_chrnum(donor);
    new->chroffset = Substring_chroffset(donor);
    new->chrhigh = Substring_chrhigh(donor);
    new->chrlength = Substring_chrlength(donor);
    new->plusp = Substring_plusp(donor);

  } else if (shortdistancep == true) {
    new->hittype = SPLICE;
    new->genestrand = Substring_genestrand(donor);
    new->chrnum = Substring_chrnum(donor);
    new->chroffset = Substring_chroffset(donor);
    new->chrhigh = Substring_chrhigh(donor);
    new->chrlength = Substring_chrlength(donor);

#if 0
    if (sensedir == SENSE_FORWARD) {
      if (first_read_p) {
	new->plusp = Substring_plusp(acceptor);
      } else {
	new->plusp = Substring_plusp(donor);
      }

    } else if (sensedir == SENSE_ANTI) {
      if (first_read_p) {
	new->plusp = Substring_plusp(donor);
      } else {
	new->plusp = Substring_plusp(acceptor);
      }

    } else {
      /* No good selection for SENSE_NULL */
      new->plusp = Substring_plusp(donor);
    }
#else
    assert(Substring_plusp(donor) == Substring_plusp(acceptor));
    new->plusp = Substring_plusp(donor);
#endif


#if 0
  } else if (merge_samechr_p == false) {
    new->hittype = DISTANT_SPLICE;
    new->sarrayp = sarrayp;
    new->improved_by_gmap_p = false;
    new->chrnum = 0;
    new->chroffset = 0;
    new->chrhigh = 0;
    new->chrlength = 0;
#endif

  } else {
    new->sarrayp = sarrayp;
    new->improved_by_gmap_p = false;
    if (Substring_chrnum(donor) == Substring_chrnum(acceptor)) {
      new->hittype = SAMECHR_SPLICE;
      new->genestrand = Substring_genestrand(donor);
      new->chrnum = Substring_chrnum(donor);
      new->chroffset = Substring_chroffset(donor);
      new->chrhigh = Substring_chrhigh(donor);
      new->chrlength = Substring_chrlength(donor);
      new->plusp = Substring_plusp(donor); /* default value, used if merge_samechr_p is true */

    } else {
      new->hittype = TRANSLOC_SPLICE;
      new->genestrand = 0;
      new->chrnum = 0;
      new->chroffset = 0;
      new->chrhigh = 0;
      new->chrlength = 0;
    }
    
    /* new->plusp assigned below */

#if 0
    if (Substring_plusp(donor) == Substring_plusp(acceptor)) {
      new->plusp = Substring_plusp(donor);
    } else {
      /* Not sure what to do here.  Probably need to have substring->dir rather than substring->plusp. */
      /* Look at ss.samechr for an example.  plusp true => pair_inversion, plusp false => pair_scramble. */
      new->plusp = true;
    }
#endif
  }

  /* printf("Making splice with shortdistancep = %d, donor chrnum %d, and acceptor chrnum %d => chrnum %d\n",
     shortdistancep,Substring_chrnum(donor),Substring_chrnum(acceptor),new->chrnum); */

  if (donor == NULL) {
    new->genomicstart = Substring_genomicstart(acceptor);
    new->genomicend = Substring_genomicend(acceptor);

  } else if (acceptor == NULL) {
    new->genomicstart = Substring_genomicstart(donor);
    new->genomicend = Substring_genomicend(donor);

  } else if (sensedir == SENSE_FORWARD) {
    new->genomicstart = Substring_genomicstart(donor);
    new->genomicend = Substring_genomicend(acceptor);

  } else if (sensedir == SENSE_ANTI) {
    new->genomicstart = Substring_genomicstart(acceptor);
    new->genomicend = Substring_genomicend(donor);

  } else {
    abort();
  }
  debug0(printf("  hittype is %s, genomicpos %u..%u\n",
		hittype_string(new->hittype),new->genomicstart,new->genomicend));

  if (new->genomicstart < new->genomicend) {
    new->low = new->genomicstart;
    new->high = new->genomicend;

    if (ambi_left != NULL) {
      new->amb_nmatches_start = amb_nmatches;
      new->amb_nmatches_end = 0;
    } else if (ambi_right != NULL) {
      new->amb_nmatches_end = amb_nmatches;
      new->amb_nmatches_start = 0;
    } else {
      new->amb_nmatches_start = new->amb_nmatches_end = 0;
    }

    new->start_ambiguous_p = (ambi_left != NULL) ? true : false;
    new->end_ambiguous_p = (ambi_right != NULL) ? true : false;
    new->start_ambi = Intlist_to_array(&new->start_nambi,ambi_left);
    new->end_ambi = Intlist_to_array(&new->end_nambi,ambi_right);
    new->start_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_left);
    new->end_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_right);

  } else {
    new->low = new->genomicend;
    new->high = new->genomicstart;

    if (ambi_right != NULL) {
      new->amb_nmatches_start = amb_nmatches;
      new->amb_nmatches_end = 0;
    } else if (ambi_left != NULL) {
      new->amb_nmatches_end = amb_nmatches;
      new->amb_nmatches_start = 0;
    } else {
      new->amb_nmatches_start = new->amb_nmatches_end = 0;
    }

    new->start_ambiguous_p = (ambi_right != NULL) ? true : false;
    new->end_ambiguous_p = (ambi_left != NULL) ? true : false;
    new->start_ambi = Intlist_to_array(&new->start_nambi,ambi_right);
    new->end_ambi = Intlist_to_array(&new->end_nambi,ambi_left);
    new->start_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_right);
    new->end_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_left);
  }
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;


  new->nchimera_known = Substring_nchimera_known(donor) + Substring_nchimera_known(acceptor);
  new->nchimera_novel = Substring_nchimera_novel(donor) + Substring_nchimera_novel(acceptor);
  if (new->start_ambiguous_p == true && favor_ambiguous_p == true) {
    new->nchimera_known++;
    /* new->nchimera_novel--; */
  }
  if (new->end_ambiguous_p == true && favor_ambiguous_p == true) {
    new->nchimera_known++;
    /* new->nchimera_novel--; */
  }

  if (new->chrnum == 0 || (donor != NULL && acceptor != NULL && shortdistancep == false && merge_samechr_p == false)) {
    /* Always want the original query end */
    if (first_read_p == true) {
      if (invert_first_p == false) {
	if (Substring_queryend(acceptor) > Substring_queryend(donor)) {
	  substring_for_concordance = acceptor;
	  new->substring_low = new->substring_high = new->substring_acceptor;
	  new->effective_chrnum = Substring_chrnum(acceptor);
	  new->other_chrnum = Substring_chrnum(donor);
	} else {
	  substring_for_concordance = donor;
	  new->substring_low = new->substring_high = new->substring_donor;
	  new->effective_chrnum = Substring_chrnum(donor);
	  new->other_chrnum = Substring_chrnum(acceptor);
	}
      } else {
	if (Substring_querystart(acceptor) < Substring_querystart(donor)) {
	  substring_for_concordance = acceptor;
	  new->substring_low = new->substring_high = new->substring_acceptor;
	  new->effective_chrnum = Substring_chrnum(acceptor);
	  new->other_chrnum = Substring_chrnum(donor);
	} else {
	  substring_for_concordance = donor;
	  new->substring_low = new->substring_high = new->substring_donor;
	  new->effective_chrnum = Substring_chrnum(donor);
	  new->other_chrnum = Substring_chrnum(acceptor);
	}
      }

    } else {
      if (invert_second_p == false) {
	if (Substring_queryend(acceptor) > Substring_queryend(donor)) {
	  substring_for_concordance = acceptor;
	  new->substring_low = new->substring_high = new->substring_acceptor;
	  new->effective_chrnum = Substring_chrnum(acceptor);
	  new->other_chrnum = Substring_chrnum(donor);
	} else {
	  substring_for_concordance = donor;
	  new->substring_low = new->substring_high = new->substring_donor;
	  new->effective_chrnum = Substring_chrnum(donor);
	  new->other_chrnum = Substring_chrnum(acceptor);
	}
      } else {
	if (Substring_querystart(acceptor) < Substring_querystart(donor)) {
	  substring_for_concordance = acceptor;
	  new->substring_low = new->substring_high = new->substring_acceptor;
	  new->effective_chrnum = Substring_chrnum(acceptor);
	  new->other_chrnum = Substring_chrnum(donor);
	} else {
	  substring_for_concordance = donor;
	  new->substring_low = new->substring_high = new->substring_donor;
	  new->effective_chrnum = Substring_chrnum(donor);
	  new->other_chrnum = Substring_chrnum(acceptor);
	}
      }
    }

    /* Redefine based on inner substring */
    new->genomicstart = Substring_genomicstart(substring_for_concordance);
    new->genomicend = Substring_genomicend(substring_for_concordance);
    new->plusp = Substring_plusp(substring_for_concordance);
    
  } else {
    new->effective_chrnum = new->chrnum;
    new->other_chrnum = 0;

    if (donor == NULL) {
      new->substring_low = new->substring_high = new->substring1;
    } else if (acceptor == NULL) {
      new->substring_low = new->substring_high = new->substring1;
    } else if (sensedir == SENSE_FORWARD) {
#ifdef DONORIS1
      if (new->plusp == true) {
	new->substring_low = new->substring1; /* donor */
	new->substring_high = new->substring2; /* acceptor */
      } else {
	new->substring_low = new->substring2; /* acceptor */
	new->substring_high = new->substring1; /* donor */
      }
#else
      if (new->plusp == true) {
	new->substring_low = new->substring_donor;
	new->substring_high = new->substring_acceptor;
      } else {
	new->substring_low = new->substring_acceptor;
	new->substring_high = new->substring_donor;
      }
#endif

    } else if (sensedir == SENSE_ANTI) {
#ifdef DONORIS1
      if (new->plusp == true) {
	new->substring_low = new->substring2; /* acceptor */
	new->substring_high = new->substring1; /* donor */
      } else {
	new->substring_low = new->substring1; /* donor */
	new->substring_high = new->substring2; /* acceptor */
      }
#else
      if (new->plusp == true) {
	new->substring_low = new->substring_acceptor;
	new->substring_high = new->substring_donor;
      } else {
	new->substring_low = new->substring_donor;
	new->substring_high = new->substring_acceptor;
      }
#endif

    } else {
      abort();
    }
  }

  new->nmismatches_whole = nmismatches_donor + nmismatches_acceptor;
  new->score = new->ntscore = splicing_penalty + new->nmismatches_whole;
  if (sensedir == SENSE_FORWARD) {
    new->score += antistranded_penalty;
  }

  if (donor == NULL) {
    /* new->mapq_loglik = Substring_mapq_loglik(acceptor); */
    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(acceptor) + nmismatches_donor;
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(acceptor) + nmismatches_donor; */
    new->nmatches = Substring_nmatches(acceptor);
    new->nmatches_posttrim = Substring_nmatches_posttrim(acceptor);
    if (favor_ambiguous_p == true) {
      new->nmatches += amb_nmatches;
    }
    new->sensedir_nonamb = SENSE_NULL;	/* Ignore sense based on ambiguous end */
    debug0(printf("New splice has acceptor %d + amb %d matches, sensedir nonamb %d\n",
		  Substring_nmatches(acceptor),amb_nmatches,new->sensedir_nonamb));
  } else if (acceptor == NULL) {
    /* new->mapq_loglik = Substring_mapq_loglik(donor); */
    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(donor) + nmismatches_acceptor;
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(donor) + nmismatches_acceptor; */
    new->nmatches = Substring_nmatches(donor);
    new->nmatches_posttrim = Substring_nmatches_posttrim(donor);
    if (favor_ambiguous_p == true) {
      new->nmatches += amb_nmatches;
    }
    new->sensedir_nonamb = SENSE_NULL;	/* Ignore sense based on ambiguous end */
    debug0(printf("New splice has donor %d + amb %d matches, sensedir nonamb %d\n",
		  Substring_nmatches(donor),amb_nmatches,new->sensedir_nonamb));
  } else {
    /* new->mapq_loglik = Substring_mapq_loglik(donor) + Substring_mapq_loglik(acceptor); */
    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(donor) + Substring_nmismatches_bothdiff(acceptor);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(donor) + Substring_nmismatches_refdiff(acceptor); */
    new->nmatches = Substring_nmatches(donor) + Substring_nmatches(acceptor);
    new->nmatches_posttrim = Substring_nmatches_posttrim(donor) + Substring_nmatches_posttrim(acceptor);
    new->sensedir_nonamb = sensedir;
    debug0(printf("New splice has donor %d + acceptor %d matches, sensedir nonamb %d\n",
		  Substring_nmatches(donor),Substring_nmatches(acceptor),new->sensedir_nonamb));
  }
  new->sensedir = sensedir;

  if (new->substring0 != NULL) {
    new->trim_left = Substring_trim_left(new->substring0);
    new->trim_left_splicep = Substring_trim_left_splicep(new->substring0);
  } else {
    new->trim_left = Substring_trim_left(new->substring1);
    new->trim_left_splicep = Substring_trim_left_splicep(new->substring1);
  }

  if (new->substring2 != NULL) {
    new->trim_right = Substring_trim_right(new->substring2);
    new->trim_right_splicep = Substring_trim_right_splicep(new->substring2);
  } else {
    new->trim_right = Substring_trim_right(new->substring1);
    new->trim_right_splicep = Substring_trim_right_splicep(new->substring1);
  }
  
  new->penalties = splicing_penalty;


  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

#if 0
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->distance = distance;
  new->acceptor_distance = new->donor_distance = 0U;

  new->paired_usedp = false;
  new->paired_seenp = false;
  new->concordantp = false;

  new->circularpos = compute_circularpos(&new->alias,new);

  assert(new->substring1 != NULL);

  debug0(printf("Returning new splice %p at genomic %u..%u\n",new,new->genomicstart,new->genomicend));
  return new;
}



/* Never returns NULL.  Never copies substrings.  Always shortdistance. */
T
Stage3end_new_shortexon (int *found_score, Substring_T donor, Substring_T acceptor, Substring_T shortexon,
			 Chrpos_T acceptor_distance, Chrpos_T donor_distance,
			 int amb_nmatches_donor, int amb_nmatches_acceptor,
			 Intlist_T ambi_left, Intlist_T ambi_right,
			 Intlist_T amb_nmismatches_left, Intlist_T amb_nmismatches_right,
			 bool copy_donor_p, bool copy_acceptor_p, bool copy_shortexon_p,
			 int splicing_penalty, int querylength, int sensedir, bool sarrayp) {
  T new;
  int ignore;
  
  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_shortexon %p\n",new));

  new->deletion = (char *) NULL;
  new->querylength_adj = new->querylength = querylength;

  new->genestrand = Substring_genestrand(shortexon);
  new->sarrayp = sarrayp;
  new->improved_by_gmap_p = false;
  if (donor == NULL && acceptor == NULL) {
    new->hittype = ONE_THIRD_SHORTEXON;
    new->substring1 = copy_shortexon_p ? Substring_copy(shortexon) : shortexon;
    new->substring2 = (Substring_T) NULL;
    new->substring0 = (Substring_T) NULL;
    new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
    new->substringD = new->substringA = (Substring_T) NULL;
    new->sensedir_nonamb = SENSE_NULL;	/* Ignore sensedir based on double ambiguous ends */

  } else {
    if (donor == NULL) {
      new->hittype = TWO_THIRDS_SHORTEXON;
    } else if (acceptor == NULL) {
      new->hittype = TWO_THIRDS_SHORTEXON;
    } else {
      new->hittype = SHORTEXON;
    }
    new->substring1 = copy_shortexon_p ? Substring_copy(shortexon) : shortexon;
    if (sensedir == SENSE_FORWARD) {
      new->substringD = new->substring0 = copy_donor_p ? Substring_copy(donor) : donor;
      new->substringA = new->substring2 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
    } else if (sensedir == SENSE_ANTI) {
      new->substringA = new->substring0 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
      new->substringD = new->substring2 = copy_donor_p ? Substring_copy(donor) : donor;
    } else {
      abort();
    }
    new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
    new->sensedir_nonamb = sensedir;
  }
  new->sensedir = sensedir;

  new->pairarray = (struct Pair_T *) NULL;

  new->nindels = 0;
  new->indel_pos = 0;
  new->indel_low = 0;

  new->chrnum = Substring_chrnum(shortexon);
  new->chroffset = Substring_chroffset(shortexon);
  new->chrhigh = Substring_chrhigh(shortexon);
  new->chrlength = Substring_chrlength(shortexon);
  new->plusp = Substring_plusp(shortexon);

  /* printf("Making splice with shortdistancep = %d, donor chrnum %d, and acceptor chrnum %d => chrnum %d\n",
     shortdistancep,Substring_chrnum(donor),Substring_chrnum(acceptor),new->chrnum); */

  if (sensedir == SENSE_FORWARD) {
    new->genomicstart = (donor != NULL ? Substring_genomicstart(donor) : Substring_genomicstart(shortexon));
    new->genomicend = (acceptor != NULL ? Substring_genomicend(acceptor) : Substring_genomicend(shortexon));
  } else if (sensedir == SENSE_ANTI) {
    new->genomicstart = (acceptor != NULL ? Substring_genomicstart(acceptor) : Substring_genomicstart(shortexon));
    new->genomicend = (donor != NULL ? Substring_genomicend(donor) : Substring_genomicend(shortexon));
  } else {
    abort();
  }

  if (new->genomicstart < new->genomicend) {
    new->low = new->genomicstart;
    new->high = new->genomicend;

    if (ambi_left != NULL) {
      new->amb_nmatches_start = (sensedir == SENSE_FORWARD) ? amb_nmatches_donor : amb_nmatches_acceptor;
    } else {
      new->amb_nmatches_start = 0;
    }

    if (ambi_right != NULL) {
      new->amb_nmatches_end = (sensedir == SENSE_FORWARD) ? amb_nmatches_acceptor : amb_nmatches_donor;
    } else {
      new->amb_nmatches_end = 0;
    }

    new->start_ambiguous_p = (ambi_left != NULL) ? true : false;
    new->end_ambiguous_p = (ambi_right != NULL) ? true : false;
    new->start_ambi = Intlist_to_array(&new->start_nambi,ambi_left);
    new->end_ambi = Intlist_to_array(&new->end_nambi,ambi_right);
    new->start_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_left);
    new->end_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_right);

  } else {
    new->low = new->genomicend;
    new->high = new->genomicstart;

    if (ambi_right != NULL) {
      new->amb_nmatches_start = (sensedir == SENSE_FORWARD) ? amb_nmatches_donor : amb_nmatches_acceptor;
    } else {
      new->amb_nmatches_start = 0;
    }

    if (ambi_left != NULL) {
      new->amb_nmatches_end = (sensedir == SENSE_FORWARD) ? amb_nmatches_acceptor : amb_nmatches_donor;
    } else {
      new->amb_nmatches_end = 0;
    }

    new->start_ambiguous_p = (ambi_right != NULL) ? true : false;
    new->end_ambiguous_p = (ambi_left != NULL) ? true : false;
    new->start_ambi = Intlist_to_array(&new->start_nambi,ambi_right);
    new->end_ambi = Intlist_to_array(&new->end_nambi,ambi_left);
    new->start_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_right);
    new->end_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_left);
  }
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;


  new->nchimera_known = Substring_nchimera_known(shortexon) + Substring_nchimera_known(donor) + Substring_nchimera_known(acceptor);
  new->nchimera_novel = Substring_nchimera_novel(shortexon) + Substring_nchimera_novel(donor) + Substring_nchimera_novel(acceptor);
  if (new->start_ambiguous_p == true && favor_ambiguous_p == true) {
    new->nchimera_known++;
    /* new->nchimera_novel--; */
  }
  if (new->end_ambiguous_p == true && favor_ambiguous_p == true) {
    new->nchimera_known++;
    /* new->nchimera_novel--; */
  }


  new->effective_chrnum = new->chrnum;
  new->other_chrnum = 0;

  /* Currently not allowing translocations on shortexons */
  /* substring_for_concordance = (Substring_T) NULL; */

#if 0
  if (sensedir == SENSE_FORWARD) {
    if (new->plusp == true) {
      new->substring_low = (new->substringD != NULL ? new->substringD : new->substring1); /* donor */
      new->substring_high = (new->substringA != NULL ? new->substringA : new->substring1); /* acceptor */
    } else {
      new->substring_low = (new->substringA != NULL ? new->substringA : new->substring1); /* acceptor */
      new->substring_high = (new->substringD != NULL ? new->substringD : new->substring1); /* donor */
    }

  } else if (sensedir == SENSE_ANTI) {
    if (new->plusp == true) {
      new->substring_low = (new->substringA != NULL ? new->substringA : new->substring1); /* acceptor */
      new->substring_high = (new->substringD != NULL ? new->substringD : new->substring1); /* donor */
    } else {
      new->substring_low = (new->substringD != NULL ? new->substringD : new->substring1); /* donor */
      new->substring_high = (new->substringA != NULL ? new->substringA : new->substring1); /* acceptor */
    }

  } else {
    abort();
  }
#else
  if (new->plusp == true) {
    new->substring_low = (new->substring0 != NULL ? new->substring0 : new->substring1);
    new->substring_high = (new->substring2 != NULL ? new->substring2 : new->substring1);
  } else {
    new->substring_low = (new->substring2 != NULL ? new->substring2 : new->substring1);
    new->substring_high = (new->substring0 != NULL ? new->substring0 : new->substring1);
  }
#endif


#if 0
  new->mapq_loglik = Substring_mapq_loglik(donor) + Substring_mapq_loglik(acceptor) + Substring_mapq_loglik(shortexon);
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  new->nmismatches_whole = Substring_nmismatches_whole(shortexon);
  if (donor != NULL) {
    new->nmismatches_whole += Substring_nmismatches_whole(donor);
  }
  if (acceptor != NULL) {
    new->nmismatches_whole += Substring_nmismatches_whole(acceptor);
  }
  new->ntscore = splicing_penalty + splicing_penalty + new->nmismatches_whole;

  new->score = new->ntscore;
  if (sensedir == SENSE_FORWARD) {
    new->score += antistranded_penalty;
  }

  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(shortexon);
  if (donor != NULL) {
    new->nmismatches_bothdiff += Substring_nmismatches_bothdiff(donor);
  }
  if (acceptor != NULL) {
    new->nmismatches_bothdiff += Substring_nmismatches_bothdiff(acceptor);
  }
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(donor) + Substring_nmismatches_refdiff(acceptor) + Substring_nmismatches_refdiff(shortexon); */

  new->amb_nmatches_donor = amb_nmatches_donor;
  new->amb_nmatches_acceptor = amb_nmatches_acceptor;

  new->nmatches = Substring_nmatches(shortexon);
  new->nmatches_posttrim = Substring_nmatches_posttrim(shortexon);
  if (donor == NULL) {
    if (favor_ambiguous_p == true) {
      new->nmatches += amb_nmatches_donor;
    }
  } else {
    assert(amb_nmatches_donor == 0);
    new->nmatches += Substring_nmatches(donor);
  }
  if (acceptor == NULL) {
    if (favor_ambiguous_p == true) {
      new->nmatches += amb_nmatches_acceptor;
    }
  } else {
    assert(amb_nmatches_acceptor == 0);
    new->nmatches += Substring_nmatches(acceptor);
  }

  if (new->substring0 != NULL) {
    new->trim_left = Substring_trim_left(new->substring0);
    new->trim_left_splicep = Substring_trim_left_splicep(new->substring0);
  } else {
    new->trim_left = Substring_trim_left(new->substring1);
    new->trim_left_splicep = Substring_trim_left_splicep(new->substring1);
  }

  if (new->substring2 != NULL) {
    new->trim_right = Substring_trim_right(new->substring2);
    new->trim_right_splicep = Substring_trim_right_splicep(new->substring2);
  } else {
    new->trim_right = Substring_trim_right(new->substring1);
    new->trim_right_splicep = Substring_trim_right_splicep(new->substring1);
  }

  new->penalties = splicing_penalty + splicing_penalty;

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->distance = acceptor_distance + donor_distance;
  new->acceptor_distance = acceptor_distance;
  new->donor_distance = donor_distance;

  new->paired_usedp = false;
  new->paired_seenp = false;
  new->concordantp = false;

  new->circularpos = compute_circularpos(&new->alias,new);
    
  return new;
}


T
Stage3end_new_terminal (int querystart, int queryend, Univcoord_T left, Compress_T query_compress,
			int querylength, bool plusp, int genestrand,
			Endtype_T start_endtype, Endtype_T end_endtype,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			int max_mismatches_allowed, bool sarrayp) {
  T new;
  Substring_T substring;
  Univcoord_T genomicstart, genomicend, alignstart, alignend, alignstart_trim, alignend_trim;
  int nmismatches_whole, minlength;
  bool trim_left_p, trim_right_p;

  debug0(printf("\nStage3end_new_terminal possible: endtypes %s and %s, left %lu, querystart %d, queryend %d\n",
		Endtype_string(start_endtype),Endtype_string(end_endtype),left,querystart,queryend));

  if (plusp == true) {
    if ((genomicend = left + querylength) > chrhigh) {
      return (T) NULL;
    } else {
      genomicstart = left;

      alignstart = genomicstart + querystart;
      alignend = genomicstart + queryend;
    }

  } else {
    if ((genomicstart = left + querylength) > chrhigh) {
      return (T) NULL;
    } else {
      genomicend = left;
    
      alignstart = genomicstart - querystart;
      alignend = genomicstart - queryend;
    }
  }

  if (start_endtype == TERM) {
    trim_left_p = true;
  } else {
    trim_left_p = false;
  }

  if (end_endtype == TERM) {
    trim_right_p = true;
  } else {
    trim_right_p = false;
  }

  /* Note: Changing querylength/3 to querylength/2 loses about 1% of concordant reads */
  minlength = querylength/3;
  if (minlength > TERMINAL_COMPUTE_MINLENGTH) {
    minlength = TERMINAL_COMPUTE_MINLENGTH;
  }

  if ((substring = Substring_new(/*nmismatches_whole*/0,chrnum,chroffset,chrhigh,chrlength,
				 left,genomicstart,genomicend,query_compress,
				 start_endtype,end_endtype,querystart,queryend,querylength,
				 alignstart,alignend,/*genomiclength*/querylength,
				 /*extraleft*/0,/*extraright*/0,/*exactp*/false,plusp,genestrand,
				 trim_left_p,trim_right_p,minlength)) == NULL) {
    debug0(printf("returning NULL\n"));
    return (T) NULL;
    
#if 0
  } else if (start_endtype == TERM) {
    if (Substring_trim_left(substring) == 0) {
      /* Not really a terminal */
      Substring_free(&substring);
      return (T) NULL;
    }

  } else if (end_endtype == TERM) {
    if (Substring_trim_right(substring) == 0) {
      /* Not really a terminal */
      Substring_free(&substring);
      return (T) NULL;
    }
#endif
  }

  /* Re-compute nmismatches_whole and nmatches for terminal alignments */
  alignstart_trim = Substring_alignstart_trim(substring);
  alignend_trim = Substring_alignend_trim(substring);

  debug0(printf("alignstart_trim = %lu, alignend_trim = %lu\n",alignstart_trim,alignend_trim));
  if (plusp == true) {
    debug0(printf("plus: pos5 = %d, pos3 = %d\n",(int) (alignstart_trim-left),(int) (alignend_trim-left)));
    nmismatches_whole =
      Genome_count_mismatches_substring(query_compress,left,/*pos5*/alignstart_trim-left,
					/*pos3*/alignend_trim-left,/*plusp*/true,genestrand);
  } else {
    debug0(printf("minus: pos5 = %d, pos3 = %d\n",(int) (alignend_trim-left),(int) (alignstart_trim-left)));
    nmismatches_whole =
      Genome_count_mismatches_substring(query_compress,left,/*pos5*/alignend_trim-left,
					/*pos3*/alignstart_trim-left,/*plusp*/false,genestrand);
  }
  debug0(printf("Recomputing nmismatches_whole as %d\n",nmismatches_whole));

  if (nmismatches_whole > max_mismatches_allowed) {
    /* This may be dependent on the trimming algorithm, but is needed to avoid bad terminal alignments */
    Substring_free(&substring);
    debug0(printf("returning NULL\n"));
    return (T) NULL;
  } else {
    /* Code in substring.c suggests that nmismatches_bothdiff would be the same value as nmismatches_whole */
    Substring_set_nmismatches_terminal(substring,nmismatches_whole,/*nmismatches_bothdiff*/nmismatches_whole);
#if 0
    if (alignstart_trim >= alignend_trim) {
      Substring_set_endtypes(substring,/*start_endtype*/TERM,/*end_endtype*/END);
    } else {
      Substring_set_endtypes(substring,/*start_endtype*/END,/*end_endtype*/TERM);
    }
#endif
  }

  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_terminal %p: endtypes %s and %s, left %lu, genomicstart/end %lu..%lu, chrhigh %lu, chrnum %d, querystart %d, queryend %d\n",
		new,Endtype_string(start_endtype),Endtype_string(end_endtype),
		left,genomicstart,genomicend,chrhigh,chrnum,querystart,queryend));

  new->substring1 = substring;
  new->substring2 = (Substring_T) NULL;
  new->substring0 = (Substring_T) NULL;
  new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
  new->substringD = new->substringA = (Substring_T) NULL;
  new->substring_low = new->substring_high = new->substring1;

  new->pairarray = (struct Pair_T *) NULL;

  new->deletion = (char *) NULL;
  new->querylength_adj = new->querylength = querylength;
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  if (genomicstart < genomicend) {
    new->low = genomicstart;
    new->high = genomicend;
  } else {
    new->low = genomicend;
    new->high = genomicstart;
  }
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;


  new->hittype = TERMINAL;
  new->genestrand = genestrand;
  new->sarrayp = sarrayp;
  new->improved_by_gmap_p = false;

  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;
  new->plusp = plusp;

#if 0
  new->mapq_loglik = Substring_mapq_loglik(substring);
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  new->nindels = 0;
  new->indel_pos = 0;
  new->indel_low = 0;
  new->nmismatches_whole = Substring_nmismatches_whole(substring); /* This value was recomputed to include non-terminal end */
  new->ntscore = /* terminal_penalty + */ nmismatches_whole;

  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(substring);
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(substring); */

#if 0
  new->score = terminal_penalty + nmismatches_whole;
#else
  new->score = /* terminal_penalty + */ Substring_nmismatches_whole(substring);
#endif

#if 0
  new->nmatches = Substring_match_length(substring) - nmismatches;
#else
  new->nmatches = Substring_nmatches(substring);
  new->nmatches_posttrim = Substring_nmatches_posttrim(substring);
#endif

  new->trim_left = Substring_trim_left(substring);
  new->trim_right = Substring_trim_right(substring);
  new->trim_left_splicep = Substring_trim_left_splicep(substring);
  new->trim_right_splicep = Substring_trim_right_splicep(substring);

  new->penalties = 0;

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

  new->amb_nmatches_start = new->amb_nmatches_end = 0;

  new->start_ambiguous_p = new->end_ambiguous_p = false;
  new->start_ambi = new->end_ambi = (int *) NULL;
  new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
  new->start_nambi = new->end_nambi = 0;
  new->nchimera_known = 0;
  new->nchimera_novel = 0;

  new->distance = 0U;
  new->acceptor_distance = new->donor_distance = 0U;
  new->sensedir = new->sensedir_nonamb = SENSE_NULL;

  new->paired_usedp = false;
  new->paired_seenp = false;
  new->concordantp = false;

  new->circularpos = compute_circularpos(&new->alias,new);

  return new;
}


T
Stage3end_new_gmap (int nmismatches_whole, int nmatches_posttrim, int max_match_length,
		    int ambig_end_length_5, int ambig_end_length_3,
		    Splicetype_T ambig_splicetype_5, Splicetype_T ambig_splicetype_3,
		    double min_splice_prob, struct Pair_T *pairarray, int npairs,
		    int nsegments, int nintrons, int nindelbreaks,
		    Univcoord_T left, int genomiclength, bool plusp, int genestrand, int querylength,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		    int cdna_direction, int sensedir) {
  T new;
  Univcoord_T genomicstart, genomicend, genomepos;
  double prob1, prob2;

  /* In 2012-12-20, removed statements to return NULL, because GMAP alignments seem
     to be okay, at least when starting before coordinate 0 */
  /* Example (when aligned to chrM at beginning of genome) (actually aligns circularly):
GGATGAGGCAGGAATCAAAGACAGATACTGCGACATAGGGTGCTCCGGCTCCAGCGTCTCGCAATGCTATCGCGTG
ATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAG
  */
  /* However, this leads to fatal bugs later, so restored these statements */

  if (Stage3_bad_stretch_p(pairarray,npairs,/*pos5*/0,/*pos3*/querylength) == true) {
    debug0(printf("Bad GMAP: bad stretch\n"));
    return (T) NULL;

  } else if (plusp == true) {
    genomicstart = left;
    if ((genomicend = left + genomiclength) > chrhigh) {
      return (T) NULL;
    }
    if (genomicstart > genomicend) {
      /* Must have started before coordinate 0 */
      debug0(printf("plusp and genomicstart %lu > genomicend %lu => started before coordinate 0\n",
		    genomicstart,genomicend));
      return (T) NULL;
    }
  } else {
    if ((genomicstart = left + genomiclength) > chrhigh) {
      return (T) NULL;
    }
    genomicend = left;
    if (genomicend > genomicstart) {
      /* Must have started before coordinate 0 */
      debug0(printf("minusp and genomicend %lu > genomicstart %lu => started before coordinate 0\n",
		    genomicend,genomicstart));
      return (T) NULL;
    }
  }

  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_gmap %p: left %lu, genomicstart/end %lu..%lu, chrhigh %lu, chrnum %d, nmismatches %d, cdna_direction %d, sensedir %d, max_match_length %d\n",
		new,left,genomicstart,genomicend,chrhigh,chrnum,nmismatches_whole,cdna_direction,sensedir,max_match_length));
  debug0(printf("  ambig_end_length_5 %d, ambig_end_length_3 %d\n",ambig_end_length_5,ambig_end_length_3));

  new->substring1 = (Substring_T) NULL;
  new->substring2 = (Substring_T) NULL;
  new->substring0 = (Substring_T) NULL;
  new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
  new->substringD = new->substringA = (Substring_T) NULL;
  new->substring_low = new->substring_high = (Substring_T) NULL;

  new->pairarray = pairarray;
  new->npairs = npairs;
  new->nsegments = nsegments;

  new->deletion = (char *) NULL;
  new->querylength_adj = new->querylength = querylength /* - nindels */;
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  if (genomicstart < genomicend) {
    new->low = genomicstart;
    new->high = genomicend;
  } else {
    new->low = genomicend;
    new->high = genomicstart;
  }
  new->genomiclength = new->high - new->low;
  new->guided_insertlength = 0U;


  new->hittype = GMAP;
  new->genestrand = genestrand;
  new->sarrayp = false;
  new->improved_by_gmap_p = false;

  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->chrlength = chrlength;
  new->plusp = plusp;
  new->gmap_nindelbreaks = nindelbreaks;
  new->gmap_cdna_direction = cdna_direction;
  new->gmap_nintrons = nintrons;
  new->sensedir = new->sensedir_nonamb = sensedir;

#if 0
  new->mapq_loglik = Substring_mapq_loglik(substring);
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  new->nindels = 0;
  new->indel_pos = 0;
  new->indel_low = 0;
  new->nmismatches_whole = nmismatches_whole;
  new->ntscore = nmismatches_whole;

#if 1
  /* This favors the trimmed results */
  new->score = nmismatches_whole;
  new->score += localsplicing_penalty * nintrons;
  new->score += indel_penalty_middle * nindelbreaks;
  debug0(printf("gmap score = %d = %d + %d*%d + %d*%d\n",
		new->score,nmismatches_whole,localsplicing_penalty,nintrons,indel_penalty_middle,nindelbreaks));
#else
  /* This is a better way to score GMAP.  Using nmatches_pretrim puts all GMAP entries on an even level. */
  new->score = querylength - nmatches_posttrim;
  if (nindels > 0) {
    /* Account for the fact that a query insertion reduces number of possible posttrim matches */
    new->score -= nindels;
  }
  new->score += localsplicing_penalty * nintrons;
  new->score += indel_penalty_middle * nindelbreaks;
  debug0(printf("gmap score = %d = querylength %d (nindels %d, subtract if pos) - posttrim %d + %d*%d + %d*%d\n",
		new->score,querylength,nindels,nmatches_posttrim,
		localsplicing_penalty,nintrons,indel_penalty_middle,nindelbreaks));
#endif


  new->nmismatches_bothdiff = nmismatches_whole;
  /* new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(new->substring1); */
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1); */

  new->nmatches = nmatches_posttrim; /* To make addition of ambiguous lengths work, we need to use posttrim, not pretrim */
  new->nmatches_posttrim = nmatches_posttrim;
  if (favor_ambiguous_p == true) {
    new->nmatches += ambig_end_length_5 + ambig_end_length_3;
  }
  debug0(printf("  nmatches %d = posttrim %d + ambig_end_length_5 %d + ambig_end_length_3 %d\n",
		new->nmatches,nmatches_posttrim,ambig_end_length_5,ambig_end_length_3));
  new->nmatches_posttrim -= localsplicing_penalty * nintrons; /* for use in goodness_cmp procedures */
  new->nmatches_posttrim -= indel_penalty_middle * nindelbreaks; /* for use in goodness_cmp procedures */

  if (new->nmatches_posttrim < querylength/2) {
    debug0(printf("  nmatches %d < querylength %d/2, so returning NULL\n",
		  new->nmatches_posttrim,querylength));
    FREE_OUT(new);
    return NULL;
  } else if (max_match_length < gmap_min_nconsecutive) {
    debug0(printf("  max_match_length %d < %d, so returning NULL\n",max_match_length,gmap_min_nconsecutive));
    FREE_OUT(new);
    return NULL;
  }

  new->gmap_max_match_length = max_match_length;
  new->gmap_min_splice_prob = min_splice_prob;


  new->trim_left = Pair_querypos(&(pairarray[0])) - ambig_end_length_5;
  if (ambig_end_length_5 > 0) {
    new->trim_left_splicep = true;
  } else if (novelsplicingp == false) {
    new->trim_left_splicep = false;
  } else {
    genomepos = chroffset + Pair_genomepos(&(pairarray[0])) + 1U;
    if (plusp == true) {
      prob1 = Maxent_hr_acceptor_prob(genomepos,chroffset);
      prob2 = Maxent_hr_antidonor_prob(genomepos,chroffset);
      /* fprintf(stderr,"At %lu, acceptor prob %f, antidonor prob %f\n",genomepos,prob1,prob2); */
    } else {
      prob1 = Maxent_hr_donor_prob(genomepos,chroffset);
      prob2 = Maxent_hr_antiacceptor_prob(genomepos,chroffset);
      /* fprintf(stderr,"At %lu, donor prob %f, antiacceptor prob %f\n",genomepos,prob1,prob2); */
    }
    if (prob1 > 0.90 || prob2 > 0.90) {
      new->trim_left_splicep = true;
    } else {
      new->trim_left_splicep = false;
    }
  }

  new->trim_right = (querylength - 1) - Pair_querypos(&(pairarray[npairs-1])) - ambig_end_length_3;
  if (ambig_end_length_3 > 0) {
    new->trim_right_splicep = true;
  } else if (novelsplicingp == false) {
    new->trim_right_splicep = false;
  } else {
    genomepos = chroffset + Pair_genomepos(&(pairarray[npairs-1])) + 1U;
    if (plusp == true) {
      prob1 = Maxent_hr_donor_prob(genomepos,chroffset);
      prob2 = Maxent_hr_antiacceptor_prob(genomepos,chroffset);
      /* fprintf(stderr,"At %lu, donor prob %f, antiacceptor prob %f\n",genomepos,prob1,prob2); */
    } else {
      prob1 = Maxent_hr_acceptor_prob(genomepos,chroffset);
      prob2 = Maxent_hr_antidonor_prob(genomepos,chroffset);
      /* fprintf(stderr,"At %lu, acceptor prob %f, antidonor prob %f\n",genomepos,prob1,prob2); */
    }
    if (prob1 > 0.90 || prob2 > 0.90) {
      new->trim_right_splicep = true;
    } else {
      new->trim_right_splicep = false;
    }
  }

  /* new->penalties not used anyway for GMAP alignments */
  new->penalties = localsplicing_penalty * nintrons + indel_penalty_middle * nindelbreaks;
  /* new->penalties += ambig_end_length_5/ambig_end_interval; */
  /* new->penalties += ambig_end_length_3/ambig_end_interval; */

  /* new->gene_overlap = NO_KNOWN_GENE; -- initialized later when resolving multimappers */
  new->tally = -1L;

  if ((new->amb_nmatches_start = ambig_end_length_5) == 0) {
    new->gmap_start_endtype = END;
  } else if (ambig_splicetype_5 == DONOR || ambig_splicetype_5 == ANTIDONOR) {
    new->gmap_start_endtype = AMB_DON;
  } else if (ambig_splicetype_5 == ACCEPTOR || ambig_splicetype_5 == ANTIACCEPTOR) {
    new->gmap_start_endtype = AMB_ACC;
  } else {
    fprintf(stderr,"Do not recognize splicetype %d for ambig_end_length_5 %d\n",
	    ambig_splicetype_5,ambig_end_length_5);
    abort();
  }
    
  if ((new->amb_nmatches_end = ambig_end_length_3) == 0) {
    new->gmap_end_endtype = END;
  } else if (ambig_splicetype_3 == DONOR || ambig_splicetype_3 == ANTIDONOR) {
    new->gmap_end_endtype = AMB_DON;
  } else if (ambig_splicetype_3 == ACCEPTOR || ambig_splicetype_3 == ANTIACCEPTOR) {
    new->gmap_end_endtype = AMB_ACC;
  } else {
    fprintf(stderr,"Do not recognize splicetype %d for ambig_end_length_3 %d\n",
	    ambig_splicetype_3,ambig_end_length_3);
    abort();
  }

  new->start_ambiguous_p = new->end_ambiguous_p = false;
  new->start_ambi = new->end_ambi = (int *) NULL;
  new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
  new->start_nambi = new->end_nambi = 0;
  new->nchimera_known = 0;
  new->nchimera_novel = 0;

  new->distance = 0U;
  new->acceptor_distance = new->donor_distance = 0U;

  new->paired_usedp = false;
  new->paired_seenp = false;
  new->concordantp = false;

  new->circularpos = compute_circularpos(&new->alias,new);

  return new;
}


static int
Stage3end_output_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else if (x->mapq_loglik > y->mapq_loglik) {
    return -1;
  } else if (y->mapq_loglik > x->mapq_loglik) {
    return +1;
  } else if (x->guided_insertlength > 0 && y->guided_insertlength == 0) {
    return -1;
  } else if (y->guided_insertlength > 0 && x->guided_insertlength == 0) {
    return +1;
  } else if (x->guided_insertlength < y->guided_insertlength) {
    return -1;
  } else if (y->guided_insertlength < x->guided_insertlength) {
    return +1;
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
#if 0
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (y->hittype < x->hittype) {
    return +1;
#endif
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (y->genomicstart < x->genomicstart) {
    return +1;
  } else if (x->plusp == true && y->plusp == false) {
    return -1;
  } else if (x->plusp == false && y->plusp == true) {
    return +1;
  } else {
    return 0;
  }
}


static int
Stage3pair_output_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;

#ifdef USE_BINGO
  if (x->absdifflength_bingo_p == true && y->absdifflength_bingo_p == false) {
    return -1;
  } else if (y->absdifflength_bingo_p == true && x->absdifflength_bingo_p == false) {
    return +1;
  }
#endif

  if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else if (x->mapq_loglik > y->mapq_loglik) {
    return -1;
  } else if (y->mapq_loglik > x->mapq_loglik) {
    return +1;
  } else if (x->insertlength > 0 && y->insertlength == 0) {
    return -1;
  } else if (y->insertlength > 0 && x->insertlength == 0) {
    return +1;
  } else if (x->insertlength < y->insertlength) {
    return -1;
  } else if (y->insertlength < x->insertlength) {
    return +1;
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else {
    return 0;
  }
}



static float
Stage3end_compute_mapq (Stage3end_T this, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			char *quality_string, bool trim_terminals_p) {

  if (this == NULL) {
    return 0.0;

  } else if (this->hittype == GMAP) {
    this->mapq_loglik = Pair_compute_mapq(this->pairarray,this->npairs,
					  this->trim_left,this->trim_right,this->querylength_adj,
					  quality_string,trim_terminals_p);

  } else if (this->plusp == true) {
    this->mapq_loglik =
      Substring_compute_mapq(this->substring1,query_compress_fwd,quality_string,trim_terminals_p);

    if (this->substring2 != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substring2,query_compress_fwd,
			       quality_string,trim_terminals_p);
    }
    if (this->substring0 != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substring0,query_compress_fwd,
			       quality_string,trim_terminals_p);
    }

  } else {
    this->mapq_loglik =
      Substring_compute_mapq(this->substring1,query_compress_rev,
			     quality_string,trim_terminals_p);

    if (this->substring2 != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substring2,query_compress_rev,
			       quality_string,trim_terminals_p);
    }
    if (this->substring0 != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substring0,query_compress_rev,
			       quality_string,trim_terminals_p);
    }
  }

  return this->mapq_loglik;
}



static void
Stage3end_display_prep (Stage3end_T this, char *query, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			Genome_T genome) {
  char *deletion_ignore;

  if (this != NULL) {
    debug0(printf("Doing a display prep of end %p\n",this));
    if (this->hittype == GMAP) {
      this->nmismatches_refdiff = this->nmismatches_bothdiff;

    } else if (this->hittype == DELETION) {
      this->nmismatches_refdiff = 
	Substring_display_prep(&this->deletion,this->substring1,query,query_compress_fwd,query_compress_rev,
			       genome,/*deletion_pos*/this->indel_pos,
			       /*deletion_length*/this->nindels);
    } else {
      this->nmismatches_refdiff = 
	Substring_display_prep(&deletion_ignore,this->substring1,query,query_compress_fwd,query_compress_rev,
			       genome,/*deletion_pos*/-1,/*deletion_length*/0);
    }

    if (this->substring2 != NULL) {
      this->nmismatches_refdiff +=
	Substring_display_prep(&deletion_ignore,this->substring2,query,query_compress_fwd,query_compress_rev,
			       genome,/*deletion_pos*/-1,/*deletion_length*/0);
    }
    if (this->substring0 != NULL) {
      this->nmismatches_refdiff +=
	Substring_display_prep(&deletion_ignore,this->substring0,query,query_compress_fwd,query_compress_rev,
			       genome,/*deletion_pos*/-1,/*deletion_length*/0);
    }
  }
  return;
}


static int
end_matches_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->nmatches > y->nmatches) {
    return -1;
  } else if (x->nmatches < y->nmatches) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Stage3end_sort_bymatches (List_T hits) {
  List_T sorted = NULL;
  T *array;
  int n, i;

  n = List_length(hits);
  if (n == 0) {
    return (List_T) NULL;
  } else {
    array = (T *) List_to_array(hits,NULL);
    List_free(&hits);

    qsort(array,n,sizeof(T),end_matches_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,(void *) array[i]);
    }
    FREE(array);

    return sorted;
  }
}


/* Need to include criteria from end_matches_cmp to work on most likely terminals */
static int
paired_seenp_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->paired_seenp > y->paired_seenp) {
    return -1;
  } else if (y->paired_seenp > x->paired_seenp) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (x->nmatches < y->nmatches) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Stage3end_sort_by_paired_seenp (List_T hits) {
  List_T sorted = NULL;
  T *array;
  int n, i;

  n = List_length(hits);
  if (n == 0) {
    return (List_T) NULL;
  } else {
    array = (T *) List_to_array(hits,NULL);
    List_free(&hits);

    qsort(array,n,sizeof(T),paired_seenp_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,(void *) array[i]);
    }
    FREE(array);

    return sorted;
  }
}



Stage3end_T *
Stage3end_eval_and_sort (int *npaths, int *first_absmq, int *second_absmq,
			 Stage3end_T *stage3array, int maxpaths, Shortread_T queryseq,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 Genome_T genome, char *quality_string, bool displayp) {
  char *query;
  float maxlik, loglik;
  float total, q;		/* For Bayesian mapq calculation */
  int compute_npaths;
  int i;
  bool non_terminal_p;

  if (*npaths == 0) {
    /* Skip */
    *first_absmq = 0;
    *second_absmq = 0;

  } else if (*npaths == 1) {
    stage3array[0]->mapq_loglik = MAPQ_MAXIMUM_SCORE;
    stage3array[0]->mapq_score = 
      MAPQ_max_quality_score(Shortread_quality_string(queryseq),Shortread_fulllength(queryseq));
    stage3array[0]->absmq_score = MAPQ_MAXIMUM_SCORE;

    if (displayp == true) {
      query = Shortread_fullpointer_uc(queryseq);
      Stage3end_display_prep(stage3array[0],query,query_compress_fwd,query_compress_rev,
			     genome);
    }
    *first_absmq = stage3array[0]->absmq_score;
    *second_absmq = 0;

  } else {
    /* Determine whether to trim terminal ends */
    non_terminal_p = false;
    for (i = 0; i < *npaths; i++) {
      if (stage3array[i]->hittype != TERMINAL) {
	non_terminal_p = true;
      }
    }

    /* Compute mapq_loglik */
    for (i = 0; i < *npaths; i++) {
      Stage3end_compute_mapq(stage3array[i],query_compress_fwd,query_compress_rev,
			     quality_string,/*trim_terminals_p*/non_terminal_p ? false : true);
    }

    /* Sort by nmatches, then mapq.  Enforce monotonicity. */
    qsort(stage3array,*npaths,sizeof(Stage3end_T),Stage3end_output_cmp);
    for (i = *npaths - 1; i > 0; i--) {
      if (stage3array[i-1]->mapq_loglik < stage3array[i]->mapq_loglik) {
	stage3array[i-1]->mapq_loglik = stage3array[i]->mapq_loglik;
      }
    }
    maxlik = stage3array[0]->mapq_loglik;
    
    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < *npaths; i++) {
      stage3array[i]->mapq_loglik -= maxlik;
    }

    /* Save on computation if possible */
    if (*npaths < maxpaths) {
      compute_npaths = *npaths;
    } else {
      compute_npaths = maxpaths;
    }
    if (compute_npaths < 2) {
      compute_npaths = 2;
    }

    /* Compute absolute mapq */
    for (i = 0; i < compute_npaths; i++) {
      loglik = stage3array[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      stage3array[i]->absmq_score = rint(loglik);
    }
    *first_absmq = stage3array[0]->absmq_score;
    *second_absmq = stage3array[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < *npaths; i++) {
      total += (stage3array[i]->mapq_loglik = fasterexp(stage3array[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < compute_npaths; i++) {
      stage3array[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < compute_npaths; i++) {
      if ((q = 1.0 - stage3array[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	stage3array[i]->mapq_score = 96;
      } else {
	stage3array[i]->mapq_score = rint(-10.0 * log10(q));
      }
    }

    if (displayp == true) {
      /* Prepare for display */
      query = Shortread_fullpointer_uc(queryseq);
      for (i = 0; i < compute_npaths; i++) {
	Stage3end_display_prep(stage3array[i],query,query_compress_fwd,query_compress_rev,
			       genome);
      }
    }

#if 0
    /* Apply filtering for mapq unique -- currently not used since mapq_unique_score is high */
    if (stage3array[0]->mapq_score >= mapq_unique_score &&
	stage3array[1]->mapq_score < mapq_unique_score) {
      for (i = 1; i < *npaths; i++) {
	Stage3end_free(&(stage3array[i]));
      }
      *npaths = 1;
    }
#endif
  }

  return stage3array;
}


static int
insertlength_expected (int insertlength) {
  if (insertlength < expected_pairlength_low) {
    return -1;
  } else if (insertlength > expected_pairlength_very_high) {
    return -1;
  } else if (insertlength > expected_pairlength_high) {
    return 0;
  } else {
    return +1;
  }
}


/* For concordant ends */
static Chrpos_T
pair_insert_length (Stage3end_T hit5, Stage3end_T hit3) {

  if (hit5->plusp != hit3->plusp) {
    debug10(printf("pair_insert_length: hit5->plusp %d != hit3->plusp %d, so returning 0\n",
		   hit5->plusp,hit3->plusp));
    return 0;
  }

  if (hit5->chrnum != 0 && hit3->chrnum != 0) {
    if (Substring_overlap_p(hit5->substring1,hit3->substring1)) {
      return Substring_insert_length(hit5->substring1,hit3->substring1);
    } else if (hit5->substring2 != NULL && Substring_overlap_p(hit5->substring2,hit3->substring1)) {
      return Substring_insert_length(hit5->substring2,hit3->substring1);
    } else if (hit5->substring0 != NULL && Substring_overlap_p(hit5->substring0,hit3->substring1)) {
      return Substring_insert_length(hit5->substring0,hit3->substring1);
    }
  
    if (hit3->substring2 != NULL) {
      if (Substring_overlap_p(hit5->substring1,hit3->substring2)) {
	return Substring_insert_length(hit5->substring1,hit3->substring2);
      } else if (hit5->substring2 != NULL && Substring_overlap_p(hit5->substring2,hit3->substring2)) {
	return Substring_insert_length(hit5->substring2,hit3->substring2);
      } else if (hit5->substring0 != NULL && Substring_overlap_p(hit5->substring0,hit3->substring2)) {
	return Substring_insert_length(hit5->substring0,hit3->substring2);
      }
    }

    if (hit3->substring0 != NULL) {
      if (Substring_overlap_p(hit5->substring1,hit3->substring0)) {
	return Substring_insert_length(hit5->substring1,hit3->substring0);
      } else if (hit5->substring2 != NULL && Substring_overlap_p(hit5->substring2,hit3->substring0)) {
	return Substring_insert_length(hit5->substring2,hit3->substring0);
      } else if (hit5->substring0 != NULL && Substring_overlap_p(hit5->substring0,hit3->substring0)) {
	return Substring_insert_length(hit5->substring0,hit3->substring0);
      }
    }
  }

  /* No overlap found between any combination of substrings */
  if (hit5->plusp == true) {
    if (hit5->genomicend > hit3->genomicstart + hit5->querylength_adj + hit3->querylength_adj) {
      debug10(printf("pair_insert_length: no overlap found, and %lu - %lu + %d + %d < 0, so returning 0\n",
		     hit3->genomicstart,hit5->genomicend,hit5->querylength_adj,hit3->querylength_adj));
      return 0;
    } else {
      debug10(printf("pair_insert_length: no overlap found, so returning %lu - %lu + %d + %d\n",
		     hit3->genomicstart,hit5->genomicend,hit5->querylength_adj,hit3->querylength_adj));
    }
    return hit3->genomicstart - hit5->genomicend + hit5->querylength_adj + hit3->querylength_adj;
  } else {
    if (hit3->genomicstart > hit5->genomicend + hit5->querylength_adj + hit3->querylength_adj) {
      debug10(printf("pair_insert_length: no overlap found, and %lu - %lu + %d + %d < 0, so returning 0\n",
		     hit5->genomicend,hit3->genomicstart,hit5->querylength_adj,hit3->querylength_adj));
      return 0;
    } else {
      debug10(printf("pair_insert_length: no overlap found, so returning %lu - %lu + %d + %d\n",
		     hit5->genomicend,hit3->genomicstart,hit5->querylength_adj,hit3->querylength_adj));
      return hit5->genomicend - hit3->genomicstart + hit5->querylength_adj + hit3->querylength_adj;
    }
  }
}



/* For unpaired ends */
static Chrpos_T
pair_insert_length_unpaired (Stage3end_T hit5, Stage3end_T hit3) {

  if (hit5->chrnum != hit3->chrnum) {
    debug10(printf("pair_insert_length: hit5->plusp %d != hit3->plusp %d, so returning 0\n",
		   hit5->plusp,hit3->plusp));
    return 0;
  } else if (hit5->high < hit3->low) {
    return hit3->low - hit5->high + hit5->querylength_adj + hit3->querylength_adj;
  } else if (hit3->high < hit5->low) {
    return hit5->low - hit3->high + hit5->querylength_adj + hit3->querylength_adj;
  } else {
    return hit5->querylength_adj + hit3->querylength_adj;
  }
}


Stage3end_T *
Stage3end_eval_and_sort_guided (int *npaths, int *first_absmq, int *second_absmq, Stage3end_T guide,
				Stage3end_T *stage3array, int maxpaths, Shortread_T queryseq,
				Compress_T query_compress_fwd, Compress_T query_compress_rev,
				Genome_T genome, char *quality_string, bool displayp) {
  char *query;
  float maxlik, loglik;
  float total, q;		/* For Bayesian mapq calculation */
  int compute_npaths;
  int i;
  bool non_terminal_p;

  if (*npaths == 0) {
    /* Skip */
    *first_absmq = 0;
    *second_absmq = 0;

  } else if (*npaths == 1) {
    stage3array[0]->mapq_loglik = MAPQ_MAXIMUM_SCORE;
    stage3array[0]->mapq_score = 
      MAPQ_max_quality_score(Shortread_quality_string(queryseq),Shortread_fulllength(queryseq));
    stage3array[0]->absmq_score = MAPQ_MAXIMUM_SCORE;

    if (displayp == true) {
      query = Shortread_fullpointer_uc(queryseq);
      Stage3end_display_prep(stage3array[0],query,query_compress_fwd,query_compress_rev,
			     genome);
    }
    *first_absmq = stage3array[0]->absmq_score;
    *second_absmq = 0;

  } else {
    /* Determine whether to trim terminal ends */
    non_terminal_p = false;
    for (i = 0; i < *npaths; i++) {
      if (stage3array[i]->hittype != TERMINAL) {
	non_terminal_p = true;
      }
    }

    /* Compute mapq_loglik */
    for (i = 0; i < *npaths; i++) {
      Stage3end_compute_mapq(stage3array[i],query_compress_fwd,query_compress_rev,
			     quality_string,/*trim_terminals_p*/non_terminal_p ? false : true);
    }

    /* Compute insert_length relative to guide */
    for (i = 0; i < *npaths; i++) {
      stage3array[i]->guided_insertlength = pair_insert_length_unpaired(stage3array[i],guide);
    }

    /* Sort by nmatches, then mapq.  Enforce monotonicity. */
    qsort(stage3array,*npaths,sizeof(Stage3end_T),Stage3end_output_cmp);
    for (i = *npaths - 1; i > 0; i--) {
      if (stage3array[i-1]->mapq_loglik < stage3array[i]->mapq_loglik) {
	stage3array[i-1]->mapq_loglik = stage3array[i]->mapq_loglik;
      }
    }
    maxlik = stage3array[0]->mapq_loglik;
    
    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < *npaths; i++) {
      stage3array[i]->mapq_loglik -= maxlik;
    }

    /* Save on computation if possible */
    if (*npaths < maxpaths) {
      compute_npaths = *npaths;
    } else {
      compute_npaths = maxpaths;
    }
    if (compute_npaths < 2) {
      compute_npaths = 2;
    }

    /* Compute absolute mapq */
    for (i = 0; i < compute_npaths; i++) {
      loglik = stage3array[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      stage3array[i]->absmq_score = rint(loglik);
    }
    *first_absmq = stage3array[0]->absmq_score;
    *second_absmq = stage3array[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < *npaths; i++) {
      total += (stage3array[i]->mapq_loglik = fasterexp(stage3array[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < compute_npaths; i++) {
      stage3array[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < compute_npaths; i++) {
      if ((q = 1.0 - stage3array[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	stage3array[i]->mapq_score = 96;
      } else {
	stage3array[i]->mapq_score = rint(-10.0 * log10(q));
      }
    }

    if (displayp == true) {
      /* Prepare for display */
      query = Shortread_fullpointer_uc(queryseq);
      for (i = 0; i < compute_npaths; i++) {
	Stage3end_display_prep(stage3array[i],query,query_compress_fwd,query_compress_rev,
			       genome);
      }
    }

#if 0
    /* Apply filtering for mapq unique -- currently not used since mapq_unique_score is high */
    if (stage3array[0]->mapq_score >= mapq_unique_score &&
	stage3array[1]->mapq_score < mapq_unique_score) {
      for (i = 1; i < *npaths; i++) {
	Stage3end_free(&(stage3array[i]));
      }
      *npaths = 1;
    }
#endif
  }

  return stage3array;
}


/* Note: single-end terminals can be present with non-terminals when
   paired-end reads are searched for concordance, which can accumulate
   terminal alignments */

/* Pre-final: max (max-terminal, min-other)
   Final: max (min-terminal, max-GMAP, min-other) */


static List_T
Stage3end_optimal_score_aux (bool *eliminatedp, List_T hitlist, int cutoff_level, int suboptimal_mismatches,
			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     bool keep_gmap_p, bool finalp) {
  List_T optimal = NULL, p;
  T hit;
  int n;
  int minscore = MAX_READLENGTH;
  int max_nmatches = 0, max_nmatches_posttrim = 0;
  int trim_left, trim_right;
  int min_trim_left = MAX_READLENGTH, min_trim_right = MAX_READLENGTH;
  int max_trim_left_terminal = 0, max_trim_right_terminal = 0;
  int nindelbreaks;

#ifdef TRANSLOC_SPECIAL
  bool non_translocation_p = false;
#endif


  *eliminatedp = false;
  n = List_length(hitlist);
  debug4(printf("\nEntered Stage3end_optimal_score with %d hits: %s\n",
		n,finalp == true ? "FINAL" : "not final"));

  if (n <= 1) {
    return hitlist;
  }

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
#ifdef TRANSLOC_SPECIAL
    if (hit->chrnum != 0) {
      non_translocation_p = true;
    }
#endif
#if 0
    if (hit->hittype == GMAP) {
      debug4(printf("Found gmap/terminal\n"));
      gmap_terminal_p = true;
    } else if (hit->hittype == TERMINAL) {
      debug4(printf("Found gmap/terminal\n"));
      gmap_terminal_p = true;
    } else {
      debug4(printf("Found a non-gmap/terminal\n"));
      non_gmap_terminal_p = true;
    }
#endif
  }


  /* Use eventrim for comparing alignments */
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    debug4(printf("hittype: %s, trim_left: %d, trim_right %d\n",
		  hittype_string(hit->hittype),hit->trim_left,hit->trim_right));
    if (hit->hittype == TERMINAL) {
      /* Don't allow terminals to set trims */
#if 0
      if (hit->trim_left > max_trim_left_terminal) {
	max_trim_left_terminal = hit->trim_left;
      }
      if (hit->trim_right > max_trim_right_terminal) {
	max_trim_right_terminal = hit->trim_right;
      }
#endif

    } else if ((hit->hittype == INSERTION || hit->hittype == DELETION) &&
	       (hit->indel_pos < 15 || hit->indel_pos > hit->querylength_adj - 15)) {
      /* Don't allow end indels to set trims */

    } else {
      if (hit->trim_left_splicep == true) {
	if (hit->trim_left > max_trim_left_terminal) {
	  max_trim_left_terminal = hit->trim_left;
	}
      } else if (hit->trim_left < min_trim_left) {
	min_trim_left = hit->trim_left;
      }
      if (hit->trim_right_splicep == true) {
	if (hit->trim_right > max_trim_right_terminal) {
	  max_trim_right_terminal = hit->trim_right;
	}
      } else if (hit->trim_right < min_trim_right) {
	min_trim_right = hit->trim_right;
      }
    }
  }

  if (min_trim_left == MAX_READLENGTH) {
    trim_left = max_trim_left_terminal;
  } else {
    trim_left = (max_trim_left_terminal > min_trim_left) ? max_trim_left_terminal : min_trim_left;
  }
  if (min_trim_right == MAX_READLENGTH) {
    trim_right = max_trim_right_terminal;
  } else {
    trim_right = (max_trim_right_terminal > min_trim_right) ? max_trim_right_terminal : min_trim_right;
  }

  debug4(printf("non-terminals: min_trim_left: %d, min_trim_right %d\n",min_trim_left,min_trim_right));
  debug4(printf("prefinal-terminals: max_trim_left: %d, max_trim_right %d\n",
		max_trim_left_terminal,max_trim_right_terminal));
  debug4(printf("trim_left: %d, trim_right %d\n",trim_left,trim_right));

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    if (hit->hittype == TERMINAL && finalp == false) {
      /* Ignore */
      hit->score_eventrim = 0;
    } else if (hit->hittype == GMAP) {
      hit->score_eventrim = 0;  /* was hit->penalties */
      debug4(printf("score GMAP:"));
#if 0
      if (Stage3end_bad_stretch_p(hit,query_compress_fwd,query_compress_rev) == true) {
	hit->score_eventrim += 2;
	debug4(printf("  bad stretch 2."));
      }
#endif

      if (0 && hit->trim_left <= 8) {
	/* Ignore small trims */
      } else if (hit->trim_left > trim_left) {
	hit->score_eventrim += hit->trim_left - trim_left;
	debug4(printf("  add trim left (%d - %d).",hit->trim_left,trim_left));
      }
      if (0 && hit->trim_right <= 8) {
	/* Ignore small trims */
      } else if (hit->trim_right > trim_right) {
	hit->score_eventrim += hit->trim_right - trim_right;
	debug4(printf("  add trim right (%d - %d).",hit->trim_right,trim_right));
      }

      hit->score_eventrim += Pair_nmismatches_region(&nindelbreaks,hit->pairarray,hit->npairs,
						     trim_left,trim_right,hit->querylength_adj);
      debug4(printf("  add nmismatches %d.",Pair_nmismatches_region(&nindelbreaks,hit->pairarray,hit->npairs,
								    trim_left,trim_right,hit->querylength_adj)));
      hit->score_eventrim += indel_penalty_middle * nindelbreaks;
      hit->score_eventrim += hit->amb_nmatches_start / ambig_end_interval;
      debug4(printf("  add amb start %d/%d.",hit->amb_nmatches_start,ambig_end_interval));
      hit->score_eventrim += hit->amb_nmatches_end / ambig_end_interval;
      debug4(printf("  add amb end %d/%d.",hit->amb_nmatches_end,ambig_end_interval));
      debug4(printf("  RESULT: %d\n",hit->score_eventrim));

    } else {
      debug4(printf("score OTHER:"));
      hit->score_eventrim = hit->penalties;
      debug4(printf("  penalties %d.",hit->penalties));

      hit->score_eventrim += Substring_count_mismatches_region(hit->substring0,trim_left,trim_right,
							      query_compress_fwd,query_compress_rev);
      debug4(printf("  substring 0 %d.",Substring_count_mismatches_region(hit->substring0,trim_left,trim_right,
									  query_compress_fwd,query_compress_rev)));

      hit->score_eventrim += Substring_count_mismatches_region(hit->substring1,trim_left,trim_right,
							       query_compress_fwd,query_compress_rev);
      debug4(printf("  substring 1 %d.",Substring_count_mismatches_region(hit->substring1,trim_left,trim_right,
									  query_compress_fwd,query_compress_rev)));

      hit->score_eventrim += Substring_count_mismatches_region(hit->substring2,trim_left,trim_right,
							       query_compress_fwd,query_compress_rev);
      debug4(printf("  substring 2 %d.",Substring_count_mismatches_region(hit->substring2,trim_left,trim_right,
									  query_compress_fwd,query_compress_rev)));

      if (hit->hittype == INSERTION || hit->hittype == DELETION) {
	hit->score_eventrim += indel_penalty_middle;
	debug4(printf("  add indel %d.",indel_penalty_middle));
      }
      debug4(printf("  RESULT: %d\n",hit->score_eventrim));
    }
  }

  /* Compute minscore */
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->hittype == TERMINAL && finalp == false) {
      /* Don't use to determine minscore */
#ifdef TRANSLOC_SPECIAL
    } else if (hit->chrnum == 0 && non_translocation_p == true) {
      /* Skip, since we will eliminate */
#endif
    } else {
      if (hit->nmatches > max_nmatches) {
	max_nmatches = hit->nmatches;
	max_nmatches_posttrim = hit->nmatches_posttrim;
      }
#ifdef TERMINAL_SECOND_CLASS
      if (non_gmap_terminal_p == true && (hit->hittype == TERMINAL || hit->hittype == GMAP)) {
	/* Skip from setting minscore */
      }
#endif
      if (hit->score_eventrim < minscore) {
	minscore = hit->score_eventrim;
      }
    }
  }

  debug4(printf("Stage3end_optimal_score over %d hits: minscore = %d + subopt:%d\n",
		n,minscore,suboptimal_mismatches));
  minscore += suboptimal_mismatches;
  max_nmatches -= suboptimal_mismatches;
  max_nmatches_posttrim -= suboptimal_mismatches;

#if 0
  if (non_gmap_terminal_p == false && minscore > cutoff_level) {
    /* If we are down to GMAP or terminal hits, keep at least one hit */
    cutoff_level = minscore;
  }
#else
  cutoff_level = minscore;
#endif

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    if (hit->hittype == TERMINAL && finalp == false) {
      debug4(printf("Keeping a hit of type TERMINAL\n"));
      optimal = List_push(optimal,hit);
      
    } else if (keep_gmap_p == true && hit->hittype == GMAP) {
      /* GMAP hits already found to be better than their corresponding terminals */
      debug4(printf("Keeping a hit of type GMAP\n"));
      optimal = List_push(optimal,hit);

#ifdef TRANSLOC_SPECIAL
    } else if (hit->chrnum == 0 && non_translocation_p == true) {
      debug4(printf("Eliminating a hit with splice translocation\n"));
      *eliminatedp = true;
      Stage3end_free(&hit);
#endif

#ifdef TERMINAL_SECOND_CLASS
    } else if ((hit->hittype == TERMINAL || hit->hittype == GMAP) &&
	       non_gmap_terminal_p == true) {
      if (hit->nmatches >= max_nmatches) {
	debug4(printf("Keeping a terminal with nmatches %d\n",hit->nmatches));
	optimal = List_push(optimal,(void *) hit);
      } else {
	debug4(printf("Eliminating a terminal where non-terminals are present\n"));
	*eliminatedp = true;
	Stage3end_free(&hit);
      }
#endif

    } else if (hit->score_eventrim > cutoff_level) {
      /* For dibasep were previously using hit->ntscore, but gives false positives */
      debug4(printf("Eliminating a hit of type %s with score_eventrim %d > cutoff_level %d\n",
		    hittype_string(hit->hittype),hit->score_eventrim,cutoff_level));
      *eliminatedp = true;
      Stage3end_free(&hit);

    } else if (hit->score_eventrim > minscore /* && hit->nmatches_posttrim < max_nmatches_posttrim */) {
      debug4(printf("Eliminating a hit with score_eventrim %d and type %s\n",
		    hit->score_eventrim,hittype_string(hit->hittype)));
      *eliminatedp = true;
      Stage3end_free(&hit);

    } else {
      debug4(printf("Keeping a hit with score_eventrim %d and type %s\n",
		    hit->score_eventrim,hittype_string(hit->hittype)));
      optimal = List_push(optimal,hit);
    }
  }
  
  List_free(&hitlist);

  debug4(printf("hitlist now has %d entries\n",List_length(optimal)));
  return optimal;
}


List_T
Stage3end_optimal_score (List_T hitlist, int cutoff_level, int suboptimal_mismatches,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 bool keep_gmap_p, bool finalp) {
  List_T optimal;
  bool eliminatedp;


  optimal = Stage3end_optimal_score_aux(&eliminatedp,hitlist,cutoff_level,suboptimal_mismatches,
					query_compress_fwd,query_compress_rev,
					keep_gmap_p,finalp);
  while (eliminatedp == true) {
    optimal = Stage3end_optimal_score_aux(&eliminatedp,optimal,cutoff_level,suboptimal_mismatches,
					  query_compress_fwd,query_compress_rev,
					  keep_gmap_p,finalp);
  }

  return optimal;
}


List_T
Stage3end_remove_circular_alias (List_T hitlist) {
  List_T newlist = NULL, p;
  T hit;

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->chrnum == 0) {
      /* Distant splice */
      newlist = List_push(newlist,(void *) hit);

    } else if (hit->low >= hit->chroffset + hit->chrlength) {
      /* All in circular alias */
      Stage3end_free(&hit);

    } else {
      newlist = List_push(newlist,(void *) hit);
    }
  }

  List_free(&hitlist);
  return newlist;
}


int
Stage3end_noptimal (List_T hitlist) {
  int noptimal;
  List_T p;
  T hit;
  int minscore = MAX_READLENGTH;

  noptimal = 0;
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->score < minscore) {
      minscore = hit->score;
      noptimal = 0;
    }
    if (hit->score == minscore) {
      noptimal++;
    }
  }

  return noptimal;
}


static int
duplicate_sort_cmp (const void *a, const void *b) {
  int cmp;
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (x->hittype > y->hittype) {
    return +1;
  } else if (x->genomicend < y->genomicend) {
    return -1;
  } else if (x->genomicend > y->genomicend) {
    return +1;
  } else if ((cmp = Substring_compare(x->substring1,y->substring1)) != 0) {
    return cmp;
  } else if ((cmp = Substring_compare(x->substring2,y->substring2)) != 0) {
    return cmp;
  } else if ((cmp = Substring_compare(x->substring0,y->substring0)) != 0) {
    return cmp;
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
  } else if (x->sarrayp == true && y->sarrayp == false) {
    return -1;
  } else if (x->sarrayp == false && y->sarrayp == true) {
    return +1;
  } else {
    return 0;
  }
}

/* Same as duplicate_sort_cmp, except for indel_low */
static int
duplicate_equiv_cmp (const void *a, const void *b) {
  int cmp;
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (x->hittype > y->hittype) {
    return +1;
  } else if (x->genomicend < y->genomicend) {
    return -1;
  } else if (x->genomicend > y->genomicend) {
    return +1;
  } else if ((cmp = Substring_compare(x->substring1,y->substring1)) != 0) {
    return cmp;
  } else if ((cmp = Substring_compare(x->substring2,y->substring2)) != 0) {
    return cmp;
  } else if ((cmp = Substring_compare(x->substring0,y->substring0)) != 0) {
    return cmp;
#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif
#if 0
    /* Want to sort by sarrayp, but still consider them equal */
  } else if (x->sarrayp == true && y->sarrayp == false) {
    return -1;
  } else if (x->sarrayp == false && y->sarrayp == true) {
    return +1;
#endif
  } else {
    return 0;
  }
}


static int
genomicstart_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else if (x->genomicend < y->genomicend) {
    return -1;
  } else if (x->genomicend > y->genomicend) {
    return +1;
  } else {
    return 0;
  }
}

static int
genomicend_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->genomicend < y->genomicend) {
    return -1;
  } else if (x->genomicend > y->genomicend) {
    return +1;
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else {
    return 0;
  }
}


#if 0

List_T
Stage3end_mark_ambiguous_splices (bool *ambiguousp, List_T hitlist) {
  T x, y, *hits;
  int n, i, j;
  Chrpos_T splice_distance_1, splice_distance_2;

#ifndef CLEAN_SINGLE_END_AMBIGUOUS
  return hitlist;
#endif

  *ambiguousp = false;
  n = List_length(hitlist);
  debug9(printf("Entered Stage3end_mark_ambiguous_splices with %d pairs\n",n));
  if (n == 0) {
    return NULL;
  } else {
    hits = (T *) List_to_array(hitlist,NULL);
  }

  /* By genomicstart */
  debug9(printf("Stage3end_mark_ambiguous_splices: checking %d hits by genomicstart\n",n));
  qsort(hits,n,sizeof(T),genomicstart_cmp);
  for (i = 0; i < n; i++) {
    x = hits[i];
    if (x->hittype == SPLICE) {
      j = i+1;
      while (j < n && hits[j]->genomicstart == x->genomicstart) {
	y = hits[j];
	if (y->hittype == SPLICE) {
	  debug9(printf("  #%d has common start with #%d\n",i,j));
	  if (SENSE_INCONSISTENT_P(x->sensedir,y->sensedir)) {
	    debug9(printf("  #%d and #%d have different sense, so skipping\n",i,j));
	  } else if (x->score != y->score) {
	    debug9(printf("  #%d overlaps #%d and score %d != %d, so skipping\n",
			  i,j,x->score,y->score));
	  } else {
	    debug9(printf("  #%d overlaps #%d and score %d == %d, so both are ambiguous\n",
			  i,j,x->score,y->score));
	    if (x->genomicstart < x->genomicend) {
	      splice_distance_1 = x->genomicend - x->genomicstart;
	    } else {
	      splice_distance_1 = x->genomicstart - x->genomicend;
	    }

	    if (y->genomicstart < y->genomicend) {
	      splice_distance_2 = y->genomicend - y->genomicstart;
	    } else {
	      splice_distance_2 = y->genomicstart - y->genomicend;
	    }

	    debug9(printf("    splice distances: %u and %u\n",splice_distance_1,splice_distance_2));
	    if (splice_distance_1 > splice_distance_2) {
	      debug9(printf("    first is ambiguous\n"));
	      x->chimera_ambiguous_p = true;
	      *ambiguousp = true;
	    } else if (splice_distance_2 > splice_distance_1) {
	      debug9(printf("    second is ambiguous\n"));
	      y->chimera_ambiguous_p = true;
	      *ambiguousp = true;
	    }
	  }
	}
	j++;
      }
      i = j;
    }
  }
    
  /* By genomicend */
  debug9(printf("Stage3end_mark_ambiguous_splices: checking %d hits by genomicend\n",n));
  qsort(hits,n,sizeof(T),genomicend_cmp);
  for (i = 0; i < n; i++) {
    x = hits[i];
    if (x->hittype == SPLICE) {
      j = i+1;
      while (j < n && hits[j]->genomicend == x->genomicend) {
	y = hits[j];
	if (y->hittype == SPLICE) {
	  debug9(printf("  #%d has common end with #%d\n",i,j));
	  if (SENSE_INCONSISTENT_P(x->sensedir,y->sensedir)) {
	    debug9(printf("  #%d and #%d have different sense, so skipping\n",i,j));
	  } else if (x->score != y->score) {
	    debug9(printf("  #%d overlaps #%d and score %d != %d, so skipping\n",
			  i,j,x->score,y->score));
	  } else {
	    debug9(printf("  #%d overlaps #%d and score %d == %d, so both are ambiguous\n",
			  i,j,x->score,y->score));
	    if (x->genomicstart < x->genomicend) {
	      splice_distance_1 = x->genomicend - x->genomicstart;
	    } else {
	      splice_distance_1 = x->genomicstart - x->genomicend;
	    }

	    if (y->genomicstart < y->genomicend) {
	      splice_distance_2 = y->genomicend - y->genomicstart;
	    } else {
	      splice_distance_2 = y->genomicstart - y->genomicend;
	    }

	    debug9(printf("    splice distances: %u and %u\n",splice_distance_1,splice_distance_2));
	    if (splice_distance_1 > splice_distance_2) {
	      debug9(printf("    first is ambiguous\n"));
	      x->chimera_ambiguous_p = true;
	      *ambiguousp = true;
	    } else if (splice_distance_2 > splice_distance_1) {
	      debug9(printf("    second is ambiguous\n"));
	      y->chimera_ambiguous_p = true;
	      *ambiguousp = true;
	    }
	  }
	}
	j++;
      }
      i = j;
    }
  }

  debug9(
	 for (i = 0; i < n; i++) {
	   x = hits[i];
	   if (x->hittype == SPLICE) {
	     printf("  %d: %lu..%lu (plusp = %d) ambiguousp:%d\n",
		    i,x->genomicstart,x->genomicend,x->plusp,x->chimera_ambiguous_p);
	   }
	 }
	 );

  FREE(hits);

  return hitlist;
}

#endif

static void
Stage3end_print_substrings (Stage3end_T hit) {
  Substring_print_ends(hit->substring1,hit->chrnum);
  Substring_print_ends(hit->substring2,hit->chrnum);
  Substring_print_ends(hit->substring0,hit->chrnum);
  return;
}


#if 0
static bool
Stage3end_equal_p (Stage3end_T hit5, Stage3end_T hit3) {

  if (Substring_equal_p(hit5->substring1,hit3->substring1) == false) {
    return false;

  } else if (Substring_equal_p(hit5->substring2,hit3->substring2) == false) {
    return false;

  } else if (Substring_equal_p(hit5->substring0,hit3->substring0) == false) {
    return false;

  } else {
    return true;
  }
}
#endif


const Except_T Duplicate_Pairing = { "Duplicates both seen in pairing" };

List_T
Stage3end_remove_duplicates (List_T hitlist, Shortread_T queryseq1, Shortread_T queryseq2) {
#ifdef DEBUG7
  List_T p;
#endif
  T x, y, *hits;
  int n, usedi, i, j, k;
  bool *eliminate, eliminatep;

  n = List_length(hitlist);
  debug7(printf("Entered Stage3end_remove_duplicates with %d hits\n",n));
  if (n == 0) {
    return NULL;
  } else {
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hits = (T *) List_to_array(hitlist,NULL);
  }


  /* By equivalence */
  debug7(printf("Stage3end_remove_duplicates: checking %d hits by equivalence class\n",n));
  qsort(hits,n,sizeof(T),duplicate_sort_cmp);

  debug7(
	 for (i = 0; i < n; i++) {
	   x = hits[i];
	   printf("  Initial %d (%s): %p #%d:%lu..%lu, nmatches %d, score %d ",
		  i,hittype_string(x->hittype),x,x->chrnum,x->genomicstart-x->chroffset,x->genomicend-x->chroffset,
		  x->nmatches,x->score);
	   Stage3end_print_substrings(x);
	   printf("\n");
	 }
	 );

  eliminatep = false;
  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && duplicate_equiv_cmp(&(hits[j]),&(hits[i])) == 0) {
      j++;
    }

    if (j > i+1) {
      debug7(printf("Equivalence class #%d through #%d.  ",i,j-1));

      x = hits[i];
      if (x->paired_usedp == true) {
	usedi = i;
      } else {
	usedi = -1;
      }
      
      for (k = i+1; k < j; k++) {
	y = hits[k];
	if (y->paired_usedp == true) {
	  if (usedi >= 0) {
	    debug7(printf("  #%d equivalent to #%d and both used (%p and %p)\n",k,usedi,hits[k],hits[usedi]));
	    fprintf(stderr,"Duplicates of Stage3end_T both seen\n");
	    Shortread_print_query_pairedend_fasta(stderr,queryseq1,queryseq2,
						  /*invert_first_p*/false,/*invert_second_p*/true);
#ifdef CHECK_ASSERTIONS
	    Except_raise(&Duplicate_Pairing, __FILE__, __LINE__);
#endif
	  } else {
	    usedi = k;
	  }
	}
      }

      if (usedi < 0) {
	debug7(printf("None used yet so eliminating #%d through #%d\n",i+1,j-1));
	for (k = i+1; k < j; k++) {
	  eliminate[k] = true;
	  eliminatep = true;
	}
      } else {
	debug7(printf("One used already so eliminating all but #%d\n",usedi));
	for (k = i; k < j; k++) {
	  if (k != usedi) {
	    eliminate[k] = true;
	    eliminatep = true;
	  }
	}
      }
    }

    i = j;
  }
    

#if 0
  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
  }
#endif

  if (eliminatep == false) {
    debug7(printf("No eliminations, so hitlist is unchanged\n"));
  } else {
    List_free(&hitlist);
    hitlist = (List_T) NULL;
    for (i = n-1; i >= 0; i--) {
      x = hits[i];
      if (eliminate[i] == false) {
	debug7(printf("  Keeping #%d:%lu..%lu, nmatches %d (nindels %d, indel_pos %d, distance %u, chrnum %d) (plusp = %d)\n",
		      x->chrnum,x->genomicstart-x->chroffset,x->genomicend-x->chroffset,
		      x->nmatches,x->nindels,x->indel_pos,x->distance,x->chrnum,x->plusp));
	hitlist = List_push(hitlist,x);
      } else {
	debug7(printf("  Eliminating #%d:%lu..%lu, nmatches %d (nindels %d, indel_pos %d, distance %u, chrnum %d) (plusp = %d)\n",
		      x->chrnum,x->genomicstart-x->chroffset,x->genomicend-x->chroffset,
		      x->nmatches,x->nindels,x->indel_pos,x->distance,x->chrnum,x->plusp));
	Stage3end_free(&x);
      }
    }
  }

  FREE(hits);
  FREE(eliminate);

  debug7(
	 for (p = hitlist, i = 0; p != NULL; p = p->rest, i++) {
	   x = (T) p->first;
	   printf("  Final %d: #%d:%lu..%lu (plusp = %d)\n",
		  i,x->chrnum,x->genomicstart-x->chroffset,x->genomicend-x->chroffset,x->plusp);
	 }
	 );

  debug7(printf("Exited Stage3end_remove_duplicates with %d hits\n",List_length(hitlist)));
  return hitlist;
}


List_T
Stage3end_filter_terminals (List_T hits) {
  List_T filtered = NULL, p;
  T hit;

  for (p = hits; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->hittype == TERMINAL && Stage3end_terminal_length(hit) < terminal_output_minlength) {
      Stage3end_free(&hit);
    } else {
      filtered = List_push(filtered,(void *) hit);
    }
  }

  List_free(&hits);
  return filtered;
}


/* Used for eliminating exact duplicates.  Also sorts secondarily by hittype. */
static int
hit_sort_cmp (const void *a, const void *b) {
  Stage3end_T x = * (Stage3end_T *) a;
  Stage3end_T y = * (Stage3end_T *) b;
  
  if (x->plusp > y->plusp) {
    return -1;
  } else if (y->plusp > x->plusp) {
    return +1;

#if 0
  } else if (x->high < y->low) {
    return -1;
  } else if (y->high < x->low) {
    return +1;
#else
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;

  } else if (x->high < y->high) {
    return -1;
  } else if (y->high < x->high) {
    return +1;
#endif

#if 1
    /* Rank terminals last, so terminals cannot win */
  } else if (x->hittype != TERMINAL && y->hittype == TERMINAL) {
    return -1;
  } else if (x->hittype == TERMINAL && y->hittype != TERMINAL) {
    return +1;
#endif

  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
  } else if (x->nchimera_novel < y->nchimera_novel) {
    return -1;
  } else if (y->nchimera_novel < x->nchimera_novel) {
    return +1;
  } else if (x->nchimera_known > y->nchimera_known) {
    return -1;
  } else if (y->nchimera_known > x->nchimera_known) {
    return +1;
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (y->hittype < x->hittype) {
    return +1;
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
  } else if (x->sarrayp == true && y->sarrayp == false) {
    return -1;
  } else if (x->sarrayp == false && y->sarrayp == true) {
    return +1;
  } else {
    return 0;
  }
}

/* Same as hit_sort_cmp, except for hittype, nmatches_posttrim, and indel_low */
static int
hit_equiv_cmp (Stage3end_T x, Stage3end_T y) {

  if (x->plusp > y->plusp) {
    return -1;
  } else if (y->plusp > x->plusp) {
    return +1;
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high < y->high) {
    return +1;
  } else if (y->high < x->high) {
    return -1;

#if 0
  } else if (x->hittype != TERMINAL && y->hittype == TERMINAL) {
    return -1;
  } else if (x->hittype == TERMINAL && y->hittype != TERMINAL) {
    return +1;
#endif

  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
#if 0
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
#endif
  } else if (x->nchimera_novel < y->nchimera_novel) {
    return -1;
  } else if (y->nchimera_novel < x->nchimera_novel) {
    return +1;
  } else if (x->nchimera_known > y->nchimera_known) {
    return -1;
  } else if (y->nchimera_known > x->nchimera_known) {
    return +1;

#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif

#if 0
  } else if (x->sarrayp == true && y->sarrayp == false) {
    return -1;
  } else if (x->sarrayp == false && y->sarrayp == true) {
    return +1;
#endif

  } else {
    return 0;
  }
}


static int
hit_goodness_cmp (bool *equalp, Stage3end_T hit,
#ifdef DEBUG7
		  int k,
#endif
		  Stage3end_T best_hit, bool finalp) {
  double prob1, prob2;

#ifdef PRE_RESOLVE_MULTIMAPPING
  if (Stage3end_tally(x) > TALLY_RATIO*Stage3end_tally(y)) {
    debug7(printf("  #%d overlaps #%d and tally %ld > %f*%ld, so marking %d for elimination\n",
		  i,j,x->tally,TALLY_RATIO,y->tally,j));
    eliminate[j] = true;
  } else if (Stage3end_tally(y) > TALLY_RATIO*Stage3end_tally(x)) {
    debug7(printf("  #%d overlaps #%d and tally %f*%ld < %ld, so marking %d for elimination\n",
		  i,j,TALLY_RATIO,x->tally,y->tally,i));
    eliminate[i] = true;
  }
#endif

  *equalp = false;

#if 1
  if (finalp == true) {
    /* Skip */
  } else if (hit->hittype == TERMINAL || best_hit->hittype == TERMINAL) {
    /* Do not allow terminal to win or lose in pre-final stages */
    debug7(printf(" => %d ties by terminal\n",k));
    return 0;
  }
#endif

#if 0
  if (hit->hittype == TERMINAL || best_hit->hittype == TERMINAL) {
    /* Skip: Don't use scores if terminal is involved */
  } else if (hit->nindels == 0 && best_hit->nindels == 0) {
    /* Skip: Use scores only if indel is involved */
  } else if (hit->score > best_hit->score) {
    debug7(printf("  => %d loses by score\n",k));
    return -1;
  } else if (hit->score < best_hit->score) {
    debug7(printf("  => %d wins by score\n",k));
    return +1;
  }
#endif

#if 0
  if (hit->nmatches < best_hit->nmatches) {
    debug7(printf("  => %d loses by nmatches\n",k));
    return -1;
  } else if (hit->nmatches > best_hit->nmatches) {
    debug7(printf("  => %d wins by nmatches\n",k));
    return +1;
  }
#endif

  if (hit->nmatches_posttrim < best_hit->nmatches_posttrim) {
    debug7(printf("  => %d loses by nmatches\n",k));
    return -1;
  } else if (hit->nmatches_posttrim > best_hit->nmatches_posttrim) {
    debug7(printf("  => %d wins by nmatches\n",k));
    return +1;

  } else if (hit->nchimera_novel > best_hit->nchimera_novel) {
    debug7(printf("  => %d loses by nchimera_novel\n",k));
    return -1;
  } else if (hit->nchimera_novel < best_hit->nchimera_novel) {
    debug7(printf("  => %d wins by nchimera_novel\n",k));
    return +1;
    
  } else if (hit->nchimera_known < best_hit->nchimera_known) {
    debug7(printf("  => %d loses by nchimera_known\n",k));
    return -1;
  } else if (hit->nchimera_known > best_hit->nchimera_known) {
    debug7(printf("  => %d wins by nchimera_known\n",k));
    return +1;

  } else if (hit->hittype > best_hit->hittype) {
    debug7(printf("  => %d loses by hittype\n",k));
    return -1;
  } else if (hit->hittype < best_hit->hittype) {
    debug7(printf("  => %d wins by hittype\n",k));
    return +1;

  } else if (hit->nindels > best_hit->nindels) {
    debug7(printf("  => %d loses by nindels\n",k));
    return -1;
  } else if (hit->nindels < best_hit->nindels) {
    debug7(printf("  => %d wins by nindels\n",k));
    return +1;

  } else if (hit->chrnum == 0 && best_hit->chrnum != 0) {
    debug7(printf("  => %d loses because distant splice\n",k));
    return -1;
  } else if (hit->chrnum > 0 && best_hit->chrnum == 0) {
    debug7(printf("  => %d wins because not distant splice\n",k));
    return +1;

  } else if (hit->end_ambiguous_p == true && best_hit->end_ambiguous_p == false) {
    debug7(printf("  => %d loses because end is ambiguous\n",k));
    return -1;
  } else if (hit->end_ambiguous_p == false && best_hit->end_ambiguous_p == true) {
    debug7(printf("  => %d wins because end is not ambiguous\n",k));
    return +1;

  } else if (finalp == false) {
    debug7(printf("  => indistinguishable\n"));
    return 0;

  } else {
    if (hit->hittype == SPLICE && best_hit->hittype == SPLICE) {
      prob1 = Substring_chimera_prob(hit->substring_donor) + Substring_chimera_prob(hit->substring_acceptor);
      prob2 = Substring_chimera_prob(best_hit->substring_donor) + Substring_chimera_prob(best_hit->substring_acceptor);
      if (prob1 < prob2) {
	debug7(printf("  => %d loses by splice prob %f vs %f\n",prob1,prob2));
	return -1;
      } else if (prob1 > prob2) {
	debug7(printf("  => %d wins by splice prob %f vs %f\n",prob1,prob2));
	return +1;
      }
    }

    if (hit->genomiclength > best_hit->genomiclength) {
      debug7(printf("  => %d loses by genomiclength\n",k));
      return -1;
    } else if (hit->genomiclength < best_hit->genomiclength) {
      debug7(printf("  => %d wins by genomiclength\n",k));
      return +1;

    } else {
      debug7(printf("  => equal\n"));
      *equalp = true;
      return 0;
    }
  }
}


static bool
hit_subsumption (Stage3end_T x, Stage3end_T y) {
  if (x->plusp != y->plusp) {
    return false;		/* Different strands */
  } else if (x->low <= y->low && x->high >= y->high) {
    return true;
  } else if (y->low <= x->low && y->high >= x->high) {
    return true;
  } else {
    return false;
  }
}

static bool
hit_bad_superstretch_p (Stage3end_T hit_k, Stage3end_T *hits, int k, int j, bool finalp) {
  int a;
  bool equalp;

  for (a = k+1; a <= j; a++) {
    if (hit_subsumption(hit_k,hits[a]) == true) {
      debug7(printf("Testing %d because stretches over %d",k,a));
      if (hit_goodness_cmp(&equalp,hits[a],
#ifdef DEBUG7
			   a,
#endif
			   hit_k,finalp) > 0 || equalp == true) {
	debug7(printf(" => eliminating\n"));
	return true;
      }
      debug7(printf("\n"));
    }
  }
  return false;
}



List_T
Stage3end_remove_overlaps (List_T hitlist, bool finalp) {
  List_T unique = NULL;
  T best_hit, hit, *hits, *prev;
  int cmp;
  int nkept, n, i, j, k, besti;
  bool *eliminate, equalp;
#ifdef PRE_RESOLVE_MULTIMAPPING
  long int best_tally;
#endif

  n = List_length(hitlist);
  debug7(printf("Entered Stage3end_remove_overlaps with %d hits: %s\n",
		n,finalp == true ? "FINAL" : "not final"));
  if (n == 0) {
    return NULL;
  } else {
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hits = (T *) List_to_array(hitlist,NULL);
    List_free(&hitlist);
  }


  /* Step 1.  Check for exact duplicates */
  debug7(printf("Checking for exact duplicates\n"));
  qsort(hits,n,sizeof(Stage3end_T),hit_sort_cmp);

  debug7(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%lu..%lu, nmatches %d, score %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->nmatches,hit->score);
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    debug7(printf(" %d,%d",i,j));
    while (j < n && hit_equiv_cmp(hits[j],hits[i]) == 0) {
      debug7(printf("  %d is identical to %d => eliminating\n",j,i));
      eliminate[j] = true;
      j++;
    }
    i = j;
  }
  debug7(printf("\n"));


  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    } else if (hits[i]->paired_usedp == true) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hits;
  hits = (Stage3end_T *) CALLOC(nkept,sizeof(Stage3end_T));

  for (i = 0, j = 0; i < n; i++) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug7(printf("  Keeping %lu..%lu, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug7(printf("  Already paired %lu..%lu, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else {
      debug7(printf("  Eliminating %lu..%lu, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      Stage3end_free(&hit);
    }
  }

  FREE(prev);


  /* Step 2: Check for superstretches */
  n = nkept;
  debug7(printf("Checking for superstretches among %d hits within subsumption clusters\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }

  debug7(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches %d, score %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->nmatches,hit->score);
	 }
	 );

  /* Find clusters */
  i = 0;
  while (i < n) {
    j = i;
    if (hits[i]->chrnum != 0) {
      while (j+1 < n && hit_subsumption(hits[i],hits[j+1]) == true) {
	j = j+1;
      }

      if (j > i) {
	debug7(printf("Cluster from %d up through %d\n",i,j));

	/* Find bad superstretches */
	for (k = i; k <= j; k++) {
	  if (hits[k]->chrnum == 0) {
	    /* Skip */
	  } else if (hit_bad_superstretch_p(hits[k],hits,k,j,finalp) == true) {
	    eliminate[k] = true;
	  }
	}
      }
    }

    i = j+1;
  }

  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    } else if (hits[i]->paired_usedp == true) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hits;
  hits = (Stage3end_T *) CALLOC(nkept,sizeof(Stage3end_T));

  for (i = 0, j = 0; i < n; i++) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug7(printf("  Keeping %lu..%lu, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug7(printf("  Already paired %lu..%lu, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else {
      debug7(printf("  Eliminating %lu..%lu, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      Stage3end_free(&hit);
    }
  }

  FREE(prev);


  /* Step 3: Check for best within subsumption clusters */
  n = nkept;
  debug7(printf("Checking for best among %d hits within subsumption clusters\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  /* qsort(hits,n,sizeof(Stage3end_T),hit_sort_cmp); -- No need since original order was kept */
  
  debug7(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches %d, score %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->nmatches,hit->score);
	 }
	 );

  /* Find clusters from left */
  i = 0;
  while (i < n) {
    j = i;
    if (hits[i]->chrnum != 0) {
      while (j+1 < n && hit_subsumption(hits[i],hits[j+1]) == true) {
	j = j+1;
      }

      if (j > i) {
	debug7(printf("Cluster from %d up through %d\n",i,j));

	best_hit = hits[i];
	besti = i;
	debug7(printf("Assume best is %d\n",besti));

	for (k = i+1; k <= j; k++) {
	  if (hits[k]->chrnum == 0) {
	    /* Skip */
	  } else {
	    cmp = hit_goodness_cmp(&equalp,hits[k],
#ifdef DEBUG7
				   k,
#endif
				   best_hit,finalp);
	    debug7(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	    if (cmp > 0) {
	      best_hit = hits[k];
	      besti = k;
	      debug7(printf("Best is now %d\n",besti));
	    }
	  }
	}

	for (k = i; k <= j; k++) {
	  if (k == besti) {
	    /* Skip */
	  } else if (hits[k]->chrnum == 0) {
	    /* Skip */
	  } else if (hit_goodness_cmp(&equalp,hits[k],
#ifdef DEBUG7
				      k,
#endif
				      best_hit,finalp) < 0 || equalp == true) {
	    debug7(printf("  Eliminating hit %d from left, because beaten by %d\n",k,besti));
	    eliminate[k] = true;
	  }
	}
      }
    }
      
    i = j+1;
  }


  /* Find clusters starting from right */
  j = n - 1;
  while (j >= 0) {
    i = j;
    if (hits[j]->chrnum != 0) {
      while (i-1 >= 0 && hit_subsumption(hits[j],hits[i-1]) == true) {
	i = i-1;
      }

      if (i < j) {
	debug7(printf("Cluster from %d down through %d\n",j,i));
	best_hit = hits[i];
	besti = i;
	debug7(printf("Assume best is %d\n",besti));

	for (k = i+1; k <= j; k++) {
	  if (hits[k]->chrnum == 0) {
	    /* Skip */
	  } else {
	    cmp = hit_goodness_cmp(&equalp,hits[k],
#ifdef DEBUG7
				   k,
#endif
				   best_hit,finalp);
	    debug7(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	    if (cmp > 0) {
	      best_hit = hits[k];
	      besti = k;
	      debug7(printf("Best is now %d\n",besti));
	    }
	  }
	}

	for (k = i; k <= j; k++) {
	  if (k == besti) {
	    /* Skip */
	  } else if (hits[k]->chrnum == 0) {
	    /* Skip */
	  } else if (hit_goodness_cmp(&equalp,hits[k],
#ifdef DEBUG7
				      k,
#endif
				      best_hit,finalp) < 0 || equalp == true) {
	    debug7(printf("  Eliminating hit %d from right, because beaten by %d\n",k,besti));
	    eliminate[k] = true;
	  }
	}
      }
    }
      
    j = i-1;
  }


  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    } else if (hits[i]->paired_usedp == true) {
      nkept++;
    }
  }
  if (nkept == 0) {
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hits;
  hits = (Stage3end_T *) CALLOC(nkept,sizeof(Stage3end_T));

  for (i = 0, j = 0; i < n; i++) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug7(printf("  Keeping %lu..%lu, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug7(printf("  Already paired %lu..%lu, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else {
      debug7(printf("  Eliminating %lu..%lu, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      Stage3end_free(&hit);
    }
  }

  FREE(prev);


  /* Step 4: Check for identity */
  n = nkept;
  debug7(printf("Checking for duplicates among %d hits by identity\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  /* qsort(hits,n,sizeof(Stage3end_T),hit_sort_cmp); -- No need since original order was kept */

  debug7(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches %d, score %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		  hit->nmatches,hit->score);
	 }
	 );

  i = 0;
  while (i < n) {
    debug7(printf("Looking at %d with score %d\n",i,hits[i]->score));
    j = i+1;
    while (j < n && hit_equiv_cmp(hits[j],hits[i]) == 0) {
      debug7(printf("  %d equal to %d\n",j,i));
      eliminate[j] = true;
      j++;
    }

    i = j;
  }

  for (i = n-1; i >= 0; i--) {
    hit = hits[i];
    if (eliminate[i] == false) {
      unique = List_push(unique,hit);
    } else if (hit->paired_usedp == true) {
      unique = List_push(unique,hit);
    } else {
      Stage3end_free(&hit);
    }
  }

  FREE(hits);
  FREE(eliminate);


#ifdef PRE_RESOLVE_MULTIMAPPING
  if (use_tally_p == true && tally_iit != NULL) {
    if ((n = List_length(unique)) > 1) {
      hits = (T *) List_to_array(unique,NULL);
      List_free(&unique);

      best_tally = 0;
      for (i = 0; i < n; i++) {
	if (hits[i]->tally < 0) {
	  hits[i]->tally = Stage3end_compute_tally(hits[i]);
	}
	if (hits[i]->tally > best_tally) {
	  best_tally = hits[i]->tally;
	}
      }

      unique = (List_T) NULL;
      for (i = 0; i < n; i++) {
	if (hits[i]->tally < best_tally) {
	  /* Stage3end_free(&(hits[i])); */
	} else {
	  unique = List_push(unique,(void *) hits[i]);
	}
      }

      FREE(hits);
    }
  }
#endif

  debug7(printf("Exited Stage3end_remove_overlaps with %d hits\n",List_length(unique)));
  return unique;
}


List_T
Stage3end_resolve_multimapping (List_T hits) {
  List_T resolve1, resolve2, resolve3, p;
  Stage3end_T hit;

  Overlap_T best_overlap;
  long int best_tally;
  double tally_threshold;
  bool runlengthp;

  if (List_length(hits) <= 1) {
    return hits;
  }

  if (genes_iit == NULL) {
    resolve1 = hits;
  } else {
    best_overlap = NO_KNOWN_GENE;
    for (p = hits; p != NULL; p = p->rest) {
      hit = (Stage3end_T) p->first;
      if ((hit->gene_overlap = Stage3end_gene_overlap(hit)) > best_overlap) {
	best_overlap = hit->gene_overlap;
      }
    }
    if (best_overlap == NO_KNOWN_GENE) {
      resolve1 = hits;
    } else {
      resolve1 = (List_T) NULL;
      for (p = hits; p != NULL; p = p->rest) {
	hit = (Stage3end_T) p->first;
	if (hit->gene_overlap < best_overlap) {
	  Stage3end_free(&hit);
	} else {
	  resolve1 = List_push(resolve1,(void *) hit);
	}
      }
      List_free(&hits);
    }
  }
      
  if (List_length(resolve1) <= 1) {
    return resolve1;
  }

  if (tally_iit == NULL) {
    resolve2 = resolve1;
  } else {
    best_tally = 0L;
    for (p = resolve1; p != NULL; p = p->rest) {
      hit = (Stage3end_T) p->first;
      if ((hit->tally = Stage3end_compute_tally(hit)) > best_tally) {
	best_tally = hit->tally;
      }
    }
    if (best_tally == 0L) {
      resolve2 = resolve1;
    } else {
      resolve2 = (List_T) NULL;
#ifdef USE_TALLY_RATIO
      tally_threshold = (double) best_tally / TALLY_RATIO;
#else
      tally_threshold = 1.0;
#endif
      for (p = resolve1; p != NULL; p = p->rest) {
	hit = (Stage3end_T) p->first;
	if ((double) hit->tally < tally_threshold) {
	  Stage3end_free(&hit);
	} else {
	  resolve2 = List_push(resolve2,(void *) hit);
	}
      }
      List_free(&resolve1);
    }
  }


  if (List_length(resolve2) <= 1) {
    return resolve2;
  }

  if (runlength_iit == NULL) {
    resolve3 = resolve2;
  } else {
    runlengthp = false;
    for (p = resolve2; p != NULL; p = p->rest) {
      hit = (Stage3end_T) p->first;
      if (Stage3end_runlength_p(hit) == true) {
	runlengthp = true;
      }
    }
    if (runlengthp == false) {
      resolve3 = resolve2;
    } else {
      resolve3 = (List_T) NULL;
      for (p = resolve2; p != NULL; p = p->rest) {
	hit = (Stage3end_T) p->first;
	if (Stage3end_runlength_p(hit) == false) {
	  Stage3end_free(&hit);
	} else {
	  resolve3 = List_push(resolve3,(void *) hit);
	}
      }
      List_free(&resolve2);
    }
  }


  return resolve3;
}


static void
print_alignment_info (FILE *fp, int nblocks, int score, int mapq_score, bool sarrayp) {
  fprintf(fp,"segs:%d,align_score:%d,mapq:%d",nblocks,score,mapq_score);
  if (sarrayp == true) {
    fprintf(fp,",method:sa");
  }
  return;
}


Pairtype_T
Stage3_determine_pairtype (T hit5, T hit3) {
  debug14(printf("Entered Stage3_determine_pairtype\n"));
  if (hit5->effective_chrnum != hit3->effective_chrnum) {
    debug14(printf("Returning unpaired\n"));
    return UNPAIRED;
  } else if (hit5->plusp != hit3->plusp) {
    debug14(printf("Returning paired_inversion\n"));
    return PAIRED_INVERSION;
  } else if (hit5->plusp == true) {
    if (hit3->genomicend < hit5->genomicstart) {
      debug14(printf("Returning paired_scramble\n"));
      return PAIRED_SCRAMBLE;
    } else if (hit3->genomicstart > hit5->genomicend + pairmax) {
      debug14(printf("Returning paired_toolong\n"));
      return PAIRED_TOOLONG;
    } else {
      debug14(printf("Returning concordant\n"));
      return CONCORDANT;
    }
  } else {
    if (hit3->genomicend > hit5->genomicstart) {
      debug14(printf("Returning paired_scramble\n"));
      return PAIRED_SCRAMBLE;
    } else if (hit3->genomicstart + pairmax < hit5->genomicend) {
      debug14(printf("Returning paired_toolong\n"));
      return PAIRED_TOOLONG;
    } else {
      debug14(printf("Returning concordant\n"));
      return CONCORDANT;
    }
  }
}


Pairtype_T
Stage3pair_pairtype (Stage3pair_T this) {
  return this->pairtype;
}

bool
Stage3pair_circularp (Stage3pair_T this) {
  return this->circularp;
}



#if 0
static char *
unpaired_type_text (T hit5, T hit3) {
  if (hit5->chrnum != hit3->chrnum) {
    return UNPAIRED_INTERCHROM_TEXT;
  } else if (hit5->plusp != hit3->plusp) {
    return PAIRED_INVERSION_TEXT;
  } else if (hit5->plusp == true) {
    if (hit3->genomicstart < hit5->genomicstart) {
      return PAIRED_SCRAMBLE_TEXT;
    } else {
      return UNPAIRED_TOOLONG_TEXT;
    }
  } else {
    if (hit5->genomicstart < hit3->genomicstart) {
      return PAIRED_SCRAMBLE_TEXT;
    } else {
      return UNPAIRED_TOOLONG_TEXT;
    }
  }
}
#endif


/* Has a copy in pair.c */
static void
print_pair_info (FILE *fp, T hit5, T hit3, int insertlength, int pairscore,
		 Pairtype_T pairtype) {

  assert(hit5->effective_chrnum == hit3->effective_chrnum); /* Same chromosomes */

#if 0
  /* Doesn't hold for paired (inversion) */
  assert(hit5->plusp == hit3->plusp);	/* Same direction */
#endif

  fprintf(fp,"pair_score:%d",pairscore);
  fprintf(fp,",insert_length:%d",insertlength);

  switch (pairtype) {
  case CONCORDANT: break;
  case PAIRED_SCRAMBLE: fprintf(fp,",pairtype:scramble"); break;
  case PAIRED_INVERSION: fprintf(fp,",pairtype:inversion"); break;
  case PAIRED_TOOLONG: fprintf(fp,",pairtype:toolong"); break;
  case CONCORDANT_TRANSLOCATIONS: break;
  case CONCORDANT_TERMINAL: break;
  case PAIRED_UNSPECIFIED: abort();
  case UNPAIRED: abort();
  case UNSPECIFIED: abort();
  }

  return;
}

static void
print_single (FILE *fp, T this, int score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
	      bool invertp, T hit5, T hit3, int insertlength, int pairscore, 
	      Pairtype_T pairtype, int mapq_score) {
  char *chr;
  bool allocp;

  chr = Univ_IIT_label(chromosome_iit,this->chrnum,&allocp);

  fprintf(fp," ");
  Substring_print_single(fp,this->substring1,queryseq,chr,invertp);

  /* Alignment info */
  fprintf(fp,"\t");
  print_alignment_info(fp,/*nblocks*/1,score,mapq_score,this->sarrayp);

  /* Pairing info */
  if (hit5 != NULL && hit3 != NULL) {
    fprintf(fp,"\t");
    print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
  }

  fprintf(fp,"\n");

  if (allocp == true) {
    FREE(chr);
  }

  return;
}


static void
print_insertion (FILE *fp, T this, int score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
		 bool invertp, T hit5, T hit3, int insertlength, int pairscore,
		 Pairtype_T pairtype, int mapq_score) {
  char *chr;
  bool allocp;

  chr = Univ_IIT_label(chromosome_iit,this->chrnum,&allocp);

  fprintf(fp," ");
  Substring_print_insertion_1(fp,this->substring1,this->substring2,this->nindels,
			      queryseq,chr,invertp);
  /* Alignment info */
  fprintf(fp,"\t");
  print_alignment_info(fp,/*nblocks*/2,score,mapq_score,this->sarrayp);

  /* Pairing info */
  if (hit5 != NULL && hit3 != NULL) {
    fprintf(fp,"\t");
    print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
  }

  fprintf(fp,"\n");


  fprintf(fp,",");
  Substring_print_insertion_2(fp,this->substring1,this->substring2,this->nindels,
			      queryseq,chr,invertp);
  fprintf(fp,"\n");

  if (allocp == true) {
    FREE(chr);
  }

  return;
}

static void
print_deletion (FILE *fp, T this, int score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
		bool invertp, T hit5, T hit3, int insertlength, int pairscore,
		Pairtype_T pairtype, int mapq_score) {
  char *chr;
  bool allocp;

  chr = Univ_IIT_label(chromosome_iit,this->chrnum,&allocp);

  fprintf(fp," ");
  Substring_print_deletion_1(fp,this->substring1,this->substring2,this->nindels,this->deletion,
			     queryseq,chr,invertp);
  /* Alignment info */
  fprintf(fp,"\t");
  print_alignment_info(fp,/*nblocks*/2,score,mapq_score,this->sarrayp);

  /* Pairing info */
  if (hit5 != NULL && hit3 != NULL) {
    fprintf(fp,"\t");
    print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
  }

  fprintf(fp,"\n");

  fprintf(fp,",");
  Substring_print_deletion_2(fp,this->substring1,this->substring2,this->nindels,
			     queryseq,chr,invertp);
  fprintf(fp,"\n");

  if (allocp == true) {
    FREE(chr);
  }
}


static void
print_splice (FILE *fp, T chimera, int score,
	      Univ_IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, T hit5, T hit3,
	      int insertlength, int pairscore, Pairtype_T pairtype, int mapq_score) {
  Substring_T donor, acceptor;
  
  if (chimera->hittype == HALFSPLICE_DONOR) {
    donor = chimera->substring_donor;
    acceptor = (Substring_T) NULL;
    Substring_assign_donor_prob(donor);

  } else if (chimera->hittype == HALFSPLICE_ACCEPTOR) {
    acceptor = chimera->substring_acceptor;
    donor = (Substring_T) NULL;
    Substring_assign_acceptor_prob(acceptor);

  } else {
    donor = chimera->substring_donor;
    acceptor = chimera->substring_acceptor;
    Substring_assign_donor_prob(donor);
    Substring_assign_acceptor_prob(acceptor);
  }

  if (donor == NULL) {
    fprintf(fp," ");
    if (chimera->sensedir == SENSE_FORWARD) {
      Substring_print_acceptor(fp,acceptor,/*sensep*/true,invertp,queryseq,
			       chromosome_iit,donor,chimera->distance);
    } else {
      Substring_print_acceptor(fp,acceptor,/*sensep*/false,invertp,queryseq,
			       chromosome_iit,donor,chimera->distance);
    }

    /* Alignment info */
    fprintf(fp,"\t");
    print_alignment_info(fp,/*nblocks*/1,score,mapq_score,chimera->sarrayp);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      fprintf(fp,"\t");
      print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
    }
    fprintf(fp,"\n");

  } else if (acceptor == NULL) {
    fprintf(fp," ");
    if (chimera->sensedir == SENSE_FORWARD) {
      Substring_print_donor(fp,donor,/*sensep*/true,invertp,
			    queryseq,chromosome_iit,acceptor,chimera->distance);
    } else {
      Substring_print_donor(fp,donor,/*sensep*/false,invertp,
			    queryseq,chromosome_iit,acceptor,chimera->distance);
    }

    /* Alignment info */
    fprintf(fp,"\t");
    print_alignment_info(fp,/*nblocks*/1,score,mapq_score,chimera->sarrayp);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      fprintf(fp,"\t");
      print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
    }
    fprintf(fp,"\n");

  } else if (chimera->sensedir == SENSE_FORWARD && invertp == false) {
    fprintf(fp," ");
    Substring_print_donor(fp,donor,/*sensep*/true,/*invertp*/false,
			  queryseq,chromosome_iit,acceptor,chimera->distance);
    /* Alignment info */
    fprintf(fp,"\t");
    print_alignment_info(fp,/*nblocks*/2,score,mapq_score,chimera->sarrayp);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      fprintf(fp,"\t");
      print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
    }
    fprintf(fp,"\n");

    fprintf(fp,",");
    Substring_print_acceptor(fp,acceptor,/*sensep*/true,/*invertp*/false,queryseq,
			     chromosome_iit,donor,chimera->distance);
    fprintf(fp,"\n");

  } else if (chimera->sensedir == SENSE_FORWARD && invertp == true) {
    fprintf(fp," ");
    Substring_print_acceptor(fp,acceptor,/*sensep*/true,/*invertp*/true,queryseq,
			     chromosome_iit,donor,chimera->distance);
    /* Alignment info */
    fprintf(fp,"\t");
    print_alignment_info(fp,/*nblocks*/2,score,mapq_score,chimera->sarrayp);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      fprintf(fp,"\t");
      print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
    }
    fprintf(fp,"\n");

    fprintf(fp,",");
    Substring_print_donor(fp,donor,/*sensep*/true,/*invertp*/true,queryseq,
			  chromosome_iit,acceptor,chimera->distance);
    fprintf(fp,"\n");

  } else if (chimera->sensedir == SENSE_ANTI && invertp == false) {
    fprintf(fp," ");
    Substring_print_acceptor(fp,acceptor,/*sensep*/false,/*invertp*/false,queryseq,
			     chromosome_iit,donor,chimera->distance);
    /* Alignment info */
    fprintf(fp,"\t");
    print_alignment_info(fp,/*nblocks*/2,score,mapq_score,chimera->sarrayp);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      fprintf(fp,"\t");
      print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
    }
    fprintf(fp,"\n");

    fprintf(fp,",");
    Substring_print_donor(fp,donor,/*sensep*/false,/*invertp*/false,queryseq,
			  chromosome_iit,acceptor,chimera->distance);
    fprintf(fp,"\n");

  } else if (chimera->sensedir == SENSE_ANTI && invertp == true) {
    fprintf(fp," ");
    Substring_print_donor(fp,donor,/*sensep*/false,/*invertp*/true,queryseq,
			  chromosome_iit,acceptor,chimera->distance);
    /* Alignment info */
    fprintf(fp,"\t");
    print_alignment_info(fp,/*nblocks*/2,score,mapq_score,chimera->sarrayp);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      fprintf(fp,"\t");
      print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
    }
    fprintf(fp,"\n");

    fprintf(fp,",");
    Substring_print_acceptor(fp,acceptor,/*sensep*/false,/*invertp*/true,queryseq,
			     chromosome_iit,donor,chimera->distance);
    fprintf(fp,"\n");
  }

  return;
}


static void
print_shortexon (FILE *fp, T chimera, int score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, T hit5, T hit3,
		 int insertlength, int pairscore, Pairtype_T pairtype, int mapq_score) {
  Substring_T donor, acceptor, shortexon;
  Chrpos_T distance1, distance2;
  bool firstp = true;
  int nblocks = 1;
  
  shortexon = chimera->substring1;
  Substring_assign_shortexon_prob(shortexon);
  if ((donor = chimera->substringD) != NULL) {
    Substring_assign_donor_prob(donor);
    nblocks++;
  }
  if ((acceptor = chimera->substringA) != NULL) {
    Substring_assign_acceptor_prob(acceptor);
    nblocks++;
  }


  if (chimera->sensedir == SENSE_FORWARD && invertp == false) {
    distance1 = chimera->acceptor_distance;
    distance2 = chimera->donor_distance;

    if (donor != NULL) {
      fprintf(fp," ");
      Substring_print_donor(fp,donor,/*sensep*/true,/*invertp*/false,
			    queryseq,chromosome_iit,acceptor,distance1);
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score,chimera->sarrayp);
      if (hit5 != NULL && hit3 != NULL) {
	fprintf(fp,"\t"); print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
      }
      firstp = false;
      fprintf(fp,"\n");
    }

    if (firstp == true) { fprintf(fp," "); } else { fprintf(fp,","); }
    Substring_print_shortexon(fp,shortexon,/*sensep*/true,/*invertp*/false,queryseq,
			      chromosome_iit,distance1,distance2);
    if (firstp == true) {
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score,chimera->sarrayp);
      if (hit5 != NULL && hit3 != NULL) {
	fprintf(fp,"\t"); print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
      }
    }
    fprintf(fp,"\n");

    if (acceptor != NULL) {
      fprintf(fp,",");
      Substring_print_acceptor(fp,acceptor,/*sensep*/true,/*invertp*/false,queryseq,
			       chromosome_iit,donor,distance2);
      fprintf(fp,"\n");
    }

  } else if (chimera->sensedir == SENSE_FORWARD && invertp == true) {
    distance1 = chimera->donor_distance;
    distance2 = chimera->acceptor_distance;

    if (acceptor != NULL) {
      fprintf(fp," ");
      Substring_print_acceptor(fp,acceptor,/*sensep*/true,/*invertp*/true,queryseq,
			       chromosome_iit,donor,distance1);
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score,chimera->sarrayp);
      if (hit5 != NULL && hit3 != NULL) {
	fprintf(fp,"\t"); print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
      }
      firstp = false;
      fprintf(fp,"\n");
    }

    if (firstp == true) { fprintf(fp," "); } else { fprintf(fp,","); }
    Substring_print_shortexon(fp,shortexon,/*sensep*/true,/*invertp*/true,queryseq,
			      chromosome_iit,distance1,distance2);
    if (firstp == true) {
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score,chimera->sarrayp);
      if (hit5 != NULL && hit3 != NULL) {
	fprintf(fp,"\t"); print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
      }
    }
    fprintf(fp,"\n");

    if (donor != NULL) {
      fprintf(fp,",");
      Substring_print_donor(fp,donor,/*sensep*/true,/*invertp*/true,queryseq,
			    chromosome_iit,acceptor,distance2);
      fprintf(fp,"\n");
    }

  } else if (chimera->sensedir == SENSE_ANTI && invertp == false) {
    distance1 = chimera->donor_distance;
    distance2 = chimera->acceptor_distance;

    if (acceptor != NULL) {
      fprintf(fp," ");
      Substring_print_acceptor(fp,acceptor,/*sensep*/false,/*invertp*/false,queryseq,
			       chromosome_iit,donor,distance1);
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score,chimera->sarrayp);
      if (hit5 != NULL && hit3 != NULL) {
	fprintf(fp,"\t"); print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
      }
      firstp = false;
      fprintf(fp,"\n");
    }

    if (firstp == true) { fprintf(fp," "); } else { fprintf(fp,","); }
    Substring_print_shortexon(fp,shortexon,/*sensep*/false,/*invertp*/false,queryseq,
			      chromosome_iit,distance1,distance2);
    if (firstp == true) {
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score,chimera->sarrayp);
      if (hit5 != NULL && hit3 != NULL) {
	fprintf(fp,"\t"); print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
      }
    }
    fprintf(fp,"\n");

    if (donor != NULL) {
      fprintf(fp,",");
      Substring_print_donor(fp,donor,/*sensep*/false,/*invertp*/false,queryseq,
			    chromosome_iit,acceptor,distance2);
      fprintf(fp,"\n");
    }

  } else if (chimera->sensedir == SENSE_ANTI && invertp == true) {
    distance2 = chimera->donor_distance;
    distance1 = chimera->acceptor_distance;

    if (donor != NULL) {
      fprintf(fp," ");
      Substring_print_donor(fp,donor,/*sensep*/false,/*invertp*/true,queryseq,
			    chromosome_iit,acceptor,distance1);
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score,chimera->sarrayp);
      if (hit5 != NULL && hit3 != NULL) {
	fprintf(fp,"\t"); print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
      }
      firstp = false;
      fprintf(fp,"\n");
    }

    if (firstp == true) { fprintf(fp," "); } else { fprintf(fp,","); }
    Substring_print_shortexon(fp,shortexon,/*sensep*/false,/*invertp*/true,queryseq,
			      chromosome_iit,distance1,distance2);
    if (firstp == true) {
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score,chimera->sarrayp);
      if (hit5 != NULL && hit3 != NULL) {
	fprintf(fp,"\t"); print_pair_info(fp,hit5,hit3,insertlength,pairscore,pairtype);
      }
    }
    fprintf(fp,"\n");

    if (acceptor != NULL) {
      fprintf(fp,",");
      Substring_print_acceptor(fp,acceptor,/*sensep*/false,/*invertp*/true,queryseq,
			       chromosome_iit,donor,distance2);
      fprintf(fp,"\n");
    }
  }

  return;
}



/* May substitute paired-end loglik for single-end loglik */
void
Stage3end_print (FILE *fp, T this, int score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, T hit5, T hit3, int insertlength,
		 int pairscore, Pairtype_T pairtype, int mapq_score) {

  if (this->hittype == EXACT || this->hittype == SUB || this->hittype == TERMINAL) {
    print_single(fp,this,score,chromosome_iit,queryseq,invertp,
		 hit5,hit3,insertlength,pairscore,pairtype,mapq_score);
  } else if (this->hittype == INSERTION) {
    print_insertion(fp,this,score,chromosome_iit,queryseq,invertp,
		    hit5,hit3,insertlength,pairscore,pairtype,mapq_score);
  } else if (this->hittype == DELETION) {
    print_deletion(fp,this,score,chromosome_iit,queryseq,invertp,
		   hit5,hit3,insertlength,pairscore,pairtype,mapq_score);
  } else if (this->hittype == HALFSPLICE_DONOR || this->hittype == HALFSPLICE_ACCEPTOR ||
	     this->hittype == SPLICE || this->hittype == SAMECHR_SPLICE || this->hittype == TRANSLOC_SPLICE) {
    print_splice(fp,this,score,
		 chromosome_iit,queryseq,invertp,hit5,hit3,insertlength,
		 pairscore,pairtype,mapq_score);
  } else if (this->hittype == ONE_THIRD_SHORTEXON || this->hittype == TWO_THIRDS_SHORTEXON || this->hittype == SHORTEXON) {
    print_shortexon(fp,this,score,
		    chromosome_iit,queryseq,invertp,hit5,hit3,insertlength,
		    pairscore,pairtype,mapq_score);
  } else if (this->hittype == GMAP) {
    if (Shortread_invertedp(queryseq) == false) {
      Substring_print_gmap(fp,this->pairarray,this->npairs,this->nsegments,/*invertedp*/false,
			   this->gmap_start_endtype,this->gmap_end_endtype,
			   this->chrnum,this->chroffset,this->chrhigh,Shortread_fulllength(queryseq),
			   this->plusp,this->gmap_cdna_direction,this->score,insertlength,pairscore,mapq_score,
			   chromosome_iit);
    } else {
      Substring_print_gmap(fp,this->pairarray,this->npairs,this->nsegments,/*invertedp*/true,
			   this->gmap_end_endtype,this->gmap_start_endtype,
			   this->chrnum,this->chroffset,this->chrhigh,Shortread_fulllength(queryseq),
			   this->plusp,this->gmap_cdna_direction,this->score,insertlength,pairscore,mapq_score,
			   chromosome_iit);
    }

  } else {
    abort();
  }

  return;
}


static void
print_query_header (FILE *fp, char initchar, Shortread_T queryseq, bool invertp) {
  fprintf(fp,"%c",initchar);
  if (invertp == false) {
    Shortread_print_oneline(fp,queryseq);
  } else {
    Shortread_print_oneline_revcomp(fp,queryseq);
  }

  return;
}



static void
print_barcode_and_quality (FILE *fp, Shortread_T queryseq, bool invertp, int quality_shift) {
  char *barcode;

  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    fprintf(fp,"\tbarcode:%s",barcode);
  }

  if (Shortread_quality_string(queryseq) != NULL) {
    fprintf(fp,"\t");
    if (invertp == false) {
      Shortread_print_quality(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			      quality_shift,/*show_chopped_p*/true);
    } else {
      Shortread_print_quality_revcomp(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				      quality_shift,/*show_chopped_p*/true);
    }
  }

  return;
}


static void
print_one_paired_end (Result_T result, Resulttype_T resulttype,
		      char initchar, bool firstp, Univ_IIT_T chromosome_iit,
		      Shortread_T queryseq, Shortread_T headerseq1, Shortread_T headerseq2,
		      int maxpaths, bool quiet_if_excessive_p,
		      bool invertp, int quality_shift,
		      FILE *fp_nomapping_1, FILE *fp_unpaired_uniq, FILE *fp_unpaired_circular,
		      FILE *fp_unpaired_transloc, FILE *fp_unpaired_mult,
		      FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_circular,
		      FILE *fp_halfmapping_transloc, FILE *fp_halfmapping_mult,
		      FILE *fp_paired_uniq_circular, FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr,
		      FILE *fp_paired_uniq_long, FILE *fp_paired_mult,
		      FILE *fp_concordant_uniq, FILE *fp_concordant_circular,
		      FILE *fp_concordant_transloc, FILE *fp_concordant_mult) {
  Stage3pair_T *stage3pairarray, stage3pair;
  T *stage3array, *stage3array_mate, this, hit5, hit3;
  int npaths, npaths_mate, pathnum, first_absmq, second_absmq;
  bool outputp, translocationp;
  FILE *fp;

  if (resulttype == PAIREDEND_NOMAPPING) {
    /* If fails_as_input_p == true, then this case is handled by calling procedure */
    print_query_header(fp_nomapping_1,initchar,queryseq,invertp);
    fprintf(fp_nomapping_1,"\t0 %s",UNPAIRED_TEXT);

    print_barcode_and_quality(fp_nomapping_1,queryseq,invertp,quality_shift);
    
    fprintf(fp_nomapping_1,"\t");
    Shortread_print_header(fp_nomapping_1,headerseq1,headerseq2);
    fprintf(fp_nomapping_1,"\n");

  } else if (resulttype == CONCORDANT_UNIQ) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
    stage3pair = stage3pairarray[0];
    hit5 = stage3pair->hit5;
    hit3 = stage3pair->hit3;

    if (stage3pair->circularp == true) {
      fp = fp_concordant_circular;
    } else {
      fp = fp_concordant_uniq;
    }

    print_query_header(fp,initchar,queryseq,invertp);
    fprintf(fp,"\t1 %s",CONCORDANT_TEXT);
    
    print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

    fprintf(fp,"\t");
    Shortread_print_header(fp,headerseq1,headerseq2);

    
    if (firstp == true) {
#if 0
      Stage3pair_eval(stage3pairarray,/*npaths*/1,maxpaths,queryseq,queryseq_mate);
#endif
      Stage3end_print(fp,hit5,hit5->score,
		      chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		      stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    } else {
      Stage3end_print(fp,hit3,hit3->score,
		      chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		      stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    }

    fprintf(fp,"\n");

  } else if (resulttype == CONCORDANT_TRANSLOC) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    print_query_header(fp_concordant_transloc,initchar,queryseq,invertp);
    fprintf(fp_concordant_transloc,"\t%d %s",npaths,CONCORDANT_TEXT);
    fprintf(fp_concordant_transloc," (transloc)");

    print_barcode_and_quality(fp_concordant_transloc,queryseq,invertp,quality_shift);

    fprintf(fp_concordant_transloc,"\t");
    Shortread_print_header(fp_concordant_transloc,headerseq1,headerseq2);

    if (quiet_if_excessive_p && npaths > maxpaths) {
      /* No further output */
    } else {
#if 0
      if (firstp == true) {
	Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq,queryseq_mate);
      }
#endif

      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;

	if (firstp == true) {
	  Stage3end_print(fp_concordant_transloc,hit5,hit5->score,
			  chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	} else {
	  Stage3end_print(fp_concordant_transloc,hit3,hit3->score,
			  chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	}
      }
    }

    fprintf(fp_concordant_transloc,"\n");

  } else if (resulttype == CONCORDANT_MULT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    print_query_header(fp_concordant_mult,initchar,queryseq,invertp);
    fprintf(fp_concordant_mult,"\t%d %s",npaths,CONCORDANT_TEXT);

    print_barcode_and_quality(fp_concordant_mult,queryseq,invertp,quality_shift);

    fprintf(fp_concordant_mult,"\t");
    Shortread_print_header(fp_concordant_mult,headerseq1,headerseq2);

    if (quiet_if_excessive_p && npaths > maxpaths) {
      /* No further output */
    } else {
#if 0
      if (firstp == true) {
	Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq,queryseq_mate);
      }
#endif

      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;

	if (firstp == true) {
	  Stage3end_print(fp_concordant_mult,hit5,hit5->score,
			  chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	} else {
	  Stage3end_print(fp_concordant_mult,hit3,hit3->score,
			  chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	}
      }
    }

    fprintf(fp_concordant_mult,"\n");

  } else if (resulttype == PAIRED_UNIQ) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
    stage3pair = stage3pairarray[0];

    if (stage3pair->circularp == true) {
      fp = fp_paired_uniq_circular;
    } else if (stage3pair->pairtype == PAIRED_INVERSION) {
      fp = fp_paired_uniq_inv;
    } else if (stage3pair->pairtype == PAIRED_SCRAMBLE) {
      fp = fp_paired_uniq_scr;
    } else if (stage3pair->pairtype == PAIRED_TOOLONG) {
      fp = fp_paired_uniq_long;
    } else {
      fprintf(stderr,"Unexpected pairtype %d\n",stage3pair->pairtype);
      abort();
    }
    
    print_query_header(fp,initchar,queryseq,invertp);
    fprintf(fp,"\t1 %s",PAIRED_TEXT);

    print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

    fprintf(fp,"\t");
    Shortread_print_header(fp,headerseq1,headerseq2);

    hit5 = stage3pair->hit5;
    hit3 = stage3pair->hit3;
    
    if (firstp == true) {
#if 0
      Stage3pair_eval(stage3pairarray,/*npaths*/1,maxpaths,queryseq,queryseq_mate);
#endif
      Stage3end_print(fp,hit5,hit5->score,
		      chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		      stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    } else {
      Stage3end_print(fp,hit3,hit3->score,
		      chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		      stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    }

    fprintf(fp,"\n");

  } else if (resulttype == PAIRED_MULT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    print_query_header(fp_paired_mult,initchar,queryseq,invertp);
    fprintf(fp_paired_mult,"\t%d %s",npaths,PAIRED_TEXT);

    print_barcode_and_quality(fp_paired_mult,queryseq,invertp,quality_shift);

    fprintf(fp_paired_mult,"\t");
    Shortread_print_header(fp_paired_mult,headerseq1,headerseq2);

    if (quiet_if_excessive_p && npaths > maxpaths) {
      /* No further output */
    } else {
#if 0
      if (firstp == true) {
	Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq,queryseq_mate);
      }
#endif

      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;

	if (firstp == true) {
	  Stage3end_print(fp_paired_mult,hit5,hit5->score,
			  chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	} else {
	  Stage3end_print(fp_paired_mult,hit3,hit3->score,
			  chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
			  stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	}
      }
    }

    fprintf(fp_paired_mult,"\n");

  } else {
    /* Print as singles */
    if (firstp == true) {
      /* Get stage3array_mate first to avoid incorrect values for npaths */
      stage3array_mate = (T *) Result_array2(&npaths_mate,&first_absmq,&second_absmq,result);
      stage3array = (T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
    } else {
      /* Get stage3array_mate first to avoid incorrect values for npaths */
      stage3array_mate = (T *) Result_array(&npaths_mate,&first_absmq,&second_absmq,result);
      stage3array = (T *) Result_array2(&npaths,&first_absmq,&second_absmq,result);
    }

    outputp = true;
    translocationp = false;
    if (resulttype == HALFMAPPING_UNIQ) {
      if (npaths > 0 && Stage3end_circularpos(stage3array[0]) > 0) {
	fp = fp_halfmapping_circular;
      } else if (npaths_mate > 0 && Stage3end_circularpos(stage3array_mate[0]) > 0) {
	fp = fp_halfmapping_circular;
      } else {
	fp = fp_halfmapping_uniq;
      }

    } else if (resulttype == HALFMAPPING_TRANSLOC) {
      fp = fp_halfmapping_transloc;
      translocationp = true;

    } else if (resulttype == HALFMAPPING_MULT) {
      fp = fp_halfmapping_mult;
      if (quiet_if_excessive_p && npaths > maxpaths) {
	outputp = false;
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      if (npaths > 0 && Stage3end_circularpos(stage3array[0]) > 0) {
	fp = fp_unpaired_circular;
      } else if (npaths_mate > 0 && Stage3end_circularpos(stage3array_mate[0]) > 0) {
	fp = fp_unpaired_circular;
      } else {
	fp = fp_unpaired_uniq;
      }

    } else if (resulttype == UNPAIRED_TRANSLOC) {
      fp = fp_unpaired_transloc;
      translocationp = true;

    } else if (resulttype == UNPAIRED_MULT) {
      fp = fp_unpaired_mult;
      if (quiet_if_excessive_p && npaths > maxpaths) {
	outputp = false;
      }

    } else {
      fprintf(stderr,"Resulttype is %s\n",Resulttype_string(resulttype));
      abort();
    }

    print_query_header(fp,initchar,queryseq,invertp);
    fprintf(fp,"\t%d %s",npaths,UNPAIRED_TEXT);
    if (translocationp == true) {
      fprintf(fp," (transloc)");
    }

#if 0
    /* Print unpaired type for unpaired_uniq results */
    if (resulttype == UNPAIRED_UNIQ) {
      stage3array = (T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
      hit5 = stage3array[0];
      stage3array = (T *) Result_array2(&npaths,&first_absmq,&second_absmq,result);
      hit3 = stage3array[0];
      fprintf(fp," (%s)",unpaired_type_text(hit5,hit3));
    }
#endif

    print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

    fprintf(fp,"\t");
    Shortread_print_header(fp,headerseq1,headerseq2);

    if (outputp == true) {
      /* Stage3end_eval_and_sort(stage3array,npaths,maxpaths,queryseq); */
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	this = stage3array[pathnum-1];
	Stage3end_print(fp,this,this->score,
			chromosome_iit,queryseq,invertp,/*hit5*/(T) NULL,/*hit3*/(T) NULL,
			/*insertlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,this->mapq_score);
      }
    }

    fprintf(fp,"\n");
  }


  return;
}


/* Gets invert_first_p and invert_second_p from global above */

void
Stage3pair_print (Result_T result, Resulttype_T resulttype,
		  Univ_IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  int maxpaths, bool quiet_if_excessive_p,
		  bool nofailsp, bool failsonlyp, bool fails_as_input_p,
		  bool fastq_format_p, int quality_shift,
		  FILE *fp_nomapping_1, FILE *fp_nomapping_2,
		  FILE *fp_unpaired_uniq, FILE *fp_unpaired_circular,
		  FILE *fp_unpaired_transloc, FILE *fp_unpaired_mult,
		  FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_circular,
		  FILE *fp_halfmapping_transloc, FILE *fp_halfmapping_mult,
		  FILE *fp_paired_uniq_circular, FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr,
		  FILE *fp_paired_uniq_long, FILE *fp_paired_mult,
		  FILE *fp_concordant_uniq, FILE *fp_concordant_circular,
		  FILE *fp_concordant_transloc, FILE *fp_concordant_mult) {

  debug1(printf("Stage3pair_print: resulttype is %s\n",Resulttype_string(resulttype)));

  if (resulttype == PAIREDEND_NOMAPPING) {
    if (nofailsp == true) {
      /* No output */
      debug1(printf("  nofailsp is true, so no output\n"));
      return;

    } else if (fails_as_input_p == true) {
      debug1(printf("  fails as input is true, so printing\n"));
      if (fastq_format_p == true) {
	Shortread_print_query_pairedend_fastq(fp_nomapping_1,fp_nomapping_2,queryseq1,queryseq2,
					      invert_first_p,invert_second_p);
      } else {
	Shortread_print_query_pairedend_fasta(fp_nomapping_1,queryseq1,queryseq2,
					      invert_first_p,invert_second_p);
      }
      return;

    } else {
      debug1(printf("  printing failure output\n"));

      /* First end */
      print_one_paired_end(result,resulttype,'>',/*firstp*/true,chromosome_iit,
			   /*queryseq*/queryseq1,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
			   maxpaths,quiet_if_excessive_p,invert_first_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_circular,
			   fp_unpaired_transloc,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_circular,
			   fp_halfmapping_transloc,fp_halfmapping_mult,
			   fp_paired_uniq_circular,fp_paired_uniq_inv,fp_paired_uniq_scr,
			   fp_paired_uniq_long,fp_paired_mult,
			   fp_concordant_uniq,fp_concordant_circular,
			   fp_concordant_transloc,fp_concordant_mult);

      /* Second end */
      print_one_paired_end(result,resulttype,'<',/*firstp*/false,chromosome_iit,
			   /*queryseq*/queryseq2,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
			   maxpaths,quiet_if_excessive_p,invert_second_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_circular,
			   fp_unpaired_transloc,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_circular,
			   fp_halfmapping_transloc,fp_halfmapping_mult,
			   fp_paired_uniq_circular,fp_paired_uniq_inv,fp_paired_uniq_scr,
			   fp_paired_uniq_long,fp_paired_mult,
			   fp_concordant_uniq,fp_concordant_circular,
			   fp_concordant_transloc,fp_concordant_mult);
    }

  } else {
    if (failsonlyp == true) {
      /* Unwanted success: skip */
      debug1(printf("  failsonlyp is true, so no output\n"));
    
    } else {
      /* First end */
      print_one_paired_end(result,resulttype,'>',/*firstp*/true,chromosome_iit,
			   /*queryseq*/queryseq1,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
			   maxpaths,quiet_if_excessive_p,invert_first_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_circular,
			   fp_unpaired_transloc,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_circular,
			   fp_halfmapping_transloc,fp_halfmapping_mult,
			   fp_paired_uniq_circular,fp_paired_uniq_inv,fp_paired_uniq_scr,
			   fp_paired_uniq_long,fp_paired_mult,
			   fp_concordant_uniq,fp_concordant_circular,
			   fp_concordant_transloc,fp_concordant_mult);

      /* Second end */
      print_one_paired_end(result,resulttype,'<',/*firstp*/false,chromosome_iit,
			   /*queryseq*/queryseq2,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
			   maxpaths,quiet_if_excessive_p,invert_second_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_circular,
			   fp_unpaired_transloc,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_circular,
			   fp_halfmapping_transloc,fp_halfmapping_mult,
			   fp_paired_uniq_circular,fp_paired_uniq_inv,fp_paired_uniq_scr,
			   fp_paired_uniq_long,fp_paired_mult,
			   fp_concordant_uniq,fp_concordant_circular,
			   fp_concordant_transloc,fp_concordant_mult);
    }
  }
    
  return;
}


#if 0
List_T
Stage3end_filter_bymatch (List_T hitlist) {
  List_T filtered = NULL, p;
  T hit;
  int min_nmismatches_whole = 1000;

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->nmismatches_whole < min_nmismatches_whole) {
      min_nmismatches_whole = hit->nmismatches_whole;
    }
  }

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->nmismatches_whole == min_nmismatches_whole) {
      filtered = List_push(filtered,hit);
    } else {
      Stage3end_free(&hit);
    }
  }
  List_free(&hitlist);

  return filtered;
}
#endif



static Chrpos_T
overlap5_gmap_plus (int *querypos, Chrpos_T *genomicstart, Chrpos_T *genomicend,
		    Stage3end_T hit5, Stage3end_T gmap) {
  debug10(printf("Entered overlap5_gmap_plus\n"));
  *genomicstart = Substring_chrstart(hit5->substring_high);
  *genomicend = Substring_chrend(hit5->substring_high);
  return Pair_binary_search_ascending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
				      *genomicstart,*genomicend);
}


static Chrpos_T
overlap3_gmap_plus (int *querypos, Chrpos_T *genomicstart, Chrpos_T *genomicend,
		    Stage3end_T hit3, Stage3end_T gmap) {
  debug10(printf("Entered overlap3_gmap_plus\n"));
  *genomicstart = Substring_chrstart(hit3->substring_low);
  *genomicend = Substring_chrend(hit3->substring_low);
  return Pair_binary_search_ascending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
				      *genomicstart,*genomicend);
}


static Chrpos_T
overlap5_gmap_minus (int *querypos, Chrpos_T *genomicstart, Chrpos_T *genomicend,
		     Stage3end_T hit5, Stage3end_T gmap) {
  debug10(printf("Entered overlap5_gmap_minus\n"));
  *genomicstart = Substring_chrstart(hit5->substring_low);
  *genomicend = Substring_chrend(hit5->substring_low);
  return Pair_binary_search_descending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
				       *genomicstart,*genomicend);
}

static Chrpos_T
overlap3_gmap_minus (int *querypos, Chrpos_T *genomicstart, Chrpos_T *genomicend,
		     Stage3end_T hit3, Stage3end_T gmap) {
  debug10(printf("Entered overlap3_gmap_minus\n"));
  *genomicstart = Substring_chrstart(hit3->substring_high);
  *genomicend = Substring_chrend(hit3->substring_high);
  return Pair_binary_search_descending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
				       *genomicstart,*genomicend);
}


static void
resolve_inside_ambiguous_splice_plus (int *unresolved_amb_nmatches, T *hit5, T *hit3, bool *private5p, bool *private3p,
				      Univcoord_T *splicesites,
				      Compress_T query5_compress_fwd, Compress_T query3_compress_fwd,
				      int localsplicing_penalty, int querylength5, int querylength3,
				      int genestrand) {
#ifdef USE_BINGO
  Chrpos_T insertlength;
#endif
  Univcoord_T genomicstart, genomicend;
  int nbingo, bingoi5, bingoi3, nbounded, boundedi5, boundedi3, nbest, besti5, besti3, i, j;
  int best_nmismatches, nmismatches;
  bool new5p = false, new3p = false;
  T old;

  Substring_T donor, acceptor, shortexon;
  Univcoord_T segment_left;
  int nmismatches_shortend;
  int donor_knowni, acceptor_knowni;
  int splice_pos;
  int ignore_found_score = 0;

  *unresolved_amb_nmatches = 0;

  debug9(printf("resolve plus: hit5 %s ambiguous %d,%d and hit3 %s ambiguous %d,%d\n",
		hittype_string((*hit5)->hittype),(*hit5)->start_ambiguous_p,(*hit5)->end_ambiguous_p,
		hittype_string((*hit3)->hittype),(*hit3)->start_ambiguous_p,(*hit3)->end_ambiguous_p));

  if ((*hit5)->end_ambiguous_p == true && (*hit3)->start_ambiguous_p == true) {
    debug9(printf("Got ambiguous at 5' and ambiguous at 3':"));
    nbest = nbounded = nbingo = 0;
    best_nmismatches = querylength5 + querylength3;
    for (i = 0; i < (*hit5)->end_nambi; i++) {
      genomicend = splicesites[(*hit5)->end_ambi[i]];
      for (j = 0; j < (*hit3)->start_nambi; j++) {
	genomicstart = splicesites[(*hit3)->start_ambi[j]];
	debug9(printf(" %lu,%lu",genomicend,genomicstart));
	if (genomicend < genomicstart) {
	  nbounded++;
	  boundedi5 = i;
	  boundedi3 = j;

#ifdef USE_BINGO
	  insertlength = genomicstart - genomicend + querylength5 + querylength3;
	  debug9(printf(" (%u)",insertlength));
	  if (insertlength < expected_pairlength) {
	    if (expected_pairlength - insertlength <= pairlength_deviation) {
	      nbingo++;
	      bingoi5 = i;
	      bingoi3 = j;
	      debug9(printf("*"));
	    }
	  } else {
	    if (insertlength - expected_pairlength <= pairlength_deviation) {
	      nbingo++;
	      bingoi5 = i;
	      bingoi3 = j;
	      debug9(printf("*"));
	    }
	  }
#endif

	  if ((nmismatches = (*hit5)->end_amb_nmismatches[i] + (*hit3)->start_amb_nmismatches[j]) < best_nmismatches) {
	    best_nmismatches = nmismatches;
	    besti5 = i;
	    besti3 = j;
	    nbest = 1;
	  } else if (nmismatches == best_nmismatches) {
	    nbest++;
	  }
	}
      }
    }

#if 0
    /* No longer holds for GMAP */
    assert((*hit5)->amb_nmatches_end > 0);
    assert((*hit3)->amb_nmatches_start > 0);
#endif

#ifdef USE_BINGO
    if (nbingo == 1) {
      new5p = true; new3p = true;
    } else if (nbounded == 1) {
      new5p = true; new3p = true; bingoi5 = boundedi5; bingoi3 = boundedi3;
    }
#endif

    if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_nmatches = %d...%d",
		    (*hit5)->amb_nmatches_end,(*hit3)->amb_nmatches_start));
      *unresolved_amb_nmatches = (*hit5)->amb_nmatches_end + (*hit3)->amb_nmatches_start;
    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      new5p = true; new3p = true; bingoi5 = besti5; bingoi3 = besti3;
    }
    debug9(printf("\n"));

  } else if ((*hit5)->end_ambiguous_p == true) {
    debug9(printf("Got ambiguous at 5' (%s):",hittype_string((*hit5)->hittype)));
    nbest = nbounded = nbingo = 0;
    best_nmismatches = querylength5;
    for (i = 0; i < (*hit5)->end_nambi; i++) {
      genomicend = splicesites[(*hit5)->end_ambi[i]];
      debug9(printf(" %lu",genomicend));
      if (genomicend < (*hit3)->genomicstart /*allow overlap*/+ querylength3) {
	nbounded++;
	boundedi5 = i;

#ifdef USE_BINGO
	insertlength = (*hit3)->genomicstart - genomicend + querylength5 + querylength3;
	debug9(printf(" (%u)",insertlength));
	if (insertlength < expected_pairlength) {
	  if (expected_pairlength - insertlength <= pairlength_deviation) {
	    nbingo++;
	    bingoi5 = i;
	    debug9(printf("*"));
	  }
	} else {
	  if (insertlength - expected_pairlength <= pairlength_deviation) {
	    nbingo++;
	    bingoi5 = i;
	    debug9(printf("*"));
	  }
	}
#endif

	if ((nmismatches = (*hit5)->end_amb_nmismatches[i]) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	  besti5 = i;
	  nbest = 1;
	} else if (nmismatches == best_nmismatches) {
	  nbest++;
	}
      }
    }

#if 0
    /* No longer holds for GMAP */
    assert((*hit5)->amb_nmatches_end > 0);
    assert((*hit3)->amb_nmatches_start == 0);
#endif

#ifdef USE_BINGO
    if (nbingo == 1) {
      new5p = true;
    } else if (nbounded == 1) {
      new5p = true; bingoi5 = boundedi5;
    }
#endif

    if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_nmatches = %d...%d",
		    (*hit5)->amb_nmatches_end,(*hit3)->amb_nmatches_start));
      *unresolved_amb_nmatches = (*hit5)->amb_nmatches_end;
    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      new5p = true; bingoi5 = besti5;
    }
    debug9(printf("\n"));

  } else if ((*hit3)->start_ambiguous_p == true) {
    debug9(printf("Got ambiguous at 3':"));
    nbest = nbounded = nbingo = 0;
    best_nmismatches = querylength3;
    for (j = 0; j < (*hit3)->start_nambi; j++) {
      genomicstart = splicesites[(*hit3)->start_ambi[j]];
      debug9(printf(" %lu",genomicstart));
      if ((*hit5)->genomicend < genomicstart /*allow overlap*/+ querylength5) {
	nbounded++;
	boundedi3 = j;

#ifdef USE_BINGO
	insertlength = genomicstart - (*hit5)->genomicend + querylength5 + querylength3;
	debug9(printf(" (%u)",insertlength));
	if (insertlength < expected_pairlength) {
	  if (expected_pairlength - insertlength <= pairlength_deviation) {
	    nbingo++;
	    bingoi3 = j;
	    debug9(printf("*"));
	  }
	} else {
	  if (insertlength - expected_pairlength <= pairlength_deviation) {
	    nbingo++;
	    bingoi3 = j;
	    debug9(printf("*"));
	  }
	}
#endif

	if ((nmismatches = (*hit3)->start_amb_nmismatches[j]) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	  besti3 = j;
	  nbest = 1;
	} else if (nmismatches == best_nmismatches) {
	  nbest++;
	}
      }
    }

#if 0
    /* No longer holds for GMAP */
    assert((*hit5)->amb_nmatches_end == 0);
    assert((*hit3)->amb_nmatches_start > 0);
#endif

#ifdef USE_BINGO
    if (nbingo == 1) {
      new3p = true;
    } else if (nbounded == 1) {
      new3p = true; bingoi3 = boundedi3;
    }
#endif

    if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_nmatches = %d...%d",
		    (*hit5)->amb_nmatches_end,(*hit3)->amb_nmatches_start));
      *unresolved_amb_nmatches = (*hit3)->amb_nmatches_start;
    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      new3p = true; bingoi3 = besti3;
    }
    debug9(printf("\n"));
  }

  if (new5p == false) {
    /* Skip */
  } else if ((*hit5)->hittype == ONE_THIRD_SHORTEXON || (*hit5)->hittype == TWO_THIRDS_SHORTEXON) {

    if ((*hit5)->sensedir == SENSE_FORWARD) {
      /* End 1 */
      shortexon = (*hit5)->substring1;
	
      donor_knowni = Substring_splicesites_i_D(shortexon);
      splice_pos = Substring_chimera_pos_D(shortexon);
      acceptor_knowni = (*hit5)->end_ambi[bingoi5];
      nmismatches_shortend = (*hit5)->end_amb_nmismatches[bingoi5];
      segment_left = splicesites[acceptor_knowni] - splice_pos;
	
      if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
					     /*prob*/2.0,/*left*/segment_left,query5_compress_fwd,
					     querylength5,/*plusp*/true,genestrand,/*sensep*/true,
					     Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					     Substring_chrhigh(shortexon),Substring_chrlength(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 1: Splice from donor #%d to acceptor #%d, with nmismatches %d\n",
		      donor_knowni,acceptor_knowni,nmismatches_shortend));
	old = *hit5;
	*hit5 = Stage3end_new_shortexon(&ignore_found_score,/*donor*/(*hit5)->substringD,acceptor,shortexon,
					/*acceptor_distance*/(*hit5)->acceptor_distance,
					/*donor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
					(*hit5)->amb_nmatches_donor,/*amb_nmatches_acceptor*/0,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength5,/*sensedir*/SENSE_FORWARD,
					/*sarrayp*/false);
	if (*private5p == true) {
	  Stage3end_free(&old);
	}
	*private5p = true;
      }

    } else if ((*hit5)->sensedir == SENSE_ANTI) {
      /* End 6 */
      shortexon = (*hit5)->substring1;

      acceptor_knowni = Substring_splicesites_i_A(shortexon);
      splice_pos = Substring_chimera_pos_A(shortexon);
      donor_knowni = (*hit5)->end_ambi[bingoi5];
      nmismatches_shortend = (*hit5)->end_amb_nmismatches[bingoi5];
      segment_left = splicesites[donor_knowni] - splice_pos;

      if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
				       /*prob*/2.0,/*left*/segment_left,query5_compress_fwd,
				       querylength5,/*plusp*/true,genestrand,/*sensep*/false,
				       Substring_chrnum(shortexon),Substring_chroffset(shortexon),
				       Substring_chrhigh(shortexon),Substring_chrlength(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 6: Splice from antiacceptor #%d to antidonor #%d, with nmismatches %d\n",
		      acceptor_knowni,donor_knowni,nmismatches_shortend));
	old = *hit5;
	*hit5 = Stage3end_new_shortexon(&ignore_found_score,donor,/*acceptor*/(*hit5)->substringA,shortexon,
					/*acceptor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
					/*donor_distance*/(*hit5)->donor_distance,
					/*amb_nmatches_donor*/0,(*hit5)->amb_nmatches_acceptor,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength5,/*sensedir*/SENSE_ANTI,
					/*sarrayp*/false);
	if (*private5p == true) {
	  Stage3end_free(&old);
	}
	*private5p = true;
      }

    } else {
      fprintf(stderr,"Shortexon hit5 has no sensedir\n");
      abort();
    }


  } else if ((*hit5)->hittype == HALFSPLICE_DONOR) {
    /* End 1 */
    assert((*hit5)->sensedir == SENSE_FORWARD);
    donor = (*hit5)->substring_donor;

    donor_knowni = Substring_splicesites_i(donor);
    splice_pos = Substring_chimera_pos(donor);
    acceptor_knowni = (*hit5)->end_ambi[bingoi5];
    nmismatches_shortend = (*hit5)->end_amb_nmismatches[bingoi5];
    segment_left = splicesites[acceptor_knowni] - splice_pos;

    if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
					   /*prob*/2.0,/*left*/segment_left,query5_compress_fwd,
					   querylength5,/*plusp*/true,genestrand,/*sensep*/true,
					   Substring_chrnum(donor),Substring_chroffset(donor),
					   Substring_chrhigh(donor),Substring_chrlength(donor))) != NULL) {
      debug9(printf("Resolved halfsplice_donor, End 1: Splice from donor #%d to acceptor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(donor),Substring_splicesites_i(acceptor),nmismatches_shortend));
      old = *hit5;
      *hit5 = Stage3end_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				   donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength5,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,/*sensedir*/SENSE_FORWARD,
				   /*sarrayp*/false);
      if (*private5p == true) {
	Stage3end_free(&old);
      }
      *private5p = true;
    }
	  
  } else if ((*hit5)->hittype == HALFSPLICE_ACCEPTOR) {
    /* End 6 */
    assert((*hit5)->sensedir == SENSE_ANTI);
    acceptor = (*hit5)->substring_acceptor;

    acceptor_knowni = Substring_splicesites_i(acceptor);
    splice_pos = Substring_chimera_pos(acceptor);
    donor_knowni = (*hit5)->end_ambi[bingoi5];
    nmismatches_shortend = (*hit5)->end_amb_nmismatches[bingoi5];
    segment_left = splicesites[donor_knowni] - splice_pos;

    if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
				     /*prob*/2.0,/*left*/segment_left,query5_compress_fwd,
				     querylength5,/*plusp*/true,genestrand,/*sensep*/false,
				     Substring_chrnum(acceptor),Substring_chroffset(acceptor),
				     Substring_chrhigh(acceptor),Substring_chrlength(acceptor))) != NULL) {
      debug9(printf("Resolved halfsplice_acceptor, End 6: Splice from antiacceptor #%d to antidonor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(acceptor),Substring_splicesites_i(donor),nmismatches_shortend));
      old = *hit5;
      *hit5 = Stage3end_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				   donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength5,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,/*sensedir*/SENSE_ANTI,
				   /*sarrayp*/false);
      if (*private5p == true) {
	Stage3end_free(&old);
      }
      *private5p = true;
    }

  } else {
    fprintf(stderr,"Unexpected hittype %d for ambiguous end\n",(*hit5)->hittype);
    abort();
  }

  if (new3p == false) {
    /* Skip */

  } else if ((*hit3)->hittype == ONE_THIRD_SHORTEXON || (*hit3)->hittype == TWO_THIRDS_SHORTEXON) {

    if ((*hit3)->sensedir == SENSE_ANTI) {
      /* End 5 */
      shortexon = (*hit3)->substring1;

      donor_knowni = Substring_splicesites_i_D(shortexon);
      splice_pos = Substring_chimera_pos_D(shortexon);
      acceptor_knowni = (*hit3)->start_ambi[bingoi3];
      nmismatches_shortend = (*hit3)->start_amb_nmismatches[bingoi3];
      segment_left = splicesites[acceptor_knowni] - splice_pos;

      if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
					     /*prob*/2.0,segment_left,query3_compress_fwd,
					     querylength3,/*plusp*/true,genestrand,/*sensep*/false,
					     Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					     Substring_chrhigh(shortexon),Substring_chrlength(shortexon))) != NULL) {
	debug9(printf("Resolved shortexonr, End 5: Splice from antidonor #%d to antiacceptor #%d, with nmismatches %d\n",
		      donor_knowni,acceptor_knowni,nmismatches_shortend));
	old = *hit3;
	*hit3 = Stage3end_new_shortexon(&ignore_found_score,/*donor*/(*hit3)->substringD,acceptor,shortexon,
					/*acceptor_distance*/(*hit3)->acceptor_distance,
					/*donor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
					(*hit3)->amb_nmatches_donor,/*amb_nmatches_acceptor*/0,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength3,/*sensedir*/SENSE_ANTI,
					/*sarrayp*/false);
	if (*private3p == true) {
	  Stage3end_free(&old);
	}
	*private3p = true;
      }

    } else if ((*hit3)->sensedir == SENSE_FORWARD) {
      /* End 2 */
      shortexon = (*hit3)->substring1;

      acceptor_knowni = Substring_splicesites_i_A(shortexon);
      splice_pos = Substring_chimera_pos_A(shortexon);
      donor_knowni = (*hit3)->start_ambi[bingoi3];
      nmismatches_shortend = (*hit3)->start_amb_nmismatches[bingoi3];
      segment_left = splicesites[donor_knowni] - splice_pos;

      if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
				       /*prob*/2.0,segment_left,query3_compress_fwd,
				       querylength3,/*plusp*/true,genestrand,/*sensep*/true,
				       Substring_chrnum(shortexon),Substring_chroffset(shortexon),
				       Substring_chrhigh(shortexon),Substring_chrlength(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 2: Splice from acceptor #%d to donor #%d, with nmismatches %d\n",
		      acceptor_knowni,donor_knowni,nmismatches_shortend));
	old = *hit3;
	*hit3 = Stage3end_new_shortexon(&ignore_found_score,donor,/*acceptor*/(*hit3)->substringA,shortexon,
					/*acceptor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
					/*donor_distance*/(*hit3)->donor_distance,
					/*amb_nmatches_donor*/0,(*hit3)->amb_nmatches_acceptor,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength3,/*sensedir*/SENSE_FORWARD,
					/*sarrayp*/false);
	if (*private3p == true) {
	  Stage3end_free(&old);
	}
	*private3p = true;
      }

    } else {
      fprintf(stderr,"Shortexon hit5 has no sensedir\n");
      abort();
    }


  } else if ((*hit3)->hittype == HALFSPLICE_DONOR) {
    /* End 5 */
    assert((*hit3)->sensedir == SENSE_ANTI);
    donor = (*hit3)->substring_donor;

    donor_knowni = Substring_splicesites_i(donor);
    splice_pos = Substring_chimera_pos(donor);
    acceptor_knowni = (*hit3)->start_ambi[bingoi3];
    nmismatches_shortend = (*hit3)->start_amb_nmismatches[bingoi3];
    segment_left = splicesites[acceptor_knowni] - splice_pos;

    if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
					   /*prob*/2.0,segment_left,query3_compress_fwd,
					   querylength3,/*plusp*/true,genestrand,/*sensep*/false,
					   Substring_chrnum(donor),Substring_chroffset(donor),
					   Substring_chrhigh(donor),Substring_chrlength(donor))) != NULL) {
      debug9(printf("Resolved halfsplice donor, End 5: Splice from antidonor #%d to antiacceptor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(donor),Substring_splicesites_i(acceptor),nmismatches_shortend));
      old = *hit3;
      *hit3 = Stage3end_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				   donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength3,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				   /*sensedir*/SENSE_ANTI,/*sarrayp*/false);
      if (*private3p == true) {
	Stage3end_free(&old);
      }
      *private3p = true;
    }

  } else if ((*hit3)->hittype == HALFSPLICE_ACCEPTOR) {
    /* End 2 */
    assert((*hit3)->sensedir == SENSE_FORWARD);
    acceptor = (*hit3)->substring_acceptor;

    acceptor_knowni = Substring_splicesites_i(acceptor);
    splice_pos = Substring_chimera_pos(acceptor);
    donor_knowni = (*hit3)->start_ambi[bingoi3];
    nmismatches_shortend = (*hit3)->start_amb_nmismatches[bingoi3];
    segment_left = splicesites[donor_knowni] - splice_pos;

    if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
				     /*prob*/2.0,segment_left,query3_compress_fwd,
				     querylength3,/*plusp*/true,genestrand,/*sensep*/true,
				     Substring_chrnum(acceptor),Substring_chroffset(acceptor),
				     Substring_chrhigh(acceptor),Substring_chrlength(acceptor))) != NULL) {
      debug9(printf("Resolved halfsplice acceptor, End 2: Splice from acceptor #%d (%lu) to donor #%d (%lu), with nmismatches %d\n",
		    Substring_splicesites_i(acceptor),splicesites[acceptor_knowni],
		    Substring_splicesites_i(donor),splicesites[donor_knowni],nmismatches_shortend));
      old = *hit3;
      *hit3 = Stage3end_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				   donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength3,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				   /*sensedir*/SENSE_FORWARD,/*sarrayp*/false);
      if (*private3p == true) {
	Stage3end_free(&old);
      }
      *private3p = true;
    }

  } else {
    fprintf(stderr,"Unexpected hittype %d for ambiguous end\n",(*hit3)->hittype);
    abort();
  }
      
  return;
}


static void
resolve_inside_ambiguous_splice_minus (int *unresolved_amb_nmatches, T *hit5, T *hit3, bool *private5p, bool *private3p,
				       Univcoord_T *splicesites,
				       Compress_T query5_compress_rev, Compress_T query3_compress_rev,
				       int localsplicing_penalty, int querylength5, int querylength3,
				       int genestrand) {
#ifdef USE_BINGO
  int insertlength;
#endif
  Univcoord_T genomicstart, genomicend;
  int nbingo, bingoi5, bingoi3, nbounded, boundedi5, boundedi3, nbest, besti5, besti3, i, j;
  int best_nmismatches, nmismatches;
  bool new5p = false, new3p = false;
  T old;

  Substring_T donor, acceptor, shortexon;
  Univcoord_T segment_left;
  int nmismatches_shortend;
  int donor_knowni, acceptor_knowni;
  int splice_pos;
  int ignore_found_score = 0;

  *unresolved_amb_nmatches = 0;

  debug9(printf("resolve minus: hit5 %s ambiguous %d,%d and hit3 %s ambiguous %d,%d\n",
		hittype_string((*hit5)->hittype),(*hit5)->start_ambiguous_p,(*hit5)->end_ambiguous_p,
		hittype_string((*hit3)->hittype),(*hit3)->start_ambiguous_p,(*hit3)->end_ambiguous_p));

  if ((*hit5)->end_ambiguous_p == true && (*hit3)->start_ambiguous_p == true) {
    debug9(printf("Got ambiguous at 5' and ambiguous at 3':"));
    nbest = nbounded = nbingo = 0;
    best_nmismatches = querylength5 + querylength3;
    for (i = 0; i < (*hit5)->end_nambi; i++) {
      genomicend = splicesites[(*hit5)->end_ambi[i]];
      for (j = 0; j < (*hit3)->start_nambi; j++) {
	genomicstart = splicesites[(*hit3)->start_ambi[j]];
	debug9(printf(" %lu,%lu",genomicend,genomicstart));
	if (genomicstart < genomicend) {
	  nbounded++;
	  boundedi5 = i;
	  boundedi3 = j;

#ifdef USE_BINGO
	  insertlength = genomicend - genomicstart + querylength5 + querylength3;
	  debug9(printf(" (%u)",insertlength));
	  if (insertlength < expected_pairlength) {
	    if (expected_pairlength - insertlength <= pairlength_deviation) {
	      nbingo++;
	      bingoi5 = i;
	      bingoi3 = j;
	      debug9(printf("*"));
	    }
	  } else {
	    if (insertlength - expected_pairlength <= pairlength_deviation) {
	      nbingo++;
	      bingoi5 = i;
	      bingoi3 = j;
	      debug9(printf("*"));
	    }
	  }
#endif

	  if ((nmismatches = (*hit5)->end_amb_nmismatches[i] + (*hit3)->start_amb_nmismatches[j]) < best_nmismatches) {
	    best_nmismatches = nmismatches;
	    besti5 = i;
	    besti3 = j;
	    nbest = 1;
	  } else if (nmismatches == best_nmismatches) {
	    nbest++;
	  }
	}
      }
    }

#if 0
    /* No longer holds for GMAP */
    assert((*hit5)->amb_nmatches_end > 0);
    assert((*hit3)->amb_nmatches_start > 0);
#endif

#ifdef USE_BINGO
    if (nbingo == 1) {
      new5p = true; new3p = true;
    } else if (nbounded == 1) {
      new5p = true; new3p = true; bingoi5 = boundedi5; bingoi3 = boundedi3;
    }
#endif

    if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_nmatches = %d...%d",
		    (*hit5)->amb_nmatches_end,(*hit3)->amb_nmatches_start));
      *unresolved_amb_nmatches = (*hit5)->amb_nmatches_end + (*hit3)->amb_nmatches_start;
    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      new5p = true; new3p = true; bingoi5 = besti5; bingoi3 = besti3;
    }
    debug9(printf("\n"));

  } else if ((*hit5)->end_ambiguous_p == true) {
    debug9(printf("Got ambiguous at 5':"));
    nbest = nbounded = nbingo = 0;
    best_nmismatches = querylength5;
    for (i = 0; i < (*hit5)->end_nambi; i++) {
      genomicend = splicesites[(*hit5)->end_ambi[i]];
      debug9(printf(" %lu",genomicend));
      if ((*hit3)->genomicstart < genomicend /*allow overlap*/+ querylength3) {
	nbounded++;
	boundedi5 = i;
	boundedi3 = j;

#ifdef USE_BINGO
	insertlength = genomicend - (*hit3)->genomicstart + querylength5 + querylength3;
	debug9(printf(" (%u)",insertlength));
	if (insertlength < expected_pairlength) {
	  if (expected_pairlength - insertlength <= pairlength_deviation) {
	    nbingo++;
	    bingoi5 = i;
	    debug9(printf("*"));
	  }
	} else {
	  if (insertlength - expected_pairlength <= pairlength_deviation) {
	    nbingo++;
	    bingoi5 = i;
	    debug9(printf("*"));
	  }
	}
#endif

	if ((nmismatches = (*hit5)->end_amb_nmismatches[i]) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	  besti5 = i;
	  nbest = 1;
	} else if (nmismatches == best_nmismatches) {
	  nbest++;
	}
      }
    }

#if 0
    /* No longer holds for GMAP */
    assert((*hit5)->amb_nmatches_end > 0);
    assert((*hit3)->amb_nmatches_start == 0);
#endif

#ifdef USE_BINGO
    if (nbingo == 1) {
      new5p = true;
    } else if (nbounded == 1) {
      new5p = true; bingoi5 = boundedi5;
    }
#endif

    if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_nmatches = %d...%d",
		    (*hit5)->amb_nmatches_end,(*hit3)->amb_nmatches_start));
      *unresolved_amb_nmatches = (*hit5)->amb_nmatches_end;
    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      new5p = true; bingoi5 = besti5;
    }
    debug9(printf("\n"));

  } else if ((*hit3)->start_ambiguous_p == true) {
    debug9(printf("Got ambiguous at 3':"));
    nbest = nbounded = nbingo = 0;
    best_nmismatches = querylength3;
    for (j = 0; j < (*hit3)->start_nambi; j++) {
      genomicstart = splicesites[(*hit3)->start_ambi[j]];
      debug9(printf(" %lu",genomicstart));
      if (genomicstart < (*hit5)->genomicend /*allow overlap*/+ querylength5) {
	nbounded++;
	boundedi3 = j;

#ifdef USE_BINGO
	insertlength = (*hit5)->genomicend - genomicstart + querylength5 + querylength3;
	debug9(printf(" (%u)",insertlength));
	if (insertlength < expected_pairlength) {
	  if (expected_pairlength - insertlength <= pairlength_deviation) {
	    nbingo++;
	    bingoi3 = j;
	    debug9(printf("*"));
	  }
	} else {
	  if (insertlength - expected_pairlength <= pairlength_deviation) {
	    nbingo++;
	    bingoi3 = j;
	    debug9(printf("*"));
	  }
	}
#endif

	if ((nmismatches = (*hit3)->start_amb_nmismatches[j]) < best_nmismatches) {
	  best_nmismatches = nmismatches;
	  besti3 = j;
	  nbest = 1;
	} else if (nmismatches == best_nmismatches) {
	  nbest++;
	}
      }
    }

#if 0
    /* No longer holds for GMAP */
    assert((*hit5)->amb_nmatches_end == 0);
    assert((*hit3)->amb_nmatches_start > 0);
#endif

#ifdef USE_BINGO
    if (nbingo == 1) {
      new3p = true;
    } else if (nbounded == 1) {
      new3p = true; bingoi3 = boundedi3;
    }
#endif

    if (nbest == 0) {
      debug9(printf("\nnbest is zero: amb_nmatches = %d...%d",
		    (*hit5)->amb_nmatches_end,(*hit3)->amb_nmatches_start));
      *unresolved_amb_nmatches = (*hit3)->amb_nmatches_start;
    } else if (nbest == 1) {
      debug9(printf("\nnbest is 1, with nmismatches %d\n",best_nmismatches));
      new3p = true; bingoi3 = besti3;
    }
    debug9(printf("\n"));
  }

  if (new5p == false) {
    /* Skip */

  } else if ((*hit5)->hittype == ONE_THIRD_SHORTEXON || (*hit5)->hittype == TWO_THIRDS_SHORTEXON) {
    if ((*hit5)->sensedir == SENSE_FORWARD) {
      /* End 3 */
      shortexon = (*hit5)->substring1;

      donor_knowni = Substring_splicesites_i_D(shortexon);
      splice_pos = Substring_chimera_pos_D(shortexon);
      acceptor_knowni = (*hit5)->end_ambi[bingoi5];
      nmismatches_shortend = (*hit5)->end_amb_nmismatches[bingoi5];
      segment_left = splicesites[acceptor_knowni] - (querylength5 - splice_pos);

      if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,querylength5 - splice_pos,
					     nmismatches_shortend,
					     /*prob*/2.0,segment_left,query5_compress_rev,
					     querylength5,/*plusp*/false,genestrand,/*sensep*/true,
					     Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					     Substring_chrhigh(shortexon),Substring_chrlength(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 3: Splice from donor #%d to acceptor #%d, with nmismatches %d\n",
		      donor_knowni,acceptor_knowni,nmismatches_shortend));
	old = *hit5;
	*hit5 = Stage3end_new_shortexon(&ignore_found_score,/*donor*/(*hit5)->substringD,acceptor,shortexon,
					/*acceptor_distance*/(*hit5)->acceptor_distance,
					/*donor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
					(*hit5)->amb_nmatches_donor,/*amb_nmatches_acceptor*/0,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength5,/*sensedir*/SENSE_FORWARD,
					/*sarrayp*/false);
	if (*private5p == true) {
	  Stage3end_free(&old);
	}
	*private5p = true;
      }

    } else if ((*hit5)->sensedir == SENSE_ANTI) {
      /* End 8 */
      shortexon = (*hit5)->substring1;

      acceptor_knowni = Substring_splicesites_i_A(shortexon);
      splice_pos = Substring_chimera_pos_A(shortexon);
      donor_knowni = (*hit5)->end_ambi[bingoi5];
      nmismatches_shortend = (*hit5)->end_amb_nmismatches[bingoi5];
      segment_left = splicesites[donor_knowni] - (querylength5 - splice_pos);

      if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,querylength5 - splice_pos,
				       nmismatches_shortend,
				       /*prob*/2.0,segment_left,query5_compress_rev,
				       querylength5,/*plusp*/false,genestrand,/*sensep*/false,
				       Substring_chrnum(shortexon),Substring_chroffset(shortexon),
				       Substring_chrhigh(shortexon),Substring_chrlength(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 8: Splice from antiacceptor #%d to antidonor #%d, with nmismatches_shortend %d\n",
		      acceptor_knowni,donor_knowni,nmismatches_shortend));
	old = *hit5;
	*hit5 = Stage3end_new_shortexon(&ignore_found_score,donor,/*acceptor*/(*hit5)->substringA,shortexon,
					/*acceptor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
					/*donor_distance*/(*hit5)->donor_distance,
					/*amb_nmatches_donor*/0,(*hit5)->amb_nmatches_acceptor,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength5,/*sensedir*/SENSE_ANTI,
					/*sarrayp*/false);
	if (*private5p == true) {
	  Stage3end_free(&old);
	}
	*private5p = true;
      }
	
    } else {
      fprintf(stderr,"Shortexon hit5 has no sensedir\n");
      abort();
    }

  } else if ((*hit5)->hittype == HALFSPLICE_DONOR) {
    /* End 3 */
    assert((*hit5)->sensedir == SENSE_FORWARD);
    donor = (*hit5)->substring_donor;

    donor_knowni = Substring_splicesites_i(donor);
    splice_pos = Substring_chimera_pos(donor);
    acceptor_knowni = (*hit5)->end_ambi[bingoi5];
    nmismatches_shortend = (*hit5)->end_amb_nmismatches[bingoi5];
    segment_left = splicesites[acceptor_knowni] - (querylength5 - splice_pos);

    if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,querylength5 - splice_pos,
					   nmismatches_shortend,
					   /*prob*/2.0,segment_left,query5_compress_rev,
					   querylength5,/*plusp*/false,genestrand,/*sensep*/true,
					   Substring_chrnum(donor),Substring_chroffset(donor),
					   Substring_chrhigh(donor),Substring_chrlength(donor))) != NULL) {
      debug9(printf("Resolved halfsplice, End 3: Splice from donor #%d to acceptor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(donor),Substring_splicesites_i(acceptor),nmismatches_shortend));
      old = *hit5;
      *hit5 = Stage3end_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				   donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength5,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,
				   /*sensedir*/SENSE_FORWARD,/*sarrayp*/false);
      if (*private5p == true) {
	Stage3end_free(&old);
      }
      *private5p = true;
    }

  } else if ((*hit5)->hittype == HALFSPLICE_ACCEPTOR) {
    /* End 8 */
    assert((*hit5)->sensedir == SENSE_ANTI);
    acceptor = (*hit5)->substring_acceptor;

    acceptor_knowni = Substring_splicesites_i(acceptor);
    splice_pos = Substring_chimera_pos(acceptor);
    donor_knowni = (*hit5)->end_ambi[bingoi5];
    nmismatches_shortend = (*hit5)->end_amb_nmismatches[bingoi5];
    segment_left = splicesites[donor_knowni] - (querylength5 - splice_pos);

    if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,querylength5 - splice_pos,
				     nmismatches_shortend,
				     /*prob*/2.0,segment_left,query5_compress_rev,
				     querylength5,/*plusp*/false,genestrand,/*sensep*/false,
				     Substring_chrnum(acceptor),Substring_chroffset(acceptor),
				     Substring_chrhigh(acceptor),Substring_chrlength(acceptor))) != NULL) {
      debug9(printf("Resolved halfsplice acceptor, End 8: Splice from antiacceptor #%d to antidonor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(acceptor),Substring_splicesites_i(donor),nmismatches_shortend));
      old = *hit5;
      *hit5 = Stage3end_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				   donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength5,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,
				   /*sensedir*/SENSE_ANTI,/*sarrayp*/false);
      if (*private5p == true) {
	Stage3end_free(&old);
      }
      *private5p = true;
    }

  } else {
    fprintf(stderr,"Unexpected hittype %d for ambiguous end\n",(*hit5)->hittype);
    abort();
  }


  if (new3p == false) {
    /* Skip */

  } else if ((*hit3)->hittype == ONE_THIRD_SHORTEXON || (*hit3)->hittype == TWO_THIRDS_SHORTEXON) {
    if ((*hit3)->sensedir == SENSE_ANTI) {
      /* End 7 */
      shortexon = (*hit3)->substring1;

      donor_knowni = Substring_splicesites_i_D(shortexon);
      splice_pos = Substring_chimera_pos_D(shortexon);
      acceptor_knowni = (*hit3)->start_ambi[bingoi3];
      nmismatches_shortend = (*hit3)->start_amb_nmismatches[bingoi3];
      segment_left = splicesites[acceptor_knowni] - (querylength3 - splice_pos);

      if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,querylength3 - splice_pos,
					     nmismatches_shortend,
					     /*prob*/2.0,segment_left,query3_compress_rev,
					     querylength3,/*plusp*/false,genestrand,/*sensep*/false,
					     Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					     Substring_chrhigh(shortexon),Substring_chrlength(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 7: Splice from antidonor #%d to antiacceptor #%d, with nmismatches %d\n",
		      donor_knowni,acceptor_knowni,nmismatches_shortend));
	old = *hit3;
	*hit3 = Stage3end_new_shortexon(&ignore_found_score,/*donor*/(*hit3)->substringD,acceptor,shortexon,
					/*acceptor_distance*/(*hit3)->acceptor_distance,
					/*donor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
					(*hit3)->amb_nmatches_donor,/*amb_nmatches_acceptor*/0,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength3,/*sensedir*/SENSE_ANTI,
					/*sarrayp*/false);
	if (*private3p == true) {
	  Stage3end_free(&old);
	}
	*private3p = true;
      }

    } else if ((*hit3)->sensedir == SENSE_FORWARD) {
      /* End 4 */
      shortexon = (*hit3)->substring1;

      acceptor_knowni = Substring_splicesites_i_A(shortexon);
      splice_pos = Substring_chimera_pos_A(shortexon);
      donor_knowni = (*hit3)->start_ambi[bingoi3];
      nmismatches_shortend = (*hit3)->start_amb_nmismatches[bingoi3];
      segment_left = splicesites[donor_knowni] - (querylength3 - splice_pos);

      if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,querylength3 - splice_pos,
				       nmismatches_shortend,
				       /*prob*/2.0,segment_left,query3_compress_rev,
				       querylength3,/*plusp*/false,genestrand,/*sensep*/true,
				       Substring_chrnum(shortexon),Substring_chroffset(shortexon),
				       Substring_chrhigh(shortexon),Substring_chrlength(shortexon))) != NULL) {
	debug9(printf("Resolved halfsplice_acceptor, End 4: Splice from acceptor #%d to #%d, with nmismatches %d\n",
		      acceptor_knowni,donor_knowni,nmismatches_shortend));
	old = *hit3;
	*hit3 = Stage3end_new_shortexon(&ignore_found_score,donor,/*acceptor*/(*hit3)->substringA,shortexon,
					/*acceptor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
					/*donor_distance*/(*hit3)->donor_distance,
					/*amb_nmatches_donor*/0,(*hit3)->amb_nmatches_acceptor,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength3,/*sensedir*/SENSE_FORWARD,
					/*sarrayp*/false);
	if (*private3p == true) {
	  Stage3end_free(&old);
	}
	*private3p = true;
      }

    } else {
      fprintf(stderr,"Shortexon hit3 has no sensedir\n");
      abort();
    }

  } else if ((*hit3)->hittype == HALFSPLICE_DONOR) {
    /* End 7 */
    assert((*hit3)->sensedir == SENSE_ANTI);
    donor = (*hit3)->substring_donor;

    donor_knowni = Substring_splicesites_i(donor);
    splice_pos = Substring_chimera_pos(donor);
    acceptor_knowni = (*hit3)->start_ambi[bingoi3];
    nmismatches_shortend = (*hit3)->start_amb_nmismatches[bingoi3];
    segment_left = splicesites[acceptor_knowni] - (querylength3 - splice_pos);

    if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,querylength3 - splice_pos,
					   nmismatches_shortend,
					   /*prob*/2.0,segment_left,query3_compress_rev,
					   querylength3,/*plusp*/false,genestrand,/*sensep*/false,
					   Substring_chrnum(donor),Substring_chroffset(donor),
					   Substring_chrhigh(donor),Substring_chrlength(donor))) != NULL) {
      debug9(printf("Resolved halfsplice_donor, End 7: Splice from antidonor #%d to antiacceptor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(donor),Substring_splicesites_i(acceptor),nmismatches_shortend));
      old = *hit3;
      *hit3 = Stage3end_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				   donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength3,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				   /*sensedir*/SENSE_ANTI,/*sarrayp*/false);
      if (*private3p == true) {
	Stage3end_free(&old);
      }
      *private3p = true;
    }

  } else if ((*hit3)->hittype == HALFSPLICE_ACCEPTOR) {
    /* End 4 */
    assert((*hit3)->sensedir == SENSE_FORWARD);
    acceptor = (*hit3)->substring_acceptor;

    acceptor_knowni = Substring_splicesites_i(acceptor);
    splice_pos = Substring_chimera_pos(acceptor);
    donor_knowni = (*hit3)->start_ambi[bingoi3];
    nmismatches_shortend = (*hit3)->start_amb_nmismatches[bingoi3];
    segment_left = splicesites[donor_knowni] - (querylength3 - splice_pos);

    if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,querylength3 - splice_pos,
				     nmismatches_shortend,
				     /*prob*/2.0,segment_left,query3_compress_rev,
				     querylength3,/*plusp*/false,genestrand,/*sensep*/true,
				     Substring_chrnum(acceptor),Substring_chroffset(acceptor),
				     Substring_chrhigh(acceptor),Substring_chrlength(acceptor))) != NULL) {
      debug9(printf("Resolved halfsplice_acceptor, End 4: Splice from acceptor #%d to #%d, with nmismatches %d\n",
		    Substring_splicesites_i(acceptor),Substring_splicesites_i(donor),nmismatches_shortend));
      old = *hit3;
      *hit3 = Stage3end_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				   donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength3,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				   /*sensedir*/SENSE_FORWARD,/*sarrayp*/false);
      if (*private3p == true) {
	Stage3end_free(&old);
      }
      *private3p = true;
    }

  } else {
    fprintf(stderr,"Unexpected hittype %d for ambiguous end\n",(*hit3)->hittype);
    abort();
  }

  return;
}



static void
alias_circular (T hit) {
  Chrpos_T chrlength = hit->chrlength;

  assert(hit->alias == -1);
  if (hit->hittype == GMAP) {
    Pair_alias_circular(hit->pairarray,hit->npairs,chrlength);

  } else {
    hit->genomicstart += chrlength;
    hit->genomicend += chrlength;
    hit->low += chrlength;
    hit->high += chrlength;
    Substring_alias_circular(hit->substring0);
    Substring_alias_circular(hit->substring1);
    Substring_alias_circular(hit->substring2);
  }
  hit->alias = +1;

  return;
}



static void
unalias_circular (T hit) {
  Chrpos_T chrlength = hit->chrlength;

  assert(hit->alias == +1);
  if (hit->hittype == GMAP) {
    Pair_unalias_circular(hit->pairarray,hit->npairs,chrlength);

  } else {
    hit->genomicstart -= chrlength;
    hit->genomicend -= chrlength;
    hit->low -= chrlength;
    hit->high -= chrlength;
    Substring_unalias_circular(hit->substring0);
    Substring_unalias_circular(hit->substring1);
    Substring_unalias_circular(hit->substring2);
  }
  hit->alias = -1;

  return;
}


List_T
Stage3end_unalias_circular (List_T hitlist) {
  List_T p;
  T hit;

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->alias == +1) {
      unalias_circular(hit);
    }
  }

  return hitlist;
}



Stage3pair_T
Stage3pair_new (T hit5, T hit3,	Univcoord_T *splicesites,
		Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		int genestrand,	Pairtype_T pairtype, int localsplicing_penalty,
		bool private5p, bool private3p, bool expect_concordant_p) {
  Stage3pair_T new;
  Stage3end_T copy;
  Chrpos_T chrstart, chrend, chrpos;
  int querypos;
  int unresolved_amb_nmatches = 0;
  int found_score = 0;
  bool overreach5p, overreach3p;

  int querylength5 = hit5->querylength_adj;
  int querylength3 = hit3->querylength_adj;

  debug10(printf("\nStage3pair_new called with pairtype %s and chrnum %d, %d (effective %d, %d)\n",
		 Pairtype_string(pairtype),hit5->chrnum,hit3->chrnum,hit5->effective_chrnum,hit3->effective_chrnum));

  if (hit5->hittype == TERMINAL && Stage3end_terminal_length(hit5) < terminal_output_minlength) {
    debug10(printf("5' terminal length %d is not long enough (< %d)\n",
		   Stage3end_terminal_length(hit5),terminal_output_minlength));
    if (private5p == true) {
      Stage3end_free(&hit5);
    }
    if (private3p == true) {
      Stage3end_free(&hit3);
    }
    return (Stage3pair_T) NULL;
  } else if (hit3->hittype == TERMINAL && Stage3end_terminal_length(hit3) < terminal_output_minlength) {
    debug10(printf("3' terminal length %d is not long enough (< %d)\n",
		   Stage3end_terminal_length(hit3),terminal_output_minlength));
    if (private5p == true) {
      Stage3end_free(&hit5);
    }
    if (private3p == true) {
      Stage3end_free(&hit3);
    }
    return (Stage3pair_T) NULL;
  } else {
    new = (Stage3pair_T) MALLOC_OUT(sizeof(*new));
  }

  if (pairtype == PAIRED_UNSPECIFIED || pairtype == UNSPECIFIED) {
    /* Can get here from running GMAP improvement on a paired result */
    pairtype = Stage3_determine_pairtype(hit5,hit3);
    debug10(printf("  Changing pairtype to %s\n",Pairtype_string(pairtype)));
    if (pairtype == CONCORDANT) {
      expect_concordant_p = true;
    }
  }
  new->pairtype = pairtype;
  new->genestrand = genestrand;

#if 0
  new->mapq_loglik = hit5->mapq_loglik + hit3->mapq_loglik;
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  if (hit5->hittype == GMAP && hit3->hittype == GMAP) {
    debug10(printf("Got hit5 and hit3 both of type GMAP\n"));

    /* Do not try to resolve ambiguity on inside of concordant ends */
    if (hit5->plusp == true && hit3->plusp == true) {
      new->dir = +1;
      new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      debug10(printf("plus, no overlap: insert length %d = start3 %lu - end5 %lu + %d + %d\n",
		     new->insertlength,hit3->genomicstart,hit5->genomicend,querylength5,querylength3));
    } else if (hit5->plusp == false && hit3->plusp == false) {
      new->dir = -1;
      new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      debug10(printf("minus, no overlap: insert length %d = end5 %lu - start3 %lu + %d + %d\n",
		     new->insertlength,hit5->genomicend,hit3->genomicstart,querylength5,querylength3));
    } else {
      new->dir = 0;
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit5->hittype == GMAP) {
    debug10(printf("Got hit5 of type GMAP\n"));
    if (hit5->plusp == true && hit3->plusp == true) {
      new->dir = +1;

      if (expect_concordant_p == true) {
	/* Try to resolve ambiguity on inside of concordant ends */
	resolve_inside_ambiguous_splice_plus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					     splicesites,query5_compress_fwd,query3_compress_fwd,
					     localsplicing_penalty,querylength5,querylength3,genestrand);
      }

      /* Have 5-start..end and 3-start..end */
      debug10(printf("plus: comparing hit5->genomicend %lu <= hit3->genomicstart %lu\n",
		     hit5->genomicend,hit3->genomicstart));

      if (hit5->genomicend <= hit3->genomicstart) {
	/* No overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("plus, no overlap: insert length %d = start3 %lu - end5 %lu + %d + %d\n",
		       new->insertlength,hit3->genomicstart,hit5->genomicend,querylength5,querylength3));
      } else if ((chrpos = overlap3_gmap_plus(&querypos,&chrstart,&chrend,/*hit*/hit3,/*gmap*/hit5)) > 0U) {
	new->insertlength = /* end3 */ chrend - /* start5 */ (chrpos - querypos);
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("plus, overlap: insert length %d = end3 %lu - start5 (%lu - %d)\n",
		       new->insertlength,chrend,chrpos,querypos));
      } else {
	/* Still no overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      }

    } else if (hit5->plusp == false && hit3->plusp == false) {
      new->dir = -1;

      if (expect_concordant_p == true) {
	/* Try to resolve ambiguity on inside of concordant ends */
	resolve_inside_ambiguous_splice_minus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					      splicesites,query5_compress_rev,query3_compress_rev,
					      localsplicing_penalty,querylength5,querylength3,genestrand);
      }

      /* Have 3-end..start and 5-end..start */
      debug10(printf("minus: comparing hit3->genomicstart %lu <= hit5->genomicend %lu\n",
		     hit3->genomicstart,hit5->genomicend));

      if (hit3->genomicstart <= hit5->genomicend) {
	/* No overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("minus, no overlap: insert length %d = end5 %lu - start3 %lu + %d + %d\n",
		       new->insertlength,hit5->genomicend,hit3->genomicstart,querylength5,querylength3));
      } else if ((chrpos = overlap3_gmap_minus(&querypos,&chrstart,&chrend,/*hit*/hit3,/*gmap*/hit5)) > 0U) {
	new->insertlength = /* start5 */ (chrpos + querypos) - /* end3 */ chrend + 1;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("minus, overlap: insert length %d = start5 (%lu + %d) - end3 %lu + 1\n",
		       new->insertlength,chrpos,querypos,chrend));
      } else {
	/* Still no overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      }
    } else {
      new->dir = 0;
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit3->hittype == GMAP) {
    debug10(printf("Got hit3 of type GMAP\n"));
    if (hit5->plusp == true && hit3->plusp == true) {
      new->dir = +1;

      if (expect_concordant_p == true) {
	/* Try to resolve ambiguity on inside of concordant ends */
	resolve_inside_ambiguous_splice_plus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					     splicesites,query5_compress_fwd,query3_compress_fwd,
					     localsplicing_penalty,querylength5,querylength3,genestrand);
      }

      /* Have 5-start..end and 3-start..end */
      debug10(printf("plus: comparing hit5->genomicend %lu <= hit3->genomicstart %lu\n",
		     hit5->genomicend,hit3->genomicstart));

      if (hit5->genomicend <= hit3->genomicstart) {
	/* No overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("plus, no overlap: insert length %d = start3 %lu - end5 %lu + %d + %d\n",
		       new->insertlength,hit3->genomicstart,hit5->genomicend,querylength5,querylength3));
      } else if ((chrpos = overlap5_gmap_plus(&querypos,&chrstart,&chrend,/*hit*/hit5,/*gmap*/hit3)) > 0U) {
	new->insertlength = /* end3 */ (chrpos - querypos + querylength3) - /* start5 */ chrstart;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("plus, overlap: insert length %d = end3 (%lu - %d + %d) - start5 %lu\n",
		       new->insertlength,chrpos,querypos,querylength3,chrstart));
      } else {
	/* Still no overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      }

    } else if (hit5->plusp == false && hit3->plusp == false) {
      new->dir = -1;

      if (expect_concordant_p == true) {
	/* Try to resolve ambiguity on inside of concordant ends */
	resolve_inside_ambiguous_splice_minus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					      splicesites,query5_compress_rev,query3_compress_rev,
					      localsplicing_penalty,querylength5,querylength3,genestrand);
      }

      /* Have 3-end..start and 5-end..start */
      debug10(printf("minus: comparing hit3->genomicstart %lu <= hit5->genomicend %lu\n",
		     hit3->genomicstart,hit5->genomicend));
      if (hit3->genomicstart <= hit5->genomicend) {
	/* No overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("minus, no overlap: insert length %d = end5 %lu - start3 %lu + %d + %d\n",
		       new->insertlength,hit5->genomicend,hit3->genomicstart,querylength5,querylength3));
      } else if ((chrpos = overlap5_gmap_minus(&querypos,&chrstart,&chrend,/*hit*/hit5,/*gmap*/hit3)) > 0U) {
	new->insertlength = /* start5 */ chrstart - /* end3 */ (chrpos + querypos - querylength3) - 1;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
	debug10(printf("minus, overlap: insert length %d = start5 %lu - end3 (%lu + %d - %d) - 1\n",
		       new->insertlength,chrstart,chrpos,querypos,querylength3));
      } else {
	/* Still no overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
	new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      }
    } else {
      new->dir = 0;
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit5->plusp == true && hit3->plusp == false) {
    new->dir = 0;
    
    /* Have 5-start..end and 3-end..start */
    /*   or 3-end..start and 5-start..end */

    if (hit5->genomicend < hit3->genomicend) {
      new->insertlength = (hit3->genomicend - hit5->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else if (hit3->genomicstart < hit5->genomicstart) {
      new->insertlength = (hit5->genomicstart - hit3->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else {
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit5->plusp == false && hit3->plusp == true) {
    new->dir = 0;
    
    /* Have 5-end..start and 3-start..end */
    /*   or 3-start..end and 5-end..start */

    if (hit5->genomicstart < hit3->genomicstart) {
      new->insertlength = (hit3->genomicstart - hit5->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else if (hit3->genomicend < hit5->genomicend) {
      new->insertlength = (hit5->genomicend - hit3->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    } else {
      new->insertlength = pair_insert_length_unpaired(hit5,hit3); /* was 0 */
      new->insertlength_expected_sign = false;
    }

  } else if (hit5->plusp == true) {
    /* Concordant directions on same chromosome (plus) */
    debug10(printf("Concordant on plus strand\n"));
    new->dir = +1;

    if (expect_concordant_p == true) {
      overreach5p = overreach3p = false;
      if (hit5->hittype == SPLICE) {
	if (Substring_alignstart(hit5->substring2) > hit3->genomicend) {
	  if (Substring_alignend(hit5->substring1) < hit3->genomicstart) {
	    overreach5p = true;
	  }
	}
      }
      if (hit3->hittype == SPLICE) {
	if (Substring_alignend(hit3->substring1) < hit5->genomicstart) {
	  if (Substring_alignstart(hit3->substring2) > hit5->genomicend) {
	    overreach3p = true;
	  }
	}
      }

      if (overreach5p == true || overreach3p == true) {
	/* Either overreach */
	debug5(printf("  Returning NULL because of dual overreach\n"));
	if (private5p == true) {
	  Stage3end_free(&hit5);
	}
	if (private3p == true) {
	  Stage3end_free(&hit3);
	}
	FREE_OUT(new);
	return (Stage3pair_T) NULL;

#if 0
      } else if (overreach5p == true) {
	/* Overreach of hit5 */
	debug9(printf("Overreach of hit5 of type SPLICE.  Removing substring2\n"));
	if (hit5->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_whole(hit5->substring1),
				      /*nmismatches_acceptor*/0,/*donor*/hit5->substring1,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,
				      /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/0,/*amb_nmismatches_right*/0,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,
				      /*sensedir*/hit5->sensedir,hit5->sarrayp);
	} else if (hit5->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_whole(hit5->substring1),/*donor*/NULL,
				      /*acceptor*/hit5->substring1,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,
				      /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/0,/*amb_nmismatches_right*/0,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,
				      /*sensedir*/hit5->sensedir,hit5->sarrayp);
	} else {
	  abort();
	}
	if (private5p == true) {
	  Stage3end_free(&hit5);
	} else {
	  private5p = true;
	}
	hit5 = copy;

      } else if (overreach3p == true) {
	/* Overreach of hit3 */
	debug9(printf("Overreach of hit3 of type SPLICE.  Removing substring1\n"));
	if (hit3->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_whole(hit3->substring2),/*donor*/NULL,
				      /*acceptor*/hit3->substring2,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,
				      /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/0,/*amb_nmismatches_right*/0,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,
				      /*sensedir*/hit3->sensedir,hit3->sarrayp);
	} else if (hit3->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_whole(hit3->substring2),
				      /*nmismatches_acceptor*/0,/*donor*/hit3->substring2,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,
				      /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/0,/*amb_nmismatches_right*/0,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,
				      /*sensedir*/hit3->sensedir,hit3->sarrayp);
	} else {
	  abort();
	}
	if (private3p == true) {
	  Stage3end_free(&hit3);
	} else {
	  private3p = true;
	}
	hit3 = copy;
#endif
      }

      /* Try to resolve ambiguity on inside of concordant ends */
      resolve_inside_ambiguous_splice_plus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					   splicesites,query5_compress_fwd,query3_compress_fwd,
					   localsplicing_penalty,querylength5,querylength3,genestrand);
    }

    /* Have 5-start..end and 3-start..end */
    if (hit5->genomicend < hit3->genomicstart) {
      /* No overlap */
      new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      debug10(printf("plus, no overlap: insert length %d = start3 %lu - end5 %lu + %d + %d\n",
		     new->insertlength,hit3->genomicstart,hit5->genomicend,querylength5,querylength3));
#if 0
    } else if (hit5->genomicend > hit3->genomicend + SUBSUMPTION_SLOP) {
      /* hit5 subsumes hit3 */
      debug10(printf("plus, subsumption %lu > %lu\n",hit5->genomicend,hit3->genomicend));
      new->insertlength = 0;
      new->insertlength_expected_sign = false;
#endif
    } else {
      new->insertlength = pair_insert_length(hit5,hit3);
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    }


  } else {
    /* Concordant directions on same chromosome (minus) */
    debug10(printf("Concordant on minus strand\n"));
    new->dir = -1;

    if (expect_concordant_p == true) {
      overreach5p = overreach3p = false;
      if (hit5->hittype == SPLICE) {
	if (Substring_alignstart(hit5->substring2) < hit3->genomicend) {
	  if (Substring_alignend(hit5->substring1) > hit3->genomicstart) {
	    overreach5p = true;
	  }
	}
      }
      if (hit3->hittype == SPLICE) {
	if (Substring_alignend(hit3->substring1) > hit5->genomicstart) {
	  if (Substring_alignstart(hit3->substring2) < hit5->genomicend) {
	    overreach3p = true;
	  }
	}
      }

      if (overreach5p == true || overreach3p == true) {
	/* Either overreach */
	debug5(printf("  Returning NULL because of dual overreach\n"));
	if (private5p == true) {
	  Stage3end_free(&hit5);
	}
	if (private3p == true) {
	  Stage3end_free(&hit3);
	}
	FREE_OUT(new);
	return (Stage3pair_T) NULL;

#if 0
      } else if (overreach5p == true) {
	/* Overreach of hit5 */
	debug9(printf("Overreach of hit5 of type SPLICE.  Removing substring2\n"));
	if (hit5->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_whole(hit5->substring1),
				      /*nmismatches_acceptor*/0,/*donor*/hit5->substring1,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,
				      /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/0,/*amb_nmismatches_right*/0,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				      /*sensedir*/hit5->sensedir,hit5->sarrayp);
	} else if (hit5->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_whole(hit5->substring1),/*donor*/NULL,
				      /*acceptor*/hit5->substring1,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit5->querylength,
				      /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/0,/*amb_nmismatches_right*/0,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				      /*sensedir*/hit5->sensedir,hit5->sarrayp);
	} else {
	  abort();
	}
	if (private5p == true) {
	  Stage3end_free(&hit5);
	} else {
	  private5p = true;
	}
	hit5 = copy;

      } else if (overreach3p == true) {
	/* Overreach of hit3 */
	debug9(printf("Overreach of hit3 of type SPLICE.  Removing substring1\n"));
	if (hit3->sensedir == SENSE_FORWARD) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/0,
				      /*nmismatches_acceptor*/Substring_nmismatches_whole(hit3->substring2),/*donor*/NULL,
				      /*acceptor*/hit3->substring2,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,
				      /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/0,/*amb_nmismatches_right*/0,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				      /*sensedir*/hit3->sensedir,hit3->sarrayp);
	} else if (hit3->sensedir == SENSE_ANTI) {
	  copy = Stage3end_new_splice(&found_score,/*nmismatches_donor*/Substring_nmismatches_whole(hit3->substring2),
				      /*nmismatches_acceptor*/0,/*donor*/hit3->substring2,/*acceptor*/NULL,/*distance*/0U,
				      /*shortdistancep*/true,localsplicing_penalty,hit3->querylength,
				      /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/0,/*amb_nmismatches_right*/0,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				      /*sensedir*/hit3->sensedir,hit3->sarrayp);
	} else {
	  abort();
	}
	if (private3p == true) {
	  Stage3end_free(&hit3);
	} else {
	  private3p = true;
	}
	hit3 = copy;
#endif
      }

      /* Try to resolve ambiguity on inside of concordant ends */
      resolve_inside_ambiguous_splice_minus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					    splicesites,query5_compress_rev,query3_compress_rev,
					    localsplicing_penalty,querylength5,querylength3,genestrand);
    }

    /* Have 3-end..start and 5-end..start */
    if (hit3->genomicstart < hit5->genomicend) {
      /* No overlap */
      new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
      debug10(printf("minus, no overlap: insert length %d = end5 %lu - start3 %lu + %d + %d\n",
		     new->insertlength,hit5->genomicend,hit3->genomicstart,querylength5,querylength3));
#if 0
    } else if (hit3->genomicstart > hit5->genomicstart + SUBSUMPTION_SLOP) {
      /* hit3 subsumes hit5 */
      debug10(printf("minus, subsumption %lu > %lu\n",hit3->genomicstart,hit5->genomicstart));
      new->insertlength = 0;
      new->insertlength_expected_sign = false;
#endif
    } else {
      new->insertlength = pair_insert_length(hit5,hit3);
      new->insertlength_expected_sign = insertlength_expected(new->insertlength);
    }

  }


  debug5(printf("\nGot insertlength of %d\n",new->insertlength));
  if (new->insertlength <= 0) {
    /* Not concordant */
#ifdef USE_BINGO
    new->absdifflength_bingo_p = false;
#endif
#ifdef USE_ABSDIFFLENGTH
    new->absdifflength = -1U;
#endif

    if (expect_concordant_p == true) {
      debug5(printf("  Returning NULL\n"));
      if (private5p == true) {
	Stage3end_free(&hit5);
      }
      if (private3p == true) {
	Stage3end_free(&hit3);
      }
      FREE_OUT(new);
      return (Stage3pair_T) NULL;
    }

  } else if (new->insertlength > pairmax && expect_concordant_p == true) {
    debug5(printf("  Returning NULL\n"));
    if (private5p == true) {
      Stage3end_free(&hit5);
    }
    if (private3p == true) {
      Stage3end_free(&hit3);
    }
    FREE_OUT(new);
    return (Stage3pair_T) NULL;

  } else {
#ifdef USE_ABSDIFFLENGTH
    if (new->insertlength < expected_pairlength) {
      new->absdifflength = expected_pairlength - new->insertlength;
    } else {
      new->absdifflength = new->insertlength - expected_pairlength;
    }
#endif
#ifdef USE_BINGO
    if (new->absdifflength <= pairlength_deviation) {
      new->absdifflength_bingo_p = true;
    } else {
      new->absdifflength_bingo_p = false;
    }
#endif
  }

  if (SENSE_CONSISTENT_P(hit5->sensedir,hit3->sensedir)) {
    debug0(printf("senses are consistent\n"));
    new->sense_consistent_p = true;
  } else {
    debug0(printf("senses are inconsistent\n"));
    new->sense_consistent_p = false;
  }

  /* Do not alter score, so the alignmnent terminates at the known splice site  */
  new->score = hit5->score + hit3->score /* + unresolved_amb_nmatches */;
  new->nmatches = hit5->nmatches + hit3->nmatches - unresolved_amb_nmatches;
  new->nmatches_posttrim = hit5->nmatches_posttrim + hit3->nmatches_posttrim;
  new->indel_low = hit5->indel_low + hit3->indel_low;
  /* new->overlap_known_gene_p = false; -- initialized later when resolving multimappers */
  new->tally = -1L;

  new->low = (hit5->low < hit3->low) ? hit5->low : hit3->low;
  new->high = (hit5->high > hit3->high) ? hit5->high : hit3->high;

#if 0
  if (new->low > new->high) {
    fprintf(stderr,"new->low %lu > new->high %lu, hit5->chrnum %d\n",new->low,new->high,hit5->chrnum);
    abort();
  }
#endif

  if (hit5->chrnum == 0 || hit3->chrnum == 0) {
    new->outerlength = querylength5 + querylength3;
  } else {
    new->outerlength = new->high - new->low;
  }

  new->hit5 = hit5;
  new->private5p = private5p;

  new->hit3 = hit3;
  new->private3p = private3p;

  if (expect_concordant_p == true) {
    hit5->paired_usedp = true;
    hit3->paired_usedp = true;
  }

  new->nchimera_known = hit5->nchimera_known + hit3->nchimera_known;
  new->nchimera_novel = hit5->nchimera_novel + hit3->nchimera_novel;

  debug0(printf("Created new pair %p from %p and %p with private %d, %d\n",new,hit5,hit3,private5p,private3p));
  debug0(printf("  hittypes %s and %s\n",hittype_string(hit5->hittype),hittype_string(hit3->hittype)));
  debug0(printf("  sensedirs %d and %d\n",hit5->sensedir,hit3->sensedir));
  debug0(printf("  genomicpos %lu..%lu and %lu..%lu\n",hit5->genomicstart,hit5->genomicend,hit3->genomicstart,hit3->genomicend));

  if (hit5->circularpos < 0 && hit3->circularpos < 0) {
    new->circularp = false;
  } else {
    new->circularp = true;
  }

  if (hit5->alias > 0) {
    debug0(printf("Unaliasing 5' end\n"));
    if (private5p == false) {
      new->hit5 = Stage3end_copy(hit5);
      new->private5p = true;
    }
    unalias_circular(new->hit5);
  }

  if (hit3->alias > 0) {
    debug0(printf("Unaliasing 3' end\n"));
    if (private3p == false) {
      new->hit3 = Stage3end_copy(hit3);
      new->private3p = true;
    }
    unalias_circular(new->hit3);
  }

  return new;
}


void
Stage3pair_privatize (Stage3pair_T *array, int npairs) {
  Stage3pair_T hitpair;
  int i;

  debug0(printf("Call to Stage3pair_privatize on %d hitpairs\n",npairs));

  for (i = 0; i < npairs; i++) {
    hitpair = array[i];
    debug0(printf("  Pair with hitpairs %p (private %d), %p (private %d)\n",
		  hitpair->hit5,hitpair->private5p,hitpair->hit3,hitpair->private3p));
  
    if (hitpair->private5p == false) {
      hitpair->hit5 = Stage3end_copy(hitpair->hit5);
      hitpair->private5p = true;
    }
    if (hitpair->private3p == false) {
      hitpair->hit3 = Stage3end_copy(hitpair->hit3);
      hitpair->private3p = true;
    }
  }

  return;
}



#if 0
static int
chimera_match_distance_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->nmismatches_whole < y->nmismatches_whole) {
    return -1;
  } else if (x->nmismatches_whole > y->nmismatches_whole) {
    return +1;
  } else if (x->distance < y->distance) {
    return -1;
  } else if (x->distance > y->distance) {
    return +1;
  } else {
    return 0;
  }
}
#endif

#if 0
List_T
Stage3end_sort_bymatchdist (List_T hitlist, int maxchimerapaths) {
#ifdef DEBUG
  T hit;
#endif
  List_T sorted = NULL, p;
  T *hits;
  int npaths, n, i;

  if ((n = List_length(hitlist)) == 0) {
    return NULL;
  }

  hits = (T *) CALLOC(n,sizeof(T));
  for (p = hitlist, i = 0; p != NULL; p = p->rest) {
    hits[i++] = (T) p->first;
  }
  List_free(&hitlist);

  qsort(hits,n,sizeof(T),chimera_match_distance_cmp);

  if (n < maxchimerapaths) {
    npaths = n;
  } else {
    npaths = maxchimerapaths;
  }
  for (i = n-1; i >= npaths; i--) {
    Stage3end_free(&(hits[i]));
  }
  for (i = npaths-1; i >= 0; i--) {
    sorted = List_push(sorted,hits[i]);
  }
  FREE(hits);

  debug(
	for (p = sorted, i = 0; p != NULL; p = p->rest, i++) {
	  hit = (T) p->first;
	  printf("  Final %d: chr %d -- %d\n",
		 i,hit->substring1->chrnum,hit->substring2->chrnum);
	}
	);

  return sorted;
}
#endif



/* Used for eliminating exact duplicates.  Also sorts secondarily by hittype. */
static int
hitpair_sort_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;
  
  if (x->dir != 0 && y->dir == 0) {
    return -1;
  } else if (x->dir == 0 && y->dir != 0) {
    return +1;
  } else if (x->dir > 0 && y->dir < 0) {
    return -1;
  } else if (x->dir < 0 && y->dir > 0) {
    return +1;

#if 0
  } else if (x->high < y->low) {
    return -1;
  } else if (y->high < x->low) {
    return +1;

  } else if (x->hit5->high < y->hit5->low) {
    return -1;
  } else if (y->hit5->high < x->hit5->low) {
    return +1;

  } else if (x->hit3->high < y->hit3->low) {
    return -1;
  } else if (y->hit3->high < x->hit3->low) {
    return +1;
#else
  } else if (x->hit5->low < y->hit5->low) {
    return -1;
  } else if (y->hit5->low < x->hit5->low) {
    return +1;

  } else if (x->hit5->high < y->hit5->high) {
    return -1;
  } else if (y->hit5->high < x->hit5->high) {
    return +1;

  } else if (x->hit3->low < y->hit3->low) {
    return -1;
  } else if (y->hit3->low < x->hit3->low) {
    return +1;

  } else if (x->hit3->high < y->hit3->high) {
    return -1;
  } else if (y->hit3->high < x->hit3->high) {
    return +1;
#endif


#if 1
    /* Rank terminals last, so terminals cannot win */
  } else if (x->hit5->hittype != TERMINAL && x->hit3->hittype != TERMINAL && 
	     (y->hit5->hittype == TERMINAL || y->hit3->hittype == TERMINAL)) {
    return -1;
  } else if ((x->hit5->hittype == TERMINAL || x->hit3->hittype == TERMINAL) &&
	     y->hit5->hittype != TERMINAL && y->hit3->hittype != TERMINAL) {
    return +1;
#endif

  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
  } else if (x->nchimera_novel < y->nchimera_novel) {
    return -1;
  } else if (y->nchimera_novel < x->nchimera_novel) {
    return +1;
  } else if (x->nchimera_known > y->nchimera_known) {
    return -1;
  } else if (y->nchimera_known > x->nchimera_known) {
    return +1;
  } else if (x->hit5->hittype < y->hit5->hittype) {
    return -1;
  } else if (y->hit5->hittype < x->hit5->hittype) {
    return +1;
  } else if (x->hit3->hittype < y->hit3->hittype) {
    return -1;
  } else if (y->hit3->hittype < x->hit3->hittype) {
    return +1;
  } else if (x->sense_consistent_p == true && y->sense_consistent_p == false) {
    return -1;
  } else if (x->sense_consistent_p == false && y->sense_consistent_p == true) {
    return +1;
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
  } else {
    return 0;
  }
}


/* Same as hitpair_sort_cmp, except for hittype, nmatches_posttrim, and indel_low */
static int
hitpair_equiv_cmp (Stage3pair_T x, Stage3pair_T y) {

  if (x->dir != 0 && y->dir == 0) {
    return -1;
  } else if (x->dir == 0 && y->dir != 0) {
    return +1;
  } else if (x->dir > 0 && y->dir < 0) {
    return -1;
  } else if (x->dir < 0 && y->dir > 0) {
    return +1;
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high < y->high) {
    return -1;
  } else if (y->high < x->high) {
    return +1;

  } else if (x->hit5->low < y->hit5->low) {
    return -1;
  } else if (y->hit5->low < x->hit5->low) {
    return +1;
  } else if (x->hit5->high < y->hit5->high) {
    return -1;
  } else if (y->hit5->high < x->hit5->high) {
    return +1;

  } else if (x->hit3->low < y->hit3->low) {
    return -1;
  } else if (y->hit3->low < x->hit3->low) {
    return +1;
  } else if (x->hit3->high < y->hit3->high) {
    return -1;
  } else if (y->hit3->high < x->hit3->high) {
    return +1;

#if 0
  } else if (x->hit5->hittype != TERMINAL && x->hit3->hittype != TERMINAL && 
	     (y->hit5->hittype == TERMINAL || y->hit3->hittype == TERMINAL)) {
    return -1;
  } else if ((x->hit5->hittype == TERMINAL || x->hit3->hittype == TERMINAL) &&
	     y->hit5->hittype != TERMINAL && y->hit3->hittype != TERMINAL) {
    return +1;
#endif

#if 0
  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else if (x->nmatches_posttrim > y->nmatches_posttrim) {
    return -1;
  } else if (y->nmatches_posttrim > x->nmatches_posttrim) {
    return +1;
#endif
  } else if (x->nchimera_novel < y->nchimera_novel) {
    return -1;
  } else if (y->nchimera_novel < x->nchimera_novel) {
    return +1;
  } else if (x->nchimera_known > y->nchimera_known) {
    return -1;
  } else if (y->nchimera_known > x->nchimera_known) {
    return +1;
  } else if (x->sense_consistent_p == true && y->sense_consistent_p == false) {
    return -1;
  } else if (x->sense_consistent_p == false && y->sense_consistent_p == true) {
    return +1;

#if 0
  } else if (x->indel_low < y->indel_low) {
    return -1;
  } else if (y->indel_low < x->indel_low) {
    return +1;
#endif

  } else {
    return 0;
  }
}


static int
hitpair_position_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;
  
  if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high < y->high) {
    return +1;
  } else if (y->high < x->high) {
    return -1;
  } else {
    return 0;
  }
}


#if 0
static bool
hitpair_equal (Stage3pair_T x, Stage3pair_T y) {
  if (x->dir != y->dir) {
    debug8(printf("=>F "));
    return false;		/* Different strands */
  } else if (x->hit5->low != y->hit5->low) {
    debug8(printf("=>F "));
    return false;
  } else if (x->hit5->high != y->hit5->high) {
    debug8(printf("=>F "));
    return false;
  } else if (x->hit3->low != y->hit3->low) {
    debug8(printf("=>F "));
    return false;
  } else if (x->hit3->high != y->hit3->high) {
    debug8(printf("=>F "));
    return false;
  } else {
    debug8(printf("=>T "));
    return true;
  }
}
#endif


#if 0
static bool
hitpair_overlap (Stage3pair_T x, Stage3pair_T y) {
#if 0
  if ((x->hit5->hittype == SPLICE || x->hit3->hittype == SPLICE) &&
      (y->hit5->hittype == SPLICE || y->hit3->hittype == SPLICE)) {
    /* Special case: pairs involving splices don't overlap */
    return false;
  }
#endif
  if (x->dir != y->dir) {
    return false;		/* Different strands */
  } else if (x->high < y->low) {
    return false;
  } else if (x->low > y->high) {
    return false;
  } else {
    return true;
  }
}
#endif


static bool
hitpair_subsumption (Stage3pair_T x, Stage3pair_T y) {
  if (x->dir != y->dir) {
    return false;		/* Different strands */
  } else if (x->low <= y->low && x->high >= y->high) {
    return true;
  } else if (y->low <= x->low && y->high >= x->high) {
    return true;
    
    /* Test each end of the pair.  Example: 1586..1512 and 1400..1468 should subsume 1586..1512 and 1564..1617 */
  } else if (x->hit5->low <= y->hit5->low && x->hit5->high >= y->hit5->high) {
    return true;
  } else if (y->hit5->low <= x->hit5->low && y->hit5->high >= x->hit5->high) {
    return true;

  } else if (x->hit3->low <= y->hit3->low && x->hit3->high >= y->hit3->high) {
    return true;
  } else if (y->hit3->low <= x->hit3->low && y->hit3->high >= x->hit3->high) {
    return true;

  } else {
    return false;
  }
}


static int
pair_matches_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;

  if (x->nmatches > y->nmatches) {
    return -1;
  } else if (x->nmatches < y->nmatches) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Stage3pair_sort_bymatches (List_T hits) {
  List_T sorted = NULL;
  Stage3pair_T *array;
  int n, i;

  n = List_length(hits);
  if (n == 0) {
    return (List_T) NULL;
  } else {
    array = (Stage3pair_T *) List_to_array(hits,NULL);
    List_free(&hits);

    qsort(array,n,sizeof(Stage3pair_T),pair_matches_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,(void *) array[i]);
    }
    FREE(array);

    return sorted;
  }
}


#if 0
static int
hitpair_distance_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;

  debug(printf("Comparing %u with %u\n",x->absdifflength,y->absdifflength));
  if (x->absdifflength < y->absdifflength) {
    return -1;
  } else if (x->absdifflength > y->absdifflength) {
    return +1;
  } else {
    return 0;
  }
}
#endif


#if 0
List_T
Stage3pair_sort_distance (List_T hitpairlist) {
#ifdef DEBUG
  Stage3pair_T hitpair;
  List_T p;
#endif
  List_T sorted = NULL;
  Stage3pair_T *hitpairs;
  int n, i;

  n = List_length(hitpairlist);
  if (n == 0) {
    return NULL;
  } else {
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    List_free(&hitpairlist);
  }

  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_distance_cmp);

  for (i = n-1; i >= 0; i--) {
    sorted = List_push(sorted,hitpairs[i]);
  }
  FREE(hitpairs);

  debug(
	for (p = sorted, i = 0; p != NULL; p = p->rest, i++) {
	  hitpair = (Stage3pair_T) p->first;
	  printf("  Final %d: %lu-%lu (dir = %d), insert length %u\n",
		 i,hitpair->low,hitpair->high,hitpair->dir,hitpair->insertlength);
	}
	);

  return sorted;
}
#endif


#if 0
List_T
Stage3pair_remove_duplicates_exact (List_T hitpairlist) {
  List_T unique = NULL;
  Stage3pair_T hitpair, *hitpairs;
  int n, i, j;
  bool *eliminate;

  n = List_length(hitpairlist);
  debug8(printf("Entered Stage3pair_remove_duplicates_exact with %d pairs\n",n));
  if (n == 0) {
    return NULL;
  } else {
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    List_free(&hitpairlist);
  }

  debug8(printf("Checking for exact duplicates\n"));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s, %s-%s): %p, %lu..%lu  %u..%u|%u..%u (dir = %d), nmatches: %d\n",
		  i,Pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
		  hittype_string(hitpair->hit3->hittype),hitpair,hitpair->low,hitpair->high,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hitpair->dir,hitpair->nmatches);
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && hitpair_equiv_cmp(hitpairs[j],hitpairs[i]) == 0) {
      debug8(printf("  %d is identical to %d => eliminating\n",j,i));
      eliminate[j] = true;
      j++;
    }
    i = j;
  }

  for (i = n-1; i >= 0; i--) {
    hitpair = hitpairs[i];
    if (eliminate[i] == false) {
      unique = List_push(unique,hitpair);
    } else {
      Stage3pair_free(&hitpair);
    }
  }

  FREE(hitpairs);
  FREE(eliminate);

  debug8(printf("Exited Stage3pair_remove_duplicates_exact with %d pairs\n",List_length(unique)));
  return unique;
}
#endif


static int
hitpair_goodness_cmp (bool *equalp, Stage3pair_T hitpair,
#ifdef DEBUG8
		      int k,
#endif
		      Stage3pair_T best_hitpair, bool finalp) {
  double prob1, prob2;

#if 0
  int hitpair_nmatches, best_hitpair_nmatches;
  int max_trim_left, max_trim_right;
  Stage3end_T hit5, besthit5, hit3, besthit3;

  if (hitpair->absdifflength_bingo_p < best_hitpair->absdifflength_bingo_p) {
    /* k is worse */
    debug8(printf(" => %d loses by absdifflength (bingo)\n",k));
    return -1;
  } else if (hitpair->absdifflength_bingo_p > best_hitpair->absdifflength_bingo_p) {
    /* k is better */
    debug8(printf(" => %d wins by absdifflength (bingo)\n",k));
    return +1;
  }
#endif

#ifdef PRE_RESOLVE_MULTIMAPPING
  if (TALLY_RATIO*Stage3pair_tally(hitpair) < Stage3pair_tally(best_hitpair)) {
    /* k is worse */
    debug8(printf(" => %d loses by tally\n",k));
    return -1;
  } else if (Stage3pair_tally(hitpair) > TALLY_RATIO*Stage3pair_tally(best_hitpair)) {
    /* k is better */
    debug8(printf(" => %d wins by tally\n",k));
    return +1;
  }
#endif

  *equalp = false;

#if 1
  if (finalp == true) {
    /* Skip */
  } else if (hitpair->hit5->hittype == TERMINAL || hitpair->hit3->hittype == TERMINAL ||
	     best_hitpair->hit5->hittype == TERMINAL || best_hitpair->hit3->hittype == TERMINAL) {
    /* Do not allow terminal to win or lose in pre-final stages */
    debug8(printf(" => %d ties by terminal\n",k));
    return 0;
  }
#endif

#if 0
  /* hitpair_nmatches = hitpair->nmatches; */
  /* best_hitpair_nmatches = best_hitpair->nmatches; */

  if (hitpair->hit5->hittype == TERMINAL || best_hitpair->hit5->hittype == TERMINAL ||
      hitpair->hit3->hittype == TERMINAL || best_hitpair->hit3->hittype == TERMINAL) {
    /* Skip: Don't use scores if terminal is involved */

#if 0
  } else if (hitpair->hit5->hittype == GMAP && best_hitpair->hit5->hittype == GMAP &&
	     hitpair->hit3->hittype == GMAP && best_hitpair->hit3->hittype == GMAP) {
    /* Dual GMAP alignments: Compare only in trimmed region */
    hit5 = hitpair->hit5;
    besthit5 = best_hitpair->hit5;
    max_trim_left = (hit5->trim_left > besthit5->trim_left) ? hit5->trim_left : besthit5->trim_left;
    max_trim_right = (hit5->trim_right > besthit5->trim_right) ? hit5->trim_right : besthit5->trim_right;
    hitpair_nmatches = Pair_array_nmatches_posttrim(hit5->pairarray,hit5->npairs,
						    /*pos5*/max_trim_left,/*pos3*/hit5->querylength_adj - max_trim_right);
    best_hitpair_nmatches = Pair_array_nmatches_posttrim(besthit5->pairarray,besthit5->npairs,
							 /*pos5*/max_trim_left,/*pos3*/besthit5->querylength_adj - max_trim_right);
    debug8(printf(" gmap/gmap on 5' end with trim %d left, %d right: %d versus %d",
		  max_trim_left,max_trim_right,hitpair_nmatches,best_hitpair_nmatches));

    hit3 = hitpair->hit3;
    besthit3 = best_hitpair->hit3;
    max_trim_left = (hit3->trim_left > besthit3->trim_left) ? hit3->trim_left : besthit3->trim_left;
    max_trim_right = (hit3->trim_right > besthit3->trim_right) ? hit3->trim_right : besthit3->trim_right;
    hitpair_nmatches += Pair_array_nmatches_posttrim(hit3->pairarray,hit3->npairs,
						     /*pos5*/max_trim_left,/*pos3*/hit3->querylength_adj - max_trim_right);
    best_hitpair_nmatches += Pair_array_nmatches_posttrim(besthit3->pairarray,besthit3->npairs,
							  /*pos5*/max_trim_left,/*pos3*/besthit3->querylength_adj - max_trim_right);
    debug8(printf(" gmap/gmap on 3' end with trim %d left, %d right: %d versus %d",
		  max_trim_left,max_trim_right,hitpair_nmatches,best_hitpair_nmatches));

  } else if (hitpair->hit5->hittype == GMAP && best_hitpair->hit5->hittype == GMAP) {
    /* 5' GMAP alignments: Compare only in trimmed region */
    hit5 = hitpair->hit5;
    besthit5 = best_hitpair->hit5;
    max_trim_left = (hit5->trim_left > besthit5->trim_left) ? hit5->trim_left : besthit5->trim_left;
    max_trim_right = (hit5->trim_right > besthit5->trim_right) ? hit5->trim_right : besthit5->trim_right;
    hitpair_nmatches = Pair_array_nmatches_posttrim(hit5->pairarray,hit5->npairs,
						    /*pos5*/max_trim_left,/*pos3*/hit5->querylength_adj - max_trim_right);
    best_hitpair_nmatches = Pair_array_nmatches_posttrim(besthit5->pairarray,besthit5->npairs,
							 /*pos5*/max_trim_left,/*pos3*/besthit5->querylength_adj - max_trim_right);
    debug8(printf(" gmap/gmap on 5' end with trim %d left, %d right: %d versus %d",
		  max_trim_left,max_trim_right,hitpair_nmatches,best_hitpair_nmatches));

    hitpair_nmatches += hitpair->hit3->nmatches;
    best_hitpair_nmatches += best_hitpair->hit3->nmatches;

  } else if (hitpair->hit3->hittype == GMAP && best_hitpair->hit3->hittype == GMAP) {
    /* 3' GMAP alignments: Compare only in trimmed region */
    hit3 = hitpair->hit3;
    besthit3 = best_hitpair->hit3;
    max_trim_left = (hit3->trim_left > besthit3->trim_left) ? hit3->trim_left : besthit3->trim_left;
    max_trim_right = (hit3->trim_right > besthit3->trim_right) ? hit3->trim_right : besthit3->trim_right;
    hitpair_nmatches = Pair_array_nmatches_posttrim(hit3->pairarray,hit3->npairs,
						     /*pos5*/max_trim_left,/*pos3*/hit3->querylength_adj - max_trim_right);
    best_hitpair_nmatches = Pair_array_nmatches_posttrim(besthit3->pairarray,besthit3->npairs,
							  /*pos5*/max_trim_left,/*pos3*/besthit3->querylength_adj - max_trim_right);
    debug8(printf(" gmap/gmap on 3' end with trim %d left, %d right: %d versus %d",
		  max_trim_left,max_trim_right,hitpair_nmatches,best_hitpair_nmatches));

    hitpair_nmatches += hitpair->hit5->nmatches;
    best_hitpair_nmatches += best_hitpair->hit5->nmatches;
#endif

  } else if (hitpair->hit5->nindels == 0 && best_hitpair->hit5->nindels == 0 &&
	     hitpair->hit3->nindels == 0 && best_hitpair->hit3->nindels == 0) {
    /* Skip: Use scores only if indel is involved */
  } else if (hitpair->score > best_hitpair->score) {
    /* k is worse */
    debug8(printf(" => %d loses by score\n",k));
    return -1;
  } else if (hitpair->score < best_hitpair->score) {
    /* k is better */
    debug8(printf(" => %d wins by score\n",k));
    return +1;
  }
#endif


#if 0
  if (hitpair->nmatches < best_hitpair->nmatches) {
    /* k is worse */
    debug8(printf(" => %d loses by nmatches\n",k));
    return -1;
  } else if (hitpair->nmatches > best_hitpair->nmatches) {
    /* k is better */
    debug8(printf(" => %d wins by nmatches\n",k));
    return +1;
  }
#endif

  if (hitpair->nmatches_posttrim < best_hitpair->nmatches_posttrim) {
    /* k is worse */
    debug8(printf(" => %d loses by nmatches_posttrim\n",k));
    return -1;
  } else if (hitpair->nmatches_posttrim > best_hitpair->nmatches_posttrim) {
    /* k is better */
    debug8(printf(" => %d wins by nmatches_posttrim\n",k));
    return +1;

  } else if (hitpair->nchimera_novel > best_hitpair->nchimera_novel) {
    /* k is worse */
    debug8(printf(" => %d loses by nchimera_novel\n",k));
    return -1;
  } else if (hitpair->nchimera_novel < best_hitpair->nchimera_novel) {
    /* k is better */
    debug8(printf(" => %d wins by nchimera_novel\n",k));
    return +1;

    /* Favoring nchimera_known helps before outerlength favors known
       splices over novel ones */
  } else if (hitpair->nchimera_known < best_hitpair->nchimera_known) {
    /* k is worse */
    debug8(printf(" => %d loses by nchimera_known\n",k));
    return -1;
  } else if (hitpair->nchimera_known > best_hitpair->nchimera_known) {
    /* k is better */
    debug8(printf(" => %d wins by nchimera_known\n",k));
    return +1;

#if 0
  } else if (hitpair->absdifflength < best_hitpair->absdifflength) {
    /* k is worse */
    debug8(printf(" => %d loses by absdifflength\n",k));
    return -1;
  } else if (hitpair->absdifflength > best_hitpair->absdifflength) {
    /* k is better */
    debug8(printf(" => %d wins by absdifflength\n",k));
    return +1;
#endif

  } else if (hitpair->hit5->hittype > best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype >= best_hitpair->hit3->hittype) {
    /* k is worse */
    debug8(printf(" => %d loses by hittype\n",k));
    return -1;

  } else if (hitpair->hit5->hittype >= best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype > best_hitpair->hit3->hittype) {
    /* k is worse */
    debug8(printf(" => %d loses by hittype\n",k));
    return -1;

  } else if (hitpair->hit5->hittype < best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype <= best_hitpair->hit3->hittype) {
    /* k is better */
    debug8(printf(" => %d wins by hittype\n",k));
    return +1;

  } else if (hitpair->hit5->hittype <= best_hitpair->hit5->hittype &&
	     hitpair->hit3->hittype < best_hitpair->hit3->hittype) {
    /* k is better */
    debug8(printf(" => %d wins by hittype\n",k));
    return +1;

  } else if (finalp == false) {
    debug8(printf("  => indistinguishable\n"));
    return 0;

#ifdef USE_ABSDIFFLENGTH
    /* If insert length is within deviation of expected pairlength, favor it */
  } else if (best_hitpair->absdifflength <= (Chrpos_T) pairlength_deviation &&
	     hitpair->absdifflength > (Chrpos_T) pairlength_deviation) {
    /* k is worse */
    debug8(printf(" => %d loses by absdifflength within deviation %d\n",k,pairlength_deviation));
    return -1;
  } else if (hitpair->absdifflength <= (Chrpos_T) pairlength_deviation &&
	     best_hitpair->absdifflength > (Chrpos_T) pairlength_deviation) {
    /* k is better */
    debug8(printf(" => %d wins by absdifflength within deviation %d\n",k,pairlength_deviation));
    return +1;
#endif

    /* Previously favored longer insert lengths to give more compact
       splices.  However, we now accept splices first that give
       expected pairlength */
  } else if (hitpair->insertlength_expected_sign == -1 && best_hitpair->insertlength_expected_sign == +1) {
    /* k is worse */
    debug8(printf(" => %d loses by insertlength_expected_sign\n",k));
    return -1;
  } else if (hitpair->insertlength_expected_sign == +1 && best_hitpair->insertlength_expected_sign == -1) {
    /* k is better */
    debug8(printf(" => %d wins by insertlength_expected_sign\n",k));
    return +1;

    /* Next we look at splice probability */
  } else {
    if (hitpair->hit5->hittype == SPLICE && best_hitpair->hit5->hittype == SPLICE &&
	hitpair->hit3->hittype == SPLICE && best_hitpair->hit3->hittype == SPLICE) {
      debug8(printf(" => dual splice"));
      prob1 = Substring_chimera_prob(hitpair->hit5->substring_donor) + Substring_chimera_prob(hitpair->hit5->substring_acceptor) +
	Substring_chimera_prob(hitpair->hit3->substring_donor) + Substring_chimera_prob(hitpair->hit3->substring_acceptor);
      prob2 = Substring_chimera_prob(best_hitpair->hit5->substring_donor) + Substring_chimera_prob(best_hitpair->hit5->substring_acceptor) +
	Substring_chimera_prob(best_hitpair->hit3->substring_donor) + Substring_chimera_prob(best_hitpair->hit3->substring_acceptor);
      if (prob1 + 0.3 < prob2) {
	/* k is worse */
	debug8(printf(" => %d loses by dual splice prob %f vs %f\n",k,prob1,prob2));
	return -1;
      } else if (prob1 > prob2 + 0.3) {
	/* k is better */
	debug8(printf(" => %d wins by dual splice prob %f vs %f\n",k,prob1,prob2));
	return +1;
      }

    } else if (hitpair->hit5->hittype == SPLICE && best_hitpair->hit5->hittype == SPLICE) {
      debug8(printf(" => splice on hit5"));
      prob1 = Substring_chimera_prob(hitpair->hit5->substring_donor) + Substring_chimera_prob(hitpair->hit5->substring_acceptor);
      prob2 = Substring_chimera_prob(best_hitpair->hit5->substring_donor) + Substring_chimera_prob(best_hitpair->hit5->substring_acceptor);
      if (prob1 + 0.3 < prob2) {
	/* k is worse */
	debug8(printf(" => %d loses by splice prob %f vs %f\n",k,prob1,prob2));
	return -1;
      } else if (prob1 > prob2 + 0.3) {
	/* k is better */
	debug8(printf(" => %d wins by splice prob %f vs %f\n",k,prob1,prob2));
	return +1;
      }

    } else if (hitpair->hit3->hittype == SPLICE && best_hitpair->hit3->hittype == SPLICE) {
      debug8(printf(" => splice on hit3"));
      prob1 = Substring_chimera_prob(hitpair->hit3->substring_donor) + Substring_chimera_prob(hitpair->hit3->substring_acceptor);
      prob2 = Substring_chimera_prob(best_hitpair->hit3->substring_donor) + Substring_chimera_prob(best_hitpair->hit3->substring_acceptor);
      if (prob1 + 0.3 < prob2) {
	/* k is worse */
	debug8(printf(" => %d loses by splice prob %f vs %f\n",k,prob1,prob2));
	return -1;
      } else if (prob1 > prob2 + 0.3) {
	/* k is better */
	debug8(printf(" => %d wins by splice prob %f vs %f\n",k,prob1,prob2));
	return +1;
      }
    }

    /* Overlapping ends worse than separate ends */
    if (hitpair->insertlength <= hitpair->hit5->querylength_adj + hitpair->hit3->querylength_adj &&
	best_hitpair->insertlength > best_hitpair->hit5->querylength_adj + best_hitpair->hit3->querylength_adj) {
      debug8(printf(" => %d loses by being overlapping\n",k));
      return -1;
    } else if (hitpair->insertlength > hitpair->hit5->querylength_adj + hitpair->hit3->querylength_adj &&
	       best_hitpair->insertlength <= best_hitpair->hit5->querylength_adj + best_hitpair->hit3->querylength_adj) {
      debug8(printf(" => %d wins by being separate\n",k));
      return +1;

      /* Next, favor shorter outerlengths to give more compact splices or closer pairs */
    } else if (hitpair->outerlength > best_hitpair->outerlength) {
      /* k is worse */
      debug8(printf(" => %d loses by outerlength\n",k));
      return -1;
    } else if (hitpair->outerlength < best_hitpair->outerlength) {
      /* k is better */
      debug8(printf(" => %d wins by outerlength\n",k));
      return +1;
      
    } else {
      if (hitpair->insertlength_expected_sign >= 0 && best_hitpair->insertlength_expected_sign >= 0) {
	/* Both insert lengths are short, so favor shorter insert length */
	debug8(printf(" => short insertlengths"));
	/* Favor shorter insert lengths */
	if (hitpair->insertlength > best_hitpair->insertlength) {
	  /* k is worse */
	  debug8(printf(" => %d loses by insertlength\n",k));
	  return -1;
	} else if (hitpair->insertlength < best_hitpair->insertlength) {
	  /* k is better */
	  debug8(printf(" => %d wins by insertlength\n",k));
	  return +1;
	}

      } else {
	/* Both insert lengths are long, so favor longer insert length to give more compact splices */
	debug8(printf(" => long insertlengths"));
	if (hitpair->insertlength < best_hitpair->insertlength) {
	  /* k is worse */
	  debug8(printf(" => %d loses by insertlength\n",k));
	  return -1;
	} else if (hitpair->insertlength > best_hitpair->insertlength) {
	  /* k is better */
	  debug8(printf(" => %d wins by insertlength\n",k));
	  return +1;
	}
      }

      debug8(printf("  => equal\n"));
      *equalp = true;
      return 0;
    }
  }
}


static bool
hitpair_bad_superstretch_p (Stage3pair_T hitpair_k, Stage3pair_T *hitpairs, int k, int j,
			    bool finalp) {
  int a;
  bool equalp;

  for (a = k+1; a <= j; a++) {
    if (hitpair_subsumption(hitpair_k,hitpairs[a]) == true) {
      debug8(printf("Testing %d because stretches over %d",k,a));
      if (hitpair_goodness_cmp(&equalp,hitpairs[a],
#ifdef DEBUG8
			       a,
#endif
			       hitpair_k,finalp) > 0 || equalp == true) {
	debug8(printf(" => eliminating\n"));
	return true;
      }
      debug8(printf("\n"));
    }
  }
  return false;
}


static List_T
pair_remove_overlaps (List_T hitpairlist, bool translocp, bool finalp) {
  List_T unique = NULL;
  Stage3pair_T best_hitpair, hitpair, *hitpairs, *prev;
  int cmp;
  int nkept, n, i, j, k, besti;
  bool *eliminate, equalp;
#ifdef PRE_RESOLVE_MULTIMAPPING
  long int best_tally;
#endif

  n = List_length(hitpairlist);
  debug8(printf("  Entering pair_remove_overlaps with %d pairs: %s\n",
		n,finalp == true ? "FINAL" : "not final"));

  if (n < 2) {
    debug8(printf("  Exiting pair_remove_overlaps with %d < 2 pairs\n",n));
    return hitpairlist;
  } else {
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    List_free(&hitpairlist);
  }

  /* Step 1.  Check for exact duplicates */
  debug8(printf("  Step 1. Checking for exact duplicates\n"));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s, %s-%s): %p, %lu..%lu  %u..%u|%u..%u (dir = %d), nmatches: %d (%d posttrim), indel_low %d and %d\n",
		  i,Pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
		  hittype_string(hitpair->hit3->hittype),hitpair,hitpair->low,hitpair->high,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hitpair->dir,hitpair->nmatches,hitpair->nmatches_posttrim,
		  hitpair->hit5->indel_low,hitpair->hit3->indel_low);
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    debug8(printf(" %d,%d",i,j));
    while (j < n && hitpair_equiv_cmp(hitpairs[j],hitpairs[i]) == 0) {
      debug8(printf("  %d is identical to %d => eliminating\n",j,i));
      eliminate[j] = true;
      j++;
    }
    i = j;
  }
  debug8(printf("\n"));

  nkept = 0;
  for (i = 0; i < n; i++) {
    if (eliminate[i] == false) {
      nkept++;
    }
  }
  if (nkept == 0) {
    /* All entries eliminated one another, so keep the first one */
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hitpairs;
  hitpairs = (Stage3pair_T *) CALLOC(nkept,sizeof(Stage3pair_T));

  for (i = 0, j = 0; i < n; i++) {
    hitpair = prev[i];
    if (eliminate[i] == false) {
      debug8(printf("  Keeping %lu..%lu  %u..%u|%u..%u, nmatches (trimmed) %d, score %d, (dir = %d)\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches,hitpair->score,hitpair->dir));
      hitpairs[j++] = hitpair;
    } else {
      debug8(printf("  Eliminating %lu..%lu  %u..%u|%u..%u, nmatches (trimmed) %d, score %d, (dir = %d)\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches,hitpair->score,hitpair->dir));
      Stage3pair_free(&hitpair);
    }
  }
    
  FREE(prev);

  if (translocp == false) {
    /* Step 2: Check for superstretches */
    n = nkept;
    debug8(printf("  Step 2.  Checking for superstretches among %d hitpairs within subsumption clusters\n",n));

    for (i = 0; i < n; i++) {
      eliminate[i] = false;
    }

    debug8(
	   for (i = 0; i < n; i++) {
	     hitpair = hitpairs[i];
	     printf("  Initial %d (%s, %s-%s): %p, %lu..%lu  %u..%u|%u..%u (dir = %d), score: %d, nmatches: %d (%d posttrim), nnovel: %d, nknown: %d, insertlength: %u, outerlength: %u\n",
		    i,Pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
		    hittype_string(hitpair->hit3->hittype),hitpair,hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->dir,hitpair->score,hitpair->nmatches,hitpair->nmatches_posttrim,
		    hitpair->nchimera_novel,hitpair->nchimera_known,
		    hitpair->insertlength,hitpair->outerlength);
	   }
	   );

    /* Find clusters */
    i = 0;
    while (i < n) {
      j = i;
      while (j+1 < n && hitpair_subsumption(hitpairs[i],hitpairs[j+1]) == true) {
	j = j+1;
      }

      if (j > i) {
	debug8(printf("Cluster from %d up through %d\n",i,j));

	/* Find bad superstretches */
	for (k = i; k <= j; k++) {
	  if (hitpair_bad_superstretch_p(hitpairs[k],hitpairs,k,j,finalp) == true) {
	    eliminate[k] = true;
	  }
	}
      }

      i = j+1;
    }

    nkept = 0;
    for (i = 0; i < n; i++) {
      if (eliminate[i] == false) {
	nkept++;
      }
    }
    if (nkept == 0) {
      /* All entries eliminated one another, so keep the first one */
      eliminate[0] = false;
      nkept = 1;
    }

    prev = hitpairs;
    hitpairs = (Stage3pair_T *) CALLOC(nkept,sizeof(Stage3pair_T));

    for (i = 0, j = 0; i < n; i++) {
      hitpair = prev[i];
      if (eliminate[i] == false) {
	debug8(printf("  Keeping %lu..%lu  %u..%u|%u..%u, nmatches (trimmed) %d (dir = %d)\n",
		      hitpair->low,hitpair->high,
		      hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		      hitpair->nmatches,hitpair->dir));
	hitpairs[j++] = hitpair;
      } else {
	debug8(printf("  Eliminating %lu..%lu  %u..%u|%u..%u, nmatches (trimmed) %d (dir = %d)\n",
		      hitpair->low,hitpair->high,
		      hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		      hitpair->nmatches,hitpair->dir));
	Stage3pair_free(&hitpair);
      }
    }

    FREE(prev);
  }


  if (translocp == false) {
    /* Step 3: Check for best within subsumption clusters */
    n = nkept;
    debug8(printf("  Step 3. Checking for best among %d hitpairs within subsumption clusters\n",n));

    for (i = 0; i < n; i++) {
      eliminate[i] = false;
    }
    /* qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp); -- No need since original order was kept */
  
    debug8(
	   for (i = 0; i < n; i++) {
	     hitpair = hitpairs[i];
	     printf("  Initial %d (%s, %s-%s): %p, %lu..%lu  %u..%u|%u..%u (dir = %d), score: %d, nmatches: %d (%d posttrim), nnovel: %d, nknown: %d, insertlength: %u, outerlength: %u\n",
		    i,Pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
		    hittype_string(hitpair->hit3->hittype),hitpair,hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->dir,hitpair->score,hitpair->nmatches,hitpair->nmatches_posttrim,
		    hitpair->nchimera_novel,hitpair->nchimera_known,
		    hitpair->insertlength,hitpair->outerlength);
	   }
	   );


    /* Find clusters from left */
    i = 0;
    while (i < n) {
      j = i;
      while (j+1 < n && hitpair_subsumption(hitpairs[i],hitpairs[j+1]) == true) {
	j = j+1;
      }

      if (j > i) {
	debug8(printf("Cluster from %d up through %d\n",i,j));

	best_hitpair = hitpairs[i];
	besti = i;
	debug8(printf("Assume best is %d\n",besti));

	for (k = i+1; k <= j; k++) {
	  cmp = hitpair_goodness_cmp(&equalp,hitpairs[k],
#ifdef DEBUG8
				     k,
#endif
				     best_hitpair,finalp);
	  debug8(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	  if (cmp > 0) {
	    best_hitpair = hitpairs[k];
	    besti = k;
	    debug8(printf("Best is now %d\n",besti));
	  }
	}

	for (k = i; k <= j; k++) {
	  if (k == besti) {
	    /* Skip */
	  } else if (hitpair_goodness_cmp(&equalp,hitpairs[k],
#ifdef DEBUG8
					  k,
#endif
					  best_hitpair,finalp) < 0 || equalp == true) {
	    debug8(printf("  Eliminating hitpair %d from left, because beaten by %d\n",k,besti));
	    eliminate[k] = true;
	  }
	}
      }
      
      i = j+1;
    }


    /* Find clusters starting from right */
    j = n - 1;
    while (j >= 0) {
      i = j;
      while (i-1 >= 0 && hitpair_subsumption(hitpairs[j],hitpairs[i-1]) == true) {
	i = i-1;
      }

      if (i < j) {
	debug8(printf("Cluster from %d down through %d\n",j,i));
	best_hitpair = hitpairs[i];
	besti = i;
	debug8(printf("Assume best is %d\n",besti));

	for (k = i+1; k <= j; k++) {
	  cmp = hitpair_goodness_cmp(&equalp,hitpairs[k],
#ifdef DEBUG8
				     k,
#endif
				     best_hitpair,finalp);
	  debug8(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	  if (cmp > 0) {
	    best_hitpair = hitpairs[k];
	    besti = k;
	    debug8(printf("Best is now %d\n",besti));
	  }
	}

	for (k = i; k <= j; k++) {
	  if (k == besti) {
	    /* Skip */
	  } else if (hitpair_goodness_cmp(&equalp,hitpairs[k],
#ifdef DEBUG8
					  k,
#endif
					  best_hitpair,finalp) < 0 || equalp == true) {
	    debug8(printf("  Eliminating hitpair %d from right, because beaten by %d\n",k,besti));
	    eliminate[k] = true;
	  }
	}
      }
      
      j = i-1;
    }


    nkept = 0;
    for (i = 0; i < n; i++) {
      if (eliminate[i] == false) {
	nkept++;
      }
    }
    if (nkept == 0) {
      eliminate[0] = false;
      nkept = 1;
    }

    prev = hitpairs;
    hitpairs = (Stage3pair_T *) CALLOC(nkept,sizeof(Stage3pair_T));

    for (i = 0, j = 0; i < n; i++) {
      hitpair = prev[i];
      if (eliminate[i] == false) {
	debug8(printf("  Keeping %lu..%lu  %u..%u|%u..%u, nmatches (trimmed) %d (dir = %d)\n",
		      hitpair->low,hitpair->high,
		      hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		      hitpair->nmatches,hitpair->dir));
	hitpairs[j++] = hitpair;
      } else {
	debug8(printf("  Eliminating %lu..%lu  %u..%u|%u..%u, nmatches (trimmed) %d (dir = %d)\n",
		      hitpair->low,hitpair->high,
		      hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		      hitpair->nmatches,hitpair->dir));
	Stage3pair_free(&hitpair);
      }
    }

    FREE(prev);
  }

  /* Step 4: Check for identity */
  n = nkept;
  debug8(printf("  Step 4.  Checking for duplicates among %d hitpairs by identity\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  /* qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp); -- No need since original order was kept */

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s, %s-%s): %p, %lu..%lu  %u..%u|%u..%u (dir = %d), score: %d, nmatches: %d (%d posttrim), insertlength: %u, outerlength: %u\n",
		  i,Pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
		  hittype_string(hitpair->hit3->hittype),hitpair,hitpair->low,hitpair->high,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hitpair->dir,hitpair->score,hitpair->nmatches,hitpair->nmatches_posttrim,
		  hitpair->insertlength,hitpair->outerlength);
	 }
	 );

#if 0
  i = 0;
  while (i < n) {
    debug8(printf("Looking at %d with score %d, absdifflength %d, and outerlength %u\n",
		  i,hitpairs[i]->score,hitpairs[i]->absdifflength,hitpairs[i]->outerlength));
    besti = i;
    best_score = hitpairs[i]->score;
    best_absdifflength = hitpairs[i]->absdifflength;
    best_outerlength = hitpairs[i]->outerlength;

    j = i+1;
    while (j < n && hitpair_overlap(hitpairs[j],hitpairs[j-1]) == true) {
      debug8(printf("  %d overlaps with %d with score %d, absdifflength %d, and outerlength %u",
		    i,j,hitpairs[j]->score,hitpairs[j]->outerlength));
      if (hitpairs[j]->score < best_score) {
	debug8(printf(" => best by score\n"));
	best_score = hitpairs[j]->score;
	best_absdifflength = hitpairs[j]->absdifflength;
	best_outerlength = hitpairs[j]->outerlength;
	besti = j;
      } else if (hitpairs[j]->score > best_score) {
	/* Skip */
      } else if (hitpairs[j]->absdifflength < best_absdifflength) {
	debug8(printf(" => best by absdifflength\n"));
	best_absdifflength = hitpairs[j]->absdifflength;
	best_outerlength = hitpairs[j]->outerlength;
	besti = j;
      } else if (hitpairs[j]->absdifflength > best_absdifflength) {
	/* Skip */
      } else if (hitpairs[j]->outerlength < best_outerlength) {
	debug8(printf(" => best by outerlength\n"));
	best_outerlength = hitpairs[j]->outerlength;
	besti = j;
      }
	
      debug8(printf("\n"));
      j++;
    }

    i = j;
  }
#else
  i = 0;
  while (i < n) {
    debug8(printf("Looking at %d with score %d, insert length %d, and outerlength %u\n",
		  i,hitpairs[i]->score,hitpairs[i]->insertlength,hitpairs[i]->outerlength));
    j = i+1;
    while (j < n && hitpair_equiv_cmp(hitpairs[j],hitpairs[i]) == 0) {
      debug8(printf("  %d equal to %d\n",j,i));
      eliminate[j] = true;
      j++;
    }

    i = j;
  }
#endif

  for (i = n-1; i >= 0; i--) {
    hitpair = hitpairs[i];
    if (eliminate[i] == false) {
      unique = List_push(unique,hitpair);
    } else {
      Stage3pair_free(&hitpair);
    }
  }

  FREE(hitpairs);
  FREE(eliminate);


#ifdef PRE_RESOLVE_MULTIMAPPING
  if (use_tally_p == true && tally_iit != NULL) {
    if ((n = List_length(unique)) > 1) {
      hitpairs = (Stage3pair_T *) List_to_array(unique,NULL);
      List_free(&unique);

      best_tally = 0;
      for (i = 0; i < n; i++) {
	if (hitpairs[i]->tally < 0) {
	  hitpairs[i]->tally = Stage3end_compute_tally(hitpairs[i]->hit5) + Stage3end_compute_tally(hitpairs[i]->hit3);
	}
	if (hitpairs[i]->tally > best_tally) {
	  best_tally = hitpairs[i]->tally;
	}
      }

      unique = (List_T) NULL;
      for (i = 0; i < n; i++) {
	if (hitpairs[i]->tally < best_tally) {
	  Stage3pair_free(&(hitpairs[i]));
	} else {
	  unique = List_push(unique,(void *) hitpairs[i]);
	}
      }

      FREE(hitpairs);
    }
  }
#endif

  debug8(printf("  Exited pair_remove_overlaps with %d pairs\n",List_length(unique)));
  return unique;
}


List_T
Stage3pair_remove_overlaps (List_T hitpairlist, bool translocp, bool finalp) {
  List_T unique_separate, unique_overlapping,
    separate = NULL, overlapping = NULL, p;
  Stage3pair_T hitpair;

  List_T indep_overlapping = NULL;
  Stage3pair_T *array_separate, *array_overlapping;
  Stage3pair_T hitpair_overlapping;
  Univcoord_T low, high;
  bool subsumedp, equalp;
  int n_separate, n_overlapping, i, j;


  for (p = hitpairlist; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    if (hitpair->insertlength <= hitpair->hit5->querylength_adj + hitpair->hit3->querylength_adj) {
      overlapping = List_push(overlapping,(void *) hitpair);
    } else {
      separate = List_push(separate,(void *) hitpair);
    }
  }
  List_free(&hitpairlist);

  debug8(printf("Calling Stage3pair_remove_overlaps for separate pair ends\n"));
  unique_separate = pair_remove_overlaps(separate,translocp,finalp);

  debug8(printf("Calling Stage3pair_remove_overlaps for overlapping pair ends\n"));
  unique_overlapping = pair_remove_overlaps(overlapping,translocp,finalp);
  
  if (unique_overlapping == NULL) {
    return unique_separate;
  } else if (unique_separate == NULL) {
    return unique_overlapping;
  } else {
    debug8(printf("Have both overlapping and separate\n"));
    n_overlapping = List_length(unique_overlapping);
    array_overlapping = (Stage3pair_T *) List_to_array(unique_overlapping,NULL);
    List_free(&unique_overlapping);
    n_separate = List_length(unique_separate);
    array_separate = (Stage3pair_T *) List_to_array(unique_separate,NULL);
    /* List_free(&unique_separate); -- save for final result */

    qsort(array_overlapping,n_overlapping,sizeof(Stage3pair_T),hitpair_position_cmp);
    qsort(array_separate,n_separate,sizeof(Stage3pair_T),hitpair_position_cmp);

    i = j = 0;
    for (i = 0; i < n_overlapping; i++) {
      hitpair_overlapping = array_overlapping[i];
      low = hitpair_overlapping->low;
      high = hitpair_overlapping->high;
      while (j >= 0 && array_separate[j]->high >= low) {
	j--;
      }
      j += 1;

      subsumedp = false;
      while (j < n_separate && subsumedp == false && array_separate[j]->low <= high) {
	if (hitpair_goodness_cmp(&equalp,array_separate[j],
#ifdef DEBUG8
				 j,
#endif
				 hitpair_overlapping,finalp) > 0) {
	  debug8(printf("separate pair %d better than overlapping pair %d\n",j,i));
	  subsumedp = hitpair_subsumption(array_separate[j],hitpair_overlapping);
	  debug8(printf("  checking if separate pair %d subsumes overlapping pair %d => %d\n",
			j,i,subsumedp));
	}
	j++;
      }
      j -= 1;

      if (subsumedp == true) {
	Stage3pair_free(&hitpair_overlapping);
      } else {
	indep_overlapping = List_push(indep_overlapping,(void *) hitpair_overlapping);
      }
    }

    FREE(array_separate);
    FREE(array_overlapping);

    return List_append(unique_separate,indep_overlapping);
  }
}


List_T
Stage3pair_resolve_multimapping (List_T hitpairs) {
  List_T resolve1, resolve2, resolve3, p;
  Stage3pair_T hitpair;

  Overlap_T best_overlap;
  long int best_tally;
  double tally_threshold;
  bool runlengthp;


  if (List_length(hitpairs) <= 1) {
    return hitpairs;
  }

  if (genes_iit == NULL) {
    resolve1 = hitpairs;
  } else {
    best_overlap = NO_KNOWN_GENE;
    for (p = hitpairs; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if ((hitpair->gene_overlap = Stage3pair_gene_overlap(hitpair)) > best_overlap) {
	best_overlap = hitpair->gene_overlap;
      }
    }
    if (best_overlap == NO_KNOWN_GENE) {
      resolve1 = hitpairs;
    } else {
      resolve1 = (List_T) NULL;
      for (p = hitpairs; p != NULL; p = p->rest) {
	hitpair = (Stage3pair_T) p->first;
	if (hitpair->gene_overlap < best_overlap) {
	  Stage3pair_free(&hitpair);
	} else {
	  resolve1 = List_push(resolve1,(void *) hitpair);
	}
      }
      List_free(&hitpairs);
    }
  }
      
  if (List_length(resolve1) <= 1) {
    return resolve1;
  }

  if (tally_iit == NULL) {
    resolve2 = resolve1;
  } else {
    best_tally = 0L;
    for (p = resolve1; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if ((hitpair->tally = Stage3end_compute_tally(hitpair->hit5) + Stage3end_compute_tally(hitpair->hit3)) > best_tally) {
	best_tally = hitpair->tally;
      }
    }
    if (best_tally == 0L) {
      resolve2 = resolve1;
    } else {
      resolve2 = (List_T) NULL;
#ifdef USE_TALLY_RATIO
      tally_threshold = (double) best_tally / TALLY_RATIO;
#else
      tally_threshold = 1.0;
#endif
      for (p = resolve1; p != NULL; p = p->rest) {
	hitpair = (Stage3pair_T) p->first;
	if ((double) hitpair->tally < tally_threshold) {
	  Stage3pair_free(&hitpair);
	} else {
	  resolve2 = List_push(resolve2,(void *) hitpair);
	}
      }
      List_free(&resolve1);
    }
  }


  if (List_length(resolve2) <= 1) {
    return resolve2;
  }

  if (runlength_iit == NULL) {
    resolve3 = resolve2;
  } else {
    runlengthp = false;
    for (p = resolve2; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (Stage3end_runlength_p(hitpair->hit5) == true || Stage3end_runlength_p(hitpair->hit3) == true) {
	runlengthp = true;
      }
    }
    if (runlengthp == false) {
      resolve3 = resolve2;
    } else {
      resolve3 = (List_T) NULL;
      for (p = resolve2; p != NULL; p = p->rest) {
	hitpair = (Stage3pair_T) p->first;
	if (Stage3end_runlength_p(hitpair->hit5) == false && Stage3end_runlength_p(hitpair->hit3) == false) {
	  Stage3pair_free(&hitpair);
	} else {
	  resolve3 = List_push(resolve3,(void *) hitpair);
	}
      }
      List_free(&resolve2);
    }
  }


  return resolve3;
}



Stage3pair_T *
Stage3pair_eval_and_sort (int *npaths, int *first_absmq, int *second_absmq,
			  Stage3pair_T *stage3pairarray, int maxpaths,
			  Shortread_T queryseq5, Shortread_T queryseq3,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev, 
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev, 
			  Genome_T genome, char *quality_string_5, char *quality_string_3) {
  char *query5, *query3;
  float maxlik, loglik;

  float total, q;
  int mapq_score;
  bool non_terminal_5p, non_terminal_3p;

  int compute_npaths;
  int i;

  if (*npaths == 0) {
    /* Skip */
    *first_absmq = 0;
    *second_absmq = 0;

  } else if (*npaths == 1) {
    stage3pairarray[0]->mapq_loglik = MAPQ_MAXIMUM_SCORE;
    stage3pairarray[0]->mapq_score = MAPQ_max_quality_score(quality_string_5,Shortread_fulllength(queryseq5));
    if ((mapq_score = MAPQ_max_quality_score(quality_string_3,Shortread_fulllength(queryseq3))) > stage3pairarray[0]->mapq_score) {
      stage3pairarray[0]->mapq_score = mapq_score;
    }
    stage3pairarray[0]->absmq_score = MAPQ_MAXIMUM_SCORE;

    query5 = Shortread_fullpointer_uc(queryseq5);
    query3 = Shortread_fullpointer_uc(queryseq3);

    assert(stage3pairarray[0]->private5p == true);
    assert(stage3pairarray[0]->private3p == true);
    Stage3end_display_prep(stage3pairarray[0]->hit5,query5,query5_compress_fwd,query5_compress_rev,
			   genome);
    Stage3end_display_prep(stage3pairarray[0]->hit3,query3,query3_compress_fwd,query3_compress_rev,
			   genome);

    *first_absmq = stage3pairarray[0]->absmq_score;
    *second_absmq = 0;

  } else {
    /* Determine whether to trim terminal ends */
    non_terminal_5p = non_terminal_3p = false;
    for (i = 0; i < *npaths; i++) {
      if (stage3pairarray[i]->hit5->hittype != TERMINAL) {
	non_terminal_5p = true;
      }
      if (stage3pairarray[i]->hit3->hittype != TERMINAL) {
	non_terminal_3p = true;
      }
    }

    /* Compute mapq_loglik */
    for (i = 0; i < *npaths; i++) {
      stage3pairarray[i]->mapq_loglik =
	Stage3end_compute_mapq(stage3pairarray[i]->hit5,query5_compress_fwd,query5_compress_rev,
			       quality_string_5,/*trim_terminals_p*/non_terminal_5p ? false : true);
      stage3pairarray[i]->mapq_loglik +=
      Stage3end_compute_mapq(stage3pairarray[i]->hit3,query3_compress_fwd,query3_compress_rev,
			     quality_string_3,/*trim_terminals_p*/non_terminal_3p ? false : true);
    }

    /* Sort by nmatches, then mapq, and then insert length.  Enforce monotonicity. */
    qsort(stage3pairarray,*npaths,sizeof(Stage3pair_T),Stage3pair_output_cmp);
    for (i = *npaths - 1; i > 0; i--) {
      if (stage3pairarray[i-1]->mapq_loglik < stage3pairarray[i]->mapq_loglik) {
	stage3pairarray[i-1]->mapq_loglik = stage3pairarray[i]->mapq_loglik;
      }
    }
    maxlik = stage3pairarray[0]->mapq_loglik;

    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < *npaths; i++) {
      stage3pairarray[i]->mapq_loglik -= maxlik;
    }

    /* Save on computation if possible */
    if (*npaths < maxpaths) {
      compute_npaths = *npaths;
    } else {
      compute_npaths = maxpaths;
    }
    if (compute_npaths < 2) {
      compute_npaths = 2;
    }

    /* Compute absolute mapq */
    for (i = 0; i < compute_npaths; i++) {
      loglik = stage3pairarray[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      stage3pairarray[i]->absmq_score = rint(loglik);
    }
    *first_absmq = stage3pairarray[0]->absmq_score;
    *second_absmq = stage3pairarray[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < *npaths; i++) {
      total += (stage3pairarray[i]->mapq_loglik = fasterexp(stage3pairarray[i]->mapq_loglik));
    }

    /* Prepare for display */
    query5 = Shortread_fullpointer_uc(queryseq5);
    query3 = Shortread_fullpointer_uc(queryseq3);
    for (i = 0; i < compute_npaths; i++) {
      assert(stage3pairarray[i]->private5p == true);
      assert(stage3pairarray[i]->private3p == true);
      Stage3end_display_prep(stage3pairarray[i]->hit5,query5,query5_compress_fwd,query5_compress_rev,
			     genome);
      Stage3end_display_prep(stage3pairarray[i]->hit3,query3,query3_compress_fwd,query3_compress_rev,
			     genome);
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < compute_npaths; i++) {
      stage3pairarray[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < compute_npaths; i++) {
      if ((q = 1.0 - stage3pairarray[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	stage3pairarray[i]->mapq_score = 96;
      } else {
	stage3pairarray[i]->mapq_score = rint(-10.0 * log10(q));
      }
    }

#if 0
    /* Apply filtering for mapq unique -- currently not used since mapq_unique_score is high */
    if (stage3pairarray[0]->mapq_score >= mapq_unique_score &&
	stage3pairarray[1]->mapq_score < mapq_unique_score) {
      for (i = 1; i < *npaths; i++) {
	Stage3pair_free(&(stage3pairarray[i]));
      }
      *npaths = 1;
    }
#endif
  }

  return stage3pairarray;
}


List_T
Stage3pair_remove_excess_terminals (List_T hitpairlist) {
  List_T cleaned = NULL, p;
  Stage3pair_T hitpair;
  int splice5p = false, splice3p = false, terminal5p = false, terminal3p = false;
  int best_splice5_score = 1000000, best_splice3_score = 1000000;

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->hit5->hittype == SPLICE || hitpair->hit5->hittype == SHORTEXON) {
      splice5p = true;
      if (hitpair->hit5->score < best_splice5_score) {
	best_splice5_score = hitpair->hit5->score;
      }
    } else if (hitpair->hit5->hittype == TERMINAL) {
      terminal5p = true;
    }

    if (hitpair->hit3->hittype == SPLICE || hitpair->hit3->hittype == SHORTEXON) {
      splice3p = true;
      if (hitpair->hit3->score < best_splice3_score) {
	best_splice3_score = hitpair->hit3->score;
      }
    } else if (hitpair->hit3->hittype == TERMINAL) {
      terminal3p = true;
    }
  }

  if ((splice5p == true && terminal5p == true) ||
      (splice3p == true && terminal3p == true)) {
    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->hit5->hittype == TERMINAL && splice5p == true && hitpair->hit5->score >= best_splice5_score) {
	Stage3pair_free(&hitpair);
      } else if (hitpair->hit3->hittype == TERMINAL && splice3p == true && hitpair->hit3->score >= best_splice3_score) {
	Stage3pair_free(&hitpair);
      } else {
	cleaned = List_push(cleaned,hitpair);
      }
    }

    List_free(&hitpairlist);
    return cleaned;

  } else {
    return hitpairlist;
  }

}



/* terminal alignments need to win on nmatches */
static List_T
Stage3pair_optimal_score_aux (bool *eliminatedp, List_T hitpairlist, int cutoff_level, int suboptimal_mismatches,
			      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			      bool keep_gmap_p, bool finalp) {
  List_T optimal = NULL, p;
  Stage3pair_T hitpair;
  T hit5, hit3;
  int cutoff_level_5, cutoff_level_3, score;
  int n;
  int minscore5 = MAX_READLENGTH, minscore3 = MAX_READLENGTH, minscore = MAX_READLENGTH + MAX_READLENGTH;
  /* int max_nmatches = 0, max_nmatches_posttrim, minscore = MAX_READLENGTH + MAX_READLENGTH; */
#ifdef USE_OPTIMAL_SCORE_BINGO
  int minscore_bingo = MAX_READLENGTH + MAX_READLENGTH;
#endif
  int trim_left_5, trim_right_5, trim_left_3, trim_right_3;
  int min_trim_left_5 = MAX_READLENGTH, min_trim_right_5 = MAX_READLENGTH,
    min_trim_left_3 = MAX_READLENGTH, min_trim_right_3 = MAX_READLENGTH;
  int max_trim_left_terminal_5 = 0, max_trim_right_terminal_5 = 0, 
    max_trim_left_terminal_3 = 0, max_trim_right_terminal_3 = 0;
  int nindelbreaks;
  bool non_double_terminal_p = false, non_terminal_5p = false, non_terminal_3p = false;

#ifdef TRANSLOC_SPECIAL
  bool non_translocation_p = false;
#endif


  *eliminatedp = false;

  n = List_length(hitpairlist);
  debug6(printf("\nEntered Stage3pair_optimal_score with %d hitpairs: %s\n",
		n,finalp == true ? "FINAL" : "not final"));
  
  if (n <= 1) {
    return hitpairlist;
  }

  p = hitpairlist;
  while (non_double_terminal_p == false && p != NULL) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->hit5->hittype != TERMINAL) {
      non_double_terminal_p = true;
    } else if (hitpair->hit3->hittype != TERMINAL) {
      non_double_terminal_p = true;
    }
    p = p->rest;
  }
  debug6(printf("non_double_terminal_p: %d\n",non_double_terminal_p));

  p = hitpairlist;
  while (non_terminal_5p == false && p != NULL) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->hit5->hittype != TERMINAL) {
      non_terminal_5p = true;
    }
    p = p->rest;
  }
  debug6(printf("non_terminal_5p: %d\n",non_terminal_5p));

  p = hitpairlist;
  while (non_terminal_3p == false && p != NULL) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->hit3->hittype != TERMINAL) {
      non_terminal_3p = true;
    }
    p = p->rest;
  }
  debug6(printf("non_terminal_3p: %d\n",non_terminal_3p));



  /* Use eventrim for comparing alignments */
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    hit5 = hitpair->hit5;
    hit3 = hitpair->hit3;

    debug6(printf("hit5 %lu..%lu type %s, trim_left: %d%s, trim_right %d%s.  hit3 %lu..%lu type %s, trim_left %d%s, trim_right %d%s.\n",
		  hit5->genomicstart,hit5->genomicend,hittype_string(hit5->hittype),
		  hit5->trim_left,hit5->trim_left_splicep ? " (splice)" : "",
		  hit5->trim_right,hit5->trim_right_splicep ? " (splice)" : "",
		  hit3->genomicstart,hit3->genomicend,hittype_string(hit3->hittype),
		  hit3->trim_left,hit3->trim_left_splicep ? " (splice)" : "",
		  hit3->trim_right,hit3->trim_right_splicep ? " (splice)" : ""));
    if (hit5->hittype == TERMINAL) {
      /* Don't allow terminals to set trims */
#if 0
      if (hit5->trim_left > max_trim_left_terminal_5) {
	max_trim_left_terminal_5 = hit5->trim_left;
      }
      if (hit5->trim_right > max_trim_right_terminal_5) {
	max_trim_right_terminal_5 = hit5->trim_right;
      }
#endif

    } else if ((hit5->hittype == INSERTION || hit5->hittype == DELETION) &&
	       (hit5->indel_pos < 15 || hit5->indel_pos > hit5->querylength_adj - 15)) {
      /* Don't allow end indels to set trims */

    } else {
      if (hit5->trim_left_splicep == true) {
	if (hit5->trim_left > max_trim_left_terminal_5) {
	  max_trim_left_terminal_5 = hit5->trim_left;
	}
      } else if (hit5->trim_left < min_trim_left_5) {
	min_trim_left_5 = hit5->trim_left;
      }
      if (hit5->trim_right_splicep == true) {
	if (hit5->trim_right > max_trim_right_terminal_5) {
	  max_trim_right_terminal_5 = hit5->trim_right;
	}
      } else if (hit5->trim_right < min_trim_right_5) {
	min_trim_right_5 = hit5->trim_right;
      }
    }

    if (hit3->hittype == TERMINAL) {
      /* Don't allow terminals to set trims */
#if 0
      if (hit3->trim_left > max_trim_left_terminal_3) {
	max_trim_left_terminal_3 = hit3->trim_left;
      }
      if (hit3->trim_right > max_trim_right_terminal_3) {
	max_trim_right_terminal_3 = hit3->trim_right;
      }
#endif

    } else if ((hit3->hittype == INSERTION || hit3->hittype == DELETION) &&
	       (hit3->indel_pos < 15 || hit3->indel_pos > hit3->querylength_adj - 15)) {
      /* Don't allow end indels to set trims */

    } else {
      if (hit3->trim_left_splicep == true) {
	if (hit3->trim_left > max_trim_left_terminal_3) {
	  max_trim_left_terminal_3 = hit3->trim_left;
	}
      } else if (hit3->trim_left < min_trim_left_3) {
	min_trim_left_3 = hit3->trim_left;
      }
      if (hit3->trim_right_splicep == true) {
	if (hit3->trim_right > max_trim_right_terminal_3) {
	  max_trim_right_terminal_3 = hit3->trim_right;
	}
      } else if (hit3->trim_right < min_trim_right_3) {
	min_trim_right_3 = hit3->trim_right;
      }
    }
  }

  if (min_trim_left_5 == MAX_READLENGTH) {
    trim_left_5 = max_trim_left_terminal_5;
  } else {
    trim_left_5 = (max_trim_left_terminal_5 > min_trim_left_5) ? max_trim_left_terminal_5 : min_trim_left_5;
  }
  if (min_trim_right_5 == MAX_READLENGTH) {
    trim_right_5 = max_trim_right_terminal_5;
  } else {
    trim_right_5 = (max_trim_right_terminal_5 > min_trim_right_5) ? max_trim_right_terminal_5 : min_trim_right_5;
  }

  if (min_trim_left_3 == MAX_READLENGTH) {
    trim_left_3 = max_trim_left_terminal_3;
  } else {
    trim_left_3 = (max_trim_left_terminal_3 > min_trim_left_3) ? max_trim_left_terminal_3 : min_trim_left_3;
  }
  if (min_trim_right_3 == MAX_READLENGTH) {
    trim_right_3 = max_trim_right_terminal_3;
  } else {
    trim_right_3 = (max_trim_right_terminal_3 > min_trim_right_3) ? max_trim_right_terminal_3 : min_trim_right_3;
  }

  debug6(printf("non-terminals: hit5 min_trim_left: %d, min_trim_right %d\n",
		min_trim_left_5,min_trim_right_5));
  debug6(printf("prefinal-terminals: hit5 max_trim_left: %d, max_trim_right %d\n",
		max_trim_left_terminal_5,max_trim_right_terminal_5));
  debug6(printf("overall: trim_left %d, trim_right %d\n",trim_left_5,trim_right_5));

  debug6(printf("non-terminals: hit3 min_trim_left: %d, min_trim_right %d\n",
		min_trim_left_3,min_trim_right_3));
  debug6(printf("prefinal-terminals: hit3 max_trim_left: %d, max_trim_right %d\n",
		max_trim_left_terminal_3,max_trim_right_terminal_3));
  debug6(printf("overall: trim_left %d, trim_right %d\n",trim_left_3,trim_right_3));


  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    hit5 = hitpair->hit5;
    hit3 = hitpair->hit3;

    if (hit5->hittype == TERMINAL && non_double_terminal_p == true && finalp == false) {
      /* Ignore */
      hit5->score_eventrim = 0;
    } else if (hit5->hittype == GMAP) {
      hit5->score_eventrim = 0;  /* was hit5->penalties */
      debug6(printf("score 5' GMAP:"));
#if 0
      if (Stage3end_bad_stretch_p(hit5,query5_compress_fwd,query5_compress_rev) == true) {
	hit5->score_eventrim += 2;
	debug6(printf("  bad stretch 2."));
      }
#endif

      if (0 && hit5->trim_left <= 8) {
	/* Ignore small trims */
      } else if (hit5->trim_left > trim_left_5) {
	hit5->score_eventrim += hit5->trim_left - trim_left_5;
	debug6(printf("  add trim left (%d - %d).",hit5->trim_left,trim_left_5));
      }
      if (0 && hit5->trim_right <= 8) {
	/* Ignore small trims */
      } else if (hit5->trim_right > trim_right_5) {
	hit5->score_eventrim += hit5->trim_right - trim_right_5;
	debug6(printf("  add trim right (%d - %d).",hit5->trim_right,trim_right_5));
      }

      hit5->score_eventrim += Pair_nmismatches_region(&nindelbreaks,hit5->pairarray,hit5->npairs,
						      trim_left_5,trim_right_5,hit5->querylength_adj);
      debug6(printf("  add nmismatches %d.",Pair_nmismatches_region(&nindelbreaks,hit5->pairarray,hit5->npairs,
								    trim_left_5,trim_right_5,hit5->querylength_adj)));
      hit5->score_eventrim += indel_penalty_middle * nindelbreaks;
      hit5->score_eventrim += hit5->amb_nmatches_start / ambig_end_interval;
      debug6(printf("  add amb start %d/%d.",hit5->amb_nmatches_start,ambig_end_interval));
      hit5->score_eventrim += hit5->amb_nmatches_end / ambig_end_interval;
      debug6(printf("  add amb end %d/%d.",hit5->amb_nmatches_end,ambig_end_interval));
      debug6(printf("  RESULT: %d\n",hit5->score_eventrim));
      
    } else {
      debug6(printf("score 5' OTHER:"));
      hit5->score_eventrim = hit5->penalties;
      debug6(printf("  penalties %d.",hit5->penalties));

      hit5->score_eventrim += Substring_count_mismatches_region(hit5->substring0,trim_left_5,trim_right_5,
								query5_compress_fwd,query5_compress_rev);
      debug6(printf("  substring 0 %d.",Substring_count_mismatches_region(hit5->substring0,trim_left_5,trim_right_5,
									  query5_compress_fwd,query5_compress_rev)));

      hit5->score_eventrim += Substring_count_mismatches_region(hit5->substring1,trim_left_5,trim_right_5,
							       query5_compress_fwd,query5_compress_rev);
      debug6(printf("  substring 1 %d.",Substring_count_mismatches_region(hit5->substring1,trim_left_5,trim_right_5,
									  query5_compress_fwd,query5_compress_rev)));

      hit5->score_eventrim += Substring_count_mismatches_region(hit5->substring2,trim_left_5,trim_right_5,
							       query5_compress_fwd,query5_compress_rev);
      debug6(printf("  substring 2 %d.",Substring_count_mismatches_region(hit5->substring2,trim_left_5,trim_right_5,
									  query5_compress_fwd,query5_compress_rev)));

      if (hit5->hittype == INSERTION || hit5->hittype == DELETION) {
	hit5->score_eventrim += indel_penalty_middle;
	debug6(printf("  add indel %d.",indel_penalty_middle));
      }
      debug6(printf("  RESULT: %d\n",hit5->score_eventrim));
    }

    if (hit3->hittype == TERMINAL && non_double_terminal_p == true && finalp == false) {
      /* Ignore */
      hit3->score_eventrim = 0;
    } else if (hit3->hittype == GMAP) {
      hit3->score_eventrim = 0;  /* was hit3->penalties */
      debug6(printf("score 3' GMAP:"));
#if 0
      if (Stage3end_bad_stretch_p(hit3,query3_compress_fwd,query3_compress_rev) == true) {
	hit3->score_eventrim += 2;
	debug6(printf("  bad stretch 2."));
      }
#endif

      if (0 && hit3->trim_left <= 8) {
	/* Ignore small trims */
      } else if (hit3->trim_left > trim_left_3) {
	hit3->score_eventrim += hit3->trim_left - trim_left_3;
	debug6(printf("  add trim left (%d - %d).",hit3->trim_left,trim_left_3));
      }
      if (0 && hit3->trim_right <= 8) {
	/* Ignore small trims */
      } else if (hit3->trim_right > trim_right_3) {
	hit3->score_eventrim += hit3->trim_right - trim_right_3;
	debug6(printf("  add trim right (%d - %d).",hit3->trim_right,trim_right_3));
      }

      hit3->score_eventrim += Pair_nmismatches_region(&nindelbreaks,hit3->pairarray,hit3->npairs,
						      trim_left_3,trim_right_3,hit3->querylength_adj);
      debug6(printf("  add nmismatches %d.",Pair_nmismatches_region(&nindelbreaks,hit3->pairarray,hit3->npairs,
								    trim_left_3,trim_right_3,hit3->querylength_adj)));
      hit3->score_eventrim += indel_penalty_middle * nindelbreaks;
      hit3->score_eventrim += hit3->amb_nmatches_start / ambig_end_interval;
      debug6(printf("  add amb start %d/%d.",hit3->amb_nmatches_start,ambig_end_interval));
      hit3->score_eventrim += hit3->amb_nmatches_end / ambig_end_interval;
      debug6(printf("  add amb end %d/%d.",hit3->amb_nmatches_end,ambig_end_interval));
      debug6(printf("  RESULT: %d\n",hit3->score_eventrim));
    } else {
      debug6(printf("score 3' OTHER:"));
      hit3->score_eventrim = hit3->penalties;
      debug6(printf("  penalties %d.",hit3->penalties));

      hit3->score_eventrim += Substring_count_mismatches_region(hit3->substring0,trim_left_3,trim_right_3,
								query3_compress_fwd,query3_compress_rev);
      debug6(printf("  substring 0 %d.",Substring_count_mismatches_region(hit3->substring0,trim_left_3,trim_right_3,
									  query3_compress_fwd,query3_compress_rev)));

      hit3->score_eventrim += Substring_count_mismatches_region(hit3->substring1,trim_left_3,trim_right_3,
							       query3_compress_fwd,query3_compress_rev);
      debug6(printf("  substring 1 %d.",Substring_count_mismatches_region(hit3->substring1,trim_left_3,trim_right_3,
									  query3_compress_fwd,query3_compress_rev)));

      hit3->score_eventrim += Substring_count_mismatches_region(hit3->substring2,trim_left_3,trim_right_3,
							       query3_compress_fwd,query3_compress_rev);
      debug6(printf("  substring 2 %d.",Substring_count_mismatches_region(hit3->substring2,trim_left_3,trim_right_3,
									  query3_compress_fwd,query3_compress_rev)));

      if (hit3->hittype == INSERTION || hit3->hittype == DELETION) {
	hit3->score_eventrim += indel_penalty_middle;
	debug6(printf("  add indel %d.",indel_penalty_middle));
      }
      debug6(printf("  RESULT: %d\n",hit3->score_eventrim));
    }

    hitpair->score_eventrim = hit5->score_eventrim + hit3->score_eventrim;
    if (hitpair->score_eventrim < minscore) {
      minscore = hitpair->score_eventrim;
    }
  }


  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    debug6(printf("%lu..%lu  %u..%u|%u..%u types %s and %s, score_eventrim %d+%d, pairlength %d, outerlength %u\n",
		  hitpair->low,hitpair->high,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),
		  hitpair->hit5->score_eventrim,hitpair->hit3->score_eventrim,
		  hitpair->insertlength,hitpair->outerlength));

    if (hitpair->hit5->hittype == TERMINAL && non_terminal_5p == true) {
      /* Don't use to determine minscore5 */
    } else if (hitpair->hit5->score_eventrim < minscore5) {
      minscore5 = hitpair->hit5->score_eventrim;
    }
    if (hitpair->hit3->hittype == TERMINAL && non_terminal_3p == true) {
      /* Don't use to determine minscore3 */
    } else if (hitpair->hit3->score_eventrim < minscore3) {
      minscore3 = hitpair->hit3->score_eventrim;
    }
  }
  debug6(printf("Stage3pair_optimal_score over %d pairs: minscore = %d and %d + subopt:%d\n",
		n,minscore5,minscore3,suboptimal_mismatches));

  if (non_double_terminal_p == true && finalp == false) {
    /* finalp == false.  Add suboptimal_mismatches to each end. */
    minscore5 += suboptimal_mismatches;
    minscore3 += suboptimal_mismatches;
    cutoff_level_5 = minscore5;
    cutoff_level_3 = minscore3;

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;

      if (hitpair->hit5->hittype == TERMINAL || hitpair->hit3->hittype == TERMINAL) {
	debug6(printf("Prefinal: Keeping a hit pair of type %s-%s with score_eventrim %d and %d, because finalp is false\n",
		      hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),
		      hitpair->hit5->score_eventrim,hitpair->hit3->score_eventrim));
	optimal = List_push(optimal,hitpair);

      } else if (keep_gmap_p == true && (hitpair->hit5->hittype == GMAP || hitpair->hit3->hittype == GMAP)) {
	/* GMAP hits already found to be better than their corresponding terminals */
	debug6(printf("Prefinal: Keeping a hit pair of type %s-%s with score_eventrim %d and %d, because keep_gmap_p is true\n",
		      hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),
		      hitpair->hit5->score_eventrim,hitpair->hit3->score_eventrim));
	optimal = List_push(optimal,hitpair);

      } else if (hitpair->hit5->score_eventrim > cutoff_level_5 && hitpair->hit3->score_eventrim > cutoff_level_3) {
	debug6(printf("Prefinal: Eliminating a hit pair at %lu..%lu  %u..%u|%u..%u with score_eventrim_5 %d > cutoff_level_5 %d and score_eventrim_3 %d > cutoff_level_3 %d (finalp %d)\n",
		      hitpair->low,hitpair->high,
		      hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		      hitpair->hit5->score_eventrim,cutoff_level_5,hitpair->hit3->score_eventrim,cutoff_level_3,finalp));
	*eliminatedp = true;
	Stage3pair_free(&hitpair);

      } else {
	debug6(printf("Prefinal: Keeping a hit pair with score_eventrim %d and %d (cutoff_level %d and %d)\n",
		      hitpair->hit5->score_eventrim,hitpair->hit3->score_eventrim,cutoff_level_5,cutoff_level_3));
	optimal = List_push(optimal,hitpair);
      }
    }

  } else {
    /* non_double_terminal_p == false (so need to prune results) or
       finalp == true.  Add suboptimal_mismatches to overall score. */
#if 0
    if (minscore5 + minscore3 < minscore) {
      cutoff_level = minscore + suboptimal_mismatches;
      debug6(printf("cutoff level %d = minscore %d + subopt %d\n",cutoff_level,minscore,suboptimal_mismatches));
    } else {
      cutoff_level = minscore5 + minscore3 + suboptimal_mismatches;
      debug6(printf("cutoff level %d = minscore5 %d + minscore3 %d + subopt %d\n",cutoff_level,minscore5,minscore3,suboptimal_mismatches));
    }
#else
    cutoff_level = minscore + suboptimal_mismatches;
    debug6(printf("cutoff level %d = minscore %d + subopt %d\n",cutoff_level,minscore,suboptimal_mismatches));
#endif


    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      debug6(printf("Final: %lu..%lu  %u..%u|%u..%u types %s and %s, score_eventrim %d (%d+%d), pairlength %d, outerlength %u\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),
		    hitpair->score_eventrim,hitpair->hit5->score_eventrim,hitpair->hit3->score_eventrim,
		    hitpair->insertlength,hitpair->outerlength));

#if 0
      if (hitpair->hit5->hittype != TERMINAL) {
	score5 = hitpair->hit5->score_eventrim;
      } else if (non_terminal_5p == true) {
	score5 = hitpair->hit5->score_eventrim;
      } else if (non_double_terminal_p == false) {
	score3 = hitpair->hit5->score_eventrim;
      } else {
	score5 = 0;
      }
      if (hitpair->hit3->hittype != TERMINAL) {
	score3 = hitpair->hit3->score_eventrim;
      } else if (non_terminal_3p == true) {
	score3 = hitpair->hit3->score_eventrim;
      } else if (non_double_terminal_p == false) {
	score3 = hitpair->hit3->score_eventrim;
      } else {
	score3 = 0;
      }
      score = score5 + score3;
#else
      score = hitpair->score_eventrim;
#endif

      if (keep_gmap_p == true && (hitpair->hit5->hittype == GMAP || hitpair->hit3->hittype == GMAP)) {
	/* GMAP hits already found to be better than their corresponding terminals */
	debug6(printf("Final: Keeping a hit pair of type %s-%s with score_eventrim %d, because keep_gmap_p is true\n",
		      hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),hitpair->score_eventrim));
	optimal = List_push(optimal,hitpair);

      } else if (score > cutoff_level) {
	debug6(printf("Final: Eliminating a hit pair at %lu..%lu  %u..%u|%u..%u with score %d > cutoff_level %d (finalp %d)\n",
		      hitpair->low,hitpair->high,
		      hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		      score,cutoff_level,finalp));
	*eliminatedp = true;
	Stage3pair_free(&hitpair);

      } else {
	debug6(printf("Final: Keeping a hit pair with score_eventrim %d (cutoff_level %d)\n",
		      hitpair->score_eventrim,cutoff_level));
	optimal = List_push(optimal,hitpair);
      }
    }
  }


  List_free(&hitpairlist);


#if 0
  /* Filter on pairlength */
  if (optimal != NULL) {
    hitpairlist = optimal;
    optimal = (List_T) NULL;

    hitpair = (Stage3pair_T) hitpairlist->first;
    best_absdifflength = hitpair->absdifflength;
    best_outerlength = hitpair->outerlength;

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->absdifflength < best_absdifflength) {
	best_absdifflength = hitpair->absdifflength;
	best_outerlength = hitpair->outerlength;
      } else if (hitpair->absdifflength > best_absdifflength) {
	/* Skip */
      } else if (hitpair->outerlength < best_outerlength) {
	best_outerlength = hitpair->outerlength;
      }
    }

    for (p = hitpairlist; p != NULL; p = p->rest) {
      hitpair = (Stage3pair_T) p->first;
      if (hitpair->absdifflength > best_absdifflength) {
	debug6(printf("Eliminating a hit pair with absdifflength %d\n",hitpair->absdifflength));
	*eliminatedp = true;
	Stage3pair_free(&hitpair);
      } else if (hitpair->outerlength > best_outerlength) {
	debug6(printf("Eliminating a hit pair with outerlength %u\n",hitpair->outerlength));
	*eliminatedp = true;
	Stage3pair_free(&hitpair);
      } else {
	debug6(printf("Keeping a hit pair with absdifflength %d and outerlength %d\n",
		     hitpair->absdifflength,hitpair->outerlength));
	optimal = List_push(optimal,hitpair);
      }
    }

    List_free(&hitpairlist);
  }
#endif

  return optimal;
}


List_T
Stage3pair_optimal_score (List_T hitpairlist, int cutoff_level, int suboptimal_mismatches,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			  bool keep_gmap_p, bool finalp) {
  List_T optimal;
  bool eliminatedp;

  optimal = Stage3pair_optimal_score_aux(&eliminatedp,hitpairlist,cutoff_level,suboptimal_mismatches,
					 query5_compress_fwd,query5_compress_rev,
					 query3_compress_fwd,query3_compress_rev,
					 keep_gmap_p,finalp);
  while (eliminatedp == true) {
    optimal = Stage3pair_optimal_score_aux(&eliminatedp,optimal,cutoff_level,suboptimal_mismatches,
					   query5_compress_fwd,query5_compress_rev,
					   query3_compress_fwd,query3_compress_rev,
					   keep_gmap_p,finalp);
  }

  return optimal;
}



bool
Stage3pair_sense_consistent_p (List_T hitpairlist) {
  Stage3pair_T hitpair;
  T hit5, hit3;
  List_T p;

  for (p = hitpairlist; p != NULL; p = List_next(p)) {
    hitpair = (Stage3pair_T) List_head(p);
    hit5 = hitpair->hit5;
    hit3 = hitpair->hit3;
    if (hit5->sensedir == hit3->sensedir) {
      return true;
    } else if (hit5->hittype == GMAP && hit5->gmap_nintrons > 0 && hit5->sensedir == 0) {
      /* false */
    } else if (hit3->hittype == GMAP && hit3->gmap_nintrons > 0 && hit3->sensedir == 0) {
      /* false */
    } else {
      return true;
    }
  }
  return false;
}


/* Want to unalias plus and alias minus */
List_T
Stage3end_linearize_5 (List_T hitlist) {
  Chrpos_T chrlength;
  T hit;
  List_T p;

  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (T) List_head(p);
    debug12(chrlength = hit->chrlength);
    debug12(printf("Looking at 5' end %lu..%lu against chrlength %u\n",
		   hit->genomicstart,hit->genomicend,chrlength));

    if (hit->alias == 0) {
      /* Skip */

    } else if (hit->alias > 0) {
      if (hit->plusp == true) {
	unalias_circular(hit);
      }

    } else {
      /* hit->alias < 0 */
      if (hit->plusp == false) {
	alias_circular(hit);
      }
    }
  }

  return hitlist;
}


/* Want to alias plus and unalias minus */
List_T
Stage3end_linearize_3 (List_T hitlist) {
  Chrpos_T chrlength;
  T hit;
  List_T p;

  for (p = hitlist; p != NULL; p = List_next(p)) {
    hit = (T) List_head(p);
    debug12(chrlength = hit->chrlength);
    debug12(printf("Looking at 3' end %lu..%lu against chrlength %u\n",
		   hit->genomicstart,hit->genomicend,chrlength));

    if (hit->alias == 0) {
      /* Skip */

    } else if (hit->alias < 0) {
      if (hit->plusp == true) {
	alias_circular(hit);
      }

    } else {
      /* hit->alias > 0 */
      if (hit->plusp == false) {
	unalias_circular(hit);
      }
    }
  }

  return hitlist;
}



List_T
Stage3pair_remove_circular_alias (List_T hitpairlist) {
  List_T newlist = NULL, p;
  Stage3pair_T hitpair;

  debug12(printf("Stage3pair_remove_circular_alias called with %d hitpairs\n",
		 List_length(hitpairlist)));
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->low >= hitpair->hit5->chroffset + hitpair->hit5->chrlength) {
      /* Both ends in circular alias */
      debug12(printf("Both ends in circular alias\n"));
      Stage3pair_free(&hitpair);

    } else {
      newlist = List_push(newlist,(void *) hitpair);
    }
  }

  List_free(&hitpairlist);
  return newlist;
}


static List_T
pair_up_concordant_aux (bool *abort_pairing_p, int *found_score, int *nconcordant, int *nsamechr,
			List_T *samechr, List_T *conc_transloc, List_T hitpairs,
			T **hits5_plus, int *nhits5_plus, T **hits5_minus, int *nhits5_minus,
			T **hits3_plus, int *nhits3_plus, T **hits3_minus, int *nhits3_minus,
			bool *sorted5p, bool *sorted3p,
			int cutoff_level_5, int cutoff_level_3, int subopt_levels,

			Univcoord_T *splicesites,
			Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			int querylength5, int querylength3, int maxpairedpaths,
			int splicing_penalty, int genestrand) {
  int new_found_score = *found_score;
  int pairscore, score5, score3, i, j;
  List_T q, prev_start;
  Stage3pair_T stage3pair;
  T *hits5, *hits3, hit5, hit3;
  int nhits5, nhits3;
  Univcoord_T insert_start;


  prev_start = hitpairs;
  pairscore = 0;
  while (*abort_pairing_p == false && pairscore <= *found_score + subopt_levels &&
	 pairscore <= cutoff_level_5 + cutoff_level_3) {
    debug5a(printf("pairscore = %d\n",pairscore));
    for (score5 = 0; score5 <= pairscore; score5++) {
      debug5a(printf("score5 = %d (cutoff %d), score3 = %d (cutoff %d)\n",
		     score5,cutoff_level_5,pairscore-score5,cutoff_level_3));

      if (score5 <= cutoff_level_5 && ((score3 = pairscore - score5) <= cutoff_level_3)) {
	/* Sort this level if necessary: 5' by genomicend, 3' by genomicstart */
	if (sorted5p[score5] == false) {
	  if (nhits5_plus[score5] > 0) {
	    qsort(hits5_plus[score5],nhits5_plus[score5],sizeof(T),genomicend_cmp);
	  }
	  if (nhits5_minus[score5] > 0) {
	    qsort(hits5_minus[score5],nhits5_minus[score5],sizeof(T),genomicend_cmp);
	  }
	  sorted5p[score5] = true;
	}
	if (sorted3p[score3] == false) {
	  if (nhits3_plus[score3] > 0) {
	    qsort(hits3_plus[score3],nhits3_plus[score3],sizeof(T),genomicstart_cmp);
	  }
	  if (nhits3_minus[score3] > 0) {
	    qsort(hits3_minus[score3],nhits3_minus[score3],sizeof(T),genomicstart_cmp);
	  }
	  sorted3p[score3] = true;
	}

	/* plus/plus: hits5_plus against hits3_plus (really on minus) */
	hits5 = hits5_plus[score5];
	hits3 = hits3_plus[score3];
	nhits5 = nhits5_plus[score5];
	nhits3 = nhits3_plus[score3];
	if (nhits5 > 0 && nhits3 > 0) {
	  debug5(printf("at score %d, nhits5_plus = %d; at score %d, nhits3_plus = %d\n",
			 score5,nhits5,score3,nhits3));

	  i = j = 0;
	  while (*abort_pairing_p == false && i < nhits5) {
	    hit5 = hits5[i];
	    insert_start = hit5->genomicend - querylength5;
	    debug5(printf("plus/plus: i=%d/%d %lu..%lu %s %s %p\n",
			  i,nhits5,hit5->genomicstart,hit5->genomicend,
			  print_sense(hit5->sensedir),hittype_string(hit5->hittype),hit5));

	    while (j >= 0 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %lu..%lu %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j--;
	    }
	    j++;		/* Finish backup */

	    while (j < nhits3 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %lu..%lu %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j++;
	    }

	    while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %lu..%lu %s %s %p",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      hit3 = hits3[j];
		
	      /* Want only pairs not previously seen */
	      if (hit5->paired_seenp == false || hit3->paired_seenp == false) {
		if (hit5->effective_chrnum != hit3->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 /* && hit5->other_chrnum != hit3->other_chrnum */) {
		  /* Could potentially miss an alignment if the two ends overlap */
		  debug5(printf(" => double splice translocations"));
		} else if (hit5->chrnum == 0 || hit3->chrnum == 0) {
		  debug5(printf(" => conc_transloc effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(hit5,hit3,splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/CONCORDANT_TRANSLOCATIONS,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {
		    *conc_transloc = List_push(*conc_transloc,(void *) stage3pair);
		  }

		} else if (SENSE_INCONSISTENT_P(hit5->sensedir,hit3->sensedir)) {
		  debug5(printf(" => sense inconsistent: %d | %d = %d",hit5->sensedir,hit3->sensedir,hit5->sensedir|hit3->sensedir));
		} else if (hit3->genomicend < hit5->genomicstart) {
		  debug5(printf(" => scramble because end3 %lu < start5 %lu\n",hit3->genomicend,hit5->genomicstart));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
		} else {
		  debug5(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(hit5,hit3,splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/CONCORDANT,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {

		    if (pairscore < new_found_score) {
		      new_found_score = pairscore;
		      debug5(printf(" => tentatively updating found_score to be %d",new_found_score));
		      hitpairs = List_push(hitpairs,(void *) stage3pair);
		      (*nconcordant)++;

		    } else {
		      hitpairs = List_push(hitpairs,(void *) stage3pair);
		      (*nconcordant)++;
		    }

		    if (0 && *nconcordant > maxpairedpaths) {
		      debug(printf(" -- %d concordant paths exceeds %d",*nconcordant,maxpairedpaths));
		      *abort_pairing_p = true;
		    }
		  }
		}
	      }
	      debug5(printf("\n"));

	      j++;
	    }
	    j--;		/* Finish advance */

	    i++;
	  }
	}

	/* minus/minus: hits3_minus (really on plus) against hits5_minus */
	hits3 = hits3_minus[score3];
	hits5 = hits5_minus[score5];
	nhits3 = nhits3_minus[score3];
	nhits5 = nhits5_minus[score5];
	if (nhits3 > 0 && nhits5 > 0) {
	  debug5(printf("at score %d, nhits5_minus = %d; at score %d, nhits3_minus = %d\n",
			score5,nhits5,score3,nhits3));

	  i = j = 0;
	  while (*abort_pairing_p == false && i < nhits3) {
	    hit3 = hits3[i];
	    insert_start = hit3->genomicstart - querylength3;
	    debug5(printf("minus/minus: i=%d/%d %lu..%lu %s %s %p\n",
			  i,nhits3,hit3->genomicstart,hit3->genomicend,
			  print_sense(hit3->sensedir),hittype_string(hit3->hittype),hit3));

	    while (j >= 0 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %lu..%lu %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j--;
	    }
	    j++;			/* Finish backup */

	    while (j < nhits5 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %lu..%lu %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j++;
	    }

	    while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %lu..%lu %s %s %p",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      hit5 = hits5[j];

	      /* Want only pairs not previously seen */
	      if (hit3->paired_seenp == false || hit5->paired_seenp == false) {
		if (hit3->effective_chrnum != hit5->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 /* && hit5->other_chrnum != hit3->other_chrnum */) {
		  /* Could potentially miss an alignment if the two ends overlap */
		  debug5(printf(" => double splice translocations"));

		} else if (hit5->chrnum == 0 || hit3->chrnum == 0) {
		  debug5(printf(" => conc_transloc effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(hit5,hit3,splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/CONCORDANT_TRANSLOCATIONS,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {
		    *conc_transloc = List_push(*conc_transloc,(void *) stage3pair);
		  }

		} else if (SENSE_INCONSISTENT_P(hit3->sensedir,hit5->sensedir)) {
		  debug5(printf(" => sense inconsistent: %d | %d = %d",hit5->sensedir,hit3->sensedir,hit5->sensedir|hit3->sensedir));
		} else if (hit5->genomicstart < hit3->genomicend) {
		  debug5(printf(" => scramble because start5 %lu < end3 %lu\n",hit5->genomicstart,hit3->genomicend));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
		} else {
		  debug5(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(hit5,hit3,splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/CONCORDANT,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {

		    if (pairscore < new_found_score) {
		      new_found_score = pairscore;
		      debug5(printf(" => updating new_found_score to be %d",new_found_score));
		      hitpairs = List_push(hitpairs,(void *) stage3pair);
		      (*nconcordant)++;

		    } else {
		      hitpairs = List_push(hitpairs,(void *) stage3pair);
		      (*nconcordant)++;
		    }

		    if (0 && *nconcordant > maxpairedpaths) {
		      debug(printf(" -- %d concordant paths exceeds %d",*nconcordant,maxpairedpaths));
		      *abort_pairing_p = true;
		    }
		  }
		}
	      }
	      debug5(printf("\n"));

	      j++;
	    }
	    j--;		/* Finish advance */

	    i++;
	  }
	}

	/* plus/minus (inversions): hits5_plus against hits3_minus */
	hits5 = hits5_plus[score5];
	hits3 = hits3_minus[score3];
	nhits5 = nhits5_plus[score5];
	nhits3 = nhits3_minus[score3];
	if (nhits5 > 0 && nhits3 > 0) {
	  debug5(printf("at score %d, nhits5_plus = %d; at score %d, nhits3_minus = %d\n",
			score5,nhits5,score3,nhits3));

	  i = j = 0;
	  while (*abort_pairing_p == false && i < nhits5) {
	    hit5 = hits5[i];
	    insert_start = hit5->genomicend - querylength5;
	    debug5(printf("plus/minus: i=%d/%d %lu..%lu %s %s %p\n",
			  i,nhits5,hit5->genomicstart,hit5->genomicend,
			  print_sense(hit5->sensedir),hittype_string(hit5->hittype),hit5));

	    while (j >= 0 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %lu..%lu %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j--;
	    }
	    j++;		/* Finish backup */

	    while (j < nhits3 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %lu..%lu %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j++;
	    }

	    while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %lu..%lu %s %s %p",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      hit3 = hits3[j];
		
	      /* Want only pairs not previously seen */
	      if (hit5->paired_seenp == false || hit3->paired_seenp == false) {
		if (hit5->effective_chrnum != hit3->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 /* && hit5->other_chrnum != hit3->other_chrnum */) {
		  debug5(printf(" => double splice translocations"));
		} else if (SENSE_INCONSISTENT_FOR_INVERSION_P(hit5->sensedir,hit3->sensedir)) {
		  debug5(printf(" => sense inconsistent for inversion"));
#if 0
		} else if (hits3[j]->genomicstart + querylength3 <= insert_start) {
		  debug5(printf(" => scramble"));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
#endif
		} else {
		  debug5(printf(" => inversion effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_INVERSION,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
		}
	      }
	      debug5(printf("\n"));

	      j++;
	    }
	    j--;		/* Finish advance */

	    i++;
	  }
	}

	/* minus/plus: hits3_plus against hits5_minus */
	hits3 = hits3_plus[score3];
	hits5 = hits5_minus[score5];
	nhits3 = nhits3_plus[score3];
	nhits5 = nhits5_minus[score5];
	if (nhits3 > 0 && nhits5 > 0) {
	  debug5(printf("at score %d, nhits5_minus = %d; at score %d, nhits3_plus = %d\n",
			score5,nhits5,score3,nhits3));

	  i = j = 0;
	  while (*abort_pairing_p == false && i < nhits3) {
	    hit3 = hits3[i];
	    insert_start = hit3->genomicstart - querylength3;
	    debug5(printf("minus/plus: i=%d/%d %lu..%lu %s %s %p\n",
			  i,nhits3,hit3->genomicstart,hit3->genomicend,
			  print_sense(hit3->sensedir),hittype_string(hit3->hittype),hit3));

	    while (j >= 0 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %lu..%lu %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j--;
	    }
	    j++;			/* Finish backup */

	    while (j < nhits5 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %lu..%lu %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j++;
	    }

	    while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %lu..%lu %s %s %p",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      hit5 = hits5[j];

	      /* Want only pairs not previously seen */
	      if (hit3->paired_seenp == false || hit5->paired_seenp == false) {
		if (hit3->effective_chrnum != hit5->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 /* && hit5->other_chrnum != hit3->other_chrnum */) {
		  debug5(printf(" => double splice translocations"));
		} else if (SENSE_INCONSISTENT_FOR_INVERSION_P(hit3->sensedir,hit5->sensedir)) {
		  debug5(printf(" => sense inconsistent for inversion"));
#if 0
		} else if (hits5[j]->genomicend + querylength5 <= insert_start) {
		  debug5(printf(" => scramble"));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
#endif
		} else {
		  debug5(printf(" => inversion effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if (*nsamechr <= maxpairedpaths &&
		      (stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_INVERSION,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		    (*nsamechr)++;
		  }
		}
	      }
	      debug5(printf("\n"));

	      j++;
	    }
	    j--;		/* Finish advance */

	    i++;
	  }
	}

      }
    }

    /* Mark all concordant pairs (and terminal and double_terminal) found at this pairscore level */
    for (q = hitpairs; q != prev_start; q = List_next(q)) {
      stage3pair = (Stage3pair_T) List_head(q);
      stage3pair->hit5->paired_seenp = true;
      stage3pair->hit3->paired_seenp = true;
    }
    prev_start = hitpairs;

    if (*abort_pairing_p == false) {
      *found_score = new_found_score;
    }

    pairscore++;
  }

  return hitpairs;
}


static int
sort_hits_by_score (T ***hits_plus, T ***hits_minus, int **nhits_plus, int **nhits_minus,
		    List_T *hitarray, int narray, int cutoff_level) {
  int score;
  int nhits, i;
  List_T q;
  T hit;

  /* Find sizes for allocating memory */
  debug5(printf("Sizes of pieces by score level:\n"));
  *nhits_plus = (int *) CALLOC(cutoff_level+1,sizeof(int));
  *nhits_minus = (int *) CALLOC(cutoff_level+1,sizeof(int));
  for (i = 0; i < narray; i++) {
    debug5(printf("  array score level %d with %d hits\n",i,List_length(hitarray[i])));
    for (q = hitarray[i]; q != NULL; q = q->rest) {
      hit = (T) q->first;
      debug5(printf(" : %p score %d, type %s\n",hit,hit->score,hittype_string(hit->hittype)));
      assert(hit->score >= 0);
      if (hit->score > cutoff_level) {
	debug5(printf("Skipping hit with score %d > cutoff level %d\n",hit->score,cutoff_level));
      } else if (hit->plusp == true) {
	(*nhits_plus)[hit->score]++;
      } else {
	(*nhits_minus)[hit->score]++;
      }
    }
  }

  debug5(
	 printf("Sizes of pieces by score level and plus/minus:\n");
	 for (score = 0; score <= cutoff_level; score++) {
	   printf("  score %d: %d plus, %d minus\n",score,(*nhits_plus)[score],(*nhits_minus)[score]);
	 }
	 );


  /* Reset cutoff_level */
  score = 0;
  nhits = (*nhits_plus)[score] + (*nhits_minus)[score];
  while (score+1 <= cutoff_level && nhits + (*nhits_plus)[score+1] + (*nhits_minus)[score+1] < MAX_HITS) {
    nhits += (*nhits_plus)[score+1] + (*nhits_minus)[score+1];
    score++;
    debug5(printf("Allowing score to go to %d, because nhits = %d\n",score,nhits));
  }
  debug5(printf("Resetting cutoff_level to be %d\n",score));
  cutoff_level = score;


  /* Store hits */
  *hits_plus = (T **) CALLOC(cutoff_level+1,sizeof(T *));
  for (score = 0; score <= cutoff_level; score++) {
    if ((*nhits_plus)[score] == 0) {
      (*hits_plus)[score] = (T *) NULL;
    } else {
      (*hits_plus)[score] = (T *) CALLOC((*nhits_plus)[score],sizeof(Stage3end_T));
    }
  }

  *hits_minus = (T **) CALLOC(cutoff_level+1,sizeof(T *));
  for (score = 0; score <= cutoff_level; score++) {
    if ((*nhits_minus)[score] == 0) {
      (*hits_minus)[score] = (T *) NULL;
    } else {
      (*hits_minus)[score] = (T *) CALLOC((*nhits_minus)[score],sizeof(Stage3end_T));
    }
  }

  for (score = 0; score <= cutoff_level; score++) {
    (*nhits_plus)[score] = 0;
    (*nhits_minus)[score] = 0;
  }

  for (i = 0; i < narray; i++) {
    for (q = hitarray[i]; q != NULL; q = q->rest) {
      hit = (T) q->first;
      if (hit->score > cutoff_level) {
	/* Skip */
      } else if (hit->plusp == true) {
	(*hits_plus)[hit->score][(*nhits_plus)[hit->score]++] = hit;
      } else {
	(*hits_minus)[hit->score][(*nhits_minus)[hit->score]++] = hit;
      }
    }
  }

  return cutoff_level;
}


/* Finds concordant pairs if nconcordant is 0 */
List_T
Stage3_pair_up_concordant (bool *abort_pairing_p, int *found_score, int *nconcordant, int *nsamechr,
			   List_T *samechr, List_T *conc_transloc, List_T *with_terminal,
			   List_T hitpairs, List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			   List_T terminals5, List_T terminals3,
			   int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			   Univcoord_T *splicesites,
			   Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			   Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			   int querylength5, int querylength3, int maxpairedpaths,
			   int splicing_penalty, int genestrand) {
  T **hits5_plus, **hits5_minus, **hits3_plus, **hits3_minus;
  int *nhits5_plus, *nhits5_minus, *nhits3_plus, *nhits3_minus;
  T **terminals5_plus, **terminals5_minus, **terminals3_plus, **terminals3_minus;
  int *nterminals5_plus, *nterminals5_minus, *nterminals3_plus, *nterminals3_minus;
  int score5, score3;
  bool *sorted_hits5_p, *sorted_hits3_p, *sorted_terminals5_p = NULL, *sorted_terminals3_p = NULL;
  int cutoff_level_hits5, cutoff_level_hits3, cutoff_level_terminals5 = 0, cutoff_level_terminals3 = 0;
  int ignore_found_score;
  
  debug5(printf("Starting Stage3_pair_up_concordant with %d concordant, narray5 %d, narray3 %d, terminals5 %p, terminals3 %p, found_score %d\n",
		*nconcordant,narray5,narray3,terminals5,terminals3,*found_score));

  /* Find sizes for allocating memory */
  debug5(printf("Sizes of 5-end pieces by score level:\n"));
  cutoff_level_hits5 = sort_hits_by_score(&hits5_plus,&hits5_minus,&nhits5_plus,&nhits5_minus,
					  hitarray5,narray5,cutoff_level_5);
  sorted_hits5_p = (bool *) CALLOC(cutoff_level_hits5+1,sizeof(bool));

  debug5(printf("Sizes of 3-end pieces by score level:\n"));
  cutoff_level_hits3 = sort_hits_by_score(&hits3_plus,&hits3_minus,&nhits3_plus,&nhits3_minus,
					  hitarray3,narray3,cutoff_level_3);
  sorted_hits3_p = (bool *) CALLOC(cutoff_level_hits3+1,sizeof(bool));

  if (terminals5 == NULL && terminals3 == NULL) {
    /* Look for concordant pairs among the non-terminals */
    hitpairs = pair_up_concordant_aux(&(*abort_pairing_p),&(*found_score),&(*nconcordant),&(*nsamechr),
				      &(*samechr),&(*conc_transloc),hitpairs,
				      hits5_plus,nhits5_plus,hits5_minus,nhits5_minus,
				      hits3_plus,nhits3_plus,hits3_minus,nhits3_minus,
				      /*sorted5p*/sorted_hits5_p,/*sorted3p*/sorted_hits3_p,
				      cutoff_level_hits5,cutoff_level_hits3,subopt_levels,
				      splicesites,query5_compress_fwd,query5_compress_rev,
				      query3_compress_fwd,query3_compress_rev,querylength5,querylength3,
				      maxpairedpaths,splicing_penalty,genestrand);

  } else {
    /* Look for single terminals */
    if (terminals3 != NULL) {
      cutoff_level_terminals3 = sort_hits_by_score(&terminals3_plus,&terminals3_minus,&nterminals3_plus,&nterminals3_minus,
						   &terminals3,/*narray3*/1,cutoff_level_3);
      sorted_terminals3_p = (bool *) CALLOC(cutoff_level_terminals3+1,sizeof(bool));

      /* Do not allow terminals to alter found_score */
      ignore_found_score = *found_score;
      *with_terminal = pair_up_concordant_aux(&(*abort_pairing_p),&ignore_found_score,&(*nconcordant),&(*nsamechr),
					      &(*samechr),&(*conc_transloc),*with_terminal,
					      hits5_plus,nhits5_plus,hits5_minus,nhits5_minus,
					      terminals3_plus,nterminals3_plus,terminals3_minus,nterminals3_minus,
					      /*sorted5p*/sorted_hits5_p,/*sorted3p*/sorted_terminals3_p,
					      cutoff_level_hits5,cutoff_level_terminals3,subopt_levels,
					      splicesites,query5_compress_fwd,query5_compress_rev,
					      query3_compress_fwd,query3_compress_rev,querylength5,querylength3,
					      maxpairedpaths,splicing_penalty,genestrand);
    }

    if (terminals5 != NULL) {
      cutoff_level_terminals5 = sort_hits_by_score(&terminals5_plus,&terminals5_minus,&nterminals5_plus,&nterminals5_minus,
						   &terminals5,/*narray5*/1,cutoff_level_5);
      sorted_terminals5_p = (bool *) CALLOC(cutoff_level_terminals5+1,sizeof(bool));

      /* Do not allow terminals to alter found_score */
      ignore_found_score = *found_score;
      *with_terminal = pair_up_concordant_aux(&(*abort_pairing_p),&ignore_found_score,&(*nconcordant),&(*nsamechr),
					      &(*samechr),&(*conc_transloc),*with_terminal,
					      terminals5_plus,nterminals5_plus,terminals5_minus,nterminals5_minus,
					      hits3_plus,nhits3_plus,hits3_minus,nhits3_minus,
					      /*sorted5p*/sorted_terminals5_p,/*sorted3p*/sorted_hits3_p,
					      cutoff_level_terminals5,cutoff_level_hits3,subopt_levels,
					      splicesites,query5_compress_fwd,query5_compress_rev,
					      query3_compress_fwd,query3_compress_rev,querylength5,querylength3,
					      maxpairedpaths,splicing_penalty,genestrand);
    }

    /* Previously required *with_terminal to be NULL also */
    if (terminals3 != NULL && terminals5 != NULL) {
      /* Look for double terminals */
      /* Do not allow terminals to alter found_score */
      ignore_found_score = *found_score;
      *with_terminal = pair_up_concordant_aux(&(*abort_pairing_p),&ignore_found_score,&(*nconcordant),&(*nsamechr),
					      &(*samechr),&(*conc_transloc),*with_terminal,
					      terminals5_plus,nterminals5_plus,terminals5_minus,nterminals5_minus,
					      terminals3_plus,nterminals3_plus,terminals3_minus,nterminals3_minus,
					      /*sorted5p*/sorted_terminals5_p,/*sorted3p*/sorted_terminals3_p,
					      cutoff_level_terminals5,cutoff_level_terminals3,subopt_levels,
					      splicesites,query5_compress_fwd,query5_compress_rev,
					      query3_compress_fwd,query3_compress_rev,querylength5,querylength3,
					      maxpairedpaths,splicing_penalty,genestrand);
    }

    if (sorted_terminals3_p != NULL) {
      FREE(sorted_terminals3_p);
      for (score3 = 0; score3 <= cutoff_level_terminals3; score3++) {
	FREE(terminals3_plus[score3]);
	FREE(terminals3_minus[score3]);
      }
      FREE(terminals3_plus);
      FREE(terminals3_minus);
      FREE(nterminals3_plus);
      FREE(nterminals3_minus);
    }

    if (sorted_terminals5_p != NULL) {
      FREE(sorted_terminals5_p);
      for (score5 = 0; score5 <= cutoff_level_terminals5; score5++) {
	FREE(terminals5_plus[score5]);
	FREE(terminals5_minus[score5]);
      }
      FREE(terminals5_plus);
      FREE(terminals5_minus);
      FREE(nterminals5_plus);
      FREE(nterminals5_minus);
    }
  }


  FREE(sorted_hits3_p);
  for (score3 = 0; score3 <= cutoff_level_hits3; score3++) {
    FREE(hits3_plus[score3]);
    FREE(hits3_minus[score3]);
  }
  FREE(hits3_plus);
  FREE(hits3_minus);
  FREE(nhits3_plus);
  FREE(nhits3_minus);

  FREE(sorted_hits5_p);
  for (score5 = 0; score5 <= cutoff_level_hits5; score5++) {
    FREE(hits5_plus[score5]);
    FREE(hits5_minus[score5]);
  }
  FREE(hits5_plus);
  FREE(hits5_minus);
  FREE(nhits5_plus);
  FREE(nhits5_minus);

  debug5(printf("Finished with Stage3_pair_up_concordant: %d concordant, %d samechr, %d conc_transloc, %d with_terminal\n",
		List_length(hitpairs),List_length(*samechr),List_length(*conc_transloc),List_length(*with_terminal)));

  return hitpairs;
}


