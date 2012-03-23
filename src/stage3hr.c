static char rcsid[] = "$Id: stage3hr.c 60000 2012-03-20 19:45:48Z twu $";
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


#define MAX_HITS 100000


#define CONCORDANT_TEXT "concordant"
#define PAIRED_TEXT "paired"
#define UNPAIRED_TEXT "unpaired"

#ifdef USE_TALLY_RATIO
#define TALLY_RATIO 2.0
#endif

#define SUBSUMPTION_SLOP 10	/* Should allow for short insert lengths */
/* #define TERMINAL_SECOND_CLASS 1 -- enabling this leads to poor results */


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

/* Stage3_optimal_score */
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



#define MAPQ_MAXIMUM_SCORE 40


static bool invert_first_p;
static bool invert_second_p;
static IIT_T genes_iit;
static int *genes_divint_crosstable;
static IIT_T tally_iit;
static int *tally_divint_crosstable;
static IIT_T runlength_iit;
static int *runlength_divint_crosstable;

static int pairmax;
static int expected_pairlength;
static int pairlength_deviation;

static int localsplicing_penalty;
static int indel_penalty_middle;
static int antistranded_penalty;
static bool favor_multiexon_p;


/* Probably not good to use in certain genomic regions, unless we also
   use known splicesites with distance information. */
/* But sometimes need to use to get correct mapping */
static bool favor_ambiguous_p;


void
Stage3hr_setup (bool invert_first_p_in, bool invert_second_p_in,
		IIT_T genes_iit_in, int *genes_divint_crosstable_in,
		IIT_T tally_iit_in, int *tally_divint_crosstable_in,
		IIT_T runlength_iit_in, int *runlength_divint_crosstable_in,
		bool distances_observed_p,
		int pairmax_in, int expected_pairlength_in, int pairlength_deviation_in,
		int localsplicing_penalty_in, int indel_penalty_middle_in,
		int antistranded_penalty_in, bool favor_multiexon_p_in) {
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

  pairmax = pairmax_in;
  expected_pairlength = expected_pairlength_in;
  pairlength_deviation = pairlength_deviation_in;

  if (distances_observed_p == true) {
    favor_ambiguous_p = false;
  } else {
    favor_ambiguous_p = true;
  }

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

  Chrnum_T chrnum; /* Needed for printing paired-end results.  A chrnum of 0 indicates a distant splice. */
  Chrnum_T effective_chrnum;	/* For determining concordance */
  Chrnum_T other_chrnum;	/* 0 for non-translocations, and other chrnum besides effective_chrnum for translocations */
  Genomicpos_T chroffset;
  Genomicpos_T chrhigh;

  int querylength_adj;		/* Adjusted for insertions */

  Genomicpos_T genomicstart;
  Genomicpos_T genomicend;
  bool plusp;

  Genomicpos_T low;
  Genomicpos_T high;

  double mapq_loglik;
  int mapq_score;
  int absmq_score;		/* Absolute MAPQ, for XQ and X2 flags */

  int score;			/* Includes colordiffs and penalties */
  int ntscore;			/* Includes penalties */
  int nmatches;
  int nmatches_posttrim;

  int score_eventrim; /* Used by Stage3end_optimal_score for comparing terminals and non-terminals */


  Overlap_T gene_overlap;
  long int tally;

  int nmismatches_whole;
  int nmismatches_bothdiff;
  int nmismatches_refdiff;	/* Set only for display */

  int nindels;			/* for indels */
  int indel_pos;		/* for indels.  Relative to querypos 0 */
  int indel_low;		/* for indels.  Relative to chromosomal low end of read, but still 0 if no indel. */
  char *deletion;		/* for deletions */

  Genomicpos_T distance;	/* for splicing or shortexon (sum of two distances) */
  Genomicpos_T acceptor_distance; /* for shortexon */
  Genomicpos_T donor_distance;	  /* for shortexon */
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
  bool gmap_triedp;
};


static char *
pairtype_string (Pairtype_T pairtype) {
  switch (pairtype) {
  case CONCORDANT: return "concordant";
  case PAIRED_UNSPECIFIED: return "paired_unknown";
  case PAIRED_INVERSION: return "inversion";
  case PAIRED_SCRAMBLE: return "scramble";
  case PAIRED_TOOLONG: return "toolong";
  case TRANSLOCATION: return "translocation";
  case UNPAIRED: return "unpaired";
  default: fprintf(stderr,"Unexpected pairtype %d\n",pairtype); abort();
  }
}


struct Stage3pair_T {
  Pairtype_T pairtype;
  int genestrand;

  T hit5;
  T hit3;
  bool private5p;			/* A private copy separate from hits5 and hits3, and not a pointer */
  bool private3p;

  Genomicpos_T low;
  Genomicpos_T high;
  int insertlength;
  Genomicpos_T outerlength;

  double mapq_loglik;
  int mapq_score;
  int absmq_score;

  int score;
  int nmatches;
  int nmatches_posttrim;
  int indel_low; /* For ranking identical indel alignments, so we pick lowest coord */

  Overlap_T gene_overlap;
  long int tally;

  Genomicpos_T absdifflength;
#ifdef USE_BINGO
  bool absdifflength_bingo_p;
#endif
  int dir;			/* -1, 0, or +1 */
  bool sense_consistent_p;

  int nchimera_known;
  int nchimera_novel;
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

Genomicpos_T
Stage3end_chroffset (T this) {
  return this->chroffset;
}

Genomicpos_T
Stage3end_chrhigh (T this) {
  return this->chrhigh;
}

Genomicpos_T
Stage3end_genomicstart (T this) {
  return this->genomicstart;
}

Genomicpos_T
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
Genomicpos_T
Stage3end_genomic_alignment_length (T this) {
  Genomicpos_T length;

  length = Substring_genomic_alignment_length(this->substring1);
  length += Substring_genomic_alignment_length(this->substring2);
  length += Substring_genomic_alignment_length(this->substring0);
  if (this->hittype == DELETION) {
    length += (Genomicpos_T) this->nindels;
  }
  return length;
}


Genomicpos_T
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


Genomicpos_T
Stage3end_distance (T this) {
  return this->distance;
}

Genomicpos_T
Stage3end_shortexon_acceptor_distance (T this) {
  return this->acceptor_distance;
}

Genomicpos_T
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
  if (this->sensedir == SENSE_FORWARD) {
    return +1;
  } else if (this->sensedir == SENSE_ANTI) {
    return -1;
  } else {
    return 0;
  }
}

bool
Stage3end_gmap_triedp (T this) {
  return this->gmap_triedp;
}

void
Stage3end_set_gmap_triedp (T this) {
  this->gmap_triedp = true;
  return;
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
  assert(this->hittype != GMAP);

  /* indel + splice => requires gmap
     doublesplice + splice => requires gmap
     other cases should have already been covered by gsnap
  */

  if (this->hittype != INSERTION && this->hittype != DELETION && this->hittype != SHORTEXON) {
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


bool
Stage3end_bad_stretch_p (T this, Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  if (this->hittype == GMAP) {
    return Stage3_bad_stretch_p(this->pairarray,this->npairs);
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
Stage3end_genomicbound_from_start (Genomicpos_T *genomicbound, T this, int overlap, Genomicpos_T chroffset) {
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
Stage3end_genomicbound_from_end (Genomicpos_T *genomicbound, T this, int overlap, Genomicpos_T chroffset) {
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

Genomicpos_T
Stage3pair_pairlength (Stage3pair_T this) {
  return this->insertlength;
}

int
Stage3pair_nmatches (Stage3pair_T this) {
  return this->nmatches;
}


int
Stage3pair_overlap (Stage3pair_T this) {
  Stage3end_T hit5, hit3;
  int totallength;

  hit5 = this->hit5;
  hit3 = this->hit3;

  if (hit5->hittype == SAMECHR_SPLICE || hit5->hittype == TRANSLOC_SPLICE) {
    return 0;
  } else if (hit3->hittype == SAMECHR_SPLICE || hit3->hittype == TRANSLOC_SPLICE) {
    return 0;
  } else if (this->insertlength <= hit5->querylength_adj) {
    return 0;
  } else if (this->insertlength <= hit3->querylength_adj) {
    return 0;
  } else {
    totallength = hit5->querylength_adj + hit3->querylength_adj;
    if (this->insertlength >= totallength) {
      return 0;
    } else {
      return (totallength - this->insertlength);
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

  new->chrnum = old->chrnum;
  new->effective_chrnum = old->effective_chrnum;
  new->other_chrnum = old->other_chrnum;
  new->chroffset = old->chroffset;
  new->chrhigh = old->chrhigh;

  new->querylength_adj = old->querylength_adj;

  new->genomicstart = old->genomicstart;
  new->genomicend = old->genomicend;
  new->plusp = old->plusp;

  new->low = old->low;
  new->high = old->high;

  new->mapq_loglik = old->mapq_loglik;
  new->mapq_score = old->mapq_score;
  new->absmq_score = old->absmq_score;

  new->score = old->score;
  new->ntscore = old->ntscore;
  new->nmatches = old->nmatches;
  new->nmatches_posttrim = old->nmatches_posttrim;

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
  new->gmap_triedp = old->gmap_triedp;

  return new;
}



T
Stage3end_new_exact (int *found_score, Genomicpos_T left, int genomiclength, Compress_T query_compress,
		     bool plusp, int genestrand, Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh) {
  T new;
  Substring_T substring;
  Genomicpos_T genomicstart, genomicend;

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + genomiclength;
  } else {
    genomicstart = left + genomiclength;
    genomicend = left;
  }

  if ((substring = Substring_new(/*nmismatches*/0,chrnum,chroffset,chrhigh,
				 left,genomicstart,genomicend,query_compress,
				 /*start_endtype*/END,/*end_endtype*/END,
				 /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				 /*alignstart*/genomicstart,/*alignend*/genomicend,
				 genomiclength,/*extraleft*/0,/*extraright*/0,/*exactp*/true,plusp,genestrand,
				 /*trim_left_p*/false,/*trim_right_p*/false,/*minlength*/0)) == NULL) {
    return (T) NULL;

  } else {
    new = (T) MALLOC_OUT(sizeof(*new));
    debug0(printf("Stage3end_new_exact %p: left %u, chrnum %d\n",new,left,chrnum));

    new->substring1 = substring;
    new->substring2 = (Substring_T) NULL;
    new->substring0 = (Substring_T) NULL;
    new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
    new->substringD = new->substringA = (Substring_T) NULL;
    new->substring_low = new->substring_high = new->substring1;

    new->pairarray = (struct Pair_T *) NULL;

    new->deletion = (char *) NULL;
    new->querylength_adj = genomiclength;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    if (genomicstart < genomicend) {
      new->low = genomicstart;
      new->high = genomicend;
    } else {
      new->low = genomicend;
      new->high = genomicstart;
    }
    debug0(printf("Assigned %u to low and %u to high\n",new->low,new->high));


    new->hittype = EXACT;
    new->genestrand = genestrand;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
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
    new->gmap_triedp = false;

    return new;
  }
}


T
Stage3end_new_substitution (int *found_score, int nmismatches_whole, Genomicpos_T left,
			    int genomiclength, Compress_T query_compress,
			    bool plusp, int genestrand, Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh) {
  T new;
  Substring_T substring;
  Genomicpos_T genomicstart, genomicend;

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + genomiclength;
  } else {
    genomicstart = left + genomiclength;
    genomicend = left;
  }

  if ((substring = Substring_new(nmismatches_whole,chrnum,chroffset,chrhigh,
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
    debug0(printf("Stage3end_new_substitution %p: left %u, chrnum %d, nmismatches %d\n",
		  new,left,chrnum,nmismatches_whole));

    new->substring1 = substring;
    new->substring2 = (Substring_T) NULL;
    new->substring0 = (Substring_T) NULL;
    new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
    new->substringD = new->substringA = (Substring_T) NULL;
    new->substring_low = new->substring_high = new->substring1;

    new->pairarray = (struct Pair_T *) NULL;

    new->deletion = (char *) NULL;
    new->querylength_adj = genomiclength;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    if (genomicstart < genomicend) {
      new->low = genomicstart;
      new->high = genomicend;
    } else {
      new->low = genomicend;
      new->high = genomicstart;
    }

    if (nmismatches_whole == 0) {
      /* Proper hittype needed so we can eliminate identical hits */
      new->hittype = EXACT;
    } else {
      new->hittype = SUB;
    }
    new->genestrand = genestrand;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
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
    new->gmap_triedp = false;

    return new;
  }
}



T
Stage3end_new_insertion (int *found_score, int nindels, int indel_pos, int nmismatches1_whole, int nmismatches2_whole,
			 Genomicpos_T left, int genomiclength, Compress_T query_compress,
			 int querylength, bool plusp, int genestrand, Chrnum_T chrnum, Genomicpos_T chroffset,
			 Genomicpos_T chrhigh, int indel_penalty) {
  T new;
  Substring_T substring1, substring2;
  int querystart1, queryend1, querystart2, queryend2;
  Genomicpos_T genomicstart, genomicend;
  Genomicpos_T alignstart1, alignend1, alignstart2, alignend2;

  debug2(printf("Entered with left %u, querylength %d, genomiclength %d, indel_pos %d\n",
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
    genomicstart = left;
    genomicend = left + genomiclength;

    alignstart1 = genomicstart;
    alignend1 = alignstart2 = genomicstart + indel_pos;
    alignend2 = genomicend/* - nindels*/;

  } else {
    genomicstart = left + genomiclength;
    genomicend = left;

    alignstart1 = genomicstart;
    alignend1 = alignstart2 = genomicstart - indel_pos;
    alignend2 = genomicend/* + nindels*/;

  }

  if ((substring1 = Substring_new(nmismatches1_whole,chrnum,chroffset,chrhigh,
				  left,genomicstart,genomicend,query_compress,
				  /*start_endtype*/END,/*end_endtype*/INS,
				  querystart1,queryend1,querylength,alignstart1,alignend1,genomiclength,
				  /*extraleft*/0,/*extraright*/0,/*exactp*/false,plusp,genestrand,
				  /*trim_left_p (previously was end1_indel_p ? false : true)*/true,
				  /*trim_right_p*/false,/*minlength*/0)) == NULL) {
    return (T) NULL;

  } else if ((substring2 = Substring_new(nmismatches2_whole,chrnum,chroffset,chrhigh,
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
    debug0(printf("Stage3end_new_insertion %p: left %u, chrnum %d, nmismatches %d+%d, indel_pos %d, nindels %d\n",
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
    new->querylength_adj = querylength /* - nindels */;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    if (genomicstart < genomicend) {
      new->low = genomicstart;
      new->high = genomicend;
    } else {
      new->low = genomicend;
      new->high = genomicstart;
    }


    new->hittype = INSERTION;
    new->genestrand = genestrand;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
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
#endif
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
    new->gmap_triedp = false;

    return new;
  }
}


T
Stage3end_new_deletion (int *found_score, int nindels, int indel_pos, int nmismatches1_whole, int nmismatches2_whole,
			Genomicpos_T left, int genomiclength, Compress_T query_compress,
			int querylength, bool plusp, int genestrand, Chrnum_T chrnum, Genomicpos_T chroffset,
			Genomicpos_T chrhigh, int indel_penalty) {
  T new;
  Substring_T substring1, substring2;
  int querystart1, queryend1, querystart2, queryend2;
  Genomicpos_T genomicstart, genomicend;
  Genomicpos_T alignstart1, alignend1, alignstart2, alignend2;
  Genomicpos_T left2;

  debug3(printf("Entered with left %u, querylength %d, genomiclength %d, indel_pos %d\n",
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
    genomicstart = left;
    genomicend = left + genomiclength;

    alignstart1 = genomicstart;
    alignend1 = genomicstart + indel_pos;
    alignstart2 = alignend1 + nindels;
    alignend2 = genomicend/* + nindels*/;

    /* left1 = left; */
    left2 = left + nindels;

    debug3(printf("plusp is true.  genomicstart %u, genomicend %u, alignstart1 %u, alignend1 %u, alignstart2 %u, alignend2 %u, left1 %u, left2 %u\n",
		  genomicstart,genomicend,alignstart1,alignend1,alignstart2,alignend2,left,left2));

  } else {
    genomicstart = left + genomiclength;
    genomicend = left;

    alignstart1 = genomicstart;
    alignend1 = genomicstart - indel_pos;
    alignstart2 = alignend1 - nindels;
    alignend2 = genomicend/* - nindels*/;

    /* left1 = left; */
    left2 = left + nindels;

    debug3(printf("plusp is false.  genomicstart %u, genomicend %u, alignstart1 %u, alignend1 %u, alignstart2 %u, alignend2 %u, left1 %u, left2 %u\n",
		  genomicstart,genomicend,alignstart1,alignend1,alignstart2,alignend2,left,left2));
  }


  if ((substring1 = Substring_new(nmismatches1_whole,chrnum,chroffset,chrhigh,
				  left,genomicstart,genomicend,query_compress,
				  /*start_endtype*/END,/*end_endtype*/DEL,
				  querystart1,queryend1,querylength,alignstart1,alignend1,genomiclength,
				  /*extraleft*/0,/*extraright*/0,/*exactp*/false,plusp,genestrand,
				  /*trim_left_p (previously was end1_indel_p ? false : true)*/true,
				  /*trim_right_p*/false,/*minlength*/0)) == NULL) {
    return (T) NULL;

  } else if ((substring2 = Substring_new(nmismatches2_whole,chrnum,chroffset,chrhigh,
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
    debug0(printf("Stage3end_new_deletion %p: left %u, chrnum %d, nmismatches %d+%d, indel_pos %d, nindels %d\n",
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

    new->hittype = DELETION;
    new->genestrand = genestrand;

    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->chrhigh = chrhigh;
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
#endif
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
    new->gmap_triedp = false;

    return new;
  }
}


/* Never returns NULL */
T
Stage3end_new_splice (int *found_score, int nmismatches_donor, int nmismatches_acceptor,
		      Substring_T donor, Substring_T acceptor, Genomicpos_T distance,
		      bool shortdistancep, int splicing_penalty, int querylength,
		      int amb_nmatches, Intlist_T ambi_left, Intlist_T ambi_right,
		      Intlist_T amb_nmismatches_left, Intlist_T amb_nmismatches_right,
		      bool copy_donor_p, bool copy_acceptor_p, bool first_read_p, int sensedir) {
  T new;
  int ignore;
  Substring_T substring_for_concordance; /* always the inner substring */
  
  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_splice %p with sensedir %d\n",new,sensedir));

  new->deletion = (char *) NULL;
  new->querylength_adj = querylength;

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

  if (donor == NULL) {
    new->hittype = HALFSPLICE_ACCEPTOR;
    new->genestrand = Substring_genestrand(acceptor);
    new->chrnum = Substring_chrnum(acceptor);
    new->chroffset = Substring_chroffset(acceptor);
    new->chrhigh = Substring_chrhigh(acceptor);
    new->plusp = Substring_plusp(acceptor);

  } else if (acceptor == NULL) {
    new->hittype = HALFSPLICE_DONOR;
    new->genestrand = Substring_genestrand(donor);
    new->chrnum = Substring_chrnum(donor);
    new->chroffset = Substring_chroffset(donor);
    new->chrhigh = Substring_chrhigh(donor);
    new->plusp = Substring_plusp(donor);

  } else if (shortdistancep == true) {
    new->hittype = SPLICE;
    new->genestrand = Substring_genestrand(donor);
    new->chrnum = Substring_chrnum(donor);
    new->chroffset = Substring_chroffset(donor);
    new->chrhigh = Substring_chrhigh(donor);
    new->plusp = Substring_plusp(donor);

#if 0
  } else if (merge_samechr_p == false) {
    new->hittype = DISTANT_SPLICE;
    new->chrnum = 0;
    new->chroffset = 0;
    new->chrhigh = 0;
#endif

  } else {
    if (Substring_chrnum(donor) == Substring_chrnum(acceptor)) {
      new->hittype = SAMECHR_SPLICE;
      new->genestrand = Substring_genestrand(donor);
      new->chrnum = Substring_chrnum(donor);
      new->chroffset = Substring_chroffset(donor);
      new->chrhigh = Substring_chrhigh(donor);
    } else {
      new->hittype = TRANSLOC_SPLICE;
      new->genestrand = 0;
      new->chrnum = 0;
      new->chroffset = 0;
      new->chrhigh = 0;
    }

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
  debug0(printf("  hittype is %s\n",hittype_string(new->hittype)));

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

  if (new->genomicstart < new->genomicend) {
    new->plusp = true;
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
    new->plusp = false;
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

  if (new->chrnum == 0) {
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
  new->gmap_triedp = false;

  assert(new->substring1 != NULL);

  return new;
}



/* Never returns NULL.  Never copies substrings.  Always shortdistance. */
T
Stage3end_new_shortexon (int *found_score, Substring_T donor, Substring_T acceptor, Substring_T shortexon,
			 Genomicpos_T acceptor_distance, Genomicpos_T donor_distance,
			 int amb_nmatches_donor, int amb_nmatches_acceptor,
			 Intlist_T ambi_left, Intlist_T ambi_right,
			 Intlist_T amb_nmismatches_left, Intlist_T amb_nmismatches_right,
			 bool copy_donor_p, bool copy_acceptor_p, bool copy_shortexon_p,
			 int splicing_penalty, int querylength, int sensedir) {
  T new;
  int ignore;
  
  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_shortexon %p\n",new));

  new->deletion = (char *) NULL;
  new->querylength_adj = querylength;

  new->genestrand = Substring_genestrand(shortexon);
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
  new->gmap_triedp = false;

  return new;
}


T
Stage3end_new_terminal (int querystart, int queryend, int nmismatches_whole,
			Genomicpos_T left, Compress_T query_compress,
			int querylength, bool plusp, int genestrand,
			Endtype_T start_endtype, Endtype_T end_endtype,
			Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh,
			int max_mismatches_allowed) {
  T new;
  Substring_T substring;
  Genomicpos_T genomicstart, genomicend, alignstart, alignend, alignstart_trim, alignend_trim;
  bool trim_left_p, trim_right_p;

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;

    alignstart = genomicstart + querystart;
    alignend = genomicstart + queryend;

  } else {
    genomicstart = left + querylength;
    genomicend = left;
    
    alignstart = genomicstart - querystart;
    alignend = genomicstart - queryend;
  }

  if (start_endtype == TERM) {
    trim_left_p = true;
    trim_right_p = false;
  } else if (end_endtype == TERM) {
    trim_left_p = false;
    trim_right_p = true;
  } else {
    abort();
  }

  if ((substring = Substring_new(nmismatches_whole,chrnum,chroffset,chrhigh,
				 left,genomicstart,genomicend,query_compress,
				 start_endtype,end_endtype,querystart,queryend,querylength,
				 alignstart,alignend,/*genomiclength*/querylength,
				 /*extraleft*/0,/*extraright*/0,/*exactp*/false,plusp,genestrand,
				 trim_left_p,trim_right_p,/*minlength*/querylength/2)) == NULL) {
    return (T) NULL;
    
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
  }

  /* Re-compute nmismatches_whole and nmatches for terminal alignments */
  debug0(printf("Recomputing nmismatches_whole from %d ",nmismatches_whole));
  alignstart_trim = Substring_alignstart_trim(substring);
  alignend_trim = Substring_alignend_trim(substring);

  if (plusp == true) {
    nmismatches_whole = 
      Genome_count_mismatches_substring(query_compress,left,/*pos5*/alignstart_trim-left,
					/*pos3*/alignend_trim-left,/*plusp*/true,genestrand);
  } else {
    nmismatches_whole = 
      Genome_count_mismatches_substring(query_compress,left,/*pos5*/alignend_trim-left,
					/*pos3*/alignstart_trim-left,/*plusp*/false,genestrand);
  }
  debug0(printf("to %d\n",nmismatches_whole));

  if (nmismatches_whole > max_mismatches_allowed) {
    Substring_free(&substring);
    return (T) NULL;
  } else {
    Substring_set_nmismatches_terminal(substring,nmismatches_whole);
  }

  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_terminal %p: left %u, genomicstart/end %u..%u, chrhigh %u, chrnum %d, querystart %d, queryend %d\n",
		new,left,genomicstart,genomicend,chrhigh,chrnum,querystart,queryend));

  new->substring1 = substring;
  new->substring2 = (Substring_T) NULL;
  new->substring0 = (Substring_T) NULL;
  new->substring_donor = new->substring_acceptor = (Substring_T) NULL;
  new->substringD = new->substringA = (Substring_T) NULL;
  new->substring_low = new->substring_high = new->substring1;

  new->pairarray = (struct Pair_T *) NULL;

  new->deletion = (char *) NULL;
  new->querylength_adj = querylength;
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  if (genomicstart < genomicend) {
    new->low = genomicstart;
    new->high = genomicend;
  } else {
    new->low = genomicend;
    new->high = genomicstart;
  }

  new->hittype = TERMINAL;
  new->genestrand = genestrand;

  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
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
  new->gmap_triedp = false;

  return new;
}


T
Stage3end_new_gmap (int nmismatches_whole, int nmatches_pretrim, int nmatches_posttrim,
		    int ambig_end_length_5, int ambig_end_length_3,
		    Splicetype_T ambig_splicetype_5, Splicetype_T ambig_splicetype_3,
		    struct Pair_T *pairarray, int npairs, int nsegments, int nintrons, int nindelbreaks,
		    Genomicpos_T left, int genomiclength, bool plusp, int genestrand, int querylength,
		    Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh, int sensedir) {
  T new;
  Genomicpos_T genomicstart, genomicend;

  /* Removed statements to return NULL, because GMAP alignments seem
     to be okay, at least when starting before coordinate 0 */
  /* Example (when aligned to chrM at beginning of genome) (actually aligns circularly):
GGATGAGGCAGGAATCAAAGACAGATACTGCGACATAGGGTGCTCCGGCTCCAGCGTCTCGCAATGCTATCGCGTG
ATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAG
  */

  if (plusp == true) {
    genomicstart = left;
    if ((genomicend = left + genomiclength) > chrhigh) {
      /* return (T) NULL; */
    }
    if (genomicstart > genomicend) {
      /* Must have started before coordinate 0 */
      debug0(printf("plusp and genomicstart %u > genomicend %u => started before coordinate 0\n",
		    genomicstart,genomicend));
      /* return (T) NULL; */
    }
  } else {
    if ((genomicstart = left + genomiclength) > chrhigh) {
      /* return (T) NULL; */
    }
    genomicend = left;
    if (genomicend > genomicstart) {
      /* Must have started before coordinate 0 */
      debug0(printf("minusp and genomicend %u > genomicstart %u => started before coordinate 0\n",
		    genomicend,genomicstart));
      /* return (T) NULL; */
    }
  }

  new = (T) MALLOC_OUT(sizeof(*new));
  debug0(printf("Stage3end_new_gmap %p: left %u, genomicstart/end %u..%u, chrhigh %u, chrnum %d, nmismatches %d, sensedir %d\n",
		new,left,genomicstart,genomicend,chrhigh,chrnum,nmismatches_whole,sensedir));

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
  new->querylength_adj = querylength /* - nindels */;
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  if (genomicstart < genomicend) {
    new->low = genomicstart;
    new->high = genomicend;
  } else {
    new->low = genomicend;
    new->high = genomicstart;
  }

  new->hittype = GMAP;
  new->genestrand = genestrand;

  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;
  new->plusp = plusp;
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

#if 0
  /* This favors the trimmed results */
  new->score = nmismatches_whole;
  new->score += localsplicing_penalty * nintrons;
  new->score += indel_penalty_middle * nindelbreaks;
  debug0(printf("gmap score = %d = %d + %d*%d + %d*%d\n",
		new->score,nmismatches_whole,localsplicing_penalty,nintrons,indel_penalty_middle,nindelbreaks));
#else
  /* This is a better way to score GMAP.  Using nmatches_pretrim puts all GMAP entries on an even level. */
  new->score = querylength - nmatches_posttrim;
  new->score += localsplicing_penalty * nintrons;
  new->score += indel_penalty_middle * nindelbreaks;
  debug0(printf("gmap score = %d = querylength %d - posttrim %d (pretrim %d) + %d*%d + %d*%d\n",
		new->score,querylength,nmatches_posttrim,nmatches_pretrim,
		localsplicing_penalty,nintrons,indel_penalty_middle,nindelbreaks));
#endif


  new->nmismatches_bothdiff = nmismatches_whole;
  /* new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(new->substring1); */
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1); */

  new->nmatches = nmatches_pretrim;
  new->nmatches_posttrim = nmatches_posttrim;
  if (favor_ambiguous_p == true) {
    new->nmatches += ambig_end_length_5 + ambig_end_length_3;
  }

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
  new->gmap_triedp = false;	/* if set to true, prevents running GMAP against GMAP */

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
  } else {
    return 0;
  }
}



static double
Stage3end_compute_mapq (Stage3end_T this, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			char *quality_string) {

  if (this == NULL) {
    return 0.0;

  } else if (this->hittype == GMAP) {
    this->mapq_loglik = Pair_compute_mapq(this->pairarray,this->npairs,quality_string);

  } else if (this->plusp == true) {
    this->mapq_loglik =
      Substring_compute_mapq(this->substring1,query_compress_fwd,quality_string);

    if (this->substring2 != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substring2,query_compress_fwd,
			       quality_string);
    }
    if (this->substring0 != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substring0,query_compress_fwd,
			       quality_string);
    }

  } else {
    this->mapq_loglik =
      Substring_compute_mapq(this->substring1,query_compress_rev,
			     quality_string);

    if (this->substring2 != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substring2,query_compress_rev,
			       quality_string);
    }
    if (this->substring0 != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substring0,query_compress_rev,
			       quality_string);
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




Stage3end_T *
Stage3end_eval_and_sort (int *npaths, int *second_absmq,
			 Stage3end_T *stage3array, int maxpaths, Shortread_T queryseq,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 Genome_T genome, char *quality_string, bool displayp) {
  char *query;
  double maxlik, loglik;
  double total, q;		/* For Bayesian mapq calculation */
  int compute_npaths;
  int i;

  if (*npaths == 0) {
    /* Skip */

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
    *second_absmq = 0;

  } else {
    /* Compute mapq_loglik */
    for (i = 0; i < *npaths; i++) {
      Stage3end_compute_mapq(stage3array[i],query_compress_fwd,query_compress_rev,
			     quality_string);
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
    *second_absmq = stage3array[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < *npaths; i++) {
      total += (stage3array[i]->mapq_loglik = exp(stage3array[i]->mapq_loglik));
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
Stage3end_trim_left (T this) {
  if (this->substring0 != NULL) {
    return Substring_trim_left(this->substring0);
  } else {
    return Substring_trim_left(this->substring1);
  }
}

static int
Stage3end_trim_right (T this) {
  if (this->substring2 != NULL) {
    return Substring_trim_right(this->substring2);
  } else {
    return Substring_trim_right(this->substring1);
  }
}



/* Note: single-end terminals can be present with non-terminals when
   paired-end reads are searched for concordance, which can accumulate
   terminal alignments */
List_T
Stage3end_optimal_score (List_T hitlist, int cutoff_level, int suboptimal_mismatches,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 bool keep_gmap_p) {
  List_T optimal = NULL, p;
  T hit;
  int n;
  int minscore = MAX_READLENGTH;
  int max_nmatches = 0, max_nmatches_posttrim = 0;
  bool non_translocation_p = false;
  bool non_terminal_p = false, terminalp = false;
  int max_trim_left = 0, max_trim_right = 0, trim;

  n = List_length(hitlist);
  if (n == 0) {
    return NULL;
  }

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->chrnum != 0) {
      non_translocation_p = true;
    }
    if (hit->hittype == GMAP) {
      /* Skip */
    } else if (hit->hittype == TERMINAL) {
      debug4(printf("Found a terminal\n"));
      terminalp = true;
    } else {
      debug4(printf("Found a non-terminal\n"));
      non_terminal_p = true;
    }
    hit->score_eventrim = hit->score;
  }


  if (terminalp == true && non_terminal_p == true) {
    /* Use eventrim for comparing terminals with non-terminals */
    for (p = hitlist; p != NULL; p = p->rest) {
      hit = (T) p->first;
      if (hit->hittype == GMAP) {
	/* Skip */
      } else if (hit->hittype == TERMINAL) {
	/* Skip: Play by non-terminal rules */
      } else {
	debug4(printf("trim_left: %d, trim_right %d\n",Stage3end_trim_left(hit),Stage3end_trim_right(hit)));
	if ((trim = Stage3end_trim_left(hit)) > max_trim_left) {
	  max_trim_left = trim;
	}
	if ((trim = Stage3end_trim_right(hit)) > max_trim_right) {
	  max_trim_right = trim;
	}
      }
    }
    debug4(printf("max_trim_left: %d, max_trim_right %d\n",max_trim_left,max_trim_right));

    for (p = hitlist; p != NULL; p = p->rest) {
      hit = (T) p->first;
      if (hit->hittype == GMAP) {
	/* Skip */
      } else {
	hit->score_eventrim = Substring_count_mismatches_region(hit->substring0,max_trim_left,max_trim_right,
								query_compress_fwd,query_compress_rev);
	hit->score_eventrim += Substring_count_mismatches_region(hit->substring1,max_trim_left,max_trim_right,
								 query_compress_fwd,query_compress_rev);
	hit->score_eventrim += Substring_count_mismatches_region(hit->substring2,max_trim_left,max_trim_right,
								 query_compress_fwd,query_compress_rev);
	debug4(printf("score_eventrim = %d\n",hit->score_eventrim));
      }
    }
  }

  /* Compute minscore */
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->chrnum == 0 && non_translocation_p == true) {
      /* Skip, since we will eliminate */
    } else {
      if (hit->nmatches > max_nmatches) {
	max_nmatches = hit->nmatches;
	max_nmatches_posttrim = hit->nmatches_posttrim;
      }
#ifdef TERMINAL_SECOND_CLASS
      if (non_terminal_p == true && hit->hittype == TERMINAL) {
	/* Skip from setting minscore */
      }
#endif
      if (hit->score_eventrim <= cutoff_level) {
	if (hit->score_eventrim < minscore) {
	  minscore = hit->score_eventrim;
	}
      }
    }
  }

  debug4(printf("Stage3end_optimal_score over %d hits: minscore = %d + subopt:%d\n",
		n,minscore,suboptimal_mismatches));
  minscore += suboptimal_mismatches;
  max_nmatches -= suboptimal_mismatches;
  max_nmatches_posttrim -= suboptimal_mismatches;

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;

    if (keep_gmap_p == true && hit->hittype == GMAP) {
      /* GMAP hits already found to be better than their corresponding terminals */
      debug4(printf("Keeping a hit of type GMAP\n"));
      optimal = List_push(optimal,hit);

    } else if (hit->chrnum == 0 && non_translocation_p == true) {
      debug4(printf("Eliminating a hit with splice translocation\n"));
      Stage3end_free(&hit);

#ifdef TERMINAL_SECOND_CLASS
    } else if (hit->hittype == TERMINAL && non_terminal_p == true) {
      if (hit->nmatches >= max_nmatches) {
	debug4(printf("Keeping a terminal with nmatches %d\n",hit->nmatches));
	optimal = List_push(optimal,(void *) hit);
      } else {
	debug4(printf("Eliminating a terminal where non-terminals are present\n"));
	Stage3end_free(&hit);
      }
#endif

    } else if (hit->score_eventrim > cutoff_level) {
      /* For dibasep were previously using hit->ntscore, but gives false positives */
      debug4(printf("Eliminating a hit of type %s with score %d > cutoff_level %d\n",
		    hittype_string(hit->hittype),hit->score_eventrim,cutoff_level));
      Stage3end_free(&hit);

    } else if (hit->score_eventrim > minscore && hit->nmatches_posttrim < max_nmatches_posttrim) {
      debug4(printf("Eliminating a hit with score %d and type %s\n",
		    hit->score_eventrim,hittype_string(hit->hittype)));
      Stage3end_free(&hit);

    } else {
      debug4(printf("Keeping a hit with score %d and type %s\n",
		    hit->score_eventrim,hittype_string(hit->hittype)));
      optimal = List_push(optimal,hit);
    }
  }
  
  List_free(&hitlist);

  debug4(printf("hitlist now has %d entries\n",List_length(optimal)));
  return optimal;
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
  Genomicpos_T splice_distance_1, splice_distance_2;

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
	     printf("  %d: %u..%u (plusp = %d) ambiguousp:%d\n",
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
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches: %d, ",
		  i,hittype_string(x->hittype),x,x->chrnum,x->genomicstart-x->chroffset,x->genomicend-x->chroffset,x->nmatches);
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
	debug7(printf("  Keeping #%d:%u..%u, nmatches %d (nindels %d, indel_pos %d, distance %u, chrnum %d) (plusp = %d)\n",
		      x->chrnum,x->genomicstart-x->chroffset,x->genomicend-x->chroffset,
		      x->nmatches,x->nindels,x->indel_pos,x->distance,x->chrnum,x->plusp));
	hitlist = List_push(hitlist,x);
      } else {
	debug7(printf("  Eliminating #%d:%u..%u, nmatches %d (nindels %d, indel_pos %d, distance %u, chrnum %d) (plusp = %d)\n",
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
	   printf("  Final %d: #%d:%u..%u (plusp = %d)\n",
		  i,x->chrnum,x->genomicstart-x->chroffset,x->genomicend-x->chroffset,x->plusp);
	 }
	 );

  debug7(printf("Exited Stage3end_remove_duplicates with %d hits\n",List_length(hitlist)));
  return hitlist;
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

  } else {
    return 0;
  }
}


static int
hit_goodness_cmp (Stage3end_T hit,
#ifdef DEBUG7
		  int k,
#endif
		  Stage3end_T best_hit) {

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

  if (hit->hittype == TERMINAL && best_hit->hittype != TERMINAL) {
    debug7(printf(" => %d loses by terminal\n",k));
    return -1;
  } else if (hit->hittype != TERMINAL && best_hit->hittype == TERMINAL) {
    debug7(printf(" => %d wins by terminal\n",k));
    return +1;

  } else if (hit->score > best_hit->score) {
    debug7(printf("  => %d loses by score\n",k));
    return -1;
  } else if (hit->score < best_hit->score) {
    debug7(printf("  => %d wins by score\n",k));
    return +1;

  } else if (hit->nmatches < best_hit->nmatches) {
    debug7(printf("  => %d loses by nmatches\n",k));
    return -1;
  } else if (hit->nmatches > best_hit->nmatches) {
    debug7(printf("  => %d wins by nmatches\n",k));
    return +1;

  } else if (hit->nmatches_posttrim < best_hit->nmatches_posttrim) {
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

  } else {
    return 0;
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
hit_bad_superstretch_p (Stage3end_T hit_k, Stage3end_T *hits, int k, int j) {
  int a;

  for (a = k+1; a <= j; a++) {
    if (hit_subsumption(hit_k,hits[a]) == true) {
      debug8(printf("Testing %d because stretches over %d",k,a));
      if (hit_goodness_cmp(hits[a],
#ifdef DEBUG7
			   a,
#endif
			   hit_k) >= 0) {
	debug7(printf(" => eliminating\n"));
	return true;
      }
      debug7(printf("\n"));
    }
  }
  return false;
}



List_T
Stage3end_remove_overlaps (List_T hitlist) {
  List_T unique = NULL;
  T best_hit, hit, *hits, *prev;
  int cmp;
  int nkept, n, i, j, k, besti;
  bool *eliminate;
#ifdef PRE_RESOLVE_MULTIMAPPING
  long int best_tally;
#endif

  n = List_length(hitlist);
  debug7(printf("Entered Stage3end_remove_overlaps with %d hits\n",n));
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
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches: %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,hit->nmatches);
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
      debug7(printf("  Keeping %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug7(printf("  Already paired %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else {
      debug7(printf("  Eliminating %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
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
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches: %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,hit->nmatches);
	 }
	 );

  /* Find clusters */
  i = 0;
  while (i < n) {
    j = i;
    while (j+1 < n && hit_subsumption(hits[i],hits[j+1]) == true) {
      j = j+1;
    }

    if (j > i) {
      debug7(printf("Cluster from %d up through %d\n",i,j));

      /* Find bad superstretches */
      for (k = i; k <= j; k++) {
	if (hit_bad_superstretch_p(hits[k],hits,k,j) == true) {
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
      debug7(printf("  Keeping %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug7(printf("  Already paired %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else {
      debug7(printf("  Eliminating %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
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
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches: %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,hit->nmatches);
	 }
	 );

  /* Find clusters from left */
  i = 0;
  while (i < n) {
    j = i;
    while (j+1 < n && hit_subsumption(hits[i],hits[j+1]) == true) {
      j = j+1;
    }

    if (j > i) {
      debug7(printf("Cluster from %d up through %d\n",i,j));

      best_hit = hits[i];
      besti = i;
      debug7(printf("Assume best is %d\n",besti));

      for (k = i+1; k <= j; k++) {
	cmp = hit_goodness_cmp(hits[k],
#ifdef DEBUG7
			       k,
#endif
			       best_hit);
	debug7(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	if (cmp > 0) {
	  best_hit = hits[k];
	  besti = k;
	  debug7(printf("Best is now %d\n",besti));
	}
      }

      for (k = i; k <= j; k++) {
	if (k == besti) {
	  /* Skip */
	} else if (hit_goodness_cmp(hits[k],
#ifdef DEBUG7
				    k,
#endif
				    best_hit) <= 0) {
	  debug7(printf("  Eliminating hit %d from left, because beaten by %d\n",k,besti));
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
    while (i-1 >= 0 && hit_subsumption(hits[j],hits[i-1]) == true) {
      i = i-1;
    }

    if (i < j) {
      debug7(printf("Cluster from %d down through %d\n",j,i));
      best_hit = hits[i];
      besti = i;
      debug7(printf("Assume best is %d\n",besti));

      for (k = i+1; k <= j; k++) {
	cmp = hit_goodness_cmp(hits[k],
#ifdef DEBUG7
			       k,
#endif
			       best_hit);
	debug7(printf("Comparison of %d with best %d yields %d\n",k,besti,cmp));
	if (cmp > 0) {
	  best_hit = hits[k];
	  besti = k;
	  debug7(printf("Best is now %d\n",besti));
	}
      }

      for (k = i; k <= j; k++) {
	if (k == besti) {
	  /* Skip */
	} else if (hit_goodness_cmp(hits[k],
#ifdef DEBUG7
				    k,
#endif
				    best_hit) <= 0) {
	  debug7(printf("  Eliminating hit %d from right, because beaten by %d\n",k,besti));
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
      debug7(printf("  Keeping %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else if (hit->paired_usedp == true) {
      debug7(printf("  Already paired %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
		    hit->low,hit->high,hit->nmatches,hit->plusp));
      hits[j++] = hit;
    } else {
      debug7(printf("  Eliminating %u..%u, nmatches (trimmed) %d (plusp = %d)\n",
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
	   printf("  Initial %d (%s): %p #%d:%u..%u, nmatches: %d\n",
		  i,hittype_string(hit->hittype),hit,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,hit->nmatches);
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
print_alignment_info (FILE *fp, int nblocks, int score, int mapq_score) {
  fprintf(fp,"segs:%d,align_score:%d,mapq:%d",nblocks,score,mapq_score);
  return;
}


Pairtype_T
Stage3_determine_pairtype (T hit5, T hit3) {
  if (hit5->chrnum != hit3->chrnum) {
    return UNPAIRED;
  } else if (hit5->plusp != hit3->plusp) {
    return PAIRED_INVERSION;
  } else if (hit5->plusp == true) {
    if (hit3->genomicend < hit5->genomicstart) {
      return PAIRED_SCRAMBLE;
    } else if (hit3->genomicstart > hit5->genomicend + pairmax) {
      return PAIRED_TOOLONG;
    } else {
      return CONCORDANT;
    }
  } else {
    if (hit3->genomicend > hit5->genomicstart) {
      return PAIRED_SCRAMBLE;
    } else if (hit3->genomicstart + pairmax < hit5->genomicend) {
      return PAIRED_TOOLONG;
    } else {
      return CONCORDANT;
    }
  }
}


Pairtype_T
Stage3pair_pairtype (Stage3pair_T this) {
  return this->pairtype;
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
  case TRANSLOCATION: break;
  case PAIRED_UNSPECIFIED: abort();
  case UNPAIRED: abort();
  }

  return;
}

static void
print_single (FILE *fp, T this, int score, IIT_T chromosome_iit, Shortread_T queryseq,
	      bool invertp, T hit5, T hit3, int insertlength, int pairscore, 
	      Pairtype_T pairtype, int mapq_score) {
  char *chr;
  bool allocp;

  chr = IIT_label(chromosome_iit,this->chrnum,&allocp);

  fprintf(fp," ");
  Substring_print_single(fp,this->substring1,queryseq,chr,invertp);

  /* Alignment info */
  fprintf(fp,"\t");
  print_alignment_info(fp,/*nblocks*/1,score,mapq_score);

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
print_insertion (FILE *fp, T this, int score, IIT_T chromosome_iit, Shortread_T queryseq,
		 bool invertp, T hit5, T hit3, int insertlength, int pairscore,
		 Pairtype_T pairtype, int mapq_score) {
  char *chr;
  bool allocp;

  chr = IIT_label(chromosome_iit,this->chrnum,&allocp);

  fprintf(fp," ");
  Substring_print_insertion_1(fp,this->substring1,this->substring2,this->nindels,
			      queryseq,chr,invertp);
  /* Alignment info */
  fprintf(fp,"\t");
  print_alignment_info(fp,/*nblocks*/2,score,mapq_score);

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
print_deletion (FILE *fp, T this, int score, IIT_T chromosome_iit, Shortread_T queryseq,
		bool invertp, T hit5, T hit3, int insertlength, int pairscore,
		Pairtype_T pairtype, int mapq_score) {
  char *chr;
  bool allocp;

  chr = IIT_label(chromosome_iit,this->chrnum,&allocp);

  fprintf(fp," ");
  Substring_print_deletion_1(fp,this->substring1,this->substring2,this->nindels,this->deletion,
			     queryseq,chr,invertp);
  /* Alignment info */
  fprintf(fp,"\t");
  print_alignment_info(fp,/*nblocks*/2,score,mapq_score);

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
	      IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, T hit5, T hit3,
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
    print_alignment_info(fp,/*nblocks*/1,score,mapq_score);

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
    print_alignment_info(fp,/*nblocks*/1,score,mapq_score);

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
    print_alignment_info(fp,/*nblocks*/2,score,mapq_score);

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
    print_alignment_info(fp,/*nblocks*/2,score,mapq_score);

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
    print_alignment_info(fp,/*nblocks*/2,score,mapq_score);

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
    print_alignment_info(fp,/*nblocks*/2,score,mapq_score);

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
		 IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, T hit5, T hit3,
		 int insertlength, int pairscore, Pairtype_T pairtype, int mapq_score) {
  Substring_T donor, acceptor, shortexon;
  Genomicpos_T distance1, distance2;
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
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score);
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
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score);
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
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score);
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
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score);
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
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score);
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
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score);
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
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score);
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
      fprintf(fp,"\t"); print_alignment_info(fp,nblocks,score,mapq_score);
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
		 IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, T hit5, T hit3, int insertlength,
		 int pairscore, Pairtype_T pairtype, int mapq_score) {
  int cdna_direction;

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
    if (this->sensedir == SENSE_FORWARD) {
      cdna_direction = +1;
    } else if (this->sensedir == SENSE_ANTI) {
      cdna_direction = -1;
    } else {
      cdna_direction = 0;
    }

    if (Shortread_invertedp(queryseq) == false) {
      Substring_print_gmap(fp,this->pairarray,this->npairs,this->nsegments,/*invertedp*/false,
			   this->gmap_start_endtype,this->gmap_end_endtype,
			   this->chrnum,this->chroffset,this->chrhigh,Shortread_fulllength(queryseq),
			   this->plusp,cdna_direction,this->score,insertlength,pairscore,mapq_score,
			   chromosome_iit);
    } else {
      Substring_print_gmap(fp,this->pairarray,this->npairs,this->nsegments,/*invertedp*/true,
			   this->gmap_end_endtype,this->gmap_start_endtype,
			   this->chrnum,this->chroffset,this->chrhigh,Shortread_fulllength(queryseq),
			   this->plusp,cdna_direction,this->score,insertlength,pairscore,mapq_score,
			   chromosome_iit);
    }

  } else {
    abort();
  }

  return;
}


#if 0
static Genomicpos_T
compute_pairlength (Genomicpos_T start1, Genomicpos_T end1, 
		    Genomicpos_T start2, Genomicpos_T end2) {
  Genomicpos_T low1, high1, low2, high2;

  if (start1 < end1) {
    low1 = start1;
    high1 = end1;
  } else {
    low1 = end1;
    high1 = start1;
  }

  if (start2 < end2) {
    low2 = start2;
    high2 = end2;
  } else {
    low2 = end2;
    high2 = start2;
  }

  if (low2 > high1) {
    return low2 - high1;
  } else if (low1 > high2) {
    return low1 - high2;
  } else {
    return 0U;
  }
}
#endif


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
		      char initchar, bool firstp, IIT_T chromosome_iit,
		      Shortread_T queryseq, Shortread_T headerseq,
		      int maxpaths, bool quiet_if_excessive_p,
		      bool invertp, int quality_shift,
		      FILE *fp_nomapping_1, FILE *fp_unpaired_uniq, FILE *fp_unpaired_transloc, FILE *fp_unpaired_mult,
		      FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_transloc, FILE *fp_halfmapping_mult,
		      FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr, FILE *fp_paired_uniq_long,
		      FILE *fp_paired_mult, FILE *fp_concordant_uniq, FILE *fp_concordant_transloc, FILE *fp_concordant_mult) {
  Stage3pair_T *stage3pairarray, stage3pair;
  T *stage3array, this, hit5, hit3;
  int npaths, pathnum, second_absmq;
  bool outputp, translocationp;
  FILE *fp;

  if (resulttype == PAIREDEND_NOMAPPING) {
    /* If fails_as_input_p == true, then this case is handled by calling procedure */
    print_query_header(fp_nomapping_1,initchar,queryseq,invertp);
    fprintf(fp_nomapping_1,"\t0 %s",UNPAIRED_TEXT);

    print_barcode_and_quality(fp_nomapping_1,queryseq,invertp,quality_shift);
    
    fprintf(fp_nomapping_1,"\t");
    Shortread_print_header(fp_nomapping_1,headerseq);
    fprintf(fp_nomapping_1,"\n");

  } else if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC) {
    if (resulttype == CONCORDANT_TRANSLOC) {
      fp = fp_concordant_transloc;
    } else {
      fp = fp_concordant_uniq;
    }

    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&second_absmq,result);

    print_query_header(fp,initchar,queryseq,invertp);
    fprintf(fp,"\t1 %s",CONCORDANT_TEXT);
    if (resulttype == CONCORDANT_TRANSLOC) {
      fprintf(fp," (transloc)");
    }
    
    print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

    fprintf(fp,"\t");
    Shortread_print_header(fp,headerseq);

    stage3pair = stage3pairarray[0];
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

  } else if (resulttype == CONCORDANT_MULT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&second_absmq,result);

    print_query_header(fp_concordant_mult,initchar,queryseq,invertp);
    fprintf(fp_concordant_mult,"\t%d %s",npaths,CONCORDANT_TEXT);

    print_barcode_and_quality(fp_concordant_mult,queryseq,invertp,quality_shift);

    fprintf(fp_concordant_mult,"\t");
    Shortread_print_header(fp_concordant_mult,headerseq);

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
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&second_absmq,result);
    stage3pair = stage3pairarray[0];

    if (stage3pair->pairtype == PAIRED_INVERSION) {
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
    Shortread_print_header(fp,headerseq);

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
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&second_absmq,result);

    print_query_header(fp_paired_mult,initchar,queryseq,invertp);
    fprintf(fp_paired_mult,"\t%d %s",npaths,PAIRED_TEXT);

    print_barcode_and_quality(fp_paired_mult,queryseq,invertp,quality_shift);

    fprintf(fp_paired_mult,"\t");
    Shortread_print_header(fp_paired_mult,headerseq);

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
      stage3array = (T *) Result_array(&npaths,&second_absmq,result);
    } else {
      stage3array = (T *) Result_array2(&npaths,&second_absmq,result);
    }

    outputp = true;
    translocationp = false;
    if (resulttype == HALFMAPPING_UNIQ) {
      fp = fp_halfmapping_uniq;

    } else if (resulttype == HALFMAPPING_TRANSLOC) {
      fp = fp_halfmapping_transloc;
      translocationp = true;

    } else if (resulttype == HALFMAPPING_MULT) {
      fp = fp_halfmapping_mult;
      if (quiet_if_excessive_p && npaths > maxpaths) {
	outputp = false;
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      fp = fp_unpaired_uniq;

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
      stage3array = (T *) Result_array(&npaths,&second_absmq,result);
      hit5 = stage3array[0];
      stage3array = (T *) Result_array2(&npaths,&second_absmq,result);
      hit3 = stage3array[0];
      fprintf(fp," (%s)",unpaired_type_text(hit5,hit3));
    }
#endif

    print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

    fprintf(fp,"\t");
    Shortread_print_header(fp,headerseq);

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
		  IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  int maxpaths, bool quiet_if_excessive_p,
		  bool nofailsp, bool failsonlyp, bool fails_as_input_p,
		  bool fastq_format_p, int quality_shift,
		  FILE *fp_nomapping_1, FILE *fp_nomapping_2,
		  FILE *fp_unpaired_uniq, FILE *fp_unpaired_transloc, FILE *fp_unpaired_mult,
		  FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_transloc, FILE *fp_halfmapping_mult,
		  FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr,
		  FILE *fp_paired_uniq_long, FILE *fp_paired_mult,
		  FILE *fp_concordant_uniq, FILE *fp_concordant_transloc, FILE *fp_concordant_mult) {

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
      print_one_paired_end(result,resulttype,'>',/*firstp*/true,
			   chromosome_iit,queryseq1,/*headerseq*/queryseq1,
			   maxpaths,quiet_if_excessive_p,invert_first_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_transloc,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_transloc,fp_halfmapping_mult,
			   fp_paired_uniq_inv,fp_paired_uniq_scr,fp_paired_uniq_long,
			   fp_paired_mult,fp_concordant_uniq,fp_concordant_transloc,fp_concordant_mult);

      /* Second end */
      print_one_paired_end(result,resulttype,'<',/*firstp*/false,
			   chromosome_iit,queryseq2,/*headerseq*/queryseq1,
			   maxpaths,quiet_if_excessive_p,invert_second_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_transloc,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_transloc,fp_halfmapping_mult,
			   fp_paired_uniq_inv,fp_paired_uniq_scr,fp_paired_uniq_long,
			   fp_paired_mult,fp_concordant_uniq,fp_concordant_transloc,fp_concordant_mult);
    }

  } else {
    if (failsonlyp == true) {
      /* Unwanted success: skip */
      debug1(printf("  failsonlyp is true, so no output\n"));
    
    } else {

      /* First end */
      print_one_paired_end(result,resulttype,'>',/*firstp*/true,
			   chromosome_iit,queryseq1,/*headerseq*/queryseq1,
			   maxpaths,quiet_if_excessive_p,invert_first_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_transloc,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_transloc,fp_halfmapping_mult,
			   fp_paired_uniq_inv,fp_paired_uniq_scr,fp_paired_uniq_long,
			   fp_paired_mult,fp_concordant_uniq,fp_concordant_transloc,fp_concordant_mult);

      /* Second end */
      print_one_paired_end(result,resulttype,'<',/*firstp*/false,
			   chromosome_iit,queryseq2,/*headerseq*/queryseq1,
			   maxpaths,quiet_if_excessive_p,invert_second_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_transloc,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_transloc,fp_halfmapping_mult,
			   fp_paired_uniq_inv,fp_paired_uniq_scr,fp_paired_uniq_long,
			   fp_paired_mult,fp_concordant_uniq,fp_concordant_transloc,fp_concordant_mult);
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



static Genomicpos_T
pair_insert_length (Stage3end_T hit5, Stage3end_T hit3) {

  if (hit5->plusp != hit3->plusp) {
    debug10(printf("pair_insert_length: hit5->plusp %d != hit3->plusp %d, so returning 0\n",
		   hit5->plusp,hit3->plusp));
    return 0;
  }

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

  /* No overlap found between any combination of substrings */
  if (hit5->plusp == true) {
    if (hit5->genomicend > hit3->genomicstart + hit5->querylength_adj + hit3->querylength_adj) {
      debug10(printf("pair_insert_length: no overlap found, and %u - %u + %d + %d < 0, so returning 0\n",
		     hit3->genomicstart,hit5->genomicend,hit5->querylength_adj,hit3->querylength_adj));
      return 0;
    } else {
      debug10(printf("pair_insert_length: no overlap found, so returning %u - %u + %d + %d\n",
		     hit3->genomicstart,hit5->genomicend,hit5->querylength_adj,hit3->querylength_adj));
    }
    return hit3->genomicstart - hit5->genomicend + hit5->querylength_adj + hit3->querylength_adj;
  } else {
    if (hit3->genomicstart > hit5->genomicend + hit5->querylength_adj + hit3->querylength_adj) {
      debug10(printf("pair_insert_length: no overlap found, and %u - %u + %d + %d < 0, so returning 0\n",
		     hit5->genomicend,hit3->genomicstart,hit5->querylength_adj,hit3->querylength_adj));
      return 0;
    } else {
      debug10(printf("pair_insert_length: no overlap found, so returning %u - %u + %d + %d\n",
		     hit5->genomicend,hit3->genomicstart,hit5->querylength_adj,hit3->querylength_adj));
      return hit5->genomicend - hit3->genomicstart + hit5->querylength_adj + hit3->querylength_adj;
    }
  }
}


#if 0
static Genomicpos_T
overlap_gmap_plus_old (int *querypos, Genomicpos_T *genomicstart, Genomicpos_T *genomicend,
		       Stage3end_T hit, Stage3end_T gmap) {
  Genomicpos_T genomicpos;

  debug10(printf("Entered overlap_gmap_plus\n"));
  *genomicstart = Substring_chrstart(hit->substring1);
  *genomicend = Substring_chrend(hit->substring1);
  if ((genomicpos = Pair_binary_search_ascending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
						 *genomicstart,*genomicend)) > 0U) {
    debug10(printf("substring1 %u..%u => genomicpos %u, querypos %d\n",
		   *genomicstart,*genomicend,genomicpos,*querypos));
    return genomicpos;
  }
  
  if (hit->substring2 != NULL) {
    *genomicstart = Substring_chrstart(hit->substring2);
    *genomicend = Substring_chrend(hit->substring2);
    if ((genomicpos = Pair_binary_search_ascending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
						   *genomicstart,*genomicend)) > 0U) {
      debug10(printf("substring2 %u..%u => genomicpos %u, querypos %d\n",
		     *genomicstart,*genomicend,genomicpos,*querypos));
      return genomicpos;
    }
  }

  if (hit->substring0 != NULL) {
    *genomicstart = Substring_chrstart(hit->substring0);
    *genomicend = Substring_chrend(hit->substring0);
    if ((genomicpos = Pair_binary_search_ascending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
						   *genomicstart,*genomicend)) > 0U) {
      debug10(printf("substring0 %u..%u => genomicpos %u, querypos %d\n",
		     *genomicstart,*genomicend,genomicpos,*querypos));
      return genomicpos;
    }
  }

  return 0U;
}


static Genomicpos_T
overlap_gmap_minus_old (int *querypos, Genomicpos_T *genomicstart, Genomicpos_T *genomicend,
			Stage3end_T hit, Stage3end_T gmap) {
  Genomicpos_T genomicpos;

  debug10(printf("Entered overlap_gmap_minus\n"));
  *genomicstart = Substring_chrstart(hit->substring1);
  *genomicend = Substring_chrend(hit->substring1);
  if ((genomicpos = Pair_binary_search_descending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
						  *genomicstart,*genomicend)) > 0U) {
    debug10(printf("substring1 %u..%u => genomicpos %u, querypos %d\n",
		   *genomicstart,*genomicend,genomicpos,*querypos));
    return genomicpos;
  }
  
  if (hit->substring2 != NULL) {
    *genomicstart = Substring_chrstart(hit->substring2);
    *genomicend = Substring_chrend(hit->substring2);
    if ((genomicpos = Pair_binary_search_descending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
						    *genomicstart,*genomicend)) > 0U) {
      debug10(printf("substring2 %u..%u => genomicpos %u, querypos %d\n",
		     *genomicstart,*genomicend,genomicpos,*querypos));
      return genomicpos;
    }
  }

  if (hit->substring0 != NULL) {
    *genomicstart = Substring_chrstart(hit->substring0);
    *genomicend = Substring_chrend(hit->substring0);
    if ((genomicpos = Pair_binary_search_descending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
						    *genomicstart,*genomicend)) > 0U) {
      debug10(printf("substring0 %u..%u => genomicpos %u, querypos %d\n",
		     *genomicstart,*genomicend,genomicpos,*querypos));
      return genomicpos;
    }
  }

  return 0U;
}
#endif


#if 0
/* Doesn't work if both hits cross same intron */
static Genomicpos_T
pair_insert_length_faulty (Stage3end_T hit5, Stage3end_T hit3) {

#if 0
  /* ? Doesn't hold for paired (inversion) */
  assert(hit5->plusp == hit3->plusp);
#endif

  if (hit5->plusp != hit3->plusp) {
    return 0;

  } else if (hit5->plusp == true) {
    if (Substring_overlap_p(hit5->substring_high,hit3->substring_low)) {
      return Substring_insert_length(hit5->substring_high,hit3->substring_low);
    } else {
      return 0;
    }
  } else {
    if (Substring_overlap_p(hit5->substring_low,hit3->substring_high)) {
      return Substring_insert_length(hit5->substring_low,hit3->substring_high);
    } else {
      return 0;
    }
  }
}
#endif

  


static Genomicpos_T
overlap5_gmap_plus (int *querypos, Genomicpos_T *genomicstart, Genomicpos_T *genomicend,
		    Stage3end_T hit5, Stage3end_T gmap) {
  debug10(printf("Entered overlap5_gmap_plus\n"));
  *genomicstart = Substring_chrstart(hit5->substring_high);
  *genomicend = Substring_chrend(hit5->substring_high);
  return Pair_binary_search_ascending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
				      *genomicstart,*genomicend);
}


static Genomicpos_T
overlap3_gmap_plus (int *querypos, Genomicpos_T *genomicstart, Genomicpos_T *genomicend,
		    Stage3end_T hit3, Stage3end_T gmap) {
  debug10(printf("Entered overlap3_gmap_plus\n"));
  *genomicstart = Substring_chrstart(hit3->substring_low);
  *genomicend = Substring_chrend(hit3->substring_low);
  return Pair_binary_search_ascending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
				      *genomicstart,*genomicend);
}


static Genomicpos_T
overlap5_gmap_minus (int *querypos, Genomicpos_T *genomicstart, Genomicpos_T *genomicend,
		     Stage3end_T hit5, Stage3end_T gmap) {
  debug10(printf("Entered overlap5_gmap_minus\n"));
  *genomicstart = Substring_chrstart(hit5->substring_low);
  *genomicend = Substring_chrend(hit5->substring_low);
  return Pair_binary_search_descending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
				       *genomicstart,*genomicend);
}

static Genomicpos_T
overlap3_gmap_minus (int *querypos, Genomicpos_T *genomicstart, Genomicpos_T *genomicend,
		     Stage3end_T hit3, Stage3end_T gmap) {
  debug10(printf("Entered overlap3_gmap_minus\n"));
  *genomicstart = Substring_chrstart(hit3->substring_high);
  *genomicend = Substring_chrend(hit3->substring_high);
  return Pair_binary_search_descending(&(*querypos),/*lowi*/0,/*highi*/gmap->npairs,gmap->pairarray,
				       *genomicstart,*genomicend);
}


static void
resolve_inside_ambiguous_splice_plus (int *unresolved_amb_nmatches, T *hit5, T *hit3, bool *private5p, bool *private3p,
				      Genomicpos_T *splicesites,
				      Compress_T query5_compress_fwd, Compress_T query3_compress_fwd,
				      int localsplicing_penalty, int querylength5, int querylength3,
				      int genestrand) {
#ifdef USE_BINGO
  Genomicpos_T insertlength;
#endif
  Genomicpos_T genomicstart, genomicend;
  int nbingo, bingoi5, bingoi3, nbounded, boundedi5, boundedi3, nbest, besti5, besti3, i, j;
  int best_nmismatches, nmismatches;
  bool new5p = false, new3p = false;

  Substring_T donor, acceptor, shortexon;
  Genomicpos_T segment_left;
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
	debug9(printf(" %u,%u",genomicend,genomicstart));
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
      debug9(printf(" %u",genomicend));
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
      debug9(printf(" %u",genomicstart));
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
					     Substring_chrhigh(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 1: Splice from donor #%d to acceptor #%d, with nmismatches %d\n",
		      donor_knowni,acceptor_knowni,nmismatches_shortend));
	*hit5 = Stage3end_new_shortexon(&ignore_found_score,/*donor*/(*hit5)->substringD,acceptor,shortexon,
					/*acceptor_distance*/(*hit5)->acceptor_distance,
					/*donor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
					(*hit5)->amb_nmatches_donor,/*amb_nmatches_acceptor*/0,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength5,/*sensedir*/SENSE_FORWARD);
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
				       Substring_chrhigh(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 6: Splice from antiacceptor #%d to antidonor #%d, with nmismatches %d\n",
		      acceptor_knowni,donor_knowni,nmismatches_shortend));
	*hit5 = Stage3end_new_shortexon(&ignore_found_score,donor,/*acceptor*/(*hit5)->substringA,shortexon,
					/*acceptor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
					/*donor_distance*/(*hit5)->donor_distance,
					/*amb_nmatches_donor*/0,(*hit5)->amb_nmatches_acceptor,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength5,/*sensedir*/SENSE_ANTI);
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
					   Substring_chrhigh(donor))) != NULL) {
      debug9(printf("Resolved halfsplice_donor, End 1: Splice from donor #%d to acceptor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(donor),Substring_splicesites_i(acceptor),nmismatches_shortend));
      *hit5 = Stage3end_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				   donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength5,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,/*sensedir*/SENSE_FORWARD);
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
				     Substring_chrhigh(acceptor))) != NULL) {
      debug9(printf("Resolved halfsplice_acceptor, End 6: Splice from antiacceptor #%d to antidonor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(acceptor),Substring_splicesites_i(donor),nmismatches_shortend));
      *hit5 = Stage3end_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				   donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength5,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,/*sensedir*/SENSE_ANTI);
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
					     Substring_chrhigh(shortexon))) != NULL) {
	debug9(printf("Resolved shortexonr, End 5: Splice from antidonor #%d to antiacceptor #%d, with nmismatches %d\n",
		      donor_knowni,acceptor_knowni,nmismatches_shortend));
	*hit3 = Stage3end_new_shortexon(&ignore_found_score,/*donor*/(*hit3)->substringD,acceptor,shortexon,
					/*acceptor_distance*/(*hit3)->acceptor_distance,
					/*donor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
					(*hit3)->amb_nmatches_donor,/*amb_nmatches_acceptor*/0,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength3,/*sensedir*/SENSE_ANTI);
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
				       Substring_chrhigh(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 2: Splice from acceptor #%d to donor #%d, with nmismatches %d\n",
		      acceptor_knowni,donor_knowni,nmismatches_shortend));
	*hit3 = Stage3end_new_shortexon(&ignore_found_score,donor,/*acceptor*/(*hit3)->substringA,shortexon,
					/*acceptor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
					/*donor_distance*/(*hit3)->donor_distance,
					/*amb_nmatches_donor*/0,(*hit3)->amb_nmatches_acceptor,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength3,/*sensedir*/SENSE_FORWARD);
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
					   Substring_chrhigh(donor))) != NULL) {
      debug9(printf("Resolved halfsplice donor, End 5: Splice from antidonor #%d to antiacceptor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(donor),Substring_splicesites_i(acceptor),nmismatches_shortend));
      *hit3 = Stage3end_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				   donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength3,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				   /*sensedir*/SENSE_ANTI);
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
				     Substring_chrhigh(acceptor))) != NULL) {
      debug9(printf("Resolved halfsplice acceptor, End 2: Splice from acceptor #%d (%u) to donor #%d (%u), with nmismatches %d\n",
		    Substring_splicesites_i(acceptor),splicesites[acceptor_knowni],
		    Substring_splicesites_i(donor),splicesites[donor_knowni],nmismatches_shortend));
      *hit3 = Stage3end_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				   donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength3,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				   /*sensedir*/SENSE_FORWARD);
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
				       Genomicpos_T *splicesites,
				       Compress_T query5_compress_rev, Compress_T query3_compress_rev,
				       int localsplicing_penalty, int querylength5, int querylength3,
				       int genestrand) {
#ifdef USE_BINGO
  Genomicpos_T insertlength;
#endif
  Genomicpos_T genomicstart, genomicend;
  int nbingo, bingoi5, bingoi3, nbounded, boundedi5, boundedi3, nbest, besti5, besti3, i, j;
  int best_nmismatches, nmismatches;
  bool new5p = false, new3p = false;

  Substring_T donor, acceptor, shortexon;
  Genomicpos_T segment_left;
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
	debug9(printf(" %u,%u",genomicend,genomicstart));
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
      debug9(printf(" %u",genomicend));
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
      debug9(printf(" %u",genomicstart));
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
					     Substring_chrhigh(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 3: Splice from donor #%d to acceptor #%d, with nmismatches %d\n",
		      donor_knowni,acceptor_knowni,nmismatches_shortend));
	*hit5 = Stage3end_new_shortexon(&ignore_found_score,/*donor*/(*hit5)->substringD,acceptor,shortexon,
					/*acceptor_distance*/(*hit5)->acceptor_distance,
					/*donor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
					(*hit5)->amb_nmatches_donor,/*amb_nmatches_acceptor*/0,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength5,/*sensedir*/SENSE_FORWARD);
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
				       Substring_chrhigh(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 8: Splice from antiacceptor #%d to antidonor #%d, with nmismatches_shortend %d\n",
		      acceptor_knowni,donor_knowni,nmismatches_shortend));
	*hit5 = Stage3end_new_shortexon(&ignore_found_score,donor,/*acceptor*/(*hit5)->substringA,shortexon,
					/*acceptor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
					/*donor_distance*/(*hit5)->donor_distance,
					/*amb_nmatches_donor*/0,(*hit5)->amb_nmatches_acceptor,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength5,/*sensedir*/SENSE_ANTI);
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
					   Substring_chrhigh(donor))) != NULL) {
      debug9(printf("Resolved halfsplice, End 3: Splice from donor #%d to acceptor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(donor),Substring_splicesites_i(acceptor),nmismatches_shortend));
      *hit5 = Stage3end_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				   donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength5,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,
				   /*sensedir*/SENSE_FORWARD);
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
				     Substring_chrhigh(acceptor))) != NULL) {
      debug9(printf("Resolved halfsplice acceptor, End 8: Splice from antiacceptor #%d to antidonor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(acceptor),Substring_splicesites_i(donor),nmismatches_shortend));
      *hit5 = Stage3end_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				   donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength5,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,
				   /*sensedir*/SENSE_ANTI);
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
					     Substring_chrhigh(shortexon))) != NULL) {
	debug9(printf("Resolved shortexon, End 7: Splice from antidonor #%d to antiacceptor #%d, with nmismatches %d\n",
		      donor_knowni,acceptor_knowni,nmismatches_shortend));
	*hit3 = Stage3end_new_shortexon(&ignore_found_score,/*donor*/(*hit3)->substringD,acceptor,shortexon,
					/*acceptor_distance*/(*hit3)->acceptor_distance,
					/*donor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
					(*hit3)->amb_nmatches_donor,/*amb_nmatches_acceptor*/0,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength3,/*sensedir*/SENSE_ANTI);
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
				       Substring_chrhigh(shortexon))) != NULL) {
	debug9(printf("Resolved halfsplice_acceptor, End 4: Splice from acceptor #%d to #%d, with nmismatches %d\n",
		      acceptor_knowni,donor_knowni,nmismatches_shortend));
	*hit3 = Stage3end_new_shortexon(&ignore_found_score,donor,/*acceptor*/(*hit3)->substringA,shortexon,
					/*acceptor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
					/*donor_distance*/(*hit3)->donor_distance,
					/*amb_nmatches_donor*/0,(*hit3)->amb_nmatches_acceptor,
					/*ambi_left*/NULL,/*ambi_right*/NULL,
					/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
					/*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
					localsplicing_penalty,querylength3,/*sensedir*/SENSE_FORWARD);
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
					   Substring_chrhigh(donor))) != NULL) {
      debug9(printf("Resolved halfsplice_donor, End 7: Splice from antidonor #%d to antiacceptor #%d, with nmismatches %d\n",
		    Substring_splicesites_i(donor),Substring_splicesites_i(acceptor),nmismatches_shortend));
      *hit3 = Stage3end_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				   donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength3,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				   /*sensedir*/SENSE_ANTI);
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
				     Substring_chrhigh(acceptor))) != NULL) {
      debug9(printf("Resolved halfsplice_acceptor, End 4: Splice from acceptor #%d to #%d, with nmismatches %d\n",
		    Substring_splicesites_i(acceptor),Substring_splicesites_i(donor),nmismatches_shortend));
      *hit3 = Stage3end_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				   donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				   /*shortdistancep*/true,localsplicing_penalty,querylength3,
				   /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
				   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				   /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				   /*sensedir*/SENSE_FORWARD);
      *private3p = true;
    }

  } else {
    fprintf(stderr,"Unexpected hittype %d for ambiguous end\n",(*hit3)->hittype);
    abort();
  }

  return;
}



Stage3pair_T
Stage3pair_new (T hit5, T hit3,	Genomicpos_T *splicesites,
		Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		int genestrand,	Pairtype_T pairtype, int localsplicing_penalty,
		bool private5p, bool private3p, bool expect_concordant_p) {
  Stage3pair_T new = (Stage3pair_T) MALLOC_OUT(sizeof(*new));
  Genomicpos_T genomicstart, genomicend, genomicpos;
  int querypos;
  int unresolved_amb_nmatches = 0;

  int querylength5 = hit5->querylength_adj;
  int querylength3 = hit3->querylength_adj;

  new->pairtype = pairtype;
  new->genestrand = genestrand;

#if 0
  new->mapq_loglik = hit5->mapq_loglik + hit3->mapq_loglik;
  new->mapq_score = 0;
  new->absmq_score = 0;
#endif

  debug10(printf("\nStage3pair_new called with pairtype %s and chrnum %d, %d\n",
		 pairtype_string(pairtype),hit5->chrnum,hit3->chrnum));

  if (hit5->chrnum == 0 || hit3->chrnum == 0) {
    new->dir = 0;
    new->insertlength = pair_insert_length(hit5,hit3);

  } else if (hit5->hittype == GMAP && hit3->hittype == GMAP) {
    debug10(printf("Got hit5 and hit3 both of type GMAP\n"));

    /* Do not try to resolve ambiguity on inside of concordant ends */
    if (hit5->plusp == true && hit3->plusp == true) {
      new->dir = +1;
      new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      debug10(printf("plus, no overlap: insert length %d = start3 %u - end5 %u + %d + %d\n",
		     new->insertlength,hit3->genomicstart,hit5->genomicend,querylength5,querylength3));
    } else if (hit5->plusp == false && hit3->plusp == false) {
      new->dir = -1;
      new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      debug10(printf("minus, no overlap: insert length %d = end5 %u - start3 %u + %d + %d\n",
		     new->insertlength,hit5->genomicend,hit3->genomicstart,querylength5,querylength3));
    } else {
      new->dir = 0;
      new->insertlength = 0;
    }

  } else if (hit5->hittype == GMAP) {
    debug10(printf("Got hit5 of type GMAP\n"));
    if (hit5->plusp == true && hit3->plusp == true) {
      new->dir = +1;

      /* Try to resolve ambiguity on inside of concordant ends */
      assert(expect_concordant_p == true);
      resolve_inside_ambiguous_splice_plus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					   splicesites,query5_compress_fwd,query3_compress_fwd,
					   localsplicing_penalty,querylength5,querylength3,genestrand);

      /* Have 5-start..end and 3-start..end */
      debug10(printf("plus: comparing hit5->genomicend %u <= hit3->genomicstart %u\n",
		     hit5->genomicend,hit3->genomicstart));

      if (hit5->genomicend <= hit3->genomicstart) {
	/* No overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
	debug10(printf("plus, no overlap: insert length %d = start3 %u - end5 %u + %d + %d\n",
		       new->insertlength,hit3->genomicstart,hit5->genomicend,querylength5,querylength3));
      } else if ((genomicpos = overlap3_gmap_plus(&querypos,&genomicstart,&genomicend,/*hit*/hit3,/*gmap*/hit5)) > 0U) {
	new->insertlength = /* end3 */ genomicend - /* start5 */ (genomicpos - querypos);
	debug10(printf("plus, overlap: insert length %d = end3 %u - start5 (%u - %d)\n",
		       new->insertlength,genomicend,genomicpos,querypos));
      } else {
	/* Still no overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      }

    } else if (hit5->plusp == false && hit3->plusp == false) {
      new->dir = -1;

      /* Try to resolve ambiguity on inside of concordant ends */
      assert(expect_concordant_p == true);
      resolve_inside_ambiguous_splice_minus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					    splicesites,query5_compress_rev,query3_compress_rev,
					    localsplicing_penalty,querylength5,querylength3,genestrand);

      /* Have 3-end..start and 5-end..start */
      debug10(printf("minus: comparing hit3->genomicstart %u <= hit5->genomicend %u\n",
		     hit3->genomicstart,hit5->genomicend));

      if (hit3->genomicstart <= hit5->genomicend) {
	/* No overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
	debug10(printf("minus, no overlap: insert length %d = end5 %u - start3 %u + %d + %d\n",
		       new->insertlength,hit5->genomicend,hit3->genomicstart,querylength5,querylength3));
      } else if ((genomicpos = overlap3_gmap_minus(&querypos,&genomicstart,&genomicend,/*hit*/hit3,/*gmap*/hit5)) > 0U) {
	new->insertlength = /* start5 */ (genomicpos + querypos) - /* end3 */ genomicend + 1;
	debug10(printf("minus, overlap: insert length %d = start5 (%u + %d) - end3 %u + 1\n",
		       new->insertlength,genomicpos,querypos,genomicend));
      } else {
	/* Still no overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      }
    } else {
      new->dir = 0;
      new->insertlength = 0;
    }

  } else if (hit3->hittype == GMAP) {
    debug10(printf("Got hit3 of type GMAP\n"));
    if (hit5->plusp == true && hit3->plusp == true) {
      new->dir = +1;

      /* Try to resolve ambiguity on inside of concordant ends */
      assert(expect_concordant_p == true);
      resolve_inside_ambiguous_splice_plus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					   splicesites,query5_compress_fwd,query3_compress_fwd,
					   localsplicing_penalty,querylength5,querylength3,genestrand);

      /* Have 5-start..end and 3-start..end */
      debug10(printf("plus: comparing hit5->genomicend %u <= hit3->genomicstart %u\n",
		     hit5->genomicend,hit3->genomicstart));

      if (hit5->genomicend <= hit3->genomicstart) {
	/* No overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
	debug10(printf("plus, no overlap: insert length %d = start3 %u - end5 %u + %d + %d\n",
		       new->insertlength,hit3->genomicstart,hit5->genomicend,querylength5,querylength3));
      } else if ((genomicpos = overlap5_gmap_plus(&querypos,&genomicstart,&genomicend,/*hit*/hit5,/*gmap*/hit3)) > 0U) {
	new->insertlength = /* end3 */ (genomicpos - querypos + querylength3) - /* start5 */ genomicstart;
	debug10(printf("plus, overlap: insert length %d = end3 (%u - %d + %d) - start5 %u\n",
		       new->insertlength,genomicpos,querypos,querylength3,genomicstart));
      } else {
	/* Still no overlap */
	new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      }

    } else if (hit5->plusp == false && hit3->plusp == false) {
      new->dir = -1;

      /* Try to resolve ambiguity on inside of concordant ends */
      assert(expect_concordant_p == true);
      resolve_inside_ambiguous_splice_minus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					    splicesites,query5_compress_rev,query3_compress_rev,
					    localsplicing_penalty,querylength5,querylength3,genestrand);

      /* Have 3-end..start and 5-end..start */
      debug10(printf("minus: comparing hit3->genomicstart %u <= hit5->genomicend %u\n",
		     hit3->genomicstart,hit5->genomicend));
      if (hit3->genomicstart <= hit5->genomicend) {
	/* No overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
	debug10(printf("minus, no overlap: insert length %d = end5 %u - start3 %u + %d + %d\n",
		       new->insertlength,hit5->genomicend,hit3->genomicstart,querylength5,querylength3));
      } else if ((genomicpos = overlap5_gmap_minus(&querypos,&genomicstart,&genomicend,/*hit*/hit5,/*gmap*/hit3)) > 0U) {
	new->insertlength = /* start5 */ genomicstart - /* end3 */ (genomicpos + querypos - querylength3) - 1;
	debug10(printf("minus, overlap: insert length %d = start5 %u - end3 (%u + %d - %d) - 1\n",
		       new->insertlength,genomicstart,genomicpos,querypos,querylength3));
      } else {
	/* Still no overlap */
	new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      }
    } else {
      new->dir = 0;
      new->insertlength = 0;
    }

  } else if (hit5->plusp == true && hit3->plusp == false) {
    new->dir = 0;
    
    /* Have 5-start..end and 3-end..start */
    /*   or 3-end..start and 5-start..end */

    if (hit5->genomicend < hit3->genomicend) {
      new->insertlength = (hit3->genomicend - hit5->genomicend) + querylength5 + querylength3;
    } else if (hit3->genomicstart < hit5->genomicstart) {
      new->insertlength = (hit5->genomicstart - hit3->genomicstart) + querylength5 + querylength3;
    } else {
      new->insertlength = 0;
    }

  } else if (hit5->plusp == false && hit3->plusp == true) {
    new->dir = 0;
    
    /* Have 5-end..start and 3-start..end */
    /*   or 3-start..end and 5-end..start */

    if (hit5->genomicstart < hit3->genomicstart) {
      new->insertlength = (hit3->genomicstart - hit5->genomicstart) + querylength5 + querylength3;
    } else if (hit3->genomicend < hit5->genomicend) {
      new->insertlength = (hit5->genomicend - hit3->genomicend) + querylength5 + querylength3;
    } else {
      new->insertlength = 0;
    }

  } else if (hit5->plusp == true) {
    /* Concordant directions on same chromosome (plus) */
    debug10(printf("Concordant on plus strand\n"));
    new->dir = +1;

    if (expect_concordant_p == true) {
      /* Try to resolve ambiguity on inside of concordant ends */
      resolve_inside_ambiguous_splice_plus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					   splicesites,query5_compress_fwd,query3_compress_fwd,
					   localsplicing_penalty,querylength5,querylength3,genestrand);
    }

    /* Have 5-start..end and 3-start..end */
    if (hit5->genomicend < hit3->genomicstart) {
      /* No overlap */
      new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      debug10(printf("plus, no overlap: insert length %d = start3 %u - end5 %u + %d + %d\n",
		     new->insertlength,hit3->genomicstart,hit5->genomicend,querylength5,querylength3));
#if 0
    } else if (hit5->genomicend > hit3->genomicend + SUBSUMPTION_SLOP) {
      /* hit5 subsumes hit3 */
      debug10(printf("plus, subsumption %u > %u\n",hit5->genomicend,hit3->genomicend));
      new->insertlength = 0;
#endif
    } else {
      new->insertlength = pair_insert_length(hit5,hit3);
    }


  } else {
    /* Concordant directions on same chromosome (minus) */
    debug10(printf("Concordant on minus strand\n"));
    new->dir = -1;

    if (expect_concordant_p == true) {
      /* Try to resolve ambiguity on inside of concordant ends */
      resolve_inside_ambiguous_splice_minus(&unresolved_amb_nmatches,&hit5,&hit3,&private5p,&private3p,
					    splicesites,query5_compress_rev,query3_compress_rev,
					    localsplicing_penalty,querylength5,querylength3,genestrand);
    }

    /* Have 3-end..start and 5-end..start */
    if (hit3->genomicstart < hit5->genomicend) {
      /* No overlap */
      new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      debug10(printf("minus, no overlap: insert length %d = end5 %u - start3 %u + %d + %d\n",
		     new->insertlength,hit5->genomicend,hit3->genomicstart,querylength5,querylength3));
#if 0
    } else if (hit3->genomicstart > hit5->genomicstart + SUBSUMPTION_SLOP) {
      /* hit3 subsumes hit5 */
      debug10(printf("minus, subsumption %u > %u\n",hit3->genomicstart,hit5->genomicstart));
      new->insertlength = 0;
#endif
    } else {
      new->insertlength = pair_insert_length(hit5,hit3);
    }

  }


  debug5(printf("\nGot insertlength of %d\n",new->insertlength));
  if (new->insertlength <= 0) {
    /* Not concordant */
#ifdef USE_BINGO
    new->absdifflength_bingo_p = false;
#endif
    new->absdifflength = -1U;

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

  } else {
    if (new->insertlength < expected_pairlength) {
      new->absdifflength = expected_pairlength - new->insertlength;
    } else {
      new->absdifflength = new->insertlength - expected_pairlength;
    }
#ifdef USE_BINGO
    if (new->absdifflength <= pairlength_deviation) {
      new->absdifflength_bingo_p = true;
    } else {
      new->absdifflength_bingo_p = false;
    }
#endif
  }

  if (SENSE_CONSISTENT_P(hit5->sensedir,hit3->sensedir)) {
    new->sense_consistent_p = true;
  } else {
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
    fprintf(stderr,"new->low %u > new->high %u, hit5->chrnum %d\n",new->low,new->high,hit5->chrnum);
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

  debug0(printf("Created new pair from %p and %p with private %d, %d\n",hit5,hit3,private5p,private3p));
  debug0(printf("  hittypes %d and %d\n",hit5->hittype,hit3->hittype));
  debug0(printf("  sensedirs %d and %d\n",hit5->sensedir,hit3->sensedir));

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
  } else if (x->low < y->low) {
    return -1;
  } else if (y->low < x->low) {
    return +1;
  } else if (x->high < y->high) {
    return +1;
  } else if (y->high < x->high) {
    return -1;

  } else if (x->hit5->low < y->hit5->low) {
    return -1;
  } else if (y->hit5->low < x->hit5->low) {
    return +1;
  } else if (x->hit5->high < y->hit5->high) {
    return +1;
  } else if (y->hit5->high < x->hit5->high) {
    return -1;

  } else if (x->hit3->low < y->hit3->low) {
    return -1;
  } else if (y->hit3->low < x->hit3->low) {
    return +1;
  } else if (x->hit3->high < y->hit3->high) {
    return +1;
  } else if (y->hit3->high < x->hit3->high) {
    return -1;

#if 0
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
	  printf("  Final %d: %u-%u (dir = %d), absdistance %u\n",
		 i,hitpair->low,hitpair->high,hitpair->dir,hitpair->absdifflength);
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
	   printf("  Initial %d (%s): %u..%u  %u..%u|%u..%u (dir = %d), nmatches: %d\n",
		  i,pairtype_string(hitpair->pairtype),hitpair->low,hitpair->high,
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
hitpair_goodness_cmp (Stage3pair_T hitpair,
#ifdef DEBUG8
		      int k,
#endif
		      Stage3pair_T best_hitpair) {

#if 0
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

  if ((hitpair->hit5->hittype == TERMINAL || hitpair->hit3->hittype == TERMINAL) &&
      best_hitpair->hit5->hittype != TERMINAL && best_hitpair->hit3->hittype != TERMINAL) {
    /* k is worse */
    debug8(printf(" => %d loses by terminal\n",k));
    return -1;
  } else if (hitpair->hit5->hittype != TERMINAL && hitpair->hit3->hittype != TERMINAL && 
	     (best_hitpair->hit5->hittype == TERMINAL || best_hitpair->hit3->hittype == TERMINAL)) {
    /* k is better */
    debug8(printf(" => %d wins by terminal\n",k));
    return +1;

  } else if (hitpair->score > best_hitpair->score) {
    /* k is worse */
    debug8(printf(" => %d loses by score\n",k));
    return -1;
  } else if (hitpair->score < best_hitpair->score) {
    /* k is better */
    debug8(printf(" => %d wins by score\n",k));
    return +1;

  } else if (hitpair->nmatches < best_hitpair->nmatches) {
    /* k is worse */
    debug8(printf(" => %d loses by nmatches\n",k));
    return -1;
  } else if (hitpair->nmatches > best_hitpair->nmatches) {
    /* k is better */
    debug8(printf(" => %d wins by nmatches\n",k));
    return +1;

  } else if (hitpair->nmatches_posttrim < best_hitpair->nmatches_posttrim) {
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

  } else if (hitpair->outerlength > best_hitpair->outerlength) {
    /* k is worse */
    debug8(printf(" => %d loses by outerlength\n",k));
    return -1;
  } else if (hitpair->outerlength < best_hitpair->outerlength) {
    /* k is better */
    debug8(printf(" => %d wins by outerlength\n",k));
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

    /* If insert length is within deviation of expected pairlength, favor it */
  } else if (best_hitpair->absdifflength <= (Genomicpos_T) pairlength_deviation &&
	     hitpair->absdifflength > (Genomicpos_T) pairlength_deviation) {
    /* k is worse */
    debug8(printf(" => %d loses by absdifflength within deviation %d\n",k,pairlength_deviation));
    return -1;
  } else if (hitpair->absdifflength <= (Genomicpos_T) pairlength_deviation &&
	     best_hitpair->absdifflength > (Genomicpos_T) pairlength_deviation) {
    /* k is better */
    debug8(printf(" => %d wins by absdifflength within deviation %d\n",k,pairlength_deviation));
    return +1;

    /* Otherwise, favor longer insert lengths to give more compact splices */
  } else if (hitpair->insertlength < best_hitpair->insertlength) {
    /* k is worse */
    debug8(printf(" => %d loses by insertlength\n",k));
    return -1;
  } else if (hitpair->insertlength > best_hitpair->insertlength) {
    /* k is better */
    debug8(printf(" => %d wins by insertlength\n",k));
    return +1;

  } else {
    return 0;
  }
}


static bool
hitpair_bad_superstretch_p (Stage3pair_T hitpair_k, Stage3pair_T *hitpairs, int k, int j) {
  int a;

  for (a = k+1; a <= j; a++) {
    if (hitpair_subsumption(hitpair_k,hitpairs[a]) == true) {
      debug8(printf("Testing %d because stretches over %d",k,a));
      if (hitpair_goodness_cmp(hitpairs[a],
#ifdef DEBUG8
			       a,
#endif
			       hitpair_k) >= 0) {
	debug8(printf(" => eliminating\n"));
	return true;
      }
      debug8(printf("\n"));
    }
  }
  return false;
}


static List_T
pair_remove_overlaps (List_T hitpairlist) {
  List_T unique = NULL;
  Stage3pair_T best_hitpair, hitpair, *hitpairs, *prev;
  int cmp;
  int nkept, n, i, j, k, besti;
  bool *eliminate;
#ifdef PRE_RESOLVE_MULTIMAPPING
  long int best_tally;
#endif

  n = List_length(hitpairlist);
  debug8(printf("Entered Stage3pair_remove_overlaps with %d pairs\n",n));
  if (n == 0) {
    return NULL;
  } else {
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    List_free(&hitpairlist);
  }

  /* Step 1.  Check for exact duplicates */
  debug8(printf("Checking for exact duplicates\n"));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s, %s-%s): %u..%u  %u..%u|%u..%u (dir = %d), nmatches: %d (%d posttrim), indel_low %d and %d\n",
		  i,pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
		  hittype_string(hitpair->hit3->hittype),hitpair->low,hitpair->high,
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
      debug8(printf("  Keeping %u..%u  %u..%u|%u..%u, nmatches (trimmed) %d (dir = %d)\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches,hitpair->dir));
      hitpairs[j++] = hitpair;
    } else {
      debug8(printf("  Eliminating %u..%u  %u..%u|%u..%u, nmatches (trimmed) %d (dir = %d)\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches,hitpair->dir));
      Stage3pair_free(&hitpair);
    }
  }

  FREE(prev);


  /* Step 2: Check for superstretches */
  n = nkept;
  debug8(printf("Checking for superstretches among %d hitpairs within subsumption clusters\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s, %s-%s): %u..%u  %u..%u|%u..%u (dir = %d), score: %d, nmatches: %d (%d posttrim), absdifflength: %d, nnovel: %d, nknown: %d, insertlength: %u, outerlength: %u\n",
		  i,pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
		  hittype_string(hitpair->hit3->hittype),hitpair->low,hitpair->high,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hitpair->dir,hitpair->score,hitpair->nmatches,hitpair->nmatches_posttrim,
		  hitpair->absdifflength,hitpair->nchimera_novel,hitpair->nchimera_known,
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
	if (hitpair_bad_superstretch_p(hitpairs[k],hitpairs,k,j) == true) {
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
      debug8(printf("  Keeping %u..%u  %u..%u|%u..%u, nmatches (trimmed) %d (dir = %d)\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches,hitpair->dir));
      hitpairs[j++] = hitpair;
    } else {
      debug8(printf("  Eliminating %u..%u  %u..%u|%u..%u, nmatches (trimmed) %d (dir = %d)\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches,hitpair->dir));
      Stage3pair_free(&hitpair);
    }
  }

  FREE(prev);


  /* Step 3: Check for best within subsumption clusters */
  n = nkept;
  debug8(printf("Checking for best among %d hitpairs within subsumption clusters\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  /* qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp); -- No need since original order was kept */
  
  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s, %s-%s): %u..%u  %u..%u|%u..%u (dir = %d), score: %d, nmatches: %d (%d posttrim), absdifflength: %d, nnovel: %d, nknown: %d, insertlength: %u, outerlength: %u\n",
		  i,pairtype_string(hitpair->pairtype),hittype_string(hitpair->hit5->hittype),
		  hittype_string(hitpair->hit3->hittype),hitpair->low,hitpair->high,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hitpair->dir,hitpair->score,hitpair->nmatches,hitpair->nmatches_posttrim,
		  hitpair->absdifflength,hitpair->nchimera_novel,hitpair->nchimera_known,
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
	cmp = hitpair_goodness_cmp(hitpairs[k],
#ifdef DEBUG8
				   k,
#endif
				   best_hitpair);
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
	} else if (hitpair_goodness_cmp(hitpairs[k],
#ifdef DEBUG8
					k,
#endif
					best_hitpair) <= 0) {
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
	cmp = hitpair_goodness_cmp(hitpairs[k],
#ifdef DEBUG8
				   k,
#endif
				   best_hitpair);
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
	} else if (hitpair_goodness_cmp(hitpairs[k],
#ifdef DEBUG8
					k,
#endif
					best_hitpair) <= 0) {
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
      debug8(printf("  Keeping %u..%u  %u..%u|%u..%u, nmatches (trimmed) %d (dir = %d)\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches,hitpair->dir));
      hitpairs[j++] = hitpair;
    } else {
      debug8(printf("  Eliminating %u..%u  %u..%u|%u..%u, nmatches (trimmed) %d (dir = %d)\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->nmatches,hitpair->dir));
      Stage3pair_free(&hitpair);
    }
  }

  FREE(prev);


  /* Step 4: Check for identity */
  n = nkept;
  debug8(printf("Checking for duplicates among %d hitpairs by identity\n",n));

  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
  /* qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_sort_cmp); -- No need since original order was kept */

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s): %u..%u  %u..%u|%u..%u (dir = %d), score: %d, nmatches: %d (%d posttrim), absdifflength: %d, insertlength: %u, outerlength: %u\n",
		  i,pairtype_string(hitpair->pairtype),hitpair->low,hitpair->high,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hitpair->dir,hitpair->score,hitpair->nmatches,hitpair->nmatches_posttrim,
		  hitpair->absdifflength,hitpair->insertlength,hitpair->outerlength);
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
    debug8(printf("Looking at %d with score %d, absdifflength %d, and outerlength %u\n",
		  i,hitpairs[i]->score,hitpairs[i]->absdifflength,hitpairs[i]->outerlength));
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

  debug8(printf("Exited Stage3pair_remove_overlaps with %d pairs\n",List_length(unique)));
  return unique;
}


List_T
Stage3pair_remove_overlaps (List_T hitpairlist) {
  List_T indep_overlapping = NULL, unique_separate, unique_overlapping,
    separate = NULL, overlapping = NULL, p;
  Stage3pair_T *array_separate, *array_overlapping;
  Stage3pair_T hitpair, hitpair_overlapping;
  Genomicpos_T low, high;
  bool subsumedp;
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

  unique_separate = pair_remove_overlaps(separate);
  unique_overlapping = pair_remove_overlaps(overlapping);
  
  if (unique_overlapping == NULL) {
    return unique_separate;
  } else if (unique_separate == NULL) {
    return unique_overlapping;
  } else {
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
	subsumedp = hitpair_subsumption(array_separate[j],hitpair_overlapping);
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
Stage3pair_eval_and_sort (int *npaths, int *second_absmq,
			  Stage3pair_T *stage3pairarray, int maxpaths,
			  Shortread_T queryseq5, Shortread_T queryseq3,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev, 
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev, 
			  Genome_T genome, char *quality_string_5, char *quality_string_3) {
  char *query5, *query3;
  double maxlik, loglik;

  double total, q;
  int mapq_score;

  int compute_npaths;
  int i;

  if (*npaths == 0) {
    /* Skip */

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

    *second_absmq = 0;

  } else {
    /* Compute mapq_loglik */
    for (i = 0; i < *npaths; i++) {
      stage3pairarray[i]->mapq_loglik =
	Stage3end_compute_mapq(stage3pairarray[i]->hit5,query5_compress_fwd,query5_compress_rev,
			       quality_string_5);
      stage3pairarray[i]->mapq_loglik +=
      Stage3end_compute_mapq(stage3pairarray[i]->hit3,query3_compress_fwd,query3_compress_rev,
			     quality_string_3);
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
    *second_absmq = stage3pairarray[1]->absmq_score;


    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < *npaths; i++) {
      total += (stage3pairarray[i]->mapq_loglik = exp(stage3pairarray[i]->mapq_loglik));
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
List_T
Stage3pair_optimal_score (List_T hitpairlist, int cutoff_level, int suboptimal_mismatches, bool keep_gmap_p) {
  List_T optimal = NULL, p;
  Stage3pair_T hitpair;
  int n;
  int minscore = MAX_READLENGTH + MAX_READLENGTH;
  int max_nmatches = 0, max_nmatches_posttrim;
#ifdef USE_OPTIMAL_SCORE_BINGO
  int minscore_bingo = MAX_READLENGTH + MAX_READLENGTH;
#endif
  bool non_translocation_p = false;
  bool non_terminal5_p = false, non_terminal3_p = false;

  n = List_length(hitpairlist);
  if (n == 0) {
    return NULL;
  }

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->hit5->chrnum != 0 && hitpair->hit3->chrnum != 0) {
      non_translocation_p = true;
    }
    if (hitpair->hit5->hittype != TERMINAL) {
      non_terminal5_p = true;
    }
    if (hitpair->hit3->hittype != TERMINAL) {
      non_terminal3_p = true;
    }
  }

  debug6(printf("non_terminal5_p = %d\n",non_terminal5_p));
  debug6(printf("non_terminal3_p = %d\n",non_terminal3_p));

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    debug6(printf("%u..%u  %u..%u|%u..%u types %s and %s, score %d (%d+%d), pairlength %d, outerlength %u\n",
		  hitpair->low,hitpair->high,
		  hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		  hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		  hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),
		  hitpair->score,hitpair->hit5->score,hitpair->hit3->score,
		  hitpair->insertlength,hitpair->outerlength));
    if ((hitpair->hit5->chrnum == 0 || hitpair->hit3->chrnum == 0) && non_translocation_p == true) {
      /* Skip, since we will eliminate anyway */
    } else {
      if (hitpair->nmatches > max_nmatches) {
	max_nmatches = hitpair->nmatches;
	max_nmatches_posttrim = hitpair->nmatches_posttrim;
      }
#ifdef TERMINAL_SECOND_CLASS
      if (non_terminal5_p == true && hitpair->hit5->hittype == TERMINAL) {
	/* Skip from setting minscore */
      } else if (non_terminal3_p == true && hitpair->hit3->hittype == TERMINAL) {
	/* Skip from setting minscore */
      }
#endif
      if (hitpair->score <= cutoff_level) {
	if (hitpair->score < minscore) {
	  minscore = hitpair->score;
	}
#ifdef USE_OPTIMAL_SCORE_BINGO
	if (hitpair->absdifflength_bingo_p == true) {
	  if (hitpair->score < minscore_bingo) {
	    minscore_bingo = hitpair->score;
	  }
	}
#endif
      }
    }
  }

#ifdef USE_OPTIMAL_SCORE_BINGO
  debug6(printf("Stage3pair_optimal_score over %d pairs: minscore = %d + subopt:%d, minscore_bingo = %d\n",
		n,minscore,suboptimal_mismatches,minscore_bingo));
#else
  debug6(printf("Stage3pair_optimal_score over %d pairs: minscore = %d + subopt:%d, max_nmatches = %d + subopt:%d\n",
		n,minscore,suboptimal_mismatches,max_nmatches,suboptimal_mismatches));
#endif
  minscore += suboptimal_mismatches;
  max_nmatches -= suboptimal_mismatches;
  max_nmatches_posttrim -= suboptimal_mismatches;

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

    if (keep_gmap_p == true && (hitpair->hit5->hittype == GMAP || hitpair->hit3->hittype == GMAP)) {
      /* GMAP hits already found to be better than their corresponding terminals */
      debug6(printf("Keeping a hit pair of type GMAP with score %d\n",hitpair->score));
      optimal = List_push(optimal,hitpair);

    } else if ((hitpair->hit5->chrnum == 0 || hitpair->hit3->chrnum == 0) && non_translocation_p == true) {
      debug6(printf("Eliminating a hit pair at %u..%u  %u..%u|%u..%u with splice translocation\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset));
      Stage3pair_free(&hitpair);

#ifdef TERMINAL_SECOND_CLASS
    } else if (hitpair->hit5->hittype == TERMINAL && non_terminal5_p == true) {
      if (hitpair->nmatches >= max_nmatches) {
	debug6(printf("Keeping a 5' terminal pair with nmatches %d\n",hitpair->nmatches));
	optimal = List_push(optimal,hitpair);
      } else {
	debug6(printf("Eliminating a hit pair at %u..%u  %u..%u|%u..%u with terminal on 5'\n",
		      hitpair->low,hitpair->high,
		      hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset));
	Stage3pair_free(&hitpair);
      }

    } else if (hitpair->hit3->hittype == TERMINAL && non_terminal3_p == true) {
      if (hitpair->nmatches >= max_nmatches) {
	debug6(printf("Keeping a 3' terminal pair with nmatches %d\n",hitpair->nmatches));
	optimal = List_push(optimal,hitpair);
      } else {
	debug6(printf("Eliminating a hit pair at %u..%u  %u..%u|%u..%u with terminal on 3'\n",
		      hitpair->low,hitpair->high,
		      hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		      hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset));
	Stage3pair_free(&hitpair);
      }
#endif

    } else if (hitpair->score > cutoff_level) {
      debug6(printf("Eliminating a hit pair at %u..%u  %u..%u|%u..%u with score %d > cutoff_level %d\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->score,cutoff_level));
      Stage3pair_free(&hitpair);

#ifdef USE_OPTIMAL_SCORE_BINGO
    } else if (hitpair->score >= minscore_bingo && hitpair->absdifflength_bingo_p == false) {
      debug6(printf("Eliminating a non-bingo hit pair at %u..%u  %u..%u|%u..%u with score %d and pairlength %u\n",
		    hitpair->low,hitpair->high,
		    hitpair->hit5->low - hitpair->hit5->chroffset,hitpair->hit5->high - hitpair->hit5->chroffset,
		    hitpair->hit3->low - hitpair->hit3->chroffset,hitpair->hit3->high - hitpair->hit3->chroffset,
		    hitpair->score,hitpair->insertlength));
      Stage3pair_free(&hitpair);
#endif

    } else if (hitpair->score > minscore && hitpair->nmatches_posttrim < max_nmatches_posttrim) {
      debug6(printf("Eliminating a hit pair with score %d\n",hitpair->score));
      Stage3pair_free(&hitpair);

    } else {
      debug6(printf("Keeping a hit pair with score %d\n",hitpair->score));
      optimal = List_push(optimal,hitpair);
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
	Stage3pair_free(&hitpair);
      } else if (hitpair->outerlength > best_outerlength) {
	debug6(printf("Eliminating a hit pair with outerlength %u\n",hitpair->outerlength));
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



/* Finds concordant pairs if nconcordant is 0 */
List_T
Stage3_pair_up_concordant (bool *abort_pairing_p, int *found_score, int *nconcordant,
			   List_T *samechr, List_T *conc_transloc, List_T hitpairs,
			   List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			   int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			   Genomicpos_T *splicesites,
			   Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			   Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			   int querylength5, int querylength3, int maxpairedpaths,
			   int splicing_penalty, int genestrand) {
#if 0
  /* No need for copies anymore, because effective_chrnum determined by inner substrings */
  List_T copies = NULL;
#endif
  int new_found_score = *found_score;
  List_T q, prev_start;
  Stage3pair_T stage3pair;
  T **hits5_plus, **hits5_minus, **hits3_plus, **hits3_minus, *hits5, *hits3, hit5, hit3;
  int *nhits5_plus, *nhits5_minus, *nhits3_plus, *nhits3_minus, nhits5, nhits3;
  int pairscore, score5, score3, i, j;
  bool *sorted5p, *sorted3p;
  Genomicpos_T insert_start;
  
  debug5(printf("Starting Stage3_pair_up_concordant with %d concordant, narray5 %d, narray3 %d, found_score %d\n",
		*nconcordant,narray5,narray3,*found_score));

  /* Find sizes for allocating memory */
  debug5(printf("Sizes of 5-end pieces by score level:\n"));
  nhits5_plus = (int *) CALLOC(cutoff_level_5+1,sizeof(int));
  nhits5_minus = (int *) CALLOC(cutoff_level_5+1,sizeof(int));
  for (i = 0; i < narray5; i++) {
    debug5(printf("  array5 score level %d with %d hits\n",i,List_length(hitarray5[i])));
    for (q = hitarray5[i]; q != NULL; q = q->rest) {
      hit5 = (T) q->first;
      debug5(printf("5': %p score %d, type %s\n",hit5,hit5->score,hittype_string(hit5->hittype)));
      assert(hit5->score >= 0);
      if (hit5->score > cutoff_level_5) {
	debug5(printf("Skipping hit with score %d > cutoff level %d\n",hit5->score,cutoff_level_5));
#if 0
      } else if (hit5->chimera_ambiguous_p == true) {
	/* Allow insert length to determine correct chimera from among the ambiguous */
	debug5(printf("Skipping hit with score %d because chimera is ambiguous\n",hit5->score));

      } else if (hit5->chrnum == 0) {
	if (allow_concordant_translocations_p == true) {
	  /* Translocation: Enter under each substring chrnum */
	  if (Substring_plusp(hit5->substring1) == true) {
	    nhits5_plus[hit5->score]++;
	  } else {
	    nhits5_minus[hit5->score]++;
	  }

	  if (Substring_plusp(hit5->substring2) == true) {
	    nhits5_plus[hit5->score]++;
	  } else {
	    nhits5_minus[hit5->score]++;
	  }
	}
#endif

      } else if (hit5->plusp == true) {
	nhits5_plus[hit5->score]++;
      } else {
	nhits5_minus[hit5->score]++;
      }
    }
  }

  debug5(
	 printf("Sizes of 5-end pieces by score level and plus/minus:\n");
	 for (score5 = 0; score5 <= cutoff_level_5; score5++) {
	   printf("  score %d: %d plus, %d minus\n",score5,nhits5_plus[score5],nhits5_minus[score5]);
	 }
	 );


  /* Reset cutoff_level_5 */
  score5 = 0;
  nhits5 = nhits5_plus[score5] + nhits5_minus[score5];
  while (score5+1 <= cutoff_level_5 && nhits5 + nhits5_plus[score5+1] + nhits5_minus[score5+1] < MAX_HITS) {
    nhits5 += nhits5_plus[score5+1] + nhits5_minus[score5+1];
    score5++;
    debug5(printf("Allowing score5 to go to %d, because nhits = %d\n",score5,nhits5));
  }
  debug5(printf("Resetting cutoff_level_5 to be %d\n",score5));
  cutoff_level_5 = score5;


  debug5(printf("Sizes of 3-end pieces by score level:\n"));
  nhits3_plus = (int *) CALLOC(cutoff_level_3+1,sizeof(int));
  nhits3_minus = (int *) CALLOC(cutoff_level_3+1,sizeof(int));
  for (i = 0; i < narray3; i++) {
    debug5(printf(" array3 score level %d with %d hits\n",i,List_length(hitarray3[i])));
    for (q = hitarray3[i]; q != NULL; q = q->rest) {
      hit3 = (T) q->first;
      debug5(printf("3': %p score %d, type %s\n",hit3,hit3->score,hittype_string(hit3->hittype)));
      assert(hit3->score >= 0);
      if (hit3->score > cutoff_level_3) {
	debug5(printf("Skipping hit with score %d > cutoff level %d\n",hit3->score,cutoff_level_3));
#if 0
      } else if (hit3->chimera_ambiguous_p == true) {
	/* Allow insert length to determine correct chimera from among the ambiguous */
	debug5(printf("Skipping hit with score %d because chimera is ambiguous\n",hit3->score));

      } else if (hit3->chrnum == 0) {
	if (allow_concordant_translocations_p == true) {
	  /* Translocation: Enter under each substring chrnum */
	  if (Substring_plusp(hit3->substring1) == true) {
	    nhits3_plus[hit3->score]++;
	  } else {
	    nhits3_minus[hit3->score]++;
	  }

	  if (Substring_plusp(hit3->substring2) == true) {
	    nhits3_plus[hit3->score]++;
	  } else {
	    nhits3_minus[hit3->score]++;
	  }
	}
#endif

      } else if (hit3->plusp == true) {
	nhits3_plus[hit3->score]++;
      } else {
	nhits3_minus[hit3->score]++;
      }
    }
  }
  debug5(
	 printf("Sizes of 3-end pieces by score level and plus/minus:\n");
	 for (score3 = 0; score3 <= cutoff_level_3; score3++) {
	   printf("  score %d: %d plus, %d minus\n",score3,nhits3_plus[score3],nhits3_minus[score3]);
	 }
	 );


  /* Reset cutoff_level_3 */
  score3 = 0;
  nhits3 = nhits3_plus[score3] + nhits3_minus[score3];
  while (score3+1 <= cutoff_level_3 && nhits3 + nhits3_plus[score3+1] + nhits3_minus[score3+1] < MAX_HITS) {
    nhits3 += nhits3_plus[score3+1] + nhits3_minus[score3+1];
    score3++;
    debug5(printf("Allowing score3 to go to %d, because nhits = %d\n",score3,nhits3));
  }
  debug5(printf("Resetting cutoff_level_3 to be %d\n",score3));
  cutoff_level_3 = score3;



  /* Store hits5 */
  hits5_plus = (T **) CALLOC(cutoff_level_5+1,sizeof(T *));
  for (score5 = 0; score5 <= cutoff_level_5; score5++) {
    if (nhits5_plus[score5] == 0) {
      hits5_plus[score5] = (T *) NULL;
    } else {
      hits5_plus[score5] = (T *) CALLOC(nhits5_plus[score5],sizeof(Stage3end_T));
    }
  }

  hits5_minus = (T **) CALLOC(cutoff_level_5+1,sizeof(T *));
  for (score5 = 0; score5 <= cutoff_level_5; score5++) {
    if (nhits5_minus[score5] == 0) {
      hits5_minus[score5] = (T *) NULL;
    } else {
      hits5_minus[score5] = (T *) CALLOC(nhits5_minus[score5],sizeof(Stage3end_T));
    }
  }

  for (score5 = 0; score5 <= cutoff_level_5; score5++) {
    nhits5_plus[score5] = 0;
    nhits5_minus[score5] = 0;
  }

  for (i = 0; i < narray5; i++) {
    for (q = hitarray5[i]; q != NULL; q = q->rest) {
      hit5 = (T) q->first;
      if (hit5->score > cutoff_level_5) {
	/* Skip */
#if 0
      } else if (hit5->chimera_ambiguous_p == true) {
	/* Allow insert length to determine correct chimera from among the ambiguous */
	/* Skip */

      } else if (hit5->chrnum == 0) {
	if (allow_concordant_translocations_p == true) {
	  /* Translocation: Enter under each substring chrnum */
	  copy = Stage3end_copy(hit5);
	  copy->effective_chrnum = Substring_chrnum(hit5->substring1);
	  copy->other_chrnum = Substring_chrnum(hit5->substring2);
	  copy->genomicstart = Substring_genomicstart(hit5->substring1);
	  copy->genomicend = Substring_genomicend(hit5->substring1);
	  copy->plusp = Substring_plusp(hit5->substring1);
	  copies = List_push(copies,(void *) copy);
	  if (copy->plusp == true) {
	    hits5_plus[hit5->score][nhits5_plus[hit5->score]++] = copy;
	  } else {
	    hits5_minus[hit5->score][nhits5_minus[hit5->score]++] = copy;
	  }

	  copy = Stage3end_copy(hit5);
	  copy->effective_chrnum = Substring_chrnum(hit5->substring2);
	  copy->other_chrnum = Substring_chrnum(hit5->substring1);
	  copy->genomicstart = Substring_genomicstart(hit5->substring2);
	  copy->genomicend = Substring_genomicend(hit5->substring2);
	  copy->plusp = Substring_plusp(hit5->substring2);
	  copies = List_push(copies,(void *) copy);
	  if (copy->plusp == true) {
	    hits5_plus[hit5->score][nhits5_plus[hit5->score]++] = copy;
	  } else {
	    hits5_minus[hit5->score][nhits5_minus[hit5->score]++] = copy;
	  }
	}
#endif

      } else if (hit5->plusp == true) {
	hits5_plus[hit5->score][nhits5_plus[hit5->score]++] = hit5;
      } else {
	hits5_minus[hit5->score][nhits5_minus[hit5->score]++] = hit5;
      }
    }
  }


  /* Store hits3 */
  hits3_plus = (T **) CALLOC(cutoff_level_3+1,sizeof(T *));
  for (score3 = 0; score3 <= cutoff_level_3; score3++) {
    if (nhits3_plus[score3] == 0) {
      hits3_plus[score3] = (T *) NULL;
    } else {
      hits3_plus[score3] = (T *) CALLOC(nhits3_plus[score3],sizeof(Stage3end_T));
    }
  }

  hits3_minus = (T **) CALLOC(cutoff_level_3+1,sizeof(T *));
  for (score3 = 0; score3 <= cutoff_level_3; score3++) {
    if (nhits3_minus[score3] == 0) {
      hits3_minus[score3] = (T *) NULL;
    } else {
      hits3_minus[score3] = (T *) CALLOC(nhits3_minus[score3],sizeof(Stage3end_T));
    }
  }

  for (score3 = 0; score3 <= cutoff_level_3; score3++) {
    nhits3_plus[score3] = 0;
    nhits3_minus[score3] = 0;
  }

  for (i = 0; i < narray3; i++) {
    for (q = hitarray3[i]; q != NULL; q = q->rest) {
      hit3 = (T) q->first;
      if (hit3->score > cutoff_level_3) {
	/* Skip */
#if 0
      } else if (hit3->chimera_ambiguous_p == true) {
	/* Allow insert length to determine correct chimera from among the ambiguous */
	/* Skip */
      } else if (hit3->chrnum == 0) {
	if (allow_concordant_translocations_p == true) {
	  /* Translocation: Enter under each substring chrnum */
	  copy = Stage3end_copy(hit3);
	  copy->effective_chrnum = Substring_chrnum(hit3->substring1);
	  copy->other_chrnum = Substring_chrnum(hit3->substring2);
	  copy->genomicstart = Substring_genomicstart(hit3->substring1);
	  copy->genomicend = Substring_genomicend(hit3->substring1);
	  copy->plusp = Substring_plusp(hit3->substring1);
	  copies = List_push(copies,(void *) copy);
	  if (copy->plusp == true) {
	    hits3_plus[hit3->score][nhits3_plus[hit3->score]++] = copy;
	  } else {
	    hits3_minus[hit3->score][nhits3_minus[hit3->score]++] = copy;
	  }

	  copy = Stage3end_copy(hit3);
	  copy->effective_chrnum = Substring_chrnum(hit3->substring2);
	  copy->other_chrnum = Substring_chrnum(hit3->substring1);
	  copy->genomicstart = Substring_genomicstart(hit3->substring2);
	  copy->genomicend = Substring_genomicend(hit3->substring2);
	  copy->plusp = Substring_plusp(hit3->substring2);
	  copies = List_push(copies,(void *) copy);
	  if (copy->plusp == true) {
	    hits3_plus[hit3->score][nhits3_plus[hit3->score]++] = copy;
	  } else {
	    hits3_minus[hit3->score][nhits3_minus[hit3->score]++] = copy;
	  }
	}
#endif

      } else if (hit3->plusp == true) {
	hits3_plus[hit3->score][nhits3_plus[hit3->score]++] = hit3;
      } else {
	hits3_minus[hit3->score][nhits3_minus[hit3->score]++] = hit3;
      }
    }
  }


  /* Look for concordant pairs */
  sorted5p = (bool *) CALLOC(cutoff_level_5+1,sizeof(bool));
  sorted3p = (bool *) CALLOC(cutoff_level_3+1,sizeof(bool));

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
	    debug5(printf("plus/plus: i=%d/%d %u..%u %s %s %p\n",
			  i,nhits5,hit5->genomicstart,hit5->genomicend,
			  print_sense(hit5->sensedir),hittype_string(hit5->hittype),hit5));

	    while (j >= 0 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j--;
	    }
	    j++;		/* Finish backup */

	    while (j < nhits3 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j++;
	    }

	    while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s %p",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      hit3 = hits3[j];
		
	      /* Want only pairs not previously seen */
	      if (hit5->paired_seenp == false || hit3->paired_seenp == false) {
#if 0
		if (hit5->hittype == TERMINAL && hit3->hittype == TERMINAL) {
		  debug5(printf(" => both terminal"));
		}
#endif
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
						   /*pairtype*/TRANSLOCATION,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {
		    *conc_transloc = List_push(*conc_transloc,(void *) stage3pair);
		  }

		} else if (SENSE_INCONSISTENT_P(hit5->sensedir,hit3->sensedir)) {
		  debug5(printf(" => sense inconsistent"));
		} else if (hit3->genomicend < hit5->genomicstart) {
		  debug5(printf(" => scramble because end3 %u < start5 %u\n",hit3->genomicend,hit5->genomicstart));
		  if ((stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		  }
		} else {
		  debug5(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(hit5,hit3,splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/CONCORDANT,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {
		    hitpairs = List_push(hitpairs,(void *) stage3pair);

		    if (pairscore < new_found_score) {
		      new_found_score = pairscore;
		      debug5(printf(" => tentatively updating found_score to be %d",new_found_score));
		    }
		    if (++(*nconcordant) > maxpairedpaths) {
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
	    debug5(printf("minus/minus: i=%d/%d %u..%u %s %s %p\n",
			  i,nhits3,hit3->genomicstart,hit3->genomicend,
			  print_sense(hit3->sensedir),hittype_string(hit3->hittype),hit3));

	    while (j >= 0 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j--;
	    }
	    j++;			/* Finish backup */

	    while (j < nhits5 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j++;
	    }

	    while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s %p",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      hit5 = hits5[j];

	      /* Want only pairs not previously seen */
	      if (hit3->paired_seenp == false || hit5->paired_seenp == false) {
#if 0
		if (hit3->hittype == TERMINAL && hit5->hittype == TERMINAL) {
		  debug5(printf(" => both terminal"));
		}
#endif
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
						   /*pairtype*/TRANSLOCATION,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {
		    *conc_transloc = List_push(*conc_transloc,(void *) stage3pair);
		  }

		} else if (SENSE_INCONSISTENT_P(hit3->sensedir,hit5->sensedir)) {
		  debug5(printf(" => sense inconsistent"));
		} else if (hit5->genomicstart < hit3->genomicend) {
		  debug5(printf(" => scramble because start5 %u < end3 %u\n",hit5->genomicstart,hit3->genomicend));
		  if ((stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		  }
		} else {
		  debug5(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(hit5,hit3,splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/CONCORDANT,splicing_penalty,
						   /*private5p*/false,/*private3p*/false,/*expect_concordant_p*/true)) != NULL) {
		    hitpairs = List_push(hitpairs,(void *) stage3pair);

		    if (pairscore < new_found_score) {
		      new_found_score = pairscore;
		      debug5(printf(" => updating new_found_score to be %d",new_found_score));
		    }
		    if (++(*nconcordant) > maxpairedpaths) {
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
	    debug5(printf("plus/minus: i=%d/%d %u..%u %s %s %p\n",
			  i,nhits5,hit5->genomicstart,hit5->genomicend,
			  print_sense(hit5->sensedir),hittype_string(hit5->hittype),hit5));

	    while (j >= 0 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j--;
	    }
	    j++;		/* Finish backup */

	    while (j < nhits3 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      j++;
	    }

	    while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s %p",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype),hits3[j]));
	      hit3 = hits3[j];
		
	      /* Want only pairs not previously seen */
	      if (hit5->paired_seenp == false || hit3->paired_seenp == false) {
#if 0
		if (hit5->hittype == TERMINAL && hit3->hittype == TERMINAL) {
		  debug5(printf(" => both terminal"));
		}
#endif
		if (hit5->effective_chrnum != hit3->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 /* && hit5->other_chrnum != hit3->other_chrnum */) {
		  debug5(printf(" => double splice translocations"));
		} else if (SENSE_INCONSISTENT_FOR_INVERSION_P(hit5->sensedir,hit3->sensedir)) {
		  debug5(printf(" => sense inconsistent"));
#if 0
		} else if (hits3[j]->genomicstart + querylength3 <= insert_start) {
		  debug5(printf(" => scramble"));
		  if ((stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		  }
#endif
		} else {
		  debug5(printf(" => inversion effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_INVERSION,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
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
	    debug5(printf("minus/plus: i=%d/%d %u..%u %s %s %p\n",
			  i,nhits3,hit3->genomicstart,hit3->genomicend,
			  print_sense(hit3->sensedir),hittype_string(hit3->hittype),hit3));

	    while (j >= 0 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j--;
	    }
	    j++;			/* Finish backup */

	    while (j < nhits5 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s %p\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      j++;
	    }

	    while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s %p",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype),hits5[j]));
	      hit5 = hits5[j];

	      /* Want only pairs not previously seen */
	      if (hit3->paired_seenp == false || hit5->paired_seenp == false) {
#if 0
		if (hit3->hittype == TERMINAL && hit5->hittype == TERMINAL) {
		  debug5(printf(" => both terminal"));
		}
#endif
		if (hit3->effective_chrnum != hit5->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 /* && hit5->other_chrnum != hit3->other_chrnum */) {
		  debug5(printf(" => double splice translocations"));
		} else if (SENSE_INCONSISTENT_FOR_INVERSION_P(hit3->sensedir,hit5->sensedir)) {
		  debug5(printf(" => sense inconsistent"));
#if 0
		} else if (hits5[j]->genomicend + querylength5 <= insert_start) {
		  debug5(printf(" => scramble"));
		  if ((stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
		  }
#endif
		} else {
		  debug5(printf(" => inversion effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  if ((stage3pair = Stage3pair_new(Stage3end_copy(hit5),Stage3end_copy(hit3),splicesites,
						   query5_compress_fwd,query5_compress_rev,
						   query3_compress_fwd,query3_compress_rev,genestrand,
						   /*pairtype*/PAIRED_INVERSION,splicing_penalty,
						   /*private5p*/true,/*private3p*/true,/*expect_concordant_p*/false)) != NULL) {
		    *samechr = List_push(*samechr,(void *) stage3pair);
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

    /* Mark all concordant pairs found at this pairscore level */
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


  FREE(sorted3p);
  FREE(sorted5p);

#if 0
  for (q = copies; q != NULL; q = List_next(q)) {
    copy = (T) List_head(q);
    Stage3end_free(&copy);
  }
  List_free(&copies);
#endif

  for (score3 = 0; score3 <= cutoff_level_3; score3++) {
    FREE(hits3_plus[score3]);
    FREE(hits3_minus[score3]);
  }
  FREE(hits3_plus);
  FREE(hits3_minus);
  FREE(nhits3_plus);
  FREE(nhits3_minus);

  for (score5 = 0; score5 <= cutoff_level_5; score5++) {
    FREE(hits5_plus[score5]);
    FREE(hits5_minus[score5]);
  }
  FREE(hits5_plus);
  FREE(hits5_minus);
  FREE(nhits5_plus);
  FREE(nhits5_minus);

  debug5(printf("Finished with Stage3_pair_up_concordant: %d concordant, %d samechr, %d conc_transloc\n",
		List_length(hitpairs),List_length(*samechr),List_length(*conc_transloc)));

  return hitpairs;
}


