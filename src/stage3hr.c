static char rcsid[] = "$Id: stage3hr.c 37257 2011-03-28 16:39:14Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage3hr.h"
#include "stage1hr.h"		/* To get MAX_QUERYLENGTH */

#include <stdlib.h>		/* For qsort */
#include <string.h>
#include <strings.h>
#include <ctype.h>		/* For islower */
#include <math.h>		/* For exp() and log10() */
#include "assert.h"
#include "mem.h"
#include "chrnum.h"
#include "complement.h"
#include "interval.h"
#include "listdef.h"
#include "substring.h"
#include "mapq.h"


#define MAX_HITS 100000


#define CONCORDANT_TEXT "concordant"
#define PAIRED_TEXT "paired"
#define UNPAIRED_TEXT "unpaired"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

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

/* insert length calculation */
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

/* Stage3pair_optimal_score */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif


/* Stage3_remove_duplicates */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Stage3pair_remove_duplicates */
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



static bool invert_first_p;
static bool invert_second_p;


void
Stage3hr_setup (bool invert_first_p_in, bool invert_second_p_in) {
  invert_first_p = invert_first_p_in;
  invert_second_p = invert_second_p_in;
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



/* Note: Substring_T has genomiclength, but not Stage3_T */
#define T Stage3_T
struct T {
  Hittype_T hittype;
  Chrnum_T chrnum; /* Needed for printing paired-end results.  A chrnum of 0 indicates a distant splice. */
  Chrnum_T effective_chrnum;	/* For determining concordance */
  Chrnum_T other_chrnum;	/* 0 for non-translocations, and other chrnum besides effective_chrnum for translocations */
  Genomicpos_T chroffset;

  int querylength_adj;		/* Adjusted for indels */

  Genomicpos_T genomicstart;
  Genomicpos_T genomicend;
  bool plusp;

  Genomicpos_T low;
  Genomicpos_T high;

  double mapq_loglik;
  int mapq_score;

  int score;			/* Includes colordiffs and penalties */
  int ntscore;			/* Includes penalties */
  int nmatches;

  int nmismatches_whole;
  int nmismatches_bothdiff;
  int nmismatches_refdiff;	/* Set only for display */

  int nindels;			/* for indels */
  int indel_pos;		/* for indels */
  char *deletion;		/* for deletions */

  Genomicpos_T distance;	/* for splicing or shortexon (sum of two distances) */
  Genomicpos_T acceptor_distance; /* for shortexon */
  Genomicpos_T donor_distance;	  /* for shortexon */
  int sensedir;			/* for splicing */

  bool start_ambiguous_p;
  bool end_ambiguous_p;

  int *start_ambi;
  int *end_ambi;
  int start_nambi;
  int end_nambi;
  int *start_amb_nmismatches;
  int *end_amb_nmismatches;

  /* Single: substring1 */
  /* Indel: substring1 + substring 2*/
  /* Halfsplice: substring1 */
  /* Splice: substring1 (donor) + substring2 (acceptor) */
  /* Shortexon: substring1 (shortexon) + substringD + substringA */

  Substring_T substring1;
  Substring_T substring2;
  Substring_T substringD;
  Substring_T substringA;

  Substring_T substring_low;	/* For SAM output */

  bool paired_seenp;   /* for paired-end.  set to true by Stage3_pair_up(). */
  bool concordantp;    /* for paired-end.  set to true by Stage3_pair_up(). */
};


static char *
pairtype_string (Pairtype_T pairtype) {
  switch (pairtype) {
  case CONCORDANT: return "concordant";
  case PAIRED_INVERSION: return "inversion";
  case PAIRED_SCRAMBLE: return "scramble";
  case PAIRED_TOOLONG: return "toolong";
  case UNPAIRED: return "unpaired";
  default: abort();
  }
}


struct Stage3pair_T {
  Pairtype_T pairtype;

  T hit5;
  T hit3;

  Genomicpos_T low;
  Genomicpos_T high;
  int insertlength;
  Genomicpos_T outerlength;

  double mapq_loglik;
  int mapq_score;

  int score;
  int nmatches;
  Genomicpos_T absdifflength;
  bool absdifflength_bingo_p;
  int dir;			/* -1, 0, or +1 */
  bool sense_consistent_p;
};



Hittype_T
Stage3_hittype (T this) {
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
  case ONE_THIRD_SHORTEXON: return "one-third-shortexon";
  case TWO_THIRDS_SHORTEXON: return "two-thirds-shortexon";
  case SHORTEXON: return "shortexon";
  case TERMINAL: return "terminal";
  default: abort();
  }
}


Chrnum_T
Stage3_chrnum (T this) {
  return this->chrnum;
}

Genomicpos_T
Stage3_chroffset (T this) {
  return this->chroffset;
}

Genomicpos_T
Stage3_genomicstart (T this) {
  return this->genomicstart;
}

Genomicpos_T
Stage3_genomicend (T this) {
  return this->genomicend;
}

/* For Goby */
int
Stage3_query_alignment_length (T this) {
  int length;

  length = Substring_match_length(this->substring1);
  length += Substring_match_length(this->substring2);
  length += Substring_match_length(this->substringD);
  length += Substring_match_length(this->substringA);
  if (this->hittype == INSERTION) {
    length += this->nindels;
  }
  return length;
}

/* For Goby */
Genomicpos_T
Stage3_genomic_alignment_length (T this) {
  Genomicpos_T length;

  length = Substring_genomic_alignment_length(this->substring1);
  length += Substring_genomic_alignment_length(this->substring2);
  length += Substring_genomic_alignment_length(this->substringD);
  length += Substring_genomic_alignment_length(this->substringA);
  if (this->hittype == DELETION) {
    length += (Genomicpos_T) this->nindels;
  }
  return length;
}


Genomicpos_T
Stage3_chrpos_low_trim (T this) {
  if (this->plusp == true) {
    return Substring_alignstart_trim(this->substring_low) - Substring_chroffset(this->substring_low);
  } else {
    return Substring_alignend_trim(this->substring_low) - Substring_chroffset(this->substring_low);
  }
}

int
Stage3_mapq_score (T this) {
  return this->mapq_score;
}

int
Stage3_score (T this) {
  return this->score;
}

int
Stage3_nmismatches_whole (T this) {
  return this->nmismatches_whole;
}

int
Stage3_nmismatches_bothdiff (T this) {
  return this->nmismatches_bothdiff;
}

int
Stage3_nmismatches_refdiff (T this) {
  return this->nmismatches_refdiff;
}

int
Stage3_nindels (T this) {
  return this->nindels;
}

int
Stage3_indel_pos (T this) {
  return this->indel_pos;
}


bool
Stage3_plusp (T this) {
  return this->plusp;
}

Substring_T
Stage3_substring1 (T this) {
  return this->substring1;
}

Substring_T
Stage3_substring2 (T this) {
  return this->substring2;
}

char *
Stage3_deletion_string (T this) {
  return this->deletion;
}

Substring_T
Stage3_substringD (T this) {
  if (this->hittype == SPLICE) {
    return this->substring1;
  } else if (this->hittype == SHORTEXON || this->hittype == ONE_THIRD_SHORTEXON || this->hittype == TWO_THIRDS_SHORTEXON) {
    return this->substringD;
  } else {
    fprintf(stderr,"Called Stage3_substringD with hittype %s\n",
	    hittype_string(this->hittype));
    abort();
  }
}

Substring_T
Stage3_substringA (T this) {
  if (this->hittype == SPLICE) {
    return this->substring2;
  } else if (this->hittype == SHORTEXON || this->hittype == ONE_THIRD_SHORTEXON || this->hittype == TWO_THIRDS_SHORTEXON) {
    return this->substringA;
  } else {
    fprintf(stderr,"Called Stage3_substringA with hittype %s\n",
	    hittype_string(this->hittype));
    abort();
  }
}

Substring_T
Stage3_substring_low (T this) {
  return this->substring_low;
}

Genomicpos_T
Stage3_distance (T this) {
  return this->distance;
}

Genomicpos_T
Stage3_shortexon_acceptor_distance (T this) {
  return this->acceptor_distance;
}

Genomicpos_T
Stage3_shortexon_donor_distance (T this) {
  return this->donor_distance;
}

int
Stage3_sensedir (T this) {
  return this->sensedir;
}



void
Stage3_free (T *old) {
  FREE((*old)->end_ambi);
  FREE((*old)->start_ambi);
  FREE((*old)->end_amb_nmismatches);
  FREE((*old)->start_amb_nmismatches);

  if ((*old)->deletion != NULL) {
    FREE((*old)->deletion);
  }
  if ((*old)->substring1 != NULL) {
    Substring_free(&(*old)->substring1);
  }
  if ((*old)->substring2 != NULL) {
    Substring_free(&(*old)->substring2);
  }
  if ((*old)->substringD != NULL) {
    Substring_free(&(*old)->substringD);
  }
  if ((*old)->substringA != NULL) {
    Substring_free(&(*old)->substringA);
  }

  FREE(*old);
  return;
}

Stage3_T
Stage3pair_hit5 (Stage3pair_T this) {
  return this->hit5;
}

Stage3_T
Stage3pair_hit3 (Stage3pair_T this) {
  return this->hit3;
}

int
Stage3pair_mapq_score (Stage3pair_T this) {
  return this->mapq_score;
}

Genomicpos_T
Stage3pair_pairlength (Stage3pair_T this) {
  return this->insertlength;
}

void
Stage3pair_free (Stage3pair_T *old) {
  if ((*old)->hit3 != NULL) {
    Stage3_free(&(*old)->hit3);
  }
  if ((*old)->hit5 != NULL) {
    Stage3_free(&(*old)->hit5);
  }
  FREE(*old);
  return;
}
  

static char complCode[128] = COMPLEMENT_LC;

static char *
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return complement;
}


static T
Stage3_copy (T old) {
  T new = (T) MALLOC(sizeof(*new));

  new->hittype = old->hittype;
  new->chrnum = old->chrnum;
  new->effective_chrnum = old->effective_chrnum;
  new->other_chrnum = old->other_chrnum;
  new->chroffset = old->chroffset;

  new->querylength_adj = old->querylength_adj;

  new->genomicstart = old->genomicstart;
  new->genomicend = old->genomicend;
  new->plusp = old->plusp;

  new->low = old->low;
  new->high = old->high;

  new->mapq_loglik = old->mapq_loglik;
  new->mapq_score = old->mapq_score;

  new->score = old->score;
  new->ntscore = old->ntscore;
  new->nmatches = old->nmatches;

  new->nmismatches_whole = old->nmismatches_whole;
  new->nmismatches_bothdiff = old->nmismatches_bothdiff;
  new->nmismatches_refdiff = old->nmismatches_refdiff;

  new->nindels = old->nindels;
  new->indel_pos = old->indel_pos;
  if (old->deletion == NULL) {
    new->deletion = (char *) NULL;
  } else {
    new->deletion = (char *) CALLOC(strlen(old->deletion)+1,sizeof(char));
    strcpy(new->deletion,old->deletion);
  }
  
  new->distance = old->distance;
  new->acceptor_distance = old->acceptor_distance;
  new->donor_distance = old->donor_distance;
  new->sensedir = old->sensedir;

  new->start_ambiguous_p = old->start_ambiguous_p;
  new->end_ambiguous_p = old->end_ambiguous_p;

  if ((new->start_nambi = old->start_nambi) == 0) {
    new->start_ambi = (int *) NULL;
    new->start_amb_nmismatches = (int *) NULL;
  } else {
    new->start_ambi = (int *) CALLOC(old->start_nambi,sizeof(int));
    memcpy(new->start_ambi,old->start_ambi,old->start_nambi*sizeof(int));
    new->start_amb_nmismatches = (int *) CALLOC(old->start_nambi,sizeof(int));
    memcpy(new->start_amb_nmismatches,old->start_amb_nmismatches,old->start_nambi*sizeof(int));
  }

  if ((new->end_nambi = old->end_nambi) == 0) {
    new->end_ambi = (int *) NULL;
    new->end_amb_nmismatches = (int *) NULL;
  } else {
    new->end_ambi = (int *) CALLOC(old->end_nambi,sizeof(int));
    memcpy(new->end_ambi,old->end_ambi,old->end_nambi*sizeof(int));
    new->end_amb_nmismatches = (int *) CALLOC(old->end_nambi,sizeof(int));
    memcpy(new->end_amb_nmismatches,old->end_amb_nmismatches,old->end_nambi*sizeof(int));
  }

  new->substring1 = Substring_copy(old->substring1);
  new->substring2 = Substring_copy(old->substring2);
  new->substringD = Substring_copy(old->substringD);
  new->substringA = Substring_copy(old->substringA);

  if (old->substring_low == NULL) {
    new->substring_low = NULL;
  } else if (old->substring_low == old->substring1) {
    new->substring_low = new->substring1;
  } else if (old->substring_low == old->substring2) {
    new->substring_low = new->substring2;
  } else if (old->substring_low == old->substringD) {
    new->substring_low = new->substringD;
  } else if (old->substring_low == old->substringA) {
    new->substring_low = new->substringA;
  } else {
    fprintf(stderr,"substring_low is not NULL, substring1, or substring2, or substringD, or substringA\n");
    abort();
  }

  new->paired_seenp = old->paired_seenp;
  new->concordantp = old->concordantp;

  return new;
}


T
Stage3_new_exact (int *found_score, Genomicpos_T left, int genomiclength,
		  Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		  bool plusp, Chrnum_T chrnum, Genomicpos_T chroffset) {
  T new;
  Substring_T substring;
  Genomicpos_T genomicstart, genomicend;

  debug0(printf("Stage3_new_exact: left %u, chrnum %d\n",left,chrnum));

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + genomiclength;
  } else {
    genomicstart = left + genomiclength;
    genomicend = left;
  }

  if ((substring = Substring_new(/*nmismatches*/0,/*ncolordiffs*/0,chrnum,
				 chroffset,left,genomicstart,genomicend,query_compress,genome_blocks,snp_blocks,
				 /*start_endtype*/END,/*end_endtype*/END,
				 /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				 /*alignstart*/genomicstart,/*alignend*/genomicend,
				 genomiclength,/*extraleft*/0,/*extraright*/0,
				 /*exactp*/true,/*query*/NULL,plusp,
				 /*trim_left_p*/false,/*trim_right_p*/false,
				 /*minlength*/0,/*dibasep*/false,/*cmetp*/false)) == NULL) {
    return (T) NULL;

  } else {
    new = (T) MALLOC(sizeof(*new));
    new->substring1 = substring;
    new->substring2 = (Substring_T) NULL;
    new->substringD = (Substring_T) NULL;
    new->substringA = (Substring_T) NULL;
    new->substring_low = new->substring1;

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
    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring);
    new->mapq_score = 0;
#endif

    new->nindels = 0;
    new->nmismatches_whole = 0;
    new->nmismatches_bothdiff = 0;
    /* new->nmismatches_refdiff = 0; */
    new->ntscore = 0;
    new->score = 0;
    new->nmatches = genomiclength;
    *found_score = 0;

    new->start_ambiguous_p = new->end_ambiguous_p = false;
    new->start_ambi = new->end_ambi = (int *) NULL;
    new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
    new->start_nambi = new->end_nambi = 0;
    new->distance = 0U;
    new->acceptor_distance = new->donor_distance = 0U;

    new->paired_seenp = false;
    new->concordantp = false;

    return new;
  }
}


T
Stage3_new_substitution (int *found_score, int nmismatches_whole, int ncolordiffs, Genomicpos_T left,
			 int genomiclength, Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			 bool plusp, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
			 bool dibasep, bool cmetp) {
  T new;
  Substring_T substring;
  Genomicpos_T genomicstart, genomicend;

  debug0(printf("Stage3_new_substitution: left %u, chrnum %d, nmismatches %d\n",
	       left,chrnum,nmismatches_whole));

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + genomiclength;
  } else {
    genomicstart = left + genomiclength;
    genomicend = left;
  }

  if ((substring = Substring_new(nmismatches_whole,ncolordiffs,chrnum,
				 chroffset,left,genomicstart,genomicend,query_compress,genome_blocks,snp_blocks,
				 /*start_endtype*/END,/*end_endtype*/END,
				 /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				 /*alignstart*/genomicstart,/*alignend*/genomicend,
				 genomiclength,/*extraleft*/0,/*extraright*/0,/*exactp*/false,query,plusp,
				 /*trim_left_p*/true,/*trim_right_p*/true,
				 /*minlength*/genomiclength/2,dibasep,cmetp)) == NULL) {
    return (T) NULL;

  } else {
    new = (T) MALLOC(sizeof(*new));
    new->substring1 = substring;
    new->substring2 = (Substring_T) NULL;
    new->substringD = (Substring_T) NULL;
    new->substringA = (Substring_T) NULL;
    new->substring_low = new->substring1;

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

    new->hittype = SUB;
    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring);
    new->mapq_score = 0;
#endif

    new->nindels = 0;
    new->nmismatches_whole = nmismatches_whole;
    new->ntscore = nmismatches_whole;
    new->score = nmismatches_whole + ncolordiffs;

    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(new->substring1);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1); */

#if 0
    /* This method was previously the only one correct for SNP-tolerant alignment */
    new->nmatches = Substring_match_length(new->substring1) - new->total_nmismatches;
#else
    /* This method is now correct for SNP-tolerant alignment */
    new->nmatches = Substring_nmatches(new->substring1);
#endif

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->start_ambiguous_p = new->end_ambiguous_p = false;
    new->start_ambi = new->end_ambi = (int *) NULL;
    new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
    new->start_nambi = new->end_nambi = 0;
    new->distance = 0U;
    new->acceptor_distance = new->donor_distance = 0U;

    new->paired_seenp = false;
    new->concordantp = false;

    return new;
  }
}



T
Stage3_new_insertion (int *found_score, int nindels, int indel_pos, int nmismatches1_whole, int nmismatches2_whole,
		      int ncolordiffs1, int ncolordiffs2, Genomicpos_T left, int genomiclength,
		      Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		      int querylength, bool plusp, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
		      int indel_penalty, bool dibasep, bool cmetp) {
  T new;
  Substring_T substring1, substring2;
  int querystart1, queryend1, querystart2, queryend2;
  Genomicpos_T genomicstart, genomicend;
  Genomicpos_T alignstart1, alignend1, alignstart2, alignend2;

  debug0(printf("Stage3_new_insertion: left %u, chrnum %d, indel_pos %d, nindels %d\n",
	       left,chrnum,indel_pos,nindels));

  debug2(printf("Entered with left %u, querylength %d, genomiclength %d, indel_pos %d\n",
		left,querylength,genomiclength,indel_pos));
  debug2(printf("q: %s\n",query));
#if 0
  debug2(printf("g: %s\n",genomicseg));
#endif

  assert(nindels > 0);

  if (dibasep) {
  } else {
    ncolordiffs1 = ncolordiffs2 = 0;
  }

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

  if ((substring1 = Substring_new(nmismatches1_whole,ncolordiffs1,chrnum,
				  chroffset,left,genomicstart,genomicend,query_compress,genome_blocks,snp_blocks,
				  /*start_endtype*/END,/*end_endtype*/INS,
				  querystart1,queryend1,querylength,alignstart1,alignend1,genomiclength,
				  /*extraleft*/0,/*extraright*/0,/*exactp*/false,query,plusp,
				  /*trim_left_p*/true,/*trim_right_p*/false,
				  /*minlength*/0,dibasep,cmetp)) == NULL) {
    return (T) NULL;

  } else if ((substring2 = Substring_new(nmismatches2_whole,ncolordiffs2,chrnum,
					 chroffset,left,genomicstart,genomicend,query_compress,genome_blocks,snp_blocks,
					 /*start_endtype*/INS,/*end_endtype*/END,
					 querystart2,queryend2,querylength,alignstart2,alignend2,genomiclength,
					 /*extraleft*/0,/*extraright*/0,/*exactp*/false,query,plusp,
					 /*trim_left_p*/false,/*trim_right_p*/true,
					 /*minlength*/0,dibasep,cmetp)) == NULL) {
    Substring_free(&substring1);
    return (T) NULL;

  } else {
    new = (T) MALLOC(sizeof(*new));
    new->substring1 = substring1;
    new->substring2 = substring2;
    new->substringD = (Substring_T) NULL;
    new->substringA = (Substring_T) NULL;
    if (plusp == true) {
      new->substring_low = new->substring1;
    } else {
      new->substring_low = new->substring2;
    }

    new->deletion = (char *) NULL;
    new->querylength_adj = querylength - nindels;
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
    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring1) + Substring_mapq_loglik(substring2) + 
      MAPQ_loglik_exact(quality_string,queryend1,querystart2);
    new->mapq_score = 0;
#endif

    new->nindels = nindels;
    new->nmismatches_whole = nmismatches1_whole + nmismatches2_whole;
    new->ntscore = indel_penalty + nmismatches1_whole + nmismatches2_whole;
    new->score = new->ntscore + ncolordiffs1 + ncolordiffs2;

    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(new->substring1) + Substring_nmismatches_bothdiff(new->substring2);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1) + Substring_nmismatches_refdiff(new->substring2); */

#if 0
    /* This method is correct for SNP-tolerant alignment */
    new->nmatches = Substring_match_length(new->substring1) + Substring_match_length(new->substring2) - new->total_nmismatches;
#else
    /* This method is now correct for SNP-tolerant alignment */
    new->nmatches = Substring_nmatches(new->substring1) + Substring_nmatches(new->substring2);
#endif

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->indel_pos = indel_pos;

    new->start_ambiguous_p = new->end_ambiguous_p = false;
    new->start_ambi = new->end_ambi = (int *) NULL;
    new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
    new->start_nambi = new->end_nambi = 0;
    new->distance = 0U;
    new->acceptor_distance = new->donor_distance = 0U;

    new->paired_seenp = false;
    new->concordantp = false;

    return new;
  }
}


T
Stage3_new_deletion (int *found_score, int nindels, int indel_pos, int nmismatches1_whole, int nmismatches2_whole,
		     int ncolordiffs1, int ncolordiffs2, Genomicpos_T left, int genomiclength,
		     Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		     int querylength, bool plusp, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
		     int indel_penalty, bool dibasep, bool cmetp) {
  T new;
  Substring_T substring1, substring2;
  int querystart1, queryend1, querystart2, queryend2;
  Genomicpos_T genomicstart, genomicend;
  Genomicpos_T alignstart1, alignend1, alignstart2, alignend2;
  Genomicpos_T left2;

  debug0(printf("Stage3_new_deletion: left %u, chrnum %d, nmismatches %d, indel_pos %d, nindels %d\n",
	       left,chrnum,indel_pos,nindels));

  debug3(printf("Entered with left %u, querylength %d, genomiclength %d, indel_pos %d\n",
		left,querylength,genomiclength,indel_pos));
  debug3(printf("q: %s\n",query));
  debug3(printf("g: %s\n",genomicseg));

  assert(nindels > 0);

  if (dibasep) {
  } else {
    ncolordiffs1 = ncolordiffs2 = 0;
  }

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


  if ((substring1 = Substring_new(nmismatches1_whole,ncolordiffs1,chrnum,chroffset,
				  left,genomicstart,genomicend,query_compress,genome_blocks,snp_blocks,
				  /*start_endtype*/END,/*end_endtype*/DEL,
				  querystart1,queryend1,querylength,alignstart1,alignend1,genomiclength,
				  /*extraleft*/0,/*extraright*/0,/*exactp*/false,query,plusp,
				  /*trim_left_p*/true,/*trim_right_p*/false,
				  /*minlength*/0,dibasep,cmetp)) == NULL) {
    return (T) NULL;

  } else if ((substring2 = Substring_new(nmismatches2_whole,ncolordiffs2,chrnum,chroffset,
					 left,genomicstart,genomicend,query_compress,genome_blocks,snp_blocks,
					 /*start_endtype*/DEL,/*end_endtype*/END,
					 querystart2,queryend2,querylength,alignstart2,alignend2,genomiclength,
					 /*extraleft*/0,/*extraright*/0,/*exactp*/false,query,plusp,
					 /*trim_left_p*/false,/*trim_right_p*/true,
					 /*minlength*/0,dibasep,cmetp)) == NULL) {
    Substring_free(&substring1);
    return (T) NULL;
    
  } else {
    new = (T) MALLOC(sizeof(*new));
    new->substring1 = substring1;
    new->substring2 = substring2;
    new->substringD = (Substring_T) NULL;
    new->substringA = (Substring_T) NULL;

#if 0
    new->deletion = (char *) CALLOC(nindels+1,sizeof(char));
    if (plusp == true) {
      strncpy(new->deletion,&(genomicseg[indel_pos]),nindels);
      new->substring_low = new->substring1;
    } else {
      make_complement_buffered(new->deletion,&(genomicseg[querylength-indel_pos]),nindels);
      new->substring_low = new->substring2;
    }
#else
    /* Initialize so Substring_free will not try to free */
    new->deletion = (char *) NULL;
    if (plusp == true) {
      new->substring_low = new->substring1;
    } else {
      new->substring_low = new->substring2;
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
    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

#if 0
    new->mapq_loglik = Substring_mapq_loglik(substring1) + Substring_mapq_loglik(substring2);
    new->mapq_score = 0;
#endif

    new->nindels = nindels;
    new->nmismatches_whole = nmismatches1_whole + nmismatches2_whole;
    new->ntscore = indel_penalty + nmismatches1_whole + nmismatches2_whole;
    new->score = new->ntscore + ncolordiffs1 + ncolordiffs2;

    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(new->substring1) + Substring_nmismatches_bothdiff(new->substring2);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(new->substring1) + Substring_nmismatches_refdiff(new->substring2); */

#if 0
    /* This method is correct for SNP-tolerant alignment */
    new->nmatches = Substring_match_length(new->substring1) + Substring_match_length(new->substring2) - new->total_nmismatches;
#else
    /* This method is now correct for SNP-tolerant alignment */
    new->nmatches = Substring_nmatches(new->substring1) + Substring_nmatches(new->substring2);
#endif

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->indel_pos = indel_pos;

    new->start_ambiguous_p = new->end_ambiguous_p = false;
    new->start_ambi = new->end_ambi = (int *) NULL;
    new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
    new->start_nambi = new->end_nambi = 0;
    new->distance = 0U;
    new->acceptor_distance = new->donor_distance = 0U;

    new->paired_seenp = false;
    new->concordantp = false;

    return new;
  }
}


/* Never returns NULL */
T
Stage3_new_splice (int *found_score, int nmismatches_donor, int nmismatches_acceptor,
		   Substring_T donor, Substring_T acceptor, Genomicpos_T distance,
		   bool shortdistancep, int splicing_penalty, int querylength,
		   Intlist_T ambi_left, Intlist_T ambi_right,
		   Intlist_T amb_nmismatches_left, Intlist_T amb_nmismatches_right,
		   bool copy_donor_p, bool copy_acceptor_p, bool first_read_p, int sensedir) {
  T new;
  int ignore;
  Substring_T substring_for_concordance;
  
  debug0(printf("Stage3_new_splice\n"));

  new = (T) MALLOC(sizeof(*new));
  new->deletion = (char *) NULL;
  new->querylength_adj = querylength;

  if (donor == NULL) {
    new->hittype = HALFSPLICE_ACCEPTOR;
    new->substring1 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
    new->substring2 = (Substring_T) NULL;
    
  } else if (acceptor == NULL) {
    new->hittype = HALFSPLICE_DONOR;
    new->substring1 = copy_donor_p ? Substring_copy(donor) : donor;
    new->substring2 = (Substring_T) NULL;

  } else {
    new->hittype = SPLICE;
    new->substring1 = copy_donor_p ? Substring_copy(donor) : donor;
    new->substring2 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
  }
  new->substringD = (Substring_T) NULL;
  new->substringA = (Substring_T) NULL;
  new->nindels = 0;


  if (donor == NULL) {
    new->chrnum = Substring_chrnum(acceptor);
    new->chroffset = Substring_chroffset(acceptor);
    new->plusp = Substring_plusp(acceptor);

  } else if (acceptor == NULL) {
    new->chrnum = Substring_chrnum(donor);
    new->chroffset = Substring_chroffset(donor);
    new->plusp = Substring_plusp(donor);

  } else if (shortdistancep == true) {
    new->chrnum = Substring_chrnum(donor);
    new->chroffset = Substring_chroffset(donor);
    new->plusp = Substring_plusp(donor);

  } else {
    if (Substring_chrnum(donor) == Substring_chrnum(acceptor)) {
      new->chrnum = Substring_chrnum(donor);
      new->chroffset = Substring_chroffset(donor);
    } else {
      new->chrnum = 0;
      new->chroffset = 0;
    }

    if (Substring_plusp(donor) == Substring_plusp(acceptor)) {
      new->plusp = Substring_plusp(donor);
    } else {
      /* Not sure what to do here.  Probably need to have substring->dir rather than substring->plusp. */
      /* Look at ss.samechr for an example.  plusp true => pair_inversion, plusp false => pair_scramble. */
      new->plusp = true;
    }
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

  if (new->genomicstart < new->genomicend) {
    new->low = new->genomicstart;
    new->high = new->genomicend;

    new->start_ambiguous_p = (ambi_left != NULL) ? true : false;
    new->end_ambiguous_p = (ambi_right != NULL) ? true : false;
    new->start_ambi = Intlist_to_array(&new->start_nambi,ambi_left);
    new->end_ambi = Intlist_to_array(&new->end_nambi,ambi_right);
    new->start_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_left);
    new->end_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_right);

  } else {
    new->low = new->genomicend;
    new->high = new->genomicstart;

    new->start_ambiguous_p = (ambi_right != NULL) ? true : false;
    new->end_ambiguous_p = (ambi_left != NULL) ? true : false;
    new->start_ambi = Intlist_to_array(&new->start_nambi,ambi_right);
    new->end_ambi = Intlist_to_array(&new->end_nambi,ambi_left);
    new->start_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_right);
    new->end_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_left);
  }


  if (new->chrnum == 0) {
    new->substring_low = (Substring_T) NULL;
    
    /* Always want the original query end */
    if (first_read_p == true) {
      if (invert_first_p == false) {
	if (Substring_queryend(acceptor) > Substring_queryend(donor)) {
	  substring_for_concordance = acceptor;
	  new->effective_chrnum = Substring_chrnum(acceptor);
	  new->other_chrnum = Substring_chrnum(donor);
	} else {
	  substring_for_concordance = donor;
	  new->effective_chrnum = Substring_chrnum(donor);
	  new->other_chrnum = Substring_chrnum(acceptor);
	}
      } else {
	if (Substring_querystart(acceptor) < Substring_querystart(donor)) {
	  substring_for_concordance = acceptor;
	  new->effective_chrnum = Substring_chrnum(acceptor);
	  new->other_chrnum = Substring_chrnum(donor);
	} else {
	  substring_for_concordance = donor;
	  new->effective_chrnum = Substring_chrnum(donor);
	  new->other_chrnum = Substring_chrnum(acceptor);
	}
      }

    } else {
      if (invert_second_p == false) {
	if (Substring_queryend(acceptor) > Substring_queryend(donor)) {
	  substring_for_concordance = acceptor;
	  new->effective_chrnum = Substring_chrnum(acceptor);
	  new->other_chrnum = Substring_chrnum(donor);
	} else {
	  substring_for_concordance = donor;
	  new->effective_chrnum = Substring_chrnum(donor);
	  new->other_chrnum = Substring_chrnum(acceptor);
	}
      } else {
	if (Substring_querystart(acceptor) < Substring_querystart(donor)) {
	  substring_for_concordance = acceptor;
	  new->effective_chrnum = Substring_chrnum(acceptor);
	  new->other_chrnum = Substring_chrnum(donor);
	} else {
	  substring_for_concordance = donor;
	  new->effective_chrnum = Substring_chrnum(donor);
	  new->other_chrnum = Substring_chrnum(acceptor);
	}
      }
    }

    new->genomicstart = Substring_genomicstart(substring_for_concordance);
    new->genomicend = Substring_genomicend(substring_for_concordance);

  } else {
    new->effective_chrnum = new->chrnum;
    new->other_chrnum = 0;

    if (donor == NULL) {
      new->substring_low = new->substring1;
    } else if (acceptor == NULL) {
      new->substring_low = new->substring1;
    } else if (sensedir == SENSE_FORWARD) {
      if (new->plusp == true) {
	new->substring_low = new->substring1; /* donor */
      } else {
	new->substring_low = new->substring2; /* acceptor */
      }

    } else if (sensedir == SENSE_ANTI) {
      if (new->plusp == true) {
	new->substring_low = new->substring2; /* acceptor */
      } else {
	new->substring_low = new->substring1; /* donor */
      }

    } else {
      abort();
    }
  }

  new->nmismatches_whole = nmismatches_donor + nmismatches_acceptor;
  new->ntscore = splicing_penalty + new->nmismatches_whole;

  if (donor == NULL) {
    /* new->mapq_loglik = Substring_mapq_loglik(acceptor); */
    new->score = new->ntscore + Substring_ncolordiffs(acceptor) + nmismatches_donor;
    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(acceptor) + nmismatches_donor;
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(acceptor) + nmismatches_donor; */
    new->nmatches = Substring_nmatches(acceptor);
  } else if (acceptor == NULL) {
    /* new->mapq_loglik = Substring_mapq_loglik(donor); */
    new->score = new->ntscore + Substring_ncolordiffs(donor);
    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(donor) + nmismatches_acceptor;
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(donor) + nmismatches_acceptor; */
    new->nmatches = Substring_nmatches(donor);
  } else {
    /* new->mapq_loglik = Substring_mapq_loglik(donor) + Substring_mapq_loglik(acceptor); */
    new->score = new->ntscore + Substring_ncolordiffs(donor) + Substring_ncolordiffs(acceptor);
    new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(donor) + Substring_nmismatches_bothdiff(acceptor);
    /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(donor) + Substring_nmismatches_refdiff(acceptor); */
    new->nmatches = Substring_nmatches(donor) + Substring_nmatches(acceptor);
  }
  new->mapq_score = 0;

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->distance = distance;
  new->acceptor_distance = new->donor_distance = 0U;
  new->sensedir = sensedir;

  new->paired_seenp = false;
  new->concordantp = false;

  return new;
}



/* Never returns NULL.  Never copies substrings.  Always shortdistance. */
T
Stage3_new_shortexon (int *found_score, Substring_T donor, Substring_T acceptor, Substring_T shortexon,
		      Genomicpos_T acceptor_distance, Genomicpos_T donor_distance,
		      Intlist_T ambi_left, Intlist_T ambi_right,
		      Intlist_T amb_nmismatches_left, Intlist_T amb_nmismatches_right,
		      bool copy_donor_p, bool copy_acceptor_p, bool copy_shortexon_p,
		      int splicing_penalty, int querylength, bool first_read_p, int sensedir) {
  T new;
  int ignore;
  
  debug0(printf("Stage3_new_shortexon\n"));

  new = (T) MALLOC(sizeof(*new));
  new->deletion = (char *) NULL;
  new->querylength_adj = querylength;

  if (donor == NULL && acceptor == NULL) {
    new->hittype = ONE_THIRD_SHORTEXON;
    new->substring1 = copy_shortexon_p ? Substring_copy(shortexon) : shortexon;
    new->substring2 = (Substring_T) NULL;
    new->substringD = (Substring_T) NULL;
    new->substringA = (Substring_T) NULL;

  } else if (donor == NULL || acceptor == NULL) {
    new->hittype = TWO_THIRDS_SHORTEXON;
    new->substring1 = copy_shortexon_p ? Substring_copy(shortexon) : shortexon;
    new->substring2 = (Substring_T) NULL;
    new->substringD = copy_donor_p ? Substring_copy(donor) : donor;
    new->substringA = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;

  } else {
    new->hittype = SHORTEXON;
    new->substring1 = copy_shortexon_p ? Substring_copy(shortexon) : shortexon;
    new->substring2 = (Substring_T) NULL;
    new->substringD = copy_donor_p ? Substring_copy(donor) : donor;
    new->substringA = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;
  }

  new->nindels = 0;


  new->chrnum = Substring_chrnum(shortexon);
  new->chroffset = Substring_chroffset(shortexon);
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

    new->start_ambiguous_p = (ambi_left != NULL) ? true : false;
    new->end_ambiguous_p = (ambi_right != NULL) ? true : false;
    new->start_ambi = Intlist_to_array(&new->start_nambi,ambi_left);
    new->end_ambi = Intlist_to_array(&new->end_nambi,ambi_right);
    new->start_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_left);
    new->end_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_right);

  } else {
    new->low = new->genomicend;
    new->high = new->genomicstart;

    new->start_ambiguous_p = (ambi_right != NULL) ? true : false;
    new->end_ambiguous_p = (ambi_left != NULL) ? true : false;
    new->start_ambi = Intlist_to_array(&new->start_nambi,ambi_right);
    new->end_ambi = Intlist_to_array(&new->end_nambi,ambi_left);
    new->start_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_right);
    new->end_amb_nmismatches = Intlist_to_array(&ignore,amb_nmismatches_left);
  }

  new->effective_chrnum = new->chrnum;
  new->other_chrnum = 0;

  /* Currently not allowing translocations on shortexons */
  /* substring_for_concordance = (Substring_T) NULL; */

  if (sensedir == SENSE_FORWARD) {
    if (new->plusp == true) {
      new->substring_low = (new->substringD != NULL ? new->substringD : new->substring1); /* donor */
    } else {
      new->substring_low = (new->substringA != NULL ? new->substringA : new->substring1); /* acceptor */
    }

  } else if (sensedir == SENSE_ANTI) {
    if (new->plusp == true) {
      new->substring_low = (new->substringA != NULL ? new->substringA : new->substring1); /* acceptor */
    } else {
      new->substring_low = (new->substringD != NULL ? new->substringD : new->substring1); /* donor */
    }

  } else {
    abort();
  }

#if 0
  new->mapq_loglik = Substring_mapq_loglik(donor) + Substring_mapq_loglik(acceptor) + Substring_mapq_loglik(shortexon);
  new->mapq_score = 0;
#endif

  new->nmismatches_whole = Substring_nmismatches_whole(shortexon);
  if (donor != NULL) {
    new->nmismatches_whole += Substring_nmismatches_whole(donor);
  }
  if (acceptor != NULL) {
    new->nmismatches_whole += Substring_nmismatches_whole(acceptor);
  }
  new->ntscore = splicing_penalty + splicing_penalty + new->nmismatches_whole;

  new->score = new->ntscore + Substring_ncolordiffs(shortexon);
  if (donor != NULL) {
    new->score += Substring_ncolordiffs(donor);
  }
  if (acceptor != NULL) {
    new->score += Substring_ncolordiffs(acceptor);
  }

  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(shortexon);
  if (donor != NULL) {
    new->nmismatches_bothdiff += Substring_nmismatches_bothdiff(donor);
  }
  if (acceptor != NULL) {
    new->nmismatches_bothdiff += Substring_nmismatches_bothdiff(acceptor);
  }
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(donor) + Substring_nmismatches_refdiff(acceptor) + Substring_nmismatches_refdiff(shortexon); */

  new->nmatches = Substring_nmatches(shortexon);
  if (donor != NULL) {
    new->nmatches += Substring_nmatches(donor);
  }
  if (acceptor != NULL) {
    new->nmatches += Substring_nmatches(acceptor);
  }

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->distance = acceptor_distance + donor_distance;
  new->acceptor_distance = acceptor_distance;
  new->donor_distance = donor_distance;
  new->sensedir = sensedir;

  new->paired_seenp = false;
  new->concordantp = false;

  return new;
}


T
Stage3_new_terminal (int querystart, int queryend, int nmismatches_whole, int ncolordiffs,
		     Genomicpos_T left, Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		     int querylength, bool plusp, Endtype_T start_endtype, Endtype_T end_endtype, char *query, 
		     Chrnum_T chrnum, Genomicpos_T chroffset,
		     int terminal_penalty, int max_mismatches_allowed,
		     bool dibasep, bool cmetp) {
  T new;
  Substring_T substring;
  Genomicpos_T genomicstart, genomicend, alignstart, alignend;

  debug0(printf("Stage3_new_terminal: left %u, chrnum %d, querystart %d, queryend %d\n",
	       left,chrnum,querystart,queryend));

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

  if ((substring = Substring_new(nmismatches_whole,ncolordiffs,chrnum,chroffset,
				 left,genomicstart,genomicend,query_compress,genome_blocks,snp_blocks,
				 start_endtype,end_endtype,querystart,queryend,querylength,
				 alignstart,alignend,/*genomiclength*/querylength,
				 /*extraleft*/0,/*extraright*/0,/*exactp*/false,query,plusp,
				 /*trim_left_p*/true,/*trim_right_p*/true,
				 /*minlength*/querylength/2,dibasep,cmetp)) == NULL) {
    return (T) NULL;
    
  } else if (Substring_nmismatches_whole(substring) > max_mismatches_allowed) {
    Substring_free(&substring);
    return (T) NULL;
  }

  if (start_endtype == TERM && end_endtype == TERM) {
    /* This is no longer allowed */
    abort();
    if (Substring_trim_left(substring) == 0 && Substring_trim_right(substring) == 0) {
      Substring_free(&substring);
      return (T) NULL;
    } else if (Substring_trim_left(substring) == 0) {
      Substring_set_start_endtype(substring,END);
    } else if (Substring_trim_right(substring) == 0) {
      Substring_set_end_endtype(substring,END);
    }

  } else if (start_endtype == TERM) {
    if (Substring_trim_left(substring) == 0) {
      Substring_free(&substring);
      return (T) NULL;
    }

  } else if (end_endtype == TERM) {
    if (Substring_trim_right(substring) == 0) {
      Substring_free(&substring);
      return (T) NULL;
    }

  } else {
    abort();
  }

  new = (T) MALLOC(sizeof(*new));
  new->substring1 = substring;
  new->substring2 = (Substring_T) NULL;
  new->substringD = (Substring_T) NULL;
  new->substringA = (Substring_T) NULL;
  new->substring_low = new->substring1;

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
  new->chrnum = new->effective_chrnum = chrnum;
  new->other_chrnum = 0;
  new->chroffset = chroffset;
  new->plusp = plusp;

#if 0
  new->mapq_loglik = Substring_mapq_loglik(substring);
  new->mapq_score = 0;
#endif

  new->nindels = 0;
  new->nmismatches_whole = nmismatches_whole;
  new->ntscore = terminal_penalty + nmismatches_whole;

  new->nmismatches_bothdiff = Substring_nmismatches_bothdiff(substring);
  /* new->nmismatches_refdiff = Substring_nmismatches_refdiff(substring); */

#if 0
  new->score = terminal_penalty + nmismatches_whole + ncolordiffs;
#else
  new->score = terminal_penalty + Substring_nmismatches_whole(substring);
#endif

#if 0
  new->nmatches = Substring_match_length(substring) - (nmismatches + ncolordiffs);
#else
  new->nmatches = Substring_nmatches(substring);
#endif
  new->start_ambiguous_p = new->end_ambiguous_p = false;
  new->start_ambi = new->end_ambi = (int *) NULL;
  new->start_amb_nmismatches = new->end_amb_nmismatches = (int *) NULL;
  new->start_nambi = new->end_nambi = 0;

  new->distance = 0U;
  new->acceptor_distance = new->donor_distance = 0U;
  new->sensedir = SENSE_NULL;

  new->paired_seenp = false;
  new->concordantp = false;

  return new;
}


static int
Stage3_output_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->mapq_loglik > y->mapq_loglik) {
    return -1;
  } else if (y->mapq_loglik > x->mapq_loglik) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
    return +1;
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (y->hittype < x->hittype) {
    return +1;
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

  if (x->absdifflength_bingo_p == true && y->absdifflength_bingo_p == false) {
    return -1;
  } else if (y->absdifflength_bingo_p == true && x->absdifflength_bingo_p == false) {
    return +1;
  } else if (x->insertlength > 0 && y->insertlength == 0) {
    return -1;
  } else if (y->insertlength > 0 && x->insertlength == 0) {
    return +1;
  } else if (x->insertlength < y->insertlength) {
    return -1;
  } else if (y->insertlength < x->insertlength) {
    return +1;
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
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
Stage3_compute_mapq (Stage3_T this, char *query, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     UINT4 *genome_blocks, UINT4 *snp_blocks, char *quality_string, bool dibasep, bool cmetp) {

  if (this == NULL) {
    return 0.0;

  } else if (this->plusp == true) {
    this->mapq_loglik =
      Substring_compute_mapq(this->substring1,query,query_compress_fwd,
			     genome_blocks,snp_blocks,quality_string,dibasep,cmetp);

    if (this->substring2 != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substring2,query,query_compress_fwd,
			       genome_blocks,snp_blocks,quality_string,dibasep,cmetp);
    }
    if (this->substringD != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substringD,query,query_compress_fwd,
			       genome_blocks,snp_blocks,quality_string,dibasep,cmetp);
    }
    if (this->substringA != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substringA,query,query_compress_fwd,
			       genome_blocks,snp_blocks,quality_string,dibasep,cmetp);
    }

  } else {
    this->mapq_loglik =
      Substring_compute_mapq(this->substring1,query,query_compress_rev,
			     genome_blocks,snp_blocks,quality_string,dibasep,cmetp);

    if (this->substring2 != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substring2,query,query_compress_rev,
			       genome_blocks,snp_blocks,quality_string,dibasep,cmetp);
    }
    if (this->substringD != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substringD,query,query_compress_rev,
			       genome_blocks,snp_blocks,quality_string,dibasep,cmetp);
    }
    if (this->substringA != NULL) {
      this->mapq_loglik +=
	Substring_compute_mapq(this->substringA,query,query_compress_rev,
			       genome_blocks,snp_blocks,quality_string,dibasep,cmetp);
    }
  }

  return this->mapq_loglik;
}



static void
Stage3_display_prep (Stage3_T this, char *query, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome, bool dibasep, bool cmetp) {
  char *deletion_ignore;

  if (this != NULL) {
    if (this->hittype == DELETION) {
      this->nmismatches_refdiff = 
	Substring_display_prep(&this->deletion,this->substring1,query,query_compress_fwd,query_compress_rev,
			       genome_blocks,snp_blocks,genome,/*deletion_pos*/this->indel_pos,
			       /*deletion_length*/this->nindels,dibasep,cmetp);
    } else {
      this->nmismatches_refdiff = 
	Substring_display_prep(&deletion_ignore,this->substring1,query,query_compress_fwd,query_compress_rev,
			       genome_blocks,snp_blocks,genome,/*deletion_pos*/-1,/*deletion_length*/0,dibasep,cmetp);
    }

    if (this->substring2 != NULL) {
      this->nmismatches_refdiff +=
	Substring_display_prep(&deletion_ignore,this->substring2,query,query_compress_fwd,query_compress_rev,
			       genome_blocks,snp_blocks,genome,/*deletion_pos*/-1,/*deletion_length*/0,dibasep,cmetp);
    }
    if (this->substringD != NULL) {
      this->nmismatches_refdiff +=
	Substring_display_prep(&deletion_ignore,this->substringD,query,query_compress_fwd,query_compress_rev,
			       genome_blocks,snp_blocks,genome,/*deletion_pos*/-1,/*deletion_length*/0,dibasep,cmetp);
    }
    if (this->substringA != NULL) {
      this->nmismatches_refdiff +=
	Substring_display_prep(&deletion_ignore,this->substringA,query,query_compress_fwd,query_compress_rev,
			       genome_blocks,snp_blocks,genome,/*deletion_pos*/-1,/*deletion_length*/0,dibasep,cmetp);
    }
  }
  return;
}


Stage3_T *
Stage3_eval_and_sort (Stage3_T *stage3array, int npaths, int maxpaths, Shortread_T queryseq,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,
		      UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome,
		      char *quality_string, bool dibasep, bool cmetp) {
  double maxlik, total;
  char *query;
  int i;

  if (npaths == 0) {
    /* Skip */

  } else if (npaths == 1) {
    stage3array[0]->mapq_score = 
      MAPQ_max_quality_score(Shortread_quality_string(queryseq),Shortread_fulllength(queryseq));
    query = Shortread_fullpointer_uc(queryseq);
    Stage3_display_prep(stage3array[0],query,query_compress_fwd,query_compress_rev,
			genome_blocks,snp_blocks,genome,dibasep,cmetp);

  } else {
    /* Compute mapq_loglik */
    query = Shortread_fullpointer_uc(queryseq);
    for (i = 0; i < npaths; i++) {
      Stage3_compute_mapq(stage3array[i],query,query_compress_fwd,query_compress_rev,
			  genome_blocks,snp_blocks,quality_string,dibasep,cmetp);
    }

    qsort(stage3array,npaths,sizeof(Stage3_T),Stage3_output_cmp);
    maxlik = stage3array[0]->mapq_loglik;
    
    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < npaths; i++) {
      stage3array[i]->mapq_loglik -= maxlik;
    }

    total = 0.0;
    for (i = 0; i < npaths; i++) {
      total += (stage3array[i]->mapq_loglik = exp(stage3array[i]->mapq_loglik));
    }

    /* Save on computation if possible */
    if (maxpaths < npaths) npaths = maxpaths;

    /* Prepare for display */
    for (i = 0; i < npaths; i++) {
      Stage3_display_prep(stage3array[i],query,query_compress_fwd,query_compress_rev,
			  genome_blocks,snp_blocks,genome,dibasep,cmetp);
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < npaths; i++) {
      stage3array[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < npaths; i++) {
      stage3array[i]->mapq_score = rint(-10.0 * log10(1.0 - stage3array[i]->mapq_loglik));
      /* printf("Converting %.2f to %d\n",stage3array[i]->mapq_loglik,stage3array[i]->mapq_score); */
    }

  }

  return stage3array;
}



List_T
Stage3_optimal_score (List_T hitlist, int cutoff_level, int suboptimal_mismatches) {
  List_T optimal = NULL, p;
  T hit;
  int n;
  int minscore = MAX_QUERYLENGTH;
  bool non_translocation_p = false;

  n = List_length(hitlist);
  if (n == 0) {
    return NULL;
  }

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->chrnum != 0) {
      non_translocation_p = true;
    }
  }

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->chrnum == 0 && non_translocation_p == true) {
      /* Skip, since we will eliminate */
    } else if (hit->hittype == TERMINAL) {
      /* For dibasep were previously using hit->ntscore, but gives false positives */
      /* Skip from setting minscore */
    } else if (hit->score <= cutoff_level) {
#ifdef CLEAN_SINGLE_END_AMBIGUOUS
      if (hit->chimera_ambiguous_p == false) {
#endif
	if (hit->score < minscore) {
	  minscore = hit->score;
	}
#ifdef CLEAN_SINGLE_END_AMBIGUOUS
      }
#endif
    }
  }

  debug(printf("Stage3_optimal_score over %d hits: minscore = %d + subopt:%d\n",
	       n,minscore,suboptimal_mismatches));
  minscore += suboptimal_mismatches;

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    if (hit->chrnum == 0 && non_translocation_p == true) {
      debug(printf("Eliminating a hit with splice translocation\n"));
      Stage3_free(&hit);
    } else if (hit->score > cutoff_level) {
      /* For dibasep were previously using hit->ntscore, but gives false positives */
      debug(printf("Eliminating a hit with ntscore %d > cutoff_level %d\n",hit->ntscore,cutoff_level));
      Stage3_free(&hit);
#ifdef CLEAN_SINGLE_END_AMBIGUOUS
    } else if (hit->chimera_ambiguous_p == true) {
      debug(printf("Eliminating a hit with ntscore %d and ambiguous splice\n",hit->ntscore));
      Stage3_free(&hit);
#endif
    } else if (hit->score <= minscore) {
      debug(printf("Keeping a hit with score %d\n",hit->score));
      optimal = List_push(optimal,hit);
    } else {
      debug(printf("Eliminating a hit with score %d\n",hit->score));
      Stage3_free(&hit);
    }
  }
  
  List_free(&hitlist);

  return optimal;
}


int
Stage3_noptimal (List_T hitlist) {
  int noptimal;
  List_T p;
  T hit;
  int minscore = MAX_QUERYLENGTH;

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
genomicstart_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
#if 0
  } else if (x->genomicend < y->genomicend) {
    return -1;
  } else if (x->genomicend > y->genomicend) {
    return +1;
#endif
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
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
#if 0
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
#endif
  } else if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else {
    return 0;
  }
}


#if 0

List_T
Stage3_mark_ambiguous_splices (bool *ambiguousp, List_T hitlist) {
  T x, y, *hits;
  int n, i, j;
  Genomicpos_T splice_distance_1, splice_distance_2;

#ifndef CLEAN_SINGLE_END_AMBIGUOUS
  return hitlist;
#endif

  *ambiguousp = false;
  n = List_length(hitlist);
  debug9(printf("Entered Stage3_mark_ambiguous_splices with %d pairs\n",n));
  if (n == 0) {
    return NULL;
  } else {
    hits = (T *) List_to_array(hitlist,NULL);
  }

  /* By genomicstart */
  debug9(printf("Stage3_mark_ambiguous_splices: checking %d hits by genomicstart\n",n));
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
  debug9(printf("Stage3_mark_ambiguous_splices: checking %d hits by genomicend\n",n));
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


List_T
Stage3_remove_duplicates (List_T hitlist) {
#ifdef DEBUG7
  List_T p;
#endif
  List_T unique = NULL;
  T x, y, hit, *hits, *prev;
  int n, i, j, nkept;
  bool *eliminate;

  n = List_length(hitlist);
  debug7(printf("Entered Stage3_remove_duplicates with %d hits\n",n));
  if (n == 0) {
    return NULL;
  } else {
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hits = (T *) List_to_array(hitlist,NULL);
    List_free(&hitlist);
  }


  /* By genomicstart */
  debug7(printf("Stage3_remove_duplicates: checking %d hits by genomicstart\n",n));
  qsort(hits,n,sizeof(T),genomicstart_cmp);

  debug7(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d (%s): #%d:%u..%u, nmatches: %d\n",
		  i,hittype_string(hit->hittype),hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,hit->nmatches);
	 }
	 );

  i = 0;
  while (i < n) {
    if (eliminate[i] == true) {
      debug7(printf("Skipping %d because already eliminated\n",i));
      i++;
    } else {
      x = hits[i];
      j = i+1;
      while (j < n && hits[j]->genomicstart == x->genomicstart) {
	debug7(printf("Looking at %d and %d with common start\n",i,j));
	y = hits[j];
	if (SENSE_INCONSISTENT_P(x->sensedir,y->sensedir)) {
	  debug7(printf("  #%d and #%d have different sense, so skipping\n",i,j));

	} else if (x->nmatches > y->nmatches) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d > %d, so marking %d for elimination\n",
			i,j,x->nmatches,y->nmatches,j));
	  eliminate[j] = true;
	} else if (x->nmatches < y->nmatches) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d < %d, so marking %d for elimination\n",
			i,j,x->nmatches,y->nmatches,i));
	  eliminate[i] = true;

	} else if (x->score < y->score) {
	  debug7(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
			i,j,x->score,y->score,j));
	  eliminate[j] = true;
	} else if (x->score > y->score) {
	  debug7(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
			i,j,x->score,y->score,i));
	  eliminate[i] = true;

	} else if (x->hittype < y->hittype) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d == %d, so marking %d for elimination\n",
			i,j,x->nmatches,y->nmatches,j));
	  eliminate[j] = true;
	} else if (x->hittype > y->hittype) {
	  debug7(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
			i,j,x->nmatches,y->nmatches,i));
	  eliminate[i] = true;
	} else if (x->nindels < y->nindels) {
	  debug7(printf("  #%d overlaps #%d and fewer indels, so marking %d for elimination\n",
			i,j,j));
	  eliminate[j] = true;
	} else if (x->nindels > y->nindels) {
	  debug7(printf("  #%d overlaps #%d and fewer indels, so marking %d for elimination\n",
			i,j,i));
	  eliminate[i] = true;
	} else if (x->chrnum == 0 && y->chrnum != 0) {
	  debug7(printf("  #%d overlaps #%d and first one is a distant splice, so marking first one for elimination\n",
			i,j));
	  eliminate[i] = true;
	} else if (x->chrnum > 0 && y->chrnum == 0) {
	  debug7(printf("  #%d overlaps #%d and second one is a distant splice, so marking second one for elimination\n",
			i,j));
	  eliminate[j] = true;

#ifdef CLEAN_SINGLE_END_AMBIGUOUS
	} else if (x->distance > y->distance) {
	  debug7(printf("  #%d overlaps #%d and first one has longer distance, so marking first one for elimination\n",
			i,j));
	  eliminate[i] = true;
	} else if (x->distance < y->distance) {
	  debug7(printf("  #%d overlaps #%d and second one has longer distance, so marking second one for elimination\n",
			i,j));
	  eliminate[j] = true;
#endif

	} else if (x->genomicend == y->genomicend) {
	  debug7(printf("  #%d overlaps #%d and equal so eliminate second one\n",
			i,j));
	  eliminate[j] = true;
	}
	j++;
      }
      i = j;
    }
  }
    
  j = 0;
  prev = hits;
  hits = (T *) CALLOC(n,sizeof(T));

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

  for (i = n-1; i >= 0; i--) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug7(printf("  Keeping #%d:%u..%u, nmatches %d (nindels %d, distance %u, chrnum %d) (plusp = %d) on basis of genomicstart\n",
		    hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		    hit->nmatches,hit->nindels,hit->distance,hit->chrnum,hit->plusp));
      hits[j++] = hit;
    } else {
      debug7(printf("  Eliminating #%d:%u..%u, nmatches %d (nindels %d, distance %u, chrnum %d) (plusp = %d) on basis of genomicstart\n",
		    hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		    hit->nmatches,hit->nindels,hit->distance,hit->chrnum,hit->plusp));
      Stage3_free(&hit);
    }
  }
  
  FREE(prev);
  FREE(eliminate);


  /* By genomicend */
  n = j;
  eliminate = (bool *) CALLOC(n,sizeof(bool));
  debug7(printf("Stage3_remove_duplicates: checking %d hits by genomicend\n",n));
  qsort(hits,n,sizeof(T),genomicend_cmp);
  i = 0;
  while (i < n) {
    if (eliminate[i] == true) {
      debug7(printf("Skipping %d because already eliminated\n",i));
      i++;
    } else {
      x = hits[i];
      j = i+1;
      while (j < n && hits[j]->genomicend == x->genomicend) {
	debug7(printf("Looking at %d and %d with common end\n",i,j));
	y = hits[j];
	if (SENSE_INCONSISTENT_P(x->sensedir,y->sensedir)) {
	  debug7(printf("  #%d and #%d have different sense, so skipping\n",i,j));

	} else if (x->nmatches > y->nmatches) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d > %d, so marking %d for elimination\n",
			i,j,x->nmatches,y->nmatches,j));
	  eliminate[j] = true;
	} else if (x->nmatches < y->nmatches) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d < %d, so marking %d for elimination\n",
			i,j,x->nmatches,y->nmatches,i));
	  eliminate[i] = true;

	} else if (x->score < y->score) {
	  debug7(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
			i,j,x->score,y->score,j));
	  eliminate[j] = true;
	} else if (x->score > y->score) {
	  debug7(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
			i,j,x->score,y->score,i));
	  eliminate[i] = true;

	} else if (x->hittype < y->hittype) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d == %d, so marking %d for elimination\n",
			i,j,x->nmatches,y->nmatches,j));
	  eliminate[j] = true;
	} else if (x->hittype > y->hittype) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d < %d, so marking %d for elimination\n",
			i,j,x->nmatches,y->nmatches,i));
	  eliminate[i] = true;
	} else if (x->nindels < y->nindels) {
	  debug7(printf("  #%d overlaps #%d and fewer indels, so marking %d for elimination\n",
			i,j,j));
	  eliminate[j] = true;
	} else if (x->nindels > y->nindels) {
	  debug7(printf("  #%d overlaps #%d and fewer indels, so marking %d for elimination\n",
			i,j,i));
	  eliminate[i] = true;
	} else if (x->chrnum == 0 && y->chrnum != 0) {
	  debug7(printf("  #%d overlaps #%d and first one is a distant splice, so marking first one for elimination\n",
			i,j));
	  eliminate[i] = true;
	} else if (x->chrnum != 0 && y->chrnum == 0) {
	  debug7(printf("  #%d overlaps #%d and second one is a distant splice, so marking second one for elimination\n",
			i,j));
	  eliminate[j] = true;

#ifdef CLEAN_SINGLE_END_AMBIGUOUS
	} else if (x->distance > y->distance) {
	  debug7(printf("  #%d overlaps #%d and first one has longer distance, so marking first one for elimination\n",
			i,j));
	  eliminate[i] = true;
	} else if (x->distance < y->distance) {
	  debug7(printf("  #%d overlaps #%d and second one has longer distance, so marking second one for elimination\n",
			i,j));
	  eliminate[j] = true;
#endif

	} else if (x->genomicstart == y->genomicstart) {
	  debug7(printf("  #%d overlaps #%d and equal so eliminate second one\n",
			i,j));
	  eliminate[j] = true;
	}
	j++;
      }
      i = j;
    }
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
  }

  for (i = n-1; i >= 0; i--) {
    hit = hits[i];
    if (eliminate[i] == false) {
      debug7(printf("  Keeping #%d:%u..%u, nmatches %d (nindels %d, distance %u, chrnum %d) (plusp = %d)\n",
		    hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		    hit->nmatches,hit->nindels,hit->distance,hit->chrnum,hit->plusp));
      unique = List_push(unique,hit);
    } else {
      debug7(printf("  Eliminating #%d:%u..%u, nmatches %d (nindels %d, distance %u, chrnum %d) (plusp = %d)\n",
		    hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,
		    hit->nmatches,hit->nindels,hit->distance,hit->chrnum,hit->plusp));
      Stage3_free(&hit);
    }
  }
  FREE(hits);
  FREE(eliminate);

  debug7(
	 for (p = unique, i = 0; p != NULL; p = p->rest, i++) {
	   hit = (T) p->first;
	   printf("  Final %d: #%d:%u..%u (plusp = %d)\n",
		  i,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,hit->plusp);
	 }
	 );

  return unique;
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
    if (hit3->genomicstart < hit5->genomicstart) {
      return PAIRED_SCRAMBLE;
    } else {
      return PAIRED_TOOLONG;
    }
  } else {
    if (hit5->genomicstart < hit3->genomicstart) {
      return PAIRED_SCRAMBLE;
    } else {
      return PAIRED_TOOLONG;
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
  case UNPAIRED: abort();
  }

  return;
}

static void
print_single (FILE *fp, T this, int score, IIT_T chromosome_iit, Shortread_T queryseq,
	      bool invertp, T hit5, T hit3, int insertlength, int pairscore, 
	      Pairtype_T pairtype, int mapq_score) {
  int querylength;
  char *chr;
  bool allocp;

  querylength = Shortread_fulllength(queryseq);
  chr = IIT_label(chromosome_iit,this->chrnum,&allocp);

  fprintf(fp," ");
  Substring_print_single(fp,this->substring1,queryseq,chr,querylength,invertp);

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
  Substring_print_deletion_2(fp,this->substring1,this->substring2,this->nindels,this->deletion,
			     queryseq,chr,invertp);
  fprintf(fp,"\n");

  if (allocp == true) {
    FREE(chr);
  }
}


static void
print_splice (FILE *fp, T chimera, int score, UINT4 *genome_blocks,
	      IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, T hit5, T hit3,
	      int insertlength, int pairscore, Pairtype_T pairtype, int mapq_score) {
  Substring_T donor, acceptor;
  
  if (chimera->hittype == HALFSPLICE_DONOR) {
    donor = chimera->substring1;
    acceptor = (Substring_T) NULL;
    Substring_assign_donor_prob(donor,genome_blocks);

  } else if (chimera->hittype == HALFSPLICE_ACCEPTOR) {
    acceptor = chimera->substring1;
    donor = (Substring_T) NULL;
    Substring_assign_acceptor_prob(acceptor,genome_blocks);

  } else {
    donor = chimera->substring1;
    acceptor = chimera->substring2;
    Substring_assign_donor_prob(donor,genome_blocks);
    Substring_assign_acceptor_prob(acceptor,genome_blocks);
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
print_shortexon (FILE *fp, T chimera, int score, UINT4 *genome_blocks,
		 IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, T hit5, T hit3,
		 int insertlength, int pairscore, Pairtype_T pairtype, int mapq_score) {
  Substring_T donor, acceptor, shortexon;
  Genomicpos_T distance1, distance2;
  bool firstp = true;
  int nblocks = 1;
  
  shortexon = chimera->substring1;
  Substring_assign_shortexon_prob(shortexon,genome_blocks);
  if ((donor = chimera->substringD) != NULL) {
    Substring_assign_donor_prob(donor,genome_blocks);
    nblocks++;
  }
  if ((acceptor = chimera->substringA) != NULL) {
    Substring_assign_acceptor_prob(acceptor,genome_blocks);
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
Stage3_print (FILE *fp, T this, int score, UINT4 *genome_blocks,
	      IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, T hit5, T hit3, int insertlength,
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
  } else if (this->hittype == HALFSPLICE_DONOR || this->hittype == HALFSPLICE_ACCEPTOR || this->hittype == SPLICE) {
    print_splice(fp,this,score,genome_blocks,
		 chromosome_iit,queryseq,invertp,hit5,hit3,insertlength,
		 pairscore,pairtype,mapq_score);
  } else if (this->hittype == ONE_THIRD_SHORTEXON || this->hittype == TWO_THIRDS_SHORTEXON || this->hittype == SHORTEXON) {
    print_shortexon(fp,this,score,genome_blocks,
		    chromosome_iit,queryseq,invertp,hit5,hit3,insertlength,
		    pairscore,pairtype,mapq_score);
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
      Shortread_print_quality(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,quality_shift,
			      /*show_chopped_p*/true);
    } else {
      Shortread_print_quality_revcomp(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,quality_shift,
				      /*show_chopped_p*/true);
    }
  }

  return;
}


static void
print_one_paired_end (Result_T result, Resulttype_T resulttype, bool translocationp,
		      char initchar, bool firstp, UINT4 *genome_blocks, IIT_T chromosome_iit,
		      Shortread_T queryseq, Shortread_T queryseq_mate, Shortread_T headerseq,
		      Genomicpos_T pairmax, int maxpaths, bool quiet_if_excessive_p,
		      bool invertp, int quality_shift,
		      FILE *fp_nomapping_1, FILE *fp_unpaired_uniq, FILE *fp_unpaired_mult,
		      FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_mult,
		      FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr, FILE *fp_paired_uniq_long,
		      FILE *fp_paired_mult, FILE *fp_concordant_uniq, FILE *fp_concordant_mult) {
  Stage3pair_T *stage3pairarray, stage3pair;
  T *stage3array, this, hit5, hit3;
  int npaths, pathnum;
  bool outputp;
  FILE *fp;

  if (resulttype == PAIREDEND_NOMAPPING) {
    /* If fails_as_input_p == true, then this case is handled by calling procedure */
    print_query_header(fp_nomapping_1,initchar,queryseq,invertp);
    fprintf(fp_nomapping_1,"\t0 %s",UNPAIRED_TEXT);
    if (translocationp == true) {
      fprintf(fp_nomapping_1," (transloc)");
    }

    print_barcode_and_quality(fp_nomapping_1,queryseq,invertp,quality_shift);
    
    fprintf(fp_nomapping_1,"\t");
    Shortread_print_header(fp_nomapping_1,headerseq);
    fprintf(fp_nomapping_1,"\n");

  } else if (resulttype == CONCORDANT_UNIQ) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);

    print_query_header(fp_concordant_uniq,initchar,queryseq,invertp);
    fprintf(fp_concordant_uniq,"\t1 %s",CONCORDANT_TEXT);
    if (translocationp == true) {
      fprintf(fp_concordant_uniq," (transloc)");
    }
    
    print_barcode_and_quality(fp_concordant_uniq,queryseq,invertp,quality_shift);

    fprintf(fp_concordant_uniq,"\t");
    Shortread_print_header(fp_concordant_uniq,headerseq);

    stage3pair = stage3pairarray[0];
    hit5 = stage3pair->hit5;
    hit3 = stage3pair->hit3;
    
    if (firstp == true) {
#if 0
      Stage3pair_eval(stage3pairarray,/*npaths*/1,maxpaths,queryseq,queryseq_mate);
#endif
      Stage3_print(fp_concordant_uniq,hit5,hit5->score,genome_blocks,
		   chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		   stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    } else {
      Stage3_print(fp_concordant_uniq,hit3,hit3->score,genome_blocks,
		   chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		   stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    }

    fprintf(fp_concordant_uniq,"\n");

  } else if (resulttype == CONCORDANT_MULT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);

    print_query_header(fp_concordant_mult,initchar,queryseq,invertp);
    fprintf(fp_concordant_mult,"\t%d %s",npaths,CONCORDANT_TEXT);
    if (translocationp == true) {
      fprintf(fp_concordant_mult," (transloc)");
    }

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
	  Stage3_print(fp_concordant_mult,hit5,hit5->score,genome_blocks,
		       chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		       stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	} else {
	  Stage3_print(fp_concordant_mult,hit3,hit3->score,genome_blocks,
		       chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		       stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	}
      }
    }

    fprintf(fp_concordant_mult,"\n");

  } else if (resulttype == PAIRED_UNIQ) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);
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
    if (translocationp == true) {
      fprintf(fp," (transloc)");
    }

    print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

    fprintf(fp,"\t");
    Shortread_print_header(fp,headerseq);

    hit5 = stage3pair->hit5;
    hit3 = stage3pair->hit3;
    
    if (firstp == true) {
#if 0
      Stage3pair_eval(stage3pairarray,/*npaths*/1,maxpaths,queryseq,queryseq_mate);
#endif
      Stage3_print(fp,hit5,hit5->score,genome_blocks,
		   chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		   stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    } else {
      Stage3_print(fp,hit3,hit3->score,genome_blocks,
		   chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		   stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
    }

    fprintf(fp,"\n");

  } else if (resulttype == PAIRED_MULT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);

    print_query_header(fp_paired_mult,initchar,queryseq,invertp);
    fprintf(fp_paired_mult,"\t%d %s",npaths,PAIRED_TEXT);
    if (translocationp == true) {
      fprintf(fp_paired_mult," (transloc)");
    }

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
	  Stage3_print(fp_paired_mult,hit5,hit5->score,genome_blocks,
		       chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		       stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	} else {
	  Stage3_print(fp_paired_mult,hit3,hit3->score,genome_blocks,
		       chromosome_iit,queryseq,invertp,hit5,hit3,stage3pair->insertlength,
		       stage3pair->score,stage3pair->pairtype,stage3pair->mapq_score);
	}
      }
    }

    fprintf(fp_paired_mult,"\n");

  } else {
    /* Print as singles */
    if (firstp == true) {
      stage3array = (T *) Result_array(&npaths,result);
    } else {
      stage3array = (T *) Result_array2(&npaths,result);
    }

    outputp = true;
    if (resulttype == HALFMAPPING_UNIQ) {
      fp = fp_halfmapping_uniq;

    } else if (resulttype == HALFMAPPING_MULT) {
      fp = fp_halfmapping_mult;
      if (quiet_if_excessive_p && npaths > maxpaths) {
	outputp = false;
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      fp = fp_unpaired_uniq;

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
      stage3array = (T *) Result_array(&npaths,result);
      hit5 = stage3array[0];
      stage3array = (T *) Result_array2(&npaths,result);
      hit3 = stage3array[0];
      fprintf(fp," (%s)",unpaired_type_text(hit5,hit3));
    }
#endif

    print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

    fprintf(fp,"\t");
    Shortread_print_header(fp,headerseq);

    if (outputp == true) {
      /* Stage3_eval_and_sort(stage3array,npaths,maxpaths,queryseq); */
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	this = stage3array[pathnum-1];
	Stage3_print(fp,this,this->score,genome_blocks,
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
Stage3_print_paired (Result_T result, Resulttype_T resulttype, bool translocationp, UINT4 *genome_blocks,
		     IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		     Genomicpos_T pairmax, int maxpaths, bool quiet_if_excessive_p,
		     bool nofailsp, bool failsonlyp, bool fails_as_input_p,
		     bool fastq_format_p, int quality_shift,
		     FILE *fp_nomapping_1, FILE *fp_nomapping_2,
		     FILE *fp_unpaired_uniq, FILE *fp_unpaired_mult,
		     FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_mult,
		     FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr,
		     FILE *fp_paired_uniq_long, FILE *fp_paired_mult,
		     FILE *fp_concordant_uniq, FILE *fp_concordant_mult) {

  debug1(printf("Stage3_print_paired: resulttype is %s\n",Resulttype_string(resulttype)));

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
      print_one_paired_end(result,resulttype,translocationp,'>',/*firstp*/true,genome_blocks,
			   chromosome_iit,queryseq1,/*queryseq_mate*/queryseq2,/*headerseq*/queryseq1,
			   pairmax,maxpaths,quiet_if_excessive_p,invert_first_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_mult,
			   fp_paired_uniq_inv,fp_paired_uniq_scr,fp_paired_uniq_long,
			   fp_paired_mult,fp_concordant_uniq,fp_concordant_mult);

      /* Second end */
      print_one_paired_end(result,resulttype,translocationp,'<',/*firstp*/false,genome_blocks,
			   chromosome_iit,queryseq2,/*queryseq_mate*/queryseq1,/*headerseq*/queryseq1,
			   pairmax,maxpaths,quiet_if_excessive_p,invert_second_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_mult,
			   fp_paired_uniq_inv,fp_paired_uniq_scr,fp_paired_uniq_long,
			   fp_paired_mult,fp_concordant_uniq,fp_concordant_mult);
    }

  } else {
    if (failsonlyp == true) {
      /* Unwanted success: skip */
      debug1(printf("  failsonlyp is true, so no output\n"));
    
    } else {

      /* First end */
      print_one_paired_end(result,resulttype,translocationp,'>',/*firstp*/true,genome_blocks,
			   chromosome_iit,queryseq1,/*queryseq_mate*/queryseq2,/*headerseq*/queryseq1,
			   pairmax,maxpaths,quiet_if_excessive_p,invert_first_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_mult,
			   fp_paired_uniq_inv,fp_paired_uniq_scr,fp_paired_uniq_long,
			   fp_paired_mult,fp_concordant_uniq,fp_concordant_mult);

      /* Second end */
      print_one_paired_end(result,resulttype,translocationp,'<',/*firstp*/false,genome_blocks,
			   chromosome_iit,queryseq2,/*queryseq_mate*/queryseq1,/*headerseq*/queryseq1,
			   pairmax,maxpaths,quiet_if_excessive_p,invert_second_p,
			   quality_shift,fp_nomapping_1,fp_unpaired_uniq,fp_unpaired_mult,
			   fp_halfmapping_uniq,fp_halfmapping_mult,
			   fp_paired_uniq_inv,fp_paired_uniq_scr,fp_paired_uniq_long,
			   fp_paired_mult,fp_concordant_uniq,fp_concordant_mult);
    }
  }
    
  return;
}


#if 0
List_T
Stage3_filter_bymatch (List_T hitlist) {
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
      Stage3_free(&hit);
    }
  }
  List_free(&hitlist);

  return filtered;
}
#endif



static Genomicpos_T
pair_insert_length (Stage3_T hit5, Stage3_T hit3) {

  if (Substring_overlap_p(hit5->substring1,hit3->substring1)) {
    return Substring_insert_length(hit5->substring1,hit3->substring1);
  } else if (hit5->substring2 != NULL && Substring_overlap_p(hit5->substring2,hit3->substring1)) {
    return Substring_insert_length(hit5->substring2,hit3->substring1);
  } else if (hit5->substringD != NULL && Substring_overlap_p(hit5->substringD,hit3->substring1)) {
    return Substring_insert_length(hit5->substringD,hit3->substring1);
  } else if (hit5->substringA != NULL && Substring_overlap_p(hit5->substringA,hit3->substring1)) {
    return Substring_insert_length(hit5->substringA,hit3->substring1);
  }
  
  if (hit3->substring2 != NULL) {
    if (Substring_overlap_p(hit5->substring1,hit3->substring2)) {
      return Substring_insert_length(hit5->substring1,hit3->substring2);
    } else if (hit5->substring2 != NULL && Substring_overlap_p(hit5->substring2,hit3->substring2)) {
      return Substring_insert_length(hit5->substring2,hit3->substring2);
    } else if (hit5->substringD != NULL && Substring_overlap_p(hit5->substringD,hit3->substring2)) {
      return Substring_insert_length(hit5->substringD,hit3->substring2);
    } else if (hit5->substringA != NULL && Substring_overlap_p(hit5->substringA,hit3->substring2)) {
      return Substring_insert_length(hit5->substringA,hit3->substring2);
    }
  }

  if (hit3->substringD != NULL) {
    if (Substring_overlap_p(hit5->substring1,hit3->substringD)) {
      return Substring_insert_length(hit5->substring1,hit3->substringD);
    } else if (hit5->substring2 != NULL && Substring_overlap_p(hit5->substring2,hit3->substringD)) {
      return Substring_insert_length(hit5->substring2,hit3->substringD);
    } else if (hit5->substringD != NULL && Substring_overlap_p(hit5->substringD,hit3->substringD)) {
      return Substring_insert_length(hit5->substringD,hit3->substringD);
    } else if (hit5->substringA != NULL && Substring_overlap_p(hit5->substringA,hit3->substringD)) {
      return Substring_insert_length(hit5->substringA,hit3->substringD);
    }
  }

  if (hit3->substringA != NULL) {
    if (Substring_overlap_p(hit5->substring1,hit3->substringA)) {
      return Substring_insert_length(hit5->substring1,hit3->substringA);
    } else if (hit5->substring2 != NULL && Substring_overlap_p(hit5->substring2,hit3->substringA)) {
      return Substring_insert_length(hit5->substring2,hit3->substringA);
    } else if (hit5->substringD != NULL && Substring_overlap_p(hit5->substringD,hit3->substringA)) {
      return Substring_insert_length(hit5->substringD,hit3->substringA);
    } else if (hit5->substringA != NULL && Substring_overlap_p(hit5->substringA,hit3->substringA)) {
      return Substring_insert_length(hit5->substringA,hit3->substringA);
    }
  }

  return 0;
}
  



Stage3pair_T
Stage3pair_new (T hit5, T hit3, char *query5, char *query3,
		Genomicpos_T *splicesites, Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		UINT4 *genome_blocks, UINT4 *snp_blocks,
		Genomicpos_T expected_pairlength, Genomicpos_T pairlength_deviation,
		Pairtype_T pairtype, int localsplicing_penalty, bool dibasep, bool cmetp) {
  Stage3pair_T new = (Stage3pair_T) MALLOC(sizeof(*new));
  Genomicpos_T insertlength, genomicstart, genomicend;
  int nbingo, bingoi5, bingoi3, nbounded, boundedi5, boundedi3, i, j;
  bool new5p = false, new3p = false;
  bool copy5p = true, copy3p = true;

  Substring_T donor, acceptor, shortexon;
  Genomicpos_T segment_left;
  int nmismatches_shortend;
  int donor_knowni, acceptor_knowni;
  int splice_pos;
  int ignore_found_score = 0;

  int querylength5 = hit5->querylength_adj;
  int querylength3 = hit3->querylength_adj;

  new->pairtype = pairtype;

#if 0
  new->mapq_loglik = hit5->mapq_loglik + hit3->mapq_loglik;
  new->mapq_score = 0;
#endif

  if (hit5->chrnum == 0 || hit3->chrnum == 0) {
    new->dir = 0;
    new->insertlength = pair_insert_length(hit5,hit3);

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
    /* Try to resolve ambiguity on inside of paired ends */

    /* Concordant directions on same chromosome (plus) */
    debug9(printf("Concordant on plus strand\n"));
    debug9(printf("hit5 %s ambiguous %d,%d and hit3 %s ambiguous %d,%d\n",
		  hittype_string(hit5->hittype),hit5->start_ambiguous_p,hit5->end_ambiguous_p,
		  hittype_string(hit3->hittype),hit3->start_ambiguous_p,hit3->end_ambiguous_p));
    new->dir = +1;

    if (hit5->end_ambiguous_p == true && hit3->start_ambiguous_p == true) {
      debug9(printf("Got ambiguous at 5' and ambiguous at 3':"));
      nbounded = nbingo = 0;
      for (i = 0; i < hit5->end_nambi; i++) {
	genomicend = splicesites[hit5->end_ambi[i]];
	for (j = 0; j < hit3->start_nambi; j++) {
	  genomicstart = splicesites[hit3->start_ambi[j]];
	  debug9(printf(" %u,%u",genomicend,genomicstart));
	  if (genomicend < genomicstart) {
	    nbounded++;
	    boundedi5 = i;
	    boundedi3 = j;
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
	  }
	}
      }
      if (nbingo == 1) {
	new5p = true; new3p = true;
      } else if (nbounded == 1) {
	new5p = true; new3p = true; bingoi5 = boundedi5; bingoi3 = boundedi3;
      }
      debug9(printf("\n"));

    } else if (hit5->end_ambiguous_p == true) {
      debug9(printf("Got ambiguous at 5' (%s):",hittype_string(hit5->hittype)));
      nbounded = nbingo = 0;
      for (i = 0; i < hit5->end_nambi; i++) {
	genomicend = splicesites[hit5->end_ambi[i]];
	debug9(printf(" %u",genomicend));
	if (genomicend < hit3->genomicstart) {
	  nbounded++;
	  boundedi5 = i;
	  insertlength = hit3->genomicstart - genomicend + querylength5 + querylength3;
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
	}
      }
      if (nbingo == 1) {
	new5p = true;
      } else if (nbounded == 1) {
	new5p = true; bingoi5 = boundedi5;
      }
      debug9(printf("\n"));

    } else if (hit3->start_ambiguous_p == true) {
      debug9(printf("Got ambiguous at 3':"));
      nbounded = nbingo = 0;
      for (j = 0; j < hit3->start_nambi; j++) {
	genomicstart = splicesites[hit3->start_ambi[j]];
	debug9(printf(" %u",genomicstart));
	if (hit5->genomicend < genomicstart) {
	  nbounded++;
	  boundedi3 = j;
	  insertlength = genomicstart - hit5->genomicend + querylength5 + querylength3;
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
	}
      }
      if (nbingo == 1) {
	new3p = true;
      } else if (nbounded == 1) {
	new3p = true; bingoi3 = boundedi3;
      }
      debug9(printf("\n"));
    }

    if (new5p == false) {
      /* Skip */
    } else if (hit5->hittype == ONE_THIRD_SHORTEXON || hit5->hittype == TWO_THIRDS_SHORTEXON) {

      if (hit5->sensedir == SENSE_FORWARD) {
	/* End 1 */
	shortexon = hit5->substring1;
	
	donor_knowni = Substring_splicesites_i_D(shortexon);
	splice_pos = Substring_chimera_pos_D(shortexon);
	acceptor_knowni = hit5->end_ambi[bingoi5];
	nmismatches_shortend = hit5->end_amb_nmismatches[bingoi5];
	segment_left = splicesites[acceptor_knowni] - splice_pos;
	
	if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
					       /*ncolordiffs*/0,/*prob*/2.0,/*left*/segment_left,query5_compress_fwd,
					       genome_blocks,snp_blocks,querylength5,/*plusp*/true,/*sensep*/true,
					       query5,Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					       dibasep,cmetp)) != NULL) {
	  debug9(printf("Resolved shortexon, End 1: Splice from donor #%d to acceptor #%d\n",
			donor_knowni,acceptor_knowni));
	  hit5 = Stage3_new_shortexon(&ignore_found_score,/*donor*/hit5->substringD,acceptor,shortexon,
				      /*acceptor_distance*/hit5->acceptor_distance,
				      /*donor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				      /*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
				      localsplicing_penalty,querylength5,/*first_read_p*/true,/*sensedir*/SENSE_FORWARD);
	  copy5p = false;
	}

      } else if (hit5->sensedir == SENSE_ANTI) {
	/* End 6 */
	shortexon = hit5->substring1;

	acceptor_knowni = Substring_splicesites_i_A(shortexon);
	splice_pos = Substring_chimera_pos_A(shortexon);
	donor_knowni = hit5->end_ambi[bingoi5];
	nmismatches_shortend = hit5->end_amb_nmismatches[bingoi5];
	segment_left = splicesites[donor_knowni] - splice_pos;

	if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
					 /*ncolordiffs*/0,/*prob*/2.0,/*left*/segment_left,query5_compress_fwd,
					 genome_blocks,snp_blocks,querylength5,/*plusp*/true,/*sensep*/false,
					 query5,Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					 dibasep,cmetp)) != NULL) {
	  debug9(printf("Resolved shortexon, End 6: Splice from antiacceptor #%d to antidonor #%d\n",
			acceptor_knowni,donor_knowni));
	  hit5 = Stage3_new_shortexon(&ignore_found_score,donor,/*acceptor*/hit5->substringA,shortexon,
				      /*acceptor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				      /*donor_distance*/hit5->donor_distance,
				      /*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
				      localsplicing_penalty,querylength5,/*first_read_p*/true,/*sensedir*/SENSE_ANTI);
	  copy5p = false;
	}

      } else {
	fprintf(stderr,"Shortexon hit5 has no sensedir\n");
	abort();
      }


    } else if (hit5->hittype == HALFSPLICE_DONOR) {
      /* End 1 */
      assert(hit5->sensedir == SENSE_FORWARD);
      donor = hit5->substring1;

      donor_knowni = Substring_splicesites_i(donor);
      splice_pos = Substring_chimera_pos(donor);
      acceptor_knowni = hit5->end_ambi[bingoi5];
      nmismatches_shortend = hit5->end_amb_nmismatches[bingoi5];
      segment_left = splicesites[acceptor_knowni] - splice_pos;

      if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
					     /*ncolordiffs*/0,/*prob*/2.0,/*left*/segment_left,query5_compress_fwd,
					     genome_blocks,snp_blocks,querylength5,
					     /*plusp*/true,/*sensep*/true,
					     query5,Substring_chrnum(donor),Substring_chroffset(donor),
					     dibasep,cmetp)) != NULL) {
	debug9(printf("Resolved halfsplice_donor, End 1: Splice from donor #%d to acceptor #%d\n",
		      Substring_splicesites_i(donor),Substring_splicesites_i(acceptor)));
	hit5 = Stage3_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				 donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				 /*shortdistancep*/true,localsplicing_penalty,querylength5,
				 /*ambi_left*/NULL,/*ambi_right*/NULL,
				 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				 /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,/*sensedir*/SENSE_FORWARD);
	copy5p = false;
      }
	  
    } else if (hit5->hittype == HALFSPLICE_ACCEPTOR) {
      /* End 6 */
      assert(hit5->sensedir == SENSE_ANTI);
      acceptor = hit5->substring1;

      acceptor_knowni = Substring_splicesites_i(acceptor);
      splice_pos = Substring_chimera_pos(acceptor);
      donor_knowni = hit5->end_ambi[bingoi5];
      nmismatches_shortend = hit5->end_amb_nmismatches[bingoi5];
      segment_left = splicesites[donor_knowni] - splice_pos;

      if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
				       /*ncolordiffs*/0,/*prob*/2.0,/*left*/segment_left,query5_compress_fwd,
				       genome_blocks,snp_blocks,querylength5,
				       /*plusp*/true,/*sensep*/false,
				       query5,Substring_chrnum(acceptor),Substring_chroffset(acceptor),
				       dibasep,cmetp)) != NULL) {
	debug9(printf("Resolved halfsplice_acceptor, End 6: Splice from antiacceptor #%d to antidonor #%d\n",
		      Substring_splicesites_i(acceptor),Substring_splicesites_i(donor)));
	hit5 = Stage3_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				 donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				 /*shortdistancep*/true,localsplicing_penalty,querylength5,
				 /*ambi_left*/NULL,/*ambi_right*/NULL,
				 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				 /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,/*sensedir*/SENSE_ANTI);
	copy5p = false;
      }

    } else {
      fprintf(stderr,"Unexpected hittype %d for ambiguous end\n",hit5->hittype);
      abort();
    }

    if (new3p == false) {
      /* Skip */

    } else if (hit3->hittype == ONE_THIRD_SHORTEXON || hit3->hittype == TWO_THIRDS_SHORTEXON) {

      if (hit3->sensedir == SENSE_ANTI) {
	/* End 5 */
	shortexon = hit3->substring1;

	donor_knowni = Substring_splicesites_i_D(shortexon);
	splice_pos = Substring_chimera_pos_D(shortexon);
	acceptor_knowni = hit3->start_ambi[bingoi3];
	nmismatches_shortend = hit3->start_amb_nmismatches[bingoi3];
	segment_left = splicesites[acceptor_knowni] - splice_pos;

	if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
					       /*ncolordiffs*/0,/*prob*/2.0,segment_left,query3_compress_fwd,
					       genome_blocks,snp_blocks,querylength3,/*plusp*/true,/*sensep*/false,
					       query3,Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					       dibasep,cmetp)) != NULL) {
	  debug9(printf("Resolved shortexonr, End 5: Splice from antidonor #%d to antiacceptor #%d\n",
			donor_knowni,acceptor_knowni));
	  hit3 = Stage3_new_shortexon(&ignore_found_score,/*donor*/hit3->substringD,acceptor,shortexon,
				      /*acceptor_distance*/hit3->acceptor_distance,
				      /*donor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				      /*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
				      localsplicing_penalty,querylength3,/*first_read_p*/false,/*sensedir*/SENSE_ANTI);
	  copy3p = false;
	}

      } else if (hit3->sensedir == SENSE_FORWARD) {
	/* End 2 */
	shortexon = hit3->substring1;

	acceptor_knowni = Substring_splicesites_i_A(shortexon);
	splice_pos = Substring_chimera_pos_A(shortexon);
	donor_knowni = hit3->start_ambi[bingoi3];
	nmismatches_shortend = hit3->start_amb_nmismatches[bingoi3];
	segment_left = splicesites[donor_knowni] - splice_pos;

	if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
					 /*ncolordiffs*/0,/*prob*/2.0,segment_left,query3_compress_fwd,
					 genome_blocks,snp_blocks,querylength3,/*plusp*/true,/*sensep*/true,
					 query3,Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					 dibasep,cmetp)) != NULL) {
	  debug9(printf("Resolved shortexon, End 2: Splice from acceptor #%d to donor #%d\n",
			acceptor_knowni,donor_knowni));
	  hit3 = Stage3_new_shortexon(&ignore_found_score,donor,/*acceptor*/hit3->substringA,shortexon,
				      /*acceptor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				      /*donor_distance*/hit3->donor_distance,
				      /*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
				      localsplicing_penalty,querylength3,/*first_read_p*/false,/*sensedir*/SENSE_FORWARD);
	  copy3p = false;
	}

      } else {
	fprintf(stderr,"Shortexon hit5 has no sensedir\n");
	abort();
      }


    } else if (hit3->hittype == HALFSPLICE_DONOR) {
      /* End 5 */
      assert(hit3->sensedir == SENSE_ANTI);
      donor = hit3->substring1;

      donor_knowni = Substring_splicesites_i(donor);
      splice_pos = Substring_chimera_pos(donor);
      acceptor_knowni = hit3->start_ambi[bingoi3];
      nmismatches_shortend = hit3->start_amb_nmismatches[bingoi3];
      segment_left = splicesites[acceptor_knowni] - splice_pos;

      if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
					     /*ncolordiffs*/0,/*prob*/2.0,segment_left,query3_compress_fwd,genome_blocks,snp_blocks,
					     querylength3,/*plusp*/true,/*sensep*/false,
					     query3,Substring_chrnum(donor),Substring_chroffset(donor),
					     dibasep,cmetp)) != NULL) {
	debug9(printf("Resolved halfsplice donor, End 5: Splice from antidonor #%d to antiacceptor #%d\n",
		      Substring_splicesites_i(donor),Substring_splicesites_i(acceptor)));
	hit3 = Stage3_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				 donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				 /*shortdistancep*/true,localsplicing_penalty,querylength3,
				 /*ambi_left*/NULL,/*ambi_right*/NULL,
				 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				 /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				 /*sensedir*/SENSE_ANTI);
	copy3p = false;
      }

    } else if (hit3->hittype == HALFSPLICE_ACCEPTOR) {
      /* End 2 */
      assert(hit3->sensedir == SENSE_FORWARD);
      acceptor = hit3->substring1;

      acceptor_knowni = Substring_splicesites_i(acceptor);
      splice_pos = Substring_chimera_pos(acceptor);
      donor_knowni = hit3->start_ambi[bingoi3];
      nmismatches_shortend = hit3->start_amb_nmismatches[bingoi3];
      segment_left = splicesites[donor_knowni] - splice_pos;

      if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,splice_pos,nmismatches_shortend,
				       /*ncolordiffs*/0,/*prob*/2.0,segment_left,query3_compress_fwd,genome_blocks,snp_blocks,
				       querylength3,/*plusp*/true,/*sensep*/true,
				       query3,Substring_chrnum(acceptor),Substring_chroffset(acceptor),
				       dibasep,cmetp)) != NULL) {
	debug9(printf("Resolved halfsplice acceptor, End 2: Splice from acceptor #%d (%u) to donor #%d (%u)\n",
		      Substring_splicesites_i(acceptor),splicesites[acceptor_knowni],
		      Substring_splicesites_i(donor),splicesites[donor_knowni]));
	hit3 = Stage3_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				 donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				 /*shortdistancep*/true,localsplicing_penalty,querylength3,
				 /*ambi_left*/NULL,/*ambi_right*/NULL,
				 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				 /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				 /*sensedir*/SENSE_FORWARD);
	copy3p = false;
      }

    } else {
      fprintf(stderr,"Unexpected hittype %d for ambiguous end\n",hit3->hittype);
      abort();
    }
      

    /* Have 5-start..end and 3-start..end */
    if (hit5->genomicend < hit3->genomicstart) {
      /* No overlap */
      new->insertlength = (hit3->genomicstart - hit5->genomicend) + querylength5 + querylength3;
      debug4(printf("plus, no overlap: insert length %d = start3 %u - end5 %u + %d + %d\n",
		    new->insertlength,hit3->genomicstart,hit5->genomicend,querylength5,querylength3));
    } else {
      new->insertlength = pair_insert_length(hit5,hit3);
    }


  } else {
    /* Try to resolve ambiguity on inside of paired ends */

    /* Concordant directions on same chromosome (minus) */
    debug9(printf("Concordant on minus strand\n"));
    debug9(printf("hit5 %s ambiguous %d,%d and hit3 %s ambiguous %d,%d\n",
		  hittype_string(hit5->hittype),hit5->start_ambiguous_p,hit5->end_ambiguous_p,
		  hittype_string(hit3->hittype),hit3->start_ambiguous_p,hit3->end_ambiguous_p));
    new->dir = -1;

    if (hit5->end_ambiguous_p == true && hit3->start_ambiguous_p == true) {
      debug9(printf("Got ambiguous at 5' and ambiguous at 3':"));
      nbounded = nbingo = 0;
      for (i = 0; i < hit5->end_nambi; i++) {
	genomicend = splicesites[hit5->end_ambi[i]];
	for (j = 0; j < hit3->start_nambi; j++) {
	  genomicstart = splicesites[hit3->start_ambi[j]];
	  debug9(printf(" %u,%u",genomicend,genomicstart));
	  if (genomicstart < genomicend) {
	    nbounded++;
	    boundedi5 = i;
	    boundedi3 = j;
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
	  }
	}
      }
      if (nbingo == 1) {
	new5p = true; new3p = true;
      } else if (nbounded == 1) {
	new5p = true; new3p = true; bingoi5 = boundedi5; bingoi3 = boundedi3;
      }
      debug9(printf("\n"));

    } else if (hit5->end_ambiguous_p == true) {
      debug9(printf("Got ambiguous at 5':"));
      nbounded = nbingo = 0;
      for (i = 0; i < hit5->end_nambi; i++) {
	genomicend = splicesites[hit5->end_ambi[i]];
	debug9(printf(" %u",genomicend));
	if (hit3->genomicstart < genomicend) {
	  nbounded++;
	  boundedi5 = i;
	  boundedi3 = j;
	  insertlength = genomicend - hit3->genomicstart + querylength5 + querylength3;
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
	}
      }
      if (nbingo == 1) {
	new5p = true;
      } else if (nbounded == 1) {
	new5p = true; bingoi5 = boundedi5;
      }
      debug9(printf("\n"));

    } else if (hit3->start_ambiguous_p == true) {
      debug9(printf("Got ambiguous at 3':"));
      nbounded = nbingo = 0;
      for (j = 0; j < hit3->start_nambi; j++) {
	genomicstart = splicesites[hit3->start_ambi[j]];
	debug9(printf(" %u",genomicstart));
	if (genomicstart < hit5->genomicend) {
	  nbounded++;
	  boundedi3 = j;
	  insertlength = hit5->genomicend - genomicstart + querylength5 + querylength3;
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
	}
      }
      if (nbingo == 1) {
	new3p = true;
      } else if (nbounded == 1) {
	new3p = true; bingoi3 = boundedi3;
      }
      debug9(printf("\n"));
    }

    if (new5p == false) {
      /* Skip */

    } else if (hit5->hittype == ONE_THIRD_SHORTEXON || hit5->hittype == TWO_THIRDS_SHORTEXON) {
      if (hit5->sensedir == SENSE_FORWARD) {
	/* End 3 */
	shortexon = hit5->substring1;

	donor_knowni = Substring_splicesites_i_D(shortexon);
	splice_pos = Substring_chimera_pos_D(shortexon);
	acceptor_knowni = hit5->end_ambi[bingoi5];
	nmismatches_shortend = hit5->end_amb_nmismatches[bingoi5];
	segment_left = splicesites[acceptor_knowni] - (querylength5 - splice_pos);

	if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,querylength5 - splice_pos,
					       nmismatches_shortend,/*ncolordiffs*/0,
					       /*prob*/2.0,segment_left,query5_compress_rev,
					       genome_blocks,snp_blocks,querylength5,/*plusp*/false,/*sensep*/true,
					       query5,Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					       dibasep,cmetp)) != NULL) {
	  debug9(printf("Resolved shortexon, End 3: Splice from donor #%d to acceptor #%d\n",
			donor_knowni,acceptor_knowni));
	  hit5 = Stage3_new_shortexon(&ignore_found_score,/*donor*/hit5->substringD,acceptor,shortexon,
				      /*acceptor_distance*/hit5->acceptor_distance,
				      /*donor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				      /*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
				      localsplicing_penalty,querylength5,/*first_read_p*/true,
				      /*sensedir*/SENSE_FORWARD);
	  copy5p = false;
	}

      } else if (hit5->sensedir == SENSE_ANTI) {
	/* End 8 */
	shortexon = hit5->substring1;

	acceptor_knowni = Substring_splicesites_i_A(shortexon);
	splice_pos = Substring_chimera_pos_A(shortexon);
	donor_knowni = hit5->end_ambi[bingoi5];
	nmismatches_shortend = hit5->end_amb_nmismatches[bingoi5];
	segment_left = splicesites[donor_knowni] - (querylength5 - splice_pos);

	if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,querylength5 - splice_pos,
					 nmismatches_shortend,/*ncolordiffs*/0,
					 /*prob*/2.0,segment_left,query5_compress_rev,
					 genome_blocks,snp_blocks,querylength5,/*plusp*/false,/*sensep*/false,
					 query5,Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					 dibasep,cmetp)) != NULL) {
	  debug9(printf("Resolved shortexon, End 8: Splice from antiacceptor #%d to antidonor #%d\n",
			acceptor_knowni,donor_knowni));
	  hit5 = Stage3_new_shortexon(&ignore_found_score,donor,/*acceptor*/hit5->substringA,shortexon,
				      /*acceptor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				      /*donor_distance*/hit5->donor_distance,
				      /*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
				      localsplicing_penalty,querylength5,/*first_read_p*/true,
				      /*sensedir*/SENSE_ANTI);
	  copy5p = false;
	}
	
      } else {
	fprintf(stderr,"Shortexon hit5 has no sensedir\n");
	abort();
      }

    } else if (hit5->hittype == HALFSPLICE_DONOR) {
      /* End 3 */
      assert(hit5->sensedir == SENSE_FORWARD);
      donor = hit5->substring1;

      donor_knowni = Substring_splicesites_i(donor);
      splice_pos = Substring_chimera_pos(donor);
      acceptor_knowni = hit5->end_ambi[bingoi5];
      nmismatches_shortend = hit5->end_amb_nmismatches[bingoi5];
      segment_left = splicesites[acceptor_knowni] - (querylength5 - splice_pos);

      if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,querylength5 - splice_pos,
					     nmismatches_shortend,/*ncolordiffs*/0,
					     /*prob*/2.0,segment_left,query5_compress_rev,genome_blocks,snp_blocks,
					     querylength5,/*plusp*/false,/*sensep*/true,
					     query5,Substring_chrnum(donor),Substring_chroffset(donor),
					     dibasep,cmetp)) != NULL) {
	debug9(printf("Resolved halfsplice, End 3: Splice from donor #%d to acceptor #%d\n",
		      Substring_splicesites_i(donor),Substring_splicesites_i(acceptor)));
	hit5 = Stage3_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				 donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				 /*shortdistancep*/true,localsplicing_penalty,querylength5,
				 /*ambi_left*/NULL,/*ambi_right*/NULL,
				 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				 /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/true,
				 /*sensedir*/SENSE_FORWARD);
	copy5p = false;
      }

    } else if (hit5->hittype == HALFSPLICE_ACCEPTOR) {
      /* End 8 */
      assert(hit5->sensedir == SENSE_ANTI);
      acceptor = hit5->substring1;

      acceptor_knowni = Substring_splicesites_i(acceptor);
      splice_pos = Substring_chimera_pos(acceptor);
      donor_knowni = hit5->end_ambi[bingoi5];
      nmismatches_shortend = hit5->end_amb_nmismatches[bingoi5];
      segment_left = splicesites[donor_knowni] - (querylength5 - splice_pos);

      if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,querylength5 - splice_pos,
				       nmismatches_shortend,/*ncolordiffs*/0,
				       /*prob*/2.0,segment_left,query5_compress_rev,genome_blocks,snp_blocks,
				       querylength5,/*plusp*/false,/*sensep*/false,
				       query5,Substring_chrnum(acceptor),Substring_chroffset(acceptor),
				       dibasep,cmetp)) != NULL) {
	debug9(printf("Resolved halfsplice acceptor, End 8: Splice from antiacceptor #%d to antidonor #%d\n",
		      Substring_splicesites_i(acceptor),Substring_splicesites_i(donor)));
	hit5 = Stage3_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				 donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				 /*shortdistancep*/true,localsplicing_penalty,querylength5,
				 /*ambi_left*/NULL,/*ambi_right*/NULL,
				 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				 /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/true,
				 /*sensedir*/SENSE_ANTI);
	copy5p = false;
      }

    } else {
      fprintf(stderr,"Unexpected hittype %d for ambiguous end\n",hit5->hittype);
      abort();
    }


    if (new3p == false) {
      /* Skip */

    } else if (hit3->hittype == ONE_THIRD_SHORTEXON || hit3->hittype == TWO_THIRDS_SHORTEXON) {
      if (hit3->sensedir == SENSE_ANTI) {
	/* End 7 */
	shortexon = hit3->substring1;

	donor_knowni = Substring_splicesites_i_D(shortexon);
	splice_pos = Substring_chimera_pos_D(shortexon);
	acceptor_knowni = hit3->start_ambi[bingoi3];
	nmismatches_shortend = hit3->start_amb_nmismatches[bingoi3];
	segment_left = splicesites[acceptor_knowni] - (querylength3 - splice_pos);

	if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,querylength3 - splice_pos,
					       nmismatches_shortend,/*ncolordiffs*/0,
					       /*prob*/2.0,segment_left,query3_compress_rev,genome_blocks,snp_blocks,
					       querylength3,/*plusp*/false,/*sensep*/false,
					       query3,Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					       dibasep,cmetp)) != NULL) {
	  debug9(printf("Resolved shortexon, End 7: Splice from antidonor #%d to antiacceptor #%d\n",
			donor_knowni,acceptor_knowni));
	  hit3 = Stage3_new_shortexon(&ignore_found_score,/*donor*/hit3->substringD,acceptor,shortexon,
				      /*acceptor_distance*/hit3->acceptor_distance,
				      /*donor_distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				      /*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				      /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*copy_shortexon_p*/true,
				      localsplicing_penalty,querylength3,/*first_read_p*/false,
				      /*sensedir*/SENSE_ANTI);
	  copy3p = false;
	}

      } else if (hit3->sensedir == SENSE_FORWARD) {
	/* End 4 */
	shortexon = hit3->substring1;

	acceptor_knowni = Substring_splicesites_i_A(shortexon);
	splice_pos = Substring_chimera_pos_A(shortexon);
	donor_knowni = hit3->start_ambi[bingoi3];
	nmismatches_shortend = hit3->start_amb_nmismatches[bingoi3];
	segment_left = splicesites[donor_knowni] - (querylength3 - splice_pos);

	if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,querylength3 - splice_pos,
					 nmismatches_shortend,/*ncolordiffs*/0,
					 /*prob*/2.0,segment_left,query3_compress_rev,genome_blocks,snp_blocks,
					 querylength3,/*plusp*/false,/*sensep*/true,
					 query3,Substring_chrnum(shortexon),Substring_chroffset(shortexon),
					 dibasep,cmetp)) != NULL) {
	  debug9(printf("Resolved halfsplice_acceptor, End 4: Splice from acceptor #%d to #%d\n",
			acceptor_knowni,donor_knowni));
	  hit3 = Stage3_new_shortexon(&ignore_found_score,donor,/*acceptor*/hit3->substringA,shortexon,
				      /*acceptor_distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				      /*donor_distance*/hit3->donor_distance,
				      /*ambi_left*/NULL,/*ambi_right*/NULL,
				      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				      /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*copy_shortexon_p*/true,
				      localsplicing_penalty,querylength3,/*first_read_p*/false,
				      /*sensedir*/SENSE_FORWARD);
	  copy3p = false;
	}

      } else {
	fprintf(stderr,"Shortexon hit3 has no sensedir\n");
	abort();
      }

    } else if (hit3->hittype == HALFSPLICE_DONOR) {
      /* End 7 */
      assert(hit3->sensedir == SENSE_ANTI);
      donor = hit3->substring1;

      donor_knowni = Substring_splicesites_i(donor);
      splice_pos = Substring_chimera_pos(donor);
      acceptor_knowni = hit3->start_ambi[bingoi3];
      nmismatches_shortend = hit3->start_amb_nmismatches[bingoi3];
      segment_left = splicesites[acceptor_knowni] - (querylength3 - splice_pos);

      if ((acceptor = Substring_new_acceptor(acceptor_knowni,/*joffset*/0,querylength3 - splice_pos,
					     nmismatches_shortend,/*ncolordiffs*/0,
					     /*prob*/2.0,segment_left,query3_compress_rev,genome_blocks,snp_blocks,
					     querylength3,/*plusp*/false,/*sensep*/false,
					     query3,Substring_chrnum(donor),Substring_chroffset(donor),
					     dibasep,cmetp)) != NULL) {
	debug9(printf("Resolved halfsplice_donor, End 7: Splice from antidonor #%d to antiacceptor #%d\n",
		      Substring_splicesites_i(donor),Substring_splicesites_i(acceptor)));
	hit3 = Stage3_new_splice(&ignore_found_score,Substring_nmismatches_whole(donor),/*nmismatches_acceptor*/nmismatches_shortend,
				 donor,acceptor,/*distance*/splicesites[acceptor_knowni] - splicesites[donor_knowni],
				 /*shortdistancep*/true,localsplicing_penalty,querylength3,
				 /*ambi_left*/NULL,/*ambi_right*/NULL,
				 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				 /*copy_donor_p*/true,/*copy_acceptor_p*/false,/*first_read_p*/false,
				 /*sensedir*/SENSE_ANTI);
	copy3p = false;
      }

    } else if (hit3->hittype == HALFSPLICE_ACCEPTOR) {
      /* End 4 */
      assert(hit3->sensedir == SENSE_FORWARD);
      acceptor = hit3->substring1;

      acceptor_knowni = Substring_splicesites_i(acceptor);
      splice_pos = Substring_chimera_pos(acceptor);
      donor_knowni = hit3->start_ambi[bingoi3];
      nmismatches_shortend = hit3->start_amb_nmismatches[bingoi3];
      segment_left = splicesites[donor_knowni] - (querylength3 - splice_pos);

      if ((donor = Substring_new_donor(donor_knowni,/*joffset*/0,querylength3 - splice_pos,
				       nmismatches_shortend,/*ncolordiffs*/0,
				       /*prob*/2.0,segment_left,query3_compress_rev,genome_blocks,snp_blocks,
				       querylength3,/*plusp*/false,/*sensep*/true,
				       query3,Substring_chrnum(acceptor),Substring_chroffset(acceptor),
				       dibasep,cmetp)) != NULL) {
	debug9(printf("Resolved halfsplice_acceptor, End 4: Splice from acceptor #%d to #%d\n",
		      Substring_splicesites_i(acceptor),Substring_splicesites_i(donor)));
	hit3 = Stage3_new_splice(&ignore_found_score,/*nmismatches_donor*/nmismatches_shortend,Substring_nmismatches_whole(acceptor),
				 donor,acceptor,/*distance*/splicesites[donor_knowni] - splicesites[acceptor_knowni],
				 /*shortdistancep*/true,localsplicing_penalty,querylength3,
				 /*ambi_left*/NULL,/*ambi_right*/NULL,
				 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
				 /*copy_donor_p*/false,/*copy_acceptor_p*/true,/*first_read_p*/false,
				 /*sensedir*/SENSE_FORWARD);
	copy3p = false;
      }

    } else {
      fprintf(stderr,"Unexpected hittype %d for ambiguous end\n",hit3->hittype);
      abort();
    }


    /* Have 3-end..start and 5-end..start */
    if (hit3->genomicstart < hit5->genomicend) {
      /* No overlap */
      new->insertlength = (hit5->genomicend - hit3->genomicstart) + querylength5 + querylength3;
      debug4(printf("minus, no overlap: insert length %d = end5 %u - start3 %u + %d + %d\n",
		    new->insertlength,hit5->genomicend,hit3->genomicstart,querylength5,querylength3));
      
    } else {
      new->insertlength = pair_insert_length(hit5,hit3);
    }

  }


  if (new->insertlength == 0) {
    /* Not concordant */
    new->absdifflength_bingo_p = false;
    new->absdifflength = -1U;
  } else {
    if (new->insertlength < expected_pairlength) {
      new->absdifflength = expected_pairlength - new->insertlength;
    } else {
      new->absdifflength = new->insertlength - expected_pairlength;
    }
    if (new->absdifflength <= pairlength_deviation) {
      new->absdifflength_bingo_p = true;
    } else {
      new->absdifflength_bingo_p = false;
    }
  }

  if (SENSE_CONSISTENT_P(hit5->sensedir,hit3->sensedir)) {
    new->sense_consistent_p = true;
  } else {
    new->sense_consistent_p = false;
  }


  new->score = hit5->score + hit3->score;
  new->nmatches = hit5->nmatches + hit3->nmatches;

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

  if (copy5p == true) {
    new->hit5 = Stage3_copy(hit5);
  } else {
    new->hit5 = hit5;
  }

  if (copy3p == true) {
    new->hit3 = Stage3_copy(hit3);
  } else {
    new->hit3 = hit3;
  }

  return new;
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
Stage3_sort_bymatchdist (List_T hitlist, int maxchimerapaths) {
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
    Stage3_free(&(hits[i]));
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



/* Also sorts secondarily by score */
static int
hitpair_position_cmp (const void *a, const void *b) {
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
    return -1;
  } else if (y->high < x->high) {
    return +1;
  } else if (x->sense_consistent_p == true && y->sense_consistent_p == false) {
    return -1;
  } else if (x->sense_consistent_p == false && y->sense_consistent_p == true) {
    return +1;
  } else {
    return 0;
  }
}

static bool
hitpair_equal (Stage3pair_T x, Stage3pair_T y) {
  if (x->dir != y->dir) {
    return false;		/* Different strands */
  } else if (x->hit5->low != y->hit5->low) {
    return false;
  } else if (x->hit5->high != y->hit5->high) {
    return false;
  } else if (x->hit3->low != y->hit3->low) {
    return false;
  } else if (x->hit3->high != y->hit3->high) {
    return false;
  } else {
    return true;
  }
}

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


List_T
Stage3pair_remove_duplicates_exact (List_T hitpairlist) {
  List_T unique = NULL;
  Stage3pair_T hitpair, *hitpairs;
  int n, i, j;
  bool *eliminate;

  n = List_length(hitpairlist);
  debug8(printf("Entered Stage3pair_remove_duplicates with %d pairs\n",n));
  if (n == 0) {
    return NULL;
  } else {
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    List_free(&hitpairlist);
  }

  debug8(printf("Checking for exact duplicates\n"));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_position_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s): %u..%u (dir = %d), nmatches: %d\n",
		  i,pairtype_string(hitpair->pairtype),hitpair->low,hitpair->high,hitpair->dir,hitpair->nmatches);
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && hitpair_equal(hitpairs[j],hitpairs[i]) == true) {
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

  return unique;
}


List_T
Stage3pair_remove_duplicates (List_T hitpairlist) {
  List_T unique = NULL;
  Stage3pair_T hitpair, *hitpairs, *prev;
  int nkept, n, i, j, k, l;
  bool *eliminate;

  n = List_length(hitpairlist);
  debug8(printf("Entered Stage3pair_remove_duplicates with %d pairs\n",n));
  if (n == 0) {
    return NULL;
  } else {
    eliminate = (bool *) CALLOC(n,sizeof(bool));
    hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
    List_free(&hitpairlist);
  }

  /* Step 1.  Check for exact duplicates */
  debug8(printf("Checking for exact duplicates\n"));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_position_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s): %u..%u (dir = %d), nmatches: %d\n",
		  i,pairtype_string(hitpair->pairtype),hitpair->low,hitpair->high,hitpair->dir,hitpair->nmatches);
	 }
	 );

  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && hitpair_equal(hitpairs[j],hitpairs[i]) == true) {
      debug8(printf("  %d is identical to %d => eliminating\n",j,i));
      eliminate[j] = true;
      j++;
    }
    i = j;
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
      debug8(printf("  Keeping %u..%u, nmatches %d (dir = %d)\n",
		    hitpair->low,hitpair->high,hitpair->nmatches,hitpair->dir));
      hitpairs[j++] = hitpair;
    } else {
      debug8(printf("  Eliminating %u..%u, nmatches %d (dir = %d)\n",
		    hitpair->low,hitpair->high,hitpair->nmatches,hitpair->dir));
      Stage3pair_free(&hitpair);
    }
  }

  FREE(prev);
#if 0
  /* Re-use previously allocated */
  FREE(eliminate);
#endif

  /* Step 2: Check for subsumption */
  n = nkept;
  debug8(printf("Checking for duplicates among %d hitpairs by subsumption\n",n));

#if 0
  /* Re-use previously allocated */
  eliminate = (bool *) CALLOC(n,sizeof(bool));
#else
  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
#endif
  /* qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_position_cmp); -- No need since original order was kept */
  
  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s): %u-%u (dir = %d), score: %d, nmatches: %d, absdifflength: %d, outerlength: %u\n",
		  i,pairtype_string(hitpair->pairtype),hitpair->low,hitpair->high,hitpair->dir,hitpair->score,hitpair->nmatches,
		  hitpair->absdifflength,hitpair->outerlength);
	 }
	 );


  i = 0;
  while (i < n) {
    j = i;
    while (j+1 < n && hitpair_overlap(hitpairs[j+1],hitpairs[j]) == true) {
      j = j+1;
    }

    debug8(printf("Cluster from %d through %d\n",i,j));

    for (k = i; k < j; k++) {
      for (l = k+1; l <= j; l++) {
	debug8(printf("%d vs %d",k,l));
	if (hitpair_subsumption(hitpairs[k],hitpairs[l]) == false) {
	  /* Skip */

	} else if (hitpairs[k]->absdifflength_bingo_p == true &&
		   hitpairs[l]->absdifflength_bingo_p == false) {
	  debug8(printf(" => eliminate %d by absdifflength (bingo)",l));
	  eliminate[l] = true;
	} else if (hitpairs[l]->absdifflength_bingo_p == true &&
		   hitpairs[k]->absdifflength_bingo_p == false) {
	  debug8(printf(" => eliminate %d by absdifflength (bingo)",k));
	  eliminate[k] = true;

	} else if (hitpairs[k]->nmatches > hitpairs[l]->nmatches) {
	  debug8(printf(" => eliminate %d by nmatches",l));
	  eliminate[l] = true;
	} else if (hitpairs[l]->nmatches > hitpairs[k]->nmatches) {
	  debug8(printf(" => eliminate %d by nmatches",k));
	  eliminate[k] = true;
#if 0
	} else if (hitpairs[k]->score < hitpairs[l]->score) {
	  debug8(printf(" => eliminate %d by score",l));
	  eliminate[l] = true;
	} else if (hitpairs[l]->score < hitpairs[k]->score) {
	  debug8(printf(" => eliminate %d by score",k));
	  eliminate[k] = true;
#endif
	} else if (hitpairs[k]->absdifflength < hitpairs[l]->absdifflength) {
	  debug8(printf(" => eliminate %d by absdifflength",l));
	  eliminate[l] = true;
	} else if (hitpairs[l]->absdifflength < hitpairs[k]->absdifflength) {
	  debug8(printf(" => eliminate %d by absdifflength",k));
	  eliminate[k] = true;
	} else if (hitpairs[k]->outerlength < hitpairs[l]->outerlength) {
	  debug8(printf(" => eliminate %d by outerlength",l));
	  eliminate[l] = true;
	} else if (hitpairs[l]->outerlength < hitpairs[k]->outerlength) {
	  debug8(printf(" => eliminate %d by outerlength",k));
	  eliminate[k] = true;
	}
	debug8(printf("\n"));
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
    eliminate[0] = false;
    nkept = 1;
  }

  prev = hitpairs;
  hitpairs = (Stage3pair_T *) CALLOC(nkept,sizeof(Stage3pair_T));

  for (i = 0, j = 0; i < n; i++) {
    hitpair = prev[i];
    if (eliminate[i] == false) {
      debug8(printf("  Keeping %u..%u, nmatches %d (dir = %d)\n",
		    hitpair->low,hitpair->high,hitpair->nmatches,hitpair->dir));
      hitpairs[j++] = hitpair;
    } else {
      debug8(printf("  Eliminating %u..%u, nmatches %d (dir = %d)\n",
		    hitpair->low,hitpair->high,hitpair->nmatches,hitpair->dir));
      Stage3pair_free(&hitpair);
    }
  }

  FREE(prev);
#if 0
  /* Re-use previously allocated */
  FREE(eliminate);
#endif

  /* Step 3: Check for identity */
  n = nkept;
  debug8(printf("Checking for duplicates among %d hitpairs by identity\n",n));

#if 0
  /* Re-use previously allocated */
  eliminate = (bool *) CALLOC(n,sizeof(bool));
#else
  for (i = 0; i < n; i++) {
    eliminate[i] = false;
  }
#endif
  /* qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_position_cmp); -- No need since original order was kept */

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d (%s): %u-%u (dir = %d), score: %d, nmatches: %d, absdifflength: %d, outerlength: %u\n",
		  i,pairtype_string(hitpair->pairtype),hitpair->low,hitpair->high,hitpair->dir,hitpair->score,hitpair->nmatches,
		  hitpair->absdifflength,hitpair->outerlength);
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
    while (j < n && hitpair_equal(hitpairs[j],hitpairs[i]) == true) {
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

  return unique;
}


void
Stage3pair_eval (Stage3pair_T *stage3pairarray, int npaths, int maxpaths,
		 Shortread_T queryseq5, Shortread_T queryseq3,
		 Compress_T query5_compress_fwd, Compress_T query5_compress_rev, 
		 Compress_T query3_compress_fwd, Compress_T query3_compress_rev, 
		 UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome,
		 char *quality_string_5, char *quality_string_3,
		 bool dibasep, bool cmetp) {
  char *query5, *query3;
  double maxlik, total, q;
  int mapq_score;
  int i;

  if (npaths == 0) {
    /* Skip */

  } else if (npaths == 1) {
    stage3pairarray[0]->mapq_score = MAPQ_max_quality_score(quality_string_5,Shortread_fulllength(queryseq5));
    if ((mapq_score = MAPQ_max_quality_score(quality_string_3,Shortread_fulllength(queryseq3))) > stage3pairarray[0]->mapq_score) {
      stage3pairarray[0]->mapq_score = mapq_score;
    }
    query5 = Shortread_fullpointer_uc(queryseq5);
    query3 = Shortread_fullpointer_uc(queryseq3);

    Stage3_display_prep(stage3pairarray[0]->hit5,query5,query5_compress_fwd,query5_compress_rev,
			genome_blocks,snp_blocks,genome,dibasep,cmetp);
    Stage3_display_prep(stage3pairarray[0]->hit3,query3,query3_compress_fwd,query3_compress_rev,
			genome_blocks,snp_blocks,genome,dibasep,cmetp);

  } else {
    /* Sort by insert length */
    qsort(stage3pairarray,npaths,sizeof(Stage3pair_T),Stage3pair_output_cmp);

    /* Compute mapq_loglik */
    query5 = Shortread_fullpointer_uc(queryseq5);
    query3 = Shortread_fullpointer_uc(queryseq3);
    for (i = 0; i < npaths; i++) {
      stage3pairarray[i]->mapq_loglik =
	Stage3_compute_mapq(stage3pairarray[i]->hit5,query5,query5_compress_fwd,query5_compress_rev,
			    genome_blocks,snp_blocks,quality_string_5,dibasep,cmetp);
      stage3pairarray[i]->mapq_loglik +=
      Stage3_compute_mapq(stage3pairarray[i]->hit3,query3,query3_compress_fwd,query3_compress_rev,
			  genome_blocks,snp_blocks,quality_string_3,dibasep,cmetp);

    }

    /* Find maximum likelihood */
    maxlik = stage3pairarray[0]->mapq_loglik;
    for (i = 1; i < npaths; i++) {
      if (stage3pairarray[i]->mapq_loglik > maxlik) {
	maxlik = stage3pairarray[i]->mapq_loglik;
      }
    }

    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < npaths; i++) {
      stage3pairarray[i]->mapq_loglik -= maxlik;
    }

    total = 0.0;
    for (i = 0; i < npaths; i++) {
      total += (stage3pairarray[i]->mapq_loglik = exp(stage3pairarray[i]->mapq_loglik));
    }

    /* Save on computation if possible */
    if (maxpaths < npaths) npaths = maxpaths;

    /* Prepare for display */
    query5 = Shortread_fullpointer_uc(queryseq5);
    query3 = Shortread_fullpointer_uc(queryseq3);
    for (i = 0; i < npaths; i++) {
      Stage3_display_prep(stage3pairarray[i]->hit5,query5,query5_compress_fwd,query5_compress_rev,
			  genome_blocks,snp_blocks,genome,dibasep,cmetp);
      Stage3_display_prep(stage3pairarray[i]->hit3,query3,query3_compress_fwd,query3_compress_rev,
			  genome_blocks,snp_blocks,genome,dibasep,cmetp);
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < npaths; i++) {
      stage3pairarray[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < npaths; i++) {
      if ((q = 1.0 - stage3pairarray[i]->mapq_loglik) < 1e-40) {
	stage3pairarray[i]->mapq_score = 40;
      } else {
	stage3pairarray[i]->mapq_score = rint(-10.0 * log10(q));
      }
    }

  }

  return;
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



List_T
Stage3pair_optimal_score (List_T hitpairlist, int cutoff_level, int suboptimal_mismatches) {
  List_T optimal = NULL, p;
  Stage3pair_T hitpair;
  int n;
  int minscore = MAX_QUERYLENGTH + MAX_QUERYLENGTH;
  int minscore_bingo = MAX_QUERYLENGTH + MAX_QUERYLENGTH;
  bool non_translocation_p = false;

  n = List_length(hitpairlist);
  if (n == 0) {
    return NULL;
  }

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->hit5->chrnum != 0 && hitpair->hit3->chrnum != 0) {
      non_translocation_p = true;
    }
  }

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    debug6(printf("%u..%u types %s and %s, score %d, pairlength %d, outerlength %u\n",
		  hitpair->low,hitpair->high,hittype_string(hitpair->hit5->hittype),hittype_string(hitpair->hit3->hittype),
		  hitpair->score,hitpair->insertlength,hitpair->outerlength));
    if ((hitpair->hit5->chrnum == 0 || hitpair->hit3->chrnum == 0) && non_translocation_p == true) {
      /* Skip, since we will eliminate anyway */
    } else if (hitpair->hit5->hittype == TERMINAL) {
      /* Skip from setting minscore */
    } else if (hitpair->hit3->hittype == TERMINAL) {
      /* Skip from setting minscore */
    } else if (hitpair->score <= cutoff_level) {
      if (hitpair->score < minscore) {
	minscore = hitpair->score;
      }
      if (hitpair->absdifflength_bingo_p == true) {
	if (hitpair->score < minscore_bingo) {
	  minscore_bingo = hitpair->score;
	}
      }
    }
  }

  debug6(printf("Stage3pair_optimal_score over %d pairs: minscore = %d + subopt:%d, minscore_bingo = %d\n",
	       n,minscore,suboptimal_mismatches,minscore_bingo));
  minscore += suboptimal_mismatches;

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;

    if ((hitpair->hit5->chrnum == 0 || hitpair->hit3->chrnum == 0) && non_translocation_p == true) {
      debug6(printf("Eliminating a hit pair at %u..%u with splice translocation\n",hitpair->low,hitpair->high));
      Stage3pair_free(&hitpair);
    } else if (hitpair->score > cutoff_level) {
      debug6(printf("Eliminating a hit pair at %u..%u with score %d > cutoff_level %d\n",
		   hitpair->low,hitpair->high,hitpair->score,cutoff_level));
      Stage3pair_free(&hitpair);
    } else if (hitpair->score >= minscore_bingo && hitpair->absdifflength_bingo_p == false) {
      debug6(printf("Eliminating a non-bingo hit pair at %u..%u with score %d and pairlength %u\n",
		   hitpair->low,hitpair->high,hitpair->score,hitpair->insertlength));
      Stage3pair_free(&hitpair);
    } else if (hitpair->score <= minscore) {
      debug6(printf("Keeping a hit pair with score %d\n",hitpair->score));
      optimal = List_push(optimal,hitpair);
    } else {
      debug6(printf("Eliminating a hit pair with score %d\n",hitpair->score));
      Stage3pair_free(&hitpair);
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
			   List_T *samechr, List_T hitpairs,
			   List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			   int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			   Genomicpos_T *splicesites, char *query5, char *query3,
			   Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			   Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			   UINT4 *genome_blocks, UINT4 *snp_blocks, Genomicpos_T pairmax,
			   Genomicpos_T expected_pairlength, Genomicpos_T pairlength_deviation,
			   int querylength5, int querylength3, int maxpairedpaths,
			   bool allow_concordant_translocations_p,
			   int splicing_penalty, bool dibasep, bool cmetp) {
#if 0
  /* No need for copies anymore, because effective_chrnum determined by inner substrings */
  List_T copies = NULL;
#endif
  List_T q, prev_start;
  Stage3pair_T stage3pair;
  T **hits5_plus, **hits5_minus, **hits3_plus, **hits3_minus, *hits5, *hits3, hit5, hit3;
  int *nhits5_plus, *nhits5_minus, *nhits3_plus, *nhits3_minus, nhits5, nhits3;
  int pairscore, score5, score3, i, j;
  bool *sorted5p, *sorted3p;
  Genomicpos_T insert_start;
  
  debug5(printf("Starting Stage3_pair_up_concordant with %d concordant, narray5 %d, narray3 %d, found_score %d, allow_concordant_translocations %d\n",
		*nconcordant,narray5,narray3,*found_score,allow_concordant_translocations_p));

  /* Find sizes for allocating memory */
  debug5(printf("Sizes of 5-end pieces by score level:\n"));
  nhits5_plus = (int *) CALLOC(cutoff_level_5+1,sizeof(int));
  nhits5_minus = (int *) CALLOC(cutoff_level_5+1,sizeof(int));
  for (i = 0; i < narray5; i++) {
    debug5(printf("  array5 score level %d with %d hits\n",i,List_length(hitarray5[i])));
    for (q = hitarray5[i]; q != NULL; q = q->rest) {
      hit5 = (T) q->first;
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
      hits5_plus[score5] = (T *) CALLOC(nhits5_plus[score5],sizeof(Stage3_T));
    }
  }

  hits5_minus = (T **) CALLOC(cutoff_level_5+1,sizeof(T *));
  for (score5 = 0; score5 <= cutoff_level_5; score5++) {
    if (nhits5_minus[score5] == 0) {
      hits5_minus[score5] = (T *) NULL;
    } else {
      hits5_minus[score5] = (T *) CALLOC(nhits5_minus[score5],sizeof(Stage3_T));
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
	  copy = Stage3_copy(hit5);
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

	  copy = Stage3_copy(hit5);
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
      hits3_plus[score3] = (T *) CALLOC(nhits3_plus[score3],sizeof(Stage3_T));
    }
  }

  hits3_minus = (T **) CALLOC(cutoff_level_3+1,sizeof(T *));
  for (score3 = 0; score3 <= cutoff_level_3; score3++) {
    if (nhits3_minus[score3] == 0) {
      hits3_minus[score3] = (T *) NULL;
    } else {
      hits3_minus[score3] = (T *) CALLOC(nhits3_minus[score3],sizeof(Stage3_T));
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
	  copy = Stage3_copy(hit3);
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

	  copy = Stage3_copy(hit3);
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
    debug5(printf("pairscore = %d\n",pairscore));
    for (score5 = 0; score5 <= pairscore; score5++) {
      debug5(printf("score5 = %d (cutoff %d), score3 = %d (cutoff %d)\n",
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
	    debug5(printf("plus/plus: i=%d/%d %u..%u %s %s\n",
			  i,nhits5,hit5->genomicstart,hit5->genomicend,
			  print_sense(hit5->sensedir),hittype_string(hit5->hittype)));

	    while (j >= 0 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype)));
	      j--;
	    }
	    j++;		/* Finish backup */

	    while (j < nhits3 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype)));
	      j++;
	    }

	    while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype)));
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
		} else if (SENSE_INCONSISTENT_P(hit5->sensedir,hit3->sensedir)) {
		  debug5(printf(" => sense inconsistent"));
		} else if (hit3->genomicend < hit5->genomicstart) {
		  debug5(printf(" => scramble because end3 %u < start5 %u\n",hit3->genomicend,hit5->genomicstart));
		  *samechr = List_push(*samechr,
				       (void *) Stage3pair_new(hit5,hit3,query5,query3,
							       splicesites,query5_compress_fwd,query5_compress_rev,
							       query3_compress_fwd,query3_compress_rev,genome_blocks,snp_blocks,
							       expected_pairlength,pairlength_deviation,
							       /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,dibasep,cmetp));
		} else {
		  debug5(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  hitpairs = List_push(hitpairs,
				       (void *) Stage3pair_new(hit5,hit3,query5,query3,
							       splicesites,query5_compress_fwd,query5_compress_rev,
							       query3_compress_fwd,query3_compress_rev,genome_blocks,snp_blocks,
							       expected_pairlength,pairlength_deviation,
							       /*pairtype*/CONCORDANT,splicing_penalty,dibasep,cmetp));
		  if (pairscore < *found_score) {
		    *found_score = pairscore;
		    debug5(printf(" => updating found_score to be %d",*found_score));
		  }
		  if (++(*nconcordant) > maxpairedpaths) {
		    debug(printf(" -- %d concordant paths exceeds %d",*nconcordant,maxpairedpaths));
		    *abort_pairing_p = true;
		  }
		  /* hit5->new_paired_seenp = true; -- Mark at end of round */
		  /* hit3->new_paired_seenp = true; -- Mark at end of round */
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
	    debug5(printf("minus/minus: i=%d/%d %u..%u %s %s\n",
			  i,nhits3,hit3->genomicstart,hit3->genomicend,
			  print_sense(hit3->sensedir),hittype_string(hit3->hittype)));

	    while (j >= 0 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype)));
	      j--;
	    }
	    j++;			/* Finish backup */

	    while (j < nhits5 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype)));
	      j++;
	    }

	    while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype)));
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
		} else if (SENSE_INCONSISTENT_P(hit3->sensedir,hit5->sensedir)) {
		  debug5(printf(" => sense inconsistent"));
		} else if (hit5->genomicstart < hit3->genomicend) {
		  debug5(printf(" => scramble because start5 %u < end3 %u\n",hit5->genomicstart,hit3->genomicend));
		  *samechr = List_push(*samechr,
				       (void *) Stage3pair_new(hit5,hit3,query5,query3,
							       splicesites,query5_compress_fwd,query5_compress_rev,
							       query3_compress_fwd,query3_compress_rev,genome_blocks,snp_blocks,
							       expected_pairlength,pairlength_deviation,
							       /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,dibasep,cmetp));
		} else {
		  debug5(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  hitpairs = List_push(hitpairs,
				       (void *) Stage3pair_new(hit5,hit3,query5,query3,
							       splicesites,query5_compress_fwd,query5_compress_rev,
							       query3_compress_fwd,query3_compress_rev,genome_blocks,snp_blocks,
							       expected_pairlength,pairlength_deviation,
							       /*pairtype*/CONCORDANT,splicing_penalty,dibasep,cmetp));
		  if (pairscore < *found_score) {
		    *found_score = pairscore;
		    debug5(printf(" => updating found_score to be %d",*found_score));
		  }
		  if (++(*nconcordant) > maxpairedpaths) {
		    debug(printf(" -- %d concordant paths exceeds %d",*nconcordant,maxpairedpaths));
		    *abort_pairing_p = true;
		  }
		  /* hit3->new_paired_seenp = true; -- Mark at end of round */
		  /* hit3->new_paired_seenp = true; -- Mark at end of round */
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
	    debug5(printf("plus/minus: i=%d/%d %u..%u %s %s\n",
			  i,nhits5,hit5->genomicstart,hit5->genomicend,
			  print_sense(hit5->sensedir),hittype_string(hit5->hittype)));

	    while (j >= 0 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype)));
	      j--;
	    }
	    j++;		/* Finish backup */

	    while (j < nhits3 && 
		   hits3[j]->genomicstart + querylength3 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s\n",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype)));
	      j++;
	    }

	    while (j < nhits3 && hits3[j]->genomicstart + querylength3 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s",
			    j,nhits3,hits3[j]->genomicstart,hits3[j]->genomicend,
			    print_sense(hits3[j]->sensedir),hittype_string(hits3[j]->hittype)));
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
		  *samechr = List_push(*samechr,
				       (void *) Stage3pair_new(hit5,hit3,query5,query3,
							       splicesites,query5_compress_fwd,query5_compress_rev,
							       query3_compress_fwd,query3_compress_rev,genome_blocks,snp_blocks,
							       expected_pairlength,pairlength_deviation,
							       /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,dibasep,cmetp));
#endif
		} else {
		  debug5(printf(" => inversion effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  *samechr = List_push(*samechr,
				       (void *) Stage3pair_new(hit5,hit3,query5,query3,
							       splicesites,query5_compress_fwd,query5_compress_rev,
							       query3_compress_fwd,query3_compress_rev,genome_blocks,snp_blocks,
							       expected_pairlength,pairlength_deviation,
							       /*pairtype*/PAIRED_INVERSION,splicing_penalty,dibasep,cmetp));
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
	    debug5(printf("minus/plus: i=%d/%d %u..%u %s %s\n",
			  i,nhits3,hit3->genomicstart,hit3->genomicend,
			  print_sense(hit3->sensedir),hittype_string(hit3->hittype)));

	    while (j >= 0 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax > insert_start) {
	      debug5(printf("  backup: j=%d/%d %u..%u %s %s\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype)));
	      j--;
	    }
	    j++;			/* Finish backup */

	    while (j < nhits5 && 
		   hits5[j]->genomicend + querylength5 /* for scramble: */ + pairmax <= insert_start) {
	      debug5(printf("  advance: j=%d/%d %u..%u %s %s\n",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype)));
	      j++;
	    }

	    while (j < nhits5 && hits5[j]->genomicend + querylength5 <= pairmax + insert_start) {
	      debug5(printf("  overlap: j=%d/%d %u..%u %s %s",
			    j,nhits5,hits5[j]->genomicstart,hits5[j]->genomicend,
			    print_sense(hits5[j]->sensedir),hittype_string(hits5[j]->hittype)));
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
		  *samechr = List_push(*samechr,
				       (void *) Stage3pair_new(hit5,hit3,query5,query3,
							       splicesites,query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
							       genome_blocks,snp_blocks,
							       expected_pairlength,pairlength_deviation,
							       /*pairtype*/PAIRED_SCRAMBLE,splicing_penalty,dibasep,cmetp));
#endif
		} else {
		  debug5(printf(" => inversion effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  *samechr = List_push(*samechr,
				       (void *) Stage3pair_new(hit5,hit3,query5,query3,
							       splicesites,query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
							       genome_blocks,snp_blocks,
							       expected_pairlength,pairlength_deviation,
							       /*pairtype*/PAIRED_INVERSION,splicing_penalty,dibasep,cmetp));
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
    pairscore++;
  }


  FREE(sorted3p);
  FREE(sorted5p);

#if 0
  for (q = copies; q != NULL; q = List_next(q)) {
    copy = (T) List_head(q);
    Stage3_free(&copy);
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

  debug5(printf("Finished with Stage3_pair_up_concordant: %d pairs\n",List_length(hitpairs)));

  return hitpairs;
}


