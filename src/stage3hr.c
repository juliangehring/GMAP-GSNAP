static char rcsid[] = "$Id: stage3hr.c,v 1.107 2010-07-27 00:02:16 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage3hr.h"
#include "stage1hr.h"		/* To get MAX_QUERYLENGTH */

#include <stdlib.h>		/* For qsort */
#include <string.h>
#include <strings.h>
#include <ctype.h>		/* For islower */
#include "mem.h"
#include "chrnum.h"
#include "complement.h"
#include "maxent.h"
#include "interval.h"
#include "listdef.h"
#include "substring.h"


#define MAX_HITS 100000

#if 0
#define TRANSLOCATION_TEXT "pair_translocation"
#define INVERSION_TEXT "pair_inversion"
#define SCRAMBLE_TEXT "pair_scramble"
#define TOOLONG_TEXT "pair_toolong"
#endif

#define CONCORDANT_TEXT "concordant"
#define SAMECHR_TEXT "samechr"
#define UNPAIRED_TEXT "unpaired"

#define SAMECHR_TOOLONG_TEXT "toolong"
#define SAMECHR_INVERSION_TEXT "inversion"
#define SAMECHR_SCRAMBLE_TEXT "scramble"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
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

/* Stage3_pair_up_concordant */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Stage3_pair_up_samechr */
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

/* Stage3_mark_ambiguous_splices */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif



#ifdef DEBUG5
static char *
print_sense (int sense) {
  if (sense == SENSE_NULL) {
    return "null";
  } else if (sense == SENSE_ANTI) {
    return "antisense";
  } else if (sense == SENSE_FORWARD) {
    return "sense";
  } else {
    abort();
  }
}
#endif


void
Stage3hr_print_nsnpdiffs (bool labelsp) {
  Substring_print_nsnpdiffs(labelsp);
  return;
}

void
Stage3hr_print_ncolordiffs () {
  Substring_print_ncolordiffs();
  return;
}


/* Note: Substring_T has genomiclength, but not Stage3_T */
#define T Stage3_T
struct T {
  Hittype_T hittype;
  Chrnum_T chrnum; /* Needed for printing paired-end results.  A chrnum of 0 indicates a distant splice. */
  Chrnum_T effective_chrnum;	/* For concordance and samechr */
  Chrnum_T other_chrnum;	/* 0 for non-translocations, and other chrnum besides effective_chrnum for translocations */
  Genomicpos_T chroffset;

  Genomicpos_T genomicstart;
  Genomicpos_T genomicend;
  bool plusp;

  int score;			/* Includes colordiffs and penalties */
  int ntscore;			/* Includes penalties */
  int nmatches;
  int total_nmismatches;
  int geneprob;

  int nindels;			/* for indels */
  int indel_pos;		/* for indels */
  char *deletion;		/* for deletions */

  double chimera_prob;		/* for splicing */
  Genomicpos_T distance;	/* for splicing */
  int sensedir;			/* for splicing */
  bool chimera_ambiguous_p;	/* for splicing */
#if 0
  double half_intron_score;	/* combination of support, nmismatches, and splice site probability */
#endif
#ifdef USE_HALFINTRON_SUPPORT
  int halfintron_support;	/* support - 3*nmismatches */
#endif

  Substring_T substring1;	/* indel or donor */
  Substring_T substring2;	/* indel or acceptor */
  Substring_T substring_low;	/* For SAM output */

  bool paired_seenp;   /* for paired-end.  set to true by Stage3_pair_up(). */
  bool concordantp;    /* for paired-end.  set to true by Stage3_pair_up(). */
};


struct Stage3pair_T {
  T hit5;
  T hit3;
  Pairtype_T pairtype;

  Genomicpos_T low;
  Genomicpos_T high;
  int pairedlength;
  int score;
  int nmatches;
  Genomicpos_T absdifflength;
  int dir;			/* -1, 0, or +1 */
  bool sense_consistent_p;
};

#if 0
struct Stage3chimera_T {
  T donor;
  T acceptor;
  Genomicpos_T distance;
  int total_nmismatches;
  double total_prob;
  bool sensep;
};
#endif


Hittype_T
Stage3_hittype (T this) {
  return this->hittype;
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


int
Stage3_score (T this) {
  return this->score;
}

int
Stage3_geneprob (T this) {
  return this->geneprob;
}

int
Stage3_nmismatches (T this) {
  return this->total_nmismatches;
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

Substring_T
Stage3_substring_low (T this) {
  return this->substring_low;
}

Genomicpos_T
Stage3_distance (T this) {
  return this->distance;
}

int
Stage3_sensedir (T this) {
  return this->sensedir;
}



void
Stage3_free (T *old) {
  if ((*old)->deletion != NULL) {
    FREE((*old)->deletion);
  }
  if ((*old)->substring1 != NULL) {
    Substring_free(&(*old)->substring1);
  }
  if ((*old)->substring2 != NULL) {
    Substring_free(&(*old)->substring2);
  }
  FREE(*old);
  return;
}

Pairtype_T
Stage3pair_pairtype (Stage3pair_T this) {
  return this->pairtype;
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
Stage3pair_hit5_score (Stage3pair_T this) {
  return this->hit5->score;
}

int
Stage3pair_hit3_score (Stage3pair_T this) {
  return this->hit3->score;
}

Genomicpos_T
Stage3pair_pairlength (Stage3pair_T this) {
  return this->pairedlength;
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
  new->genomicstart = old->genomicstart;
  new->genomicend = old->genomicend;
  new->plusp = old->plusp;

  new->score = old->score;
  new->ntscore = old->ntscore;
  new->total_nmismatches = old->total_nmismatches;
  new->geneprob = old->geneprob;

  new->nindels = old->nindels;
  new->indel_pos = old->indel_pos;
  if (old->deletion == NULL) {
    new->deletion = (char *) NULL;
  } else {
    new->deletion = (char *) CALLOC(strlen(old->deletion)+1,sizeof(char));
    strcpy(new->deletion,old->deletion);
  }
  
  new->chimera_prob = old->chimera_prob;
  new->chimera_ambiguous_p = old->chimera_ambiguous_p;
  new->distance = old->distance;
#if 0
  new->half_intron_score = old->half_intron_score;
#endif
#ifdef USE_HALFINTRON_SUPPORT
  new->halfintron_support = old->halfintron_support;
#endif
  new->sensedir = old->sensedir;

  new->substring1 = Substring_copy(old->substring1);
  new->substring2 = Substring_copy(old->substring2);
  if (old->substring_low == NULL) {
    new->substring_low = NULL;
  } else if (old->substring_low == old->substring1) {
    new->substring_low = new->substring1;
  } else if (old->substring_low == old->substring2) {
    new->substring_low = new->substring2;
  } else {
    fprintf(stderr,"substring_low is not NULL, substring1, or substring2\n");
    abort();
  }

  new->paired_seenp = old->paired_seenp;
  new->concordantp = old->concordantp;

  return new;
}


T
Stage3_new_exact (int *found_score, Genomicpos_T left, int genomiclength, bool plusp, Chrnum_T chrnum, Genomicpos_T chroffset) {
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

  if ((substring = Substring_new(/*nmismatches*/0,/*ncolordiffs*/0,chrnum,
				 chroffset,genomicstart,genomicend,
				 /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				 /*alignstart*/genomicstart,/*alignend*/genomicend,
				 genomiclength,/*extraleft*/0,/*extraright*/0,
				 /*genomicseg*/NULL,/*query*/NULL,plusp,
				 /*trim_left_p*/false,/*trim_right_p*/false,/*trim_maxlength*/0,
				 /*minlength*/0,/*dibasep*/false,/*cmetp*/false)) == NULL) {
    return (T) NULL;

  } else {
    new = (T) MALLOC(sizeof(*new));
    new->substring1 = substring;
    new->substring2 = (Substring_T) NULL;
    new->substring_low = new->substring1;

    new->deletion = (char *) NULL;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    new->hittype = EXACT;
    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

    new->nindels = 0;
    new->total_nmismatches = 0;
    new->ntscore = 0;
    new->score = 0;
    new->nmatches = genomiclength;
    new->geneprob = -1;
    *found_score = 0;

    new->chimera_prob = 2.0;
    new->chimera_ambiguous_p = false;
    new->distance = 0U;
#if 0
    new->half_intron_score = 0.0;
#endif
#ifdef USE_HALFINTRON_SUPPORT
    new->halfintron_support = 0;
#endif

    new->paired_seenp = false;
    new->concordantp = false;

    return new;
  }
}


T
Stage3_new_substitution (int *found_score, int nmismatches, int ncolordiffs, Genomicpos_T left,
			 int genomiclength, bool plusp, char *genomicseg, char *query,
			 Chrnum_T chrnum, Genomicpos_T chroffset, 
			 int trim_maxlength, bool dibasep, bool cmetp) {
  T new;
  Substring_T substring;
  Genomicpos_T genomicstart, genomicend;

  debug(printf("Entered new_substitution with nmismatches %d\n",nmismatches));

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + genomiclength;
  } else {
    genomicstart = left + genomiclength;
    genomicend = left;
  }

  if ((substring = Substring_new(nmismatches,ncolordiffs,chrnum,
				 chroffset,genomicstart,genomicend,
				 /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				 /*alignstart*/genomicstart,/*alignend*/genomicend,
				 genomiclength,/*extraleft*/0,/*extraright*/0,
				 genomicseg,query,plusp,
				 /*trim_left_p*/true,/*trim_right_p*/true,trim_maxlength,
				 /*minlength*/genomiclength/2,dibasep,cmetp)) == NULL) {
    return (T) NULL;

  } else {
    new = (T) MALLOC(sizeof(*new));
    new->substring1 = substring;
    new->substring2 = (Substring_T) NULL;
    new->substring_low = new->substring1;

    new->deletion = (char *) NULL;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    new->hittype = SUB;
    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

    new->nindels = 0;
    new->total_nmismatches = nmismatches;
    new->ntscore = nmismatches;
    new->score = nmismatches + ncolordiffs;
#if 0
    new->nmatches = Substring_match_length(new->substring1) - new->score;
#else
    new->nmatches = Substring_nmatches(new->substring1);
#endif
    new->geneprob = -1;

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->chimera_prob = 2.0;
    new->chimera_ambiguous_p = false;
    new->distance = 0U;
#if 0
    new->half_intron_score = 0.0;
#endif
#ifdef USE_HALFINTRON_SUPPORT
    new->halfintron_support = 0;
#endif

    new->paired_seenp = false;
    new->concordantp = false;

    return new;
  }
}



T
Stage3_new_insertion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
		      int ncolordiffs1, int ncolordiffs2,
		      Genomicpos_T left, int genomiclength, int querylength, bool plusp,
		      char *genomicseg, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
		      int indel_penalty, int trim_maxlength, bool dibasep, bool cmetp) {
  T new;
  Substring_T substring1, substring2;
  int querystart1, queryend1, querystart2, queryend2;
  Genomicpos_T genomicstart, genomicend;
  Genomicpos_T alignstart1, alignend1, alignstart2, alignend2;

  debug2(printf("Entered with left %u, querylength %d, genomiclength %d, indel_pos %d\n",
		left,querylength,genomiclength,indel_pos));
  debug2(printf("q: %s\n",query));
  debug2(printf("g: %s\n",genomicseg));

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

  if ((substring1 = Substring_new(nmismatches1,ncolordiffs1,chrnum,
				  chroffset,genomicstart,genomicend,
				  querystart1,queryend1,querylength,
				  alignstart1,alignend1,genomiclength,
				  /*extraleft*/0,/*extraright*/0,
				  genomicseg,query,plusp,
				  /*trim_left_p*/true,/*trim_right_p*/false,
				  trim_maxlength,/*minlength*/0,dibasep,cmetp)) == NULL) {
    return (T) NULL;

  } else if ((substring2 = Substring_new(nmismatches2,ncolordiffs2,chrnum,
					 chroffset,genomicstart,genomicend,
					 querystart2,queryend2,querylength,
					 alignstart2,alignend2,genomiclength,
					 /*extraleft*/0,/*extraright*/0,
					 genomicseg,query,plusp,
					 /*trim_left_p*/false,/*trim_right_p*/true,
					 trim_maxlength,/*minlength*/0,dibasep,cmetp)) == NULL) {
    Substring_free(&substring1);
    return (T) NULL;

  } else {
    new = (T) MALLOC(sizeof(*new));
    new->substring1 = substring1;
    new->substring2 = substring2;
    if (plusp == true) {
      new->substring_low = new->substring1;
    } else {
      new->substring_low = new->substring2;
    }

    new->deletion = (char *) NULL;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    new->hittype = INS;
    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

    new->nindels = nindels;
    new->total_nmismatches = nmismatches1 + nmismatches2;
    new->ntscore = indel_penalty + nmismatches1 + nmismatches2;
    new->score = new->ntscore + ncolordiffs1 + ncolordiffs2;
#if 0
    new->nmatches = Substring_match_length(new->substring1) + Substring_match_length(new->substring2) - new->score;
#else
    new->nmatches = Substring_nmatches(new->substring1) + Substring_nmatches(new->substring2);
#endif
    new->geneprob = -1;

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->indel_pos = indel_pos;

    new->chimera_prob = 2.0;
    new->chimera_ambiguous_p = false;
    new->distance = 0U;
#if 0
    new->half_intron_score = 0.0;
#endif
#ifdef USE_HALFINTRON_SUPPORT
    new->halfintron_support = 0;
#endif

    new->paired_seenp = false;
    new->concordantp = false;

    return new;
  }
}


T
Stage3_new_deletion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
		     int ncolordiffs1, int ncolordiffs2,
		     Genomicpos_T left, int genomiclength, int querylength, bool plusp,
		     char *genomicseg, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
		     int indel_penalty, int trim_maxlength, bool dibasep, bool cmetp) {
  T new;
  Substring_T substring1, substring2;
  int querystart1, queryend1, querystart2, queryend2;
  Genomicpos_T genomicstart, genomicend;
  Genomicpos_T alignstart1, alignend1, alignstart2, alignend2;

  debug3(printf("Entered with left %u, querylength %d, genomiclength %d, indel_pos %d\n",
		left,querylength,genomiclength,indel_pos));
  debug3(printf("q: %s\n",query));
  debug3(printf("g: %s\n",genomicseg));

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

  } else {
    genomicstart = left + genomiclength;
    genomicend = left;

    alignstart1 = genomicstart;
    alignend1 = genomicstart - indel_pos;
    alignstart2 = alignend1 - nindels;
    alignend2 = genomicend/* - nindels*/;
  }

  if ((substring1 = Substring_new(nmismatches1,ncolordiffs1,chrnum,chroffset,
				  genomicstart,genomicend,
				  querystart1,queryend1,querylength,
				  alignstart1,alignend1,genomiclength,
				  /*extraleft*/0,/*extraright*/0,	
				  genomicseg,query,plusp,
				  /*trim_left_p*/true,/*trim_right_p*/false,
				  trim_maxlength,/*minlength*/0,dibasep,cmetp)) == NULL) {
    return (T) NULL;

  } else if ((substring2 = Substring_new(nmismatches2,ncolordiffs2,chrnum,chroffset,
					 genomicstart,genomicend,
					 querystart2,queryend2,querylength,
					 alignstart2,alignend2,genomiclength,
					 /*extraleft*/0,/*extraright*/0,
					 genomicseg,query,plusp,
					 /*trim_left_p*/false,/*trim_right_p*/true,
					 trim_maxlength,/*minlength*/0,dibasep,cmetp)) == NULL) {
    Substring_free(&substring1);
    return (T) NULL;
    
  } else {
    new = (T) MALLOC(sizeof(*new));
    new->substring1 = substring1;
    new->substring2 = substring2;

    new->deletion = (char *) CALLOC(nindels+1,sizeof(char));
    if (plusp == true) {
      strncpy(new->deletion,&(genomicseg[indel_pos]),nindels);
      new->substring_low = new->substring1;
    } else {
      make_complement_buffered(new->deletion,&(genomicseg[querylength-indel_pos]),nindels);
      new->substring_low = new->substring2;
    }
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    new->hittype = DEL;
    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->plusp = plusp;
    new->sensedir = SENSE_NULL;

    new->nindels = nindels;
    new->total_nmismatches = nmismatches1 + nmismatches2;
    new->ntscore = indel_penalty + nmismatches1 + nmismatches2;
    new->score = new->ntscore + ncolordiffs1 + ncolordiffs2;
#if 0
    new->nmatches = Substring_match_length(new->substring1) + Substring_match_length(new->substring2) - new->score;
#else
    new->nmatches = Substring_nmatches(new->substring1) + Substring_nmatches(new->substring2);
#endif
    new->geneprob = -1;

    if (new->score < *found_score) {
      *found_score = new->score;
    }

    new->indel_pos = indel_pos;

    new->chimera_prob = 2.0;
    new->chimera_ambiguous_p = false;
    new->distance = 0U;
#if 0
    new->half_intron_score = 0.0;
#endif
#ifdef USE_HALFINTRON_SUPPORT
    new->halfintron_support = 0;
#endif

    new->paired_seenp = false;
    new->concordantp = false;

    return new;
  }
}


T
Stage3_new_terminal (int querystart, int queryend, int nmismatches, int ncolordiffs,
		     Genomicpos_T left, int querylength,
		     bool plusp, char *genomicseg, char *query, 
		     Chrnum_T chrnum, Genomicpos_T chroffset,
		     int terminal_penalty, int trim_maxlength, bool dibasep, bool cmetp) {
  T new;
  Substring_T substring;
  Genomicpos_T genomicstart, genomicend, alignstart, alignend;

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

  if ((substring = Substring_new(nmismatches,ncolordiffs,chrnum,chroffset,
				 genomicstart,genomicend,querystart,queryend,querylength,
				 alignstart,alignend,/*genomiclength*/querylength,
				 /*extraleft*/0,/*extraright*/0,genomicseg,query,plusp,
				 /*trim_left_p*/true,/*trim_right_p*/true,trim_maxlength,
				 /*minlength*/querylength/2,dibasep,cmetp)) == NULL) {
    return (T) NULL;
    
  } else {
    new = (T) MALLOC(sizeof(*new));
    new->substring1 = substring;
    new->substring2 = (Substring_T) NULL;
    new->substring_low = new->substring1;

    new->deletion = (char *) NULL;
    new->genomicstart = genomicstart;
    new->genomicend = genomicend;

    new->hittype = TERMINAL;
    new->chrnum = new->effective_chrnum = chrnum;
    new->other_chrnum = 0;
    new->chroffset = chroffset;
    new->plusp = plusp;
    new->total_nmismatches = nmismatches;
    new->ntscore = terminal_penalty + nmismatches;
    new->score = terminal_penalty + nmismatches + ncolordiffs;
#if 0
    new->nmatches = Substring_match_length(new->substring1) - (nmismatches + ncolordiffs);
#else
    new->nmatches = Substring_nmatches(new->substring1);
#endif
    /* new->chimera_prob = Substring_chimera_prob(acceptor); */
    new->chimera_ambiguous_p = false;

    new->distance = 0U;
    new->sensedir = SENSE_NULL;
    new->geneprob = -1;

    new->paired_seenp = false;
    new->concordantp = false;

    return new;
  }
}


int
Stage3_output_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->nmatches > y->nmatches) {
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


List_T
Stage3_optimal_score (List_T hitlist, int cutoff_level, int suboptimal_mismatches) {
  List_T optimal = NULL, p;
  T hit;
  int n;
  int minscore = MAX_QUERYLENGTH;
  int best_terminal_nmatches = 0;

  n = List_length(hitlist);
  if (n == 0) {
    return NULL;
  }

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    /* For dibasep were previously using hit->ntscore, but gives false positives */
    if (hit->score <= cutoff_level) {
      if (hit->hittype == TERMINAL) {
	if (hit->nmatches > best_terminal_nmatches) {
	  best_terminal_nmatches = hit->nmatches;
	}
      } else if (hit->chimera_ambiguous_p == false) {
	if (hit->score < minscore) {
	  minscore = hit->score;
	}
      }
    }
  }

  debug(printf("Stage3_optimal_score over %d hits: minscore = %d + subopt:%d, best_terminal_nmatches %d\n",
	       n,minscore,suboptimal_mismatches,best_terminal_nmatches));
  minscore += suboptimal_mismatches;

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    /* For dibasep were previously using hit->ntscore, but gives false positives */
    if (hit->score > cutoff_level) {
      debug(printf("Eliminating a hit with ntscore %d > cutoff_level %d\n",
		   hit->ntscore,cutoff_level));
      Stage3_free(&hit);
    } else if (hit->chimera_ambiguous_p == true) {
      debug(printf("Eliminating a hit with ntscore %d and ambiguous splice\n",
		   hit->ntscore));
      Stage3_free(&hit);
    } else if (hit->score <= minscore) {
      if (hit->hittype == TERMINAL) {
	if (hit->nmatches < best_terminal_nmatches) {
	  debug(printf("Eliminating a terminal hit with nmatches %d\n",hit->nmatches));
	  Stage3_free(&hit);
	} else {
	  debug(printf("Keeping a terminal hit with nmatches %d\n",hit->nmatches));
	  optimal = List_push(optimal,hit);
	}
      } else {
	debug(printf("Keeping a hit with score %d\n",hit->score));
	optimal = List_push(optimal,hit);
      }
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
  } else {
    return 0;
  }
}


List_T
Stage3_mark_ambiguous_splices (bool *ambiguousp, List_T hitlist) {
  T x, y, *hits;
  int n, i, j;

  *ambiguousp = false;
  n = List_length(hitlist);
  debug9(printf("Entered Stage3_mark_ambiguous_splices with %d pairs\n",n));
  if (n == 0) {
    return NULL;
  }

  hits = (T *) List_to_array(hitlist,NULL);

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
	    x->chimera_ambiguous_p = true;
	    y->chimera_ambiguous_p = true;
	    *ambiguousp = true;
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
	    x->chimera_ambiguous_p = true;
	    y->chimera_ambiguous_p = true;
	    *ambiguousp = true;
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
  }

  eliminate = (bool *) CALLOC(n,sizeof(bool));

  hits = (T *) List_to_array(hitlist,NULL);
  List_free(&hitlist);


  /* By genomicstart */
  debug7(printf("Stage3_remove_duplicates: checking %d hits by genomicstart\n",n));
  qsort(hits,n,sizeof(T),genomicstart_cmp);

  debug7(
	 for (i = 0; i < n; i++) {
	   hit = hits[i];
	   printf("  Initial %d: #%d:%u..%u, nmatches: %d\n",
		  i,hit->chrnum,hit->genomicstart-hit->chroffset,hit->genomicend-hit->chroffset,hit->nmatches);
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
#if 0
	} else if (x->score < y->score) {
	  debug7(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
		       i,j,x->score,y->score,j));
	  eliminate[j] = true;
	} else if (x->score > y->score) {
	  debug7(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
		       i,j,x->score,y->score,i));
	  eliminate[i] = true;
#else
	} else if (x->nmatches > y->nmatches) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d > %d, so marking %d for elimination\n",
		       i,j,x->nmatches,y->nmatches,j));
	  eliminate[j] = true;
	} else if (x->nmatches < y->nmatches) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d < %d, so marking %d for elimination\n",
		       i,j,x->nmatches,y->nmatches,i));
	  eliminate[i] = true;
#endif
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

	} else if (x->distance > y->distance) {
	  debug7(printf("  #%d overlaps #%d and first one has longer distance, so marking first one for elimination\n",
		       i,j));
	  eliminate[i] = true;
	} else if (x->distance < y->distance) {
	  debug7(printf("  #%d overlaps #%d and second one has longer distance, so marking second one for elimination\n",
		       i,j));
	  eliminate[j] = true;

	} else if (x->genomicend == y->genomicend) {
#if 0
	  if (x->half_intron_score != 0.0 && y->half_intron_score != 0.0) {
	    if (x->half_intron_score > y->half_intron_score) {
	      debug7(printf("  #%d overlaps #%d and equal and first has better support score %f > %f, so marking second one for elimination\n",
			   i,j,x->half_intron_score,y->half_intron_score));
	      eliminate[j] = true;
	    } else if (y->half_intron_score > x->half_intron_score) {
	      debug7(printf("  #%d overlaps #%d and equal and second has better support score %f > %f, so marking first one for elimination\n",
			   i,j,y->half_intron_score,x->half_intron_score));
	      eliminate[i] = true;
	    } else {
	      debug7(printf("  #%d overlaps #%d and equal so eliminate second one\n",
			   i,j));
	      eliminate[j] = true;
	    }
	  } else {
#endif
	    debug7(printf("  #%d overlaps #%d and equal so eliminate second one\n",
			 i,j));
	    eliminate[j] = true;
#if 0
	  }
#endif
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
#if 0
	} else if (x->score < y->score) {
	  debug7(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
		       i,j,x->score,y->score,j));
	  eliminate[j] = true;
	} else if (x->score > y->score) {
	  debug7(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
		       i,j,x->score,y->score,i));
	  eliminate[i] = true;
#else
	} else if (x->nmatches > y->nmatches) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d > %d, so marking %d for elimination\n",
		       i,j,x->nmatches,y->nmatches,j));
	  eliminate[j] = true;
	} else if (x->nmatches < y->nmatches) {
	  debug7(printf("  #%d overlaps #%d and nmatches %d < %d, so marking %d for elimination\n",
		       i,j,x->nmatches,y->nmatches,i));
	  eliminate[i] = true;
#endif
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

	} else if (x->distance > y->distance) {
	  debug7(printf("  #%d overlaps #%d and first one has longer distance, so marking first one for elimination\n",
		       i,j));
	  eliminate[i] = true;
	} else if (x->distance < y->distance) {
	  debug7(printf("  #%d overlaps #%d and second one has longer distance, so marking second one for elimination\n",
		       i,j));
	  eliminate[j] = true;

	} else if (x->genomicstart == y->genomicstart) {
#if 0
	  if (x->half_intron_score != 0.0 && y->half_intron_score != 0.0) {
	    if (x->half_intron_score > y->half_intron_score) {
	      debug7(printf("  #%d overlaps #%d and equal and first has better support score %f > %f, so marking second one for elimination\n",
			   i,j,x->half_intron_score,y->half_intron_score));
	      eliminate[j] = true;
	    } else if (y->half_intron_score > x->half_intron_score) {
	      debug7(printf("  #%d overlaps #%d and equal and second has better support score %f > %f, so marking first one for elimination\n",
			   i,j,y->half_intron_score,x->half_intron_score));
	      eliminate[i] = true;
	    } else {
	      debug7(printf("  #%d overlaps #%d and equal so eliminate second one\n",
			   i,j));
	      eliminate[j] = true;
	    }
	  } else {
#endif
	    debug7(printf("  #%d overlaps #%d and equal so eliminate second one\n",
			 i,j));
	    eliminate[j] = true;
#if 0
	  }
#endif
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
print_alignment_info (int nblocks, int score, int geneprob) {
  printf("segs:%d,align_score:%d",nblocks,score);
  if (geneprob >= 0) {
    printf(",geneprob:%d",geneprob);
  }
  return;
}


static char *
samechr_type_text (Genomicpos_T genomicstart5, bool plus5p,
		   Genomicpos_T genomicstart3, bool plus3p, Genomicpos_T pairedlength,
		   Genomicpos_T pairmax) {
  if (plus5p != plus3p) {
    return SAMECHR_INVERSION_TEXT;
  } else if (plus5p == true) {
    if (genomicstart3 < genomicstart5) {
      return SAMECHR_SCRAMBLE_TEXT;
    } else if (pairedlength < pairmax) {
      return CONCORDANT_TEXT;
    } else {
      return SAMECHR_TOOLONG_TEXT;
    }
  } else {
    if (genomicstart5 < genomicstart3) {
      return SAMECHR_SCRAMBLE_TEXT;
    } else if (pairedlength < pairmax) {
      return CONCORDANT_TEXT;
    } else {
      return SAMECHR_TOOLONG_TEXT;
    }
  }
}


static void
print_pair_info (T hit5, T hit3, Pairtype_T pairtype, int pairedlength, int pairscore) {

  printf("pair_score:%d",pairscore);
  if (pairtype == CONCORDANT) {
    printf(",pair_length:%d",pairedlength);
  } else if (hit5->chrnum != hit3->chrnum) {
    /* Different chromosomes */
  } else if (hit5->plusp != hit3->plusp) {
    printf(",pair_length:%d",pairedlength);
    /* printf(",%s",INVERSION_TEXT); */
  } else if (hit5->plusp == true) {
    printf(",pair_length:%d",pairedlength);
    if (hit3->genomicstart < hit5->genomicstart) {
      /* printf(",%s",SCRAMBLE_TEXT); */
    } else {
      /* printf(",%s",TOOLONG_TEXT); */
    }
  } else {
    printf(",pair_length:%d",pairedlength);
    if (hit5->genomicstart < hit3->genomicstart) {
      /* printf(",%s",SCRAMBLE_TEXT); */
    } else {
      /* printf(",%s",TOOLONG_TEXT); */
    }
  }

  return;
}

static void
print_single (T this, int score, int geneprob, IIT_T chromosome_iit, Sequence_T queryseq,
	      bool invertp, IIT_T snps_iit, int *snps_divint_crosstable,
	      T hit5, T hit3, Pairtype_T pairtype, int pairedlength, int pairscore) {
  int querylength;
  char *chr;
  bool allocp;

  querylength = Sequence_fulllength(queryseq);
  chr = IIT_label(chromosome_iit,this->chrnum,&allocp);

  printf(" ");
  Substring_print_single(this->substring1,this->hittype,queryseq,chr,querylength,invertp,
			 snps_iit,snps_divint_crosstable);

  /* Alignment info */
  printf("\t");
  print_alignment_info(/*nblocks*/1,score,geneprob);

  /* Pairing info */
  if (hit5 != NULL && hit3 != NULL) {
    printf("\t");
    print_pair_info(hit5,hit3,pairtype,pairedlength,pairscore);
  }

  printf("\n");

  if (allocp == true) {
    FREE(chr);
  }

  return;
}


static void
print_insertion (T this, int score, int geneprob, IIT_T chromosome_iit, Sequence_T queryseq,
		 bool invertp, IIT_T snps_iit, int *snps_divint_crosstable,
		 T hit5, T hit3, Pairtype_T pairtype, int pairedlength, int pairscore) {
  char *chr;
  bool allocp;

  chr = IIT_label(chromosome_iit,this->chrnum,&allocp);

  printf(" ");
  Substring_print_insertion_1(this->substring1,this->substring2,this->nindels,
			      queryseq,chr,invertp,snps_iit,snps_divint_crosstable);
  /* Alignment info */
  printf("\t");
  print_alignment_info(/*nblocks*/2,score,geneprob);

  /* Pairing info */
  if (hit5 != NULL && hit3 != NULL) {
    printf("\t");
    print_pair_info(hit5,hit3,pairtype,pairedlength,pairscore);
  }

  printf("\n");


  printf(",");
  Substring_print_insertion_2(this->substring1,this->substring2,this->nindels,
			      queryseq,chr,invertp,snps_iit,snps_divint_crosstable);
  printf("\n");

  if (allocp == true) {
    FREE(chr);
  }

  return;
}

static void
print_deletion (T this, int score, int geneprob, IIT_T chromosome_iit, Sequence_T queryseq,
		bool invertp, IIT_T snps_iit, int *snps_divint_crosstable,
		T hit5, T hit3, Pairtype_T pairtype, int pairedlength, int pairscore) {
  char *chr;
  bool allocp;

  chr = IIT_label(chromosome_iit,this->chrnum,&allocp);

  printf(" ");
  Substring_print_deletion_1(this->substring1,this->substring2,this->nindels,this->deletion,
			     queryseq,chr,invertp,snps_iit,snps_divint_crosstable);
  /* Alignment info */
  printf("\t");
  print_alignment_info(/*nblocks*/2,score,geneprob);

  /* Pairing info */
  if (hit5 != NULL && hit3 != NULL) {
    printf("\t");
    print_pair_info(hit5,hit3,pairtype,pairedlength,pairscore);
  }

  printf("\n");

  printf(",");
  Substring_print_deletion_2(this->substring1,this->substring2,this->nindels,this->deletion,
			     queryseq,chr,invertp,snps_iit,snps_divint_crosstable);
  printf("\n");

  if (allocp == true) {
    FREE(chr);
  }
}


static void
print_splice (T chimera, int score, int geneprob, Genome_T genome, IIT_T chromosome_iit,
	      Sequence_T queryseq, bool invertp, IIT_T snps_iit, int *snps_divint_crosstable,
	      IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
	      int *splicesites_divint_crosstable, T hit5, T hit3,
	      Pairtype_T pairtype, int pairedlength, int pairscore) {
  Substring_T donor, acceptor;
  
  donor = chimera->substring1;
  acceptor = chimera->substring2;

  Substring_assign_donor_prob(donor,genome,chromosome_iit);
  Substring_assign_acceptor_prob(acceptor,genome,chromosome_iit);

  if (acceptor == NULL) {
    abort();			/* Using terminals now instead of halfdonor */
    /* Single sequence */
    printf(" ");
    Substring_print_donor(donor,/*sensep*/chimera->sensedir==SENSE_FORWARD,invertp,
			  queryseq,chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			  donor_typeint,splicesites_divint_crosstable,/*endp*/true,
			  /*acceptor*/NULL,/*chimera_distance*/0);
    /* Alignment info */
    printf("\t");
    print_alignment_info(/*nblocks*/1,score,geneprob);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      printf("\t");
      print_pair_info(hit5,hit3,pairtype,pairedlength,pairscore);
    }
    printf("\n");

  } else if (donor == NULL) {
    abort();			/* Using terminals now instead of halfacceptor */
    /* Single sequence */
    printf(" ");
    Substring_print_acceptor(acceptor,/*sensep*/chimera->sensedir==SENSE_FORWARD,invertp,
			     queryseq,chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			     acceptor_typeint,splicesites_divint_crosstable,/*endp*/true,
			     /*donor*/NULL,/*chimera_distance*/0);
    /* Alignment info */
    printf("\t");
    print_alignment_info(/*nblocks*/1,score,geneprob);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      printf("\t");
      print_pair_info(hit5,hit3,pairtype,pairedlength,pairscore);
    }
    printf("\n");

  } else if (chimera->sensedir == SENSE_FORWARD && invertp == false) {
    printf(" ");
    Substring_print_donor(donor,/*sensep*/true,/*invertp*/false,
			  queryseq,chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			  donor_typeint,splicesites_divint_crosstable,/*endp*/false,
			  acceptor,chimera->distance);
    /* Alignment info */
    printf("\t");
    print_alignment_info(/*nblocks*/2,score,geneprob);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      printf("\t");
      print_pair_info(hit5,hit3,pairtype,pairedlength,pairscore);
    }
    printf("\n");

    printf(",");
    Substring_print_acceptor(acceptor,/*sensep*/true,/*invertp*/false,queryseq,
			     chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			     acceptor_typeint,splicesites_divint_crosstable,/*endp*/false,
			     donor,chimera->distance);
    printf("\n");

  } else if (chimera->sensedir == SENSE_FORWARD && invertp == true) {
    printf(" ");
    Substring_print_acceptor(acceptor,/*sensep*/true,/*invertp*/true,queryseq,
			     chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			     acceptor_typeint,splicesites_divint_crosstable,/*endp*/false,
			     donor,chimera->distance);
    /* Alignment info */
    printf("\t");
    print_alignment_info(/*nblocks*/2,score,geneprob);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      printf("\t");
      print_pair_info(hit5,hit3,pairtype,pairedlength,pairscore);
    }
    printf("\n");

    printf(",");
    Substring_print_donor(donor,/*sensep*/true,/*invertp*/true,queryseq,
			  chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			  donor_typeint,splicesites_divint_crosstable,/*endp*/false,
			  acceptor,chimera->distance);
    printf("\n");

  } else if (chimera->sensedir == SENSE_ANTI && invertp == false) {
    printf(" ");
    Substring_print_acceptor(acceptor,/*sensep*/false,/*invertp*/false,queryseq,
			     chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			     acceptor_typeint,splicesites_divint_crosstable,/*endp*/false,
			     donor,chimera->distance);
    /* Alignment info */
    printf("\t");
    print_alignment_info(/*nblocks*/2,score,geneprob);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      printf("\t");
      print_pair_info(hit5,hit3,pairtype,pairedlength,pairscore);
    }
    printf("\n");

    printf(",");
    Substring_print_donor(donor,/*sensep*/false,/*invertp*/false,queryseq,
			  chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			  donor_typeint,splicesites_divint_crosstable,/*endp*/false,
			  acceptor,chimera->distance);
    printf("\n");

  } else if (chimera->sensedir == SENSE_ANTI && invertp == true) {
    printf(" ");
    Substring_print_donor(donor,/*sensep*/false,/*invertp*/true,queryseq,
			  chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			  donor_typeint,splicesites_divint_crosstable,/*endp*/false,
			  acceptor,chimera->distance);
    /* Alignment info */
    printf("\t");
    print_alignment_info(/*nblocks*/2,score,geneprob);

    /* Pairing info */
    if (hit5 != NULL && hit3 != NULL) {
      printf("\t");
      print_pair_info(hit5,hit3,pairtype,pairedlength,pairscore);
    }
    printf("\n");

    printf(",");
    Substring_print_acceptor(acceptor,/*sensep*/false,/*invertp*/true,queryseq,
			     chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			     acceptor_typeint,splicesites_divint_crosstable,/*endp*/false,
			     donor,chimera->distance);
    printf("\n");
  }

  return;
}


static void
print_terminal (T this, int score, int geneprob, IIT_T chromosome_iit, Sequence_T queryseq,
		bool invertp, IIT_T snps_iit, int *snps_divint_crosstable,
		T hit5, T hit3, Pairtype_T pairtype, int pairedlength, int pairscore) {
  char *chr;
  bool allocp;

  chr = IIT_label(chromosome_iit,this->chrnum,&allocp);

  printf(" ");
  Substring_print_terminal(this->substring1,invertp,queryseq,chr,
			   snps_iit,snps_divint_crosstable);

  /* Alignment info */
  printf("\t");
  print_alignment_info(/*nblocks*/1,score,geneprob);

  /* Pairing info */
  if (hit5 != NULL && hit3 != NULL) {
    printf("\t");
    print_pair_info(hit5,hit3,pairtype,pairedlength,pairscore);
  }

  printf("\n");

  if (allocp == true) {
    FREE(chr);
  }

  return;
}



void
Stage3_print (T this, int score, int geneprob, Genome_T genome, IIT_T chromosome_iit, Sequence_T queryseq,
	      IIT_T snps_iit, int *snps_divint_crosstable,
	      IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
	      int *splicesites_divint_crosstable, bool invertp,
	      T hit5, T hit3, Pairtype_T pairtype, int pairedlength, int pairscore) {

  if (this->hittype == EXACT || this->hittype == SUB) {
    print_single(this,score,geneprob,chromosome_iit,queryseq,invertp,
		 snps_iit,snps_divint_crosstable,hit5,hit3,pairtype,pairedlength,pairscore);
  } else if (this->hittype == INS) {
    print_insertion(this,score,geneprob,chromosome_iit,queryseq,invertp,
		    snps_iit,snps_divint_crosstable,hit5,hit3,pairtype,pairedlength,pairscore);
  } else if (this->hittype == DEL) {
    print_deletion(this,score,geneprob,chromosome_iit,queryseq,invertp,
		   snps_iit,snps_divint_crosstable,hit5,hit3,pairtype,pairedlength,pairscore);
  } else if (this->hittype == SPLICE) {
    print_splice(this,score,geneprob,genome,chromosome_iit,queryseq,invertp,
		 snps_iit,snps_divint_crosstable,splicesites_iit,
		 donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		 hit5,hit3,pairtype,pairedlength,pairscore);
  } else if (this->hittype == TERMINAL) {
    print_terminal(this,score,geneprob,chromosome_iit,queryseq,invertp,
		   snps_iit,snps_divint_crosstable,hit5,hit3,pairtype,pairedlength,pairscore);
  }
  return;
}


void
Stage3_geneprob_eval (List_T list, IIT_T geneprob_iit, IIT_T chromosome_iit, Sequence_T queryseq) {
  List_T p;
  T this;

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    if (this->hittype == EXACT || this->hittype == SUB) {
      this->geneprob = Substring_geneprob_eval_single(this->substring1,geneprob_iit,chromosome_iit);
    } else if (this->hittype == INS) {
      this->geneprob = Substring_geneprob_eval_insertion(this->substring1,this->substring2,this->indel_pos,this->nindels,
							 geneprob_iit,chromosome_iit,queryseq);
    } else if (this->hittype == DEL) {
      this->geneprob = Substring_geneprob_eval_deletion(this->substring1,this->substring2,this->indel_pos,this->nindels,
							geneprob_iit,chromosome_iit,queryseq);
    } else if (this->hittype == SPLICE) {
      this->geneprob = Substring_geneprob_eval_splice(/*donor*/this->substring1,/*acceptor*/this->substring2,
						      geneprob_iit,chromosome_iit,Sequence_fulllength(queryseq),
						      /*sensep*/this->sensedir==SENSE_FORWARD);
    }
  }

  return;
}

void
Stage3pair_geneprob_eval (List_T list, IIT_T geneprob_iit, IIT_T chromosome_iit,
			  Sequence_T queryseq5, Sequence_T queryseq3) {
  List_T p;
  Stage3pair_T pair;
  T this;

  for (p = list; p != NULL; p = List_next(p)) {
    pair = (Stage3pair_T) List_head(p);
    this = pair->hit5;
    if (this->hittype == EXACT || this->hittype == SUB) {
      this->geneprob = Substring_geneprob_eval_single(this->substring1,geneprob_iit,chromosome_iit);
    } else if (this->hittype == INS) {
      this->geneprob = Substring_geneprob_eval_insertion(this->substring1,this->substring2,this->indel_pos,this->nindels,
							 geneprob_iit,chromosome_iit,queryseq5);
    } else if (this->hittype == DEL) {
      this->geneprob = Substring_geneprob_eval_deletion(this->substring1,this->substring2,this->indel_pos,this->nindels,
							geneprob_iit,chromosome_iit,queryseq5);
    } else if (this->hittype == SPLICE) {
      this->geneprob = Substring_geneprob_eval_splice(/*donor*/this->substring1,/*acceptor*/this->substring2,
						      geneprob_iit,chromosome_iit,Sequence_fulllength(queryseq5),
						      /*sensep*/this->sensedir==SENSE_FORWARD);
    }

    this = pair->hit3;
    if (this->hittype == EXACT || this->hittype == SUB) {
      this->geneprob = Substring_geneprob_eval_single(this->substring1,geneprob_iit,chromosome_iit);
    } else if (this->hittype == INS) {
      this->geneprob = Substring_geneprob_eval_insertion(this->substring1,this->substring2,this->indel_pos,this->nindels,
							 geneprob_iit,chromosome_iit,queryseq3);
    } else if (this->hittype == DEL) {
      this->geneprob = Substring_geneprob_eval_deletion(this->substring1,this->substring2,this->indel_pos,this->nindels,
							geneprob_iit,chromosome_iit,queryseq3);
    } else if (this->hittype == SPLICE) {
      this->geneprob = Substring_geneprob_eval_splice(/*donor*/this->substring1,/*acceptor*/this->substring2,
						      geneprob_iit,chromosome_iit,Sequence_fulllength(queryseq3),
						      /*sensep*/this->sensedir==SENSE_FORWARD);
    }
  }

  return;
}


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



static bool
result_failed_p (Result_T result, bool fivep, int maxpaths, bool quiet_if_excessive_p) {
  Resulttype_T resulttype;
  int npaths;

  resulttype = Result_resulttype(result);
  if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR_SINGLE || resulttype == PAIREDEND_SAMECHR_MULTIPLE) {
    /* stage3pairarray = (Stage3pair_T *) */ Result_array(&npaths,result);

    if (npaths == 0 || (quiet_if_excessive_p && npaths > maxpaths)) {
      return true;
    } else {
      return false;
    }

  } else {
    /* PAIREDEND_AS_SINGLES or PAIREDEND_AS_SINGLES_UNIQUE */
    if (fivep == true) {
      /* stage3array = (T *) */ Result_array(&npaths,result);
    } else {
      /* stage3array = (T *) */ Result_array2(&npaths,result);
    }
      
    if (npaths == 0 || (quiet_if_excessive_p && npaths > maxpaths)) {
      return true;
    } else {
      return false;
    }
  }

}


static void
print_paired_hits (Result_T result, char initchar, bool fivep, 
		   Genome_T genome, IIT_T chromosome_iit,
		   Sequence_T queryseq, Sequence_T headerseq, Genomicpos_T pairmax,
		   IIT_T snps_iit, int *snps_divint_crosstable,
		   IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		   int *splicesites_divint_crosstable, int maxpaths, bool quiet_if_excessive_p,
		   bool sequence_inverted_p, bool invertp) {
  Stage3pair_T *stage3pairarray, stage3pair;
  Resulttype_T resulttype;
  T *stage3array, this, hit5, hit3;
  int npaths, pathnum;
  bool double_transloc_p, single_transloc_p;

  printf("%c",initchar);
  if (sequence_inverted_p == true) {
    invertp = (invertp == true) ? false : true;
  }
  if (invertp == false) {
    Sequence_print_oneline(stdout,queryseq);
  } else {
    Sequence_print_oneline_revcomp(stdout,queryseq);
  }

  resulttype = Result_resulttype(result);
  if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR_SINGLE || resulttype == PAIREDEND_SAMECHR_MULTIPLE) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);
    if (resulttype == PAIREDEND_CONCORDANT) {
      printf("\t%d %s\t",npaths,CONCORDANT_TEXT);
    } else if (resulttype == PAIREDEND_SAMECHR_SINGLE) {
      stage3pair = stage3pairarray[0];
      hit5 = stage3pair->hit5;
      hit3 = stage3pair->hit3;
      printf("\t%d %s\t",npaths,samechr_type_text(hit5->genomicstart,hit5->plusp,hit3->genomicstart,hit3->plusp,
						  stage3pair->pairedlength,pairmax));
    } else {
      printf("\t%d %s\t",npaths,UNPAIRED_TEXT);
    }
    Sequence_print_header(headerseq,/*checksump*/false);
    /* printf("\n"); -- included in header */

    if (npaths == 0 || (quiet_if_excessive_p && npaths > maxpaths)) {
      /* No output */
    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;

	if (fivep == true) {
	  Stage3_print(hit5,hit5->score,hit5->geneprob+hit3->geneprob,
		       genome,chromosome_iit,queryseq,
		       snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       invertp,hit5,hit3,stage3pair->pairtype,stage3pair->pairedlength,
		       stage3pair->score);
	} else {
	  Stage3_print(hit3,hit3->score,hit5->geneprob+hit3->geneprob,
		       genome,chromosome_iit,queryseq,
		       snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       invertp,hit5,hit3,stage3pair->pairtype,stage3pair->pairedlength,
		       stage3pair->score);
	}
      }
    }

  } else {
    /* PAIREDEND_AS_SINGLES or PAIREDEND_AS_SINGLES_UNIQUE */
    if (fivep == true) {
      stage3array = (T *) Result_array(&npaths,result);
    } else {
      stage3array = (T *) Result_array2(&npaths,result);
    }
    printf("\t%d %s\t",npaths,UNPAIRED_TEXT);
    Sequence_print_header(headerseq,/*checksump*/false);
    /* printf("\n"); -- included in header */
      
    if (npaths == 0 || (quiet_if_excessive_p && npaths > maxpaths)) {
      /* No output */
    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	this = stage3array[pathnum-1];
	Stage3_print(this,this->score,this->geneprob,genome,chromosome_iit,queryseq,
		     snps_iit,snps_divint_crosstable,
		     splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		     invertp,/*hit5*/(T) NULL,/*hit3*/(T) NULL,/*pairtype*/UNPAIRED,/*pairedlength*/0,
		     /*pairscore*/0);
      }
    }
  }

  printf("\n");

  return;
}


#if 0
static void
print_paired_hits_old (Result_T result, char initchar, bool fivep, 
		       Genome_T genome, IIT_T chromosome_iit,
		       Sequence_T queryseq, Sequence_T headerseq, Genomicpos_T pairmax,
		       IIT_T snps_iit, int *snps_divint_crosstable,
		       IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		       int *splicesites_divint_crosstable, int maxpaths, bool quiet_if_excessive_p,
		       bool sequence_inverted_p, bool invertp) {
  Stage3pair_T *stage3pairarray, stage3pair;
  Resulttype_T resulttype;
  Pairtype_T pairtype;
  T *stage3array, *stage3array1, *stage3array2, this, hit5, hit3;
  int npaths, pathnum;
  bool samechrp, double_transloc_p, single_transloc_p;
  Genomicpos_T pairedlength, low5, high5, low3, high3;

  printf("%c",initchar);
  if (sequence_inverted_p == true) {
    invertp = (invertp == true) ? false : true;
  }
  if (invertp == false) {
    Sequence_print_oneline(stdout,queryseq);
  } else {
    Sequence_print_oneline_revcomp(stdout,queryseq);
  }

  resulttype = Result_resulttype(result);
  if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);

    if (resulttype == PAIREDEND_CONCORDANT) {
      printf("\t%d %s\t",npaths,CONCORDANT_TEXT);
      samechrp = false;
    } else {
      printf("\t%d %s\t",npaths,SAMECHR_TEXT);
      samechrp = true;
    }
    Sequence_print_header(headerseq,/*checksump*/false);
    /* printf("\n"); -- included in header */

    if (npaths == 0 || (quiet_if_excessive_p && npaths > maxpaths)) {
      /* No output */
    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	hit5 = stage3pair->hit5;
	hit3 = stage3pair->hit3;

	if (fivep == true) {
	  Stage3_print(hit5,hit5->score,hit5->geneprob+hit3->geneprob,
		       genome,chromosome_iit,queryseq,
		       snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       invertp,hit5,hit3,stage3pair->pairtype,stage3pair->pairedlength,
		       stage3pair->score);
	} else {
	  Stage3_print(hit3,hit3->score,hit5->geneprob+hit3->geneprob,
		       genome,chromosome_iit,queryseq,
		       snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       invertp,hit5,hit3,stage3pair->pairtype,stage3pair->pairedlength,
		       stage3pair->score);
	}
      }
    }

  } else if (resulttype == PAIREDEND_AS_SINGLES || resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
    samechrp = double_transloc_p = single_transloc_p = false;
    if (resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
      stage3array1 = (T *) Result_array(&npaths,result);
      stage3array2 = (T *) Result_array2(&npaths,result);
      hit5 = stage3array1[0];
      hit3 = stage3array2[0];

      if (hit5->chrnum == hit3->chrnum) {
	if (hit5->chrnum == 0) {
	  double_transloc_p = true;
	} else {
	  samechrp = true;
	}
      } else if (hit5->chrnum == 0 || hit3->chrnum == 0) {
	single_transloc_p = true;
      }
    }


    if (double_transloc_p == true) {
      /* Unique and double translocation */

      if (Substring_chrnum(hit5->substring1) != Substring_chrnum(hit3->substring1) ||
	  Substring_chrnum(hit5->substring2) != Substring_chrnum(hit3->substring2)) {
	pairedlength = 0U;
	pairtype = UNPAIRED;
	printf("\t%d unpaired\t",npaths);

      } else {
	if ((pairedlength = compute_pairlength(Substring_genomicstart(hit5->substring1),
					       Substring_genomicend(hit5->substring1),
					       Substring_genomicstart(hit3->substring1),
					       Substring_genomicend(hit3->substring1))) < pairmax ||
	    (pairedlength = compute_pairlength(Substring_genomicstart(hit5->substring2),
					       Substring_genomicend(hit5->substring2),
					       Substring_genomicstart(hit3->substring2),
					       Substring_genomicend(hit3->substring2))) < pairmax) {
	  pairtype = CONCORDANT;
	  printf("\t%d concordant\t",npaths);
	} else {
	  pairtype = SAMECHR;
	  printf("\t%d %s\t",npaths,samechr_type_text(Substring_genomicstart(hit5->substring1),Substring_plusp(hit5->substring1),
						      Substring_genomicstart(hit3->substring1),Substring_plusp(hit3->substring1)));
	}
      }

      Sequence_print_header(headerseq,/*checksump*/false);

      if (fivep == true) {
	Stage3_print(hit5,hit5->score,hit5->geneprob+hit3->geneprob,
		     genome,chromosome_iit,queryseq,
		     snps_iit,snps_divint_crosstable,
		     splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		     invertp,hit5,hit3,pairtype,pairedlength,
		     hit5->score+hit3->score);
      } else {
	Stage3_print(hit3,hit3->score,hit5->geneprob+hit3->geneprob,
		     genome,chromosome_iit,queryseq,
		     snps_iit,snps_divint_crosstable,
		     splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		     invertp,hit5,hit3,pairtype,pairedlength,
		     hit5->score+hit3->score);
      }

    } else if (single_transloc_p == true) {
      /* Unique and single translocation */

      if (hit5->chrnum == 0) {
	if (Substring_chrnum(hit5->substring1) == hit3->chrnum) {
	  if ((pairedlength = compute_pairlength(Substring_genomicstart(hit5->substring1),Substring_genomicend(hit5->substring1),
						 hit3->genomicstart,hit3->genomicend)) < pairmax) {
	    pairtype = CONCORDANT;
	    printf("\t%d concordant\t",npaths);
	  } else {
	    pairtype = SAMECHR;
	    printf("\t%d %s\t",npaths,samechr_type_text(Substring_genomicstart(hit5->substring1),Substring_plusp(hit5->substring1),
							hit3->genomicstart,hit3->plusp));
	  }

	} else if (Substring_chrnum(hit5->substring2) == hit3->chrnum) {
	  if ((pairedlength = compute_pairlength(Substring_genomicstart(hit5->substring2),Substring_genomicend(hit5->substring2),
						 hit3->genomicstart,hit3->genomicend)) < pairmax) {
	    pairtype = CONCORDANT;
	    printf("\t%d concordant\t",npaths);
	  } else {
	    pairtype = SAMECHR;
	    printf("\t%d %s\t",npaths,samechr_type_text(Substring_genomicstart(hit5->substring2),Substring_plusp(hit5->substring2),
							hit3->genomicstart,hit3->plusp));
	  }

	} else {
	  pairedlength = 0U;
	  pairtype = UNPAIRED;
	  printf("\t%d unpaired\t",npaths);
	}

      } else {
	if (Substring_chrnum(hit3->substring1) == hit5->chrnum) {
	  if ((pairedlength = compute_pairlength(hit5->genomicstart,hit5->genomicend,
						 Substring_genomicstart(hit3->substring1),
						 Substring_genomicend(hit3->substring1))) < pairmax) {
	    pairtype = CONCORDANT;
	    printf("\t%d concordant\t",npaths);
	  } else {
	    pairtype = SAMECHR;
	    printf("\t%d %s\t",npaths,samechr_type_text(hit5->genomicstart,hit5->plusp,
							Substring_genomicstart(hit3->substring1),Substring_plusp(hit3->substring1)));
	  }

	} else if (Substring_chrnum(hit3->substring2) == hit5->chrnum) {
	  if ((pairedlength = compute_pairlength(hit5->genomicstart,hit5->genomicend,
						 Substring_genomicstart(hit3->substring2),
						 Substring_genomicend(hit3->substring2))) < pairmax) {
	    printf("\t%d concordant\t",npaths);
	    pairtype = CONCORDANT;
	  } else {
	    pairtype = SAMECHR;
	    printf("\t%d %s\t",npaths,samechr_type_text(hit5->genomicstart,hit5->plusp,
							Substring_genomicstart(hit3->substring2),Substring_plusp(hit3->substring2)));
	  }

	} else {
	  pairedlength = 0U;
	  pairtype = UNPAIRED;
	  printf("\t%d unpaired\t",npaths);
	}
      }

      Sequence_print_header(headerseq,/*checksump*/false);

      if (fivep == true) {
	Stage3_print(hit5,hit5->score,hit5->geneprob+hit3->geneprob,
		     genome,chromosome_iit,queryseq,
		     snps_iit,snps_divint_crosstable,
		     splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		     invertp,hit5,hit3,pairtype,pairedlength,
		     hit5->score+hit3->score);
      } else {
	Stage3_print(hit3,hit3->score,hit5->geneprob+hit3->geneprob,
		     genome,chromosome_iit,queryseq,
		     snps_iit,snps_divint_crosstable,
		     splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		     invertp,hit5,hit3,pairtype,pairedlength,
		     hit5->score+hit3->score);
      }

    } else if (samechrp == true) {
      /* Unique and same chromosome */
      printf("\t%d %s\t",npaths,samechr_type_text(hit5->genomicstart,hit5->plusp,hit3->genomicstart,hit3->plusp));
      Sequence_print_header(headerseq,/*checksump*/false);
      /* printf("\n"); -- included in header */

      pairedlength = compute_pairlength(hit5->genomicstart,hit5->genomicend,
					hit3->genomicstart,hit3->genomicend);
      if (fivep == true) {
	Stage3_print(hit5,hit5->score,hit5->geneprob+hit3->geneprob,
		     genome,chromosome_iit,queryseq,
		     snps_iit,snps_divint_crosstable,
		     splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		     invertp,hit5,hit3,/*pairtype*/SAMECHR,pairedlength,
		     hit5->score+hit3->score);
      } else {
	Stage3_print(hit3,hit3->score,hit5->geneprob+hit3->geneprob,
		     genome,chromosome_iit,queryseq,
		     snps_iit,snps_divint_crosstable,
		     splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		     invertp,hit5,hit3,/*pairtype*/SAMECHR,pairedlength,
		     hit5->score+hit3->score);
      }

    } else {
      /* Unpaired */
      if (fivep == true) {
	stage3array = (T *) Result_array(&npaths,result);
      } else {
	stage3array = (T *) Result_array2(&npaths,result);
      }
      printf("\t%d %s\t",npaths,UNPAIRED_TEXT);
      Sequence_print_header(headerseq,/*checksump*/false);
      /* printf("\n"); -- included in header */
      
      if (npaths == 0 || (quiet_if_excessive_p && npaths > maxpaths)) {
	/* No output */
      } else {
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	  this = stage3array[pathnum-1];
	  Stage3_print(this,this->score,this->geneprob,genome,chromosome_iit,queryseq,
		       snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       invertp,/*hit5*/(T) NULL,/*hit3*/(T) NULL,/*pairtype*/UNPAIRED,/*pairedlength*/0,
		       /*pairscore*/0);
	}
      }

    }

  } else {
    /* abort() */
  }

   printf("\n");

   return;
 }
 #endif



 /* Cannot put back output in fastq format */
 static void
 print_query_pairedend (Sequence_T queryseq1, Sequence_T queryseq2) {
   printf(">");
   Sequence_print_header(queryseq1,/*checksump*/false);
   /* printf("\n"); -- included in header */
   Sequence_print_oneline(stdout,queryseq1);
   printf("\n");
   Sequence_print_oneline_revcomp(stdout,queryseq2);
   printf("\n");

   return;
 }



 void
 Stage3_print_paired (Result_T result, Genome_T genome, IIT_T chromosome_iit, Sequence_T queryseq1, Sequence_T queryseq2,
		      Genomicpos_T pairmax, IIT_T snps_iit, int *snps_divint_crosstable,
		      IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		      int *splicesites_divint_crosstable, int maxpaths, bool quiet_if_excessive_p,
		      bool circularp, bool invertp, bool nofailsp, bool failsonlyp) {
   if (result_failed_p(result,/*fivep*/true,maxpaths,quiet_if_excessive_p) == true &&
       result_failed_p(result,/*fivep*/false,maxpaths,quiet_if_excessive_p) == true) {
     if (nofailsp == true) {
       /* No output */
       return;
     } else if (failsonlyp == true) {
       print_query_pairedend(queryseq1,queryseq2);
       return;
     }
   } else {
     if (failsonlyp == true) {
       /* No output */
       return;
     }
   }

   if (circularp == false) {
     /* First sequence */
     print_paired_hits(result,'>',/*fivep*/true,
		       genome,chromosome_iit,queryseq1,/*headerseq*/queryseq1,
		       pairmax,snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       maxpaths,quiet_if_excessive_p,/*sequence_inverted_p*/false,/*invertp*/false);

     /* Second sequence.  query2 has already been revcomp'd */
     print_paired_hits(result,'<',/*fivep*/false,
		       genome,chromosome_iit,queryseq2,/*headerseq*/queryseq1,
		       pairmax,snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       maxpaths,quiet_if_excessive_p,/*sequence_inverted_p*/true,/*invertp*/invertp);

   } else {
     /* First sequence.  query2 has already been revcomp'd */
     print_paired_hits(result,'>',/*fivep*/false,
		       genome,chromosome_iit,queryseq2,/*headerseq*/queryseq2,
		       pairmax,snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       maxpaths,quiet_if_excessive_p,/*sequence_inverted_p*/true,/*invertp*/(invertp == true) ? false : true);

     /* Second sequence */
     print_paired_hits(result,'<',/*fivep*/true,
		       genome,chromosome_iit,queryseq1,/*headerseq*/queryseq2,
		       pairmax,snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       maxpaths,quiet_if_excessive_p,/*sequence_inverted_p*/false,/*invertp*/true);
   }

   return;
 }


 #if 0
 /* Note: doesn't handle translocations, which are handled in gsnap.c as individual runs on each single-end */
 bool
 Stage3_pairable_p (T hit5, T hit3, Genomicpos_T pairmax, Pairedresult_T pairedresult) {

   debug(printf("Comparing %s%u (chr %d) with %s%u (chr %d): ",
		hit3->plusp ? "+" : "-",hit5->genomicstart,hit5->chrnum,
		hit3->plusp ? "+" : "-",hit3->genomicstart,hit3->chrnum));

   if (hit5->chrnum != hit3->chrnum) {
     debug(printf("No, different chromosomes\n\n"));
     return false;
   } else if (hit5->plusp != hit3->plusp) {
     if (pairedresult == INVERSION) {
       debug(printf("Yes, saw inversion and pairedresult is inversion\n\n"));
       return true;
     } else {
       debug(printf("No, saw inversion and pairedresult not inversion\n\n"));
       return false;
     }
   } else if (pairedresult == INVERSION) {
     debug(printf("Yes, didn't see inversion and pairedresult is inversion\n\n"));
     return false;
   } else if (pairedresult == SCRAMBLE) {
     debug(printf("No, didn't see inversion and pairedresult is scramble\n\n"));
     return true;
   } else if (hit5->plusp == true) {
     if (hit3->genomicstart < hit5->genomicstart) {
       debug(printf("No, positions wrong\n"));
       return false;
     } else if (pairedresult == INSERTION_DELETION) {
       debug(printf("Yes, pairedresult is insertion_deletion\n\n"));
       return true;
     } else if (hit3->genomicstart > hit5->genomicstart + pairmax) {
       debug(printf("No, length %u too long\n",hit3->genomicstart - hit5->genomicstart + 1));
       return false;
     } else {
       debug(printf("Yes, length %u is okay\n",hit3->genomicstart - hit5->genomicstart + 1));
       return true;
     }
   } else {
     if (hit5->genomicstart < hit3->genomicstart) {
       debug(printf("No, positions wrong\n"));
       return false;
     } else if (pairedresult == INSERTION_DELETION) {
       debug(printf("Yes, pairedresult is insertion_deletion\n\n"));
       return true;
     } else if (hit5->genomicstart > hit3->genomicstart + pairmax) {
       debug(printf("No, length %u too long\n",hit5->genomicstart - hit3->genomicstart + 1));
       return false;
     } else {
       debug(printf("Yes, length %u is okay\n",hit5->genomicstart - hit3->genomicstart + 1));
       return true;
     }
   }
 }
 #endif

 static Stage3pair_T
 Stage3pair_new (T hit5, T hit3, int querylength5, int querylength3, int expected_pairlength,
		 Pairtype_T pairtype) {
   Stage3pair_T new = (Stage3pair_T) MALLOC(sizeof(*new));
   Genomicpos_T pairlength5, pairlength3;

   new->hit5 = Stage3_copy(hit5);
   new->hit3 = Stage3_copy(hit3);
   new->pairtype = pairtype;

   new->score = hit5->score + hit3->score;
   new->nmatches = hit5->nmatches + hit3->nmatches;

   if (pairtype == SAMECHR) {
     new->dir = 0;
     new->low = hit5->genomicstart;
     if (hit5->genomicend < new->low) {
       new->low = hit5->genomicend;
     }
     if (hit3->genomicstart < new->low) {
       new->low = hit3->genomicstart;
     }
     if (hit3->genomicend < new->low) {
       new->low = hit3->genomicend;
     }

     new->high = hit5->genomicstart;
     if (hit5->genomicend > new->high) {
       new->high = hit5->genomicend;
     }
     if (hit3->genomicstart > new->high) {
       new->high = hit3->genomicstart;
     }
     if (hit3->genomicend > new->high) {
       new->high = hit3->genomicend;
     }

     new->pairedlength = new->high - new->low;

   } else if (pairtype == CONCORDANT) {
     if (hit5->chrnum == 0 && hit3->chrnum == 0) {
       new->dir = 0;

       new->low = Substring_genomicstart(hit5->substring1);
       if (Substring_genomicend(hit5->substring1) < new->low) {
	 new->low = Substring_genomicend(hit5->substring1);
       }
       if (Substring_genomicstart(hit5->substring2) < new->low) {
	 new->low = Substring_genomicstart(hit5->substring2);
       }
       if (Substring_genomicend(hit5->substring2) < new->low) {
	 new->low = Substring_genomicend(hit5->substring2);
       }
       if (Substring_genomicstart(hit3->substring1) < new->low) {
	 new->low = Substring_genomicstart(hit3->substring1);
       }
       if (Substring_genomicend(hit3->substring1) < new->low) {
	 new->low = Substring_genomicend(hit3->substring1);
       }
       if (Substring_genomicstart(hit3->substring2) < new->low) {
	 new->low = Substring_genomicstart(hit3->substring2);
       }
       if (Substring_genomicend(hit3->substring2) < new->low) {
	 new->low = Substring_genomicend(hit3->substring2);
       }

       new->high = Substring_genomicstart(hit5->substring1);
       if (Substring_genomicend(hit5->substring1) > new->high) {
	 new->high = Substring_genomicend(hit5->substring1);
       }
       if (Substring_genomicstart(hit5->substring2) > new->high) {
	 new->high = Substring_genomicstart(hit5->substring2);
       }
       if (Substring_genomicend(hit5->substring2) > new->high) {
	 new->high = Substring_genomicend(hit5->substring2);
       }
       if (Substring_genomicstart(hit3->substring1) > new->high) {
	 new->high = Substring_genomicstart(hit3->substring1);
       }
       if (Substring_genomicend(hit3->substring1) > new->high) {
	 new->high = Substring_genomicend(hit3->substring1);
       }
       if (Substring_genomicstart(hit3->substring2) > new->high) {
	 new->high = Substring_genomicstart(hit3->substring2);
       }
       if (Substring_genomicend(hit3->substring2) > new->high) {
	 new->high = Substring_genomicend(hit3->substring2);
       }

       if (Substring_chrnum(hit5->substring1) == Substring_chrnum(hit3->substring1) &&
	   Substring_plusp(hit5->substring1) == Substring_plusp(hit3->substring1)) {
	 if (Substring_plusp(hit5->substring1) == true) {
	   pairlength5 = Substring_genomicend(hit3->substring1) - (Substring_genomicend(hit5->substring1) - querylength5);
	   pairlength3 = (Substring_genomicstart(hit3->substring1) + querylength3) - Substring_genomicstart(hit5->substring1);
	 } else {
	   pairlength5 = (Substring_genomicend(hit5->substring1) + querylength5) - Substring_genomicend(hit3->substring1);
	   pairlength3 = Substring_genomicstart(hit5->substring1) - (Substring_genomicstart(hit3->substring1) - querylength3);
	 }

       } else if (Substring_chrnum(hit5->substring1) == Substring_chrnum(hit3->substring2) &&
		  Substring_plusp(hit5->substring1) == Substring_plusp(hit3->substring2)) {
	 if (Substring_plusp(hit5->substring1) == true) {
	   pairlength5 = Substring_genomicend(hit3->substring2) - (Substring_genomicend(hit5->substring1) - querylength5);
	   pairlength3 = (Substring_genomicstart(hit3->substring2) + querylength3) - Substring_genomicstart(hit5->substring1);
	 } else {
	   pairlength5 = (Substring_genomicend(hit5->substring1) + querylength5) - Substring_genomicend(hit3->substring2);
	   pairlength3 = Substring_genomicstart(hit5->substring1) - (Substring_genomicstart(hit3->substring2) - querylength3);
	 }

       } else if (Substring_chrnum(hit5->substring2) == Substring_chrnum(hit3->substring1) &&
		  Substring_plusp(hit5->substring2) == Substring_plusp(hit3->substring1)) {
	 if (Substring_plusp(hit5->substring2) == true) {
	   pairlength5 = Substring_genomicend(hit3->substring1) - (Substring_genomicend(hit5->substring2) - querylength5);
	   pairlength3 = (Substring_genomicstart(hit3->substring1) + querylength3) - Substring_genomicstart(hit5->substring2);
	 } else {
	   pairlength5 = (Substring_genomicend(hit5->substring2) + querylength5) - Substring_genomicend(hit3->substring1);
	   pairlength3 = Substring_genomicstart(hit5->substring2) - (Substring_genomicstart(hit3->substring1) - querylength3);
	 }

       } else if (Substring_chrnum(hit5->substring2) == Substring_chrnum(hit3->substring2) &&
		  Substring_plusp(hit5->substring2) == Substring_plusp(hit3->substring2)) {
	 if (Substring_plusp(hit5->substring2) == true) {
	   pairlength5 = Substring_genomicend(hit3->substring2) - (Substring_genomicend(hit5->substring2) - querylength5);
	   pairlength3 = (Substring_genomicstart(hit3->substring2) + querylength3) - Substring_genomicstart(hit5->substring2);
	 } else {
	   pairlength5 = (Substring_genomicend(hit5->substring2) + querylength5) - Substring_genomicend(hit3->substring2);
	   pairlength3 = Substring_genomicstart(hit5->substring2) - (Substring_genomicstart(hit3->substring2) - querylength3);
	 }

       } else {
	 pairlength5 = pairlength3 = 0;
       }

       if (pairlength5 < pairlength3) {
	 new->pairedlength = pairlength5;
       } else {
	 new->pairedlength = pairlength3;
       }

     } else if (hit5->chrnum == 0) {
       /* Use hit3 */
       if (hit3->plusp == true /* implies hit5->plusp is true */) {
	 new->dir = +1;
	 new->low = hit5->genomicstart;
	 new->high = hit3->genomicend;
       } else {
	 new->dir = -1;
	 new->low = hit3->genomicend;
	 new->high = hit5->genomicstart;
       }
       pairlength5 = hit3->genomicend - (hit5->genomicend - querylength5);
       pairlength3 = (hit3->genomicstart + querylength3) - hit5->genomicstart;
       if (pairlength5 < pairlength3) {
	 new->pairedlength = pairlength5;
       } else {
	 new->pairedlength = pairlength3;
       }

     } else if (hit5->plusp == true) {
       new->dir = +1;
       new->low = hit5->genomicstart;
       new->high = hit3->genomicend;

       pairlength5 = hit3->genomicend - (hit5->genomicend - querylength5);
       pairlength3 = (hit3->genomicstart + querylength3) - hit5->genomicstart;
       if (pairlength5 < pairlength3) {
	 new->pairedlength = pairlength5;
       } else {
	 new->pairedlength = pairlength3;
       }

     } else {
       new->dir = -1;
       new->low = hit3->genomicend;
       new->high = hit5->genomicstart;

       pairlength5 = (hit5->genomicend + querylength5) - hit3->genomicend;
       pairlength3 = hit5->genomicstart - (hit3->genomicstart - querylength3);
       if (pairlength5 < pairlength3) {
	 new->pairedlength = pairlength5;
       } else {
	 new->pairedlength = pairlength3;
       }
     }

   } else {
     abort();
   }

   if (new->pairedlength < expected_pairlength) {
     new->absdifflength = expected_pairlength - new->pairedlength;
   } else {
     new->absdifflength = new->pairedlength - expected_pairlength;
   }

   if (SENSE_CONSISTENT_P(hit5->sensedir,hit3->sensedir)) {
     new->sense_consistent_p = true;
   } else {
     new->sense_consistent_p = false;
   }

   return new;
 }

 List_T
 Stage3_filter_bymatch (List_T hitlist) {
   List_T filtered = NULL, p;
   T hit;
   int min_total_nmismatches = 1000;

   for (p = hitlist; p != NULL; p = p->rest) {
     hit = (T) p->first;
     if (hit->total_nmismatches < min_total_nmismatches) {
       min_total_nmismatches = hit->total_nmismatches;
     }
   }

   for (p = hitlist; p != NULL; p = p->rest) {
     hit = (T) p->first;
     if (hit->total_nmismatches == min_total_nmismatches) {
       filtered = List_push(filtered,hit);
     } else {
       Stage3_free(&hit);
     }
   }
   List_free(&hitlist);

   return filtered;
 }

 #if 0
 static int
 chimera_match_distance_cmp (const void *a, const void *b) {
   T x = * (T *) a;
   T y = * (T *) b;

   if (x->total_nmismatches < y->total_nmismatches) {
     return -1;
   } else if (x->total_nmismatches > y->total_nmismatches) {
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


/* Never returns NULL */
T
Stage3_new_splice (int *found_score, Substring_T donor, Substring_T acceptor, Genomicpos_T distance,
		   bool shortdistancep, int splicing_penalty, int querylength,
		   bool copy_donor_p, bool copy_acceptor_p, int sensedir) {
  T new;
  
  new = (T) MALLOC(sizeof(*new));
  new->deletion = (char *) NULL;
  new->hittype = SPLICE;
  new->nindels = 0;

  if (donor == NULL) {
    abort();

  } else if (acceptor == NULL) {
    abort();

  } else {
    new->substring1 = copy_donor_p ? Substring_copy(donor) : donor;
    new->substring2 = copy_acceptor_p ? Substring_copy(acceptor) : acceptor;

    if (shortdistancep == true) {
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

    if (sensedir == SENSE_FORWARD) {
      new->genomicstart = Substring_genomicstart(donor);
      new->genomicend = Substring_genomicend(acceptor);
    } else if (sensedir == SENSE_ANTI) {
      new->genomicstart = Substring_genomicstart(acceptor);
      new->genomicend = Substring_genomicend(donor);
    } else {
      abort();
    }


    if (new->chrnum == 0) {
      new->effective_chrnum = -1; /* not assigning here */
      new->other_chrnum = -1;	  /* not assigning here */
      new->substring_low = (Substring_T) NULL;
    
    } else {
      new->effective_chrnum = new->chrnum;
      new->other_chrnum = 0;
      if (sensedir == SENSE_FORWARD) {
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

    new->total_nmismatches = Substring_nmismatches(donor) + Substring_nmismatches(acceptor);
    new->ntscore = splicing_penalty + Substring_nmismatches(donor) + Substring_nmismatches(acceptor);
    new->score = new->ntscore + Substring_ncolordiffs(donor) + Substring_ncolordiffs(acceptor);
#if 0
    new->nmatches = Substring_match_length(donor) + Substring_match_length(acceptor) - 
      (Substring_nmismatches(donor) + Substring_ncolordiffs(donor) + Substring_nmismatches(acceptor) + Substring_ncolordiffs(acceptor));
#else
    new->nmatches = Substring_nmatches(donor) + Substring_nmatches(acceptor);
#endif
    new->chimera_prob = Substring_chimera_prob(donor) + Substring_chimera_prob(acceptor);
    new->chimera_ambiguous_p = false;
#if 0
    new->half_intron_score = 0.0;
#endif
  }
  new->geneprob = -1;

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->distance = distance;
  new->sensedir = sensedir;

  new->paired_seenp = false;
  new->concordantp = false;

  return new;
}


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
  } else if (x->low != y->low) {
    return false;
  } else if (x->high != y->high) {
    return false;
  } else {
    return true;
  }
}

static bool
hitpair_overlap (Stage3pair_T x, Stage3pair_T y) {
  if ((x->hit5->hittype == SPLICE || x->hit3->hittype == SPLICE) &&
      (y->hit5->hittype == SPLICE || y->hit3->hittype == SPLICE)) {
    /* Special case: pairs involving splices don't overlap */
    return false;
  } else if (x->dir != y->dir) {
    return false;		/* Different strands */
  } else if (x->high < y->low) {
    return false;
  } else if (x->low > y->high) {
    return false;
  } else {
    return true;
  }
}


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
  }

  hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
  List_free(&hitpairlist);

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


List_T
Stage3pair_remove_samechr (List_T hitpairlist) {
  List_T concordant = NULL, p;
  Stage3pair_T hitpair;

#if 0
  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->pairtype == SAMECHR) {
      Stage3pair_free(&hitpair);
    } else {
      concordant = List_push(concordant,hitpair);
    }
  }

  List_free(&hitpairlist);

  return concordant;
#else
  return hitpairlist;
#endif
}


List_T
Stage3pair_remove_duplicates (List_T hitpairlist) {
#ifdef DEBUG8
  List_T p;
#endif
  List_T unique = NULL;
  Stage3pair_T hitpair, *hitpairs, *prev;
  int best_nmatches;
  int nkept, n, i, j;
  bool *eliminate, used_equal_p;

  n = List_length(hitpairlist);
  debug8(printf("Entered Stage3pair_remove_duplicates with %d pairs\n",n));
  if (n == 0) {
    return NULL;
  }

  eliminate = (bool *) CALLOC(n,sizeof(bool));

  hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
  List_free(&hitpairlist);

  /* Check for exact duplicates */
  debug8(printf("Checking for exact duplicates\n"));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_position_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d: %u..%u (dir = %d), nmatches: %d\n",
		  i,hitpair->low,hitpair->high,hitpair->dir,hitpair->nmatches);
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

  j = 0;
  prev = hitpairs;
  hitpairs = (Stage3pair_T *) CALLOC(n,sizeof(Stage3pair_T));

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
  FREE(eliminate);


  /* Check for duplicates by equivalence according to overlap, and keep best nmatches */
  n = j;
  eliminate = (bool *) CALLOC(n,sizeof(bool));
  debug8(printf("Checking for duplicates among %d hitpairs by overlap\n",n));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_position_cmp);

  debug8(
	 for (i = 0; i < n; i++) {
	   hitpair = hitpairs[i];
	   printf("  Initial %d: %u-%u (dir = %d), nmatches: %d\n",
		  i,hitpair->low,hitpair->high,hitpair->dir,hitpair->nmatches);
	 }
	 );

  i = 0;
  while (i < n) {
    debug8(printf("Looking at %d\n",i));
    best_nmatches = hitpairs[i]->nmatches;

    j = i+1;
    while (j < n && hitpair_overlap(hitpairs[j],hitpairs[i]) == true) {
      debug8(printf("  %d overlaps with %d",j,i));
      if (hitpairs[j]->nmatches > best_nmatches) {
	best_nmatches = hitpairs[j]->nmatches;
	debug8(printf("=> best_nmatches %d\n",best_nmatches));
      }
      debug8(printf("\n"));
      j++;
    }

    used_equal_p = false;
    if (hitpairs[i]->nmatches == best_nmatches) {
      used_equal_p = true;
    } else if (hitpairs[i]->nmatches < best_nmatches) {
      eliminate[i] = true;
      debug8(printf("  %d: nmatches %d < nmatches %d, so marking for elimination\n",
		   i,hitpairs[i]->nmatches,best_nmatches));
    }

    j = i+1;
    while (j < n && hitpair_overlap(hitpairs[j],hitpairs[i]) == true) {
      if (hitpairs[j]->nmatches == best_nmatches) {
	if (used_equal_p == true) {
	  eliminate[j] = true;
	  debug8(printf("  %d: nmatches %d == nmatches %d, so marking for elimination\n",
		       j,hitpairs[j]->nmatches,best_nmatches));
	} else {
	  used_equal_p = true;
	}

      } else if (hitpairs[j]->nmatches < best_nmatches) {
	eliminate[j] = true;
	debug8(printf("  %d: nmatches %d < nmatches %d, so marking for elimination\n",
		     j,hitpairs[i]->nmatches,best_nmatches));
	
      }
      j++;
    }

    i = j;
  }

  for (i = 0; i < n; i++) {
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
Stage3pair_optimal_score (List_T hitpairlist, int cutoff_level, int suboptimal_mismatches) {
  List_T optimal = NULL, p;
  Stage3pair_T hitpair;
  int n;
  int minscore = MAX_QUERYLENGTH + MAX_QUERYLENGTH;
  int best_terminal_nmatches = 0;

  n = List_length(hitpairlist);
  if (n == 0) {
    return NULL;
  }

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->score <= cutoff_level) {
      if (hitpair->hit5->hittype == TERMINAL || hitpair->hit3->hittype == TERMINAL) {
	if (hitpair->nmatches > best_terminal_nmatches) {
	  best_terminal_nmatches = hitpair->nmatches;
	}
      } else if (hitpair->score < minscore) {
	minscore = hitpair->score;
      }
    }
  }

  debug(printf("Stage3pair_optimal_score over %d pairs: minscore = %d + subopt:%d, best_terminal_nmatches = %d\n",
	       n,minscore,suboptimal_mismatches,best_terminal_nmatches));
  minscore += suboptimal_mismatches;

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->score > cutoff_level) {
      debug(printf("Eliminating a hit pair with score %d > cutoff_level %d\n",
		   hitpair->score,cutoff_level));
      Stage3pair_free(&hitpair);
    } else if (hitpair->score <= minscore) {
      if (hitpair->hit5->hittype == TERMINAL || hitpair->hit3->hittype == TERMINAL) {
	if (hitpair->nmatches < best_terminal_nmatches) {
	  debug(printf("Eliminating a terminal hitpair with nmatches %d < best %d\n",
		       hitpair->nmatches,best_terminal_nmatches));
	  Stage3pair_free(&hitpair);
	} else {
	  debug(printf("Keeping a terminal hitpair with nmatches %d >= best %d\n",
		       hitpair->nmatches,best_terminal_nmatches));
	  optimal = List_push(optimal,hitpair);
	}
      } else {
	debug(printf("Keeping a hit pair with score %d\n",hitpair->score));
	optimal = List_push(optimal,hitpair);
      }
    } else {
      debug(printf("Eliminating a hit pair with score %d\n",hitpair->score));
      Stage3pair_free(&hitpair);
    }
  }
  
  List_free(&hitpairlist);

  return optimal;
}



/* Finds concordant pairs if nconcordant is 0 */
List_T
Stage3_pair_up_concordant (bool *abort_pairing_p, int *found_score, int *nconcordant, List_T hitpairs,
			   List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			   int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			   Genomicpos_T pairmax, Genomicpos_T expected_pairlength,
			   int querylength5, int querylength3, int maxpairedpaths) {
  List_T copies = NULL, q;
  T **hits5_plus, **hits5_minus, **hits3_plus, **hits3_minus, *hits5, *hits3, hit5, hit3, copy;
  int *nhits5_plus, *nhits5_minus, *nhits3_plus, *nhits3_minus, nhits5, nhits3;
  int pairscore, score5, score3, i, j, jj;
  bool *sorted5p, *sorted3p;
  
  debug5(printf("Starting Stage3_pair_up_concordant with %d concordant, narray5 %d, narray3 %d, found_score %d\n",
		*nconcordant,narray5,narray3,*found_score));

  /* Find sizes for allocating memory */
  debug5(printf("Sizes of 5-end pieces by level:\n"));
  nhits5_plus = (int *) CALLOC(cutoff_level_5+1,sizeof(int));
  nhits5_minus = (int *) CALLOC(cutoff_level_5+1,sizeof(int));
  for (i = 0; i < narray5; i++) {
    debug5(printf("  array5 level %d with %d hits\n",i,List_length(hitarray5[i])));
    for (q = hitarray5[i]; q != NULL; q = q->rest) {
      hit5 = (T) q->first;
      if (hit5->score > cutoff_level_5) {
	debug5(printf("Skipping hit with score %d > cutoff level %d\n",hit5->score,cutoff_level_5));
      } else if (hit5->chimera_ambiguous_p == true) {
	debug5(printf("Skipping hit with score %d because chimera is ambiguous\n",hit5->score));
      } else if (hit5->chrnum == 0) {
	/* Translocation: Enter under each substring chrnum */
	if (Substring_plusp(hit5->substring1)) {
	  nhits5_plus[hit5->score]++;
	} else {
	  nhits5_minus[hit5->score]++;
	}

	if (Substring_plusp(hit5->substring2)) {
	  nhits5_plus[hit5->score]++;
	} else {
	  nhits5_minus[hit5->score]++;
	}

      } else if (hit5->plusp) {
	nhits5_plus[hit5->score]++;
      } else {
	nhits5_minus[hit5->score]++;
      }
    }
  }

  debug5(
	 printf("Sizes of 5-end pieces by score:\n");
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


  debug5(printf("Sizes of 3-end pieces by level:\n"));
  nhits3_plus = (int *) CALLOC(cutoff_level_3+1,sizeof(int));
  nhits3_minus = (int *) CALLOC(cutoff_level_3+1,sizeof(int));
  for (i = 0; i < narray3; i++) {
    debug5(printf(" array3 level %d with %d hits\n",i,List_length(hitarray3[i])));
    for (q = hitarray3[i]; q != NULL; q = q->rest) {
      hit3 = (T) q->first;
      if (hit3->score > cutoff_level_3) {
	debug5(printf("Skipping hit with score %d > cutoff level %d\n",hit3->score,cutoff_level_3));
      } else if (hit3->chimera_ambiguous_p == true) {
	debug5(printf("Skipping hit with score %d because chimera is ambiguous\n",hit3->score));
      } else if (hit3->chrnum == 0) {
	/* Translocation: Enter under each substring chrnum */
	if (Substring_plusp(hit3->substring1)) {
	  nhits3_plus[hit3->score]++;
	} else {
	  nhits3_minus[hit3->score]++;
	}

	if (Substring_plusp(hit3->substring2)) {
	  nhits3_plus[hit3->score]++;
	} else {
	  nhits3_minus[hit3->score]++;
	}

      } else if (hit3->plusp) {
	nhits3_plus[hit3->score]++;
      } else {
	nhits3_minus[hit3->score]++;
      }
    }
  }
  debug5(
	 printf("Sizes of 3-end pieces by score:\n");
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
      } else if (hit5->chimera_ambiguous_p == true) {
	/* Skip */
      } else if (hit5->chrnum == 0) {
	/* Translocation: Enter under each substring chrnum */
	copy = Stage3_copy(hit5);
	copy->effective_chrnum = Substring_chrnum(hit5->substring1);
	copy->other_chrnum = Substring_chrnum(hit5->substring2);
	copy->genomicstart = Substring_genomicstart(hit5->substring1);
	copy->genomicend = Substring_genomicend(hit5->substring1);
	copy->plusp = Substring_plusp(hit5->substring1);
	copies = List_push(copies,(void *) copy);
	if (copy->plusp) {
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
	if (copy->plusp) {
	  hits5_plus[hit5->score][nhits5_plus[hit5->score]++] = copy;
	} else {
	  hits5_minus[hit5->score][nhits5_minus[hit5->score]++] = copy;
	}

      } else if (hit5->plusp) {
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
      } else if (hit3->chimera_ambiguous_p == true) {
	/* Skip */
      } else if (hit3->chrnum == 0) {
	/* Translocation: Enter under each substring chrnum */
	copy = Stage3_copy(hit3);
	copy->effective_chrnum = Substring_chrnum(hit3->substring1);
	copy->other_chrnum = Substring_chrnum(hit3->substring2);
	copy->genomicstart = Substring_genomicstart(hit3->substring1);
	copy->genomicend = Substring_genomicend(hit3->substring1);
	copy->plusp = Substring_plusp(hit3->substring1);
	copies = List_push(copies,(void *) copy);
	if (copy->plusp) {
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
	if (copy->plusp) {
	  hits3_plus[hit3->score][nhits3_plus[hit3->score]++] = copy;
	} else {
	  hits3_minus[hit3->score][nhits3_minus[hit3->score]++] = copy;
	}

      } else if (hit3->plusp) {
	hits3_plus[hit3->score][nhits3_plus[hit3->score]++] = hit3;
      } else {
	hits3_minus[hit3->score][nhits3_minus[hit3->score]++] = hit3;
      }
    }
  }


  /* Look for concordant pairs */
  sorted5p = (bool *) CALLOC(cutoff_level_5+1,sizeof(bool));
  sorted3p = (bool *) CALLOC(cutoff_level_3+1,sizeof(bool));

  pairscore = 0;
  while (*abort_pairing_p == false && pairscore <= *found_score + subopt_levels &&
	 pairscore <= cutoff_level_5 + cutoff_level_3) {
    debug5(printf("pairscore = %d\n",pairscore));
    for (score5 = 0; score5 <= pairscore; score5++) {
      debug5(printf("score5 = %d, score3 = %d\n",score5,pairscore-score5));

      if (score5 <= cutoff_level_5 && ((score3 = pairscore - score5) <= cutoff_level_3)) {
	/* Sort if necessary */
	if (sorted5p[score5] == false) {
	  if (nhits5_plus[score5] > 0) {
	    qsort(hits5_plus[score5],nhits5_plus[score5],sizeof(T),genomicstart_cmp);
	  }
	  if (nhits5_minus[score5] > 0) {
	    qsort(hits5_minus[score5],nhits5_minus[score5],sizeof(T),genomicstart_cmp);
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

	/* hits5_plus against hits3_plus (really on minus) */
	i = j = 0;
	hits5 = hits5_plus[score5];
	hits3 = hits3_plus[score3];
	nhits5 = nhits5_plus[score5];
	nhits3 = nhits3_plus[score3];

	while (*abort_pairing_p == false && i < nhits5 && j < nhits3) {
	  hit5 = hits5[i];
	  hit3 = hits3[j];

	  if (hit3->genomicstart + querylength3 < hit5->genomicstart) {
	    debug5(printf("plus: i=%d %u %s, j=%d %u %s, advance j\n",
			  i,hit5->genomicstart,print_sense(hit5->sensedir),
			  j,hit3->genomicstart,print_sense(hit3->sensedir)));
	    j++;
	  } else if (hit5->genomicstart + pairmax < hit3->genomicstart) {
	    debug5(printf("plus: i=%d %u %s, j=%d %u %s, advance i\n",
			  i,hit5->genomicstart,print_sense(hit5->sensedir),
			  j,hit3->genomicstart,print_sense(hit3->sensedir)));
	    i++;
	  } else {
	    debug5(printf("plus: i=%d %u %s, j=%d %u %s, overlap\n",
			  i,hit5->genomicstart,print_sense(hit5->sensedir),
			  j,hit3->genomicstart,print_sense(hit3->sensedir)));
	    for (jj = j; jj >= 0 && hits3[jj]->genomicstart + querylength3 >= hit5->genomicstart; jj--) {
	      debug5(printf("plus: i=%d %u %s, jj=%d %u %s, backup jj\n",
			    i,hit5->genomicstart,print_sense(hit5->sensedir),
			    jj,hits3[jj]->genomicstart,print_sense(hit3->sensedir)));
	    }
	    for (jj++; jj < nhits3 && hits3[jj]->genomicstart <= hit5->genomicstart + pairmax; jj++) {
	      debug5(printf("plus: i=%d %u %s, jj=%d %u %s",
			    i,hit5->genomicstart,print_sense(hit5->sensedir),
			    jj,hits3[jj]->genomicstart,print_sense(hit3->sensedir)));
	      hit3 = hits3[jj];

	      /* Only want pairs not previously seen */
	      if (hit5->paired_seenp == false || hit3->paired_seenp == false) {
		if (hit5->effective_chrnum != hit3->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (SENSE_INCONSISTENT_P(hit5->sensedir,hit3->sensedir)) {
		  debug5(printf(" => sense inconsistent"));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 && hit5->other_chrnum != hit3->other_chrnum) {
		  debug5(printf(" => translocations not reciprocal"));
		} else {
		  debug5(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
				hit5->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  hitpairs = List_push(hitpairs,Stage3pair_new(hit5,hit3,querylength5,querylength3,
							       expected_pairlength,/*pairtype*/CONCORDANT));
		  if (pairscore < *found_score) {
		    *found_score = pairscore;
		    debug5(printf(" => updating found_score to be %d",*found_score));
		  }
		  if (++(*nconcordant) > maxpairedpaths) {
		    debug(printf("%d concordant paths exceeds %d\n",*nconcordant,maxpairedpaths));
		    *abort_pairing_p = true;
		  }
		  hit5->paired_seenp = true;
		  hit3->paired_seenp = true;
		}
	      }
	      debug5(printf("\n"));
	    }
	    i++;
	  }
	}

	/* hits3_minus (really on plus) against hits5_minus */
	i = j = 0;
	hits3 = hits3_minus[score3];
	hits5 = hits5_minus[score5];
	nhits3 = nhits3_minus[score3];
	nhits5 = nhits5_minus[score5];

	while (*abort_pairing_p == false && j < nhits5 && i < nhits3) {
	  hit5 = hits5[j];
	  hit3 = hits3[i];

	  if (hit5->genomicstart + querylength5 < hit3->genomicstart) {
	    debug5(printf("minus: j=%d %u %s, i=%d %u %s, advance j\n",
			  j,hit5->genomicstart,print_sense(hit5->sensedir),
			  i,hit3->genomicstart,print_sense(hit3->sensedir)));
	    j++;
	  } else if (hit3->genomicstart + pairmax < hit5->genomicstart) {
	    debug5(printf("minus: j=%d %u %s, i=%d %u %s, advance i\n",
			  j,hit5->genomicstart,print_sense(hit5->sensedir),
			  i,hit3->genomicstart,print_sense(hit3->sensedir)));
	    i++;
	  } else {
	    debug5(printf("minus: j=%d %u %s, i=%d %u %s, overlap\n",
			  j,hit5->genomicstart,print_sense(hit5->sensedir),
			  i,hit3->genomicstart,print_sense(hit3->sensedir)));
	    for (jj = j; jj >= 0 && hits5[jj]->genomicstart + querylength5 >= hit3->genomicstart; jj--) {
	      debug5(printf("minus: jj=%d %u %s, i=%d %u %s, backup jj",
			    jj,hits5[jj]->genomicstart,print_sense(hit5->sensedir),
			    i,hit3->genomicstart,print_sense(hit3->sensedir)));
	    }
	    for (jj++; jj < nhits5 && hits5[jj]->genomicstart <= hit3->genomicstart + pairmax; jj++) {
	      debug5(printf("minus: jj=%d %u %s, i=%d %u %s",
			    jj,hits5[jj]->genomicstart,print_sense(hit5->sensedir),
			    i,hit3->genomicstart,print_sense(hit3->sensedir)));
	      hit5 = hits5[jj];

	      /* Only want pairs not previously seen */
	      if (hit3->paired_seenp == false || hit5->paired_seenp == false) {
		if (hit3->effective_chrnum != hit5->effective_chrnum) {
		  debug5(printf(" => diff chrs %d and %d",hit5->effective_chrnum,hit3->effective_chrnum));
		} else if (SENSE_INCONSISTENT_P(hit3->sensedir,hit5->sensedir)) {
		  debug5(printf(" => sense inconsistent"));
		} else if (hit5->chrnum == 0 && hit3->chrnum == 0 && hit5->other_chrnum != hit3->other_chrnum) {
		  debug5(printf(" => translocations not reciprocal"));
		} else {
		  debug5(printf(" => concordant effchr %d (chrnum5 %d, chrnum3 %d)",
				hit3->effective_chrnum,hit5->chrnum,hit3->chrnum));
		  hitpairs = List_push(hitpairs,Stage3pair_new(hit5,hit3,querylength5,querylength3,
							       expected_pairlength,/*pairtype*/CONCORDANT));
		  if (pairscore < *found_score) {
		    *found_score = pairscore;
		    debug5(printf(" => updating found_score to be %d",*found_score));
		  }
		  if (++(*nconcordant) > maxpairedpaths) {
		    debug(printf("%d concordant paths exceeds %d\n",*nconcordant,maxpairedpaths));
		    *abort_pairing_p = true;
		  }
		  hit3->paired_seenp = true;
		  hit5->paired_seenp = true;
		}
	      }
	      debug5(printf("\n"));
	    }
	    i++;
	  }
	}
      }
    }
    pairscore++;
  }

  FREE(sorted3p);
  FREE(sorted5p);

  for (q = copies; q != NULL; q = List_next(q)) {
    copy = (T) List_head(q);
    Stage3_free(&copy);
  }
  List_free(&copies);

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


/* Want one samechr, but could have duplicates or hitarray5 and
   hitarray3, so need to check at time of creating a result (after
   duplicates removed). */
List_T
Stage3_pair_up_samechr (int *found_score, List_T hitpairs,
			List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
			int cutoff_level_5, int cutoff_level_3, int subopt_levels,
			Genomicpos_T expected_pairlength, int querylength5, int querylength3,
			int maxpairedpaths) {
  List_T copies = NULL, q;
  T **hits5_all, **hits3_all, *hits5, *hits3, hit5, hit3, copy;
  Stage3pair_T stage3pair;
  int *nhits5_all, *nhits3_all, nhits5, nhits3;
  int pairscore, score5, score3, i, j, jj;
  bool abort_samechr_p = false;
  bool *sorted5p, *sorted3p;
  int nsamechr = 0;
  
  debug6(printf("Starting Stage3_pair_up_samechr\n"));

  /* Find sizes for allocating memory */
  debug6(printf("Sizes of 5-end pieces:\n"));
  nhits5_all = (int *) CALLOC(cutoff_level_5+1,sizeof(int));
  for (i = 0; i < narray5; i++) {
    for (q = hitarray5[i]; q != NULL; q = q->rest) {
      hit5 = (T) q->first;
      if (hit5->score > cutoff_level_5) {
      } else if (hit5->chrnum == 0) {
	/* Translocation: Enter under each substring chrnum */
	nhits5_all[hit5->score]++;
	nhits5_all[hit5->score]++;
      } else {
	nhits5_all[hit5->score]++;
      }
    }
  }

  debug6(
	 printf("Sizes of 5-end pieces:\n");
	 for (score5 = 0; score5 <= cutoff_level_5; score5++) {
	   printf("  score %d: %d\n",score5,nhits5_all[score5]);
	 }
	 );

  debug6(printf("Sizes of 3-end pieces:\n"));
  nhits3_all = (int *) CALLOC(cutoff_level_3+1,sizeof(int));
  for (i = 0; i < narray3; i++) {
    for (q = hitarray3[i]; q != NULL; q = q->rest) {
      hit3 = (T) q->first;
      if (hit3->score > cutoff_level_3) {
      } else if (hit3->chrnum == 0) {
	/* Translocation: Enter under each substring chrnum */
	nhits3_all[hit3->score]++;
	nhits3_all[hit3->score]++;
      } else {
	nhits3_all[hit3->score]++;
      }
    }
  }
  debug6(
	 printf("Sizes of 3-end pieces:\n");
	 for (score3 = 0; score3 <= cutoff_level_3; score3++) {
	   printf("  score %d: %d\n",score3,nhits3_all[score3]);
	 }
	 );

  /* Store hits5 */
  hits5_all = (T **) CALLOC(cutoff_level_5+1,sizeof(T *));
  for (score5 = 0; score5 <= cutoff_level_5; score5++) {
    if (nhits5_all[score5] == 0) {
      hits5_all[score5] = (T *) NULL;
    } else {
      hits5_all[score5] = (T *) CALLOC(nhits5_all[score5],sizeof(Stage3_T));
    }
  }

  for (score5 = 0; score5 <= cutoff_level_5; score5++) {
    nhits5_all[score5] = 0;
  }

  for (i = 0; i < narray5; i++) {
    for (q = hitarray5[i]; q != NULL; q = q->rest) {
      hit5 = (T) q->first;
      if (hit5->score > cutoff_level_5) {
      } else if (hit5->chrnum == 0) {
	/* Translocation: Enter under each substring chrnum */
	copy = Stage3_copy(hit5);
	copy->effective_chrnum = Substring_chrnum(hit5->substring1);
	copy->genomicstart = Substring_genomicstart(hit5->substring1);
	copy->genomicend = Substring_genomicend(hit5->substring1);
	copy->plusp = Substring_plusp(hit5->substring1);
	copies = List_push(copies,(void *) copy);
	hits5_all[hit5->score][nhits5_all[hit5->score]++] = copy;

	copy = Stage3_copy(hit5);
	copy->effective_chrnum = Substring_chrnum(hit5->substring2);
	copy->genomicstart = Substring_genomicstart(hit5->substring2);
	copy->genomicend = Substring_genomicend(hit5->substring2);
	copy->plusp = Substring_plusp(hit5->substring2);
	copies = List_push(copies,(void *) copy);
	hits5_all[hit5->score][nhits5_all[hit5->score]++] = copy;
      } else {
	hits5_all[hit5->score][nhits5_all[hit5->score]++] = hit5;
      }
    }
  }

  /* Store hits3 */
  hits3_all = (T **) CALLOC(cutoff_level_3+1,sizeof(T *));
  for (score3 = 0; score3 <= cutoff_level_3; score3++) {
    if (nhits3_all[score3] == 0) {
      hits3_all[score3] = (T *) NULL;
    } else {
      hits3_all[score3] = (T *) CALLOC(nhits3_all[score3],sizeof(Stage3_T));
    }
  }

  for (score3 = 0; score3 <= cutoff_level_3; score3++) {
    nhits3_all[score3] = 0;
  }

  for (i = 0; i < narray3; i++) {
    for (q = hitarray3[i]; q != NULL; q = q->rest) {
      hit3 = (T) q->first;
      if (hit3->score > cutoff_level_3) {
      } else if (hit3->chrnum == 0) {
	/* Translocation: Enter under each substring chrnum */
	copy = Stage3_copy(hit3);
	copy->effective_chrnum = Substring_chrnum(hit3->substring1);
	copy->genomicstart = Substring_genomicstart(hit3->substring1);
	copy->genomicend = Substring_genomicend(hit3->substring1);
	copy->plusp = Substring_plusp(hit3->substring1);
	copies = List_push(copies,(void *) copy);
	hits3_all[hit3->score][nhits3_all[hit3->score]++] = copy;

	copy = Stage3_copy(hit3);
	copy->effective_chrnum = Substring_chrnum(hit3->substring2);
	copy->genomicstart = Substring_genomicstart(hit3->substring2);
	copy->genomicend = Substring_genomicend(hit3->substring2);
	copy->plusp = Substring_plusp(hit3->substring2);
	copies = List_push(copies,(void *) copy);
	hits3_all[hit3->score][nhits3_all[hit3->score]++] = copy;
      } else {
	hits3_all[hit3->score][nhits3_all[hit3->score]++] = hit3;
      }
    }
  }


  /* Look for samechr pairs */
  sorted5p = (bool *) CALLOC(cutoff_level_5+1,sizeof(bool));
  sorted3p = (bool *) CALLOC(cutoff_level_3+1,sizeof(bool));

  pairscore = 0;
  while (abort_samechr_p == false && pairscore <= *found_score + subopt_levels) {
    for (score5 = 0; score5 <= pairscore; score5++) {
      if (score5 <= cutoff_level_5 && ((score3 = pairscore - score5) <= cutoff_level_3)) {
	/* Sort if necessary, at least by chrnum */
	if (sorted5p[score5] == false) {
	  if (nhits5_all[score5] > 0) {
	    qsort(hits5_all[score5],nhits5_all[score5],sizeof(T),genomicstart_cmp);
	  }
	  sorted5p[score5] = true;
	}
	if (sorted3p[score3] == false) {
	  if (nhits3_all[score3] > 0) {
	    qsort(hits3_all[score3],nhits3_all[score3],sizeof(T),genomicstart_cmp);
	  }
	  sorted3p[score3] = true;
	}

	/* hits5 against hits3 */
	i = j = 0;
	hits5 = hits5_all[score5];
	hits3 = hits3_all[score3];
	nhits5 = nhits5_all[score5];
	nhits3 = nhits3_all[score3];

	while (abort_samechr_p == false && i < nhits5 && j < nhits3) {
	  hit5 = hits5[i];
	  hit3 = hits3[j];
	  debug6(printf("either: i=%d %u (chr %d), j=%d %u (chr %d)",
			i,hit5->genomicstart,hit5->effective_chrnum,j,hit3->genomicstart,hit3->effective_chrnum));

	  if (hit3->effective_chrnum < hit5->effective_chrnum) {
	    j++;
	  } else if (hit5->effective_chrnum < hit3->effective_chrnum) {
	    i++;
	  } else if (hit5->concordantp == true) {
	    /* Disregard a hit involved in a concordant pair */
	    abort();		/* Should be calling this procedure only if nconcordant == 0 */
	    debug6(printf(" => hit5 concordant"));
	    i++;
	  } else {
	    for (jj = j; jj < nhits3 && hits3[jj]->effective_chrnum == hit5->effective_chrnum; jj++) {
	      hit3 = hits3[jj];
	      if (hit3->concordantp == true) {
		/* Disregard a hit involved in a concordant pair */
		abort();		/* Should be calling this procedure only if nconcordant == 0 */
		debug6(printf(" => hit3 concordant"));
	      } else /* if (hit5->paired_seenp == false || hit3->paired_seenp == false) */ {
		debug6(printf(" => samechr chr %d",hit5->effective_chrnum));
		hitpairs = List_push(hitpairs,Stage3pair_new(hit5,hit3,querylength5,querylength3,
							     expected_pairlength,/*pairtype*/SAMECHR));
		if (pairscore < *found_score) {
		  *found_score = pairscore;
		}
		if (++nsamechr >= maxpairedpaths) {
		  debug(printf("%d samechr paths exceeds %d\n",nsamechr,maxpairedpaths));
		  abort_samechr_p = true;
		  for (q = hitpairs; q != NULL; q = List_next(q)) {
		    stage3pair = (Stage3pair_T) List_head(q);
		    Stage3pair_free(&stage3pair);
		  }
		  List_free(&hitpairs);
		  hitpairs = (List_T) NULL;
		}
	      }
	    }
	    i++;
	  }
	  debug6(printf("\n"));
	}
      }
    }
    pairscore++;
  }

  FREE(sorted3p);
  FREE(sorted5p);

  for (q = copies; q != NULL; q = List_next(q)) {
    copy = (T) List_head(q);
    Stage3_free(&copy);
  }
  List_free(&copies);

  for (score3 = 0; score3 <= cutoff_level_3; score3++) {
    FREE(hits3_all[score3]);
  }
  FREE(hits3_all);
  FREE(nhits3_all);

  for (score5 = 0; score5 <= cutoff_level_5; score5++) {
    FREE(hits5_all[score5]);
  }
  FREE(hits5_all);
  FREE(nhits5_all);

  debug6(printf("Finished with Stage3_pair_up_samechr\n"));

  return hitpairs;
}

