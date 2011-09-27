static char rcsid[] = "$Id: stage3hr.c,v 1.88 2010/03/10 01:34:52 twu Exp $";
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


#if 0
#define TRANSLOCATION_TEXT "pair_translocation"
#define INVERSION_TEXT "pair_inversion"
#define SCRAMBLE_TEXT "pair_scramble"
#define TOOLONG_TEXT "pair_toolong"
#endif

#define CONCORDANT_TEXT "concordant"
#define SAMECHR_TOOLONG_TEXT "toolong"
#define SAMECHR_INVERSION_TEXT "inversion"
#define SAMECHR_SCRAMBLE_TEXT "scramble"
#define UNPAIRED_TEXT "unpaired"


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

/* Stage3_pair_up */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
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
  Genomicpos_T chroffset;

  Genomicpos_T genomicstart;
  Genomicpos_T genomicend;
  bool plusp;

  int score;			/* Includes colordiffs and penalties */
  int ntscore;			/* Includes penalties */
  int total_nmismatches;
  int geneprob;

  int nindels;			/* for indels */
  int indel_pos;		/* for indels */
  char *deletion;		/* for deletions */

  double chimera_prob;		/* for splicing */
  Genomicpos_T distance;	/* for splicing */
  bool sensep;			/* for splicing */
  double half_intron_score;	/* combination of support, nmismatches, and splice site probability */

  Substring_T substring1;	/* indel or donor */
  Substring_T substring2;	/* indel or acceptor */

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
  Genomicpos_T absdifflength;
  int dir;			/* -1, 0, or +1 */
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

Genomicpos_T
Stage3_distance (T this) {
  return this->distance;
}

bool
Stage3_sensep (T this) {
  return this->sensep;
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
  new->distance = old->distance;
  new->half_intron_score = old->half_intron_score;
  new->sensep = old->sensep;

  new->substring1 = Substring_copy(old->substring1);
  new->substring2 = Substring_copy(old->substring2);

  new->paired_seenp = old->paired_seenp;
  new->concordantp = old->concordantp;

  return new;
}


T
Stage3_new_exact (int *found_score, Genomicpos_T left, int genomiclength, bool plusp, Chrnum_T chrnum, Genomicpos_T chroffset) {
  T new = (T) MALLOC(sizeof(*new));

  new->deletion = (char *) NULL;
  if (plusp == true) {
    new->genomicstart = left;
    new->genomicend = left + genomiclength;
  } else {
    new->genomicstart = left + genomiclength;
    new->genomicend = left;
  }
  new->substring1 = Substring_new(/*nmismatches*/0,/*ncolordiffs*/0,chrnum,
				  chroffset,new->genomicstart,new->genomicend,
				  /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				  /*alignstart*/new->genomicstart,/*alignend*/new->genomicend,
				  genomiclength,/*extraleft*/0,/*extraright*/0,
				  /*genomicseg*/NULL,/*query*/NULL,plusp,
				  /*trim_ends_p*/false,/*dibasep*/false,/*cmetp*/false);

  new->substring2 = (Substring_T) NULL;

  new->hittype = EXACT;
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->plusp = plusp;

  new->nindels = 0;
  new->total_nmismatches = 0;
  new->ntscore = 0;
  new->score = 0;
  new->geneprob = -1;
  *found_score = 0;

  new->chimera_prob = 2.0;
  new->distance = 0U;
  new->half_intron_score = 0.0;

  new->paired_seenp = false;
  new->concordantp = false;

  return new;
}


T
Stage3_new_substitution (int *found_score, int nmismatches, int ncolordiffs, Genomicpos_T left,
			 int genomiclength, bool plusp, char *genomicseg, char *query,
			 Chrnum_T chrnum, Genomicpos_T chroffset, 
			 bool trim_ends_p, bool dibasep, bool cmetp) {
  T new = (T) MALLOC(sizeof(*new));

  debug(printf("Entered new_substitution with nmismatches %d\n",nmismatches));

  new->deletion = (char *) NULL;

  if (plusp == true) {
    new->genomicstart = left;
    new->genomicend = left + genomiclength;
  } else {
    new->genomicstart = left + genomiclength;
    new->genomicend = left;
  }

  new->substring1 = Substring_new(nmismatches,ncolordiffs,chrnum,
				  chroffset,new->genomicstart,new->genomicend,
				  /*querystart*/0,/*queryend*/genomiclength,/*querylength*/genomiclength,
				  /*alignstart*/new->genomicstart,/*alignend*/new->genomicend,
				  genomiclength,/*extraleft*/0,/*extraright*/0,
				  genomicseg,query,plusp,trim_ends_p,dibasep,cmetp);

  new->substring2 = (Substring_T) NULL;

  new->hittype = SUB;
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->plusp = plusp;

  new->nindels = 0;
  new->total_nmismatches = nmismatches;
  new->ntscore = nmismatches;
  new->score = nmismatches + ncolordiffs;
  new->geneprob = -1;

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->chimera_prob = 2.0;
  new->distance = 0U;
  new->half_intron_score = 0.0;

  new->paired_seenp = false;
  new->concordantp = false;

  return new;
}



T
Stage3_new_insertion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
		      int ncolordiffs1, int ncolordiffs2,
		      Genomicpos_T left, int genomiclength, int querylength, bool plusp,
		      char *genomicseg, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
		      int indel_penalty, bool trim_ends_p, bool dibasep, bool cmetp) {
  T new = (T) MALLOC(sizeof(*new));
  int querystart1, queryend1, querystart2, queryend2;
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
    new->genomicstart = left;
    new->genomicend = left + genomiclength;

    alignstart1 = new->genomicstart;
    alignend1 = alignstart2 = new->genomicstart + indel_pos;
    alignend2 = new->genomicend/* - nindels*/;

  } else {
    new->genomicstart = left + genomiclength;
    new->genomicend = left;

    alignstart1 = new->genomicstart;
    alignend1 = alignstart2 = new->genomicstart - indel_pos;
    alignend2 = new->genomicend/* + nindels*/;

  }

  new->deletion = (char *) NULL;
  new->substring1 = Substring_new(nmismatches1,ncolordiffs1,chrnum,
				  chroffset,new->genomicstart,new->genomicend,
				  querystart1,queryend1,querylength,
				  alignstart1,alignend1,genomiclength,
				  /*extraleft*/0,/*extraright*/0,
				  genomicseg,query,plusp,trim_ends_p,dibasep,cmetp);

  new->substring2 = Substring_new(nmismatches2,ncolordiffs2,chrnum,
				  chroffset,new->genomicstart,new->genomicend,
				  querystart2,queryend2,querylength,
				  alignstart2,alignend2,genomiclength,
				  /*extraleft*/0,/*extraright*/0,
				  genomicseg,query,plusp,trim_ends_p,dibasep,cmetp);

  new->hittype = INS;
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->plusp = plusp;

  new->nindels = nindels;
  new->total_nmismatches = nmismatches1 + nmismatches2;
  new->ntscore = indel_penalty + nmismatches1 + nmismatches2;
  new->score = new->ntscore + ncolordiffs1 + ncolordiffs2;
  new->geneprob = -1;

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->indel_pos = indel_pos;

  new->chimera_prob = 2.0;
  new->distance = 0U;
  new->half_intron_score = 0.0;

  new->paired_seenp = false;
  new->concordantp = false;

  return new;
}


T
Stage3_new_deletion (int *found_score, int nindels, int indel_pos, int nmismatches1, int nmismatches2,
		     int ncolordiffs1, int ncolordiffs2,
		     Genomicpos_T left, int genomiclength, int querylength, bool plusp,
		     char *genomicseg, char *query, Chrnum_T chrnum, Genomicpos_T chroffset,
		     int indel_penalty, bool trim_ends_p, bool dibasep, bool cmetp) {
  T new = (T) MALLOC(sizeof(*new));
  int querystart1, queryend1, querystart2, queryend2;
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
    new->genomicstart = left;
    new->genomicend = left + genomiclength;
    alignstart1 = new->genomicstart;
    alignend1 = new->genomicstart + indel_pos;
    alignstart2 = alignend1 + nindels;
    alignend2 = new->genomicend/* + nindels*/;

  } else {
    new->genomicstart = left + genomiclength;
    new->genomicend = left;
    alignstart1 = new->genomicstart;
    alignend1 = new->genomicstart - indel_pos;
    alignstart2 = alignend1 - nindels;
    alignend2 = new->genomicend/* - nindels*/;
  }

  new->substring1 = Substring_new(nmismatches1,ncolordiffs1,chrnum,chroffset,
				  new->genomicstart,new->genomicend,
				  querystart1,queryend1,querylength,
				  alignstart1,alignend1,genomiclength,
				  /*extraleft*/0,/*extraright*/0,	
				  genomicseg,query,plusp,trim_ends_p,dibasep,cmetp);

  new->deletion = (char *) CALLOC(nindels+1,sizeof(char));
  if (plusp == true) {
    strncpy(new->deletion,&(genomicseg[indel_pos]),nindels);
  } else {
    make_complement_buffered(new->deletion,&(genomicseg[querylength-indel_pos]),nindels);
  }

  new->substring2 = Substring_new(nmismatches2,ncolordiffs2,chrnum,chroffset,
				  new->genomicstart,new->genomicend,
				  querystart2,queryend2,querylength,
				  alignstart2,alignend2,genomiclength,
				  /*extraleft*/0,/*extraright*/0,
				  genomicseg,query,plusp,trim_ends_p,dibasep,cmetp);

  new->hittype = DEL;
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->plusp = plusp;

  new->nindels = nindels;
  new->total_nmismatches = nmismatches1 + nmismatches2;
  new->ntscore = indel_penalty + nmismatches1 + nmismatches2;
  new->score = new->ntscore + ncolordiffs1 + ncolordiffs2;
  new->geneprob = -1;

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->indel_pos = indel_pos;

  new->chimera_prob = 2.0;
  new->distance = 0U;
  new->half_intron_score = 0.0;

  new->paired_seenp = false;
  new->concordantp = false;

  return new;
}


int
Stage3_output_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->score < y->score) {
    return -1;
  } else if (x->score > y->score) {
    return +1;
  } else if (x->hittype < y->hittype) {
    return -1;
  } else if (x->hittype > y->hittype) {
    return +1;
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
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

  n = List_length(hitlist);
  if (n == 0) {
    return NULL;
  }

  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    /* For dibasep were previously using hit->ntscore, but gives false positives */
    if (hit->score <= cutoff_level) {
      if (hit->score < minscore) {
	minscore = hit->score;
      }
    }
  }

  debug(printf("Stage3_optimal_score: minscore = %d + subopt:%d, nhits = %d\n",
	       minscore,suboptimal_mismatches,n));
  minscore += suboptimal_mismatches;
  for (p = hitlist; p != NULL; p = p->rest) {
    hit = (T) p->first;
    /* For dibasep were previously using hit->ntscore, but gives false positives */
    if (hit->score > cutoff_level) {
      debug(printf("Eliminating a hit with ntscore %d > cutoff_level %d\n",
		   hit->ntscore,cutoff_level));
      Stage3_free(&hit);
    } else if (hit->score <= minscore) {
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
Stage3_remove_duplicates (List_T hitlist) {
#ifdef DEBUG
  List_T p;
#endif
  List_T unique = NULL;
  T x, y, hit, *hits, *prev;
  int n, i, j;
  bool *eliminate;

  n = List_length(hitlist);
  if (n == 0) {
    return NULL;
  }

  eliminate = (bool *) CALLOC(n,sizeof(bool));

  hits = (T *) List_to_array(hitlist,NULL);
  List_free(&hitlist);


  /* By genomicstart */
  debug(printf("Checking %d hits for duplicates by genomicstart\n",n));
  qsort(hits,n,sizeof(T),genomicstart_cmp);
  for (i = 0; i < n; i++) {
    x = hits[i];
    j = i+1;
    while (j < n && hits[j]->genomicstart == x->genomicstart) {
      y = hits[j];
      debug(printf("  #%d has common start with #%d\n",i,j));
      if ((x->half_intron_score == 0.0 || y->half_intron_score == 0.0) && x->score < y->score) {
	debug(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
		     i,j,x->score,y->score,j));
	eliminate[j] = true;
      } else if ((x->half_intron_score == 0.0 || y->half_intron_score == 0.0) && x->score > y->score) {
	debug(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
		     i,j,x->score,y->score,i));
	eliminate[i] = true;
      } else if (x->hittype < y->hittype) {
	debug(printf("  #%d overlaps #%d and score %d == %d, so marking %d for elimination\n",
		     i,j,x->score,y->score,j));
	eliminate[j] = true;
      } else if (x->hittype > y->hittype) {
	debug(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
		     i,j,x->score,y->score,i));
	eliminate[i] = true;
      } else if (x->nindels < y->nindels) {
	debug(printf("  #%d overlaps #%d and fewer indels, so marking %d for elimination\n",
		     i,j,j));
	eliminate[j] = true;
      } else if (x->nindels > y->nindels) {
	debug(printf("  #%d overlaps #%d and fewer indels, so marking %d for elimination\n",
		     i,j,i));
	eliminate[i] = true;
      } else if (x->chrnum == 0 && y->chrnum != 0) {
	debug(printf("  #%d overlaps #%d and first one is a distant splice, so marking first one for elimination\n",
		     i,j));
	eliminate[i] = true;
      } else if (x->chrnum > 0 && y->chrnum == 0) {
	debug(printf("  #%d overlaps #%d and second one is a distant splice, so marking second one for elimination\n",
		     i,j));
	eliminate[j] = true;

      } else if (x->distance > y->distance) {
	debug(printf("  #%d overlaps #%d and first one has longer distance, so marking first one for elimination\n",
		     i,j));
	eliminate[i] = true;
      } else if (x->distance < y->distance) {
	debug(printf("  #%d overlaps #%d and second one has longer distance, so marking second one for elimination\n",
		     i,j));
	eliminate[j] = true;

      } else if (x->genomicend == y->genomicend) {
	if (x->half_intron_score != 0.0 && y->half_intron_score != 0.0) {
	  if (x->half_intron_score > y->half_intron_score) {
	    debug(printf("  #%d overlaps #%d and equal and first has better support score %f > %f, so marking second one for elimination\n",
			 i,j,x->half_intron_score,y->half_intron_score));
	    eliminate[j] = true;
	  } else if (y->half_intron_score > x->half_intron_score) {
	    debug(printf("  #%d overlaps #%d and equal and second has better support score %f > %f, so marking first one for elimination\n",
			 i,j,y->half_intron_score,x->half_intron_score));
	    eliminate[i] = true;
	  } else {
	    debug(printf("  #%d overlaps #%d and equal so eliminate second one\n",
			 i,j));
	    eliminate[j] = true;
	  }
	} else {
	  debug(printf("  #%d overlaps #%d and equal so eliminate second one\n",
		       i,j));
	  eliminate[j] = true;
	}
      }
      j++;
    }
  }
    
  j = 0;
  prev = hits;
  hits = (T *) CALLOC(n,sizeof(T));
  for (i = n-1; i >= 0; i--) {
    hit = prev[i];
    if (eliminate[i] == false) {
      debug(printf("  Keeping %u..%u (nindels %d, distance %u, chrnum %d) (plusp = %d) on basis of genomicstart\n",
		   hit->genomicstart,hit->genomicend,hit->nindels,hit->distance,hit->chrnum,hit->plusp));
      hits[j++] = hit;
    } else {
      debug(printf("  Eliminating %u..%u (nindels %d, distance %u, chrnum %d) (plusp = %d) on basis of genomicstart\n",
		   hit->genomicstart,hit->genomicend,hit->nindels,hit->distance,hit->chrnum,hit->plusp));
      Stage3_free(&hit);
    }
  }
  FREE(prev);
  FREE(eliminate);

  /* By genomicend */
  n = j;
  eliminate = (bool *) CALLOC(n,sizeof(bool));
  debug(printf("Checking %d hits for duplicates by genomicend\n",n));
  qsort(hits,n,sizeof(T),genomicend_cmp);
  for (i = 0; i < n; i++) {
    x = hits[i];
    j = i+1;
    while (j < n && hits[j]->genomicend == x->genomicend) {
      y = hits[j];
      debug(printf("  #%d has common end with #%d\n",i,j));
      if ((x->half_intron_score == 0.0 || y->half_intron_score == 0.0) && x->score < y->score) {
	debug(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
		     i,j,x->score,y->score,j));
	eliminate[j] = true;
      } else if ((x->half_intron_score == 0.0 || y->half_intron_score == 0.0) && x->score > y->score) {
	debug(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
		     i,j,x->score,y->score,i));
	eliminate[i] = true;
      } else if (x->hittype < y->hittype) {
	debug(printf("  #%d overlaps #%d and score %d == %d, so marking %d for elimination\n",
		     i,j,x->score,y->score,j));
	eliminate[j] = true;
      } else if (x->hittype > y->hittype) {
	debug(printf("  #%d overlaps #%d and score %d < %d, so marking %d for elimination\n",
		     i,j,x->score,y->score,i));
	eliminate[i] = true;
      } else if (x->nindels < y->nindels) {
	debug(printf("  #%d overlaps #%d and fewer indels, so marking %d for elimination\n",
		     i,j,j));
	eliminate[j] = true;
      } else if (x->nindels > y->nindels) {
	debug(printf("  #%d overlaps #%d and fewer indels, so marking %d for elimination\n",
		     i,j,i));
	eliminate[i] = true;
      } else if (x->chrnum == 0 && y->chrnum != 0) {
	debug(printf("  #%d overlaps #%d and first one is a distant splice, so marking first one for elimination\n",
		     i,j));
	eliminate[i] = true;
      } else if (x->chrnum != 0 && y->chrnum == 0) {
	debug(printf("  #%d overlaps #%d and second one is a distant splice, so marking second one for elimination\n",
		     i,j));
	eliminate[j] = true;

      } else if (x->distance > y->distance) {
	debug(printf("  #%d overlaps #%d and first one has longer distance, so marking first one for elimination\n",
		     i,j));
	eliminate[i] = true;
      } else if (x->distance < y->distance) {
	debug(printf("  #%d overlaps #%d and second one has longer distance, so marking second one for elimination\n",
		     i,j));
	eliminate[j] = true;

      } else if (x->genomicstart == y->genomicstart) {
	if (x->half_intron_score != 0.0 && y->half_intron_score != 0.0) {
	  if (x->half_intron_score > y->half_intron_score) {
	    debug(printf("  #%d overlaps #%d and equal and first has better support score %f > %f, so marking second one for elimination\n",
			 i,j,x->half_intron_score,y->half_intron_score));
	    eliminate[j] = true;
	  } else if (y->half_intron_score > x->half_intron_score) {
	    debug(printf("  #%d overlaps #%d and equal and second has better support score %f > %f, so marking first one for elimination\n",
			 i,j,y->half_intron_score,x->half_intron_score));
	    eliminate[i] = true;
	  } else {
	    debug(printf("  #%d overlaps #%d and equal so eliminate second one\n",
			 i,j));
	    eliminate[j] = true;
	  }
	} else {
	  debug(printf("  #%d overlaps #%d and equal so eliminate second one\n",
		       i,j));
	  eliminate[j] = true;
	}
      }
      j++;
    }
  }

  for (i = n-1; i >= 0; i--) {
    hit = hits[i];
    if (eliminate[i] == false) {
      debug(printf("  Keeping %u..%u (nindels %d, distance %u, chrnum %d) (plusp = %d)\n",
		   hit->genomicstart,hit->genomicend,hit->nindels,hit->distance,hit->chrnum,hit->plusp));
      unique = List_push(unique,hit);
    } else {
      debug(printf("  Eliminating %u..%u (nindels %d, distance %u, chrnum %d) (plusp = %d)\n",
		   hit->genomicstart,hit->genomicend,hit->nindels,hit->distance,hit->chrnum,hit->plusp));
      Stage3_free(&hit);
    }
  }
  FREE(hits);
  FREE(eliminate);

  debug(
	for (p = unique, i = 0; p != NULL; p = p->rest, i++) {
	  hit = (T) p->first;
	  printf("  Final %d: %u..%u (plusp = %d)\n",i,hit->genomicstart,hit->genomicend,hit->plusp);
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
samechr_text (T hit5, T hit3) {
  if (hit5->chrnum != hit3->chrnum) {
    return (char *) NULL;
  } else if (hit5->plusp != hit3->plusp) {
    return SAMECHR_INVERSION_TEXT;
  } else if (hit5->plusp == true) {
    if (hit3->genomicstart < hit5->genomicstart) {
      return SAMECHR_SCRAMBLE_TEXT;
    } else {
      return SAMECHR_TOOLONG_TEXT;
    }
  } else {
    if (hit5->genomicstart < hit3->genomicstart) {
      return SAMECHR_SCRAMBLE_TEXT;
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
print_insertion (T this, int score, int geneprob, IIT_T chromosome_iit,
		 bool invertp, IIT_T snps_iit, int *snps_divint_crosstable,
		 T hit5, T hit3, Pairtype_T pairtype, int pairedlength, int pairscore) {
  char *chr;
  bool allocp;

  chr = IIT_label(chromosome_iit,this->chrnum,&allocp);

  printf(" ");
  Substring_print_insertion_1(this->substring1,this->substring2,this->nindels,
			      chr,invertp,snps_iit,snps_divint_crosstable);
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
			      chr,invertp,snps_iit,snps_divint_crosstable);
  printf("\n");

  if (allocp == true) {
    FREE(chr);
  }

  return;
}

static void
print_deletion (T this, int score, int geneprob, IIT_T chromosome_iit,
		bool invertp, IIT_T snps_iit, int *snps_divint_crosstable,
		T hit5, T hit3, Pairtype_T pairtype, int pairedlength, int pairscore) {
  char *chr;
  bool allocp;

  chr = IIT_label(chromosome_iit,this->chrnum,&allocp);

  printf(" ");
  Substring_print_deletion_1(this->substring1,this->substring2,this->nindels,this->deletion,
			     chr,invertp,snps_iit,snps_divint_crosstable);
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
			     chr,invertp,snps_iit,snps_divint_crosstable);
  printf("\n");

  if (allocp == true) {
    FREE(chr);
  }
}


static void
print_splice (T chimera, int score, int geneprob, Genome_T genome, IIT_T chromosome_iit,
	      bool invertp, IIT_T snps_iit, int *snps_divint_crosstable, IIT_T splicesites_iit,
	      int donor_typeint, int acceptor_typeint, int *splicesites_divint_crosstable,
	      T hit5, T hit3, Pairtype_T pairtype, int pairedlength, int pairscore) {
  Substring_T donor, acceptor;
  
  donor = chimera->substring1;
  acceptor = chimera->substring2;

  Substring_assign_donor_prob(donor,genome,chromosome_iit);
  Substring_assign_acceptor_prob(acceptor,genome,chromosome_iit);

  if (acceptor == NULL) {
    /* Single sequence */
    printf(" ");
    Substring_print_donor(donor,chimera->sensep,invertp,
			  chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
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
    /* Single sequence */
    printf(" ");
    Substring_print_acceptor(acceptor,chimera->sensep,invertp,
			     chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
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

  } else if (chimera->sensep == true && invertp == false) {
    printf(" ");
    Substring_print_donor(donor,chimera->sensep,invertp,
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
    Substring_print_acceptor(acceptor,chimera->sensep,invertp,
			     chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			     acceptor_typeint,splicesites_divint_crosstable,/*endp*/false,
			     donor,chimera->distance);
    printf("\n");

  } else if (chimera->sensep == true && invertp == true) {
    printf(" ");
    Substring_print_acceptor(acceptor,chimera->sensep,invertp,
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
    Substring_print_donor(donor,chimera->sensep,invertp,
			  chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			  donor_typeint,splicesites_divint_crosstable,/*endp*/false,
			  acceptor,chimera->distance);
    printf("\n");

  } else if (chimera->sensep == false && invertp == false) {
    printf(" ");
    Substring_print_acceptor(acceptor,chimera->sensep,invertp,
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
    Substring_print_donor(donor,chimera->sensep,invertp,
			  chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			  donor_typeint,splicesites_divint_crosstable,/*endp*/false,
			  acceptor,chimera->distance);
    printf("\n");

  } else if (chimera->sensep == false && invertp == true) {
    printf(" ");
    Substring_print_donor(donor,chimera->sensep,invertp,
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
    Substring_print_acceptor(acceptor,chimera->sensep,invertp,
			     chromosome_iit,snps_iit,snps_divint_crosstable,splicesites_iit,
			     acceptor_typeint,splicesites_divint_crosstable,/*endp*/false,
			     donor,chimera->distance);
    printf("\n");
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
    print_insertion(this,score,geneprob,chromosome_iit,invertp,
		    snps_iit,snps_divint_crosstable,hit5,hit3,pairtype,pairedlength,pairscore);
  } else if (this->hittype == DEL) {
    print_deletion(this,score,geneprob,chromosome_iit,invertp,
		   snps_iit,snps_divint_crosstable,hit5,hit3,pairtype,pairedlength,pairscore);
  } else if (this->hittype == SPLICE) {
    print_splice(this,score,geneprob,genome,chromosome_iit,invertp,
		 snps_iit,snps_divint_crosstable,splicesites_iit,
		 donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		 hit5,hit3,pairtype,pairedlength,pairscore);
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
						      geneprob_iit,chromosome_iit,Sequence_fulllength(queryseq),this->sensep);
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
						      geneprob_iit,chromosome_iit,Sequence_fulllength(queryseq5),this->sensep);
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
						      geneprob_iit,chromosome_iit,Sequence_fulllength(queryseq3),this->sensep);
    }
  }

  return;
}



static void
print_paired_hits (Result_T result, char initchar, bool fivep, 
		   Genome_T genome, IIT_T chromosome_iit,
		   Sequence_T queryseq, Sequence_T headerseq, 
		   IIT_T snps_iit, int *snps_divint_crosstable,
		   IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		   int *splicesites_divint_crosstable, int maxpaths, bool quiet_if_excessive_p,
		   bool sequence_inverted_p, bool invertp) {
  Stage3pair_T *stage3pairarray, stage3pair;
  Resulttype_T resulttype;
  T *stage3array, *stage3array1, *stage3array2, this, hit5, hit3;
  int npaths, pathnum;
  bool samechrp;
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
  if (resulttype == PAIREDEND_SAMECHR) {
    /* Should be discovered only within PAIREDEND_AS_SINGLES_UNIQUE */
    abort();

  } else if (resulttype == PAIREDEND_CONCORDANT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);

    printf("\t%d %s\t",npaths,CONCORDANT_TEXT);
    Sequence_print_header(headerseq,/*checksump*/false);
    /* printf("\n"); -- included in header */

    if (npaths == 0 || (quiet_if_excessive_p && npaths > maxpaths)) {
    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];
	if (fivep == true) {
	  Stage3_print(stage3pair->hit5,stage3pair->hit5->score,
		       stage3pair->hit5->geneprob+stage3pair->hit3->geneprob,
		       genome,chromosome_iit,queryseq,
		       snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       invertp,
		       stage3pair->hit5,stage3pair->hit3,stage3pair->pairtype,stage3pair->pairedlength,
		       stage3pair->score);
	} else {
	  Stage3_print(stage3pair->hit3,stage3pair->hit3->score,
		       stage3pair->hit5->geneprob+stage3pair->hit3->geneprob,
		       genome,chromosome_iit,queryseq,
		       snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		       invertp,
		       stage3pair->hit5,stage3pair->hit3,stage3pair->pairtype,stage3pair->pairedlength,
		       stage3pair->score);
	}
      }
    }

  } else if (resulttype == PAIREDEND_AS_SINGLES || resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
    samechrp = false;
    if (resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
      stage3array1 = (T *) Result_array(&npaths,result);
      stage3array2 = (T *) Result_array2(&npaths,result);
      hit5 = stage3array1[0];
      hit3 = stage3array2[0];

      if (hit5->chrnum == hit3->chrnum) {
	samechrp = true;
      }
    }

    if (samechrp == false) {
      if (fivep == true) {
	stage3array = (T *) Result_array(&npaths,result);
      } else {
	stage3array = (T *) Result_array2(&npaths,result);
      }
      printf("\t%d %s\t",npaths,UNPAIRED_TEXT);
      Sequence_print_header(headerseq,/*checksump*/false);
      /* printf("\n"); -- included in header */
      
      if (npaths == 0 || (quiet_if_excessive_p && npaths > maxpaths)) {
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

    } else {
      printf("\t%d %s\t",npaths,samechr_text(hit5,hit3));
      Sequence_print_header(headerseq,/*checksump*/false);
      /* printf("\n"); -- included in header */

      if (hit5->plusp == true) {
	high5 = hit5->genomicend;
	low5 = hit5->genomicstart;
      } else {
	high5 = hit5->genomicstart;
	low5 = hit5->genomicend;
      }

      if (hit3->plusp == true) {
	high3 = hit3->genomicend;
	low3 = hit3->genomicstart;
      } else {
	high3 = hit3->genomicstart;
	low3 = hit3->genomicend;
      }

      if (low3 > high5) {
	pairedlength = low3 - high5;
      } else if (low5 > high3) {
	pairedlength = low5 - high3;
      } else {
	pairedlength = 0U;
      }

      if (fivep == true) {
	Stage3_print(hit5,hit5->score,hit5->geneprob+hit3->geneprob,
		     genome,chromosome_iit,queryseq,
		     snps_iit,snps_divint_crosstable,
		     splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		     invertp,
		     hit5,hit3,/*pairtype*/SAMECHR,pairedlength,
		     hit5->score+hit3->score);
      } else {
	Stage3_print(hit3,hit3->score,hit5->geneprob+hit3->geneprob,
		     genome,chromosome_iit,queryseq,
		     snps_iit,snps_divint_crosstable,
		     splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		     invertp,
		     hit5,hit3,/*pairtype*/SAMECHR,pairedlength,
		     hit5->score+hit3->score);
      }

    }

  } else {
    /* abort() */
  }

  printf("\n");

  return;
}


void
Stage3_print_paired (Result_T result, Genome_T genome, IIT_T chromosome_iit, Sequence_T queryseq1, Sequence_T queryseq2,
		     IIT_T snps_iit, int *snps_divint_crosstable,
		     IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		     int *splicesites_divint_crosstable, int maxpaths, bool quiet_if_excessive_p,
		     bool circularp, bool invertp) {
  
  if (circularp == false) {
    /* First sequence */
    print_paired_hits(result,'>',/*fivep*/true,
		      genome,chromosome_iit,queryseq1,/*headerseq*/queryseq1,
		      snps_iit,snps_divint_crosstable,
		      splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		      maxpaths,quiet_if_excessive_p,/*sequence_inverted_p*/false,/*invertp*/false);

    /* Second sequence.  query2 has already been revcomp'd */
    print_paired_hits(result,'<',/*fivep*/false,
		      genome,chromosome_iit,queryseq2,/*headerseq*/queryseq1,
		      snps_iit,snps_divint_crosstable,
		      splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		      maxpaths,quiet_if_excessive_p,/*sequence_inverted_p*/true,/*invertp*/invertp);

  } else {
    /* First sequence.  query2 has already been revcomp'd */
    print_paired_hits(result,'>',/*fivep*/false,
		      genome,chromosome_iit,queryseq2,/*headerseq*/queryseq2,
		      snps_iit,snps_divint_crosstable,
		      splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		      maxpaths,quiet_if_excessive_p,/*sequence_inverted_p*/true,/*invertp*/(invertp == true) ? false : true);

    /* Second sequence */
    print_paired_hits(result,'<',/*fivep*/true,
		      genome,chromosome_iit,queryseq1,/*headerseq*/queryseq2,
		      snps_iit,snps_divint_crosstable,
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

  } else {
    if (hit5->plusp == true) {
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
  }

  if (new->pairedlength < expected_pairlength) {
    new->absdifflength = expected_pairlength - new->pairedlength;
  } else {
    new->absdifflength = new->pairedlength - expected_pairlength;
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


T
Stage3_new_splice (int *found_score, Substring_T donor, Substring_T acceptor, Genomicpos_T distance,
		   bool shortdistancep, int splicing_penalty, int support, int querylength,
		   bool copyp, bool sensep) {
  T new;
  
  new = (T) MALLOC(sizeof(*new));
  new->deletion = (char *) NULL;
  new->hittype = SPLICE;
  new->nindels = 0;

  if (donor == NULL) {
    new->substring1 = (Substring_T) NULL;
    new->substring2 = copyp ? Substring_copy(acceptor) : acceptor;
    new->chrnum = Substring_chrnum(acceptor);
    new->chroffset = Substring_chroffset(acceptor);
    new->genomicstart = Substring_genomicstart(acceptor);
    new->genomicend = Substring_genomicend(acceptor);
    new->plusp = Substring_plusp(acceptor);
    new->total_nmismatches = Substring_nmismatches(acceptor);
    new->ntscore = splicing_penalty + Substring_nmismatches(acceptor);
    new->score = splicing_penalty + Substring_nmismatches(acceptor) + Substring_ncolordiffs(acceptor) + 
      (querylength - support + 8) / 8;
    new->chimera_prob = Substring_chimera_prob(acceptor);
    new->half_intron_score = support*0.25 - Substring_nmismatches(acceptor) + new->chimera_prob;
    debug(printf("%d*0.25 - %d + %f => half_intron_score %f\n",
		 support,Substring_nmismatches(acceptor),new->chimera_prob,new->half_intron_score));

  } else if (acceptor == NULL) {
    new->substring1 = copyp ? Substring_copy(donor) : donor;
    new->substring2 = (Substring_T) NULL;
    new->chrnum = Substring_chrnum(donor);
    new->chroffset = Substring_chroffset(donor);
    new->genomicstart = Substring_genomicstart(donor);
    new->genomicend = Substring_genomicend(donor);
    new->plusp = Substring_plusp(donor);
    new->total_nmismatches = Substring_nmismatches(donor);
    new->ntscore = splicing_penalty + Substring_nmismatches(donor);
    new->score = splicing_penalty + Substring_nmismatches(donor) + Substring_ncolordiffs(donor) +
      (querylength - support + 8) / 8;
    new->chimera_prob = Substring_chimera_prob(donor);
    new->half_intron_score = support*0.25 + new->chimera_prob;
    debug(printf("%d*0.25 - %d + %f => half_intron_score %f\n",
		 support,Substring_nmismatches(donor),new->chimera_prob,new->half_intron_score));

  } else {
    new->substring1 = copyp ? Substring_copy(donor) : donor;
    new->substring2 = copyp ? Substring_copy(acceptor) : acceptor;
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

    if (new->chrnum > 0) {
      if (sensep == true) {
#if 0
	printf("Case 1. sense: donor:%u..%u, acceptor:%u..%u\n",
	       Substring_genomicstart(donor),Substring_genomicend(donor),
	       Substring_genomicstart(acceptor),Substring_genomicend(acceptor));
#endif
	new->genomicstart = Substring_genomicstart(donor);
	new->genomicend = Substring_genomicend(acceptor);
      } else {
#if 0
	printf("Case 2. antisense: donor:%u..%u, acceptor:%u..%u\n",
	       Substring_genomicstart(donor),Substring_genomicend(donor),
	       Substring_genomicstart(acceptor),Substring_genomicend(acceptor));
#endif
	new->genomicstart = Substring_genomicstart(acceptor);
	new->genomicend = Substring_genomicend(donor);
      }
#if 0
      printf("Picking %u..%u\n",new->genomicstart,new->genomicend);
#endif
    }

    new->total_nmismatches = Substring_nmismatches(donor) + Substring_nmismatches(acceptor);
    new->ntscore = splicing_penalty + Substring_nmismatches(donor) + Substring_nmismatches(acceptor);
    new->score = new->ntscore + Substring_ncolordiffs(donor) + Substring_ncolordiffs(acceptor);
    new->chimera_prob = Substring_chimera_prob(donor) + Substring_chimera_prob(acceptor);
    new->half_intron_score = 0.0;
  }
  new->geneprob = -1;

  if (new->score < *found_score) {
    *found_score = new->score;
  }

  new->distance = distance;
  new->sensep = sensep;

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
  } else if (x->score < y->score) {
    return -1;
  } else if (y->score < x->score) {
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

#if 0
static int
hitpair_score_cmp (const void *a, const void *b) {
  Stage3pair_T x = * (Stage3pair_T *) a;
  Stage3pair_T y = * (Stage3pair_T *) b;

  if (x->score < y->score) {
    return -1;
  } else if (x->score > y->score) {
    return +1;
  } else if (x->absdifflength < y->absdifflength) {
    return -1;
  } else if (x->absdifflength > y->absdifflength) {
    return +1;
  } else {
    return 0;
  }
}
#endif


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
}


List_T
Stage3pair_remove_duplicates (List_T hitpairlist) {
#ifdef DEBUG
  List_T p;
#endif
  List_T unique = NULL;
  Stage3pair_T hitpair, *hitpairs;
  int n, i, j;
  bool *eliminate;

  n = List_length(hitpairlist);
  if (n == 0) {
    return NULL;
  }

  debug(
	for (p = hitpairlist, i = 0; p != NULL; p = p->rest, i++) {
	  hitpair = (Stage3pair_T) p->first;
	  printf("  Initial %d: %u-%u (dir = %d)\n",
		 i,hitpair->low,hitpair->high,hitpair->dir);
	}
	);

  eliminate = (bool *) CALLOC(n,sizeof(bool));

  hitpairs = (Stage3pair_T *) List_to_array(hitpairlist,NULL);
  List_free(&hitpairlist);

  /* Check for duplicates */
  debug(printf("Checking %d hits for duplicates\n",n));
  qsort(hitpairs,n,sizeof(Stage3pair_T),hitpair_position_cmp);
  i = 0;
  while (i < n) {
    if (eliminate[i] == true) {
      i++;
    } else {
      j = i+1;
      while (j < n && hitpair_equal(hitpairs[j],hitpairs[i]) == true) {
	debug(printf("  %d (score %d) == %d (score %d), so marking %d for elimination\n",
		     i,hitpairs[i]->score,j,hitpairs[j]->score,j));
	eliminate[j] = true;
	j++;
      }
      i = j;
    }
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

  n = List_length(hitpairlist);
  if (n == 0) {
    return NULL;
  }

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->score <= cutoff_level) {
      if (hitpair->score < minscore) {
	minscore = hitpair->score;
      }
    }
  }

  debug(printf("Stage3pair_optimal_score: minscore = %d + subopt:%d, nhits = %d\n",
	       minscore,suboptimal_mismatches,n));
  minscore += suboptimal_mismatches;

  for (p = hitpairlist; p != NULL; p = p->rest) {
    hitpair = (Stage3pair_T) p->first;
    if (hitpair->score > cutoff_level) {
      debug(printf("Eliminating a hit pair with score %d > cutoff_level %d\n",
		   hitpair->score,cutoff_level));
      Stage3pair_free(&hitpair);
    } else if (hitpair->score <= minscore) {
      optimal = List_push(optimal,hitpair);
    } else {
      debug(printf("Eliminating a hit pair with score %d\n",hitpair->score));
      Stage3pair_free(&hitpair);
    }
  }
  
  List_free(&hitpairlist);

  return optimal;
}



/* Finds concordant pairs and samechr pairs if nconcordant is 0 */
List_T
Stage3_pair_up (bool *abort_pairing_p, int *found_score, int *nconcordant, int *nsamechr, List_T hitpairs,
		List_T *hitarray5, int narray5, List_T *hitarray3, int narray3,
		int cutoff_level_5, int cutoff_level_3, int subopt_levels,
		Genomicpos_T pairmax, Genomicpos_T expected_pairlength,
		int querylength5, int querylength3, int maxpairedpaths) {
  List_T q;
  T **hits5_plus, **hits5_minus, **hits3_plus, **hits3_minus, *hits5, *hits3, hit5, hit3;
  int *nhits5_plus = 0, *nhits5_minus = 0, *nhits3_plus = 0, *nhits3_minus = 0, nhits5, nhits3;
  int pairscore, score5, score3, i, j, jj;
  bool *sorted5p, *sorted3p;
  
  debug5(printf("Starting Stage3_pair_up_concordant with %d concordant and %d samechr\n",
		*nconcordant,*nsamechr));

  /* Find sizes for allocating memory */
  debug5(printf("Sizes of 5-end pieces:\n"));
  nhits5_plus = (int *) CALLOC(cutoff_level_5+1,sizeof(int));
  nhits5_minus = (int *) CALLOC(cutoff_level_5+1,sizeof(int));
  for (i = 0; i < narray5; i++) {
    for (q = hitarray5[i]; q != NULL; q = q->rest) {
      hit5 = (T) q->first;
      if (hit5->chrnum == 0) {
	/* Cannot be paired */
      } else if (hit5->score > cutoff_level_5) {
      } else if (hit5->plusp) {
	nhits5_plus[hit5->score]++;
      } else {
	nhits5_minus[hit5->score]++;
      }
    }
  }

  debug5(
	 printf("Sizes of 5-end pieces:\n");
	 for (score5 = 0; score5 <= cutoff_level_5; score5++) {
	   printf("  score %d: %d plus, %d minus\n",score5,nhits5_plus[score5],nhits5_minus[score5]);
	 }
	 );

  debug5(printf("Sizes of 3-end pieces:\n"));
  nhits3_plus = (int *) CALLOC(cutoff_level_3+1,sizeof(int));
  nhits3_minus = (int *) CALLOC(cutoff_level_3+1,sizeof(int));
  for (i = 0; i < narray3; i++) {
    for (q = hitarray3[i]; q != NULL; q = q->rest) {
      hit3 = (T) q->first;
      if (hit3->chrnum == 0) {
	/* Cannot be paired */
      } else if (hit3->score > cutoff_level_3) {
      } else if (hit3->plusp) {
	nhits3_plus[hit3->score]++;
      } else {
	nhits3_minus[hit3->score]++;
      }
    }
  }
  debug5(
	 printf("Sizes of 3-end pieces:\n");
	 for (score3 = 0; score3 <= cutoff_level_3; score3++) {
	   printf("  score %d: %d plus, %d minus\n",score3,nhits3_plus[score3],nhits3_minus[score3]);
	 }
	 );

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
      if (hit5->chrnum == 0) {
	/* Cannot be paired */
      } else if (hit5->score > cutoff_level_5) {
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
      if (hit3->chrnum == 0) {
	/* Cannot be paired */
      } else if (hit3->score > cutoff_level_3) {
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
  while (*abort_pairing_p == false && pairscore <= *found_score + subopt_levels) {
    for (score5 = 0; score5 <= pairscore; score5++) {
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
	  debug5(printf("plus: i=%d %u, j=%d %u",i,hit5->genomicstart,j,hit3->genomicstart));

	  if (hit3->genomicstart + querylength3 < hit5->genomicstart) {
	    j++;
	  } else if (hit5->genomicstart + pairmax < hit3->genomicstart) {
	    i++;
	  } else {
	    for (jj = j; jj >= 0 && hits3[jj]->genomicstart + querylength3 >= hit5->genomicstart; jj--) ;
	    for (jj++; jj < nhits3 && hits3[jj]->genomicstart <= hit5->genomicstart + pairmax; jj++) {
	      hit3 = hits3[jj];

	      /* Only want pairs not previously seen */
	      if (hit5->paired_seenp == false || hit3->paired_seenp == false) {
		if (hit5->chrnum == hit3->chrnum) {
		  debug5(printf(" => concordant %d",hit5->chrnum));
		  hitpairs = List_push(hitpairs,Stage3pair_new(hit5,hit3,querylength5,querylength3,
							       expected_pairlength,/*pairtype*/CONCORDANT));
		  if (pairscore < *found_score) {
		    *found_score = pairscore;
		  }
		  if (++(*nconcordant) > maxpairedpaths) {
		    debug(printf("%d concordant paths exceeds %d\n",*nconcordant,maxpairedpaths));
		    *abort_pairing_p = true;
		  }
		  hit5->paired_seenp = true;
		  hit3->paired_seenp = true;
		}
	      }

	    }
	    i++;
	  }
	  debug5(printf("\n"));
	}

	/* hits3_minus (really on plus) against hits5_minus */
	i = j = 0;
	hits3 = hits3_minus[score3];
	hits5 = hits5_minus[score5];
	nhits3 = nhits3_minus[score3];
	nhits5 = nhits5_minus[score5];

	while (*abort_pairing_p == false && i < nhits3 && j < nhits5) {
	  hit3 = hits3[i];
	  hit5 = hits5[j];
	  debug5(printf("minus: i=%d %u, j=%d %u",i,hit3->genomicstart,j,hit5->genomicstart));

	  if (hit5->genomicstart + querylength5 < hit3->genomicstart) {
	    j++;
	  } else if (hit3->genomicstart + pairmax < hit5->genomicstart) {
	    i++;
	  } else {
	    for (jj = j; jj >= 0 && hits5[jj]->genomicstart + querylength5 >= hit3->genomicstart; jj--) ;
	    for (jj++; jj < nhits5 && hits5[jj]->genomicstart <= hit3->genomicstart + pairmax; jj++) {
	      hit5 = hits5[jj];

	      /* Only want pairs not previously seen */
	      if (hit3->paired_seenp == false || hit5->paired_seenp == false) {
		if (hit3->chrnum == hit5->chrnum) {
		  debug5(printf(" => concordant %d",hit3->chrnum));
		  hitpairs = List_push(hitpairs,Stage3pair_new(hit5,hit3,querylength5,querylength3,
							       expected_pairlength,/*pairtype*/CONCORDANT));
		  if (pairscore < *found_score) {
		    *found_score = pairscore;
		  }
		  if (++(*nconcordant) > maxpairedpaths) {
		    debug(printf("%d concordant paths exceeds %d\n",*nconcordant,maxpairedpaths));
		    *abort_pairing_p = true;
		  }
		  hit3->paired_seenp = true;
		  hit5->paired_seenp = true;
		}
	      }
	    }
	    i++;
	  }
	  debug5(printf("\n"));
	}
      }
    }
    pairscore++;
  }

  FREE(sorted3p);
  FREE(sorted5p);

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

  debug5(printf("Finished with Stage3_pair_up\n"));

  return hitpairs;
}


#if 0
/* Code to look for samechr */

  if (*abort_pairing_p == false && *nconcordant == 0) {
    /* Look for samechr pairs */

    nhits5 = nhits5_plus + nhits5_minus;
    nhits3 = nhits3_plus + nhits3_minus;
    debug5(printf("Looking for samechr in %d hits5 against %d hits3\n",nhits5,nhits3));

    /* Sort */
    if (nhits5 == 0) {
      hits5 = (T *) NULL;
    } else {
      hits5 = (T *) CALLOC(nhits5,sizeof(Stage3_T));
      j = 0;
      for (i = 0; i < nhits5_plus; i++) {
	hits5[j++] = hits5_plus[i];
      }
      for (i = 0; i < nhits5_minus; i++) {
	hits5[j++] = hits5_minus[i];
      }
      qsort(hits5,nhits5,sizeof(T),genomicstart_cmp);
    }

    if (nhits3 == 0) {
      hits3 = (T *) NULL;
    } else {
      hits3 = (T *) CALLOC(nhits3,sizeof(Stage3_T));
      j = 0;
      for (i = 0; i < nhits3_plus; i++) {
	hits3[j++] = hits3_plus[i];
      }
      for (i = 0; i < nhits3_minus; i++) {
	hits3[j++] = hits3_minus[i];
      }
      qsort(hits3,nhits3,sizeof(T),genomicstart_cmp);
    }

    /* hits5 against hits3 */
    i = j = 0;
    while (*nsamechr <= maxpairedpaths && i < nhits5 && j < nhits3) {
      hit5 = hits5[i];
      hit3 = hits3[j];
      debug5(printf("either: i=%d %u (chr %d, seenp %d), j=%d %u (chr %d, seenp %d)",
		    i,hit5->genomicstart,hit5->chrnum,hit5->paired_seenp,
		    j,hit3->genomicstart,hit3->chrnum,hit3->paired_seenp));
      
      if (hit3->chrnum < hit5->chrnum) {
	j++;
      } else if (hit5->chrnum < hit3->chrnum) {
	i++;
      } else if (hit5->concordantp == true) {
	/* Disregard a hit involved in a concordant pair */
	debug5(printf(" => hit5 concordant"));
	i++;
      } else {
	for (jj = j; jj < nhits3 && hits3[jj]->chrnum == hit5->chrnum; jj++) {
	  hit3 = hits3[jj];
	  if (hit3->concordantp == true) {
	    /* Disregard a hit involved in a concordant pair */
	    debug5(printf(" => hit3 concordant"));
	  } else if (hit5->paired_seenp == false || hit3->paired_seenp == false) {
	    /* Want only pairs not previously seen */
	    debug5(printf(" => samechr"));
	    hitpairs = List_push(hitpairs,Stage3pair_new(hit5,hit3,querylength5,querylength3,
							 expected_pairlength,/*pairtype*/SAMECHR));
	    if (++(*nsamechr) >= maxpairedpaths) {
	      debug(printf("%d samechr paths exceeds %d\n",*nsamechr,maxpairedpaths));
	    }
	  }
	}
	i++;
      }
      debug5(printf("\n"));
    }

    FREE(hits5);
    FREE(hits3);
  }

#endif

