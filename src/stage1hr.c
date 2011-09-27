static char rcsid[] = "$Id: stage1hr.c,v 1.248 2010/03/10 01:33:32 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
#define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "stage1hr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memset() */
#include <math.h>
#include <ctype.h>		/* for tolower() */
#include "mem.h"
#include "reader.h"
#include "oligo.h"
#include "indexdb.h"
#include "indexdb_hr.h"
#include "listdef.h"
#include "intlist.h"
#include "intlistdef.h"
#include "stage3hr.h"
#include "substring.h"
#include "complement.h"
#include "compress.h"
#include "genome_hr.h"
#include "maxent.h"
#include "iitdef.h"
#include "interval.h"
#include "spanningelt.h"
#include "cmet.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#endif


#define NREQUIRED_FAST 2	/* For candidate generation using
				   multimiss.  A value of 2 implies 
				   specificty of a 24-mer, which
				   should be low for a human-sized
				   genome */
#define NREQUIRED_STRETCH 1	/* For multimiss, avoiding multiple_mm */

#if INDEX1PART == 12
#define LEFTREADSHIFT 8		/* 32 - 2*INDEX1PART.  chars are shifted onto left of 32-bit word */
#define OLIGOBASE_MASK 0x00FFFFFF
#else
#define LEFTREADSHIFT (32 - INDEX1PART - INDEX1PART)
#define OLIGOBASE_MASK ~(~0U << 2*INDEX1PART)
#endif


#define NMOD 3			/* Sample q3 so 12-mers q3 -> 10-mers q1 */

/* On 5' end, x = querypos.  On 3' end, x = (query_lastpos - querypos). */
#define FLOOR_END(x) ((x < 3) ? 0 : (x + 9)/12)

/* Here, x = (querypos - last_querypos).  Was (x-3)/12, but the new formula handles indels. */
#define FLOOR_MIDDLE(x) ((x < 6) ? 0 : (x + 6)/12)

#define MIN_REPEAT 3

#define OPEN_SCORE -3
#define EXTEND_SCORE -1


#define HALF_INTRON_MAX_ENDLENGTH 14 /* 12-mer plus possible shift of 2 */

#define MAX_MIDDLE_INDELS_FOUND 100
#define MAX_LOCALSPLICING_HITS 10000
#define MAX_LOCALSPLICING_ATTEMPTS 50


/* Overall flow */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Diagonals for slow multiple mm, indels, splicing */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Indels */ 
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Indels, end */ 
#ifdef DEBUG2E
#define debug2e(x) x
#else
#define debug2e(x)
#endif

/* Floors */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Spliceends, known */ 
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* Spliceends, local pair */ 
#ifdef DEBUG4P
#define debug4p(x) x
#else
#define debug4p(x)
#endif

/* Splicepairs local (short) */ 
#ifdef DEBUG4S
#define debug4s(x) x
#else
#define debug4s(x)
#endif

/* Splicepairs novel distant (long) */ 
#ifdef DEBUG4L
#define debug4l(x) x
#else
#define debug4l(x)
#endif

/* Spliceends nove distant (long) details */
#ifdef DEBUG4LD
#define debug4ld(x) x
#else
#define debug4ld(x)
#endif

/* Spliceends for distant */
#ifdef DEBUG4E
#define debug4e(x) x
#else
#define debug4e(x)
#endif

/* Half introns */ 
#ifdef DEBUG4H
#define debug4h(x) x
#else
#define debug4h(x)
#endif

/* Heapify */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif

/* Identify exact/onemiss/multimiss matches */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Identify onemiss matches, list contents */
#ifdef DEBUG7A
#define debug7a(x) x
#else
#define debug7a(x)
#endif

/* Trim ends */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* straddling at beginning of genome.  May want to turn on DEBUG11 in indexdb_hr.c */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif

/* dual_search for known splice sites */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif

/* Tournament */ 
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif



typedef struct Indelsegment_T *Indelsegment_T;
struct Indelsegment_T {
  Genomicpos_T diagonal;
  Genomicpos_T chroffset;
  Genomicpos_T chrhigh;
  Chrnum_T chrnum;

  int querypos5;
  int querypos3;

  int floor_xfirst;
  int floor_xlast;
};

typedef struct Splicesegment_T *Splicesegment_T;
struct Splicesegment_T {
  Genomicpos_T diagonal;
  Genomicpos_T chroffset;
  Chrnum_T chrnum;

  int querypos5;
  int querypos3;
};



struct Floors_T {
  int *allocated0;
  int *prev_omitted;

  int **allocated2;
  int *allocated1;
  int **scorefrom;		/* [from][to] */

  int **allocated4;
  int *allocated3;
  int **scoreto;		/* [to][from] */
};


void
Floors_free (Floors_T *old) {
  FREE((*old)->allocated1);
  FREE((*old)->allocated2);

  FREE((*old)->allocated3);
  FREE((*old)->allocated4);

  if ((*old)->allocated0) {
    FREE((*old)->allocated0);
  }
  FREE(*old);
  return;
}

#ifdef DEBUG3
static void
Floors_print (Floors_T floors, int query_lastpos) {
  int from, to;

  if (floors->prev_omitted) {
    for (to = -3; to <= query_lastpos+3; to++) {
      printf("querypos %d, prev_omitted %d\n",to,floors->prev_omitted[to]);
    }
  }

  for (from = -3; from <= query_lastpos+3; from++) {
    for (to = from+1; to <= query_lastpos+3; to++) {
      printf("from %d to %d, floor_score %d or %d",
	     from,to,floors->scorefrom[from][to],floors->scoreto[to][from]);
      if (floors->prev_omitted) {
	printf(" (prev %d)",floors->prev_omitted[to]);
      }
      printf("\n");
    }
  }

  return;
}
#endif


static Floors_T
Floors_new_standard (int querylength, int max_end_insertions) {
  Floors_T new = (Floors_T) MALLOC(sizeof(*new));
  int query_lastpos, pos, from, to;
  int extra = 1 + max_end_insertions + max_end_insertions;

  query_lastpos = querylength - INDEX1PART;
  new->allocated0 = (int *) NULL;
  new->prev_omitted = (int *) NULL;

  new->allocated2 = (int **) CALLOC(query_lastpos+extra,sizeof(int *));
  new->allocated1 = (int *) CALLOC((query_lastpos+extra)*(query_lastpos+extra),sizeof(int));
  new->allocated2[0] = &(new->allocated1[max_end_insertions]);
  for (pos = 1; pos < query_lastpos+extra; pos++) {
    new->allocated2[pos] = &(new->allocated2[pos-1][query_lastpos+extra]);
  }
  new->scorefrom = &(new->allocated2[max_end_insertions]);


  new->allocated4 = (int **) CALLOC(query_lastpos+extra,sizeof(int *));
  new->allocated3 = (int *) CALLOC((query_lastpos+extra)*(query_lastpos+extra),sizeof(int));
  new->allocated4[0] = &(new->allocated3[max_end_insertions]);
  for (pos = 1; pos < query_lastpos+extra; pos++) {
    new->allocated4[pos] = &(new->allocated4[pos-1][query_lastpos+extra]);
  }
  new->scoreto = &(new->allocated4[max_end_insertions]);


  for (to = -max_end_insertions; to <= query_lastpos+max_end_insertions; to++) {
    for (from = -max_end_insertions; from < to; from++) {
      new->scorefrom[from][to] = new->scoreto[to][from] = FLOOR_MIDDLE(to - from);
    }
  }

  debug3(printf("Floors standard:\n"));
  debug3(Floors_print(new,query_lastpos));
  return new;
}


static Floors_T
Floors_new_omitted (int querylength, int max_end_insertions, bool *omitted) {
  Floors_T new = (Floors_T) MALLOC(sizeof(*new));
  int query_lastpos, querypos, pos, from, to;
  int prev;
  int extra = 1 + max_end_insertions + max_end_insertions;
  
  query_lastpos = querylength - INDEX1PART;
  new->allocated0 = (int *) CALLOC(query_lastpos+extra,sizeof(int));
  new->prev_omitted = &(new->allocated0[max_end_insertions]);

  new->allocated2 = (int **) CALLOC(query_lastpos+extra,sizeof(int *));
  new->allocated1 = (int *) CALLOC((query_lastpos+extra)*(query_lastpos+extra),sizeof(int));
  new->allocated2[0] = &(new->allocated1[max_end_insertions]);
  for (pos = 1; pos < query_lastpos+extra; pos++) {
    new->allocated2[pos] = &(new->allocated2[pos-1][query_lastpos+extra]);
  }
  new->scorefrom = &(new->allocated2[max_end_insertions]);


  new->allocated4 = (int **) CALLOC(query_lastpos+extra,sizeof(int *));
  new->allocated3 = (int *) CALLOC((query_lastpos+extra)*(query_lastpos+extra),sizeof(int));
  new->allocated4[0] = &(new->allocated3[max_end_insertions]);
  for (pos = 1; pos < query_lastpos+extra; pos++) {
    new->allocated4[pos] = &(new->allocated4[pos-1][query_lastpos+extra]);
  }
  new->scoreto = &(new->allocated4[max_end_insertions]);


  /* Set up omitted.  Save for middle_indels computation. */
  prev = -1;
  for (querypos = -max_end_insertions; querypos < 0; querypos++) {
    new->prev_omitted[querypos] = -1;
  }
  for ( ; querypos <= query_lastpos; querypos++) {
    new->prev_omitted[querypos] = prev;
    if (omitted[querypos] == true) {
      prev = querypos;
    }
  }
  for ( ; querypos <= query_lastpos+max_end_insertions; querypos++) {
    new->prev_omitted[querypos] = prev;
  }


  for (to = -max_end_insertions; to <= query_lastpos+max_end_insertions; to++) {
    prev = new->prev_omitted[to];
    for (from = -max_end_insertions; from < prev; from++) {
      new->scorefrom[from][to] = new->scorefrom[from][prev] + FLOOR_MIDDLE(to - prev);
    }
    for ( ; from < to; from++) {
      new->scorefrom[from][to] = FLOOR_MIDDLE(to - from);
    }
  }

  for (to = -max_end_insertions; to <= query_lastpos+max_end_insertions; to++) {
    for (from = -max_end_insertions; from < to; from++) {
      new->scoreto[to][from] = new->scorefrom[from][to];
    }
  }

  debug3(
	 printf("Floors omitted:");
	 for (pos = 0; pos <= query_lastpos; pos++) {
	   if (omitted[pos] == true) {
	     printf(" %d",pos);
	   }
	 }
	 printf("\n");
	 )
  debug3(Floors_print(new,query_lastpos));
  return new;
}


/************************************************************************/


#define T Stage1_T
struct T {
  List_T plus_spanningset[3];
  List_T minus_spanningset[3];

  Genomicpos_T **plus_positions_allocated;
  Genomicpos_T **plus_positions; /* points to above[2] */
  Genomicpos_T **minus_positions_allocated;
  Genomicpos_T **minus_positions; /* points to above[2] */

  int *plus_npositions_allocated;
  int *plus_npositions;		/* points to above[2] */

  int *minus_npositions_allocated;
  int *minus_npositions;	/* points to above[2] */

  bool *plus_retrievedp_allocated;
  bool *plus_retrievedp;	/* points to above[2] */
  bool *minus_retrievedp_allocated;
  bool *minus_retrievedp;	/* points to above[2] */

  bool *plus_allocp_allocated;
  bool *plus_allocp;		/* points to above[2] */
  bool *minus_allocp_allocated;
  bool *minus_allocp;		/* points to above[2] */

  bool *validp;
  bool *omitted;
  
  Storedoligomer_T *forward_oligos_allocated;
  Storedoligomer_T *forward_oligos; /* points to above[2] */
  Storedoligomer_T *revcomp_oligos_allocated;
  Storedoligomer_T *revcomp_oligos; /* points to above[2] */

  Compress_T query_compress_fwd;
  Compress_T query_compress_rev;

  bool all_positions_fetched_p;
};

static void
stage3list_gc (List_T *old) {
  List_T p;
  Stage3_T hit;

  for (p = *old; p != NULL; p = p->rest) {
    hit = (Stage3_T) p->first;
    Stage3_free(&hit);
  }
  List_free(&(*old));
  return;
}

static void
substringlist_gc (List_T *old) {
  List_T p;
  Substring_T hit;

  for (p = *old; p != NULL; p = p->rest) {
    hit = (Substring_T) p->first;
    Substring_free(&hit);
  }
  List_free(&(*old));
  return;
}
    
void
Stage1_free (T *old, int querylength) {
  List_T p;
  Spanningelt_T spanningelt;
  int mod, i;

  /* Stage1hr_check(*old); */

  if (*old) {
    if ((*old)->query_compress_fwd) {
      Compress_free(&(*old)->query_compress_fwd);
      Compress_free(&(*old)->query_compress_rev);
    }

    for (mod = 0; mod < 3; mod++) {
      for (p = (*old)->plus_spanningset[mod]; p; p = p->rest) {
	spanningelt = (Spanningelt_T) p->first;
	Spanningelt_free(&spanningelt);
      }
      List_free(&((*old)->plus_spanningset[mod]));
    
      for (p = (*old)->minus_spanningset[mod]; p; p = p->rest) {
	spanningelt = (Spanningelt_T) p->first;
	Spanningelt_free(&spanningelt);
      }
      List_free(&((*old)->minus_spanningset[mod]));
    }

    for (i = -2; i < querylength; i++) {
      if ((*old)->plus_positions[i] != NULL) {
	if ((*old)->plus_allocp[i] == true) {
	  FREE((*old)->plus_positions[i]);
	}
      }
      if ((*old)->minus_positions[i] != NULL) {
	if ((*old)->minus_allocp[i] == true) {
	  FREE((*old)->minus_positions[i]);
	}
      }
    }

    FREE((*old)->revcomp_oligos_allocated);
    FREE((*old)->forward_oligos_allocated);
    FREE((*old)->omitted);
    FREE((*old)->validp);
    FREE((*old)->plus_positions_allocated);
    FREE((*old)->minus_positions_allocated);
    FREE((*old)->plus_npositions_allocated);
    FREE((*old)->minus_npositions_allocated);
    FREE((*old)->plus_allocp_allocated);
    FREE((*old)->minus_allocp_allocated);
    FREE((*old)->plus_retrievedp_allocated);
    FREE((*old)->minus_retrievedp_allocated);

    FREE(*old);
  }

  return;
}


/************************************************************************/

static int
read_oligos (bool *allvalidp, T this, char *queryuc_ptr, int querylength,
	     int query_lastpos, bool dibasep, bool cmetp) {
  Reader_T reader;
  int querypos, noligos = 0;
  Oligostate_T last_state = INIT;
  Storedoligomer_T forward = 0U, revcomp = 0U;

  /* This estimate may be too high */
  /* this->maxfloor = 1 + querylength/oligobase * 2; */

  reader = Reader_new(queryuc_ptr,/*querystart*/0,/*queryend*/querylength,
		      /*reader_overlap*/INDEX1PART,dibasep);

  /* Prevents us from processing invalid query 12-mers */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    this->plus_retrievedp[querypos] = true;
    this->minus_retrievedp[querypos] = true;
  }

  /* Note: leftshifting is done here, rather than in Oligo_lookup */
  debug(printf("oligobase_mask: %08X\n",OLIGOBASE_MASK));
#if 0
  *any_omitted_p = false;
  *all_omitted_p = true;
#endif
  while ((last_state = Oligo_next(last_state,&querypos,&forward,&revcomp,INDEX1PART,
				  reader,/*cdnaend*/FIVE)) != DONE) {
    this->plus_positions[querypos] = (Genomicpos_T *) NULL;
    this->minus_positions[querypos] = (Genomicpos_T *) NULL;
    this->plus_npositions[querypos] = 0;
    this->minus_npositions[querypos] = 0;

    if (last_state == VALID) {
      this->validp[querypos] = true;

      this->plus_retrievedp[querypos] = false;
      this->minus_retrievedp[querypos] = false;

      if (dibasep) {
	this->forward_oligos[querypos] = forward & OLIGOBASE_MASK;
	this->revcomp_oligos[querypos] = (revcomp >> LEFTREADSHIFT) & OLIGOBASE_MASK;
      } else if (cmetp) {
	this->forward_oligos[querypos] = Cmet_reduce_ct(forward & OLIGOBASE_MASK);
	this->revcomp_oligos[querypos] = Cmet_reduce_ga((revcomp >> LEFTREADSHIFT) & OLIGOBASE_MASK);
      } else {
	this->forward_oligos[querypos] = forward & OLIGOBASE_MASK;
	this->revcomp_oligos[querypos] = (revcomp >> LEFTREADSHIFT) & OLIGOBASE_MASK;
      }

      debug(printf("At querypos %d, read oligo = %06X\n",querypos,this->forward_oligos[querypos]));
      noligos++;
    }
  }
  if (noligos < query_lastpos + 1) {
    debug(printf("Read only %d oligos due to non-ACGT; expected %d\n",noligos,query_lastpos + 1));
    *allvalidp = false;
  } else {
    *allvalidp = true;
  }

  Reader_free(&reader);

  return noligos;
}


/************************************************************************
 *   Omitted:
 *   In all cases, want to omit poly-AT.
 *   For purposes of finding mismatches, may want to omit frequent oligomers also
 *   For purposes of finding indels, may want to omit repetitive oligomers at ends also
 ************************************************************************/

static void
omit_oligos_clear (T this, int query_lastpos) {
  int querypos;

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    this->omitted[querypos] = false;
  }
  return;
}


static void
omit_oligos_polyat (bool *all_omitted_p, bool *any_omitted_p, T this, int query_lastpos) {
  int querypos;

  *all_omitted_p = true;
  *any_omitted_p = false;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->forward_oligos[querypos] == 0U || this->revcomp_oligos[querypos] == 0U) {
      this->omitted[querypos] = true;
      *any_omitted_p = true;
    } else {
      this->omitted[querypos] = false;
      *all_omitted_p = false;
    }
  }

  return;
}


static void
omit_oligos (bool *all_omitted_p, bool *any_omitted_p, T this, int query_lastpos,
	     int indexdb_size_threshold, bool frequentp, bool repetitivep) {
  int querypos;
  bool still_repetitive_p;

  *any_omitted_p = false;

  /* Always omit poly-AT */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->forward_oligos[querypos] == 0U || this->revcomp_oligos[querypos] == 0U) {
      this->omitted[querypos] = true;
      *any_omitted_p = true;
    } else {
      this->omitted[querypos] = false;
    }
  }

  if (frequentp == true) {
    /* Omit frequent oligos, but only in the middle */
    for (querypos = 3; querypos <= query_lastpos-3; querypos++) {
      if (this->plus_npositions[querypos] > indexdb_size_threshold && 
	  this->minus_npositions[querypos] > indexdb_size_threshold) {
	this->omitted[querypos] = true;
	*any_omitted_p = true;
      }
    }
  }

  if (repetitivep == true) {
    /* Omit repetitive oligos at the ends */
    still_repetitive_p = true;
    querypos = 0;
    while (querypos <= query_lastpos && still_repetitive_p == true) {
      if (Oligo_repetitive_p(this->forward_oligos[querypos])) {
	debug(printf("Querypos %d is repetitive\n",querypos));
	this->omitted[querypos] = true;
	*any_omitted_p = true;
      } else {
	still_repetitive_p = false;
      }
      querypos++;
    }

    still_repetitive_p = true;
    querypos = query_lastpos;
    while (querypos >= 0 && still_repetitive_p == true) {
      if (Oligo_repetitive_p(this->forward_oligos[querypos])) {
	debug(printf("Querypos %d is repetitive\n",querypos));
	this->omitted[querypos] = true;
	*any_omitted_p = true;
      } else {
	still_repetitive_p = false;
      }
      querypos--;
    }
  }

  if (*any_omitted_p == false) {
    debug(printf("No oligos are omitted\n"));
    *all_omitted_p = false;
  } else {
    debug(
	  printf("Omitted oligos:");
	  for (querypos = 0; querypos <= query_lastpos; querypos++) {
	    if (this->omitted[querypos] == true) {
	      printf(" %d",querypos);
	    }
	  }
	  printf("\n"));

    *all_omitted_p = true;
    for (querypos = 0; querypos <= query_lastpos; querypos++) {
      if (this->omitted[querypos] == false) {
	*all_omitted_p = false;
      }
    }
  }

  return;
}


static void
omit_oligos_repetitive (bool *all_omitted_p, bool *any_omitted_p, T this, int query_lastpos) {
  int querypos;

  *all_omitted_p = true;
  *any_omitted_p = false;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (Oligo_repetitive_p(this->forward_oligos[querypos])) {
      this->omitted[querypos] = true;
      *any_omitted_p = true;
    } else {
      this->omitted[querypos] = false;
      *all_omitted_p = false;
    }
  }

  return;
}


static T
Stage1_new (int querylength) {
  T new = (T) MALLOC(sizeof(*new));
  int querypos, mod;
  
  for (mod = 0; mod < 3; mod++) {
    new->plus_spanningset[mod] = (List_T) NULL;
    new->minus_spanningset[mod] = (List_T) NULL;
  }

  new->plus_positions_allocated = (Genomicpos_T **) CALLOC(querylength+2,sizeof(Genomicpos_T *));
  new->plus_positions = &(new->plus_positions_allocated[2]);
  new->minus_positions_allocated = (Genomicpos_T **) CALLOC(querylength+2,sizeof(Genomicpos_T *));
  new->minus_positions = &(new->minus_positions_allocated[2]);

  new->plus_npositions_allocated = (int *) CALLOC(querylength+2,sizeof(int));
  new->plus_npositions = &(new->plus_npositions_allocated[2]);
  new->minus_npositions_allocated = (int *) CALLOC(querylength+2,sizeof(int));
  new->minus_npositions = &(new->minus_npositions_allocated[2]);

  for (querypos = -2; querypos < querylength; querypos++) {
    new->plus_positions[querypos] = (Genomicpos_T *) NULL;
    new->plus_npositions[querypos] = 0;
    new->minus_positions[querypos] = (Genomicpos_T *) NULL;
    new->minus_npositions[querypos] = 0;
  }

  new->plus_retrievedp_allocated = (bool *) CALLOC(querylength+2,sizeof(bool));
  new->minus_retrievedp_allocated = (bool *) CALLOC(querylength+2,sizeof(bool));
  new->plus_retrievedp = &(new->plus_retrievedp_allocated[2]);
  new->minus_retrievedp = &(new->minus_retrievedp_allocated[2]);

  new->plus_allocp_allocated = (bool *) CALLOC(querylength+2,sizeof(bool));
  new->minus_allocp_allocated = (bool *) CALLOC(querylength+2,sizeof(bool));
  new->plus_allocp = &(new->plus_allocp_allocated[2]);
  new->minus_allocp = &(new->minus_allocp_allocated[2]);

  new->validp = (bool *) CALLOC(querylength,sizeof(bool));
  new->omitted = (bool *) CALLOC(querylength,sizeof(bool));

  new->forward_oligos_allocated = (Storedoligomer_T *) CALLOC(querylength+2,sizeof(Storedoligomer_T));
  new->forward_oligos = &(new->forward_oligos_allocated[2]);
  new->revcomp_oligos_allocated = (Storedoligomer_T *) CALLOC(querylength+2,sizeof(Storedoligomer_T));
  new->revcomp_oligos = &(new->revcomp_oligos_allocated[2]);

  new->query_compress_fwd = (Compress_T) NULL;
  new->query_compress_rev = (Compress_T) NULL;

  new->all_positions_fetched_p = false;

  return new;
}


/************************************************************************/

static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}

#if 0
static int
count_mismatches_met (int *nunknowns, int *nmetdiffs, int start, int end, char *gbuffer, char *query,
		      bool plusp) {
  int nmismatches = 0, i;
  
  debug12(printf("query:  %s\n",query));
  debug12(printf("genome: %s\n",gbuffer));
  debug12(printf("count:  "));
  debug12(
	 for (i = 0; i < start; i++) {
	   printf(" ");
	 }
	 );

  *nunknowns = 0;
  *nmetdiffs = 0;

  if (plusp) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'C' && query[i] == 'T') {
	debug12(printf("."));
	gbuffer[i] = '.';
	(*nmetdiffs)++;
      } else if (query[i] != gbuffer[i]) {
	debug12(printf("x"));
	if (gbuffer[i] == OUTOFBOUNDS) {
	  (*nunknowns)++;
	} else {
	  nmismatches++;
	}
	gbuffer[i] = tolower(gbuffer[i]);
      } else {
	debug12(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'G' && query[i] == 'A') {
	debug12(printf("."));
	gbuffer[i] = '.';
	(*nmetdiffs)++;
      } else if (query[i] != gbuffer[i]) {
	debug12(printf("x"));
	if (gbuffer[i] == OUTOFBOUNDS) {
	  (*nunknowns)++;
	} else {
	  nmismatches++;
	}
	gbuffer[i] = tolower(gbuffer[i]);
      } else {
	debug12(printf("*"));
      }
    }
  }

  debug12(printf("\n"));

  return nmismatches;
}

static void
mark_mismatches (int start, int end, char *gbuffer, char *query) {
  int i;
  
  debug12(printf("query:  %s\n",query));
  debug12(printf("genome: %s\n",gbuffer));
  debug12(printf("count:  "));
  debug12(
	 for (i = 0; i < start; i++) {
	   printf(" ");
	 }
	 );

  for (i = start; i < end; i++) {
    debug12(printf("*"));
    if (query[i] != gbuffer[i]) {
      gbuffer[i] = tolower(gbuffer[i]);
    }
  }

  debug12(printf("\n"));

  return;
}
#endif


/************************************************************************/

#define PARENT(i) (i >> 1)
#define LEFT(i) (i << 1)
#define RIGHT(i) ((i << 1) | 1)


typedef struct Batch_T *Batch_T;

struct Batch_T {
  int querypos;
  int diagterm;
  int npositions;

  Genomicpos_T diagonal;
  Genomicpos_T *positions;
#ifdef HAVE_64_BIT
  UINT8 diagonal_add_querypos;
#endif
};


static void
Batch_init (Batch_T batch, int querypos, int diagterm, Genomicpos_T *positions, int npositions, int querylength) {

  batch->querypos = querypos;
  batch->diagterm = diagterm;
  batch->positions = positions;
#ifdef WORDS_BIGENDIAN
  batch->diagonal = Bigendian_convert_uint(*positions) + diagterm;
#else
  batch->diagonal = *positions + diagterm;
#endif
  batch->npositions = npositions;

  while (batch->npositions > 0 && batch->diagonal < querylength) {
    debug11(printf("Eliminating diagonal %u as straddling beginning of genome (Batch_init)\n",batch->diagonal));
#ifdef WORDS_BIGENDIAN
    batch->diagonal = Bigendian_convert_uint(*(++batch->positions)) + diagterm;
#else
    batch->diagonal = *(++batch->positions) + diagterm;
#endif
    batch->npositions--;
  }

#ifdef HAVE_64_BIT
  batch->diagonal_add_querypos = (UINT8) batch->diagonal;
  batch->diagonal_add_querypos <<= 32;
  batch->diagonal_add_querypos |= querypos /* Previously added 2 because querypos was -2: + 2*/;
#endif

  return;
}


static void
Batch_init_simple (Batch_T batch, Genomicpos_T *diagonals, int ndiagonals, int querylength, int querypos) {

  batch->querypos = querypos;
  batch->positions = diagonals;
  batch->diagonal = *diagonals;
  batch->npositions = ndiagonals;

  while (batch->npositions > 0 && batch->diagonal < querylength) {
    debug11(printf("Eliminating diagonal %u as straddling beginning of genome (Batch_init)\n",batch->diagonal));
    batch->diagonal = *(++batch->positions);
    batch->npositions--;
  }

  return;
}


static void
min_heap_insert (Batch_T *heap, int *heapsize, Batch_T batch) {
  int i;
#ifdef HAVE_64_BIT
  UINT8 diagonal_add_querypos;
#else
  int querypos;
  Genomicpos_T diagonal;
#endif

  /* debug0(printf("Inserting into heap: diagonal %u at querypos %d\n",diagonal,querypos)); */

  i = ++(*heapsize);
#ifdef HAVE_64_BIT
  diagonal_add_querypos = batch->diagonal_add_querypos;
  while (i > 1 && (heap[PARENT(i)]->diagonal_add_querypos > diagonal_add_querypos)) {
    heap[i] = heap[PARENT(i)];
    i = PARENT(i);
  }
#else
  querypos = batch->querypos;
  diagonal = batch->diagonal;
  /* sort primarily by diagonal, then by querypos */
  while (i > 1 && (heap[PARENT(i)]->diagonal > diagonal || 
		   (heap[PARENT(i)]->diagonal == diagonal && heap[PARENT(i)]->querypos > querypos))) {
    heap[i] = heap[PARENT(i)];
    i = PARENT(i);
  }
#endif
  heap[i] = batch;

  return;
}


static void
min_heap_insert_simple (Batch_T *heap, int *heapsize, Batch_T batch) {
  int i;
  Genomicpos_T diagonal;

  /* debug0(printf("Inserting into heap: diagonal %u at querypos %d\n",diagonal,querypos)); */

  i = ++(*heapsize);
  diagonal = batch->diagonal;
  while (i > 1 && (heap[PARENT(i)]->diagonal > diagonal)) {
    heap[i] = heap[PARENT(i)];
    i = PARENT(i);
  }
  heap[i] = batch;

  return;
}



/* Note FORMULA: formulas for querypos <-> diagonal (diagterm in call to Indexdb_read) are:

plus: diagonal = position + querylength - querypos
minus: diagonal = position + querypos + INDEX1PART

For minus, the INDEX1PART is needed in call to Indexdb_read because
position is stored at beginning of plus oligomer, which corresponds to
end of minus oligomer.  As a result, we have the following formulas:

high genomic position = diagonal (corresponds to querypos =
querylength for plus, and querypos = 0 for minus)

low genomic position = diagonal - querylength (corresponds to querypos
= 0 for plus, and querypos = querylength for minus)

*/


static List_T
report_perfect_segment (int *found_score, List_T hits, Genomicpos_T diagonal,
			IIT_T chromosome_iit, char *query, char *queryptr, int querylength,
			char *gbuffer, Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			Genome_T genome, Genome_T genomealt,
			int nmisses_allowed, bool trim_ends_p, bool dibasep, bool cmetp, bool plusp) {
  Genomicpos_T left, chroffset, chrhigh;
  Chrnum_T chrnum;
  int nmismatches, ncolordiffs;

  left = diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
  chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,left,left);
  IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
  debug(printf("Reporting perfect segment at left %u and diagonal %u, with chroffset %u and chrhigh %u\n",
	       left,diagonal,chroffset,chrhigh));

  if (diagonal > chrhigh + 1) {
    /* Query goes over end of chromosome */
    debug(printf("  Ignore: goes over end of chromosome\n"));
    return hits;

  } else if (cmetp) {
    /* Count actual number of mismatches.  May not be a perfect segment. */
    nmismatches = Genome_count_mismatches_limit(&ncolordiffs,queryptr,query_compress,genome_blocks,snp_blocks,
						left,/*pos5*/0,/*pos3*/querylength,
						/*max_mismatches_allowed*/nmisses_allowed,
						dibasep,/*cmetp*/true,plusp);
    if (nmismatches > nmisses_allowed) {
      return hits;
    } else {
      /* Don't use Stage3_new_exact, because need to mark mismatches */
      Genome_fill_buffer_simple(genome,left,querylength,gbuffer);
      return List_push(hits,(void *) Stage3_new_substitution(&(*found_score),nmismatches,/*ncolordiffs*/0,
							     left,/*genomiclength*/querylength,
							     plusp,gbuffer,query,chrnum,chroffset,
							     trim_ends_p,dibasep,/*cmetp*/true));
    }

  } else if (genomealt) {
    Genome_fill_buffer_simple(genome,left,querylength,gbuffer);
    return List_push(hits,(void *) Stage3_new_substitution(&(*found_score),/*nmismatches*/0,/*ncolordiffs*/0,
							   left,/*genomiclength*/querylength,
							   plusp,gbuffer,query,chrnum,chroffset,
							   trim_ends_p,dibasep,/*cmetp*/false));
  } else if (dibasep) {
    Genome_count_mismatches_substring(&ncolordiffs,queryptr,query_compress,genome_blocks,snp_blocks,
				      left,/*pos5*/0,/*pos3*/querylength,dibasep,/*cmetp*/false,plusp);

    /* Need to fill buffer with nucleotide genome anyway */
    Genome_fill_buffer_simple(genome,left,querylength,gbuffer);
    return List_push(hits,(void *) Stage3_new_substitution(&(*found_score),/*nmismatches*/0,ncolordiffs,
							   left,/*genomiclength*/querylength,
							   plusp,gbuffer,query,chrnum,chroffset,
							   trim_ends_p,dibasep,/*cmetp*/false));

  } else {
    return List_push(hits,(void *) Stage3_new_exact(&(*found_score),left,/*genomiclength*/querylength,plusp,
						    chrnum,chroffset));
  }
}


/* Called only by exact/sub:1 procedures, so need to do Bigendian conversion */
#ifdef WORDS_BIGENDIAN
static int
binary_search (int lowi, int highi, Genomicpos_T *positions, Genomicpos_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = (lowi+highi)/2;
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,Bigendian_convert_uint(positions[lowi]),middlei,Bigendian_convert_uint(positions[middlei]),
		   highi,Bigendian_convert_uint(positions[highi]),goal));
    if (goal < Bigendian_convert_uint(positions[middlei])) {
      highi = middlei;
    } else if (goal > Bigendian_convert_uint(positions[middlei])) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}
#else
static int
binary_search (int lowi, int highi, Genomicpos_T *positions, Genomicpos_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = (lowi+highi)/2;
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi,positions[highi],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}
#endif


/* Generalization of identify_exact_iter and identify_onemiss_iter */
static List_T
identify_multimiss_iter (int *found_score, Chrnum_T *chrnum, Genomicpos_T *chroffset, Genomicpos_T *chrhigh,
			 List_T hits, Genomicpos_T goal, List_T prev, int *nempty, 
			 int *global_miss_querypos5, int *global_miss_querypos3,
			 char *query, char *queryptr, int querylength, char *gbuffer,
			 Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			 Genome_T genome, Genome_T genomealt, IIT_T chromosome_iit, bool plusp,
			 int nmisses_allowed, int nmisses_seen, int miss_querypos5, int miss_querypos3,
			 bool trim_ends_p, bool dibasep, bool cmetp) {
  List_T spanningset;
  void *ignore;
  Spanningelt_T elt;
  Compoundpos_T compoundpos;
  Genomicpos_T local_goal, left;
  Genomicpos_T position;
  int nmismatches, ncolordiffs, j;


  debug7(printf("identify_multimiss_iter on diagonal %u with %d misses seen initially\n",
		goal,nmisses_seen));

  if (nmisses_seen > nmisses_allowed) {
    debug7(printf("Result: skipping because %d misses seen > %d allowed\n",nmisses_seen,nmisses_allowed));
    return hits;
  }

  for (spanningset = prev->rest; spanningset /* != NULL */; prev = spanningset, spanningset = spanningset->rest) {
    elt = (Spanningelt_T) spanningset->first;
    debug7(printf("nmisses seen %d, allowed %d, remaining %d, goal %u: ",
		  nmisses_seen,nmisses_allowed,List_length(prev->rest),goal));

    if (elt->intersection_diagonals != NULL) {
      /* Intersection diagonals already computed */
      if (elt->intersection_ndiagonals > 0 && *elt->intersection_diagonals < goal) {
	debug7(printf("  (%d>>",elt->intersection_diagonals));
	j = 1;
	while (j < elt->intersection_ndiagonals && elt->intersection_diagonals[j] < goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= elt->intersection_ndiagonals) {
	  j = binary_search(j >> 1,elt->intersection_ndiagonals,elt->intersection_diagonals,goal);
	} else {
	  j = binary_search(j >> 1,j,elt->intersection_diagonals,goal);
	}
	elt->intersection_diagonals += j;
	elt->intersection_ndiagonals -= j;
	debug7(printf("  >>%d)",elt->intersection_ndiagonals));
      }
      
      if (elt->intersection_ndiagonals <= 0) {
	/* List is empty, so modify spanningset and continue with one more miss seen. */
	prev->rest = List_pop(spanningset,&ignore);
	spanningset = prev;
	*nempty += 1;
	if (elt->miss_querypos5 < *global_miss_querypos5) *global_miss_querypos5 = elt->miss_querypos5;
	if (elt->miss_querypos3 > *global_miss_querypos3) *global_miss_querypos3 = elt->miss_querypos3;

	debug7(printf(" intersection empty, counts as one miss --"));
	if (++nmisses_seen > nmisses_allowed) {
	  debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	  return hits;
	} else {
	  debug7(printf("  nmisses seen %d <= allowed %d, so continuing\n",nmisses_seen,nmisses_allowed));
	  if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	  if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	  /* continue; -- naturally falls to end of loop */
	}
      } else if (*elt->intersection_diagonals > local_goal) {
	/* Already advanced past goal, so continue with one more miss seen. */
	debug7(printf(" one miss --"));
	if (++nmisses_seen > nmisses_allowed) {
	  debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	  return hits;
	} else {
	  debug7(printf("  nmisses seen %d <= allowed %d, so continuing\n",nmisses_seen,nmisses_allowed));
	  if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	  if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	  /* continue; -- naturally falls to end of loop */
	}
      } else {
	/* Found goal.  Advance past goal and continue with loop. */
	debug7(printf(" advancing\n"));
	++elt->intersection_diagonals;
	--elt->intersection_ndiagonals;
	/* continue; -- naturally falls to end of loop */
      }

    } else {
      if (elt->partnerp == true) {
	/* Partner is guaranteed to be atomic */
	local_goal = goal - elt->partner_diagterm;

#ifdef WORDS_BIGENDIAN
	if (elt->partner_npositions > 0 && Bigendian_convert_uint(*elt->partner_positions) < local_goal) {
	  debug7(printf("  (%d>>",elt->partner_npositions));
	  j = 1;
	  while (j < elt->partner_npositions && Bigendian_convert_uint(elt->partner_positions[j]) < local_goal) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= elt->partner_npositions) {
	    j = binary_search(j >> 1,elt->partner_npositions,elt->partner_positions,local_goal);
	  } else {
	    j = binary_search(j >> 1,j,elt->partner_positions,local_goal);
	  }
	  elt->partner_positions += j;
	  elt->partner_npositions -= j;
	  debug7(printf("  >>%d)",elt->partner_npositions));
	}
#else
	if (elt->partner_npositions > 0 && *elt->partner_positions < local_goal) {
	  debug7(printf("  (%d>>",elt->partner_npositions));
	  j = 1;
	  while (j < elt->partner_npositions && elt->partner_positions[j] < local_goal) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= elt->partner_npositions) {
	    j = binary_search(j >> 1,elt->partner_npositions,elt->partner_positions,local_goal);
	  } else {
	    j = binary_search(j >> 1,j,elt->partner_positions,local_goal);
	  }
	  elt->partner_positions += j;
	  elt->partner_npositions -= j;
	  debug7(printf("  >>%d)",elt->partner_npositions));
	}
#endif

	if (elt->partner_npositions <= 0) {
	  /* Empty, so modify spanningset and continue with one more miss seen. */
	  prev->rest = List_pop(spanningset,&ignore);
	  spanningset = prev;
	  *nempty += 1;
	  if (elt->miss_querypos5 < *global_miss_querypos5) *global_miss_querypos5 = elt->miss_querypos5;
	  if (elt->miss_querypos3 > *global_miss_querypos3) *global_miss_querypos3 = elt->miss_querypos3;

	  debug7(printf(" partner empty --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    continue;		/* Don't need to check main list below */
	  }
#ifdef WORDS_BIGENDIAN
	} else if (Bigendian_convert_uint(*elt->partner_positions) > local_goal) {
	  /* Advanced past local_goal, so continue with one more miss seen. */
	  debug7(printf(" not in partner --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    continue;		/* Don't need to check main list below */
	  }
#else
	} else if (*elt->partner_positions > local_goal) {
	  /* Advanced past local_goal, so continue with one more miss seen. */
	  debug7(printf(" not in partner --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    continue;		/* Don't need to check main list below */
	  }
#endif
	} else {
	  /* Found local_goal.  Advance past local_goal and continue with rest of compound querypos */
	  debug7(printf(" found in partner, so continue with rest of compound querypos\n"));
	  ++elt->partner_positions;
	  --elt->partner_npositions;
	  /* Continue below with main list */
	}
      }

      if ((compoundpos = elt->compoundpos) != NULL) {
	local_goal = goal - elt->compoundpos_diagterm;
	if (Compoundpos_search(&position,compoundpos,local_goal) <= 0) {
	  /* Empty, so modify spanningset and continue with one more miss seen. */
	  prev->rest = List_pop(spanningset,&ignore);
	  spanningset = prev;
	  *nempty += 1;
	  if (elt->miss_querypos5 < *global_miss_querypos5) *global_miss_querypos5 = elt->miss_querypos5;
	  if (elt->miss_querypos3 > *global_miss_querypos3) *global_miss_querypos3 = elt->miss_querypos3;

	  debug7(printf("  compoundpos empty --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
	} else if (position > local_goal) {
	  /* Advanced past goal.  Continue with one more miss seen. */
	  debug7(printf("  compoundpos failed %u > %u --",position,local_goal));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
	} else {
	  /* Found goal.  Advance past goal and continue with loop.  */
	  debug7(printf("  found %u, advancing...",local_goal));
	  /* continue; -- Naturally falls to end of loop */
	}
      
      } else {
	/* Ordinary querypos */
	local_goal = goal - elt->diagterm;

#ifdef WORDS_BIGENDIAN
	if (elt->npositions > 0 && Bigendian_convert_uint(*elt->positions) < local_goal) {
	  debug7(printf("  (%d>>",elt->npositions));
	  j = 1;
	  while (j < elt->npositions && Bigendian_convert_uint(elt->positions[j]) < local_goal) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= elt->npositions) {
	    j = binary_search(j >> 1,elt->npositions,elt->positions,local_goal);
	  } else {
	    j = binary_search(j >> 1,j,elt->positions,local_goal);
	  }
	  elt->positions += j;
	  elt->npositions -= j;
	  debug7(printf("  >>%d)",elt->npositions));
	}
#else
	if (elt->npositions > 0 && *elt->positions < local_goal) {
	  debug7(printf("  (%d>>",elt->npositions));
	  j = 1;
	  while (j < elt->npositions && elt->positions[j] < local_goal) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= elt->npositions) {
	    j = binary_search(j >> 1,elt->npositions,elt->positions,local_goal);
	  } else {
	    j = binary_search(j >> 1,j,elt->positions,local_goal);
	  }
	  elt->positions += j;
	  elt->npositions -= j;
	  debug7(printf("  >>%d)",elt->npositions));
	}
#endif

	if (elt->npositions <= 0) {
	  /* List is empty, so continue with one more miss seen. */
	  prev->rest = List_pop(spanningset,&ignore);
	  spanningset = prev;
	  *nempty += 1;
	  if (elt->miss_querypos5 < *global_miss_querypos5) *global_miss_querypos5 = elt->miss_querypos5;
	  if (elt->miss_querypos3 > *global_miss_querypos3) *global_miss_querypos3 = elt->miss_querypos3;

	  debug7(printf(" positions empty, counts as one miss --"));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
#ifdef WORDS_BIGENDIAN
	} else if (Bigendian_convert_uint(*elt->positions) > local_goal) {
	  /* Already advanced past goal, so continue with one more miss seen. */
	  debug7(printf(" one miss %u > %u --",Bigendian_convert_uint(*elt->positions),local_goal));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
#else
	} else if (*elt->positions > local_goal) {
	  /* Already advanced past goal, so continue with one more miss seen. */
	  debug7(printf(" one miss %u > %u --",*elt->positions,local_goal));
	  if (++nmisses_seen > nmisses_allowed) {
	    debug7(printf(" nmisses seen %d > allowed %d, so returning\n",nmisses_seen,nmisses_allowed));
	    return hits;
	  } else {
	    debug7(printf("  nmisses seen %d <= allowed %d, so continuing\n",nmisses_seen,nmisses_allowed));
	    if (elt->miss_querypos5 < miss_querypos5) miss_querypos5 = elt->miss_querypos5;
	    if (elt->miss_querypos3 > miss_querypos3) miss_querypos3 = elt->miss_querypos3;
	    /* continue; -- Naturally falls to end of loop */
	  }
#endif
	} else {
	  /* Found goal.  Advance past goal and continue with loop. */
	  debug7(printf(" advancing\n"));
	  ++elt->positions;
	  --elt->npositions;
	  /* continue; -- Naturally falls to end of loop */
	}
      }
    }
    /* End of loop */
  }
    
  /* success */
  debug7(printf("  successful candidate found, with >= %d misses, %d allowed\n",nmisses_seen,nmisses_allowed));
  if (nmisses_seen == 0) {
    return report_perfect_segment(&(*found_score),hits,/*diagonal*/goal,
				  chromosome_iit,query,queryptr,querylength,
				  gbuffer,query_compress,genome_blocks,snp_blocks,
				  genome,genomealt,nmisses_allowed,trim_ends_p,dibasep,cmetp,plusp);
  } else {
    if (goal < querylength) {
      debug7(printf("  Goes over beginning of chromosome\n"));
      return hits;
    } else {
      left = goal - querylength;	/* goal here is diagonal */
    }

    if (goal > *chrhigh) {
      *chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,left,left);
      IIT_interval_bounds(&(*chroffset),&(*chrhigh),chromosome_iit,*chrnum);
      *chrhigh += 1U;
    }
    if (goal > *chrhigh) {
      debug7(printf("  Goes over end of chromosome\n"));
      return hits;

    } else {
      if (cmetp || snp_blocks) {
	debug7(printf("  Testing in entire query\n"));
	nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryptr,query_compress,genome_blocks,snp_blocks,
							left,/*pos5*/0,/*pos3*/querylength,dibasep,cmetp,plusp);
      } else {
	debug7(printf("  Testing in query bounds %d..%d\n",miss_querypos5,miss_querypos3));
	nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryptr,query_compress,genome_blocks,snp_blocks,
							left,/*pos5*/miss_querypos5,/*pos3*/miss_querypos3,
							dibasep,cmetp,plusp);

      }
      debug7(printf("nmismatches = %d (vs %d misses allowed)\n",nmismatches,nmisses_allowed));

      if (nmismatches > nmisses_allowed) {
	debug7(printf("Result: too many mismatches\n"));
	return hits;
      } else {
	debug7(printf("Result: successful hit saved\n"));
	Genome_fill_buffer_simple(genome,left,querylength,gbuffer);
	debug(printf("Reporting hit with %d mismatches\n",nmismatches));
	return List_push(hits,(void *) Stage3_new_substitution(&(*found_score),nmismatches,ncolordiffs,
							       left,/*genomiclength*/querylength,
							       plusp,gbuffer,query,*chrnum,*chroffset,
							       trim_ends_p,dibasep,cmetp));
      }
    }
  }
}


/* Since querypos -2, -1, query_lastpos+1, and query_lastpos+2 are now
   stored as compoundpos, we no longer want to use them for boosting */
static void
most_specific_oligomer (int *best_plus_querypos, int *best_minus_querypos, T this,
			int query_lastpos, Indexdb_T indexdb, Indexdb_T indexdb2) {
  int querypos, mod;
  int best_plus_count[3], best_minus_count[3];

  /* Not needed, since this is the first procedure called */
  /* Block_restore(this->block5); */

  best_plus_querypos[0] = -3;
  best_plus_querypos[1] = -3;
  best_plus_querypos[2] = -3;
  best_minus_querypos[0] = -3;
  best_minus_querypos[1] = -3;
  best_minus_querypos[2] = -3;

  best_plus_count[0] = best_plus_count[1] = best_plus_count[2] = 0;
  best_minus_count[0] = best_minus_count[1] = best_minus_count[2] = 0;

#if 0
  if (this->validp[0] == false) {
    debug(printf("Not counting at querypos 0, neg 2 or neg 1 because validp is false at querypos 0\n"));
    this->plus_npositions[-2] = 0;
    this->plus_npositions[-1] = 0;
    this->minus_npositions[-2] = 0;
    this->minus_npositions[-1] = 0;
  } else {
#endif
    this->plus_npositions[-2] = Indexdb_count_left_subst_2(indexdb,this->forward_oligos[0]);
    this->minus_npositions[-2] = Indexdb_count_right_subst_2(indexdb2,this->revcomp_oligos[0]);
    debug(printf("Counting at querypos 0, neg 2, plus_npositions = %d (oligo %06X), minus_npositions = %d (oligo %06X)\n",
		 this->plus_npositions[-2],this->forward_oligos[0],this->minus_npositions[-2],this->revcomp_oligos[0]));
    best_plus_count[1] = this->plus_npositions[-2];
    best_minus_count[1] = this->minus_npositions[-2];

    this->plus_npositions[-1] = Indexdb_count_left_subst_1(indexdb,this->forward_oligos[0]);
    this->minus_npositions[-1] = Indexdb_count_right_subst_1(indexdb2,this->revcomp_oligos[0]);
    debug(printf("Counting at querypos 0, neg 1, plus_npositions = %d (oligo %06X), minus_npositions = %d (oligo %06X)\n",
		 this->plus_npositions[-1],this->forward_oligos[0],this->minus_npositions[-1],this->revcomp_oligos[0]));
    best_plus_count[2] = this->plus_npositions[-1];
    best_minus_count[2] = this->minus_npositions[-1];
#if 0
  }
#endif

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
#if 0
    if (this->validp[querypos] == true) {
#endif
      mod = querypos % 3;
      this->plus_npositions[querypos] = Indexdb_count_no_subst(indexdb,this->forward_oligos[querypos]);
      this->minus_npositions[querypos] = Indexdb_count_no_subst(indexdb2,this->revcomp_oligos[querypos]);
      debug(printf("Counting at querypos %d, plus_npositions = %d (oligo %06X), minus_npositions = %d (oligo %06X)\n",
		   querypos,this->plus_npositions[querypos],this->forward_oligos[querypos],
		   this->minus_npositions[querypos],this->revcomp_oligos[querypos]));

      if (best_plus_querypos[mod] < 0 || this->plus_npositions[querypos] < best_plus_count[mod]) {
	best_plus_querypos[mod] = querypos;
	best_plus_count[mod] = this->plus_npositions[querypos];
      }
      if (best_minus_querypos[mod] < 0 || this->minus_npositions[querypos] < best_minus_count[mod]) {
	best_minus_querypos[mod] = querypos;
	best_minus_count[mod] = this->minus_npositions[querypos];
      }
#if 0
    }
#endif
  }

  querypos = query_lastpos;
#if 0
  if (this->validp[querypos] == false) {
    debug(printf("Not counting at querypos %d (pos 1) or %d (pos 2) because validp is false at querypos %d\n",
		 querypos+1,querypos+2,querypos));
    this->plus_npositions[querypos+1] = 0;
    this->plus_npositions[querypos+2] = 0;
    this->minus_npositions[querypos+1] = 0;
    this->minus_npositions[querypos+2] = 0;
  } else {
#endif
    mod = (querypos+1) % 3;
    this->plus_npositions[querypos+1] = Indexdb_count_right_subst_1(indexdb,this->forward_oligos[querypos]);
    this->minus_npositions[querypos+1] = Indexdb_count_left_subst_1(indexdb2,this->revcomp_oligos[querypos]);
    debug(printf("Counting at querypos %d, pos 1, plus_npositions = %d (oligo %06X), minus_npositions = %d (oligo %06X)\n",
		 querypos,this->plus_npositions[querypos+1],this->forward_oligos[querypos],
		 this->minus_npositions[querypos+1],this->revcomp_oligos[querypos]));

#if 0
    /* Don't want boostpos to be a compoundpos */
    if (best_plus_querypos[mod] < 0 || this->plus_npositions[querypos+1] < best_plus_count[mod]) {
      best_plus_querypos[mod] = querypos+1;
      best_plus_count[mod] = this->plus_npositions[querypos+1];
    }
    if (best_minus_querypos[mod] < 0 || this->minus_npositions[querypos+1] < best_minus_count[mod]) {
      best_minus_querypos[mod] = querypos+1;
      best_minus_count[mod] = this->minus_npositions[querypos+1];
    }
#endif

    mod = (querypos+2) % 3;
    this->plus_npositions[querypos+2] = Indexdb_count_right_subst_2(indexdb,this->forward_oligos[querypos]);
    this->minus_npositions[querypos+2] = Indexdb_count_left_subst_2(indexdb2,this->revcomp_oligos[querypos]);
    debug(printf("Counting at querypos %d, pos 2, plus_npositions = %d (oligo %06X), minus_npositions = %d (oligo %06X)\n",
		 querypos,this->plus_npositions[querypos+2],this->forward_oligos[querypos],
		 this->minus_npositions[querypos+2],this->revcomp_oligos[querypos]));

#if 0
    /* Don't want boostpos to be a compoundpos */
    if (best_plus_querypos[mod] < 0 || this->plus_npositions[querypos+2] < best_plus_count[mod]) {
      best_plus_querypos[mod] = querypos+2;
      best_plus_count[mod] = this->plus_npositions[querypos+2];
    }
    if (best_minus_querypos[mod] < 0 || this->minus_npositions[querypos+2] < best_minus_count[mod]) {
      best_minus_querypos[mod] = querypos+2;
      best_minus_count[mod] = this->minus_npositions[querypos+2];
    }
#endif

#if 0
  }
#endif

  return;
}


static List_T
find_spanning_exact_matches (int *found_score, T this,
			     char *queryuc_ptr, char *queryrc, int querylength, int query_lastpos,
			     Indexdb_T indexdb, Indexdb_T indexdb2, char *gbuffer,
			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome, Genome_T genomealt,
			     IIT_T chromosome_iit, bool trim_ends_p, bool dibasep, bool cmetp) {
  List_T hits = NULL;
  List_T spanningset, sorted;
  Spanningelt_T *array;
  int best_plus_querypos[3], best_minus_querypos[3];
  Genomicpos_T *diagonals0, *positions0, diagonal0;
  int diagterm0, ndiagonals0, npositions0;
  int boostpos, mod, nelts, minscore, i;
  int global_miss_querypos5, global_miss_querypos3, elt_miss_querypos5, elt_miss_querypos3;
  int nempty;
  Chrnum_T chrnum;
  Genomicpos_T chroffset, chrhigh;

  debug(printf("Starting find_spanning_exact_matches\n"));

  /* Use shortest list for candidate generation */
  most_specific_oligomer(best_plus_querypos,best_minus_querypos,this,query_lastpos,indexdb,indexdb2);

  /* Plus */
  for (mod = 0; mod < 3; mod++) {
    chrhigh = 0U;
    spanningset = Spanningelt_set(&minscore,this->forward_oligos,&this->plus_retrievedp,&this->plus_positions,
				  this->plus_npositions,indexdb,query_lastpos,querylength,mod,/*plusp*/true);
    nelts = List_length(spanningset);
    array = (Spanningelt_T *) List_to_array(spanningset,NULL);
    List_free(&spanningset);

    boostpos = best_plus_querypos[mod];
    debug(printf("exact_matches, plus mod %d: proposed boostpos is %d\n",mod,boostpos));
    if (this->plus_npositions[boostpos] < minscore &&
	this->plus_retrievedp[boostpos] == false) {
      /* Boost */
      qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
      sorted = (List_T) NULL;
      for (i = nelts-1; i >= 0; --i) {
	sorted = List_push(sorted,array[i]);
      }
      FREE(array);
      this->plus_spanningset[mod] = sorted;

      /* Get boost positions */
      this->plus_positions[boostpos] = 
	Indexdb_read_inplace(&(this->plus_npositions[boostpos]),indexdb,this->forward_oligos[boostpos]);
      this->plus_retrievedp[boostpos] = true;
      positions0 = this->plus_positions[boostpos];
      npositions0 = this->plus_npositions[boostpos];
      diagterm0 = querylength - boostpos; /* FORMULA */

      debug(printf("*** find_spanning_exact_matches, plus mod %d, with boost @ %d (%d positions)\n",
		   mod,boostpos,npositions0);
	    Spanningelt_print_set(sorted));

      spanningset = List_push(List_copy(sorted),(void **) NULL); /* Add a dummy list elt to front */
      nempty = 0;
      global_miss_querypos5 = querylength;
      global_miss_querypos3 = 0;

      while (--npositions0 >= 0 && nempty == 0) {
#ifdef WORDS_BIGENDIAN
	diagonal0 = Bigendian_convert_uint(*positions0++) + diagterm0;
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,Bigendian_convert_uint(*positions0),diagterm0));
#else
	diagonal0 = (*positions0++) + diagterm0;
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,(*positions0),diagterm0));
#endif
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/true,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3,
				       trim_ends_p,dibasep,cmetp);
      }
      List_free(&spanningset);

    } else {
      qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
      if (nelts > 1) {
	qsort(&(array[1]),nelts-1,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
      }
      sorted = (List_T) NULL;
      for (i = nelts-1; i >= 0; --i) {
	sorted = List_push(sorted,array[i]);
      }
      FREE(array);
      this->plus_spanningset[mod] = sorted;

      debug(printf("*** find_spanning_exact_matches, plus mod %d, no boosting\n",mod));
      debug(Spanningelt_print_set(this->plus_spanningset[mod]));

      diagonals0 = Spanningelt_diagonals(&ndiagonals0,(Spanningelt_T) sorted->first,&elt_miss_querypos5,&elt_miss_querypos3);
      spanningset = List_push(List_copy(sorted->rest),(void **) NULL); /* Add a dummy list elt to front */
      nempty = 0;
      global_miss_querypos5 = querylength;
      global_miss_querypos3 = 0;

      while (--ndiagonals0 >= 0 && nempty == 0) {
	debug7(printf("diag0 %d:%u advancing\n",ndiagonals0,(*diagonals0)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,*diagonals0++,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/true,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3,
				       trim_ends_p,dibasep,cmetp);
      }
      List_free(&spanningset);
    }
  }

  /* Minus */
  for (mod = 0; mod < 3; mod++) {
    chrhigh = 0U;
    spanningset = Spanningelt_set(&minscore,this->revcomp_oligos,&this->minus_retrievedp,&this->minus_positions,
				  this->minus_npositions,indexdb2,query_lastpos,querylength,mod,/*plusp*/false);
    nelts = List_length(spanningset);
    array = (Spanningelt_T *) List_to_array(spanningset,NULL);
    List_free(&spanningset);

    boostpos = best_minus_querypos[mod];
    debug(printf("exact_matches, minus mod %d: proposed boostpos is %d\n",mod,boostpos));
    if (this->minus_npositions[boostpos] < minscore &&
	this->minus_retrievedp[boostpos] == false) {
      /* Boost */
      qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
      sorted = (List_T) NULL;
      for (i = nelts-1; i >= 0; --i) {
	sorted = List_push(sorted,array[i]);
      }
      FREE(array);
      this->minus_spanningset[mod] = sorted;

      /* Get boost positions */
      this->minus_positions[boostpos] = 
	Indexdb_read_inplace(&(this->minus_npositions[boostpos]),indexdb2,this->revcomp_oligos[boostpos]);
      this->minus_retrievedp[boostpos] = true;
      positions0 = this->minus_positions[boostpos];
      npositions0 = this->minus_npositions[boostpos];
      diagterm0 = boostpos + INDEX1PART; /* FORMULA */

      debug(printf("*** find_spanning_exact_matches, minus mod %d, with boost @ %d (%d positions)\n",
		   mod,boostpos,npositions0);
	    Spanningelt_print_set(sorted));

      spanningset = List_push(List_copy(sorted),(void **) NULL);/* Add a dummy list elt to front */
      nempty = 0;
      global_miss_querypos5 = querylength;
      global_miss_querypos3 = 0;

      while (--npositions0 >= 0 && nempty == 0) {
#ifdef WORDS_BIGENDIAN
	diagonal0 = Bigendian_convert_uint(*positions0++) + diagterm0;
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,Bigendian_convert_uint(*positions0),diagterm0));
#else
	diagonal0 = (*positions0++) + diagterm0;
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,(*positions0),diagterm0));
#endif
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/false,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3,
				       trim_ends_p,dibasep,cmetp);
      }
      List_free(&spanningset);

    } else {
      qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
      if (nelts > 1) {
	qsort(&(array[1]),nelts-1,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
      }
      sorted = (List_T) NULL;
      for (i = nelts-1; i >= 0; --i) {
	sorted = List_push(sorted,array[i]);
      }
      FREE(array);
      this->minus_spanningset[mod] = sorted;

      debug(printf("*** find_spanning_exact_matches, minus mod %d, no boosting\n",mod));
      debug(Spanningelt_print_set(this->minus_spanningset[mod]));

      diagonals0 = Spanningelt_diagonals(&ndiagonals0,(Spanningelt_T) sorted->first,&elt_miss_querypos5,&elt_miss_querypos3);
      spanningset = List_push(List_copy(sorted->rest),(void **) NULL); /* Add a dummy list elt to front */
      nempty = 0;
      global_miss_querypos5 = querylength;
      global_miss_querypos3 = 0;

      while (--ndiagonals0 >= 0 && nempty == 0) {
	debug7(printf("diag0 %d:%u advancing\n",ndiagonals0,(*diagonals0)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,*diagonals0++,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/false,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3,
				       trim_ends_p,dibasep,cmetp);
      }
      List_free(&spanningset);
    }
  }
    
  return hits;
}


static List_T
find_spanning_onemiss_matches (int *found_score, List_T hits, T this,
			       char *queryuc_ptr, char *queryrc, int querylength,
			       char *gbuffer, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			       UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome, Genome_T genomealt,
			       IIT_T chromosome_iit, bool trim_ends_p, bool dibasep, bool cmetp) {
  List_T spanningset, sorted;
  Spanningelt_T *array;
  Genomicpos_T *diagonals0, *diagonals1, diagonal0, diagonal1;
  int global_miss_querypos5, global_miss_querypos3;
  int miss0_querypos5, miss0_querypos3, miss1_querypos5, miss1_querypos3;
  int mod, nelts, i;
  int ndiagonals0, ndiagonals1;
  int nempty;
  Chrnum_T chrnum;
  Genomicpos_T chroffset, chrhigh;

  debug(printf("Starting find_spanning_onemiss_matches\n"));

  /* Plus */
  for (mod = 0; mod < 3; mod++) {
    debug(printf("Onemiss plus mod %d\n",mod));

    spanningset = this->plus_spanningset[mod];
    nelts = List_length(spanningset);
    array = (Spanningelt_T *) List_to_array(spanningset,NULL);
    /* List_free(&spanningset); */

    qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
    if (nelts > 2) {
      qsort(&(array[2]),nelts-2,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
    }
    sorted = (List_T) NULL;
    for (i = nelts-1; i >= 0; --i) {
      sorted = List_push(sorted,Spanningelt_reset(array[i]));
    }
    FREE(array);

    debug(printf("*** find_spanning_onemiss_matches, plus mod %d\n",mod));
    debug(Spanningelt_print_set(sorted));

    diagonals0 = Spanningelt_diagonals(&ndiagonals0,(Spanningelt_T) sorted->first,&miss0_querypos5,&miss0_querypos3);
    diagonals1 = Spanningelt_diagonals(&ndiagonals1,(Spanningelt_T) sorted->rest->first,&miss1_querypos5,&miss1_querypos3);
    spanningset = List_push(List_copy(sorted->rest->rest),(void **) NULL); /* Add a dummy list elt to front */
    nempty = 0;
    global_miss_querypos5 = querylength;
    global_miss_querypos3 = 0;
    List_free(&sorted);
    chrhigh = 0U;

    while (ndiagonals0 > 0 && ndiagonals1 > 0 && nempty <= 1) {
      if ((diagonal0 = (*diagonals0)) < (diagonal1 = (*diagonals1))) {
	debug7(printf("diag0 %d:%u advancing\n",ndiagonals0,(*diagonals0)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/true,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3,
				       trim_ends_p,dibasep,cmetp);
	++diagonals0;
	--ndiagonals0;

      } else if (diagonal1 < diagonal0) {
	debug7(printf("diag1 %d:%u advancing\n",ndiagonals1,(*diagonals1)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal1,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/true,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3,
				       trim_ends_p,dibasep,cmetp);
	++diagonals1;
	--ndiagonals1;

      } else {
	debug7(printf("diag0&1 %d:%u == %d:%u advancing\n",ndiagonals0,diagonal0,ndiagonals1,diagonal1));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/true,/*nmisses_allowed*/1,
				       /*nmisses_seen*/nempty,global_miss_querypos5,global_miss_querypos3,
				       trim_ends_p,dibasep,cmetp);
	++diagonals0;
	--ndiagonals0;
	++diagonals1;
	--ndiagonals1;
      }
    }

    while (--ndiagonals0 >= 0 && nempty == 0) {
      debug7(printf("diag0 %d:%u advancing\n",ndiagonals0+1,(*diagonals0)));
      hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,*diagonals0++,
				     /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				     /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
				     /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				     genome,genomealt,chromosome_iit,/*plusp*/true,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3,
				     trim_ends_p,dibasep,cmetp);
    }

    while (--ndiagonals1 >= 0 && nempty == 0) {
      debug7(printf("diag1 %d:%u advancing\n",ndiagonals1+1,(*diagonals1)));
      hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,*diagonals1++,
				     /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				     /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
				     /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				     genome,genomealt,chromosome_iit,/*plusp*/true,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3,
				     trim_ends_p,dibasep,cmetp);
    }
    
    List_free(&spanningset);
  }

  /* Minus */
  for (mod = 0; mod < 3; mod++) {
    debug(printf("Onemiss minus mod %d\n",mod));

    spanningset = this->minus_spanningset[mod];
    nelts = List_length(spanningset);
    array = (Spanningelt_T *) List_to_array(spanningset,NULL);
    /* List_free(&spanningset); */

    qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
    if (nelts > 2) {
      qsort(&(array[2]),nelts-2,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
    }
    sorted = (List_T) NULL;
    for (i = nelts-1; i >= 0; --i) {
      sorted = List_push(sorted,Spanningelt_reset(array[i]));
    }
    FREE(array);

    debug(printf("*** find_spanning_onemiss_matches, minus mod %d\n",mod));
    debug(Spanningelt_print_set(sorted));

    diagonals0 = Spanningelt_diagonals(&ndiagonals0,(Spanningelt_T) sorted->first,&miss0_querypos5,&miss0_querypos3);
    diagonals1 = Spanningelt_diagonals(&ndiagonals1,(Spanningelt_T) sorted->rest->first,&miss1_querypos5,&miss1_querypos3);
    spanningset = List_push(List_copy(sorted->rest->rest),(void **) NULL); /* Add a dummy list to front */
    nempty = 0;
    global_miss_querypos5 = querylength;
    global_miss_querypos3 = 0;
    List_free(&sorted);
    chrhigh = 0U;

    while (ndiagonals0 > 0 && ndiagonals1 > 0 && nempty <= 1) {
      if ((diagonal0 = (*diagonals0)) < (diagonal1 = (*diagonals1))) {
	debug7(printf("diag0 %d:%u advancing\n",ndiagonals0,(*diagonals0)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/false,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3,
				       trim_ends_p,dibasep,cmetp);
	++diagonals0;
	--ndiagonals0;

      } else if (diagonal1 < diagonal0) {
	debug7(printf("diag1 %d:%u advancing\n",ndiagonals1,(*diagonals1)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal1,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/false,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3,
				       trim_ends_p,dibasep,cmetp);
	++diagonals1;
	--ndiagonals1;

      } else {
	debug7(printf("diag0&1 %d:%u == %d:%u advancing\n",ndiagonals0,diagonal0,ndiagonals1,diagonal1));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/false,/*nmisses_allowed*/1,
				       /*nmisses_seen*/nempty,global_miss_querypos5,global_miss_querypos3,
				       trim_ends_p,dibasep,cmetp);
	++diagonals0;
	--ndiagonals0;
	++diagonals1;
	--ndiagonals1;
      }
    }

    while (--ndiagonals0 >= 0 && nempty == 0) {
      debug7(printf("diag0 %d:%u advancing\n",ndiagonals0+1,(*diagonals0)));
      hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,*diagonals0++,
				     /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				     /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
				     /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				     genome,genomealt,chromosome_iit,/*plusp*/false,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3,
				     trim_ends_p,dibasep,cmetp);
    }

    while (--ndiagonals1 >= 0 && nempty == 0) {
      debug7(printf("diag1 %d:%u advancing\n",ndiagonals1+1,(*diagonals1)));
      hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,*diagonals1++,
				     /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				     /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
				     /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				     genome,genomealt,chromosome_iit,/*plusp*/false,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3,
				     trim_ends_p,dibasep,cmetp);
    }

    List_free(&spanningset);
  }

  return hits;
}


static List_T
find_spanning_multimiss_matches (int *found_score, List_T hits, T this, int nrequired,
				 char *queryuc_ptr, char *queryrc, int querylength,
				 char *gbuffer, Compress_T query_compress_fwd, Compress_T query_compress_rev,
				 UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome, Genome_T genomealt,
				 IIT_T chromosome_iit, int nmisses_allowed, bool trim_ends_p, bool dibasep, bool cmetp) {
  Genomicpos_T *diagonals, diagonal;
  List_T spanningset, sorted;
  Spanningelt_T *array;
  int nunion = nmisses_allowed + nrequired, nelts;
  int heapsize, count, mod, i;
  int ndiagonals, nempty;
  int parenti, smallesti, righti;
  int global_miss_querypos5, global_miss_querypos3;
  int elt_miss_querypos5, elt_miss_querypos3;
  struct Batch_T *batchpool, sentinel_struct;
  Batch_T *heap, batch, sentinel;
  Genomicpos_T chroffset, chrhigh;
  Chrnum_T chrnum;

  debug(printf("Starting find_spanning_multimiss_matches with %d misses allowed\n",nmisses_allowed));

  sentinel_struct.diagonal = (Genomicpos_T) -1U; /* infinity */
  sentinel = &sentinel_struct;

  batchpool = (struct Batch_T *) CALLOC(nunion,sizeof(struct Batch_T));
  heap = (Batch_T *) CALLOC(2*(nunion+1)+1+1,sizeof(Batch_T)); /* being liberal with allocation */

  /* Plus */
  for (mod = 0; mod < 3; mod++) {
    debug(printf("Multimiss plus mod %d\n",mod));

    spanningset = this->plus_spanningset[mod];
    nelts = List_length(spanningset);
    array = (Spanningelt_T *) List_to_array(spanningset,NULL);
    /* List_free(&spanningset); */

    qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
    if (nelts > nunion) {
      qsort(&(array[nunion]),nelts-nunion,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
    }
    sorted = (List_T) NULL;
    for (i = nelts-1; i >= 0; --i) {
      sorted = List_push(sorted,Spanningelt_reset(array[i]));
    }
    FREE(array);

    debug(printf("*** find_spanning_multimiss_matches, %d misses allowed, plus mod %d\n",nmisses_allowed,mod));
    debug(Spanningelt_print_set(sorted));

    /* Put first few pointers into heap */
    heapsize = 0;
    spanningset = sorted;
    global_miss_querypos5 = querylength;
    global_miss_querypos3 = 0;
    for (i = 0; i < nunion && spanningset; i++, spanningset = spanningset->rest) {
      /* Get list as a special one, and perform conversion if necessary */
      diagonals = Spanningelt_diagonals(&ndiagonals,(Spanningelt_T) spanningset->first,&elt_miss_querypos5,&elt_miss_querypos3);
      if (elt_miss_querypos5 < global_miss_querypos5) global_miss_querypos5 = elt_miss_querypos5;
      if (elt_miss_querypos3 > global_miss_querypos3) global_miss_querypos3 = elt_miss_querypos3;

      batch = &(batchpool[i]);
      debug(printf("Adding batch %d of size %d...",i,ndiagonals));
      if (ndiagonals > 0) {
	Batch_init_simple(batch,diagonals,ndiagonals,querylength,/*querypos*/i);
	if (batch->npositions > 0) {
	  debug(printf("inserting into heap"));
	  min_heap_insert_simple(heap,&heapsize,batch);
	}
      }
      debug(printf("\n"));
    }
    debug(printf("heapsize is %d\n",heapsize));
    if (heapsize == 0) {
      List_free(&sorted);
    } else {
      spanningset = List_push(List_copy(spanningset),(void **) NULL); /* Add a dummy list elt to front */
      nempty = 0;
      List_free(&sorted);
      
      /* Set up rest of heap */
      for (i = heapsize+1; i <= 2*heapsize+1; i++) {
	heap[i] = sentinel;
      }

      debug7(printf("*** multimiss mod %d plus:\n",mod));

      /* Initialize loop */
      batch = heap[1];
      diagonal = batch->diagonal;
      count = 1;
      debug7(printf("at #%d, initial diagonal is %u\n",batch->querypos,diagonal));
      
      /* Update batch */
      if (--batch->npositions <= 0) {
	/* Use last entry in heap for heapify */
	batch = heap[heapsize];
	heap[heapsize--] = sentinel;
      } else {
	/* Use this batch for heapify */
	batch->diagonal = *(++batch->positions);
      }

      /* Heapify down */
      debug6(printf("Starting heapify with %u\n",diagonal));
      parenti = 1;
      smallesti = (heap[3]->diagonal < heap[2]->diagonal) ? 3 : 2;
      debug6(printf("Comparing left %d/right %d: %u and %u\n",2,3,heap[2]->diagonal,heap[3]->diagonal));
      while (batch->diagonal > heap[smallesti]->diagonal) {
	heap[parenti] = heap[smallesti];
	parenti = smallesti;
	smallesti = LEFT(parenti);
	righti = smallesti+1;
	debug6(printf("Comparing left %d/right %d: %u and %u\n",
		      smallesti,righti,heap[smallesti]->diagonal,heap[righti]->diagonal));
	if (heap[righti]->diagonal < heap[smallesti]->diagonal) {
	  smallesti = righti;
	}
      }
      heap[parenti] = batch;
      debug6(printf("Inserting at %d\n\n",parenti));

      /* Iterate through heap */
      chrhigh = 0U;
      while (heapsize > 0) {
	batch = heap[1];

	if (batch->diagonal == diagonal) {
	  count++;
	  debug7(printf("at #%d, incrementing diagonal %u to count %d\n",batch->querypos,diagonal,count));
	} else {
	  /* End of diagonal */
	  if (count >= nrequired) {
	    /* printf("Testing %d..%d\n",miss_querypos5,miss_querypos3); */
	    hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal,
					   /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
					   /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
					   /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
					   genome,genomealt,chromosome_iit,/*plusp*/true,nmisses_allowed,
					   /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3,
					   trim_ends_p,dibasep,cmetp);
	  }
	  diagonal = batch->diagonal;
	  count = 1;
	  debug7(printf("at #%d, next diagonal is %u\n",batch->querypos,diagonal));
	}

	/* Update batch */
	if (--batch->npositions <= 0) {
	  /* Use last entry in heap for heapify */
	  batch = heap[heapsize];
	  heap[heapsize--] = sentinel;
	} else {
	  /* Use this batch for heapify */
	  batch->diagonal = *(++batch->positions);
	}
      
	/* Heapify down */
	debug6(printf("Starting heapify with %u\n",diagonal));
	parenti = 1;
	smallesti = (heap[3]->diagonal < heap[2]->diagonal) ? 3 : 2;
	debug6(printf("Comparing left %d/right %d: %u and %u\n",2,3,heap[2]->diagonal,heap[3]->diagonal));
	while (batch->diagonal > heap[smallesti]->diagonal) {
	  heap[parenti] = heap[smallesti];
	  parenti = smallesti;
	  smallesti = LEFT(parenti);
	  righti = smallesti+1;
	  debug6(printf("Comparing left %d/right %d: %u and %u\n",
			smallesti,righti,heap[smallesti]->diagonal,heap[righti]->diagonal));
	  if (heap[righti]->diagonal < heap[smallesti]->diagonal) {
	    smallesti = righti;
	  }
	}
	heap[parenti] = batch;
	debug6(printf("Inserting at %d\n\n",parenti));
      }

      /* Terminate loop */
      if (count >= nrequired) {
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/true,nmisses_allowed,
				       /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3,
				       trim_ends_p,dibasep,cmetp);
      }
      List_free(&spanningset);
    }
  }

  /* Minus */
  for (mod = 0; mod < 3; mod++) {
    debug(printf("Multimiss minus mod %d\n",mod));

    spanningset = this->minus_spanningset[mod];
    nelts = List_length(spanningset);
    array = (Spanningelt_T *) List_to_array(spanningset,NULL);
    /* List_free(&spanningset); */

    qsort(array,nelts,sizeof(Spanningelt_T),Spanningelt_candidates_cmp);
    if (nelts > nunion) {
      qsort(&(array[nunion]),nelts-nunion,sizeof(Spanningelt_T),Spanningelt_pruning_cmp);
    }
    sorted = (List_T) NULL;
    for (i = nelts-1; i >= 0; --i) {
      sorted = List_push(sorted,Spanningelt_reset(array[i]));
    }
    FREE(array);

    debug(printf("*** find_spanning_multimiss_matches, %d misses_allowed, minus mod %d\n",nmisses_allowed,mod));
    debug(Spanningelt_print_set(sorted));

    /* Put first few pointers into heap */
    heapsize = 0;
    spanningset = sorted;
    global_miss_querypos5 = querylength;
    global_miss_querypos3 = 0;
    for (i = 0; i < nunion && spanningset; i++, spanningset = spanningset->rest) {
      /* Get list as a special one, and perform conversion if necessary */
      diagonals = Spanningelt_diagonals(&ndiagonals,(Spanningelt_T) spanningset->first,&elt_miss_querypos5,&elt_miss_querypos3);
      if (elt_miss_querypos5 < global_miss_querypos5) global_miss_querypos5 = elt_miss_querypos5;
      if (elt_miss_querypos3 > global_miss_querypos3) global_miss_querypos3 = elt_miss_querypos3;

      batch = &(batchpool[i]);
      debug(printf("Adding batch %d of size %d...",i,ndiagonals));
      if (ndiagonals > 0) {
	Batch_init_simple(batch,diagonals,ndiagonals,querylength,/*querypos*/i);
	if (batch->npositions > 0) {
	  debug(printf("inserting into heap"));
	  min_heap_insert_simple(heap,&heapsize,batch);
	}
      }
      debug(printf("\n"));
    }
    debug(printf("heapsize is %d\n",heapsize));

    if (heapsize == 0) {
      List_free(&sorted);
    } else {
      spanningset = List_push(List_copy(spanningset),(void **) NULL); /* Add a dummy list elt to front */
      nempty = 0;
      List_free(&sorted);

      /* Set up rest of heap */
      for (i = heapsize+1; i <= 2*heapsize+1; i++) {
	heap[i] = sentinel;
      }

      debug7(printf("*** multimiss mod %d minus:\n",mod));

      /* Initialize loop */
      batch = heap[1];
      diagonal = batch->diagonal;
      count = 1;
      debug7(printf("at #%d, initial diagonal is %u\n",batch->querypos,diagonal));

      /* Update batch */
      if (--batch->npositions <= 0) {
	/* Use last entry in heap for heapify */
	batch = heap[heapsize];
	heap[heapsize--] = sentinel;
      } else {
	/* Use this batch for heapify */
	batch->diagonal = *(++batch->positions);
      }

      /* Heapify down */
      debug6(printf("Starting heapify with %u\n",diagonal));
      parenti = 1;
      smallesti = (heap[3]->diagonal < heap[2]->diagonal) ? 3 : 2;
      debug6(printf("Comparing left %d/right %d: %u and %u\n",2,3,heap[2]->diagonal,heap[3]->diagonal));
      while (batch->diagonal > heap[smallesti]->diagonal) {
	heap[parenti] = heap[smallesti];
	parenti = smallesti;
	smallesti = LEFT(parenti);
	righti = smallesti+1;
	debug6(printf("Comparing left %d/right %d: %u and %u\n",
		      smallesti,righti,heap[smallesti]->diagonal,heap[righti]->diagonal));
	if (heap[righti]->diagonal < heap[smallesti]->diagonal) {
	  smallesti = righti;
	}
      }
      heap[parenti] = batch;
      debug6(printf("Inserting at %d\n\n",parenti));

      /* Iterate through heap */
      chrhigh = 0U;
      while (heapsize > 0) {
	batch = heap[1];

	if (batch->diagonal == diagonal) {
	  count++;
	  debug7(printf("at #%d, incrementing diagonal %u to count %d\n",batch->querypos,diagonal,count));
	} else {
	  /* End of diagonal */
	  if (count >= nrequired) {
	    hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal,
					   /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
					   /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
					   /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
					   genome,genomealt,chromosome_iit,/*plusp*/false,nmisses_allowed,
					   /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3,
					   trim_ends_p,dibasep,cmetp);
	  }
	  diagonal = batch->diagonal;
	  count = 1;
	  debug7(printf("at #%d, next diagonal is %u\n",batch->querypos,diagonal));
	}

	/* Update batch */
	if (--batch->npositions <= 0) {
	  /* Use last entry in heap for heapify */
	  batch = heap[heapsize];
	  heap[heapsize--] = sentinel;
	} else {
	  /* Use this batch for heapify */
	  batch->diagonal = *(++batch->positions);
	}

	/* Heapify down */
	debug6(printf("Starting heapify with %u\n",diagonal));
	parenti = 1;
	smallesti = (heap[3]->diagonal < heap[2]->diagonal) ? 3 : 2;
	debug6(printf("Comparing left %d/right %d: %u and %u\n",2,3,heap[2]->diagonal,heap[3]->diagonal));
	while (batch->diagonal > heap[smallesti]->diagonal) {
	  heap[parenti] = heap[smallesti];
	  parenti = smallesti;
	  smallesti = LEFT(parenti);
	  righti = smallesti+1;
	  debug6(printf("Comparing left %d/right %d: %u and %u\n",
			smallesti,righti,heap[smallesti]->diagonal,heap[righti]->diagonal));
	  if (heap[righti]->diagonal < heap[smallesti]->diagonal) {
	    smallesti = righti;
	  }
	}
	heap[parenti] = batch;
	debug6(printf("Inserting at %d\n\n",parenti));
      }

      /* Terminate loop */
      if (count >= nrequired) {
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       genome,genomealt,chromosome_iit,/*plusp*/false,nmisses_allowed,
				       /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3,
				       trim_ends_p,dibasep,cmetp);
      }
      List_free(&spanningset);
    }
  }

  FREE(heap);
  FREE(batchpool);
  return hits;
}


/************************************************************************/


#if 0
static void
trim_ends_unknowns_only (int *trim5, int *trim3, char *sequence1, char *sequence2, int length) {
  int pos;

  pos = 0;
  while (pos < length && sequence2[pos] == OUTOFBOUNDS) {
    pos++;
  }
  debug8(printf("outofbounds: trim 5': at %d: %c != %c\n",pos,sequence2[pos],OUTOFBOUNDS));
  *trim5 = pos;

  pos = length-1;
  debug8(printf("outofbounds: trim 3': %d:%c\n",pos,sequence2[pos]));
  while (pos >= 0 && sequence2[pos] == OUTOFBOUNDS) {
    pos--;
  }
  *trim3 = pos+1;
  debug8(printf("outofbounds: trim 3': %d - %d\n",length,*trim3));
  *trim3 = length - (*trim3);

  debug8(
	 printf("At query ->: %.*s\n",length,sequence1);
	 printf("At genome->: %.*s\n",length,sequence2);
	 printf("%02d %02d    ->: ",*trim5,*trim3);
	 for (pos = 0; pos < *trim5; pos++) {
	   printf(" ");
	 }
	 for ( ; pos < length - (*trim3); pos++) {
	   printf("*");
	 }
	 for ( ; pos < length; pos++) {
	   printf(" ");
	 }
	 printf("\n");
	 );

  return;
}
#endif


/************************************************************************/


/* Returns a master pointer (segments) to the block of segments */
/* If end_indel_mismatches_allowed set to 0, won't save any segments for end indels. */
static List_T
find_complete_mm_and_indelsegments (int *found_score, List_T hits, 
				    struct Indelsegment_T **segments, int *nsegments, Genomicpos_T **positions, int *npositions,
				    bool *omitted, int querylength, int query_lastpos,
				    char *query, char *queryptr, char *gbuffer, Compress_T query_compress,
				    UINT4 *genome_blocks, UINT4 *snp_blocks,
				    Genome_T genome, IIT_T chromosome_iit, Floors_T floors,
				    int max_mismatches_allowed, int indel_mismatches_allowed, int max_indel_sep, int max_end_insertions,
				    int firstbound, int lastbound, bool trim_ends_p, bool dibasep, bool cmetp, bool plusp) {
  Batch_T batch, sentinel;
  struct Batch_T sentinel_struct, *batchpool;
  Batch_T *heap;
  int heapsize = 0;
  int parenti, smallesti, righti, i;
  int querypos, first_querypos, last_querypos, prev_first_querypos, prev_last_querypos;
  int floor, floor_xfirst, floor_xlast, prev_floor_xfirst, prev_floor_xlast,
    floor_incr;
  int *floors_from_neg3, *floors_from_xfirst, *floors_to_xlast, *floors_to_pos3;
  /* int exclude_xfirst, exclude_xlast; */
  Genomicpos_T diagonal, last_diagonal, prev_diagonal, chroffset = 0U, chrhigh = 0U,
    prev_chroffset, prev_chrhigh;
  Chrnum_T chrnum, prev_chrnum;
#ifdef HAVE_64_BIT
  UINT8 diagonal_add_querypos;
#endif
  int nmismatches, ncolordiffs, total_npositions = 0;
  Indelsegment_T ptr;
  Genomicpos_T left;

  Genomicpos_T middle_indel_bound = 0U;	/* Prevents first segment from being saved */
  bool prev_savedp;		/* Not necessary to initialize, because of above */

  debug(printf("*** Starting find_complete_mm_and_indelsegments with max_mismatches_allowed %d, firstbound %d, lastbound %d ***\n",
	       max_mismatches_allowed,firstbound,lastbound));

  /* Create sentinel */
#ifdef HAVE_64_BIT
  sentinel_struct.diagonal_add_querypos = (UINT8) -1U; /* infinity */
  sentinel_struct.diagonal_add_querypos <<= 32;
#else
  sentinel_struct.querypos = querylength; /* essentially infinity */
  sentinel_struct.diagonal = (Genomicpos_T) -1U; /* infinity */
#endif
  sentinel = &sentinel_struct;

  /* Set up batches */
  batchpool = (struct Batch_T *) CALLOC(query_lastpos+1,sizeof(struct Batch_T));
  heap = (Batch_T *) CALLOC(2*(query_lastpos+1)+1+1,sizeof(Batch_T));

  /* Don't add entries for compoundpos positions (skip querypos -2, -1, lastpos+1, lastpos+2) */
  if (plusp) {
    for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
      if (omitted[querypos] == true) {
	debug0(printf("Not adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
      } else if (npositions[querypos] > 0) {
	debug0(printf("Adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
	batch = &(batchpool[i]);
	Batch_init(batch,querypos,/*diagterm*/querylength - querypos,positions[querypos],npositions[querypos],querylength);
	total_npositions += npositions[querypos];
	if (batch->npositions > 0) {
	  min_heap_insert(heap,&heapsize,batch);
	  i++;
	}
      } else {
	debug0(printf("Not adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
      }
    }
  } else {
    for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
      if (omitted[querypos] == true) {
	debug0(printf("Not adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
      } else if (npositions[querypos] > 0) {
	debug0(printf("Adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
	batch = &(batchpool[i]);
	Batch_init(batch,querypos,/*diagterm*/querypos + INDEX1PART,positions[querypos],npositions[querypos],querylength);
	total_npositions += npositions[querypos];
	if (batch->npositions > 0) {
	  min_heap_insert(heap,&heapsize,batch);
	  i++;
	}
      } else {
	debug0(printf("Not adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
      }
    }
  }


  if (i == 0) {
    FREE(heap);
    FREE(batchpool);
    *nsegments = 0;
    return hits;
  }

  /* Set up rest of heap */
  for (i = heapsize+1; i <= 2*heapsize+1; i++) {
    heap[i] = sentinel;
  }

  *segments = (struct Indelsegment_T *) CALLOC(total_npositions,sizeof(struct Indelsegment_T));
  ptr = &((*segments)[0]);

  /*
  if ((exclude_xfirst = firstbound-2-INDEX1PART-max_end_insertions) < 3) {
    exclude_xfirst = 3;
  }
  if ((exclude_xlast = lastbound+1+max_end_insertions) > query_lastpos-3) {
    exclude_xlast = query_lastpos-3;
  }
  */

#if 0
  /* Should account for firstbound and lastbound */
  floors_from_xfirst = floors->scorefrom[/* xfirst_from = */ firstbound-3+max_end_insertions];
  floors_to_xlast = floors->scoreto[/* xlast_to = */ lastbound+4-INDEX1PART-max_end_insertions];
#else
  if (INDEX1PART+max_end_insertions > query_lastpos + 3) {
    floors_from_xfirst = floors->scorefrom[query_lastpos+3];
  } else {
    floors_from_xfirst = floors->scorefrom[INDEX1PART+max_end_insertions];
  }
  if (query_lastpos-INDEX1PART-max_end_insertions < -3) {
    floors_to_xlast = floors->scoreto[-3];
  } else {
    floors_to_xlast = floors->scoreto[query_lastpos-INDEX1PART-max_end_insertions];
  }
#endif
  floors_from_neg3 = floors->scorefrom[-3];
  floors_to_pos3 = floors->scoreto[query_lastpos+3];

  /* Initialize loop */
  batch = heap[1];
  first_querypos = last_querypos = querypos = batch->querypos;
  last_diagonal = diagonal = batch->diagonal;

  floor = floor_xlast = floors_from_neg3[first_querypos];
#ifdef PREBOUNDS
  floor_xfirst = (first_querypos < 3) ? querylength : floors->scorefrom[INDEX1PART+max_end_insertions][first_querypos];
  debug(printf("first_xfirst, checking if first_querypos < 3 and taking score from %d+%d\n",INDEX1PART,max_end_insertions));
  debug(printf("first_xlast, checking if last_querypos > %d and taking score to %d-%d\n",query_lastpos-3,query_lastpos-INDEX1PART,max_end_insertions));
#else
  /* floor_xfirst = (first_querypos < exclude_xfirst) ? querylength : floors->score[firstbound-3+max_end_insertions][first_querypos]; */

  floor_xfirst = floors_from_xfirst[first_querypos] /* floors->scorefrom[xfirst_from][first_querypos] */;
  debug(printf("first_xfirst, checking if first_querypos < %d and taking score from %d+%d\n",firstbound-INDEX1PART,firstbound-3,max_end_insertions));
  debug(printf("first_xlast, checking if last_querypos > %d and taking score to %d-%d\n",lastbound+1,lastbound+3+1,max_end_insertions));
#endif
  debug0(printf("multiple_mm_%s, diagonal %u, querypos %d\n",plusp ? "plus" : "minus",diagonal,querypos));
  debug0(printf("first_querypos = %d => initialize floor %d, floor_xfirst %d, floor_xlast %d\n",
	        first_querypos,floor,floor_xfirst,floor_xlast));

  if (--batch->npositions <= 0) {
    /* Use last entry in heap for insertion */
    batch = heap[heapsize];
    querypos = batch->querypos;
    heap[heapsize--] = sentinel;

  } else {
    /* Use this batch for insertion (same querypos) */
#ifdef WORDS_BIGENDIAN
    batch->diagonal = Bigendian_convert_uint(*(++batch->positions)) + batch->diagterm;
#else
    batch->diagonal = *(++batch->positions) + batch->diagterm;
#endif
#ifdef HAVE_64_BIT
    batch->diagonal_add_querypos = (UINT8) batch->diagonal;
    batch->diagonal_add_querypos <<= 32;
    batch->diagonal_add_querypos |= querypos /* Previously added 2 because querypos was -2: + 2*/;
#endif
  }

  /* heapify */
  parenti = 1;
#ifdef HAVE_64_BIT
  diagonal_add_querypos = batch->diagonal_add_querypos;
  smallesti = (heap[3]->diagonal_add_querypos < heap[2]->diagonal_add_querypos) ? 3 : 2;
  while (diagonal_add_querypos > heap[smallesti]->diagonal_add_querypos) {
    heap[parenti] = heap[smallesti];
    parenti = smallesti;
    smallesti = LEFT(parenti);
    righti = smallesti+1;
    if (heap[righti]->diagonal_add_querypos < heap[smallesti]->diagonal_add_querypos) {
      smallesti = righti;
    }
  }
#else
  diagonal = batch->diagonal;
  smallesti = ((heap[3]->diagonal < heap[2]->diagonal) ||
	       ((heap[3]->diagonal == heap[2]->diagonal) &&
		(heap[3]->querypos < heap[2]->querypos))) ? 3 : 2;
  /* Note that diagonal/querypos will never exceed a sentinel diagonal/querypos */
  while (diagonal > heap[smallesti]->diagonal ||
	 (diagonal == heap[smallesti]->diagonal &&
	  querypos > heap[smallesti]->querypos)) {
    heap[parenti] = heap[smallesti];
    parenti = smallesti;
    smallesti = LEFT(parenti);
    righti = smallesti+1;
    if ((heap[righti]->diagonal < heap[smallesti]->diagonal) ||
		  ((heap[righti]->diagonal == heap[smallesti]->diagonal) &&
		   (heap[righti]->querypos < heap[smallesti]->querypos))) {
      smallesti = righti;
    }
  }
#endif
  heap[parenti] = batch;


  /* Continue after initialization */ 
  while (heapsize > 0) {
    batch = heap[1];
    querypos = batch->querypos;
    diagonal = batch->diagonal;

    if (diagonal == last_diagonal) {
      /* Continuing exact match or substitution */
      floor_incr = floors->scorefrom[last_querypos][querypos];
      floor += floor_incr;
      floor_xfirst += floor_incr;
      floor_xlast += floor_incr;
      debug0(printf("diagonal %u unchanged: last_querypos = %d, querypos = %d => floor increments by %d\n",
		    diagonal,last_querypos,querypos,floor_incr));
      debug0(printf("*multiple_mm_%s, diagonal %u, querypos %d, floor %d, floor_xfirst %d, floor_xlast %d\n",
		    plusp ? "plus" : "minus",diagonal,querypos,floor,floor_xfirst,floor_xlast));
    } else {
      /* End of diagonal */
      floor += floors_to_pos3[last_querypos]  /* floors->score[last_querypos][query_lastpos+3] */;
      floor_xfirst += floors_to_pos3[last_querypos]  /* floors->score[last_querypos][query_lastpos+3] */;
#ifdef PREBOUNDS
      floor_xlast += (last_querypos > query_lastpos-3) ? querylength : floors->scorefrom[last_querypos][query_lastpos-INDEX1PART-max_end_insertions];
#else
      /* floor_xlast += (last_querypos > exclude_xlast) ? querylength : floors->score[last_querypos][lastbound+3+1-INDEX1PART-max_end_insertions]; */
      floor_xlast += floors_to_xlast[last_querypos];  /* floors->score[last_querypos][xlast_to]; */
#endif
      debug0(printf("new diagonal %u > last diagonal %u: last_querypos = %d => terminate floor %d, floor_xfirst %d, floor_xlast %d\n",
		    diagonal,last_diagonal,last_querypos,floor,floor_xfirst,floor_xlast));

      if (last_diagonal > chrhigh) {
	/* update chromosome bounds, based on low end */
	chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,last_diagonal-querylength,last_diagonal-querylength);
	IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
	chrhigh += 1U;
      }
      if (last_diagonal <= chrhigh) { /* FORMULA for high position */
	/* position of high end is within current chromosome */
	debug0(printf("  => multiple_mm, floor %d (vs max_mismatches_allowed %d): diagonal %u, query %d..%d, chrbounds %u..%u, floor_xfirst %d, floor_xlast %d\n",
		      floor,max_mismatches_allowed,last_diagonal,first_querypos,last_querypos,chroffset,chrhigh,floor_xfirst,floor_xlast));
	if (floor <= max_mismatches_allowed) {
	  left = last_diagonal - querylength;
	  nmismatches = Genome_count_mismatches_limit(&ncolordiffs,queryptr,query_compress,genome_blocks,snp_blocks,
						      left,/*pos5*/0,/*pos3*/querylength,
						      max_mismatches_allowed,dibasep,cmetp,plusp);
	  debug0(printf("nmismatches = %d (vs report at %d)\n",nmismatches,max_mismatches_allowed));
	  if (nmismatches <= max_mismatches_allowed) {
	    Genome_fill_buffer_simple(genome,left,querylength,gbuffer);
	    hits = List_push(hits,(void *) Stage3_new_substitution(&(*found_score),nmismatches,ncolordiffs,
								   left,/*genomiclength*/querylength,
								   plusp,gbuffer,query,chrnum,chroffset,
								   trim_ends_p,dibasep,cmetp));
	  }
	}

	if (last_diagonal < middle_indel_bound) {
	  /* Save previous and this one */
	  if (prev_savedp == false) {
	    ptr->diagonal = prev_diagonal;
	    ptr->chrnum = prev_chrnum;
	    ptr->chroffset = prev_chroffset;
	    ptr->chrhigh = prev_chrhigh;
	    ptr->querypos5 = prev_first_querypos;
	    ptr->querypos3 = prev_last_querypos;
	    ptr->floor_xfirst = prev_floor_xfirst;
	    ptr->floor_xlast = prev_floor_xlast;
	    ptr++;
	  }
	  ptr->diagonal = last_diagonal;
	  ptr->chrnum = chrnum;
	  ptr->chroffset = chroffset;
	  ptr->chrhigh = chrhigh;
	  ptr->querypos5 = first_querypos;
	  ptr->querypos3 = last_querypos;
	  ptr->floor_xfirst = floor_xfirst;
	  ptr->floor_xlast = floor_xlast;
	  ptr++;
	  debug0(printf("Save previous diagonal %u and this diagonal %u < %u for middle indels\n",
			prev_diagonal,last_diagonal,middle_indel_bound));
	  prev_savedp = true;
	} else if (floor_xfirst /* hack until xfirst fixed: - 0 */ <= indel_mismatches_allowed || 
		   floor_xlast /* hack until xlast fixed: - 0 */ <= indel_mismatches_allowed) {
	  /* Save this one */
	  ptr->diagonal = last_diagonal;
	  ptr->chrnum = chrnum;
	  ptr->chroffset = chroffset;
	  ptr->chrhigh = chrhigh;
	  ptr->querypos5 = first_querypos;
	  ptr->querypos3 = last_querypos;
	  ptr->floor_xfirst = floor_xfirst;
	  ptr->floor_xlast = floor_xlast;
	  ptr++;
	  debug0(printf("Save this diagonal for %u for end indels (%u >= middle_indel_bound %u)\n",
			last_diagonal,last_diagonal,middle_indel_bound));
	  prev_savedp = true;
	} else {
	  debug0(printf("Don't save this diagonal %u (%u >= middle_indel_bound %u)\n",
			last_diagonal,last_diagonal,middle_indel_bound));
	  prev_savedp = false;
	}
	middle_indel_bound = last_diagonal + max_indel_sep;
	debug0(printf("Set middle_indel_bound to be %u + %u = %u\n",last_diagonal,max_indel_sep,middle_indel_bound));
	prev_diagonal = last_diagonal;
	prev_chrnum = chrnum;
	prev_chroffset = chroffset;
	prev_chrhigh = chrhigh;
	prev_first_querypos = first_querypos;
	prev_last_querypos = last_querypos;
	prev_floor_xfirst = floor_xfirst;
	prev_floor_xlast = floor_xlast;
      }
      
      /* Prepare next diagonal */
      first_querypos = querypos;
      last_diagonal = diagonal;
      floor = floor_xlast = floors_from_neg3[first_querypos] /* floors->score[-3][first_querypos] */;
#ifdef PREBOUNDS
      floor_xfirst = (first_querypos < 3) ? querylength : floors->scorefrom[INDEX1PART+max_end_insertions][first_querypos];
#else
      /* floor_xfirst = (first_querypos < exclude_xfirst) ? querylength : floors->score[firstbound-3+max_end_insertions][first_querypos]; */
      floor_xfirst = floors_from_xfirst[first_querypos];  /* floors->score[xfirst_from][first_querypos]; */
#endif
      debug0(printf("*multiple_mm_%s, diagonal %u, querypos %d, floor %d, floor_xfirst %d, floor_xlast %d\n",
		    plusp ? "plus" : "minus",diagonal,querypos,floor,floor_xfirst,floor_xlast));
      debug0(printf("start of diagonal %u, first_querypos = %d => initialize floor %d, floor_xfirst %d, floor_xlast %d\n",
		    diagonal,first_querypos,floor,floor_xfirst,floor_xlast));
    }
    last_querypos = querypos;


    if (--batch->npositions <= 0) {
      /* Use last entry in heap for insertion */
      batch = heap[heapsize];
      querypos = batch->querypos;
      heap[heapsize--] = sentinel;

    } else {
      /* Use this batch for insertion (same querypos) */
#ifdef WORDS_BIGENDIAN
      batch->diagonal = Bigendian_convert_uint(*(++batch->positions)) + batch->diagterm;
#else
      batch->diagonal = *(++batch->positions) + batch->diagterm;
#endif
#ifdef HAVE_64_BIT
      batch->diagonal_add_querypos = (UINT8) batch->diagonal;
      batch->diagonal_add_querypos <<= 32;
      batch->diagonal_add_querypos |= querypos /* Previously added 2 because querypos was -2: + 2*/;
#endif
    }

    /* heapify */
    parenti = 1;
#ifdef HAVE_64_BIT
    diagonal_add_querypos = batch->diagonal_add_querypos;
    smallesti = (heap[3]->diagonal_add_querypos < heap[2]->diagonal_add_querypos) ? 3 : 2;
    while (diagonal_add_querypos > heap[smallesti]->diagonal_add_querypos) {
      heap[parenti] = heap[smallesti];
      parenti = smallesti;
      smallesti = LEFT(parenti);
      righti = smallesti+1;
      if (heap[righti]->diagonal_add_querypos < heap[smallesti]->diagonal_add_querypos) {
	smallesti = righti;
      }
    }
#else
    diagonal = batch->diagonal;
    smallesti = ((heap[3]->diagonal < heap[2]->diagonal) ||
		 ((heap[3]->diagonal == heap[2]->diagonal) &&
		  (heap[3]->querypos < heap[2]->querypos))) ? 3 : 2;
    /* Note that diagonal/querypos will never exceed a sentinel diagonal/querypos */
    while (diagonal > heap[smallesti]->diagonal ||
	   (diagonal == heap[smallesti]->diagonal &&
	    querypos > heap[smallesti]->querypos)) {
      heap[parenti] = heap[smallesti];
      parenti = smallesti;
      smallesti = LEFT(parenti);
      righti = smallesti+1;
      if ((heap[righti]->diagonal < heap[smallesti]->diagonal) ||
	  ((heap[righti]->diagonal == heap[smallesti]->diagonal) &&
	   (heap[righti]->querypos < heap[smallesti]->querypos))) {
	smallesti = righti;
      }
    }
#endif
    heap[parenti] = batch;
  }

  /* Terminate loop. */
  floor += floors_to_pos3[last_querypos];   /* floors->score[last_querypos][query_lastpos+3]; */
  floor_xfirst += floors_to_pos3[last_querypos];  /* floors->score[last_querypos][query_lastpos+3]; */
#ifdef PREBOUNDS
  floor_xlast += (last_querypos > query_lastpos-3) ? querylength : floors->score[last_querypos][query_lastpos-INDEX1PART-max_end_insertions];
#else
  /* floor_xlast += (last_querypos > exclude_xlast) ? querylength : floors->score[last_querypos][lastbound+3+1-INDEX1PART-max_end_insertions]; */
  floor_xlast += floors_to_xlast[last_querypos];  /* floors->score[last_querypos][xlast_to]; */
#endif
  debug0(printf("no more diagonals: last_querypos = %d => terminate floor %d, floor_xfirst %d, floor_xlast %d\n",
		last_querypos,floor,floor_xfirst,floor_xlast));

  if (last_diagonal > chrhigh) {
    /* update chromosome bounds, based on low end */
    chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,last_diagonal-querylength,last_diagonal-querylength);
    IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
    chrhigh += 1U;
  }
  if (last_diagonal <= chrhigh) { /* FORMULA for high position */
    /* position of high end is within current chromosome */
    debug0(printf("  => multiple_mm, floor %d (vs max_mismatches_allowed %d): diagonal %u, query %d..%d, chrbounds %u..%u, floor_xfirst %d, floor_xlast %d\n",
		  floor,max_mismatches_allowed,last_diagonal,first_querypos,last_querypos,chroffset,chrhigh,floor_xfirst,floor_xlast));
    if (floor <= max_mismatches_allowed) {
      left = last_diagonal - querylength;
      nmismatches = Genome_count_mismatches_limit(&ncolordiffs,queryptr,query_compress,genome_blocks,snp_blocks,
						  left,/*pos5*/0,/*pos3*/querylength,
						  max_mismatches_allowed,dibasep,cmetp,plusp);
      debug0(printf("nmismatches = %d (vs report at %d)\n",nmismatches,max_mismatches_allowed));
      if (nmismatches <= max_mismatches_allowed) {
	Genome_fill_buffer_simple(genome,left,querylength,gbuffer);
	hits = List_push(hits,(void *) Stage3_new_substitution(&(*found_score),nmismatches,ncolordiffs,
							       left,/*genomiclength*/querylength,
							       plusp,gbuffer,query,chrnum,chroffset,
							       trim_ends_p,dibasep,cmetp));
      }
    }

    if (last_diagonal < middle_indel_bound) {
      if (prev_savedp == false) {
	ptr->diagonal = prev_diagonal;
	ptr->chrnum = prev_chrnum;
	ptr->chroffset = prev_chroffset;
	ptr->chrhigh = prev_chrhigh;
	ptr->querypos5 = prev_first_querypos;
	ptr->querypos3 = prev_last_querypos;
	ptr->floor_xfirst = prev_floor_xfirst;
	ptr->floor_xlast = prev_floor_xlast;
	ptr++;
      }
      ptr->diagonal = last_diagonal;
      ptr->chrnum = chrnum;
      ptr->chroffset = chroffset;
      ptr->chrhigh = chrhigh;
      ptr->querypos5 = first_querypos;
      ptr->querypos3 = last_querypos;
      ptr->floor_xfirst = floor_xfirst;
      ptr->floor_xlast = floor_xlast;
      ptr++;
      debug0(printf("Save previous diagonal %u and this diagonal %u\n",prev_diagonal,last_diagonal));
    } else if (floor_xfirst /* hack until xfirst fixed: - 0 */ <= indel_mismatches_allowed || 
	       floor_xlast /* hack until xlast fixed: - 0 */ <= indel_mismatches_allowed) {
      /* Save this one */
      ptr->diagonal = last_diagonal;
      ptr->chrnum = chrnum;
      ptr->chroffset = chroffset;
      ptr->chrhigh = chrhigh;
      ptr->querypos5 = first_querypos;
      ptr->querypos3 = last_querypos;
      ptr->floor_xfirst = floor_xfirst;
      ptr->floor_xlast = floor_xlast;
      ptr++;
      debug0(printf("Save previous diagonal %u and this diagonal %u < %u for middle indels\n",
		    prev_diagonal,last_diagonal,middle_indel_bound));
    } else {
      debug0(printf("Don't save this diagonal %u (%u >= middle_indel_bound %u)\n",
		    last_diagonal,last_diagonal,middle_indel_bound));
    }
  }

  FREE(heap);
  FREE(batchpool);

  /* Note: *segments is in descending diagonal order.  Will need to
     reverse before solving middle deletions */

  *nsegments = ptr - (*segments);
  debug2(printf("nsegments = %d\n",*nsegments));
  debug(printf("nsegments = %d (total_npositions = %d)\n",*nsegments,total_npositions));

  return hits;
}


/* Needed for novelsplicing */
static struct Splicesegment_T *
identify_segments_all (int *nsegments, Genomicpos_T **positions, int *npositions,
		       Storedoligomer_T *oligos, bool *omitted, int querylength, int query_lastpos,
		       IIT_T chromosome_iit, bool plusp) {
  struct Splicesegment_T *segments = NULL;
  Splicesegment_T ptr;
  int total_npositions = 0;
  Batch_T batch, sentinel;
  struct Batch_T sentinel_struct, *batchpool;
  Batch_T *heap;
  int heapsize = 0;
  int i;
  int parenti, smallesti, righti;
  int querypos;
  int first_querypos, last_querypos;
  Genomicpos_T diagonal, last_diagonal, chroffset = 0U, chrhigh = 0U;
  Chrnum_T chrnum;
#ifdef HAVE_64_BIT
  UINT8 diagonal_add_querypos;
#endif

  debug(printf("*** Starting identify_segments_all ***\n"));

  /* Create sentinel */
#ifdef HAVE_64_BIT
  sentinel_struct.diagonal_add_querypos = (UINT8) -1U; /* infinity */
  sentinel_struct.diagonal_add_querypos <<= 32;
#else
  sentinel_struct.querypos = querylength; /* essentially infinity */
  sentinel_struct.diagonal = (Genomicpos_T) -1U; /* infinity */
#endif
  sentinel = &sentinel_struct;

  /* Set up batches */
  batchpool = (struct Batch_T *) CALLOC(query_lastpos+1,sizeof(struct Batch_T));
  heap = (Batch_T *) CALLOC(2*(query_lastpos+1)+1+1,sizeof(Batch_T));

  /* Don't add entries for compoundpos positions (skip querypos -2, -1, lastpos+1, lastpos+2) */
  if (plusp) {
    for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
      if (omitted[querypos] == true) {
	debug0(printf("Not adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
      } else if (npositions[querypos] > 0) {
	debug0(printf("Adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
	batch = &(batchpool[i]);
	Batch_init(batch,querypos,/*diagterm*/querylength - querypos,positions[querypos],npositions[querypos],querylength);
	total_npositions += npositions[querypos];
	if (batch->npositions > 0) {
	  min_heap_insert(heap,&heapsize,batch);
	  i++;
	}
      } else {
	debug0(printf("Not adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
      }
    }
  } else {
    for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
      if (omitted[querypos] == true) {
	debug0(printf("Not adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
      } else if (npositions[querypos] > 0) {
	debug0(printf("Adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
	batch = &(batchpool[i]);
	Batch_init(batch,querypos,/*diagterm*/querypos + INDEX1PART,positions[querypos],npositions[querypos],querylength);
	total_npositions += npositions[querypos];
	if (batch->npositions > 0) {
	  min_heap_insert(heap,&heapsize,batch);
	  i++;
	}
      } else {
	debug0(printf("Not adding batch for querypos %d with %d positions, omitted %d, and oligo %06X\n",
		      querypos,npositions[querypos],omitted[querypos],oligos[querypos]));
      }
    }
  }


  if (i == 0) {
    FREE(heap);
    FREE(batchpool);
    *nsegments = 0;
    return (struct Splicesegment_T *) NULL;
  }

  /* Set up rest of heap */
  for (i = heapsize+1; i <= 2*heapsize+1; i++) {
    heap[i] = sentinel;
  }


  segments = (struct Splicesegment_T *) CALLOC(total_npositions,sizeof(struct Splicesegment_T));
  ptr = &(segments[0]);

  /* Initialize loop */
  batch = heap[1];
  first_querypos = last_querypos = querypos = batch->querypos;
  last_diagonal = diagonal = batch->diagonal;
  debug0(printf("multiple_mm_%s, diagonal %u, querypos %d\n",plusp ? "plus" : "minus",diagonal,querypos));

  if (--batch->npositions <= 0) {
    /* Use last entry in heap for insertion */
    batch = heap[heapsize];
    querypos = batch->querypos;
    heap[heapsize--] = sentinel;

  } else {
    /* Use this batch for insertion (same querypos) */
#ifdef WORDS_BIGENDIAN
    batch->diagonal = Bigendian_convert_uint(*(++batch->positions)) + batch->diagterm;
#else
    batch->diagonal = *(++batch->positions) + batch->diagterm;
#endif
#ifdef HAVE_64_BIT
    batch->diagonal_add_querypos = (UINT8) batch->diagonal;
    batch->diagonal_add_querypos <<= 32;
    batch->diagonal_add_querypos |= querypos /* Previously added 2 because querypos was -2: + 2*/;
#endif
  }

  /* heapify */
  parenti = 1;
#ifdef HAVE_64_BIT
  diagonal_add_querypos = batch->diagonal_add_querypos;
  smallesti = (heap[3]->diagonal_add_querypos < heap[2]->diagonal_add_querypos) ? 3 : 2;
  while (diagonal_add_querypos > heap[smallesti]->diagonal_add_querypos) {
    heap[parenti] = heap[smallesti];
    parenti = smallesti;
    smallesti = LEFT(parenti);
    righti = smallesti+1;
    if (heap[righti]->diagonal_add_querypos < heap[smallesti]->diagonal_add_querypos) {
      smallesti = righti;
    }
  }
#else
  diagonal = batch->diagonal;
  smallesti = ((heap[3]->diagonal < heap[2]->diagonal) ||
	       ((heap[3]->diagonal == heap[2]->diagonal) &&
		(heap[3]->querypos < heap[2]->querypos))) ? 3 : 2;
  /* Note that diagonal/querypos will never exceed a sentinel diagonal/querypos */
  while (diagonal > heap[smallesti]->diagonal ||
	 (diagonal == heap[smallesti]->diagonal &&
	  querypos > heap[smallesti]->querypos)) {
    heap[parenti] = heap[smallesti];
    parenti = smallesti;
    smallesti = LEFT(parenti);
    righti = smallesti+1;
    if ((heap[righti]->diagonal < heap[smallesti]->diagonal) ||
		  ((heap[righti]->diagonal == heap[smallesti]->diagonal) &&
		   (heap[righti]->querypos < heap[smallesti]->querypos))) {
      smallesti = righti;
    }
  }
#endif
  heap[parenti] = batch;

  /* Continue after initialization */ 
  while (heapsize > 0) {
    batch = heap[1];
    querypos = batch->querypos;
    diagonal = batch->diagonal;

    if (diagonal == last_diagonal) {
      /* Continuing exact match or substitution */
      debug0(printf("diagonal %u unchanged: last_querypos = %d, querypos = %d\n",
		    diagonal,last_querypos,querypos));
      debug0(printf("*multiple_mm_%s, diagonal %u, querypos %d\n",plusp ? "plus" : "minus",diagonal,querypos));
    } else {
      /* End of diagonal */
      debug0(printf("new diagonal %u > last diagonal %u\n",diagonal,last_diagonal));

      if (last_diagonal > chrhigh) {
	/* update chromosome bounds, based on low end */
	chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,last_diagonal-querylength,last_diagonal-querylength);
	IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
	chrhigh += 1U;
      }
      if (last_diagonal <= chrhigh) { /* FORMULA for high position */
	/* position of high end is within current chromosome */
	debug0(printf("  => identify_segments_all: diagonal %u, querypos %d..%d, chrbounds %u..%u\n",
		      last_diagonal,first_querypos,last_querypos,chroffset,chrhigh));
	/* Save all segments */
	ptr->diagonal = last_diagonal;
	ptr->querypos5 = first_querypos;
	ptr->querypos3 = last_querypos;
	ptr->chrnum = chrnum;
	ptr->chroffset = chroffset;
	ptr++;
      }

      /* Prepare next diagonal */
      first_querypos = querypos;
      last_diagonal = diagonal;
    }
    last_querypos = querypos;


    if (--batch->npositions <= 0) {
      /* Use last entry in heap for insertion */
      batch = heap[heapsize];
      querypos = batch->querypos;
      heap[heapsize--] = sentinel;

    } else {
      /* Use this batch for insertion (same querypos) */
#ifdef WORDS_BIGENDIAN
      batch->diagonal = Bigendian_convert_uint(*(++batch->positions)) + batch->diagterm;
#else
      batch->diagonal = *(++batch->positions) + batch->diagterm;
#endif
#ifdef HAVE_64_BIT
      batch->diagonal_add_querypos = (UINT8) batch->diagonal;
      batch->diagonal_add_querypos <<= 32;
      batch->diagonal_add_querypos |= querypos /* Previously added 2 because querypos was -2: + 2*/;
#endif
    }

    /* heapify */
    parenti = 1;
#ifdef HAVE_64_BIT
    diagonal_add_querypos = batch->diagonal_add_querypos;
    smallesti = (heap[3]->diagonal_add_querypos < heap[2]->diagonal_add_querypos) ? 3 : 2;
    while (diagonal_add_querypos > heap[smallesti]->diagonal_add_querypos) {
      heap[parenti] = heap[smallesti];
      parenti = smallesti;
      smallesti = LEFT(parenti);
      righti = smallesti+1;
      if (heap[righti]->diagonal_add_querypos < heap[smallesti]->diagonal_add_querypos) {
	smallesti = righti;
      }
    }
#else
    diagonal = batch->diagonal;
    smallesti = ((heap[3]->diagonal < heap[2]->diagonal) ||
		 ((heap[3]->diagonal == heap[2]->diagonal) &&
		  (heap[3]->querypos < heap[2]->querypos))) ? 3 : 2;
    /* Note that diagonal/querypos will never exceed a sentinel diagonal/querypos */
    while (diagonal > heap[smallesti]->diagonal ||
	   (diagonal == heap[smallesti]->diagonal &&
	    querypos > heap[smallesti]->querypos)) {
      heap[parenti] = heap[smallesti];
      parenti = smallesti;
      smallesti = LEFT(parenti);
      righti = smallesti+1;
      if ((heap[righti]->diagonal < heap[smallesti]->diagonal) ||
	  ((heap[righti]->diagonal == heap[smallesti]->diagonal) &&
	   (heap[righti]->querypos < heap[smallesti]->querypos))) {
	smallesti = righti;
      }
    }
#endif
    heap[parenti] = batch;
  }

  /* Terminate loop. */
  debug0(printf("no more diagonals: last_querypos = %d\n",last_querypos));

  if (last_diagonal > chrhigh) {
    /* update chromosome bounds, based on low end */
    chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,last_diagonal-querylength,last_diagonal-querylength);
    IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
    chrhigh += 1U;
  }
  if (last_diagonal <= chrhigh) { /* FORMULA for high position */
    /* position of high end is within current chromosome */
    debug0(printf("  => identify_segments_all: diagonal %u, querypos %d..%d, chrbounds %u..%u\n",
		  last_diagonal,first_querypos,last_querypos,chroffset,chrhigh));
    /* Save all segments */
    ptr->diagonal = last_diagonal;
    ptr->querypos5 = first_querypos;
    ptr->querypos3 = last_querypos;
    ptr->chrnum = chrnum;
    ptr->chroffset = chroffset;
    ptr++;
  }

  FREE(heap);
  FREE(batchpool);

  /* Note: *segments is in descending diagonal order.  Will need to
     reverse before find novel splice pairs */

  *nsegments = ptr - segments;
  debug(printf("nsegments = %d (total_npositions = %d)\n",*nsegments,total_npositions));

  return segments;
}


/* indels is positive here */
static List_T
solve_middle_insertion (bool *foundp, int *found_score, List_T hits,
			Indelsegment_T ptr, int indels,
			Compress_T query_compress, Genome_T genome,
			UINT4 *genome_blocks, UINT4 *snp_blocks,
			char *query, char *queryptr, int querylength,
			char *gbuffer, int min_indel_end_matches, int indel_penalty, int max_mismatches_allowed,
			bool trim_ends_p, bool dibasep, bool cmetp, bool plusp) {
#ifdef DEBUG2
  int i;
#endif
  Genomicpos_T left;
  int indel_pos;
  int mismatch_positions_left[MAX_QUERYLENGTH], mismatch_positions_right[MAX_QUERYLENGTH];
  int colordiffs_left[MAX_QUERYLENGTH], colordiffs_right[MAX_QUERYLENGTH];
  int nmismatches_left, nmismatches_right;
  int sum, lefti, righti;

  /* query has insertion.  Get |indels| less from genome; trim from left. */
  left = ptr->diagonal - querylength;

  debug2(Genome_fill_buffer_simple(genome,left+indels,querylength-indels,gbuffer));
  debug2(printf("solve_middle_indel, plus, insertion: Getting genome at diagonal %u - querylength %d + indels %d = %u\n",
		ptr->diagonal,querylength,indels,left+indels));
  debug2(printf("g1: %s\n",gbuffer));
  debug2(printf("q:  %s\n",queryptr));
  debug2(printf("g2: %s\n",&(gbuffer[indels])));
  
    /* No need to check chromosome bounds */
  nmismatches_left = Genome_mismatches_left(mismatch_positions_left,colordiffs_left,max_mismatches_allowed,
					    queryptr,query_compress,genome_blocks,snp_blocks,
					    left+indels,/*pos5*/0,/*pos3*/querylength,dibasep,cmetp,plusp);
  debug2(
	 printf("%d mismatches on left at:",nmismatches_left);
	 for (i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");

	 if (dibasep) {
	   printf("colordiffs on left:");
	   for (i = 0; i < nmismatches_left; i++) {
	     printf(" %d",colordiffs_left[i]);
	   }
	   printf("\n");
	 });


  /* No need to check chromosome bounds */
  nmismatches_right = Genome_mismatches_right(mismatch_positions_right,colordiffs_right,max_mismatches_allowed,
					      queryptr,query_compress,genome_blocks,snp_blocks,
					      left,/*pos5*/0,/*pos3*/querylength,dibasep,cmetp,plusp);
  debug2(
	 printf("%d mismatches on right at:",nmismatches_right);
	 for (i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
  
	 if (dibasep) {
	   printf("colordiffs on right:");
	   for (i = 0; i < nmismatches_right; i++) {
	     printf(" %d",colordiffs_right[i]);
	   }
	   printf("\n");
	 });


  for (sum = 0; sum <= max_mismatches_allowed; sum++) {
    for (lefti = 0; lefti <= sum && lefti < nmismatches_left; lefti++) {
      if ((righti = sum - lefti) < nmismatches_right && 
	  mismatch_positions_left[lefti] + indels > mismatch_positions_right[righti]) {
	indel_pos = mismatch_positions_left[lefti]; /* This appears to be correct */
	/* indel_pos = mismatch_positions_right[righti]; */
	/* printf("indel_pos using left: %d\n",mismatch_positions_left[lefti]); */
	/* printf("indel_pos using right: %d\n",mismatch_positions_right[righti]); */
	debug2(printf("indel_pos: %d\n",indel_pos));
	if (indel_pos >= min_indel_end_matches && indel_pos + indels <= querylength - min_indel_end_matches) {
	  *foundp = true;
	  Genome_fill_buffer_simple(genome,left+indels,querylength-indels,gbuffer);
	  debug2(printf("successful insertion with %d mismatches and indel_pos at %d\n",sum,indel_pos));
	  if (plusp == true) {
	    return List_push(hits,
			     (void *) Stage3_new_insertion(&(*found_score),indels,indel_pos,
							   /*nmismatches1*/lefti,/*nmismatches2*/righti,
							   /*ncolordiffs1: colordiffs_left[lefti]*/0,
							   /*ncolordiffs2: colordiffs_right[righti]*/0,
							   /*left*/left+indels,/*genomiclength*/querylength-indels,
							   querylength,/*plusp*/true,gbuffer,query,ptr->chrnum,
							   ptr->chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	  } else {
	    return List_push(hits,
			     (void *) Stage3_new_insertion(&(*found_score),indels,/*indel_pos*/querylength-indel_pos-indels,
							   /*nmismatches1*/righti,/*nmismatches2*/lefti,
							   /*ncolordiffs1: colordiffs_right[righti]*/0,
							   /*ncolordiffs2: colordiffs_left[lefti]*/0,
							   /*left*/left+indels,/*genomiclength*/querylength-indels,
							   querylength,/*plusp*/false,gbuffer,query,ptr->chrnum,
							   ptr->chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	  }
	}
      }
    }
  }

  *foundp = false;
  return hits;
}


/* indels is negative here */
static List_T
solve_middle_deletion (bool *foundp, int *found_score, List_T hits,
		       Indelsegment_T ptr, int indels,
		       Compress_T query_compress, Genome_T genome,
		       UINT4 *genome_blocks, UINT4 *snp_blocks,
		       char *query, char *queryptr, int querylength,
		       char *gbuffer, int min_indel_end_matches, int indel_penalty, int max_mismatches_allowed,
		       bool trim_ends_p, bool dibasep, bool cmetp, bool plusp) {
#ifdef DEBUG2
  int i;
#endif
  Genomicpos_T left;
  int indel_pos;
  int mismatch_positions_left[MAX_QUERYLENGTH], mismatch_positions_right[MAX_QUERYLENGTH];
  int colordiffs_left[MAX_QUERYLENGTH], colordiffs_right[MAX_QUERYLENGTH];
  int nmismatches_left, nmismatches_right;
  int sum, lefti, righti;

  /* query has deletion.  Get |indels| more from genome; add to right. */
  left = ptr->diagonal - querylength;
  debug2(Genome_fill_buffer_simple(genome,left,querylength-indels,gbuffer));
  debug2(printf("solve_middle_indel, plus, deletion: Getting genome at diagonal %u - querylength %d = %u\n",
		ptr->diagonal,querylength,left));
  debug2(printf("g1: %s\n",gbuffer));
  debug2(printf("q:  %s\n",queryptr));
  debug2(printf("g2: %s\n",&(gbuffer[-indels])));
  
  /* No need to check chromosome bounds */
  nmismatches_left = Genome_mismatches_left(mismatch_positions_left,colordiffs_left,max_mismatches_allowed,
					    queryptr,query_compress,genome_blocks,snp_blocks,
					    left,/*pos5*/0,/*pos3*/querylength,dibasep,cmetp,plusp);
  
  debug2(
	 printf("%d mismatches on left at:",nmismatches_left);
	 for (i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");

	 if (dibasep) {
	   printf("colordiffs on left:");
	   for (i = 0; i < nmismatches_left; i++) {
	     printf(" %d",colordiffs_left[i]);
	   }
	   printf("\n");
	 });
  
  /* No need to check chromosome bounds */
  nmismatches_right = Genome_mismatches_right(mismatch_positions_right,colordiffs_right,max_mismatches_allowed,
					      queryptr,query_compress,genome_blocks,snp_blocks,
					      left-indels,/*pos5*/0,/*pos3*/querylength,dibasep,cmetp,plusp);
  
  debug2(
	 printf("%d mismatches on right at:",nmismatches_right);
	 for (i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");

	 if (dibasep) {
	   printf("colordiffs on right:");
	   for (i = 0; i < nmismatches_right; i++) {
	     printf(" %d",colordiffs_right[i]);
	   }
	   printf("\n");
	 });
  
  for (sum = 0; sum <= max_mismatches_allowed; sum++) {
    for (lefti = 0; lefti <= sum && lefti < nmismatches_left; lefti++) {
      if ((righti = sum - lefti) < nmismatches_right &&
	  mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	/* indel_pos = mismatch_positions_left[lefti]; */
	indel_pos = mismatch_positions_right[righti] + 1; /* This is leftmost position in righti+1 .. lefti */
	/* printf("indel_pos using left: %d\n",mismatch_positions_left[lefti]); */
	/* printf("indel_pos using right: %d\n",mismatch_positions_right[righti]+1); */
	debug2(printf("indel_pos: %d\n",indel_pos));
	if (indel_pos >= min_indel_end_matches && indel_pos <= querylength - min_indel_end_matches) {
	  *foundp = true;
	  Genome_fill_buffer_simple(genome,left,querylength-indels,gbuffer);
	  debug2(printf("successful deletion with %d mismatches and indel_pos at %d\n",
			sum,indel_pos));
	  if (plusp == true) {
	    return List_push(hits,
			     (void *) Stage3_new_deletion(&(*found_score),-indels,indel_pos,
							  /*nmismatches1*/lefti,/*nmismatches2*/righti,
							  /*ncolordiffs1: colordiffs_left[lefti]*/0,
							  /*ncolordiffs2: colordiffs_right[righti]*/0,
							  left,/*genomiclength*/querylength-indels,
							  querylength,/*plusp*/true,gbuffer,query,ptr->chrnum,
							  ptr->chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	  } else {
	    return List_push(hits,
			     (void *) Stage3_new_deletion(&(*found_score),-indels,/*indel_pos*/querylength-indel_pos,
							  /*nmismatches1*/righti,/*nmismatches2*/lefti,
							  /*ncolordiffs1: colordiffs_right[righti]*/0,
							  /*ncolordiffs2: colordiffs_left[lefti]*/0,
							  left,/*genomiclength*/querylength-indels,
							  querylength,/*plusp*/false,gbuffer,query,ptr->chrnum,
							  ptr->chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	  }
	}
      }
    }
  }
  
  *foundp = false;
  return hits;
}


/* 

The pattern below is a middle insertion on plus strand, or middle deletion on minus strand:

diagonal 2354, querypos 18
diagonal 2354, querypos 19
diagonal 2354, querypos 20
diagonal 2354, querypos 21
diagonal 2354, querypos 22
diagonal 2354, querypos 23
diagonal 2354, querypos 24
diagonal 2356, querypos 0
diagonal 2356, querypos 1
diagonal 2356, querypos 2
diagonal 2356, querypos 3
diagonal 2356, querypos 4
diagonal 2356, querypos 5


The pattern below is a middle deletion on plus strand, or middle insertion on minus strand:

diagonal 2354, querypos 0
diagonal 2354, querypos 1
diagonal 2354, querypos 2
diagonal 2354, querypos 3
diagonal 2354, querypos 4
diagonal 2354, querypos 5
diagonal 2354, querypos 6
diagonal 2354, querypos 7
diagonal 2356, querypos 18
diagonal 2356, querypos 19
diagonal 2356, querypos 20
diagonal 2356, querypos 21
diagonal 2356, querypos 22
diagonal 2356, querypos 23
diagonal 2356, querypos 24

*/


static List_T
find_middle_indels (int *found_score, List_T hits, 
		    struct Indelsegment_T *plus_segments, struct Indelsegment_T *minus_segments,
		    int plus_nsegments, int minus_nsegments,
		    Genome_T genome, UINT4 *genome_blocks, UINT4 *snp_blocks,
		    char *queryuc_ptr, char *queryrc,
		    Floors_T floors, int querylength, int query_lastpos,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    int max_middle_insertions, int max_middle_deletions, int min_indel_end_matches,
		    int indel_penalty, int max_mismatches_allowed, bool trim_ends_p, bool dibasep, bool cmetp) {
  int indels, floor, pos, prev, middle;
  int *floors_from_neg3, *floors_to_pos3;
  int nfound = 0;
  Indelsegment_T segmenti, segmentj;
  char *gbuffer;
  bool foundp;

  debug(printf("*** find_middle_indels with max_mismatches_allowed %d ***\n",
	       max_mismatches_allowed));

  gbuffer = (char *) CALLOC(MAX_QUERYLENGTH + max_middle_deletions + 1,sizeof(char));

  floors_from_neg3 = floors->scorefrom[-3];
  floors_to_pos3 = floors->scoreto[query_lastpos+3];

  if (plus_nsegments > 1) {
    for (segmenti = plus_segments; segmenti < &(plus_segments[plus_nsegments-1]); segmenti++) {
      for (segmentj = segmenti+1; segmentj < &(plus_segments[plus_nsegments]) && 
	     segmentj->diagonal <= segmenti->diagonal + max_middle_insertions &&
	     segmentj->chrnum == segmenti->chrnum; segmentj++) {

	debug2(printf("plus insertion?  diagonal %u, querypos %d..%d => diagonal %u, querypos %d..%d => ",
		      segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
		      segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	/* j5 j3 i5 i3 */
	if (segmentj->querypos3 < segmenti->querypos5) {
	  indels = segmentj->diagonal - segmenti->diagonal; /* positive */
	  floor = floors_from_neg3[segmentj->querypos5] + floors_to_pos3[segmenti->querypos3]
	    /* floors->score[-3][segmentj->querypos5] + floors->score[segmenti->querypos3][query_lastpos+3] */ ;
	  if (floors->prev_omitted == NULL) {
	    if ((middle = FLOOR_MIDDLE(segmenti->querypos5 - segmentj->querypos3 - indels)) > 0) {
	      middle--;	/* for insertion, which looks like a mismatch */
	    }
	    debug2(printf("\nmiddle (no omission): %d\n",middle));
	    floor += middle;	
	  } else {
	    pos = segmenti->querypos5;
	    debug2(printf("\nmiddle (omission):"));
	    while (pos > segmentj->querypos3) {
	      if ((prev = floors->prev_omitted[pos]) < segmentj->querypos3) {
		prev = segmentj->querypos3;
	      }
	      if ((middle = FLOOR_MIDDLE(pos - prev - indels)) > 0) {
		middle--;	/* for insertion, which looks like a mismatch */
	      }
	      floor += middle;
	      debug2(printf("(%d..%d)+%d,",prev,pos,middle));
	      pos = prev;
	    }
	    debug2(printf("\n"));
	  }
	  if (floor <= max_mismatches_allowed) {
	    debug2(printf("insertion, floor = %d+middle+%d=%d, indels = %d\n",
			  floors->scorefrom[-3][segmentj->querypos5],
			  floors->scorefrom[segmenti->querypos3][query_lastpos+3],
			  floor,indels));
	    hits = solve_middle_insertion(&foundp,&(*found_score),hits,segmenti,indels,
					  /*query_compress*/query_compress_fwd,genome,genome_blocks,snp_blocks,
					  /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
					  min_indel_end_matches,indel_penalty,max_mismatches_allowed,
					  trim_ends_p,dibasep,cmetp,/*plusp*/true);
	    if (foundp == true) {
	      nfound++;
	    }
	  } else {
	    debug2(printf("too many mismatches, because floor %d+middle+%d=%d > %d\n",
			  floors->scorefrom[-3][segmentj->querypos5],
			  floors->scorefrom[segmenti->querypos3][query_lastpos+3],
			  floor,max_mismatches_allowed));
	  }
	} else {
	  debug2(printf("garbage, because querypos3 %d >= querypos5 %d\n",
			segmentj->querypos3,segmenti->querypos5));
	}
      }

      for (segmentj = segmenti+1; segmentj < &(plus_segments[plus_nsegments]) && 
	     segmentj->diagonal <= segmenti->diagonal + max_middle_deletions &&
	     segmentj->chrnum == segmenti->chrnum; segmentj++) {
	debug2(printf("plus deletion?  diagonal %u, querypos %d..%d => diagonal %u, querypos %d..%d => ",
		      segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
		      segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	/* i5 i3 j5 j3 */
	if (segmenti->querypos3 < segmentj->querypos5) {
	  indels = segmenti->diagonal - segmentj->diagonal; /* negative */
	  floor = floors_from_neg3[segmenti->querypos5] + floors_to_pos3[segmentj->querypos3]
	    /* floors->score[-3][segmenti->querypos5] + floors->score[segmentj->querypos3][query_lastpos+3] */;
	  if (floors->prev_omitted == NULL) {
	    if ((middle = FLOOR_MIDDLE(segmentj->querypos5 - segmenti->querypos3 /*- indels*/)) > 0) {
	      middle--;	/* for deletion, which looks like a mismatch */
	    }
	    debug2(printf("\nmiddle (no omission): %d\n",middle));
	    floor += middle;	
	  } else {
	    pos = segmentj->querypos5;
	    debug2(printf("\nmiddle (omission):"));
	    while (pos > segmenti->querypos3) {
	      if ((prev = floors->prev_omitted[pos]) < segmenti->querypos3) {
		prev = segmenti->querypos3;
	      }
	      if ((middle = FLOOR_MIDDLE(pos - prev /*- indels*/)) > 0) {
		middle--;	/* for deletion, which looks like a mismatch */
	      }
	      floor += middle;
	      debug2(printf("(%d..%d)+%d,",prev,pos,middle));
	      pos = prev;
	    }
	    debug2(printf("\n"));
	  }
	  if (floor <= max_mismatches_allowed) {
	    debug2(printf("deletion, floor = %d+middle+%d=%d, indels = %d\n",
			  floors->scorefrom[-3][segmenti->querypos5],
			  floors->scorefrom[segmentj->querypos3][query_lastpos+3],
			  floor,indels));
	    hits = solve_middle_deletion(&foundp,&(*found_score),hits,segmenti,indels,
					 /*query_compress*/query_compress_fwd,genome,genome_blocks,snp_blocks,
					 /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,gbuffer,
					 min_indel_end_matches,indel_penalty,max_mismatches_allowed,
					 trim_ends_p,dibasep,cmetp,/*plusp*/true);
	    if (foundp == true) {
	      nfound++;
	    }
	  } else {
	    debug2(printf("too many mismatches, because floor = %d+middle+%d=%d > %d\n",
			  floors->scorefrom[-3][segmenti->querypos5],
			  floors->scorefrom[segmentj->querypos3][query_lastpos+3],
			  floor,max_mismatches_allowed));
	  }
	} else {
	  debug2(printf("garbage, because querypos3 %d >= querypos5 %d\n",
			segmenti->querypos3,segmentj->querypos5));
	}
      }
    }
  }

  if (minus_nsegments > 1) {
    for (segmenti = minus_segments; segmenti < &(minus_segments[minus_nsegments-1]); segmenti++) {
      for (segmentj = segmenti+1; segmentj < &(minus_segments[minus_nsegments]) &&
	     segmentj->diagonal <= segmenti->diagonal + max_middle_deletions &&
	     segmentj->chrnum == segmenti->chrnum; segmentj++) {
	debug2(printf("minus deletion?  diagonal %u, querypos %d..%d => diagonal %u, querypos %d..%d => ",
		      segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
		      segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	/* j5 j3 i5 i3 */
	if (segmentj->querypos3 < segmenti->querypos5) {
	  indels = segmenti->diagonal - segmentj->diagonal; /* negative */
	  floor = floors_from_neg3[segmentj->querypos5] + floors_to_pos3[segmenti->querypos3]
	    /* floors->score[-3][segmentj->querypos5] + floors->score[segmenti->querypos3][query_lastpos+3] */;
	  if (floors->prev_omitted == NULL) {
	    if ((middle = FLOOR_MIDDLE(segmenti->querypos5 - segmentj->querypos3 /*- indels*/)) > 0) {
	      middle--;	/* for deletion, which looks like a mismatch */
	    }
	    debug2(printf("\nmiddle (no omission): %d\n",middle));
	    floor += middle;	
	  } else {
	    pos = segmenti->querypos5;
	    debug2(printf("\nmiddle (omission):"));
	    while (pos > segmentj->querypos3) {
	      if ((prev = floors->prev_omitted[pos]) < segmentj->querypos3) {
		prev = segmentj->querypos3;
	      }
	      if ((middle = FLOOR_MIDDLE(pos - prev /*- indels*/)) > 0) {
		middle--; /* for deletion, which looks like a mismatch */
	      }
	      floor += middle;
	      debug2(printf("(%d..%d)+%d,",prev,pos,middle));
	      pos = prev;
	    }
	    debug2(printf("\n"));
	  }
	  if (floor <= max_mismatches_allowed) {
	    debug2(printf("deletion, floor = %d+middle+%d=%d, indels = %d\n",
			  floors->scorefrom[-3][segmentj->querypos5],
			  floors->scorefrom[segmenti->querypos3][query_lastpos+3],
			  floor,indels));
	    hits = solve_middle_deletion(&foundp,&(*found_score),hits,segmenti,indels,
					 /*query_compress*/query_compress_rev,genome,genome_blocks,snp_blocks,
					 /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
					 min_indel_end_matches,indel_penalty,max_mismatches_allowed,
					 trim_ends_p,dibasep,cmetp,/*plusp*/false);
	    if (foundp == true) {
	      nfound++;
	    }
	  } else {
	    debug2(printf("too many mismatches, because floor = %d+middle+%d=%d > %d\n",
			  floors->scorefrom[-3][segmentj->querypos5],
			  floors->scorefrom[segmenti->querypos3][query_lastpos+3],
			  floor,max_mismatches_allowed));
	    debug2(printf("too many mismatches, because floor %d > %d\n",floor,max_mismatches_allowed));
	  }
	} else {
	  debug2(printf("garbage, because querypos3 %d >= querypos5 %d\n",
			segmentj->querypos3,segmenti->querypos5));
	}
      }

      for (segmentj = segmenti+1; segmentj < &(minus_segments[minus_nsegments]) &&
	     segmentj->diagonal <= segmenti->diagonal + max_middle_insertions &&
	     segmentj->chrnum == segmenti->chrnum; segmentj++) {
	debug2(printf("minus insertion?  diagonal %u, querypos %d..%d => diagonal %u, querypos %d..%d => ",
		      segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
		      segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	/* i5 i3 j5 j3 */
	if (segmenti->querypos3 < segmentj->querypos5) {
	  indels = segmentj->diagonal - segmenti->diagonal; /* positive */
	  floor = floors_from_neg3[segmenti->querypos5] + floors_to_pos3[segmentj->querypos3]
	    /* floors->score[-3][segmenti->querypos5] + floors->score[segmentj->querypos3][query_lastpos+3] */;
	  if (floors->prev_omitted == NULL) {
	    if ((middle = FLOOR_MIDDLE(segmentj->querypos5 - segmenti->querypos3 - indels)) > 0) {
	      middle--;	/* for insertion, which looks like a mismatch */
	    }
	    debug2(printf("\nmiddle (no omission): %d\n",middle));
	    floor += middle;	
	  } else {
	    pos = segmentj->querypos5;
	    debug2(printf("\nmiddle (omission):"));
	    while (pos > segmenti->querypos3) {
	      if ((prev = floors->prev_omitted[pos]) < segmenti->querypos3) {
		prev = segmenti->querypos3;
	      }
	      if ((middle = FLOOR_MIDDLE(pos - prev - indels)) > 0) {
		middle--;	/* for insertion, which looks like a mismatch */
	      }
	      floor += middle;
	      debug2(printf("(%d..%d)+%d,",prev,pos,middle));
	      pos = prev;
	    }
	    debug2(printf("\n"));
	  }
	  if (floor <= max_mismatches_allowed) {
	    debug2(printf("insertion, floor = %d+middle+%d=%d, indels = %d\n",
			  floors->scorefrom[-3][segmenti->querypos5],
			  floors->scorefrom[segmentj->querypos3][query_lastpos+3],
			  floor,indels));
	    hits = solve_middle_insertion(&foundp,&(*found_score),hits,segmenti,indels,
					  /*query_compress*/query_compress_rev,genome,genome_blocks,snp_blocks,
					  /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,gbuffer,
					  min_indel_end_matches,indel_penalty,max_mismatches_allowed,
					  trim_ends_p,dibasep,cmetp,/*plusp*/false);
	    if (foundp == true) {
	      nfound++;
	    }
	  } else {
	    debug2(printf("too many mismatches, because floor %d+middle+%d=%d > %d\n",
			  floors->scorefrom[-3][segmenti->querypos5],
			  floors->scorefrom[segmentj->querypos3][query_lastpos+3],
			  floor,max_mismatches_allowed));
	  }
	} else {
	  debug2(printf("garbage, because querypos3 %d >= querypos5 %d\n",
			segmenti->querypos3,segmentj->querypos5));
	}
      }
    }
  }

  FREE(gbuffer);

  return hits;
}


/************************************************************************/

static bool
compute_end_indels_right (int *indels, int *ncolordiffs, int *indel_pos, int length1, int querylength, 
			  Genomicpos_T left, char *queryuc_ptr, Compress_T query_compress,
			  UINT4 *genome_blocks, UINT4 *snp_blocks, int min_indel_end_matches,
			  int max_end_insertions, int max_end_deletions,
			  bool dibasep, bool cmetp, bool plusp) {
#ifdef DEBUG2E
  int i;
#endif
  int sep;
  int nmismatches_right;
  int mismatch_positions_right[MAX_QUERYLENGTH+1], colordiffs_right[MAX_QUERYLENGTH];

  debug2e(printf("Entered compute_end_indels_right with length1 = %d\n",length1));
  if (max_end_deletions > length1 - min_indel_end_matches) {
    max_end_deletions = length1 - min_indel_end_matches;
  }
  if (max_end_insertions > length1 - min_indel_end_matches) {
    max_end_insertions = length1 - min_indel_end_matches;
  }

  for (sep = 1; sep <= max_end_deletions && sep <= max_end_insertions; sep++) {
    *indels = -sep;
    debug2e(printf("Trying deletion of %d\n",*indels));
    if ((nmismatches_right = Genome_mismatches_right(mismatch_positions_right,colordiffs_right,/*max_mismatches_allowed*/0,
						     queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						     left-(*indels),/*pos5*/0,/*pos3*/querylength,
						     dibasep,cmetp,plusp)) > 0) {
      debug2e(
	      printf("%d mismatches on right at:",nmismatches_right);
	      for (i = 0; i <= nmismatches_right; i++) {
		printf(" %d",mismatch_positions_right[i]);
	      }
	      printf("\n");

	      if (dibasep) {
		printf("colordiffs on right:");
		for (i = 0; i < nmismatches_right; i++) {
		  printf(" %d",colordiffs_right[i]);
		}
		printf("\n");
	      });

      if (querylength - mismatch_positions_right[0] > length1 /* && querylength - mismatch_positions_right[0] >= querylength */) {
	if ((*indel_pos = mismatch_positions_right[0] + 1) > 0 && *indel_pos < querylength - 1) {
	  /* *nmismatches = 0; */
	  /* *ncolordiffs = colordiffs_right[0]; */
	  *ncolordiffs = 0;
	  return true;
	}
      }
    }

    *indels = +sep;
    debug2e(printf("Trying insertion of %d\n",*indels));
    if ((nmismatches_right = Genome_mismatches_right(mismatch_positions_right,colordiffs_right,/*max_mismatches_allowed*/0,
						     queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						     left-(*indels),/*pos5*/0,/*pos3*/querylength,
						     dibasep,cmetp,plusp)) > 0) {
      debug2e(
	      printf("%d mismatches on right at:",nmismatches_right);
	      for (i = 0; i <= nmismatches_right; i++) {
		printf(" %d",mismatch_positions_right[i]);
	      }
	      printf("\n");
	      
	      if (dibasep) {
		printf("colordiffs on right:");
		for (i = 0; i < nmismatches_right; i++) {
		  printf(" %d",colordiffs_right[i]);
		}
		printf("\n");
	      });

      /* Second clause derived from: querylength - mismatch_positions_right[0] + (*indels) >= querylength.  Not clear if necessary. */
      if (querylength - mismatch_positions_right[0] + (*indels) > length1 && (*indels) >= mismatch_positions_right[0]) {
	if ((*indel_pos = mismatch_positions_right[0] + 1 - (*indels)) > 0 && *indel_pos < querylength - 1) {
	  /* *nmismatches = 0; */
	  /* *ncolordiffs = colordiffs_right[0]; */
	  *ncolordiffs = 0;
	  return true;
	}
      }
    }
  }

  for ( ; sep <= max_end_deletions; sep++) {
    *indels = -sep;
    debug2e(printf("Trying deletion of %d\n",*indels));
    if ((nmismatches_right = Genome_mismatches_right(mismatch_positions_right,colordiffs_right,/*max_mismatches_allowed*/0,
						     queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						     left-(*indels),/*pos5*/0,/*pos3*/querylength,
						     dibasep,cmetp,plusp)) > 0) {
      debug2e(
	      printf("%d mismatches on right at:",nmismatches_right);
	      for (i = 0; i <= nmismatches_right; i++) {
		printf(" %d",mismatch_positions_right[i]);
	      }
	      printf("\n");

	      if (dibasep) {
		printf("colordiffs on right:");
		for (i = 0; i < nmismatches_right; i++) {
		  printf(" %d",colordiffs_right[i]);
		}
		printf("\n");
	      });

      if (querylength - mismatch_positions_right[0] > length1 /* && querylength - mismatch_positions_right[0] >= querylength */) {
	if ((*indel_pos = mismatch_positions_right[0] + 1) > 0 && *indel_pos < querylength - 1) {
	  /* *nmismatches = 0; */
	  /* *ncolordiffs = colordiffs_right[0]; */
	  *ncolordiffs = 0;
	  return true;
	}
      }
    }
  }

  for ( ; sep <= max_end_insertions; sep++) {
    *indels = +sep;
    debug2e(printf("Trying insertion of %d\n",*indels));
    if ((nmismatches_right = Genome_mismatches_right(mismatch_positions_right,colordiffs_right,/*max_mismatches_allowed*/0,
						     queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						     left-(*indels),/*pos5*/0,/*pos3*/querylength,
						     dibasep,cmetp,plusp)) > 0) {

      debug2e(
	      printf("%d mismatches on right at:",nmismatches_right);
	      for (i = 0; i <= nmismatches_right; i++) {
		printf(" %d",mismatch_positions_right[i]);
	      }
	      printf("\n");

	      if (dibasep) {
		printf("colordiffs on right:");
		for (i = 0; i < nmismatches_right; i++) {
		  printf(" %d",colordiffs_right[i]);
		}
		printf("\n");
	      });

      /* Second clause derived from: querylength - mismatch_positions_right[0] + (*indels) >= querylength.  Not clear if necessary. */
      if (querylength - mismatch_positions_right[0] + (*indels) > length1 && (*indels) >= mismatch_positions_right[0]) {
	if ((*indel_pos = mismatch_positions_right[0] + 1 - (*indels)) > 0 && *indel_pos < querylength - 1) {
	  /* *nmismatches = 0; */
	  /* *ncolordiffs = colordiffs_right[0]; */
	  *ncolordiffs = 0;
	  return true;
	}
      }
    }
  }

  return false;
}


static bool
compute_end_indels_left (int *indels, int *ncolordiffs, int *indel_pos, int length1, int querylength,
			 Genomicpos_T left, char *queryuc_ptr, Compress_T query_compress,
			 UINT4 *genome_blocks, UINT4 *snp_blocks, int min_indel_end_matches,
			 int max_end_insertions, int max_end_deletions,
			 bool dibasep, bool cmetp, bool plusp) {
#ifdef DEBUG2E
  int i;
#endif
  int sep;
  int nmismatches_left;
  int mismatch_positions_left[MAX_QUERYLENGTH+1];
  int colordiffs_left[MAX_QUERYLENGTH];

  debug2e(printf("Entered compute_end_indels_left with length1 = %d\n",length1));
  if (max_end_deletions > length1 - min_indel_end_matches) {
    max_end_deletions = length1 - min_indel_end_matches;
  }
  if (max_end_insertions > length1 - min_indel_end_matches) {
    max_end_insertions = length1 - min_indel_end_matches;
  }

  for (sep = 1; sep <= max_end_deletions && sep <= max_end_insertions; sep++) {
    if (left >= sep) {
      *indels = -sep;
      debug2e(printf("Trying deletion of %d\n",*indels));
      if ((nmismatches_left = Genome_mismatches_left(mismatch_positions_left,colordiffs_left,/*max_mismatches_allowed*/0,
						     queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						     left+(*indels),/*pos5*/0,/*pos3*/querylength,dibasep,cmetp,plusp)) > 0) {
	debug2e(
		printf("%d mismatches on left at:",nmismatches_left);
		for (i = 0; i <= nmismatches_left; i++) {
		  printf(" %d",mismatch_positions_left[i]);
		}
		printf("\n");

		if (dibasep) {
		  printf("colordiffs on left:");
		  for (i = 0; i < nmismatches_left; i++) {
		    printf(" %d",colordiffs_left[i]);
		  }
		  printf("\n");
		});

	if (mismatch_positions_left[0] > length1 /* && mismatch_positions_left[0] <= querylength */) {
	  if ((*indel_pos = mismatch_positions_left[0]) > 0 && *indel_pos < querylength - 1) {
	    /* *nmismatches = 0; */
	    /* *ncolordiffs = colordiffs_left[0]; */
	    *ncolordiffs = 0;
	    return true;
	  }
	}
      }
    }
      
    *indels = +sep;
    debug2e(printf("Trying insertion of %d\n",*indels));
    if ((nmismatches_left = Genome_mismatches_left(mismatch_positions_left,colordiffs_left,/*max_mismatches_allowed*/0,
						   queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						   left+(*indels),/*pos5*/0,/*pos3*/querylength,
						   dibasep,cmetp,plusp)) > 0) {
      debug2e(
	      printf("%d mismatches on left at:",nmismatches_left);
	      for (i = 0; i <= nmismatches_left; i++) {
		printf(" %d",mismatch_positions_left[i]);
	      }
	      printf("\n");

	      if (dibasep) {
		printf("colordiffs on left:");
		for (i = 0; i < nmismatches_left; i++) {
		  printf(" %d",colordiffs_left[i]);
		}
		printf("\n");
	      });

      if (mismatch_positions_left[0] + (*indels) > length1 && mismatch_positions_left[0] + (*indels) <= querylength) {
	if ((*indel_pos = mismatch_positions_left[0] + (*indels)) > 0 && *indel_pos < querylength - 1) {
	  /* *nmismatches = 0; */
	  /* *ncolordiffs = colordiffs_left[0]; */
	  *ncolordiffs = 0;
	  return true;
	}
      }
    }
  }

  for ( ; sep <= max_end_deletions; sep++) {
    if (left >= sep) {
      *indels = -sep;
      debug2e(printf("Trying deletion of %d\n",*indels));
      if ((nmismatches_left = Genome_mismatches_left(mismatch_positions_left,colordiffs_left,/*max_mismatches_allowed*/0,
						     queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						     left+(*indels),/*pos5*/0,/*pos3*/querylength,dibasep,cmetp,plusp)) > 0) {
	debug2e(
		printf("%d mismatches on left at:",nmismatches_left);
		for (i = 0; i <= nmismatches_left; i++) {
		  printf(" %d",mismatch_positions_left[i]);
		}
		printf("\n");

		if (dibasep) {
		  printf("colordiffs on left:");
		  for (i = 0; i < nmismatches_left; i++) {
		    printf(" %d",colordiffs_left[i]);
		  }
		  printf("\n");
		});
      
	if (mismatch_positions_left[0] > length1 /* && mismatch_positions_left <= querylength */) {
	  if ((*indel_pos = mismatch_positions_left[0]) > 0 && *indel_pos < querylength - 1) {
	    /* *nmismatches = 0; */
	    /* *ncolordiffs = colordiffs_left[0]; */
	    *ncolordiffs = 0;
	    return true;
	  }
	}
      }
    }
  }
    
  for ( ; sep <= max_end_insertions; sep++) {
    *indels = +sep;
    debug2e(printf("Trying insertion of %d\n",*indels));
    if ((nmismatches_left = Genome_mismatches_left(mismatch_positions_left,colordiffs_left,/*max_mismatches_allowed*/0,
						   queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						   left+(*indels),/*pos5*/0,/*pos3*/querylength,
						   dibasep,cmetp,plusp)) > 0) {
      debug2e(
	      printf("%d mismatches on left at:",nmismatches_left);
	      for (i = 0; i <= nmismatches_left; i++) {
		printf(" %d",mismatch_positions_left[i]);
	      }
	      printf("\n");

	      if (dibasep) {
		printf("colordiffs on left:");
		for (i = 0; i < nmismatches_left; i++) {
		  printf(" %d",colordiffs_left[i]);
		}
		printf("\n");
	      });

      if (mismatch_positions_left[0] + (*indels) > length1 && mismatch_positions_left[0] + (*indels) <= querylength) {
	if ((*indel_pos = mismatch_positions_left[0] + (*indels)) > 0 && *indel_pos < querylength - 1) {
	  /* *nmismatches = 0; */
	  /* *ncolordiffs = colordiffs_left[0]; */
	  *ncolordiffs = 0;
	  return true;
	}
      }
    }
  }

  return false;
}


/************************************************************************/

/* solve 5' end: indel in first oligomer */
static List_T
solve_first_indel_plus (int *found_score, List_T hits, Genomicpos_T diagonal, int firstbound,
			Chrnum_T chrnum, Genomicpos_T chroffset,
			Genome_T genome, UINT4 *genome_blocks, UINT4 *snp_blocks,
			char *query, char *queryptr, int querylength, Compress_T query_compress_fwd,
			char *gbuffer, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
			int indel_penalty, int max_mismatches, bool trim_ends_p, bool dibasep, bool cmetp) {
  int i;
  Genomicpos_T left;
  int indels, indel_pos, endgoal, nmismatches_long, ncolordiffs_long, ncolordiffs_short;
  int mismatch_positions[MAX_QUERYLENGTH];
  int colordiffs[MAX_QUERYLENGTH];
  int nmismatches;
  
  left = diagonal - querylength;
  if (max_end_deletions > left - chroffset) {
    max_end_deletions = left - chroffset;
    /* diagonal - querylength guaranteed to be >= chroffset, so max_end_deletions >= 0 */
  }

  debug2e(printf("\nsolve_first_indel, plus: Getting genome at diagonal %u - querylength %d - max_end_deletions %d = %u\n",
		diagonal,querylength,max_end_deletions,left-max_end_deletions));
  debug2e(Genome_fill_buffer_simple(genome,left-max_end_deletions,querylength+max_end_deletions,gbuffer));
  debug2e(printf("q: %s\ng: %s\n",queryptr,gbuffer));

  /* No need to check chromosome bounds */
  nmismatches_long = Genome_mismatches_right(mismatch_positions,colordiffs,max_mismatches,
					     queryptr,query_compress_fwd,genome_blocks,snp_blocks,
					     left,/*pos5*/firstbound,/*pos3*/querylength,dibasep,cmetp,/*plusp*/true);

  debug2e(
	 printf("first_indel_plus: %d mismatches to right of firstbound %d:",nmismatches_long,firstbound);
	 for (i = 0; i <= nmismatches_long; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");

	 if (dibasep) {
	   printf("colordiffs on right:");
	   for (i = 0; i < nmismatches_long; i++) {
	     printf(" %d",colordiffs[i]);
	   }
	   printf("\n");
	 });

  if (nmismatches_long <= max_mismatches) {
    endgoal = firstbound;
  } else {
    endgoal = mismatch_positions[max_mismatches];
  }
  if (endgoal <= firstbound /* && endgoal > min_indel_end_matches */) {
    if (compute_end_indels_left(&indels,&ncolordiffs_short,&indel_pos,/*length1*/endgoal,querylength,
				left,queryptr,/*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				min_indel_end_matches,max_end_insertions,max_end_deletions,
#if 0
				/*max_mismatches_allowed*/max_mismatches-nmismatches_long,
#endif
				dibasep,cmetp,/*plusp*/true) == true) {
      debug2e(printf("Got indel_pos %d.\n",indel_pos));

      Genome_fill_buffer_simple(genome,left-max_end_deletions,querylength+max_end_deletions,gbuffer);
      if (indels > 0) {
	/*
	  solve_first_indel, plus: Getting genome at diagonal 2353 - querylength 36 - max_end_deletions 6 = 2311
	  q: AAAANNCATTCTCCTCCGCATAAGCCTGCGTCAGAT
	  g: GTAACATGAAAACATTCTCCTCCGCATAAGCCTGCGTCAGAT
	  first_indel_plus: 1 mismatches to right of firstbound 20: 5 -1
	  Trying insertion of 2
	  1 mismatches on left at: 4 36
	  Got indel_pos 6.
	  >AAAAnnCATTCTCCTCCGCATAAGCCTGCGTCAGAT	1	34.ins4 
	  .AAAA--------------------------------	1..4	ins:2,sub:0	+MT:2320..2323
	  ,------CATTCTCCTCCGCATAAGCCTGCGTCAGAT	7..36	ins:2,sub:0	+MT:2324..2353
          .0123456
	*/

	nmismatches = 0;
	ncolordiffs_long = 0;
	for (i = nmismatches_long-1; i >= 0; i--) {
	  if (indel_pos <= mismatch_positions[i]) { /* <= because indel_pos is in long segment */
	    nmismatches++;
	    /* ncolordiffs_long = colordiffs[i]; */
	  }
	}

	hits = List_push(hits,
			 (void *) Stage3_new_insertion(&(*found_score),indels,/*indel_pos*/indel_pos-indels,
						       /*nmismatches1*/0,/*nmismatches2*/nmismatches,
						       /*ncolordiffs1*/ncolordiffs_short,/*ncolordiffs2*/ncolordiffs_long,
						       left+indels,/*genomiclength*/querylength-indels,querylength,
						       /*plusp*/true,&(gbuffer[max_end_deletions+indels]),query,
						       chrnum,chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	debug2e(printf("successful insertion at %d with %d long + %d short mismatches\n",
		       indel_pos,nmismatches,0));
      } else {
	/*
	  solve_first_indel, plus: Getting genome at diagonal 1752395434 - querylength 36 - max_end_deletions 6 = 1752395392
	  q: GATACCTTTCTCCAAAATGGATCAAAATGTTTGTTT
	  g: CAGAAGATACCTTTTTTCCAAAATGGATCAAAATGTTTGTTT
	  first_indel_plus: 2 mismatches to right of firstbound 17: 9 5 -1
	  Entered compute_end_indels_left with length1 = 5
	  Trying deletion of -1
	  1 mismatches on left at: 9 36
	  Got indel_pos 9.
	  >GATACCTTTCTCCAAAATGGATCAAAATGTTTGTTT	1	1 0 0
	  .GATACCTTTt--------------------------	1..9	del:1,sub:0	+10:74630280..74630288
	  ,---------tTCCAAAATGGATCAAAATGTTTGTTT	10..36	del:1,sub:0	+10:74630290..74630316
          .0123456789
	*/

	nmismatches = 0;
	ncolordiffs_long = 0;
	for (i = nmismatches_long-1; i >= 0; i--) {
	  if (indel_pos <= mismatch_positions[i]) { /* <= because indel_pos is in long segment */
	    nmismatches++;
	    /* ncolordiffs_long = colordiffs[i]; */
	  }
	}

	hits = List_push(hits,
			 (void *) Stage3_new_deletion(&(*found_score),-indels,indel_pos,
						      /*nmismatches1*/0,/*nmismatches2*/nmismatches,
						      /*ncolordiffs1*/ncolordiffs_short,/*ncolordiffs2*/ncolordiffs_long,
						      left+indels,/*genomiclength*/querylength-indels,querylength,
						      /*plusp*/true,&(gbuffer[max_end_deletions+indels]),query,
						      chrnum,chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	debug2e(printf("successful deletion at %d with %d long + %d short mismatches\n",
		       indel_pos,nmismatches,0));
      }
      debug2e(if (dibasep) {
	printf("colordiffs: %d short, %d long\n",ncolordiffs_short,ncolordiffs_long);
      });
    }
  }

  return hits;
}

static List_T
solve_last_indel_plus (int *found_score, List_T hits, Genomicpos_T diagonal, int lastbound,
		       Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh,
		       Genome_T genome, UINT4 *genome_blocks, UINT4 *snp_blocks,
		       char *query, char *queryptr, int querylength,
		       Compress_T query_compress_fwd, char *gbuffer, int max_end_insertions, int max_end_deletions,
		       int min_indel_end_matches, int indel_penalty, int max_mismatches,
		       bool trim_ends_p, bool dibasep, bool cmetp) {
  int i;
  Genomicpos_T left;
  int indels, indel_pos, endgoal, nmismatches_long, ncolordiffs_short, ncolordiffs_long;
  int mismatch_positions[MAX_QUERYLENGTH];
  int colordiffs[MAX_QUERYLENGTH];
  int nmismatches;
  
  left = diagonal - querylength;
#if 0
  if (max_end_deletions > chrhigh - diagonal) {
    max_end_deletions = chrhigh - diagonal;
    /* diagonal guaranteed to be <= chrhigh, so max_end_deletions >= 0 */
  }
#else
  if (max_end_deletions > left - chroffset) {
    max_end_deletions = left - chroffset;
    /* left guaranteed to be >= chroffset, so max_end_deletions >= 0 */
  }
#endif

  debug2e(printf("\nsolve_last_indel, plus: Getting genome at diagonal %u - querylength %d + %d = %u\n",
		diagonal,querylength,max_end_deletions,left));
  debug2e(Genome_fill_buffer_simple(genome,left,querylength+max_end_deletions,gbuffer));
  debug2e(printf("q: %s\ng: %s\n",queryptr,gbuffer));

  /* No need to check chromosome bounds */
  nmismatches_long = Genome_mismatches_left(mismatch_positions,colordiffs,max_mismatches,
					    queryptr,query_compress_fwd,genome_blocks,snp_blocks,
					    left,/*pos5*/0,/*pos3*/lastbound,dibasep,cmetp,/*plusp*/true);

  debug2e(
	 printf("last_indel_plus: %d mismatches to left of lastbound %d:",nmismatches_long,lastbound);
	 for (i = 0; i <= nmismatches_long; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");

	 if (dibasep) {
	   printf("colordiffs on left:");
	   for (i = 0; i < nmismatches_long; i++) {
	     printf(" %d",colordiffs[i]);
	   }
	   printf("\n");
	 });

  if (nmismatches_long <= max_mismatches) {
    endgoal = lastbound;
  } else {
    endgoal = mismatch_positions[max_mismatches];
  }
  if (endgoal >= lastbound /* && querylength - endgoal > min_indel_end_matches */) {
    if (compute_end_indels_right(&indels,&ncolordiffs_short,&indel_pos,
				 /*length1*/querylength-endgoal,querylength,
				 left,queryptr,/*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				 min_indel_end_matches,max_end_insertions,max_end_deletions,
#if 0
				 /*max_mismatches_allowed*/max_mismatches-nmismatches_long,
#endif
				 dibasep,cmetp,/*plusp*/true) == true) {
      debug2e(printf("Got indel_pos %d\n",indel_pos));

      Genome_fill_buffer_simple(genome,left,querylength+max_end_deletions,gbuffer);
      if (indels > 0) {
	/*
	  solve_last_indel, plus: Getting genome at diagonal 2355 - querylength 36 + 6 = 2319
	  q: AAAACATTCTCCTCCGCATAAGCCTGCGTCNNAGAT
	  g: AAAACATTCTCCTCCGCATAAGCCTGCGTCAGATCAAAACAC
	  last_indel_plus: 2 mismatches to left of lastbound 16: 30 31 36
	  Trying insertion of 2
	  1 mismatches on right at: 31 -1
	  Got indel_pos 30
	  >AAAACATTCTCCTCCGCATAAGCCTGCGTCnnAGAT	1	34.ins30 
	  .AAAACATTCTCCTCCGCATAAGCCTGCGTC------	1..30	ins:2,sub:1+-1=0	+MT:2320..2349
	  ,--------------------------------AGAT	33..36	ins:2,sub:0+0=0	+MT:2350..2353
	  .0123456789012345678901234567890
	*/

	nmismatches = 0;
	ncolordiffs_long = 0;
	for (i = nmismatches_long-1; i >= 0; i--) {
	  if (indel_pos > mismatch_positions[i]) { /* > because indel_pos is in insert */
	    nmismatches++;
	    /* ncolordiffs_long = colordiffs[i]; */
	  }
	}

	hits = List_push(hits,
			 (void *) Stage3_new_insertion(&(*found_score),indels,indel_pos,
						       /*nmismatches1*/nmismatches,/*nmismatches2*/0,
						       /*ncolordiffs1*/ncolordiffs_long,
						       /*ncolordiffs2*/ncolordiffs_short,
						       left,querylength-indels,querylength,
						       /*plusp*/true,gbuffer,query,chrnum,
						       chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	debug2e(printf("successful insertion at %d with %d long + %d short mismatches\n",
		       indel_pos,nmismatches,0));
      } else {
	/*
	  solve_last_indel, plus: Getting genome at diagonal 2351 - querylength 32 + 6 = 2319
	  q: AAAACATTCTCCTCCGCATAAGCCTGTCAGAT
	  g: AAAACATTCTCCTCCGCATAAGCCTGCGTCAGATCAAA
	  last_indel_plus: 1 mismatches to left of lastbound 15: 26 32
	  Trying deletion of -2
	  1 mismatches on right at: 24 -1
	  Got indel_pos 25
	  >AAAACATTCTCCTCCGCATAAGCCTGTCAGAT	1	34.del26 
	  .AAAACATTCTCCTCCGCATAAGCCTgc-----	1..25	del:2,sub:0	+MT:2320..2344
	  ,-------------------------GTCAGAT	26..32	del:2,sub:0	+MT:2347..2353
	  .01234567890123456789012345
	*/

	nmismatches = 0;
	ncolordiffs_long = 0;
	for (i = nmismatches_long-1; i >= 0; i--) {
	  if (indel_pos > mismatch_positions[i]) { /* > because indel_pos is in deletion */
	    nmismatches++;
	    /* ncolordiffs_long = colordiffs[i]; */
	  }
	}

	hits = List_push(hits,
			 (void *) Stage3_new_deletion(&(*found_score),-indels,indel_pos,
						      /*nmismatches1*/nmismatches,/*nmismatches2*/0,
						      /*ncolordiffs1*/ncolordiffs_long,
						      /*ncolordiffs2*/ncolordiffs_short,
						      left,querylength-indels,querylength,
						      /*plusp*/true,gbuffer,query,chrnum,
						      chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	debug2e(printf("successful deletion at %d with %d long + %d short mismatches\n",
		       indel_pos,nmismatches,0));
      }
      debug2e(if (dibasep) {
	printf("colordiffs: %d long, %d short\n",ncolordiffs_long,ncolordiffs_short);
      });
    }
  }
    
  return hits;
}


static List_T
solve_first_indel_minus (int *found_score, List_T hits, Genomicpos_T diagonal, int lastbound, 
			 Chrnum_T chrnum, Genomicpos_T chroffset,
			 Genome_T genome, UINT4 *genome_blocks, UINT4 *snp_blocks,
			 char *query, char *queryptr, int querylength,
			 Compress_T query_compress_rev, char *gbuffer, int max_end_insertions, int max_end_deletions,
			 int min_indel_end_matches, int indel_penalty, int max_mismatches,
			 bool trim_ends_p, bool dibasep, bool cmetp) {
  int i;
  Genomicpos_T left;
  int indels, indel_pos, endgoal, nmismatches_long, ncolordiffs_short, ncolordiffs_long;
  int mismatch_positions[MAX_QUERYLENGTH];
  int colordiffs[MAX_QUERYLENGTH];
  int nmismatches;
  
  left = diagonal - querylength;
  if (max_end_deletions > left - chroffset) {
    max_end_deletions = left - chroffset;
    /* left guaranteed to be >= chroffset, so max_end_deletions >= 0 */
  }

  debug2e(printf("\nsolve_first_indel, minus: Getting genome at diagonal %u + 12 - querylength %d = %u.\n",
		diagonal,querylength,left));
  debug2e(Genome_fill_buffer_simple(genome,left,querylength+max_end_deletions,gbuffer));
  debug2e(printf("q.rc: %s\ng:    %s\n",queryptr,gbuffer));

  /* No need to check chromosome bounds */
  nmismatches_long = Genome_mismatches_left(mismatch_positions,colordiffs,max_mismatches,
					    queryptr,query_compress_rev,genome_blocks,snp_blocks,
					    left,/*pos5*/0,/*pos3*/lastbound,dibasep,cmetp,/*plusp*/false);

  debug2e(
	 printf("first_indel_minus: %d mismatches to left of lastbound %d:",nmismatches_long,lastbound);
	 for (i = 0; i < nmismatches_long; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");

	 if (dibasep) {
	   printf("colordiffs on left:");
	   for (i = 0; i < nmismatches_long; i++) {
	     printf(" %d",colordiffs[i]);
	   }
	   printf("\n");
	 });

  if (nmismatches_long <= max_mismatches) {
    endgoal = lastbound;
  } else {
    endgoal = mismatch_positions[max_mismatches];
  }
  if (endgoal >= lastbound /* && querylength - endgoal > min_indel_end_matches */) {
    if (compute_end_indels_right(&indels,&ncolordiffs_short,&indel_pos,
				 /*length1*/querylength-endgoal,querylength,
				 left,queryptr,/*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				 min_indel_end_matches,max_end_insertions,max_end_deletions,
#if 0
				 /*max_mismatches_allowed*/max_mismatches-nmismatches_long,
#endif
				 dibasep,cmetp,/*plusp*/false) == true) {
      debug2e(printf("Got indel_pos %d\n",indel_pos));

      Genome_fill_buffer_simple(genome,left,querylength+max_end_deletions,gbuffer);
      if (indels > 0) {
	/*
	  solve_first_indel, minus: Getting genome at diagonal 2355 + 12 - querylength 36 = 2319.
	  q.rc: AAAACATTCTCCTCCGCATAAGCCTGCGTCNNAGAT
	  g:    AAAACATTCTCCTCCGCATAAGCCTGCGTCAGATCAAAACAC
	  first_indel_minus: 1 mismatches to left of lastbound 19: 30
	  Trying insertion of 2
	  1 mismatches on right at: 31 -1
	  Got indel_pos 30
	  >ATCTnnGACGCAGGCTTATGCGGAGGAGAATGTTTT	1	34.ins30  REVCOMP
	  .ATCT--------------------------------	1..4	ins:2,sub:0	-MT:2353..2350
	  ,------GACGCAGGCTTATGCGGAGGAGAATGTTTT	7..36	ins:2,sub:1	-MT:2349..2320
	  ......0987654321098765432109876543210
	*/

	nmismatches = 0;
	ncolordiffs_long = 0;
	for (i = nmismatches_long-1; i >= 0; i--) {
	  if (indel_pos > mismatch_positions[i]) { /* > because indel_pos is in insert */
	    nmismatches++;
	    /* ncolordiffs_long = colordiffs[i]; */
	  }
	}

	hits = List_push(hits,
			 (void *) Stage3_new_insertion(&(*found_score),indels,/*indel_pos*/querylength-indel_pos-indels,
						       /*nmismatches1*/0,/*nmismatches2*/nmismatches,
						       /*ncolordiffs1*/ncolordiffs_short,
						       /*ncolordiffs2*/ncolordiffs_long,
						       left,querylength-indels,querylength,
						       /*plusp*/false,gbuffer,query,chrnum,
						       chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	debug2e(printf("successful insertion at %d with %d long + %d short mismatches\n",
		       indel_pos,nmismatches,0));
      } else {
	/*
	  solve_first_indel, minus: Getting genome at diagonal 2351 + 12 - querylength 32 = 2319.
	  q.rc: AAAACATTCTCCTCCGCATAAGCCTGTCAGAT
	  g:    AAAACATTCTCCTCCGCATAAGCCTGCGTCAGATCAAA
	  first_indel_minus: 1 mismatches to left of lastbound 15: 26
	  Trying deletion of -2
	  1 mismatches on right at: 24 -1
	  Got indel_pos 25
	  successful deletion at 25 with 0 long + 0 short mismatches
	  >ATCTGACAGGCTTATGCGGAGGAGAATGTTTT	1	34.del26  REVCOMP
	  .ATCTGACgc-----------------------	1..7	del:2,sub:0	-MT:2353..2347
	  ,-------AGGCTTATGCGGAGGAGAATGTTTT	8..32	del:2,sub:0	-MT:2344..2320
	  .......54321098765432109876543210
	*/

	nmismatches = 0;
	ncolordiffs_long = 0;
	for (i = nmismatches_long-1; i >= 0; i--) {
	  if (indel_pos > mismatch_positions[i]) { /* > because mismatch is in the deleted part */
	    nmismatches++;
	    /* ncolordiffs_long = colordiffs[i]; */
	  }
	}

	hits = List_push(hits,
			 (void *) Stage3_new_deletion(&(*found_score),-indels,querylength-indel_pos,
						      /*nmismatches1*/0,/*nmismatches2*/nmismatches,
						      /*ncolordiffs1*/ncolordiffs_short,
						      /*ncolordiffs2*/ncolordiffs_long,
						      left,querylength-indels,querylength,
						      /*plusp*/false,gbuffer,query,chrnum,
						      chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	debug2e(printf("successful deletion at %d with %d long + %d short mismatches\n",
		       indel_pos,nmismatches,0));
      }
      debug2e(if (dibasep) {
	printf("colordiffs: %d short, %d long\n",ncolordiffs_short,ncolordiffs_long);
      });
    }
  }

  debug2e(printf("done\n"));
  return hits;
}

static List_T
solve_last_indel_minus (int *found_score, List_T hits, Genomicpos_T diagonal, int firstbound,
			Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh,
			Genome_T genome, UINT4 *genome_blocks, UINT4 *snp_blocks,
			char *query, char *queryptr, int querylength,
			Compress_T query_compress_rev, char *gbuffer,
			int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
			int indel_penalty,int max_mismatches, bool trim_ends_p, bool dibasep, bool cmetp) {
  int i;
  Genomicpos_T left;
  int indels, indel_pos, endgoal, nmismatches_long, ncolordiffs_short, ncolordiffs_long;
  int mismatch_positions[MAX_QUERYLENGTH];
  int colordiffs[MAX_QUERYLENGTH];
  int nmismatches;
  
  left = diagonal - querylength;
#if 0
  /* Doesn't work for 
     GTTAGGGTTAGGGTTAGAGGTTAGGGGTTAGGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG
     CTAACCCTAACCCTAACCCTAACCCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCC
  */

  if (max_end_deletions > chrhigh - diagonal) {
    max_end_deletions = chrhigh - diagonal;
    /* diagonal guaranteed to be <= chrhigh, so max_end_deletions >= 0 */
  }
#else
  if (max_end_deletions > left - chroffset) {
    max_end_deletions = left - chroffset;
    /* left guaranteed to be >= chroffset, so max_end_deletions >= 0 */
  }
#endif

  
  debug2e(printf("\nsolve_last_indel, minus: Getting genome at diagonal %u + 12 - querylength %d = %u, max_end_deletions = %d.\n",
		 diagonal,querylength,left,max_end_deletions));
  debug2e(Genome_fill_buffer_simple(genome,left-max_end_deletions,querylength+max_end_deletions,gbuffer));
  debug2e(printf("q.rc: %s\ng:    %s\n",queryptr,gbuffer));

  /* No need to check chromosome bounds */
  nmismatches_long = Genome_mismatches_right(mismatch_positions,colordiffs,max_mismatches,
					     queryptr,query_compress_rev,genome_blocks,snp_blocks,
					     left,/*pos5*/firstbound,/*pos3*/querylength,
					     dibasep,cmetp,/*plusp*/false);

  debug2e(
	 printf("last_indel_minus: %d mismatches to right of firstbound %d:",nmismatches_long,firstbound);
	 for (i = 0; i < nmismatches_long; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");

	 if (dibasep) {
	   printf("colordiffs on right:");
	   for (i = 0; i < nmismatches_long; i++) {
	     printf(" %d",colordiffs[i]);
	   }
	   printf("\n");
	 });

  if (nmismatches_long <= max_mismatches) {
    endgoal = firstbound;
  } else {
    endgoal = mismatch_positions[max_mismatches];
  }
  if (endgoal <= firstbound /* && endgoal > min_indel_end_matches */) {
    if (compute_end_indels_left(&indels,&ncolordiffs_short,&indel_pos,/*length1*/endgoal,querylength,
				left,queryptr,/*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				min_indel_end_matches,max_end_insertions,max_end_deletions,
#if 0
				/*max_mismatches_allowed*/max_mismatches-nmismatches_long,
#endif
				dibasep,cmetp,/*plusp*/false) == true) {
      debug2e(printf("Got indel_pos %d\n",indel_pos));

      Genome_fill_buffer_simple(genome,left-max_end_deletions,querylength+max_end_deletions,gbuffer);
      if (indels > 0) {
	/*
	  solve_last_indel, minus: Getting genome at diagonal 2353 + 12 - querylength 36 = 2317.
	  q.rc: AAAANNCATTCTCCTCCGCATAAGCCTGCGTCAGAT
	  g:    GTAACATGAAAACATTCTCCTCCGCATAAGCCTGCGTCAGAT
	  last_indel_minus: 1 mismatches to right of firstbound 17: 5
	  Trying insertion of 2
	  1 mismatches on left at: 4 36
	  Got indel_pos 6
	  >ATCTGACGCAGGCTTATGCGGAGGAGAATGnnTTTT	1	34.ins4  REVCOMP
	  .ATCTGACGCAGGCTTATGCGGAGGAGAATG------	1..30	ins:2,sub:0	-MT:2353..2324
	  ,--------------------------------TTTT	33..36	ins:2,sub:0	-MT:2323..2320
	  ..............................6543210
	*/

	nmismatches = 0;
	ncolordiffs_long = 0;
	for (i = nmismatches_long-1; i >= 0; i--) {
	  if (indel_pos <= mismatch_positions[i]) { /* <= because indel_pos is in long segment */
	    nmismatches++;
	    /* ncolordiffs_long = colordiffs[i]; */
	  }
	}

	hits = List_push(hits,
			 (void *) Stage3_new_insertion(&(*found_score),indels,querylength-indel_pos,
						       /*nmismatches1*/nmismatches,/*nmismatches2*/0,
						       /*ncolordiffs1*/ncolordiffs_long,
						       /*ncolordiffs2*/ncolordiffs_short,
						       left+indels,querylength-indels,querylength,
						       /*plusp*/false,&(gbuffer[max_end_deletions+indels]),query,
						       chrnum,chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	debug2e(printf("successful insertion at %d with %d long + %d short mismatches\n",
		       indel_pos,nmismatches,0));
      } else {
	/*
	  solve_last_indel, minus: Getting genome at diagonal 1752395434 + 12 - querylength 36 = 1752395398.
	  q.rc: GATACCTTTCTCCAAAATGGATCAAAATGTTTGTTT
	  g:    CAGAAGATACCTTTTTTCCAAAATGGATCAAAATGTTTGTTT
	  last_indel_minus: 2 mismatches to right of firstbound 17: 9 5
	  Entered compute_end_indels_left with length1 = 5
	  Trying deletion of -1
	  1 mismatches on left at: 9 36
	  Got indel_pos 9
	  >AAACAAACATTTTGATCCATTTTGGAGAAAGGTATC	1	1 0 0 REVCOMP
	  .AAACAAACATTTTGATCCATTTTGGAaa--------	1..27	del:1,sub:0	-10:74630316..74630290
	  ,---------------------------AAAGGTATC	28..36	del:1,sub:0	-10:74630288..74630280
	  ...........................9876543210
	*/

	nmismatches = 0;
	ncolordiffs_long = 0;
	for (i = nmismatches_long-1; i >= 0; i--) {
	  if (indel_pos <= mismatch_positions[i]) { /* <= because indel_pos is in long segment */
	    nmismatches++;
	    /* ncolordiffs_long = colordiffs[i]; */
	  }
	}

	hits = List_push(hits,(void *) Stage3_new_deletion(&(*found_score),-indels,querylength-indel_pos,
							   /*nmismatches1*/nmismatches,/*nmismatches2*/0,
							   /*ncolordiffs1*/ncolordiffs_long,
							   /*ncolordiffs2*/ncolordiffs_short,
							   left+indels,querylength-indels,querylength,
							   /*plusp*/false,&(gbuffer[max_end_deletions+indels]),query,
							   chrnum,chroffset,indel_penalty,trim_ends_p,dibasep,cmetp));
	debug2e(printf("successful deletion at %d with %d long + %d short mismatches\n",
		       indel_pos,nmismatches,0));
      }
      debug2e(if (dibasep) {
	printf("colordiffs: %d long, %d short\n",ncolordiffs_long,ncolordiffs_short);
      });
    }
  }
    
  return hits;
}


static List_T
find_end_indels (int *found_score, List_T hits,
		 struct Indelsegment_T *plus_segments, struct Indelsegment_T *minus_segments,
		 int plus_nsegments, int minus_nsegments,
		 Genome_T genome, UINT4 *genome_blocks, UINT4 *snp_blocks,
		 char *queryuc_ptr, char *queryrc,
		 int querylength, int firstbound, int lastbound,
		 Compress_T query_compress_fwd, Compress_T query_compress_rev,
		 int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		 int indel_penalty, int max_mismatches_allowed, bool trim_ends_p, bool dibasep, bool cmetp) {
  char *gbuffer;
  Indelsegment_T ptr;

  debug(printf("*** find_end_indels with max_mismatches_allowed %d ***\n",
	       max_mismatches_allowed));

  gbuffer = (char *) CALLOC(MAX_QUERYLENGTH + max_end_deletions + 1,sizeof(char));

  for (ptr = plus_segments; ptr < &(plus_segments[plus_nsegments]); ptr++) {

    if (ptr->floor_xfirst <= max_mismatches_allowed) {
      /* First indel, plus */
      debug2e(printf("floor_xfirst %d <= mismatches allowed %d\n",ptr->floor_xfirst,max_mismatches_allowed));
      hits = solve_first_indel_plus(&(*found_score),hits,ptr->diagonal,firstbound,
				    ptr->chrnum,ptr->chroffset,
				    genome,genome_blocks,snp_blocks,
				    /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
				    query_compress_fwd,gbuffer,
				    max_end_insertions,max_end_deletions,min_indel_end_matches,
				    indel_penalty,max_mismatches_allowed,trim_ends_p,dibasep,cmetp);
    }

    if (ptr->floor_xlast <= max_mismatches_allowed) {
      /* Last indel, plus */
      debug2e(printf("floor_xlast %d <= mismatches allowed %d\n",ptr->floor_xlast,max_mismatches_allowed));
      hits = solve_last_indel_plus(&(*found_score),hits,ptr->diagonal,lastbound,
				   ptr->chrnum,ptr->chroffset,ptr->chrhigh,
				   genome,genome_blocks,snp_blocks,
				   /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
				   query_compress_fwd,gbuffer,
				   max_end_insertions,max_end_deletions,min_indel_end_matches,
				   indel_penalty,max_mismatches_allowed,trim_ends_p,dibasep,cmetp);
    }
  }

  for (ptr = minus_segments; ptr < &(minus_segments[minus_nsegments]); ptr++) {

    if (ptr->floor_xfirst <= max_mismatches_allowed) { 
      /* First indel, minus */
      debug2e(printf("floor_xfirst %d <= mismatches allowed %d\n",ptr->floor_xfirst,max_mismatches_allowed));
      hits = solve_first_indel_minus(&(*found_score),hits,ptr->diagonal,lastbound,
				     ptr->chrnum,ptr->chroffset,
				     genome,genome_blocks,snp_blocks,
				     /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
				     query_compress_rev,gbuffer,
				     max_end_insertions,max_end_deletions,min_indel_end_matches,
				     indel_penalty,max_mismatches_allowed,trim_ends_p,dibasep,cmetp);
    }

    if (ptr->floor_xlast <= max_mismatches_allowed) {
      /* Last indel, minus */
      debug2e(printf("floor_xlast %d <= mismatches allowed %d\n",ptr->floor_xlast,max_mismatches_allowed));
      hits = solve_last_indel_minus(&(*found_score),hits,ptr->diagonal,firstbound,
				    ptr->chrnum,ptr->chroffset,ptr->chrhigh,
				    genome,genome_blocks,snp_blocks,
				    /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
				    query_compress_rev,gbuffer,
				    max_end_insertions,max_end_deletions,min_indel_end_matches,
				    indel_penalty,max_mismatches_allowed,trim_ends_p,dibasep,cmetp);
    }
  }

  FREE(gbuffer);

  return hits;
}



/************************************************************************
 *   Splicing
 ************************************************************************/



Genomicpos_T *
Stage1_retrieve_splicesites (Splicetype_T **splicetypes, int *nsplicesites,
			     IIT_T splicesites_iit, int *splicesites_divint_crosstable,
			     int donor_typeint, int acceptor_typeint, IIT_T chromosome_iit) {
  Genomicpos_T *splicesites, chrlength, chroffset, position;
  Genomicpos_T last_donor, last_antidonor, last_acceptor, last_antiacceptor;
  int *splicesites1;
  int divno, nsplicesites1, space, i, k;
  Chrnum_T chrnum;
  Interval_T *intervals, interval;

  space = 0;
  for (chrnum = 1; chrnum <= IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicesites_divint_crosstable[chrnum]) > 0) {
      space += IIT_nintervals(splicesites_iit,divno);
      debug(printf("divno %d, nsplicesites %d\n",divno,IIT_nintervals(splicesites_iit,divno)));
    }
  }
  debug(printf("total splicesites with duplicates: %d\n",space));

  if (space == 0) {
    splicesites = (Genomicpos_T *) NULL;
    *splicetypes = (Splicetype_T *) NULL;
  } else {
    splicesites = (Genomicpos_T *) CALLOC(space,sizeof(Genomicpos_T));
    *splicetypes = (Splicetype_T *) CALLOC(space,sizeof(Splicetype_T));
  }

  k = 0;
  for (chrnum = 1; chrnum <= IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicesites_divint_crosstable[chrnum]) > 0) {
      chroffset = IIT_interval_low(chromosome_iit,chrnum);
      chrlength = IIT_length(chromosome_iit,chrnum);
      splicesites1 = IIT_get_with_divno(&nsplicesites1,splicesites_iit,divno,
					0U,chrlength,/*sortp*/false);
      if (nsplicesites1 > 0) {
	intervals = (Interval_T *) CALLOC(nsplicesites1,sizeof(Interval_T));
	for (i = 0; i < nsplicesites1; i++) {
	  intervals[i] = &(splicesites_iit->intervals[divno][i]);
	}
	qsort(intervals,nsplicesites1,sizeof(Interval_T),Interval_cmp);
	  
	last_donor = last_antidonor = last_acceptor = last_antiacceptor = 0U;
	for (i = 0; i < nsplicesites1; i++) {
	  interval = intervals[i];
	  position = Interval_low(intervals[i]) + chroffset;
	  if (Interval_type(interval) == donor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position != last_donor) {
		last_donor = splicesites[k] = position;
		(*splicetypes)[k++] = DONOR;
	      }
	    } else {
	      if (position != last_antidonor) {
		last_antidonor = splicesites[k] = position;
		(*splicetypes)[k++] = ANTIDONOR;
	      }
	    }
	  } else if (Interval_type(interval) == acceptor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position != last_acceptor) {
		last_acceptor = splicesites[k] = position;
		(*splicetypes)[k++] = ACCEPTOR;
	      }
	    } else {
	      if (position != last_antiacceptor) {
		last_antiacceptor = splicesites[k] = position;
		(*splicetypes)[k++] = ANTIACCEPTOR;
	      }
	    }
	  }
	}
	FREE(intervals);
	FREE(splicesites1);
      }
    }
  }

  *nsplicesites = k;
  debug(printf("total splicesites without duplicates: %d\n",*nsplicesites));

  return splicesites;
}

/************************************************************************
 *  Solving known splice sites using dual search
 ************************************************************************/

#if 0

typedef struct Splice_T *Splice_T;
struct Splice_T {
  Genomicpos_T left;
  Genomicpos_T splicesite;
  int splicepos;
  int querypos;
};

static Splice_T
Splice_new (Genomicpos_T left, Genomicpos_T splicesite) {
  Splice_T new = (Splice_T) MALLOC(sizeof(*new));

  if (splicesite < left) {
    printf("PROBLEM! splicesite %u < left %u\n",splicesite,left);
    abort();
  }
  new->left = left;
  new->splicesite = splicesite;
#if 0
  /* Now assigned as a batch for each querypos */
  if (plusp) {
    new->splicepos = splicesite - left;
  } else {
    new->splicepos = querylength - (splicesite - left);
  }
  new->querypos = querypos;
#endif
  return new;
}

static void
Splice_free (Splice_T *old) {
  FREE(*old);
  return;
}

static int
Splice_cmp (const void *a, const void *b) {
  Splice_T x = * (Splice_T *) a;
  Splice_T y = * (Splice_T *) b;

  if (x->left < y->left) {
    return -1;
  } else if (y->left < x->left) {
    return +1;
  } else if (x->splicesite < y->splicesite) {
    return -1;
  } else if (y->splicesite < x->splicesite) {
    return +1;
  } else if (x->querypos < y->querypos) {
    return -1;
  } else if (y->querypos < x->querypos) {
    return +1;
  } else {
    return 0;
  }
}


/* Same as littleendian version.  Assumes positions already converted. */
static int
binary_search_noendian (int lowi, int highi, Genomicpos_T *positions, Genomicpos_T goal) {
  int middlei;

  debug10(printf("entered binary_search_noendian with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = (lowi+highi)/2;
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi,positions[highi],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}


/* positions1 has splicesites, positions2 has positions */
static void
dual_search (List_T *splices, Genomicpos_T *positions1, int highi,
	     Genomicpos_T *positions2, int highj,
	     Splicetype_T *splicetypes, int querylength, int diagterm) {
  int middlei, middlej, starti, startj, endi, endj, i, j;
  bool changep;
  Genomicpos_T splicesite, position, left, right, splicesite_adjleft, splicesite_adjright;
  Splicetype_T splicetype;

  if (highi <= /*lowi*/ 0) {
    return;
  } else if (highj <= /*lowj*/ 0) {
    return;
  } else {
    if (highi /* - lowi*/ < highj /* - lowj*/) {
      /* Given site, search against positions */
      middlei = highi/2;	/* (lowi+highi)/2 */
      splicesite = positions1[middlei];

      /* Search splicesite against right of intervals: position -
	 diagterm + querylength > splicesite, so position > splicesite + diagterm - querylength */
      middlej = binary_search(/*lowj*/0,highj,positions2,splicesite+diagterm-querylength);
#ifdef WORDS_BIGENDIAN
      debug12(printf("  Binary search of splicesite (#%d..%d) #%d:%u (adjusted to %u) returns position #%d:%u\n",
		     /*lowi*/0,highi,middlei,splicesite,splicesite+diagterm-querylength,middlej,
		     Bigendian_convert_uint(positions2[middlej])));
#else
      debug12(printf("  Binary search of splicesite (#%d..%d) #%d:%u (adjusted to %u) returns position #%d:%u\n",
		     /*lowi*/0,highi,middlei,splicesite,splicesite+diagterm-querylength,middlej,positions2[middlej]));
#endif

    } else {
      /* Given position, search against sites */
      middlej = highj/2;	/* (lowj+highj)/2; */

#ifdef WORDS_BIGENDIAN
      position = Bigendian_convert_uint(positions2[middlej]);
#else
      position = positions2[middlej];
#endif

      /* Search left of interval against splicesites: splicesite > position - diagterm */
      middlei = binary_search_noendian(/*lowi*/0,highi,positions1,/*left*/position-diagterm);
      debug12(printf("  Binary search of position (#%d..%d) #%d:%u (adjusted to %u) returns splicesite #%d:%u\n",
		     /*lowj*/0,highj,middlej,position,position-diagterm,middlei,positions1[middlei]));
    }

    /* Find starts */
    starti = middlei;
    startj = middlej;
    changep = true;
    while (changep) {
      changep = false;

      /* Scan backward on position as long as right > splicesite */
      /* positions2 - diagterm + querylength > positions1 => positions2 > positions1 + diagterm - querylength */
      splicesite_adjright = positions1[starti] + diagterm - querylength;
#ifdef WORDS_BIGENDIAN
      while (startj > 0 /*startj-1 >= lowj*/ && Bigendian_convert_uint(positions2[startj-1]) > splicesite_adjright) {
	startj--;
	changep = true;
      }
#else
      while (startj > 0 /*startj-1 >= lowj*/ && positions2[startj-1] > splicesite_adjright) {
	startj--;
	changep = true;
      }
#endif

      /* Scan backward on site as long as splicesite > left */
#ifdef WORDS_BIGENDIAN
      left = Bigendian_convert_uint(positions2[startj]) - diagterm;
#else
      left = positions2[startj] - diagterm;
#endif
      while (starti > 0 /*starti-1 >= lowi*/ && positions1[starti-1] > left) {
	starti--;
	changep = true;
      }
    }


    /* Find ends */
    endi = middlei;
    endj = middlej;
    changep = true;
    while (changep) {
      changep = false;

      /* Scan forward on site as long as splicesite < right */
#ifdef WORDS_BIGENDIAN
      right = Bigendian_convert_uint(positions2[endj]) - diagterm + querylength;
#else
      right = positions2[endj] - diagterm + querylength;
#endif
      while (endi+1 < highi && positions1[endi+1] < right) {
	endi++;
	changep = true;
      }

      /* Scan forward on position as long as left < splicesite */
      /* positions2 - diagterm < positions1 => positions2 < positions1 + diagterm */
      splicesite_adjleft = positions1[endi] + diagterm;
#ifdef WORDS_BIGENDIAN
      while (endj+1 < highj && Bigendian_convert_uint(positions2[endj+1]) < splicesite_adjleft) {
	endj++;
	changep = true;
      }
#else
      while (endj+1 < highj && positions2[endj+1] < splicesite_adjleft) {
	endj++;
	changep = true;
      }
#endif
    }


    /* Recurse on low part */
    /* dual_search(splices,positions1,lowi,starti,positions2,lowj,startj,splicetypes,
       querylength,diagterm); */
    dual_search(splices,positions1,starti,positions2,startj,splicetypes,querylength,diagterm);


    /* Process middle part */
    for (i = starti; i <= endi; i++) {
      splicesite = positions1[i];
      splicesite_adjleft = splicesite + diagterm;
      splicesite_adjright = splicesite_adjleft - querylength;
      for (j = startj; j <= endj; j++) {
	/* left < splicesite => positions2 - diagterm < splicesite => positions2 < splicesite + diagterm */
	/* splicesite < right => right > splicesite => positions2 - diagterm + querylength > splicesite =>
	   positions2 > splicesite + diagterm - querylength */
	if (positions2[j] < splicesite_adjleft && positions2[j] > splicesite_adjright) {
	  debug12(printf("splicesite #%d:%u, left #%d:%u type:%d\n",
			 i,splicesite,j,positions2[j]-diagterm,splicetypes[i]));
	  splicetype = splicetypes[i];
	  splices[splicetype] = 
	    List_push(splices[splicetype],Splice_new(/*left*/positions2[j]-diagterm,splicesite));
	}
      }
    }

    /* Recurse on high part */
    /* dual_search(splices,positions1,endi+1,highi,positions2,endj+1,highj,splicetypes,
       querylength,diagterm); */
    dual_search(splices,&(positions1[endi+1]),highi-(endi+1),
		&(positions2[endj+1]),highj-(endj+1),&(splicetypes[endi+1]),
		querylength,diagterm);

    return;
  }
}




static void
gather_spliceends_known (List_T **spliceends, List_T splices,
			 UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome,
			 IIT_T chromosome_iit, char *queryuc_ptr, char *queryrc, int querylength,
			 int query_lastpos, char *gbuffer, int max_mismatches_allowed, 
			 Compress_T query_compress, Floors_T floors,
			 bool dibasep, bool cmetp, bool donorp, bool plusp, bool sensep, bool left_floor_p, bool left_limit_p) {
  List_T filtered = NULL, p;
  Splice_T *array, splice;
  Substring_T hit;
  Genomicpos_T left, splicesite, prev_left, prev_splicesite;
  Genomicpos_T chroffset, chrhigh = 0U;
  Chrnum_T chrnum;
  int nmismatches, ncolordiffs, nunknowns, n, i;
  int querypos, splicepos, prev_querypos, prev_splicepos;
  int floor;
  int *floors_from_neg3, *floors_to_pos3;

  filtered = (List_T) NULL;
  if (left_floor_p) {
    prev_left = 0U;
    prev_querypos = 0;
    for (p = splices; p != NULL; p = p->rest) {
      splice = (Splice_T) p->first;
      if (splice->splicepos < INDEX1PART) {
	debug4(printf("Removing splice at left %u, splicesite %u, splicepos %d because too close to edge\n",
		      splice->left,splice->splicesite,splice->splicepos));
	Splice_free(&splice);
      } else if (splice->querypos > splice->splicepos - INDEX1PART) {
	debug4(printf("Removing splice at left %u, splicesite %u, splicepos %d because of querypos %d\n",
		      splice->left,splice->splicesite,splice->splicepos,splice->querypos));
	Splice_free(&splice);
      } else {
	filtered = List_push(filtered,splice);
      }
    }
  } else {
    prev_left = 0U;
    prev_querypos = 0;
    for (p = splices; p != NULL; p = p->rest) {
      splice = (Splice_T) p->first;
      if (splice->splicepos > query_lastpos) {
	debug4(printf("Removing splice at left %u, splicesite %u, splicepos %d because too close to edge\n",
		      splice->left,splice->splicesite,splice->splicepos));
	Splice_free(&splice);
      } else if (splice->querypos < splice->splicepos) {
	debug4(printf("Removing splice at left %u, splicesite %u, splicepos %d because of querypos %d\n",
		      splice->left,splice->splicesite,splice->splicepos,splice->querypos));
	Splice_free(&splice);
      } else {
	filtered = List_push(filtered,splice);
      }
    }
  }

  if ((n = List_length(filtered)) == 0) {
    return;
  } else {
    array = (Splice_T *) List_to_array(filtered,NULL);
    qsort(array,n,sizeof(Splice_T),Splice_cmp);
    List_free(&filtered);
  }

  /* Start of splicesite */
  splice = array[0];
  prev_left = splice->left;
  prev_splicesite = splice->splicesite;
  prev_splicepos = splice->splicepos;
  prev_querypos = splice->querypos;
  Splice_free(&splice);
  if (left_floor_p) {
    floor = floors_from_neg3[prev_querypos] /* floors->score[-3][prev_querypos] */;
  } else {
    floor = floors->scorefrom[prev_splicepos][prev_querypos];
  }
  debug4(printf("Start %s, %s, %s: left %u, splicesite %u, querypos %d, floor %d\n",
		donorp ? "donor" : "acceptor",plusp ? "plus" : "minus",
		sensep ? "sense" : "antisense",prev_left,prev_splicesite,prev_querypos,floor));

  for (i = 1; i < n; i++) {
    splice = array[i];
    left = splice->left;
    splicesite = splice->splicesite;
    splicepos = splice->splicepos;
    querypos = splice->querypos;

    if (left == prev_left && splicesite == prev_splicesite) {
	/* Continuing splice site */
      floor += floors->scorefrom[prev_querypos][querypos];
      debug4(printf("Cont. %s, %s, %s: left %u, splicesite %u, querypos %d->%d, floor %d\n",
		    donorp ? "donor" : "acceptor",plusp ? "plus" : "minus",
		    sensep ? "sense" : "antisense",prev_left,prev_splicesite,prev_querypos,querypos,floor));
      prev_querypos = querypos;
    } else {
      /* End of splicesite */
      if (left_floor_p) {
	floor += floors->scorefrom[prev_querypos][prev_splicepos-INDEX1PART+3];
      } else {
	floor += floors_to_pos3[prev_querypos] /*floors->score[prev_querypos][query_lastpos+3] */;
      }
      debug4(printf("End:  %s, %s, %s: splicesite %u => splicepos %d (floor %d)\n",
		    donorp ? "donor" : "acceptor",plusp ? "plus" : "minus",
		    sensep ? "sense" : "antisense",prev_splicesite,prev_splicepos,floor));
      if (floor <= max_mismatches_allowed) {
	if (left_limit_p) {
	  if (plusp) {
	    debug4(printf("Testing from 0 to %d\n",prev_splicepos));
	    nmismatches = Genome_count_mismatches_limit(&ncolordiffs,snp_blocks,query_compress,genome_blocks,
							queryuc_ptr,prev_left,/*pos5*/0,/*pos3*/prev_splicepos,max_mismatches_allowed,
							dibasep,cmetp,plusp);
	  } else {
	    debug4(printf("Testing from 0 to %d\n",querylength-prev_splicepos));
	    nmismatches = Genome_count_mismatches_limit(&ncolordiffs,snp_blocks,query_compress,genome_blocks,
							queryuc_ptr,prev_left,/*pos5*/0,/*pos3*/querylength-prev_splicepos,max_mismatches_allowed,
							dibasep,cmetp,plusp);
	  }
	} else {
	  if (plusp) {
	    debug4(printf("Testing from %d to %d\n",prev_splicepos,querylength));
	    nmismatches = Genome_count_mismatches_limit(&ncolordiffs,snp_blocks,query_compress,genome_blocks,
							queryuc_ptr,prev_left,/*pos5*/prev_splicepos,/*pos3*/querylength,max_mismatches_allowed,
							dibasep,cmetp,plusp);
	  } else {
	    debug4(printf("Testing from %d to %d\n",querylength-prev_splicepos,querylength));
	    nmismatches = Genome_count_mismatches_limit(&ncolordiffs,snp_blocks,query_compress,genome_blocks,
							queryuc_ptr,prev_left,/*pos5*/querylength-prev_splicepos,/*pos3*/querylength,max_mismatches_allowed,
							dibasep,cmetp,plusp);
	  }
	}
	debug4(Genome_fill_buffer_simple(genome,prev_left,querylength,gbuffer));
	debug4(if (plusp) {
		 printf("q: %s\ng: %s\n",queryuc_ptr,gbuffer);
	       } else {
		 printf("q.rc: %s\ng:    %s\n",queryrc,gbuffer);
	       });
	debug4(printf("Test: %s, %s, %s: splicesite %u => splicepos %d (%d mismatches)\n",
		      donorp ? "donor" : "acceptor",plusp ? "plus" : "minus",
		      sensep ? "sense" : "antisense",prev_splicesite,prev_splicepos,nmismatches));
	if (nmismatches <= max_mismatches_allowed) {
	  if (prev_left + querylength > chrhigh) {
	    /* update chromosome bounds, based on low end */
	    chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,prev_left,prev_left);
	    IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
	    chrhigh += 1U;
	  }
	  Genome_fill_buffer_simple(genome,prev_left,querylength,gbuffer);
	  if (donorp) {
	    hit = Substring_new_donor(prev_splicepos,nmismatches,/*prob*/2.0,prev_left,
				      querylength,plusp,sensep,
				      gbuffer,queryuc_ptr,chrnum,chroffset,/*knownp*/true,
				      trim_ends_p,dibasep,cmetp);
	  } else {
	    hit = Substring_new_acceptor(prev_splicepos,nmismatches,/*prob*/2.0,prev_left,
					 querylength,plusp,sensep,
					 gbuffer,queryuc_ptr,chrnum,chroffset,/*knownp*/true,
					 trim_ends_p,dibasep,cmetp);
	  }
	  (*spliceends)[nmismatches] = List_push((*spliceends)[nmismatches],(void *) hit);
	}
      }

      /* Start of splicesite */
      prev_left = left;
      prev_splicesite = splicesite;
      prev_splicepos = splicepos;
      prev_querypos = querypos;
      if (left_floor_p) {
	floor = floors_from_neg3[prev_querypos] /* floors->score[-3][prev_querypos] */;
      } else {
	floor = floors->scorefrom[prev_splicepos][prev_querypos];
      }
      debug4(printf("Start %s, %s, %s: left %u, splicesite %u, querypos %d, floor %d\n",
		    donorp ? "donor" : "acceptor",plusp ? "plus" : "minus",
		    sensep ? "sense" : "antisense",prev_left,prev_splicesite,prev_querypos,floor));
    }

    Splice_free(&splice);
  }

  /* End of splicesite */
  if (left_floor_p) {
    floor += floors->scorefrom[prev_querypos][prev_splicepos-INDEX1PART+3];
  } else {
    floor += floors_to_pos3[prev_querypos]  /* floors->score[prev_querypos][query_lastpos+3] */;
  }
  debug4(printf("End:  %s, %s, %s: splicesite %u => splicepos %d (floor %d)\n",
		donorp ? "donor" : "acceptor",plusp ? "plus" : "minus",
		sensep ? "sense" : "antisense",prev_splicesite,prev_splicepos,floor));
  if (floor <= max_mismatches_allowed) {
    if (left_limit_p) {
      if (plusp) {
	debug4(printf("Testing from 0 to %d\n",prev_splicepos));
	nmismatches = Genome_count_mismatches_limit(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						    prev_left,/*pos5*/0,/*pos3*/prev_splicepos,max_mismatches_allowed,
						    dibasep,cmetp,plusp);
      } else {
	debug4(printf("Testing from 0 to %d\n",querylength-prev_splicepos));
	nmismatches = Genome_count_mismatches_limit(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						    prev_left,/*pos5*/0,/*pos3*/querylength-prev_splicepos,max_mismatches_allowed,
						    dibasep,cmetp,plusp);
      }
    } else {
      if (plusp) {
	debug4(printf("Testing from %d to %d\n",prev_splicepos,querylength));
	nmismatches = Genome_count_mismatches_limit(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						    prev_left,/*pos5*/prev_splicepos,/*pos3*/querylength,max_mismatches_allowed,
						    dibasep,cmetp,plusp);
      } else {
	debug4(printf("Testing from %d to %d\n",querylength-prev_splicepos,querylength));
	nmismatches = Genome_count_mismatches_limit(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						    prev_left,/*pos5*/querylength-prev_splicepos,/*pos3*/querylength,max_mismatches_allowed,
						    dibasep,cmetp,plusp);
      }
    }
    debug4(Genome_fill_buffer_simple(genome,prev_left,querylength,gbuffer));
    debug4(
	   if (plusp) {
	     printf("q: %s\ng: %s\n",queryuc_ptr,gbuffer);
	   } else {
	     printf("q.rc: %s\ng:    %s\n",queryrc,gbuffer);
	   });
    debug4(printf("Test: %s, %s, %s: splicesite %u => splicepos %d (%d mismatches)\n",
		  donorp ? "donor" : "acceptor",plusp ? "plus" : "minus",
		  sensep ? "sense" : "antisense",prev_splicesite,prev_splicepos,nmismatches));
    if (nmismatches <= max_mismatches_allowed) {
      if (prev_left + querylength > chrhigh) {
	/* update chromosome bounds, based on low end */
	chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,prev_left,prev_left);
	IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
	chrhigh += 1U;
      }
      Genome_fill_buffer_simple(genome,prev_left,querylength,gbuffer);
      if (donorp) {
	hit = Substring_new_donor(prev_splicepos,nmismatches,ncolordiffs,/*prob*/2.0,prev_left,
				  querylength,plusp,sensep,
				  gbuffer,queryuc_ptr,chrnum,chroffset,/*knownp*/true,
				  trim_ends_p,dibasep,cmetp);
      } else {
	hit = Substring_new_acceptor(prev_splicepos,nmismatches,ncolordiffs,/*prob*/2.0,prev_left,
				     querylength,plusp,sensep,
				     gbuffer,queryuc_ptr,chrnum,chroffset,/*knownp*/true,
				     trim_ends_p,dibasep,cmetp);
      }
      (*spliceends)[nmismatches] = List_push((*spliceends)[nmismatches],(void *) hit);
    }
  }

  FREE(array);
  return;
}



static void
identify_spliceends_known (List_T **donors_plus, List_T **antidonors_plus,
			   List_T **acceptors_plus, List_T **antiacceptors_plus,
			   List_T **donors_minus, List_T **antidonors_minus,
			   List_T **acceptors_minus, List_T **antiacceptors_minus,
			   T this, UINT4 *genome_blocks, UINT4 *snp_blocks,
			   Genome_T genome,
			   IIT_T chromosome_iit, char *queryuc_ptr, char *queryrc, int querylength,
			   int query_lastpos, char *gbuffer, int max_mismatches_allowed, 
			   Genomicpos_T *splicesites, Splicetype_T *splicetypes, int nsplicesites,
			   Compress_T query_compress_fwd, Compress_T query_compress_rev, 
			   Floors_T floors, bool dibasep, bool cmetp) {
  List_T splices[4], lastp0, lastp1, lastp2, lastp3, p;
  Splice_T splice;
  int querypos, nmismatches, i;

  /* Plus */
  splices[DONOR] = splices[ANTIDONOR] = splices[ACCEPTOR] = splices[ANTIACCEPTOR] = (List_T) NULL;
  lastp0 = lastp1 = lastp2 = lastp3 = (List_T) NULL;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->omitted[querypos] == false) {
      dual_search(splices,splicesites,nsplicesites,
		  this->plus_positions[querypos],this->plus_npositions[querypos],
		  splicetypes,querylength,/*diagterm*/querypos);
      for (p = splices[DONOR]; p != lastp0; p = p->rest) {
	splice = (Splice_T) p->first;
	splice->querypos = querypos;
	splice->splicepos = splice->splicesite - splice->left; /* for plus */
      }
      for (p = splices[ANTIDONOR]; p != lastp1; p = p->rest) {
	splice = (Splice_T) p->first;
	splice->querypos = querypos;
	splice->splicepos = splice->splicesite - splice->left; /* for plus */
      }
      for (p = splices[ACCEPTOR]; p != lastp2; p = p->rest) {
	splice = (Splice_T) p->first;
	splice->querypos = querypos;
	splice->splicepos = splice->splicesite - splice->left; /* for plus */
      }
      for (p = splices[ANTIACCEPTOR]; p != lastp3; p = p->rest) {
	splice = (Splice_T) p->first;
	splice->querypos = querypos;
	splice->splicepos = splice->splicesite - splice->left; /* for plus */
      }
      lastp0 = splices[DONOR];
      lastp1 = splices[ANTIDONOR];
      lastp2 = splices[ACCEPTOR];
      lastp3 = splices[ANTIACCEPTOR];
    }
  }

  /* End 1. donor, plus, sense */
  if (splices[DONOR]) {
    gather_spliceends_known(&(*donors_plus),splices[DONOR],genome_blocks,snp_blocks,genome,
			    chromosome_iit,queryuc_ptr,queryrc,querylength,query_lastpos,gbuffer,
			    max_mismatches_allowed,/*query_compress*/query_compress_fwd,
			    floors,dibasep,cmetp,/*donorp*/true,/*plusp*/true,/*sensep*/true,
			    /*left_floor_p*/true,/*left_limit_p*/true);
    List_free(&splices[DONOR]);
  }

  /* End 5. donor, plus, antisense */
  if (splices[ANTIDONOR]) {
    gather_spliceends_known(&(*antidonors_plus),splices[ANTIDONOR],genome_blocks,snp_blocks,genome,
			    chromosome_iit,queryuc_ptr,queryrc,querylength,query_lastpos,gbuffer,
			    max_mismatches_allowed,/*query_compress*/query_compress_fwd,
			    floors,dibasep,cmetp,/*donorp*/true,/*plusp*/true,/*sensep*/false,
			    /*left_floor_p*/false,/*left_limit_p*/false);
    List_free(&splices[ANTIDONOR]);
  }

  /* End 2. acceptor, plus, sense */
  if (splices[ACCEPTOR]) {
    gather_spliceends_known(&(*acceptors_plus),splices[ACCEPTOR],genome_blocks,snp_blocks,genome,
			    chromosome_iit,queryuc_ptr,queryrc,querylength,query_lastpos,gbuffer,
			    max_mismatches_allowed,/*query_compress*/query_compress_fwd,
			    floors,dibasep,cmetp,/*donorp*/false,/*plusp*/true,/*sensep*/true,
			    /*left_floor_p*/false,/*left_limit_p*/false);
    List_free(&splices[ACCEPTOR]);
  }

  /* End 6. acceptor, plus, antisense */
  if (splices[ANTIACCEPTOR]) {
    gather_spliceends_known(&(*antiacceptors_plus),splices[ANTIACCEPTOR],genome_blocks,snp_blocks,genome,
			    chromosome_iit,queryuc_ptr,queryrc,querylength,query_lastpos,gbuffer,
			    max_mismatches_allowed,/*query_compress*/query_compress_fwd,
			    floors,dibasep,cmetp,/*donorp*/false,/*plusp*/true,/*sensep*/false,
			    /*left_floor_p*/true,/*left_limit_p*/true);
    List_free(&splices[ANTIACCEPTOR]);
  }


  /* Minus */
  splices[DONOR] = splices[ANTIDONOR] = splices[ACCEPTOR] = splices[ANTIACCEPTOR] = (List_T) NULL;
  lastp0 = lastp1 = lastp2 = lastp3 = (List_T) NULL;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->omitted[querypos] == false) {
      dual_search(splices,splicesites,nsplicesites,
		  this->minus_positions[querypos],this->minus_npositions[querypos],
		  splicetypes,querylength,/*diagterm*/query_lastpos-querypos);
      for (p = splices[DONOR]; p != lastp0; p = p->rest) {
	splice = (Splice_T) p->first;
	splice->querypos = querypos;
	splice->splicepos = querylength - (splice->splicesite - splice->left); /* for minus */
      }
      for (p = splices[ANTIDONOR]; p != lastp1; p = p->rest) {
	splice = (Splice_T) p->first;
	splice->querypos = querypos;
	splice->splicepos = querylength - (splice->splicesite - splice->left); /* for minus */
      }
      for (p = splices[ACCEPTOR]; p != lastp2; p = p->rest) {
	splice = (Splice_T) p->first;
	splice->querypos = querypos;
	splice->splicepos = querylength - (splice->splicesite - splice->left); /* for minus */
      }
      for (p = splices[ANTIACCEPTOR]; p != lastp3; p = p->rest) {
	splice = (Splice_T) p->first;
	splice->querypos = querypos;
	splice->splicepos = querylength - (splice->splicesite - splice->left); /* for minus */
      }
      lastp0 = splices[DONOR];
      lastp1 = splices[ANTIDONOR];
      lastp2 = splices[ACCEPTOR];
      lastp3 = splices[ANTIACCEPTOR];
    }
  }

  /* End 7. donor, minus, antisense */
  if (splices[DONOR]) {
    gather_spliceends_known(&(*antidonors_minus),splices[DONOR],genome_blocks,snp_blocks,genome,
			    chromosome_iit,queryuc_ptr,queryrc,querylength,query_lastpos,gbuffer,
			    max_mismatches_allowed,/*query_compress*/query_compress_rev,
			    floors,dibasep,cmetp,/*donorp*/true,/*plusp*/false,/*sensep*/false,
			    /*left_floor_p*/false,/*left_limit_p*/true);
    List_free(&splices[DONOR]);
  }


  /* End 3. donor, minus, sense */
  if (splices[ANTIDONOR]) {
    gather_spliceends_known(&(*donors_minus),splices[ANTIDONOR],genome_blocks,snp_blocks,genome,
			    chromosome_iit,queryuc_ptr,queryrc,querylength,query_lastpos,gbuffer,
			    max_mismatches_allowed,/*query_compress*/query_compress_rev,
			    floors,dibasep,cmetp,/*donorp*/true,/*plusp*/false,/*sensep*/true,
			    /*left_floor_p*/true,/*left_limit_p*/false);
    List_free(&splices[ANTIDONOR]);
  }


  /* End 8. acceptor, minus, antisense */
  if (splices[ACCEPTOR]) {
    gather_spliceends_known(&(*antiacceptors_minus),splices[ACCEPTOR],genome_blocks,snp_blocks,genome,
			    chromosome_iit,queryuc_ptr,queryrc,querylength,query_lastpos,gbuffer,
			    max_mismatches_allowed,/*query_compress*/query_compress_rev,
			    floors,dibasep,cmetp,/*donorp*/false,/*plusp*/false,/*sensep*/false,
			    /*left_floor_p*/true,/*left_limit_p*/false);
    List_free(&splices[ACCEPTOR]);
  }


  /* End 4. acceptor, minus, sense */
  if (splices[ANTIACCEPTOR]) {
    gather_spliceends_known(&(*acceptors_minus),splices[ANTIACCEPTOR],genome_blocks,snp_blocks,genome,
			    chromosome_iit,queryuc_ptr,queryrc,querylength,query_lastpos,gbuffer,
			    max_mismatches_allowed,/*query_compress*/query_compress_rev,
			    floors,dibasep,cmetp,/*donorp*/false,/*plusp*/false,/*sensep*/true,
			    /*left_floor_p*/false,/*left_limit_p*/true);
    List_free(&splices[ANTIACCEPTOR]);
  }


  return;
}

#endif /* Solving known splice sites using dual search */


/* Do not compare against true or false */
/* Loosest criterion */
static int
sufficient_splice_prob_local (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support < 14) {
    return (spliceprob > 0.95);
  } else if (support < 20) {
    return (spliceprob > 0.90);
  } else if (support < 26) {
    return (spliceprob > 0.85);
  } else if (support < 32) {
    return (spliceprob > 0.80);
  } else if (support < 38) {
    return (spliceprob > 0.50);
  } else {
    return 1;
  }
}

#define MIN_SPLICE_SUPPORT_DISTANT 14

/* Do not compare against true or false */
/* Moderate criterion */
static int
sufficient_splice_prob_distant (int support, int nmismatches, double spliceprob) {
  /* printf("support %d, prob %.2f\n",support,spliceprob); */
  support -= 3*nmismatches;
  if (support < MIN_SPLICE_SUPPORT_DISTANT) {
    return 0;
  } else if (support < 20) {
    return (spliceprob > 0.95);
  } else if (support < 26) {
    return (spliceprob > 0.90);
  } else if (support < 32) {
    return (spliceprob > 0.85);
  } else if (support < 38) {
    return (spliceprob > 0.80);
  } else if (support < 44) {
    return (spliceprob > 0.50);
  } else {
    return 1;
  }
}

/* Do not compare against true or false */
/* Strictest criterion */
static int
sufficient_splice_prob_halfintron (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support < 14) {
    return 0;
  } else if (support < 20) {
    return 0;
  } else if (support < 26) {
    return (spliceprob > 0.95);
  } else if (support < 32) {
    return (spliceprob > 0.90);
  } else if (support < 38) {
    return (spliceprob > 0.85);
  } else if (support < 44) {
    return (spliceprob > 0.80);
  } else if (support < 50) {
    return (spliceprob > 0.50);
  } else {
    return 1;
  }
}


static List_T
solve_splicepair_local_plus (int *found_score, int *nhits, List_T hits,
			     Splicesegment_T segmenti, Splicesegment_T segmentj,
			     int *mismatch_positions_left, int nmismatches_left,
			     int *mismatch_positions_right, int nmismatches_right,
			     UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome,
			     char *queryuc_ptr, int querylength, Compress_T query_compress, 
			     bool *segmenti_donor_knownp, bool *segmentj_acceptor_knownp,
			     bool *segmentj_donor_knownp, bool *segmenti_acceptor_knownp,
			     bool novelsplicingp, int splicing_penalty,
			     int min_localsplicing_end_matches, int max_mismatches_allowed,
			     bool trim_ends_p, bool dibasep, bool cmetp) {
  Substring_T donor, acceptor;
  Genomicpos_T segmenti_left, segmentj_left, model_left;
  int best_splice_pos, splice_pos_start, splice_pos_end, splice_pos, model_length, i;
  char gbuffer[MAX_QUERYLENGTH+1];
  char gbuffer1[MAX_QUERYLENGTH+MAXENT_MAXLENGTH+1], gbuffer2[MAX_QUERYLENGTH+MAXENT_MAXLENGTH+1];
  double segmenti_donor_prob[MAX_QUERYLENGTH+1], segmentj_acceptor_prob[MAX_QUERYLENGTH+1],
    segmentj_donor_prob[MAX_QUERYLENGTH+1], segmenti_acceptor_prob[MAX_QUERYLENGTH+1];
  int best_donor_nmismatches, best_acceptor_nmismatches, donor_nmismatches, acceptor_nmismatches, ncolordiffs;
  int donor_support, acceptor_support;
  int sum, lefti, righti;
  double best_prob, prob;
  bool sensep, novelcheckp;

  segmenti_left = segmenti->diagonal - querylength;
  segmentj_left = segmentj->diagonal - querylength;
  debug4p(printf("solve_splicepair_local_plus: Getting genome at lefti %u and leftj %u\n",
		 segmenti_left,segmentj_left));
  
  splice_pos_start = querylength;
  splice_pos_end = 0;
  for (sum = 0; sum <= max_mismatches_allowed; sum++) {
    for (lefti = 0; lefti <= sum && lefti < nmismatches_left; lefti++) {
      if ((righti = sum - lefti) < nmismatches_right &&
	  mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	debug4p(printf("At %d+%d mismatches, splice_pos using right: %d\n",lefti,righti,mismatch_positions_right[righti]+1));
	debug4p(printf("At %d+%d mismatches, splice_pos using left: %d\n",lefti,righti,mismatch_positions_left[lefti]));
	if (mismatch_positions_right[righti] + 1 < splice_pos_start) {
	  splice_pos_start = mismatch_positions_right[righti] + 1;	/* This is leftmost position in righti+1 .. lefti */
	}
	if (mismatch_positions_left[lefti] > splice_pos_end) {
	  splice_pos_end = mismatch_positions_left[lefti];	/* This is rightmost position in righti+1 .. lefti */
	}
      }
    }
  }

  /* Exclude ends */
  if (splice_pos_start < min_localsplicing_end_matches) {
    splice_pos_start = min_localsplicing_end_matches;
  }
  if (splice_pos_end > querylength - min_localsplicing_end_matches) {
    splice_pos_end = querylength - min_localsplicing_end_matches;
  }

  if (splice_pos_start <= splice_pos_end) {
    /* End 1 to End 2: sense, originally from plus strand.  No complement. */
    novelcheckp = false;
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      model_left = segmenti_left + splice_pos_start - DONOR_MODEL_LEFT_MARGIN;
      model_length = DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+(splice_pos_end-splice_pos_start)+1;
      Genome_fill_buffer_simple(genome,model_left,model_length,gbuffer1);
      novelcheckp = true;
    }
    for (splice_pos = splice_pos_start, i = 0; splice_pos <= splice_pos_end; splice_pos++, i++) {
      if (segmenti_donor_knownp[splice_pos]) {
	debug4p(printf("plus sense donor already known at %u+%d\n",segmenti_left,splice_pos));
	segmenti_donor_prob[splice_pos] = 1.0;
      } else if (novelcheckp == false) {
	segmenti_donor_prob[splice_pos] = -1.0;
      } else {
	debug4p(printf("%s\n",&(gbuffer1[i])));
	debug4p(printf("%*s^^\n",DONOR_MODEL_LEFT_MARGIN,""));
	debug4p(printf("%*s%c%c\n",DONOR_MODEL_LEFT_MARGIN,"",
		       gbuffer1[i+DONOR_MODEL_LEFT_MARGIN],gbuffer1[i+DONOR_MODEL_LEFT_MARGIN+1]));
	if (gbuffer1[i+DONOR_MODEL_LEFT_MARGIN] == 'G' && 
	    (gbuffer1[i+DONOR_MODEL_LEFT_MARGIN+1] == 'T' ||
	     gbuffer1[i+DONOR_MODEL_LEFT_MARGIN+1] == 'C')) {
	  debug4p(printf("plus sense donor at %u: %f\n",
			 segmenti_left+splice_pos,Maxent_donor_prob(&(gbuffer1[i]))));
	  segmenti_donor_prob[splice_pos] = Maxent_donor_prob(&(gbuffer1[i]));
	} else {
	  segmenti_donor_prob[splice_pos] = -1.0;
	}
      }
    }


    novelcheckp = false;
    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      model_left = segmentj_left + splice_pos_start - ACCEPTOR_MODEL_LEFT_MARGIN;
      model_length = ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+(splice_pos_end-splice_pos_start)+1;
      Genome_fill_buffer_simple(genome,model_left,model_length,gbuffer1);
      novelcheckp = true;
    }
    for (splice_pos = splice_pos_start, i = 0; splice_pos <= splice_pos_end; splice_pos++, i++) {
      if (segmentj_acceptor_knownp[splice_pos]) {
	debug4p(printf("plus sense acceptor already known at %u+%d\n",segmentj_left,splice_pos));
	segmentj_acceptor_prob[splice_pos] = 1.0;
      } else if (novelcheckp == false) {
	segmentj_acceptor_prob[splice_pos] = -1.0;
      } else {
	debug4p(printf("%s\n",&(gbuffer1[i])));
	debug4p(printf("%*s^^\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,""));
	debug4p(printf("%*s%c%c\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,"",
		       gbuffer1[i+ACCEPTOR_MODEL_LEFT_MARGIN-2],gbuffer1[i+ACCEPTOR_MODEL_LEFT_MARGIN-1]));
	if (gbuffer1[i+ACCEPTOR_MODEL_LEFT_MARGIN-2] == 'A' &&  gbuffer1[i+ACCEPTOR_MODEL_LEFT_MARGIN-1] == 'G') {
	  debug4p(printf("plus sense acceptor at %u: %f\n",
			 segmentj_left+splice_pos,Maxent_acceptor_prob(&(gbuffer1[i]))));
	  segmentj_acceptor_prob[splice_pos] = Maxent_acceptor_prob(&(gbuffer1[i]));
	} else {
	  segmentj_acceptor_prob[splice_pos] = -1.0;
	}
      }
    }
    

    /* End 7 to End 8: antisense, originally from minus strand.  Complement. */
    novelcheckp = false;
    if (novelsplicingp && segmentj_left + splice_pos_start > DONOR_MODEL_RIGHT_MARGIN) {
      model_left = segmentj_left + splice_pos_start - DONOR_MODEL_RIGHT_MARGIN - 1;
      model_length = DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+(splice_pos_end-splice_pos_start)+1;
      Genome_fill_buffer_simple(genome,model_left,model_length,gbuffer1);
      make_complement_buffered(gbuffer2,gbuffer1,model_length);
      novelcheckp = true;
    }
    for (splice_pos = splice_pos_start, i = (splice_pos_end-splice_pos_start); splice_pos <= splice_pos_end; splice_pos++, i--) {
      if (segmentj_donor_knownp[splice_pos]) {
	debug4p(printf("plus antisense donor already known at %u+%d\n",segmentj_left,splice_pos));
	segmentj_donor_prob[splice_pos] = 1.0;
      } else if (novelcheckp == false) {
	segmentj_donor_prob[splice_pos] = -1.0;
      } else {
	debug4p(printf("%s\n",&(gbuffer2[i])));
	debug4p(printf("%*s^^\n",DONOR_MODEL_LEFT_MARGIN,""));
	debug4p(printf("%*s%c%c\n",DONOR_MODEL_LEFT_MARGIN,"",
		       gbuffer2[i+DONOR_MODEL_LEFT_MARGIN],gbuffer2[i+DONOR_MODEL_LEFT_MARGIN+1]));
	if (gbuffer2[i+DONOR_MODEL_LEFT_MARGIN] == 'G' && 
	    (gbuffer2[i+DONOR_MODEL_LEFT_MARGIN+1] == 'T' ||
	     gbuffer2[i+DONOR_MODEL_LEFT_MARGIN+1] == 'C')) {
	  debug4p(printf("plus antisense donor at %u: %f\n",
			 segmentj_left+splice_pos,Maxent_donor_prob(&(gbuffer2[i]))));
	  segmentj_donor_prob[splice_pos] = Maxent_donor_prob(&(gbuffer2[i]));
	} else {
	  segmentj_donor_prob[splice_pos] = -1.0;
	}
      }
    }
    

    novelcheckp = false;
    if (novelsplicingp && segmenti_left + splice_pos_start > ACCEPTOR_MODEL_RIGHT_MARGIN) {
      model_left = segmenti_left + splice_pos_start - ACCEPTOR_MODEL_RIGHT_MARGIN - 1;
      model_length = ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+(splice_pos_end-splice_pos_start)+1;
      Genome_fill_buffer_simple(genome,model_left,model_length,gbuffer1);
      make_complement_buffered(gbuffer2,gbuffer1,model_length);
      novelcheckp = true;
    }
    for (splice_pos = splice_pos_start, i = (splice_pos_end-splice_pos_start); splice_pos <= splice_pos_end; splice_pos++, i--) {
      if (segmenti_acceptor_knownp[splice_pos]) {
	debug4p(printf("plus antisense acceptor already known at %u+%d\n",segmenti_left,splice_pos));
	segmenti_acceptor_prob[splice_pos] = 1.0;
      } else if (novelcheckp == false) {
	segmenti_acceptor_prob[splice_pos] = -1.0;
      } else {
	debug4p(printf("%s\n",&(gbuffer2[i])));
	debug4p(printf("%*s^^\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,""));
	debug4p(printf("%*s%c%c\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,"",
		       gbuffer2[i+ACCEPTOR_MODEL_LEFT_MARGIN-2],gbuffer2[i+ACCEPTOR_MODEL_LEFT_MARGIN-1]));
	if (gbuffer2[i+ACCEPTOR_MODEL_LEFT_MARGIN-2] == 'A' && gbuffer2[i+ACCEPTOR_MODEL_LEFT_MARGIN-1] == 'G') {
	  debug4p(printf("plus antisense acceptor at %u: %f\n",
			 segmenti_left+splice_pos,Maxent_acceptor_prob(&(gbuffer2[i]))));
	  segmenti_acceptor_prob[splice_pos] = Maxent_acceptor_prob(&(gbuffer2[i]));
	} else {
	  segmenti_acceptor_prob[splice_pos] = -1.0;
	}
      }
    }
    
    best_prob = 0.0;
    sensep = true;

    /* End 1 to End 2: sense, originally from plus strand.  No complement. */
    /* Donor support is chimera_pos = splice_pos.  Acceptor support is querylength - chimera_pos = querylength - splice_pos. */
    for (splice_pos = splice_pos_start; splice_pos <= splice_pos_end; splice_pos++) {
      debug4p(printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",
		     splice_pos,segmenti_donor_prob[splice_pos],segmentj_acceptor_prob[splice_pos]));
      donor_support = splice_pos;
      acceptor_support = querylength - splice_pos;
      if (sufficient_splice_prob_local(donor_support,/*nmismatches*/0,segmenti_donor_prob[splice_pos]) &&
	  sufficient_splice_prob_local(acceptor_support,/*nmismatches*/0,segmentj_acceptor_prob[splice_pos]) &&
	  (prob = segmenti_donor_prob[splice_pos] + segmentj_acceptor_prob[splice_pos]) > best_prob) {
	/* Recheck support including mismatches */
	Genome_fill_buffer_simple(genome,/*left*/segmenti_left,querylength,gbuffer);
	donor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
							      dibasep,cmetp,/*plusp*/true);
	if (sufficient_splice_prob_local(donor_support,donor_nmismatches,segmenti_donor_prob[splice_pos])) {
	  Genome_fill_buffer_simple(genome,/*left*/segmentj_left,querylength,gbuffer);
	  acceptor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
								   /*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   dibasep,cmetp,/*plusp*/true);
	  if (sufficient_splice_prob_local(acceptor_support,acceptor_nmismatches,segmentj_acceptor_prob[splice_pos])) {
	    /* Success */
	    best_prob = prob;
	    best_splice_pos = splice_pos;
	    best_donor_nmismatches = donor_nmismatches;
	    best_acceptor_nmismatches = acceptor_nmismatches;
	  }
	}
      }
    }
    
    /* End 7 to End 8: antisense, originally from minus strand.  Complement. */
    /* Donor support is querylength - chimera_pos = querylength - splice_pos.  Acceptor support is chimera_pos = splice_pos. */
    for (splice_pos = splice_pos_start; splice_pos <= splice_pos_end; splice_pos++) {
      debug4p(printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",
		     splice_pos,segmentj_donor_prob[splice_pos],segmenti_acceptor_prob[splice_pos]));
      donor_support = querylength - splice_pos;
      acceptor_support = splice_pos;
      if (sufficient_splice_prob_local(donor_support,/*nmismatches*/0,segmentj_donor_prob[splice_pos]) &&
	  sufficient_splice_prob_local(acceptor_support,/*nmismatches*/0,segmenti_acceptor_prob[splice_pos]) &&
	  (prob = segmentj_donor_prob[splice_pos] + segmenti_acceptor_prob[splice_pos]) > best_prob) {
	/* Recheck support including mismatches */
	Genome_fill_buffer_simple(genome,/*left*/segmentj_left,querylength,gbuffer);
	donor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
							      dibasep,cmetp,/*plusp*/true);
	if (sufficient_splice_prob_local(donor_support,donor_nmismatches,segmentj_donor_prob[splice_pos])) {
	  Genome_fill_buffer_simple(genome,/*left*/segmenti_left,querylength,gbuffer);
	  acceptor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
								   /*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   dibasep,cmetp,/*plusp*/true);
	  if (sufficient_splice_prob_local(acceptor_support,acceptor_nmismatches,segmenti_acceptor_prob[splice_pos])) {
	    /* Success */
	    best_prob = prob;
	    best_splice_pos = splice_pos;
	    best_donor_nmismatches = donor_nmismatches;
	    best_acceptor_nmismatches = acceptor_nmismatches;
	    sensep = false;
	  }
	}
      }
    }
    
    if (best_prob > 0.0) {
      debug4p(printf("best_prob = %f at splice_pos %d\n",best_prob,best_splice_pos));
      if (sensep) {
	/* End 1 to End 2: sense, originally from plus strand.  No complement. */
	Genome_fill_buffer_simple(genome,/*left*/segmenti_left,querylength,gbuffer);
#if 0
	donor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segmenti_left,/*pos5*/0,/*pos3*/best_splice_pos,
							      dibasep,cmetp,/*plusp*/true);
#endif
	donor = Substring_new_donor(best_splice_pos,best_donor_nmismatches,ncolordiffs,
				    /*prob*/segmenti_donor_prob[best_splice_pos],
				    /*left*/segmenti_left,querylength,/*plusp*/true,/*sensep*/true,
				    gbuffer,queryuc_ptr,segmenti->chrnum,segmenti->chroffset,
				    /*knownp*/segmenti_donor_knownp[best_splice_pos],
				    trim_ends_p,dibasep,cmetp);
      
	Genome_fill_buffer_simple(genome,/*left*/segmentj_left,querylength,gbuffer);
#if 0
	acceptor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
								 /*left*/segmentj_left,/*pos5*/best_splice_pos,/*pos3*/querylength,
								 dibasep,cmetp,/*plusp*/true);
#endif
	acceptor = Substring_new_acceptor(best_splice_pos,best_acceptor_nmismatches,ncolordiffs,
					  /*prob*/segmentj_acceptor_prob[best_splice_pos],
					  /*left*/segmentj_left,querylength,/*plusp*/true,/*sensep*/true,
					  gbuffer,queryuc_ptr,segmentj->chrnum,segmentj->chroffset,
					  /*knownp*/segmentj_acceptor_knownp[best_splice_pos],
					  trim_ends_p,dibasep,cmetp);
	*nhits += 1;
	return List_push(hits,(void *) Stage3_new_splice(&(*found_score),donor,acceptor,/*distance*/segmentj_left - segmenti_left,
							 /*shortdistancep*/true,splicing_penalty,/*support*/querylength,querylength,
							 /*copyp*/false,/*sensep*/true));
      
      } else {
	/* End 7 to End 8: antisense, originally from minus strand.  Complement. */
	Genome_fill_buffer_simple(genome,/*left*/segmentj_left,querylength,gbuffer);
#if 0
	donor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segmentj_left,/*pos5*/best_splice_pos,/*pos3*/querylength,
							      dibasep,cmetp,/*plusp*/true);
#endif
	donor = Substring_new_donor(best_splice_pos,best_donor_nmismatches,ncolordiffs,
				    /*prob*/segmentj_donor_prob[best_splice_pos],
				    /*left*/segmentj_left,querylength,/*plusp*/true,/*sensep*/false,
				    gbuffer,queryuc_ptr,segmentj->chrnum,segmentj->chroffset,
				    /*knownp*/segmentj_donor_knownp[best_splice_pos],
				    trim_ends_p,dibasep,cmetp);
      
	Genome_fill_buffer_simple(genome,/*left*/segmenti_left,querylength,gbuffer);
#if 0
	acceptor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
								 /*left*/segmenti_left,/*pos5*/0,/*pos3*/best_splice_pos,
								 dibasep,cmetp,/*plusp*/true);
#endif
	acceptor = Substring_new_acceptor(best_splice_pos,best_acceptor_nmismatches,ncolordiffs,
					  /*prob*/segmenti_acceptor_prob[best_splice_pos],
					  /*left*/segmenti_left,querylength,/*plusp*/true,/*sensep*/false,
					  gbuffer,queryuc_ptr,segmenti->chrnum,segmenti->chroffset,
					  /*knownp*/segmenti_acceptor_knownp[best_splice_pos],
					  trim_ends_p,dibasep,cmetp);
	*nhits += 1;
	return List_push(hits,(void *) Stage3_new_splice(&(*found_score),donor,acceptor,/*distance*/segmentj_left - segmenti_left,
							 /*shortdistancep*/true,splicing_penalty,/*support*/querylength,querylength,
							 /*copyp*/false,/*sensep*/false));
      }
    }
  }

  return hits;
}


static List_T
find_splicepairs_local_plus (int *found_score, List_T hits, struct Splicesegment_T *plus_segments, int plus_nsegments,
			     UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome, char *queryuc_ptr, 
			     Floors_T floors, int querylength, int query_lastpos, Compress_T query_compress,
			     Genomicpos_T *splicesites, Splicetype_T *splicetypes, int nsplicesites,
			     bool novelsplicingp, Genomicpos_T max_distance,
			     int splicing_penalty, int min_localsplicing_end_matches, int max_mismatches_allowed,
			     bool trim_ends_p, bool dibasep, bool cmetp) {
#ifdef DEBUG4S
  int i;
#endif
  int j, startj;
  Splicesegment_T segmenti, segmentj;
  Genomicpos_T segmenti_left, segmentj_left, first_segmentj_left;
  int mismatch_positions_left[MAX_QUERYLENGTH], mismatch_positions_right[MAX_QUERYLENGTH], leftmost, rightmost;
  int colordiffs_left[MAX_QUERYLENGTH], colordiffs_right[MAX_QUERYLENGTH];
  int nmismatches_left, nmismatches_right;
  bool segmenti_donor_knownp[MAX_QUERYLENGTH+1], segmentj_acceptor_knownp[MAX_QUERYLENGTH+1],
    segmentj_donor_knownp[MAX_QUERYLENGTH+1], segmenti_acceptor_knownp[MAX_QUERYLENGTH+1];
  bool firstp, advancep;
  int floor, middle, pos, prev;
  int *floors_from_neg3, *floors_to_pos3;
  int nhits, nattempts;


  debug(printf("*** Starting find_splicepairs_local_plus on %d segments with max_distance %u ***\n",
	       plus_nsegments,max_distance));
  debug(printf("Initially have %d hits\n",List_length(hits)));

  nhits = List_length(hits);

  if (plus_nsegments > 1) {
    floors_from_neg3 = floors->scorefrom[-3];
    floors_to_pos3 = floors->scoreto[query_lastpos+3];

    for (segmenti = plus_segments; segmenti < &(plus_segments[plus_nsegments-1]) && nhits < MAX_LOCALSPLICING_HITS; segmenti++) {
      advancep = false;
      firstp = true;

      nattempts = 0;
      for (segmentj = segmenti+1; segmentj < &(plus_segments[plus_nsegments]) &&
	     segmentj->diagonal <= segmenti->diagonal + max_distance &&
	     segmentj->chrnum == segmenti->chrnum &&
	     ++nattempts < MAX_LOCALSPLICING_ATTEMPTS; segmentj++) {
	debug4s(printf("plus local?  diagonal %u, querypos %d..%d => diagonal %u, querypos %d..%d => ",
		       segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
		       segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	/* i5 i3 j5 j3 */
	if (segmenti->querypos3 >= segmentj->querypos5) {
	  debug4s(printf("Bad querypos\n"));
	} else {
	  floor = floors_from_neg3[segmenti->querypos5] + floors_to_pos3[segmentj->querypos3]
	    /* floors->score[-3][segmenti->querypos5] + floors->score[segmentj->querypos3][query_lastpos+3] */;
	  if (floors->prev_omitted == NULL) {
	    if ((middle = FLOOR_MIDDLE(segmentj->querypos5 - segmenti->querypos3)) > 0) {
	      middle--;	/* for splice, which looks like a mismatch */
	    }
	    debug4s(printf("\nmiddle (no omission): %d\n",middle));
	    floor += middle;	
	  } else {
	    pos = segmentj->querypos5;
	    debug4s(printf("\nmiddle (omission):"));
	    while (pos > segmenti->querypos3) {
	      if ((prev = floors->prev_omitted[pos]) < segmenti->querypos3) {
		prev = segmenti->querypos3;
	      }
	      if ((middle = FLOOR_MIDDLE(pos - prev)) > 0) {
		middle--;	/* for splice, which looks like a mismatch */
	      }
	      floor += middle;
	      debug4s(printf("(%d..%d)+%d,",prev,pos,middle));
	      pos = prev;
	    }
	    debug4s(printf("\n"));
	  }
	  if (floor > max_mismatches_allowed) {
	    debug4s(printf("too many mismatches, floor = %d+middle+%d=%d > %d\n",
			   floors->scorefrom[-3][segmenti->querypos5],
			   floors->scorefrom[segmentj->querypos3][query_lastpos+3],
			   floor,max_mismatches_allowed));
	  } else {
	    debug4s(printf("okay, floor = %d+middle+%d=%d\n",
			   floors->scorefrom[-3][segmenti->querypos5],
			   floors->scorefrom[segmentj->querypos3][query_lastpos+3],
			   floor));

	    segmentj_left = segmentj->diagonal - querylength;

	    if (firstp == true) {
	      first_segmentj_left = segmentj_left;
	      segmenti_left = segmenti->diagonal - querylength;
	      nmismatches_left = Genome_mismatches_left(mismatch_positions_left,colordiffs_left,max_mismatches_allowed,
							queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							/*left*/segmenti_left,/*pos5*/0,/*pos3*/querylength,
							dibasep,cmetp,/*plusp*/true);
	      firstp = false;
	      debug4s(printf("%d mismatches on left at:",nmismatches_left);
		      for (i = 0; i <= nmismatches_left; i++) {
			printf(" %d",mismatch_positions_left[i]);
		      }
		      printf("\n"));
	    }
	    
	    nmismatches_right = Genome_mismatches_right(mismatch_positions_right,colordiffs_right,max_mismatches_allowed,
							queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							/*left*/segmentj_left,/*pos5*/0,/*pos3*/querylength,
							dibasep,cmetp,/*plusp*/true);
	    debug4s(printf("%d mismatches on right at:",nmismatches_right);
		    for (i = 0; i <= nmismatches_right; i++) {
		      printf(" %d",mismatch_positions_right[i]);
		    }
		    printf("\n"));
	    
	    debug4s(printf("Want leftmost %d > rightmost %d (and leftmost < %d-%d and rightmost >= %d ?)\n",
			   mismatch_positions_left[nmismatches_left-1],mismatch_positions_right[nmismatches_right-1],
			   querylength,min_localsplicing_end_matches,min_localsplicing_end_matches));

	    if (nmismatches_left > 0 && nmismatches_right > 0 && 
		(leftmost = mismatch_positions_left[nmismatches_left-1]) > (rightmost = mismatch_positions_right[nmismatches_right-1])
#if 0
		&& leftmost <= querylength - min_localsplicing_end_matches && rightmost >= min_localsplicing_end_matches
#endif
		) {
	      
	      if (advancep == false) {
		/* Advance through known splice sites */
		if (nsplicesites > 0 && *splicesites < segmenti_left) {
		  j = 1;
		  while (j < nsplicesites && splicesites[j] < segmenti_left) {
		    j <<= 1;		/* gallop by 2 */
		  }
		  if (j >= nsplicesites) {
		    j = binary_search(j >> 1,nsplicesites,splicesites,segmenti_left);
		  } else {
		    j = binary_search(j >> 1,j,splicesites,segmenti_left);
		  }
		  splicesites += j;
		  splicetypes += j;
		  nsplicesites -= j;
		}
		
		/* Ends 1 and 8: set */
		memset(segmenti_donor_knownp,0,querylength*sizeof(bool));
		memset(segmenti_acceptor_knownp,0,querylength*sizeof(bool));
		
		j = 0;
		while (j < nsplicesites && splicesites[j] < segmenti->diagonal) {
		  if (splicetypes[j] == DONOR) {
		    debug4s(printf("Setting known donor for segmenti at %u\n",splicesites[j]));
		    segmenti_donor_knownp[splicesites[j] - segmenti_left] = true;
		  } else if (splicetypes[j] == ANTIACCEPTOR) {
		    debug4s(printf("Setting known acceptor for segmenti at %u\n",splicesites[j]));
		    segmenti_acceptor_knownp[splicesites[j] - segmenti_left] = true;
		  }
		  j++;
		}
		
		/* Ends 2 and 7: advance through known splice sites */
		startj = 0;
		if (nsplicesites > 0 && *splicesites < first_segmentj_left) {
		  startj = 1;
		  while (startj < nsplicesites && splicesites[startj] < first_segmentj_left) {
		    startj <<= 1;	/* gallop by 2 */
		  }
		  if (startj >= nsplicesites) {
		    startj = binary_search(startj >> 1,nsplicesites,splicesites,first_segmentj_left);
		  } else {
		    startj = binary_search(startj >> 1,startj,splicesites,first_segmentj_left);
		  }
		}
		
		advancep = true;
	      }
	      
	      /* Ends 2 and 7: set */
	      memset(segmentj_acceptor_knownp,0,querylength*sizeof(bool));
	      memset(segmentj_donor_knownp,0,querylength*sizeof(bool));
	      
	      j = startj;
	      while (j < nsplicesites && splicesites[j] < segmentj_left) {
		j++;			/* Advance to this segmentj */
	      }
	      while (j < nsplicesites && splicesites[j] < segmentj->diagonal) {
		if (splicetypes[j] == ACCEPTOR) {
		  debug4s(printf("Setting known acceptor for segmentj at %u\n",splicesites[j]));
		  segmentj_acceptor_knownp[splicesites[j] - segmentj_left] = true;
		} else if (splicetypes[j] == ANTIDONOR) {
		  debug4s(printf("Setting known donor for segmentj at %u\n",splicesites[j]));
		  segmentj_donor_knownp[splicesites[j] - segmentj_left] = true;
		}
		j++;
	      }
	      
	      debug4s(printf("success: solve_splicepair_local_plus\n"));
	      hits = solve_splicepair_local_plus(&(*found_score),&nhits,hits,segmenti,segmentj,
						 mismatch_positions_left,nmismatches_left,
						 mismatch_positions_right,nmismatches_right,
						 genome_blocks,snp_blocks,genome,
						 queryuc_ptr,querylength,query_compress,
						 segmenti_donor_knownp,segmentj_acceptor_knownp,
						 segmentj_donor_knownp,segmenti_acceptor_knownp,
						 novelsplicingp,
						 splicing_penalty,min_localsplicing_end_matches,max_mismatches_allowed,
						 trim_ends_p,dibasep,cmetp);
	    }
	  }
	}
      }
    }
  }

  debug(printf("Finished solve_splicepair_local_plus with %d hits\n",List_length(hits)));

  return hits;
}


static List_T
solve_splicepair_local_minus (int *found_score, int *nhits, List_T hits,
			      Splicesegment_T segmenti, Splicesegment_T segmentj,
			      int *mismatch_positions_left, int nmismatches_left,
			      int *mismatch_positions_right, int nmismatches_right,
			      UINT4 *genome_blocks, UINT4 *snp_blocks,
			      Genome_T genome, char *queryuc_ptr, int querylength, Compress_T query_compress,
			      bool *segmenti_donor_knownp, bool *segmentj_acceptor_knownp,
			      bool *segmentj_donor_knownp, bool *segmenti_acceptor_knownp,
			      bool novelsplicingp, int splicing_penalty, int min_localsplicing_end_matches, int max_mismatches_allowed,
			      bool trim_ends_p, bool dibasep, bool cmetp) {
  Substring_T donor, acceptor;
  Genomicpos_T segmenti_left, segmentj_left, model_left;
  int best_splice_pos, splice_pos_start, splice_pos_end, splice_pos, model_length, i;
  char gbuffer[MAX_QUERYLENGTH+1];
  char gbuffer1[MAX_QUERYLENGTH+MAXENT_MAXLENGTH+1], gbuffer2[MAX_QUERYLENGTH+MAXENT_MAXLENGTH+1];
  int best_donor_nmismatches, best_acceptor_nmismatches, donor_nmismatches, acceptor_nmismatches, ncolordiffs;
  int donor_support, acceptor_support;
  int sum, lefti, righti;
  double segmenti_donor_prob[MAX_QUERYLENGTH+1], segmentj_acceptor_prob[MAX_QUERYLENGTH+1],
    segmentj_donor_prob[MAX_QUERYLENGTH+1], segmenti_acceptor_prob[MAX_QUERYLENGTH+1];
  double best_prob, prob;
  bool sensep, novelcheckp;

  segmenti_left = segmenti->diagonal - querylength;
  segmentj_left = segmentj->diagonal - querylength;
  debug4p(printf("solve_splicepair_local_minus: Getting genome at lefti %u and leftj %u\n",
		 segmenti_left,segmentj_left));
  
  splice_pos_start = querylength;
  splice_pos_end = 0;
  for (sum = 0; sum <= max_mismatches_allowed; sum++) {
    for (lefti = 0; lefti <= sum && lefti < nmismatches_left; lefti++) {
      if ((righti = sum - lefti) < nmismatches_right &&
	  mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	debug4p(printf("At %d+%d mismatches, splice_pos using right: %d\n",lefti,righti,mismatch_positions_right[righti]+1));
	debug4p(printf("At %d+%d mismatches, splice_pos using left: %d\n",lefti,righti,mismatch_positions_left[lefti]));
	if (mismatch_positions_right[righti] + 1 < splice_pos_start) {
	  splice_pos_start = mismatch_positions_right[righti] + 1;	/* This is leftmost position in righti+1 .. lefti */
	}
	if (mismatch_positions_left[lefti] > splice_pos_end) {
	  splice_pos_end = mismatch_positions_left[lefti];	/* This is rightmost position in righti+1 .. lefti */
	}
      }
    }
  }
	
  /* Exclude ends */
  if (splice_pos_start < min_localsplicing_end_matches) {
    splice_pos_start = min_localsplicing_end_matches;
  }
  if (splice_pos_end > querylength - min_localsplicing_end_matches) {
    splice_pos_end = querylength - min_localsplicing_end_matches;
  }

  if (splice_pos_start <= splice_pos_end) {
    /* End 3 to End 4: sense, originally from minus strand.  Complement. */
    novelcheckp = false;
    if (novelsplicingp && segmentj_left + splice_pos_start > DONOR_MODEL_RIGHT_MARGIN) {
      model_left = segmentj_left + splice_pos_start - DONOR_MODEL_RIGHT_MARGIN - 1;
      model_length = DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+(splice_pos_end-splice_pos_start)+1;
      Genome_fill_buffer_simple(genome,model_left,model_length,gbuffer1);
      make_complement_buffered(gbuffer2,gbuffer1,model_length);
      novelcheckp = true;
    }
    for (splice_pos = splice_pos_start, i = (splice_pos_end-splice_pos_start); splice_pos <= splice_pos_end; splice_pos++, i--) {
      if (segmentj_donor_knownp[splice_pos]) {
	debug4p(printf("minus sense donor already known at %u+%d\n",segmentj_left,splice_pos));
	segmentj_donor_prob[splice_pos] = 1.0;
      } else if (novelcheckp == false) {
	segmentj_donor_prob[splice_pos] = -1.0;
      } else {
	debug4p(printf("%s\n",&(gbuffer2[i])));
	debug4p(printf("%*s^^\n",DONOR_MODEL_LEFT_MARGIN,""));
	debug4p(printf("%*s%c%c\n",DONOR_MODEL_LEFT_MARGIN,"",
		       gbuffer2[i+DONOR_MODEL_LEFT_MARGIN],gbuffer2[i+DONOR_MODEL_LEFT_MARGIN+1]));
	if (gbuffer2[i+DONOR_MODEL_LEFT_MARGIN] == 'G' && 
	    (gbuffer2[i+DONOR_MODEL_LEFT_MARGIN+1] == 'T' ||
	     gbuffer2[i+DONOR_MODEL_LEFT_MARGIN+1] == 'C')) {
	  debug4p(printf("minus sense donor at %u: %f\n",
			 segmentj_left+splice_pos,Maxent_donor_prob(&(gbuffer2[i]))));
	  segmentj_donor_prob[splice_pos] = Maxent_donor_prob(&(gbuffer2[i]));
	} else {
	  segmentj_donor_prob[splice_pos] = -1.0;
	}
      }
    }

    novelcheckp = false;
    if (novelsplicingp && segmenti_left + splice_pos_start > ACCEPTOR_MODEL_RIGHT_MARGIN) {
      model_left = segmenti_left + splice_pos_start - ACCEPTOR_MODEL_RIGHT_MARGIN - 1;
      model_length = ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+(splice_pos_end-splice_pos_start)+1;
      Genome_fill_buffer_simple(genome,model_left,model_length,gbuffer1);
      make_complement_buffered(gbuffer2,gbuffer1,model_length);
      novelcheckp = true;
    }
    for (splice_pos = splice_pos_start, i = (splice_pos_end-splice_pos_start); splice_pos <= splice_pos_end; splice_pos++, i--) {
      if (segmenti_acceptor_knownp[splice_pos]) {
	debug4p(printf("minus sense acceptor already known at %u+%d\n",segmenti_left,splice_pos));
	segmenti_acceptor_prob[splice_pos] = 1.0;
      } else if (novelcheckp == false) {
	segmenti_acceptor_prob[splice_pos] = -1.0;
      } else {
	debug4p(printf("%s\n",&(gbuffer2[i])));
	debug4p(printf("%*s^^\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,""));
	debug4p(printf("%*s%c%c\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,"",
		       gbuffer2[i+ACCEPTOR_MODEL_LEFT_MARGIN-2],gbuffer2[i+ACCEPTOR_MODEL_LEFT_MARGIN-1]));
	if (gbuffer2[i+ACCEPTOR_MODEL_LEFT_MARGIN-2] == 'A' && gbuffer2[i+ACCEPTOR_MODEL_LEFT_MARGIN-1] == 'G') {
	  debug4p(printf("minus sense acceptor at %u: %f\n",
			 segmenti_left+splice_pos,Maxent_acceptor_prob(&(gbuffer2[i]))));
	  segmenti_acceptor_prob[splice_pos] = Maxent_acceptor_prob(&(gbuffer2[i]));
	} else {
	  segmenti_acceptor_prob[splice_pos] = -1.0;
	}
      }
    }


    /* End 5 to End 6: antisense, originally from plus strand.  No complement. */
    novelcheckp = false;
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      model_left = segmenti_left + splice_pos_start - DONOR_MODEL_LEFT_MARGIN;
      model_length = DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+(splice_pos_end-splice_pos_start)+1;
      Genome_fill_buffer_simple(genome,model_left,model_length,gbuffer1);
      novelcheckp = true;
    }
    for (splice_pos = splice_pos_start, i = 0; splice_pos <= splice_pos_end; splice_pos++, i++) {
      if (segmenti_donor_knownp[splice_pos]) {
	debug4p(printf("minus antisense donor already known at %u+%d\n",segmenti_left,splice_pos));
	segmenti_donor_prob[splice_pos] = 1.0;
      } else if (novelcheckp == false) {
	segmenti_donor_prob[splice_pos] = -1.0;
      } else {
	debug4p(printf("%s\n",&(gbuffer1[i])));
	debug4p(printf("%*s^^\n",DONOR_MODEL_LEFT_MARGIN,""));
	debug4p(printf("%*s%c%c\n",DONOR_MODEL_LEFT_MARGIN,"",
		       gbuffer1[i+DONOR_MODEL_LEFT_MARGIN],gbuffer1[i+DONOR_MODEL_LEFT_MARGIN+1]));
	if (gbuffer1[i+DONOR_MODEL_LEFT_MARGIN] == 'G' && 
	    (gbuffer1[i+DONOR_MODEL_LEFT_MARGIN+1] == 'T' ||
	     gbuffer1[i+DONOR_MODEL_LEFT_MARGIN+1] == 'C')) {
	  debug4p(printf("minus antisense donor at %u: %f\n",
			 segmenti_left+splice_pos,Maxent_donor_prob(&(gbuffer1[i]))));
	  segmenti_donor_prob[splice_pos] = Maxent_donor_prob(&(gbuffer1[i]));
	} else {
	  segmenti_donor_prob[splice_pos] = -1.0;
	}
      }
    }

    novelcheckp = false;
    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      model_left = segmentj_left + splice_pos_start - ACCEPTOR_MODEL_LEFT_MARGIN;
      model_length = ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+(splice_pos_end-splice_pos_start)+1;
      Genome_fill_buffer_simple(genome,model_left,model_length,gbuffer1);
      novelcheckp = true;
    }
    for (splice_pos = splice_pos_start, i = 0; splice_pos <= splice_pos_end; splice_pos++, i++) {
      if (segmentj_acceptor_knownp[splice_pos]) {
	debug4p(printf("minus antisense acceptor already known at %u+%d\n",segmentj_left,splice_pos));
	segmentj_acceptor_prob[splice_pos] = 1.0;
      } else if (novelcheckp == false) {
	segmentj_acceptor_prob[splice_pos] = -1.0;
      } else {
	debug4p(printf("%s\n",&(gbuffer1[i])));
	debug4p(printf("%*s^^\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,""));
	debug4p(printf("%*s%c%c\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,"",
		       gbuffer1[i+ACCEPTOR_MODEL_LEFT_MARGIN-2],gbuffer1[i+ACCEPTOR_MODEL_LEFT_MARGIN-1]));
	if (gbuffer1[i+ACCEPTOR_MODEL_LEFT_MARGIN-2] == 'A' && gbuffer1[i+ACCEPTOR_MODEL_LEFT_MARGIN-1] == 'G') {
	  debug4p(printf("minus antisense acceptor at %u: %f\n",
			 segmentj_left+splice_pos+1,Maxent_acceptor_prob(&(gbuffer1[i]))));
	  segmentj_acceptor_prob[splice_pos] = Maxent_acceptor_prob(&(gbuffer1[i]));
	} else {
	  segmentj_acceptor_prob[splice_pos] = -1.0;
	}
      }
    }


    best_prob = 0.0;
    sensep = true;
    /* End 3 to End 4: sense, originally from minus strand.  Complement. */
    /* Donor support is chimera_pos = querylength - splice_pos.  Acceptor support is querylength - chimera_pos = querylength - (querylength - splice_pos) = splice_pos. */
    for (splice_pos = splice_pos_start; splice_pos <= splice_pos_end; splice_pos++) {
      debug4p(printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",
		     splice_pos,segmentj_donor_prob[splice_pos],segmenti_acceptor_prob[splice_pos]));
      donor_support = querylength - splice_pos;
      acceptor_support = splice_pos;
      if (sufficient_splice_prob_local(donor_support,/*nmismatches*/0,segmentj_donor_prob[splice_pos]) &&
	  sufficient_splice_prob_local(acceptor_support,/*nmismatches*/0,segmenti_acceptor_prob[splice_pos]) &&
	  (prob = segmentj_donor_prob[splice_pos] + segmenti_acceptor_prob[splice_pos]) > best_prob) {
	/* Recheck support including mismatches */
	Genome_fill_buffer_simple(genome,/*left*/segmentj_left,querylength,gbuffer);
	donor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
							      dibasep,cmetp,/*plusp*/false);
	if (sufficient_splice_prob_local(donor_support,donor_nmismatches,segmentj_donor_prob[splice_pos])) {
	  Genome_fill_buffer_simple(genome,/*left*/segmenti_left,querylength,gbuffer);
	  acceptor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
								   /*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   dibasep,cmetp,/*plusp*/false);
	  if (sufficient_splice_prob_local(acceptor_support,acceptor_nmismatches,segmenti_acceptor_prob[splice_pos])) {
	    /* Success */
	    best_prob = prob;
	    best_splice_pos = splice_pos;
	    best_donor_nmismatches = donor_nmismatches;
	    best_acceptor_nmismatches = acceptor_nmismatches;
	  }
	}
      }
    }

    /* End 5 to End 6: antisense, originally from plus strand.  No complement. */
    /* Donor support is querylength - chimera_pos = querylength - (querylength - splice_pos) = splice_pos.  Acceptor support is chimera_pos = querylength - splice_pos. */
    for (splice_pos = splice_pos_start; splice_pos <= splice_pos_end; splice_pos++) {
      debug4p(printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",
		     splice_pos,segmenti_donor_prob[splice_pos],segmentj_acceptor_prob[splice_pos]));
      donor_support = splice_pos;
      acceptor_support = querylength - splice_pos;
      if (sufficient_splice_prob_local(donor_support,/*nmismatches*/0,segmenti_donor_prob[splice_pos]) &&
	  sufficient_splice_prob_local(acceptor_support,/*nmismatches*/0,segmentj_acceptor_prob[splice_pos]) &&
	  (prob = segmenti_donor_prob[splice_pos] + segmentj_acceptor_prob[splice_pos]) > best_prob) {
	/* Recheck support including mismatches */
	Genome_fill_buffer_simple(genome,/*left*/segmenti_left,querylength,gbuffer);
	donor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
							      dibasep,cmetp,/*plusp*/false);
	if (sufficient_splice_prob_local(donor_support,donor_nmismatches,segmenti_donor_prob[splice_pos])) {
	  Genome_fill_buffer_simple(genome,/*left*/segmentj_left,querylength,gbuffer);
	  acceptor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
								   /*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   dibasep,cmetp,/*plusp*/false);
	  if (sufficient_splice_prob_local(acceptor_support,acceptor_nmismatches,segmentj_acceptor_prob[splice_pos])) {
	    /* Success */
	    best_prob = prob;
	    best_splice_pos = splice_pos;
	    best_donor_nmismatches = donor_nmismatches;
	    best_acceptor_nmismatches = acceptor_nmismatches;
	    sensep = false;
	  }
	}
      }
    }

    if (best_prob > 0.0) {
      debug4p(printf("best_prob = %f at splice_pos %d\n",best_prob,best_splice_pos));

      if (sensep) {
	/* End 3 to End 4: sense, originally from minus strand.  Complement. */
	Genome_fill_buffer_simple(genome,/*left*/segmentj_left,querylength,gbuffer);
#if 0
	donor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segmentj_left,/*pos5*/best_splice_pos,/*pos3*/querylength,
							      dibasep,cmetp,/*plusp*/false);
#endif
	donor = Substring_new_donor(/*donor_pos*/querylength-best_splice_pos,
				    best_donor_nmismatches,ncolordiffs,
				    /*prob*/segmentj_donor_prob[best_splice_pos],/*left*/segmentj_left,
				    querylength,/*plusp*/false,/*sensep*/true,
				    gbuffer,queryuc_ptr,segmentj->chrnum,segmentj->chroffset,
				    /*knownp*/segmentj_donor_knownp[best_splice_pos],
				    trim_ends_p,dibasep,cmetp);
	      
	Genome_fill_buffer_simple(genome,/*left*/segmenti_left,querylength,gbuffer);
#if 0
	acceptor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
								 /*left*/segmenti_left,/*pos5*/0,/*pos3*/best_splice_pos,
								 dibasep,cmetp,/*plusp*/false);
#endif
	acceptor = Substring_new_acceptor(/*acceptor_pos*/querylength-best_splice_pos,
					  best_acceptor_nmismatches,ncolordiffs,
					  /*prob*/segmenti_acceptor_prob[best_splice_pos],/*left*/segmenti_left,
					  querylength,/*plusp*/false,/*sensep*/true,
					  gbuffer,queryuc_ptr,segmenti->chrnum,segmenti->chroffset,
					  /*knownp*/segmenti_acceptor_knownp[best_splice_pos],
					  trim_ends_p,dibasep,cmetp);
	*nhits += 1;
	return List_push(hits,(void *) Stage3_new_splice(&(*found_score),donor,acceptor,/*distance*/segmentj_left - segmenti_left,
							 /*shortdistancep*/true,splicing_penalty,/*support*/querylength,querylength,
							 /*copyp*/false,/*sensep*/true));
      } else {
	/* End 5 to End 6: antisense, originally from plus strand.  No complement. */
	Genome_fill_buffer_simple(genome,/*left*/segmenti_left,querylength,gbuffer);
#if 0
	donor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segmenti_left,/*pos5*/0,/*pos3*/best_splice_pos,
							      dibasep,cmetp,/*plusp*/false);
#endif
	donor = Substring_new_donor(/*donor_pos*/querylength-best_splice_pos,
				    best_donor_nmismatches,ncolordiffs,
				    /*prob*/segmenti_donor_prob[best_splice_pos],/*left*/segmenti_left,
				    querylength,/*plusp*/false,/*sensep*/false,
				    gbuffer,queryuc_ptr,segmenti->chrnum,segmenti->chroffset,
				    /*knownp*/segmenti_donor_knownp[best_splice_pos],
				    trim_ends_p,dibasep,cmetp);
	
	Genome_fill_buffer_simple(genome,/*left*/segmentj_left,querylength,gbuffer);
#if 0
	acceptor_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
								 /*left*/segmentj_left,/*pos5*/best_splice_pos,/*pos3*/querylength,
								 dibasep,cmetp,/*plusp*/false);
#endif
	acceptor = Substring_new_acceptor(/*acceptor_pos*/querylength-best_splice_pos,
					  best_acceptor_nmismatches,ncolordiffs,
					  /*prob*/segmentj_acceptor_prob[best_splice_pos],/*left*/segmentj_left,
					  querylength,/*plusp*/false,/*sensep*/false,
					  gbuffer,queryuc_ptr,segmentj->chrnum,segmentj->chroffset,
					  /*knownp*/segmentj_acceptor_knownp[best_splice_pos],
					  trim_ends_p,dibasep,cmetp);
	*nhits += 1;
	return List_push(hits,(void *) Stage3_new_splice(&(*found_score),donor,acceptor,/*distance*/segmentj_left - segmenti_left,
							 /*shortdistancep*/true,splicing_penalty,/*support*/querylength,querylength,
							 /*copyp*/false,/*sensep*/false));
      }
    }
  }

  return hits;
}

static List_T
find_splicepairs_local_minus (int *found_score, List_T hits, struct Splicesegment_T *minus_segments, int minus_nsegments,
			      UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome, char *queryuc_ptr, 
			      Floors_T floors, int querylength, int query_lastpos, Compress_T query_compress,
			      Genomicpos_T *splicesites, Splicetype_T *splicetypes, int nsplicesites,
			      bool novelsplicingp, Genomicpos_T max_distance,
			      int splicing_penalty, int min_localsplicing_end_matches, int max_mismatches_allowed,
			      bool trim_ends_p, bool dibasep, bool cmetp) {
#ifdef DEBUG4S
  int i;
#endif
  int j, startj;
  Splicesegment_T segmenti, segmentj;
  Genomicpos_T segmenti_left, segmentj_left, first_segmentj_left;
  int mismatch_positions_left[MAX_QUERYLENGTH], mismatch_positions_right[MAX_QUERYLENGTH], leftmost, rightmost;
  int colordiffs_left[MAX_QUERYLENGTH], colordiffs_right[MAX_QUERYLENGTH];
  int nmismatches_left, nmismatches_right;
  bool segmenti_donor_knownp[MAX_QUERYLENGTH+1], segmentj_acceptor_knownp[MAX_QUERYLENGTH+1],
    segmentj_donor_knownp[MAX_QUERYLENGTH+1], segmenti_acceptor_knownp[MAX_QUERYLENGTH+1];
  bool firstp, advancep;
  int floor, middle, pos, prev;
  int *floors_from_neg3, *floors_to_pos3;
  int nhits, nattempts;

  debug(printf("*** Starting find_splicepairs_local_minus on %d segments with max_distance %u ***\n",
	       minus_nsegments,max_distance));
  debug(printf("Initially have %d hits\n",List_length(hits)));

  nhits = List_length(hits);

  if (minus_nsegments > 1) {
    floors_from_neg3 = floors->scorefrom[-3];
    floors_to_pos3 = floors->scoreto[query_lastpos+3];

    for (segmenti = minus_segments; segmenti < &(minus_segments[minus_nsegments-1]) && nhits < MAX_LOCALSPLICING_HITS; segmenti++) {
      advancep = false;
      firstp = true;

      nattempts = 0;
      for (segmentj = segmenti+1; segmentj < &(minus_segments[minus_nsegments]) &&
	     segmentj->diagonal <= segmenti->diagonal + max_distance &&
	     segmentj->chrnum == segmenti->chrnum &&
	     ++nattempts < MAX_LOCALSPLICING_ATTEMPTS; segmentj++) {
	debug4s(printf("minus local?  diagonal %u, querypos %d..%d => diagonal %u, querypos %d..%d => ",
		       segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
		       segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	/* j5 j3 i5 i3 */
	if (segmentj->querypos3 >= segmenti->querypos5) {
	  debug4s(printf("Bad querypos\n"));
	} else {
	  floor = floors_from_neg3[segmentj->querypos5] + floors_to_pos3[segmenti->querypos3]
	    /* floors->score[-3][segmentj->querypos5] + floors->score[segmenti->querypos3][query_lastpos+3] */;
	  if (floors->prev_omitted == NULL) {
	    if ((middle = FLOOR_MIDDLE(segmenti->querypos5 - segmentj->querypos3 /*- indels*/)) > 0) {
	      middle--;	/* for deletion, which looks like a mismatch */
	    }
	    debug4s(printf("\nmiddle (no omission): %d\n",middle));
	    floor += middle;	
	  } else {
	    pos = segmenti->querypos5;
	    debug4s(printf("\nmiddle (omission):"));
	    while (pos > segmentj->querypos3) {
	      if ((prev = floors->prev_omitted[pos]) < segmentj->querypos3) {
		prev = segmentj->querypos3;
	      }
	      if ((middle = FLOOR_MIDDLE(pos - prev /*- indels*/)) > 0) {
		middle--; /* for deletion, which looks like a mismatch */
	      }
	      floor += middle;
	      debug4s(printf("(%d..%d)+%d,",prev,pos,middle));
	      pos = prev;
	    }
	    debug4s(printf("\n"));
	  }
	  if (floor > max_mismatches_allowed) {
	    debug4s(printf("too many mismatches, floor = %d+middle+%d=%d > %d\n",
			   floors->scorefrom[-3][segmentj->querypos5],
			   floors->scorefrom[segmenti->querypos3][query_lastpos+3],
			   floor,max_mismatches_allowed));
	  } else {
	    debug4s(printf("okay, floor = %d+middle+%d=%d\n",
			   floors->scorefrom[-3][segmentj->querypos5],
			   floors->scorefrom[segmenti->querypos3][query_lastpos+3],
			   floor));

	    segmentj_left = segmentj->diagonal - querylength;

	    if (firstp == true) {
	      first_segmentj_left = segmentj_left;
	      segmenti_left = segmenti->diagonal - querylength;
	      nmismatches_left = Genome_mismatches_left(mismatch_positions_left,colordiffs_left,max_mismatches_allowed,
							queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							/*left*/segmenti_left,/*pos5*/0,/*pos3*/querylength,
							dibasep,cmetp,/*plusp*/true);
	      firstp = false;
	      debug4s(printf("%d mismatches on left at:",nmismatches_left);
		      for (i = 0; i <= nmismatches_left; i++) {
			printf(" %d",mismatch_positions_left[i]);
		      }
		      printf("\n"));
	    }

	    nmismatches_right = Genome_mismatches_right(mismatch_positions_right,colordiffs_right,max_mismatches_allowed,
							queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							/*left*/segmentj_left,/*pos5*/0,/*pos3*/querylength,
							dibasep,cmetp,/*plusp*/true);
	    debug4s(printf("%d mismatches on right at:",nmismatches_right);
		    for (i = 0; i <= nmismatches_right; i++) {
		      printf(" %d",mismatch_positions_right[i]);
		    }
		    printf("\n"));

	    debug4s(printf("Want leftmost %d > rightmost %d (and leftmost <= %d-%d and rightmost >= %d ?)\n",
			   mismatch_positions_left[nmismatches_left-1],mismatch_positions_right[nmismatches_right-1],
			   querylength,min_localsplicing_end_matches,min_localsplicing_end_matches));

	    if (nmismatches_left > 0 && nmismatches_right > 0 && 
		(leftmost = mismatch_positions_left[nmismatches_left-1]) > (rightmost = mismatch_positions_right[nmismatches_right-1]) 
#if 0
		&& leftmost <= querylength - min_localsplicing_end_matches && rightmost >= min_localsplicing_end_matches
#endif
		) {
	      
	      if (advancep == false) {
		/* Advance through known splice sites */
		if (nsplicesites > 0 && *splicesites < segmenti_left) {
		  j = 1;
		  while (j < nsplicesites && splicesites[j] < segmenti_left) {
		    j <<= 1;		/* gallop by 2 */
		  }
		  if (j >= nsplicesites) {
		    j = binary_search(j >> 1,nsplicesites,splicesites,segmenti_left);
		  } else {
		    j = binary_search(j >> 1,j,splicesites,segmenti_left);
		  }
		  splicesites += j;
		  splicetypes += j;
		  nsplicesites -= j;
		}

		/* Ends 4 and 5: set */
		memset(segmenti_acceptor_knownp,0,querylength*sizeof(bool));
		memset(segmenti_donor_knownp,0,querylength*sizeof(bool));

		j = 0;
		while (j < nsplicesites && splicesites[j] < segmenti->diagonal) {
		  if (splicetypes[j] == ANTIACCEPTOR) {
		    debug4s(printf("Setting known acceptor for segmenti at %u\n",splicesites[j]));
		    segmenti_acceptor_knownp[splicesites[j] - segmenti_left] = true;
		  } else if (splicetypes[j] == DONOR) {
		    debug4s(printf("Setting known donor for segmenti at %u\n",splicesites[j]));
		    segmenti_donor_knownp[splicesites[j] - segmenti_left] = true;
		  }
		  j++;
		}
		
		/* Ends 3 and 6: advance through known splice sites */
		startj = 0;
		if (nsplicesites > 0 && *splicesites < first_segmentj_left) {
		  startj = 1;
		  while (startj < nsplicesites && splicesites[startj] < first_segmentj_left) {
		    startj <<= 1;	/* gallop by 2 */
		  }
		  if (startj >= nsplicesites) {
		    startj = binary_search(startj >> 1,nsplicesites,splicesites,first_segmentj_left);
		  } else {
		    startj = binary_search(startj >> 1,startj,splicesites,first_segmentj_left);
		  }
		}
		
		advancep = true;
	      }
	      
	      /* Ends 3 and 6: set */
	      memset(segmentj_donor_knownp,0,querylength*sizeof(bool));
	      memset(segmentj_acceptor_knownp,0,querylength*sizeof(bool));
	      
	      j = startj;
	      while (j < nsplicesites && splicesites[j] < segmentj_left) {
		j++;			/* Advance to this segmentj */
	      }
	      while (j < nsplicesites && splicesites[j] < segmentj->diagonal) {
		if (splicetypes[j] == ANTIDONOR) {
		  debug4s(printf("Setting known donor for segmentj at %u\n",splicesites[j]));
		  segmentj_donor_knownp[splicesites[j] - segmentj_left] = true;
		} else if (splicetypes[j] == ACCEPTOR) {
		  debug4s(printf("Setting known acceptor for segmentj at %u\n",splicesites[j]));
		  segmentj_acceptor_knownp[splicesites[j] - segmentj_left] = true;
		}
		j++;
	      }
	      
	      debug4s(printf("success: solve_splicepair_local_minus\n"));
	      hits = solve_splicepair_local_minus(&(*found_score),&nhits,hits,segmenti,segmentj,
						  mismatch_positions_left,nmismatches_left,
						  mismatch_positions_right,nmismatches_right,
						  genome_blocks,snp_blocks,genome,
						  queryuc_ptr,querylength,query_compress,
						  segmenti_donor_knownp,segmentj_acceptor_knownp,
						  segmentj_donor_knownp,segmenti_acceptor_knownp,
						  novelsplicingp,splicing_penalty,min_localsplicing_end_matches,max_mismatches_allowed,
						  trim_ends_p,dibasep,cmetp);
	    }
	  }
	}
      }
    }
  }

  debug(printf("Finished solve_splicepair_local_minus with %d hits\n",List_length(hits)));

  return hits;
}


static void
identify_spliceends_plus (List_T **donors_plus, List_T **antidonors_plus,
			  List_T **acceptors_plus, List_T **antiacceptors_plus,
			  struct Splicesegment_T *plus_segments, int plus_nsegments,
			  UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome,
			  char *queryuc_ptr, Floors_T floors, int querylength, int query_lastpos,
			  char *gbuffer, Compress_T query_compress,
			  Genomicpos_T *splicesites, Splicetype_T *splicetypes, int nsplicesites,
			  bool novelsplicingp, bool canonicalp, 
			  int max_mismatches_allowed, bool trim_ends_p, bool dibasep, bool cmetp) {
  Splicesegment_T segment;
  Substring_T hit;
  Genomicpos_T segment_left, model_left;
  int nmismatches, stopi, j;
  int splice_pos;
  char gbuffer1[MAXENT_MAXLENGTH+1], gbuffer2[MAXENT_MAXLENGTH+1];
  double prob;

  int trimpos;
  int mismatch_positions[MAX_QUERYLENGTH+1], colordiffs[MAX_QUERYLENGTH];
  int nmismatches_left, nmismatches_right;
  int ncolordiffs = 0;
  int *floors_from_neg3, *floors_to_pos3;
  bool advancep, found_known_p;

#ifdef DEBUG4E
  int i;
#endif


  floors_from_neg3 = floors->scorefrom[-3];
  floors_to_pos3 = floors->scoreto[query_lastpos+3];

  for (segment = plus_segments; segment < &(plus_segments[plus_nsegments]); segment++) {
    segment_left = segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
    debug4e(printf("identify_spliceends_plus: Checking up to %d mismatches at diagonal %u (querypos %d..%d) - querylength %d = %u\n",
		    max_mismatches_allowed,segment->diagonal,segment->querypos5,segment->querypos3,querylength,segment_left));

    debug4e(
	     Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	     printf("genome 0..: %s\n",gbuffer);
	     printf("query  0..: %s\n",queryuc_ptr);
	     );

    gbuffer[0] = '\0';		/* Signals that we haven't filled it yet */
    advancep = false;

    if (floors_from_neg3[segment->querypos5] <= max_mismatches_allowed) {
      trimpos = Genome_trim_right(queryuc_ptr,query_compress,genome_blocks,snp_blocks,
				  /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
				  dibasep,cmetp,/*plusp*/true);

      nmismatches_left = Genome_mismatches_left(mismatch_positions,colordiffs,max_mismatches_allowed,
						queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						/*left*/segment_left,/*pos5*/0,/*pos3*/trimpos,
						dibasep,cmetp,/*plusp*/true);
    
      debug4e(
	       printf("For %d up to trim %d, %d mismatches on left at:",0,trimpos,nmismatches_left);
	       for (i = 0; i <= nmismatches_left; i++) {
		 printf(" %d",mismatch_positions[i]);
	       }
	       printf("\n");
	       );
      
      if (nmismatches_left == 0) {
	stopi = trimpos;
      } else {
	stopi = mismatch_positions[nmismatches_left-1];
      }
      if (stopi > querylength - 2) {
	stopi = querylength - 2;
      }
      debug4e(printf("Search from %d up to %d\n",INDEX1PART,stopi));
      
      if (stopi >= MIN_SPLICE_SUPPORT_DISTANT) {
	/* Advance through known splice sites */
	if (/* advancep == false && */ nsplicesites > 0 && *splicesites < segment_left + MIN_SPLICE_SUPPORT_DISTANT) {
	  j = 1;
	  while (j < nsplicesites && splicesites[j] < segment_left + MIN_SPLICE_SUPPORT_DISTANT) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= nsplicesites) {
	    j = binary_search(j >> 1,nsplicesites,splicesites,segment_left + MIN_SPLICE_SUPPORT_DISTANT);
	  } else {
	    j = binary_search(j >> 1,j,splicesites,segment_left + MIN_SPLICE_SUPPORT_DISTANT);
	  }
	  splicesites += j;
	  splicetypes += j;
	  nsplicesites -= j;
	  advancep = true;
	}

	/* Known splicing */
	found_known_p = false;
	for (j = 0; j < nsplicesites && splicesites[j] <= segment_left + stopi; j++) {
	  if (splicetypes[j] == DONOR) {
	    splice_pos = splicesites[j] - segment_left;
	    nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							    /*left*/segment_left,/*pos5*/0,/*pos3*/splice_pos,
							    dibasep,cmetp,/*plusp*/true);
	    if (gbuffer[0] == '\0') {
	      Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	    }
	    debug4e(printf("Known donor for segmenti at %u (%d mismatches)\n",splicesites[j],nmismatches));
	    hit = Substring_new_donor(splice_pos,nmismatches,ncolordiffs,/*prob*/2.0,/*left*/segment_left,
				      querylength,/*plusp*/true,/*sensep*/true,
				      gbuffer,queryuc_ptr,segment->chrnum,segment->chroffset,/*knownp*/true,
				      trim_ends_p,dibasep,cmetp);
	    (*donors_plus)[nmismatches] = List_push((*donors_plus)[nmismatches],(void *) hit);
	    found_known_p = true;
	  } else if (splicetypes[j] == ANTIACCEPTOR) {
	    splice_pos = splicesites[j] - segment_left;
	    nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							    /*left*/segment_left,/*pos5*/0,/*pos3*/splice_pos,
							    dibasep,cmetp,/*plusp*/true);
	    if (gbuffer[0] == '\0') {
	      Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	    }
	    debug4e(printf("Known antiacceptor for segmenti at %u (%d mismatches)\n",splicesites[j],nmismatches));
	    hit = Substring_new_acceptor(splice_pos,nmismatches,ncolordiffs,/*prob*/2.0,/*left*/segment_left,
					 querylength,/*plusp*/true,/*sensep*/false,
					 gbuffer,queryuc_ptr,segment->chrnum,segment->chroffset,/*knownp*/true,
					 trim_ends_p,dibasep,cmetp);
	    (*antiacceptors_plus)[nmismatches] = List_push((*antiacceptors_plus)[nmismatches],(void *) hit);
	    found_known_p = true;
	  }
	}

	if (found_known_p == false && novelsplicingp == true) {
	  Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	  splice_pos = INDEX1PART;
	  nmismatches = 0;
	  while (mismatch_positions[nmismatches] <= splice_pos) { /* Must be <= */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }

	  while (splice_pos < stopi) {
	    if (mismatch_positions[nmismatches] <= splice_pos) { /* Must be <= */
	      debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	      nmismatches++;
	    }

	    if (gbuffer[splice_pos+2] == 'T' || gbuffer[splice_pos+2] == 'C') {
	      if (gbuffer[splice_pos+1] == 'G') {
		/* Looking for GT/GC */
		/* End 1.  sense direction, plus strand, donor.  Support is chimera_pos = splice_pos. */
		if (segment_left + splice_pos + 1 >= DONOR_MODEL_LEFT_MARGIN) {
		  model_left = segment_left + splice_pos + 1 - DONOR_MODEL_LEFT_MARGIN;
		  Genome_fill_buffer_simple(genome,model_left,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1,gbuffer1);
		  debug4e(printf("splice_pos %d, support %d, donor model: %s (prob = %f)\n",
				 splice_pos,splice_pos,gbuffer1,Maxent_donor_prob(gbuffer1)));
		  if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob = Maxent_donor_prob(gbuffer1))) {
		    debug4e(printf("=>sense, plus, donor: %f at %d (%d mismatches)\n",prob,splice_pos+1,nmismatches));
		    debug4e(printf("q: %s\ng: %s\n",queryuc_ptr,gbuffer));
		    hit = Substring_new_donor(splice_pos+1,nmismatches,ncolordiffs,prob,/*left*/segment_left,
					      querylength,/*plusp*/true,/*sensep*/true,
					      gbuffer,queryuc_ptr,segment->chrnum,segment->chroffset,/*knownp*/false,
					      trim_ends_p,dibasep,cmetp);
		    (*donors_plus)[nmismatches] = List_push((*donors_plus)[nmismatches],(void *) hit);
		  }
		}
#if 0
		/* Cannot advance 2 because C at splice_pos+2 could match [GC] at splice_pos+1 */
		splice_pos++;
#endif
	    
	      } else if (gbuffer[splice_pos+2] == 'T' && gbuffer[splice_pos+1] == 'C') {
		/* Looking for CT (revcomp of AG) */
		/* End 6.  antisense direction, plus strand, acceptor.  Support is chimera_pos = splice_pos. */
		if (segment_left + splice_pos >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
		  model_left = segment_left + splice_pos - ACCEPTOR_MODEL_RIGHT_MARGIN;
		  Genome_fill_buffer_simple(genome,model_left,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1,gbuffer1);
		  make_complement_buffered(gbuffer2,gbuffer1,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1);
		  debug4e(printf("splice_pos %d, support %d, acceptor model rc: %s (prob = %f)\n",
				 splice_pos,splice_pos,gbuffer2,Maxent_acceptor_prob(gbuffer2)));
		  if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob = Maxent_acceptor_prob(gbuffer2))) {
		    debug4e(printf("=>antisense, plus, antiacceptor: %f at %d (%d mismatches)\n",prob,splice_pos+1,nmismatches));
		    debug4e(printf("q: %s\ng: %s\n",queryuc_ptr,gbuffer));
		    hit = Substring_new_acceptor(splice_pos+1,nmismatches,ncolordiffs,prob,/*left*/segment_left,
						 querylength,/*plusp*/true,/*sensep*/false,
						 gbuffer,queryuc_ptr,segment->chrnum,segment->chroffset,/*knownp*/false,
						 trim_ends_p,dibasep,cmetp);
		    (*antiacceptors_plus)[nmismatches] = List_push((*antiacceptors_plus)[nmismatches],(void *) hit);
		  }
		}
		/* Can advance 2 because T at splice_pos+2 will not match [GC] at splice_pos+1 */
		splice_pos++;
	      }
	    }
	    splice_pos++;
	  }
	}
	/* End of novelsplicing */

      }
    }

    if (floors_to_pos3[segment->querypos3] <= max_mismatches_allowed) {
      trimpos = Genome_trim_left(queryuc_ptr,query_compress,genome_blocks,snp_blocks,
				 /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
				 dibasep,cmetp,/*plusp*/true);

      nmismatches_right = Genome_mismatches_right(mismatch_positions,colordiffs,max_mismatches_allowed,
						  queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						  /*left*/segment_left,/*pos5*/trimpos+1,/*pos3*/querylength,
						  dibasep,cmetp,/*plusp*/true);

      debug4e(
	       printf("For %d down to trim %d, %d mismatches on right at:",querylength,trimpos,nmismatches_right);
	       for (i = 0; i <= nmismatches_right; i++) {
		 printf(" %d",mismatch_positions[i]);
	       }
	       printf("\n");
	       );

      if (nmismatches_right == 0) {
	stopi = trimpos;
      } else {
	stopi = mismatch_positions[nmismatches_right-1];
      }
      if (stopi < 2) {
	stopi = 2;
      }
      debug4e(printf("Search from %d down to %d\n",query_lastpos-1,stopi));

      if (stopi < querylength - MIN_SPLICE_SUPPORT_DISTANT) {
	/* Advance through known splice sites */
	if (advancep == false && nsplicesites > 0 && *splicesites < segment_left + MIN_SPLICE_SUPPORT_DISTANT) {
	  j = 1;
	  while (j < nsplicesites && splicesites[j] < segment_left + MIN_SPLICE_SUPPORT_DISTANT) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= nsplicesites) {
	    j = binary_search(j >> 1,nsplicesites,splicesites,segment_left + MIN_SPLICE_SUPPORT_DISTANT);
	  } else {
	    j = binary_search(j >> 1,j,splicesites,segment_left + MIN_SPLICE_SUPPORT_DISTANT);
	  }
	  splicesites += j;
	  splicetypes += j;
	  nsplicesites -= j;
	  /* advancep = true; */
	}

	/* Known splicing */
	found_known_p = false;
	for (j = 0; j < nsplicesites && splicesites[j] < segment->diagonal; j++) {
	  if (splicetypes[j] == ACCEPTOR) {
	    if ((splice_pos = splicesites[j] - segment_left) > stopi) {
	      nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
							      dibasep,cmetp,/*plusp*/true);
	      if (gbuffer[0] == '\0') {
		Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	      }
	      debug4e(printf("Known acceptor for segmenti at %u (%d mismatches)\n",splicesites[j],nmismatches));
	      hit = Substring_new_acceptor(splice_pos,nmismatches,ncolordiffs,/*prob*/2.0,/*left*/segment_left,
					   querylength,/*plusp*/true,/*sensep*/true,
					   gbuffer,queryuc_ptr,segment->chrnum,segment->chroffset,/*knownp*/true,
					   trim_ends_p,dibasep,cmetp);
	      (*acceptors_plus)[nmismatches] = List_push((*acceptors_plus)[nmismatches],(void *) hit);
	      found_known_p = true;
	    }
	  } else if (splicetypes[j] == ANTIDONOR) {
	    if ((splice_pos = splicesites[j] - segment_left) > stopi) {
	      nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
							      dibasep,cmetp,/*plusp*/true);
	      if (gbuffer[0] == '\0') {
		Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	      }
	      debug4e(printf("Known antidonor for segmenti at %u (%d mismatches)\n",splicesites[j],nmismatches));
	      hit = Substring_new_donor(splice_pos,nmismatches,ncolordiffs,
					/*prob*/2.0,/*left*/segment_left,
					querylength,/*plusp*/true,/*sensep*/false,
					gbuffer,queryuc_ptr,segment->chrnum,segment->chroffset,/*knownp*/true,
					trim_ends_p,dibasep,cmetp);
	      (*antidonors_plus)[nmismatches] = List_push((*antidonors_plus)[nmismatches],(void *) hit);
	      found_known_p = true;
	    }
	  }
	}

	if (found_known_p == false && novelsplicingp == true) {
	  if (gbuffer[0] == '\0') {
	    Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	  }
	  splice_pos = query_lastpos - 1;
	  nmismatches = 0;
	  while (mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }

	  while (splice_pos > stopi) {
	    if (mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	      debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	      nmismatches++;
	    }

	    if (gbuffer[splice_pos-2] == 'A' || gbuffer[splice_pos-2] == 'G') {
	      if (gbuffer[splice_pos-2] == 'A' && gbuffer[splice_pos-1] == 'G') {
		/* Looking for AG */
		/* End 2. sense direction, plus strand, acceptor.  Support is querylength - chimera_pos = querylength - splice_pos. */
		if (segment_left + splice_pos >= ACCEPTOR_MODEL_LEFT_MARGIN) {
		  model_left = segment_left + splice_pos - ACCEPTOR_MODEL_LEFT_MARGIN;
		  Genome_fill_buffer_simple(genome,model_left,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1,gbuffer1);
		  debug4e(printf("splice_pos %d, support %d, acceptor model: %s (prob = %f)\n",
				 splice_pos,querylength-splice_pos,gbuffer1,Maxent_acceptor_prob(gbuffer1)));
		  if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob = Maxent_acceptor_prob(gbuffer1))) {
		    debug4e(printf("=>sense, plus, acceptor: %f at %d (%d mismatches)\n",prob,splice_pos,nmismatches));
		    debug4e(printf("q: %s\ng: %s\n",queryuc_ptr,gbuffer));
		    hit = Substring_new_acceptor(splice_pos,nmismatches,ncolordiffs,prob,/*left*/segment_left,
						 querylength,/*plusp*/true,/*sensep*/true,
						 gbuffer,queryuc_ptr,segment->chrnum,segment->chroffset,/*knownp*/false,
						 trim_ends_p,dibasep,cmetp);
		    (*acceptors_plus)[nmismatches] = List_push((*acceptors_plus)[nmismatches],(void *) hit);
		  }
		}
		/* Can advance 2, since A at splice_pos-2 will not match [GC] at splice_pos-1 */
		splice_pos--;

	      } else if (gbuffer[splice_pos-1] == 'C') {
		/* Looking for AC/GC (revcomp of GT/GC) */
		/* End 5. antisense direction, plus strand, donor.  Support is querylength - chimera_pos = querylength - splice_pos. */
		if (segment_left + splice_pos - 1 >= DONOR_MODEL_RIGHT_MARGIN) {
		  model_left = segment_left + splice_pos - DONOR_MODEL_RIGHT_MARGIN - 1;
		  Genome_fill_buffer_simple(genome,model_left,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1,gbuffer1);
		  make_complement_buffered(gbuffer2,gbuffer1,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1);
		  debug4e(printf("splice_pos %d, support %d, donor model rc: %s (prob = %f)\n",
				 splice_pos,querylength-splice_pos,gbuffer2,Maxent_donor_prob(gbuffer2)));
		  if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob = Maxent_donor_prob(gbuffer2))) {
		    debug4e(printf("=>antisense, plus, antidonor: %f at %d (%d mismatches)\n",prob,splice_pos,nmismatches));
		    debug4e(printf("q: %s\ng: %s\n",queryuc_ptr,gbuffer));
		    hit = Substring_new_donor(splice_pos,nmismatches,ncolordiffs,prob,/*left*/segment_left,
					      querylength,/*plusp*/true,/*sensep*/false,
					      gbuffer,queryuc_ptr,segment->chrnum,segment->chroffset,/*knownp*/false,
					      trim_ends_p,dibasep,cmetp);
		    (*antidonors_plus)[nmismatches] = List_push((*antidonors_plus)[nmismatches],(void *) hit);
		  }
		}
#if 0
		/* Cannot advance 2, since G at splice_pos-2 could match [GC] at splice_pos-1 */
		splice_pos--;
#endif
	      }
	    }
	    splice_pos--;
	  }
	}
	/* End of novelsplicing */

      }
    }
  }
    
  return;
}


static void
identify_spliceends_minus (List_T **donors_minus, List_T **antidonors_minus,
			   List_T **acceptors_minus, List_T **antiacceptors_minus,
			   struct Splicesegment_T *minus_segments, int minus_nsegments,
			   UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome,
			   char *queryuc_ptr, /* for debugging */ char *queryrc,
			   Floors_T floors, int querylength, int query_lastpos,
			   char *gbuffer, Compress_T query_compress,
			   Genomicpos_T *splicesites, Splicetype_T *splicetypes, int nsplicesites,
			   bool novelsplicingp, bool canonicalp, 
			   int max_mismatches_allowed, bool trim_ends_p, bool dibasep, bool cmetp) {
  Splicesegment_T segment;
  Substring_T hit;
  Genomicpos_T segment_left, model_left;
  int nmismatches, stopi, j;
  int splice_pos;
  char gbuffer1[MAXENT_MAXLENGTH+1], gbuffer2[MAXENT_MAXLENGTH+1];
  double prob;

  int trimpos;
  int mismatch_positions[MAX_QUERYLENGTH+1], colordiffs[MAX_QUERYLENGTH];
  int nmismatches_left, nmismatches_right;
  int ncolordiffs = 0;
  int *floors_from_neg3, *floors_to_pos3;
  bool advancep, found_known_p;

#ifdef DEBUG4E
  int i;
#endif


  floors_from_neg3 = floors->scorefrom[-3];
  floors_to_pos3 = floors->scoreto[query_lastpos+3];

  for (segment = minus_segments; segment < &(minus_segments[minus_nsegments]); segment++) {
    segment_left = segment->diagonal - querylength;
    debug4e(printf("identify_spliceends_minus: Getting genome at diagonal %u (querypos %d..%d) + 12 - querylength %d = %u\n",
		   segment->diagonal,segment->querypos5,segment->querypos3,querylength,segment_left));
    debug4e(
	     Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	     printf("genome   0..: %s\n",gbuffer);
	     printf("query.rc 0..: %s\n",queryrc);
	     );

    gbuffer[0] = '\0';		/* Signals that we haven't filled it yet */
    advancep = false;

    if (floors_to_pos3[segment->querypos3] <= max_mismatches_allowed) {
      trimpos = Genome_trim_right(queryuc_ptr,query_compress,genome_blocks,snp_blocks,
				  /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
				  dibasep,cmetp,/*plusp*/false);

      nmismatches_left = Genome_mismatches_left(mismatch_positions,colordiffs,max_mismatches_allowed,
						queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						/*left*/segment_left,/*pos5*/0,/*pos3*/trimpos,
						dibasep,cmetp,/*plusp*/false);

      debug4e(
	       printf("For %d up to trim %d, %d mismatches on left at:",0,trimpos,nmismatches_left);
	       for (i = 0; i <= nmismatches_left; i++) {
		 printf(" %d",mismatch_positions[i]);
	       }
	       printf("\n");
	       );

      if (nmismatches_left == 0) {
	stopi = trimpos;
      } else {
	stopi = mismatch_positions[nmismatches_left-1];
      }
      if (stopi > querylength - 2) {
	stopi = querylength - 2;
      }
      debug4e(printf("Search from %d up to %d\n",INDEX1PART,stopi));

      if (stopi >= MIN_SPLICE_SUPPORT_DISTANT) {
	/* Advance through known splice sites */
	if (/* advancep == false && */ nsplicesites > 0 && *splicesites < segment_left + MIN_SPLICE_SUPPORT_DISTANT) {
	  j = 1;
	  while (j < nsplicesites && splicesites[j] < segment_left + MIN_SPLICE_SUPPORT_DISTANT) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= nsplicesites) {
	    j = binary_search(j >> 1,nsplicesites,splicesites,segment_left + MIN_SPLICE_SUPPORT_DISTANT);
	  } else {
	    j = binary_search(j >> 1,j,splicesites,segment_left + MIN_SPLICE_SUPPORT_DISTANT);
	  }
	  splicesites += j;
	  splicetypes += j;
	  nsplicesites -= j;
	  advancep = true;
	}

	/* Known splicing */
	found_known_p = false;
	for (j = 0; j < nsplicesites && splicesites[j] <= segment_left + stopi; j++) {
	  if (splicetypes[j] == ANTIACCEPTOR) {
	    splice_pos = splicesites[j] - segment_left;
	    nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							    /*left*/segment_left,/*pos5*/0,/*pos3*/splice_pos,
							    dibasep,cmetp,/*plusp*/false);
	    if (gbuffer[0] == '\0') {
	      Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	    }
	    debug4e(printf("Known acceptor for segmenti at %u (%d mismatches)\n",splicesites[j],nmismatches));
	    hit = Substring_new_acceptor(/*acceptor_pos*/querylength-splice_pos,nmismatches,ncolordiffs,
					 /*prob*/2.0,/*left*/segment_left,querylength,
					 /*plusp*/false,/*sensep*/true,gbuffer,queryuc_ptr,segment->chrnum,
					 segment->chroffset,/*knownp*/true,trim_ends_p,dibasep,cmetp);
	    (*acceptors_minus)[nmismatches] = List_push((*acceptors_minus)[nmismatches],(void *) hit);
	    found_known_p = true;

	  } else if (splicetypes[j] == DONOR) {
	    splice_pos = splicesites[j] - segment_left;
	    nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							    /*left*/segment_left,/*pos5*/0,/*pos3*/splice_pos,
							    dibasep,cmetp,/*plusp*/false);
	    if (gbuffer[0] == '\0') {
	      Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	    }
	    debug4e(printf("Known antidonor for segmenti at %u (%d mismatches)\n",splicesites[j],nmismatches));
	    hit = Substring_new_donor(/*donor_pos*/querylength-splice_pos,nmismatches,ncolordiffs,
				      /*prob*/2.0,/*left*/segment_left,querylength,
				      /*plusp*/false,/*sensep*/false,gbuffer,queryuc_ptr,segment->chrnum,
				      segment->chroffset,/*knownp*/true,trim_ends_p,dibasep,cmetp);
	    (*antidonors_minus)[nmismatches] = List_push((*antidonors_minus)[nmismatches],(void *) hit);
	    found_known_p = true;
	  }
	}

	if (found_known_p == false && novelsplicingp == true) {
	  Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	  splice_pos = INDEX1PART;
	  nmismatches = 0;
	  while (mismatch_positions[nmismatches] <= splice_pos) { /* Must be <= */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }

	  while (splice_pos < stopi) {
	    if (mismatch_positions[nmismatches] <= splice_pos) { /* Must be <= */
	      debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	      nmismatches++;
	    }

	    if (gbuffer[splice_pos+2] == 'T' || gbuffer[splice_pos+2] == 'C') {
	      if (gbuffer[splice_pos+2] == 'T' && gbuffer[splice_pos+1] == 'C') {
		/* Looking for CT (revcomp of AG) */
		/* End 4. sense direction, minus strand, acceptor.  Support is querylength - chimera_pos = querylength - (querylength - splice_pos) = splice_pos. */
		if (segment_left + splice_pos >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
		  model_left = segment_left + splice_pos - ACCEPTOR_MODEL_RIGHT_MARGIN;
		  Genome_fill_buffer_simple(genome,model_left,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1,gbuffer1);
		  make_complement_buffered(gbuffer2,gbuffer1,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1);
		  debug4e(printf("splice_pos %d, support %d, acceptor model rc: %s (prob = %f)\n",
				 splice_pos,splice_pos,gbuffer2,Maxent_acceptor_prob(gbuffer2)));
		  if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob = Maxent_acceptor_prob(gbuffer2))) {
		    debug4e(printf("=>sense, minus, acceptor: %f at %d (%d mismatches)\n",prob,querylength-1-splice_pos,nmismatches));
		    debug4e(printf("q.rc: %s\ng:    %s\n",queryrc,gbuffer));
		    hit = Substring_new_acceptor(/*acceptor_pos*/querylength-1-splice_pos,nmismatches,ncolordiffs,
						 prob,/*left*/segment_left,querylength,
						 /*plusp*/false,/*sensep*/true,gbuffer,queryuc_ptr,segment->chrnum,
						 segment->chroffset,/*knownp*/false,trim_ends_p,dibasep,cmetp);
		    (*acceptors_minus)[nmismatches] = List_push((*acceptors_minus)[nmismatches],(void *) hit);
		  }
		}
		/* Can advance 2 because T at splice_pos+2 will not match [CG] at splice_pos+1 */
		splice_pos++;
		
	      } else if (gbuffer[splice_pos+1] == 'G') {
		/* Looking for GT/GC */
		/* End 7. antisense direction, minus strand, donor.  Support is querylength - chimera_pos = querylength - (querylength - splice_pos) = splice_pos. */
		if (segment_left + splice_pos + 1 >= DONOR_MODEL_LEFT_MARGIN) {
		  model_left = segment_left + splice_pos + 1 - DONOR_MODEL_LEFT_MARGIN;
		  Genome_fill_buffer_simple(genome,model_left,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1,gbuffer1);
		  debug4e(printf("splice_pos %d, support %d, donor model: %s (prob = %f)\n",
				 splice_pos,splice_pos,gbuffer1,Maxent_donor_prob(gbuffer1)));
		  if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob = Maxent_donor_prob(gbuffer1))) {
		    debug4e(printf("=>antisense, minus, antidonor: %f at %d (%d mismatches)\n",prob,querylength-1-splice_pos,nmismatches));
		    debug4e(printf("q.rc: %s\ng.rc: %s\n",queryrc,gbuffer));
		    hit = Substring_new_donor(/*donor_pos*/querylength-1-splice_pos,nmismatches,ncolordiffs,
					      prob,/*left*/segment_left,querylength,
					      /*plusp*/false,/*sensep*/false,gbuffer,queryuc_ptr,segment->chrnum,
					      segment->chroffset,/*knownp*/false,trim_ends_p,dibasep,cmetp);
		    (*antidonors_minus)[nmismatches] = List_push((*antidonors_minus)[nmismatches],(void *) hit);
		  }
		}
#if 0
		/* Cannot advance 2 because C at splice_pos+2 could match [CG] at splice_pos+1 */
		splice_pos++;
#endif
	      }
	    }
	    splice_pos++;
	  }
	}
	/* End of novelsplicing */

      }
    }
      
    if (floors_from_neg3[segment->querypos5] <= max_mismatches_allowed) {
      trimpos = Genome_trim_left(queryuc_ptr,query_compress,genome_blocks,snp_blocks,
				 /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
				 dibasep,cmetp,/*plusp*/false);

      nmismatches_right = Genome_mismatches_right(mismatch_positions,colordiffs,max_mismatches_allowed,
						  queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						  /*left*/segment_left,/*pos5*/trimpos+1,/*pos3*/querylength,
						  dibasep,cmetp,/*plusp*/false);

      debug4e(
	       printf("For %d down to trim %d, %d mismatches on right at:",querylength,trimpos,nmismatches_right);
	       for (i = 0; i <= nmismatches_right; i++) {
		 printf(" %d",mismatch_positions[i]);
	       }
	       printf("\n");
	       );

      if (nmismatches_right == 0) {
	stopi = trimpos;
      } else {
	stopi = mismatch_positions[nmismatches_right-1];
      }
      if (stopi < 2) {
	stopi = 2;
      }
      debug4e(printf("Search from %d down to %d\n",query_lastpos-1,stopi));

      if (stopi < querylength - MIN_SPLICE_SUPPORT_DISTANT) {
	/* Advance through known splice sites */
	if (advancep == false && nsplicesites > 0 && *splicesites < segment_left + MIN_SPLICE_SUPPORT_DISTANT) {
	  j = 1;
	  while (j < nsplicesites && splicesites[j] < segment_left + MIN_SPLICE_SUPPORT_DISTANT) {
	    j <<= 1;		/* gallop by 2 */
	  }
	  if (j >= nsplicesites) {
	    j = binary_search(j >> 1,nsplicesites,splicesites,segment_left + MIN_SPLICE_SUPPORT_DISTANT);
	  } else {
	    j = binary_search(j >> 1,j,splicesites,segment_left + MIN_SPLICE_SUPPORT_DISTANT);
	  }
	  splicesites += j;
	  splicetypes += j;
	  nsplicesites -= j;
	  /* advancep = true; */
	}

	/* Known splicing */
	found_known_p = false;
	for (j = 0; j < nsplicesites && splicesites[j] < segment->diagonal; j++) {
	  if (splicetypes[j] == ANTIDONOR) {
	    if ((splice_pos = splicesites[j] - segment_left) > stopi) {
	      nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
							      dibasep,cmetp,/*plusp*/false);
	      if (gbuffer[0] == '\0') {
		Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	      }
	      debug4e(printf("Known donor for segmenti at %u (%d mismatches)\n",splicesites[j],nmismatches));
	      hit = Substring_new_donor(/*donor_pos*/querylength-splice_pos,nmismatches,ncolordiffs,
					/*prob*/2.0,/*left*/segment_left,querylength,
					/*plusp*/false,/*sensep*/true,gbuffer,queryuc_ptr,segment->chrnum,
					segment->chroffset,/*knownp*/true,trim_ends_p,dibasep,cmetp);
	      (*donors_minus)[nmismatches] = List_push((*donors_minus)[nmismatches],(void *) hit);
	      found_known_p = true;
	    }

	  } else if (splicetypes[j] == ACCEPTOR) {
	    if ((splice_pos = splicesites[j] - segment_left) > stopi) {
	      nmismatches = Genome_count_mismatches_substring(&ncolordiffs,queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							      /*left*/segment_left,/*pos5*/splice_pos,/*pos3*/querylength,
							      dibasep,cmetp,/*plusp*/false);
	      if (gbuffer[0] == '\0') {
		Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	      }
	      debug4e(printf("Known antiacceptor for segmenti at %u (%d mismatches)\n",splicesites[j],nmismatches));
	      hit = Substring_new_acceptor(/*acceptor_pos*/querylength-splice_pos,nmismatches,ncolordiffs,
					   /*prob*/2.0,/*left*/segment_left,querylength,
					   /*plusp*/false,/*sensep*/false,gbuffer,queryuc_ptr,segment->chrnum,
					   segment->chroffset,/*knownp*/true,trim_ends_p,dibasep,cmetp);
	      (*antiacceptors_minus)[nmismatches] = List_push((*antiacceptors_minus)[nmismatches],(void *) hit);
	      found_known_p = true;
	    }
	  }
	}

	if (found_known_p == false && novelsplicingp == true) {
	  if (gbuffer[0] == '\0') {
	    Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	  }
	  splice_pos = query_lastpos - 1;
	  nmismatches = 0;
	  while (mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	    debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	    nmismatches++;
	  }

	  while (splice_pos > stopi) {
	    if (mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	      debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	      nmismatches++;
	    }

	    if (gbuffer[splice_pos-2] == 'A' || gbuffer[splice_pos-2] == 'G') {
	      if (gbuffer[splice_pos-1] == 'C') {
		/* Looking for AC/GC (revcomp of GT/GC) */
		/* End 3. sense direction, minus strand, donor.  Support is chimera_pos = querylength - splice_pos. */
		if (segment_left + splice_pos - 1 >= DONOR_MODEL_RIGHT_MARGIN) {
		  model_left = segment_left + splice_pos - DONOR_MODEL_RIGHT_MARGIN - 1;
		  Genome_fill_buffer_simple(genome,model_left,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1,gbuffer1);
		  make_complement_buffered(gbuffer2,gbuffer1,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1);
		  debug4e(printf("splice_pos %d, support %d, donor model rc: %s (prob = %f)\n",
				 splice_pos,querylength-splice_pos,gbuffer2,Maxent_donor_prob(gbuffer2)));
		  if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob = Maxent_donor_prob(gbuffer2))) {
		    debug4e(printf("=>sense, minus, donor: %f at %d (%d mismatches)\n",prob,querylength-splice_pos,nmismatches));
		    debug4e(printf("q.rc: %s\ng.rc: %s\n",queryrc,gbuffer));
		    hit = Substring_new_donor(/*donor_pos*/querylength-splice_pos,nmismatches,ncolordiffs,
					      prob,/*left*/segment_left,querylength,
					      /*plusp*/false,/*sensep*/true,gbuffer,queryuc_ptr,segment->chrnum,
					      segment->chroffset,/*knownp*/false,trim_ends_p,dibasep,cmetp);
		    (*donors_minus)[nmismatches] = List_push((*donors_minus)[nmismatches],(void *) hit);
		  }
		}
#if 0
		/* Cannot advance 2 because G at splice_pos-2 could match [CG] at splice_pos-1 */
		splice_pos--;
#endif

	      } else if (gbuffer[splice_pos-2] == 'A' && gbuffer[splice_pos-1] == 'G') {
		/* Looking for AG */
		/* End 8. antisense direction, minus strand, acceptor.  Support is chimera_pos = querylength - splice_pos. */
		if (segment_left + splice_pos >= ACCEPTOR_MODEL_LEFT_MARGIN) {
		  model_left = segment_left + splice_pos - ACCEPTOR_MODEL_LEFT_MARGIN;
		  Genome_fill_buffer_simple(genome,model_left,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1,gbuffer1);
		  debug4e(printf("splice_pos %d, support %d, acceptor model: %s (prob = %f)\n",
				 splice_pos,querylength-splice_pos,gbuffer1,Maxent_acceptor_prob(gbuffer1)));
		  if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob = Maxent_acceptor_prob(gbuffer1))) {
		    debug4e(printf("=>antisense, minus, antiacceptor: %f at %d (%d mismatches)\n",prob,querylength-splice_pos,nmismatches));
		    debug4e(printf("q.rc: %s\ng.rc: %s\n",queryrc,gbuffer));
		    hit = Substring_new_acceptor(/*acceptor_pos*/querylength-splice_pos,nmismatches,ncolordiffs,
						 prob,/*left*/segment_left,querylength,
						 /*plusp*/false,/*sensep*/false,gbuffer,queryuc_ptr,segment->chrnum,
						 segment->chroffset,/*knownp*/false,trim_ends_p,dibasep,cmetp);
		    (*antiacceptors_minus)[nmismatches] = List_push((*antiacceptors_minus)[nmismatches],(void *) hit);
		  }
		}
		/* Can advance 2 because A at splice_pos-2 will not match [CG] at splice_pos-1 */
		splice_pos--;
	      }
	    }
	    splice_pos--;
	  }
	}
	/* End of novelsplicing */

      }
    }
  }

  return;
}



static void
fetch_positions_for_all_12mers (T this, Indexdb_T indexdb, Indexdb_T indexdb2, int query_lastpos) {
  int querypos;
  
  /* querypos -2, -1, query_lastpos+1, and query_lastpos+2 are special cases */
  /* if allvalidp is true, then 0 and query_lastpos should have been done already */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->plus_retrievedp[querypos] == false) {
      /* FORMULA */
      this->plus_positions[querypos] = 
	Indexdb_read_inplace(&(this->plus_npositions[querypos]),indexdb,this->forward_oligos[querypos]);
      debug(printf("Retrieving at querypos %d, plus_npositions = %d\n",
		   querypos,this->plus_npositions[querypos]));
      this->plus_retrievedp[querypos] = true;
      this->plus_allocp[querypos] = false;
    }
    if (this->minus_retrievedp[querypos] == false) {
      /* FORMULA */
      this->minus_positions[querypos] = 
	Indexdb_read_inplace(&(this->minus_npositions[querypos]),indexdb2,this->revcomp_oligos[querypos]);
      debug(printf("Retrieving at querypos %d, minus_npositions = %d\n",
		   querypos,this->minus_npositions[querypos]));
      this->minus_retrievedp[querypos] = true;
      this->minus_allocp[querypos] = false;
    }
  }

  this->all_positions_fetched_p = true;
  
  return;
}


static List_T
find_splicepairs_distant (int *found_score, int *nsplicepairs, List_T distantsplicing,
			  List_T *donors_plus, List_T *antidonors_plus,
			  List_T *acceptors_plus, List_T *antiacceptors_plus,
			  List_T *donors_minus, List_T *antidonors_minus,
			  List_T *acceptors_minus, List_T *antiacceptors_minus,
			  Genomicpos_T shortsplicedist, int localsplicing_penalty, int distantsplicing_penalty,
			  int min_distantsplicing_end_matches, double min_distantsplicing_identity,
			  int querylength, int nmismatches_allowed, int maxchimerapaths) {
  List_T p, q, qsave;
  Substring_T donor, acceptor;
  int min_endlength_1, min_endlength_2, nmismatches1, nmismatches2, pos;
  int splicing_penalty;
  Genomicpos_T distance;
  bool shortdistancep;
  double nonidentity = 1.0 - min_distantsplicing_identity;
  
  debug(printf("Starting find_splicepairs_distant with nonidentity %f\n",nonidentity));
  debug4l(printf("Starting find_splicepairs_distant with nonidentity %f\n",nonidentity));
  
  if (nonidentity == 0.0) {
    nmismatches_allowed = 0;
  }

  for (nmismatches1 = 0; nmismatches1 <= nmismatches_allowed; nmismatches1++) {
    nmismatches2 = nmismatches_allowed - nmismatches1;

    if (nonidentity == 0.0) {
      min_endlength_1 = min_endlength_2 = min_distantsplicing_end_matches;
    } else {
      min_endlength_1 = rint((double) nmismatches1/nonidentity);
      if (min_endlength_1 < min_distantsplicing_end_matches) {
	min_endlength_1 = min_distantsplicing_end_matches;
      }
      min_endlength_2 = rint((double) nmismatches2/nonidentity);
      if (min_endlength_2 < min_distantsplicing_end_matches) {
	min_endlength_2 = min_distantsplicing_end_matches;
      }
    }

    debug4l(printf("  nmismatches1 = %d, nmismatches2 = %d, min_endlength_1 = %d, min_endlength_2 = %d\n",
		   nmismatches1,nmismatches2,min_endlength_1,min_endlength_2));

    /************************************************************************
     *   Same strands
     ************************************************************************/

    /* 1.  End 1 to End 2.  Same strands. */
    p = donors_plus[nmismatches1];
    q = acceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_known_distant (%d+%d mismatches): donors+ (%d) to acceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    debug4l(printf("nsplicepairs %d, maxchimerapaths %d\n",*nsplicepairs,maxchimerapaths));
    while (p != NULL && q != NULL && *nsplicepairs < maxchimerapaths) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("donor at %u and acceptor at %u\n",Substring_genomicstart(donor),Substring_genomicstart(acceptor)));

      if ((pos = Substring_chimera_pos(donor)) < min_endlength_1) {
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	p = p->rest;
      } else if (pos < Substring_chimera_pos(acceptor)) {
	p = p->rest;
      } else if (pos > Substring_chimera_pos(acceptor)) {
	q = q->rest;
      } else {
	/* Generate all pairs at this chimera_pos */
	qsave = q;
	while (p != NULL && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	      shortdistancep = false;
	      splicing_penalty = distantsplicing_penalty;
	    } else if (Substring_genomicstart(acceptor) > Substring_genomicstart(donor)) {
	      distance = Substring_genomicstart(acceptor) - Substring_genomicstart(donor);
	      if (distance <= shortsplicedist) {
		shortdistancep = true;
		splicing_penalty = localsplicing_penalty;
	      } else {
		shortdistancep = false;
		splicing_penalty = distantsplicing_penalty;
	      }
	    } else {
	      distance = Substring_genomicstart(donor) - Substring_genomicstart(acceptor);
	      shortdistancep = false; /* scramble */
	      splicing_penalty = distantsplicing_penalty;
	    }
	    debug4ld(printf("1-2. Pushing a candidate at chimera_pos %d (%d..%d), donor %u to acceptor %u.  shortdistancep = %d\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor),shortdistancep));
#if 0
	    if (shortdistancep) {
	      *localsplicing = List_push(*localsplicing,
					 (void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								    shortdistancep,splicing_penalty,/*support*/querylength,
								    querylength,/*copyp*/true,/*sensep*/true));
	    } else {
#endif
	      distantsplicing = List_push(distantsplicing,
					  (void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								     shortdistancep,splicing_penalty,/*support*/querylength,
								     querylength,/*copyp*/true,/*sensep*/true));
#if 0
	    }
#endif
	    (*nsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }
  
    /* 4. End 3 to End 4.  Same strands. */
    p = donors_minus[nmismatches1];
    q = acceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_known_distant (%d+%d mismatches): donors- (%d) to acceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    debug4l(printf("nsplicepairs %d, maxchimerapaths %d\n",*nsplicepairs,maxchimerapaths));
    while (p != NULL && q != NULL && *nsplicepairs < maxchimerapaths) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("donor at %u and acceptor at %u\n",Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	
      if ((pos = Substring_chimera_pos(donor)) < min_endlength_1) {
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	p = p->rest;
      } else if (pos < Substring_chimera_pos(acceptor)) {
	p = p->rest;
      } else if (pos > Substring_chimera_pos(acceptor)) {
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	      shortdistancep = false;
	      splicing_penalty = distantsplicing_penalty;
	    } else if (Substring_genomicstart(acceptor) > Substring_genomicstart(donor)) {
	      distance = Substring_genomicstart(acceptor) - Substring_genomicstart(donor);
	      shortdistancep = false; /* scramble */
	      splicing_penalty = distantsplicing_penalty;
	    } else {
	      distance = Substring_genomicstart(donor) - Substring_genomicstart(acceptor);
	      if (distance <= shortsplicedist) {
		shortdistancep = true;
		splicing_penalty = localsplicing_penalty;
	      } else {
		shortdistancep = false;
		splicing_penalty = distantsplicing_penalty;
	      }
	    }
	    debug4ld(printf("3-4. Pushing a candidate at chimera_pos %d (%d..%d), donor %u to acceptor %u.  shortdistancep = %d.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor),shortdistancep));
#if 0
	    if (shortdistancep) {
	      *localsplicing = List_push(*localsplicing,
					 (void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								    shortdistancep,distantsplicing_penalty,/*support*/querylength,
								    querylength,/*copyp*/true,/*sensep*/true));
	    } else {
#endif
	      distantsplicing = List_push(distantsplicing,
					  (void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								     shortdistancep,splicing_penalty,/*support*/querylength,
								     querylength,/*copyp*/true,/*sensep*/true));
#if 0
	    }
#endif
	    (*nsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }
  
    /* 5. End 5 to End 6.  Same strands. */
    p = antidonors_plus[nmismatches1];
    q = antiacceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_known_distant (%d+%d mismatches): antidonors+ (%d) to antiacceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    debug4l(printf("nsplicepairs %d, maxchimerapaths %d\n",*nsplicepairs,maxchimerapaths));
    while (p != NULL && q != NULL && *nsplicepairs < maxchimerapaths) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("donor at %u and acceptor at %u\n",Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	
      if ((pos = Substring_chimera_pos(donor)) < min_endlength_2) {
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	p = p->rest;
      } else if (pos < Substring_chimera_pos(acceptor)) {
	p = p->rest;
      } else if (pos > Substring_chimera_pos(acceptor)) {
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	      shortdistancep = false;
	      splicing_penalty = distantsplicing_penalty;
	    } else if (Substring_genomicstart(acceptor) > Substring_genomicstart(donor)) {
	      distance = Substring_genomicstart(acceptor) - Substring_genomicstart(donor);
	      shortdistancep = false; /* scramble */
	      splicing_penalty = distantsplicing_penalty;
	    } else {
	      distance = Substring_genomicstart(donor) - Substring_genomicstart(acceptor);
	      if (distance <= shortsplicedist) {
		shortdistancep = true;
		splicing_penalty = localsplicing_penalty;
	      } else {
		shortdistancep = false;
		splicing_penalty = distantsplicing_penalty;
	      }
	    }
	    debug4ld(printf("5-6. Pushing a candidate at chimera_pos %d (%d..%d), donor %u to acceptor %u.  shortdistancep = %d\n",
			    pos,min_endlength_2,querylength-min_endlength_1,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor),shortdistancep));
#if 0
	    if (shortdistancep) {
	      *localsplicing = List_push(*localsplicing,
					 (void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								    shortdistancep,distantsplicing_penalty,/*support*/querylength,
								    querylength,/*copyp*/true,/*sensep*/false));
	    } else {
#endif
	      distantsplicing = List_push(distantsplicing,
					  (void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								     shortdistancep,splicing_penalty,/*support*/querylength,
								     querylength,/*copyp*/true,/*sensep*/false));
#if 0
	    }
#endif
	    (*nsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }
  
    /* 8. End 7 to End 8.  Same strands. */
    p = antidonors_minus[nmismatches1];
    q = antiacceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_known_distant (%d+%d mismatches): antidonors- (%d) to antiacceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    debug4l(printf("nsplicepairs %d, maxchimerapaths %d\n",*nsplicepairs,maxchimerapaths));
    while (p != NULL && q != NULL && *nsplicepairs < maxchimerapaths) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("donor at %u and acceptor at %u\n",Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	
      if ((pos = Substring_chimera_pos(donor)) < min_endlength_2) {
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	p = p->rest;
      } else if (pos < Substring_chimera_pos(acceptor)) {
	p = p->rest;
      } else if (pos > Substring_chimera_pos(acceptor)) {
	q = q->rest;
      } else {
	qsave = q;
	  
	while (p != NULL && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	      shortdistancep = false;
	      splicing_penalty = distantsplicing_penalty;
	    } else if (Substring_genomicstart(acceptor) > Substring_genomicstart(donor)) {
	      distance = Substring_genomicstart(acceptor) - Substring_genomicstart(donor);
	      if (distance <= shortsplicedist) {
		shortdistancep = true;
		splicing_penalty = localsplicing_penalty;
	      } else {
		shortdistancep = false;
		splicing_penalty = distantsplicing_penalty;
	      }
	    } else {
	      distance = Substring_genomicstart(donor) - Substring_genomicstart(acceptor);
	      shortdistancep = false; /* scramble */
	      splicing_penalty = distantsplicing_penalty;
	    }
	    debug4ld(printf("7-8. Pushing a candidate at chimera_pos %d (%d..%d), donor %u to acceptor %u.  shortdistancep = %d.\n",
			    pos,min_endlength_2,querylength-min_endlength_1,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor),shortdistancep));
#if 0
	    if (shortdistancep) {
	      *localsplicing = List_push(*localsplicing,
					 (void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								    shortdistancep,distantsplicing_penalty,/*support*/querylength,
								    querylength,/*copyp*/true,/*sensep*/false));
	    } else {
#endif
	      distantsplicing = List_push(distantsplicing,
					  (void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								     shortdistancep,distantsplicing_penalty,/*support*/querylength,
								     querylength,/*copyp*/true,/*sensep*/false));
#if 0
	    }
#endif
	    (*nsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }


    /************************************************************************
     *   Different strands
     ************************************************************************/

    /* 2. End 1 to End 4.  Different strands. */
    p = donors_plus[nmismatches1];
    q = acceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_known_distant (%d+%d mismatches): donors+ (%d) to acceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    debug4l(printf("nsplicepairs %d, maxchimerapaths %d\n",*nsplicepairs,maxchimerapaths));
    while (p != NULL && q != NULL && *nsplicepairs < maxchimerapaths) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("donor at %u and acceptor at %u\n",Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	
      if ((pos = Substring_chimera_pos(donor)) < min_endlength_1) {
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	p = p->rest;
      } else if (pos < Substring_chimera_pos(acceptor)) {
	p = p->rest;
      } else if (pos > Substring_chimera_pos(acceptor)) {
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) - pos) > (Substring_genomicstart(donor) + pos)) {
	      distance = (Substring_genomicstart(acceptor) - pos) - (Substring_genomicstart(donor) + pos);
	    } else {
	      distance = (Substring_genomicstart(donor) + pos) - (Substring_genomicstart(acceptor) - pos);
	    }
	    debug4ld(printf("1-4. Pushing a candidate at chimera_pos %d (%d..%d), donor %u to acceptor %u.  Different strands, so not shortdistance.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								   /*shortdistancep*/false,distantsplicing_penalty,/*support*/querylength,
								   querylength,/*copyp*/true,/*sensep*/true));
	    (*nsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }
  
    /* 3. End 3 to End 2.  Different strands. */
    p = donors_minus[nmismatches1];
    q = acceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_known_distant (%d+%d mismatches): donors- (%d) to acceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    debug4l(printf("nsplicepairs %d, maxchimerapaths %d\n",*nsplicepairs,maxchimerapaths));
    while (p != NULL && q != NULL && *nsplicepairs < maxchimerapaths) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("donor at %u and acceptor at %u\n",Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	
      if ((pos = Substring_chimera_pos(donor)) < min_endlength_1) {
	p = p->rest;
      } else if (pos > querylength - min_endlength_2) {
	p = p->rest;
      } else if (pos < Substring_chimera_pos(acceptor)) {
	p = p->rest;
      } else if (pos > Substring_chimera_pos(acceptor)) {
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if (Substring_genomicstart(acceptor) > Substring_genomicstart(donor)) {
	      distance = (Substring_genomicstart(acceptor) + pos) - (Substring_genomicstart(donor) - pos);
	    } else {
	      distance = (Substring_genomicstart(donor) - pos) - (Substring_genomicstart(acceptor) + pos);
	    }
	    debug4ld(printf("3-2. Pushing a candidate at chimera_pos %d (%d..%d), donor %u to acceptor %u.  Different strands so not shortdistance.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								   /*shortdistancep*/false,distantsplicing_penalty,/*support*/querylength,
								   querylength,/*copyp*/true,/*sensep*/true));
	    (*nsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }
  
  
    /* 6. End 5 to End 8.  Different strands. */
    p = antidonors_plus[nmismatches1];
    q = antiacceptors_minus[nmismatches2];
    debug4l(printf("find_splicepairs_known_distant (%d+%d mismatches): antidonors+ (%d) to antiacceptors- (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    debug4l(printf("nsplicepairs %d, maxchimerapaths %d\n",*nsplicepairs,maxchimerapaths));
    while (p != NULL && q != NULL && *nsplicepairs < maxchimerapaths) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("donor at %u and acceptor at %u\n",Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	
      if ((pos = Substring_chimera_pos(donor)) < min_endlength_2) {
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	p = p->rest;
      } else if (pos < Substring_chimera_pos(acceptor)) {
	p = p->rest;
      } else if (pos > Substring_chimera_pos(acceptor)) {
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) - pos) > (Substring_genomicstart(donor) + pos)) {
	      distance = (Substring_genomicstart(acceptor) - pos) - (Substring_genomicstart(donor) + pos);
	    } else {
	      distance = (Substring_genomicstart(donor) + pos) - (Substring_genomicstart(acceptor) - pos);
	    }
	    debug4ld(printf("5-8. Pushing a candidate at chimera_pos %d (%d..%d), donor %u to acceptor %u.  Different strands so not shortdistance.\n",
			    pos,min_endlength_2,querylength-min_endlength_1,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								   /*shortdistancep*/false,distantsplicing_penalty,/*support*/querylength,
								   querylength,/*copyp*/true,/*sensep*/false));
	    (*nsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }
  
    /* 7. End 7 to End 6.  Different strands. */
    p = antidonors_minus[nmismatches1];
    q = antiacceptors_plus[nmismatches2];
    debug4l(printf("find_splicepairs_known_distant (%d+%d mismatches): antidonors- (%d) to antiacceptors+ (%d)\n",
		   nmismatches1,nmismatches2,List_length(p),List_length(q)));
    debug4l(printf("nsplicepairs %d, maxchimerapaths %d\n",*nsplicepairs,maxchimerapaths));
    while (p != NULL && q != NULL && *nsplicepairs < maxchimerapaths) {
      donor = (Substring_T) p->first;
      acceptor = (Substring_T) q->first;
      debug4ld(printf("donor at %u and acceptor at %u\n",Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	
      if ((pos = Substring_chimera_pos(donor)) < min_endlength_2) {
	p = p->rest;
      } else if (pos > querylength - min_endlength_1) {
	p = p->rest;
      } else if (pos < Substring_chimera_pos(acceptor)) {
	p = p->rest;
      } else if (pos > Substring_chimera_pos(acceptor)) {
	q = q->rest;
      } else {
	qsave = q;
	while (p != NULL && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) + pos) > (Substring_genomicstart(donor) - pos)) {
	      distance = (Substring_genomicstart(acceptor) + pos) - (Substring_genomicstart(donor) - pos);
	    } else {
	      distance = (Substring_genomicstart(donor) - pos) - (Substring_genomicstart(acceptor) + pos);
	    }
	    debug4ld(printf("7-6. Pushing a candidate at chimera_pos %d (%d..%d), donor %u to acceptor %u.  Different strands so not shortdistance.\n",
			    pos,min_endlength_2,querylength-min_endlength_1,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3_new_splice(&(*found_score),donor,acceptor,distance,
								   /*shortdistancep*/false,distantsplicing_penalty,/*support*/querylength,
								   querylength,/*copyp*/true,/*sensep*/false));
	    (*nsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }
  }
  
  return distantsplicing;
}


static List_T
find_splices_halfintrons (int *found_score,
			  List_T *donors_plus, List_T *antidonors_plus,
			  List_T *acceptors_plus, List_T *antiacceptors_plus,
			  List_T *donors_minus, List_T *antidonors_minus,
			  List_T *acceptors_minus, List_T *antiacceptors_minus, 
			  int max_mismatches_allowed, int splicing_penalty,
			  int querylength) {
  List_T hits = NULL, p;
  Substring_T donor, acceptor;
  int nmismatches, support, endlength;
  
  debug(printf("Starting find_splices_halfintrons\n"));
  debug(
	for (nmismatches = 0; nmismatches <= max_mismatches_allowed; nmismatches++) {
	  printf("At %d nmismatches: +donors/acceptors %d/%d, +antidonors/antiacceptors %d/%d, -donors/acceptors %d/%d, -antidonors/antiacceptors %d/%d\n",
		 nmismatches,
		 List_length((donors_plus)[nmismatches]),List_length((acceptors_plus)[nmismatches]),
		 List_length((antidonors_plus)[nmismatches]),List_length((antiacceptors_plus)[nmismatches]),
		 List_length((donors_minus)[nmismatches]),List_length((acceptors_minus)[nmismatches]),
		 List_length((antidonors_minus)[nmismatches]),List_length((antiacceptors_minus)[nmismatches]));
	});
  
  /* Donors and antiacceptors => Want chimera_pos to be at end */
  /* Acceptors and antidonors => Want chimera_pos to be at beginning */
  
  nmismatches = 0;
  /* Don't want to end when first set of hits found */
  while (/* hits == NULL && */ nmismatches <= max_mismatches_allowed) {
    /* End 1 */
    for (p = donors_plus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      
      support = Substring_chimera_pos(donor);
      /* endlength = querylength - support; */
#if 0
      /* Use same standard as sufficient_splice_prob_distant */
      if (sufficient_splice_prob_halfintron(support,Substring_nmismatches(donor),Substring_chimera_prob(donor))) {
#endif
	debug4h(printf("halfintron donor_plus: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(donor),Substring_chimera_pos(donor),Substring_nmismatches(donor)));
	hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),donor,/*acceptor*/NULL,/*distance*/0U,
							 /*shortdistancep*/false,splicing_penalty,support,
							 querylength,/*copyp*/true,/*sensep*/true));
#if 0
      } else {
	debug4h(printf("halfintron donor_plus, insufficient support: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(donor),Substring_chimera_pos(donor),Substring_nmismatches(donor)));
      }
#endif
    }
  
    /* End 2 */
    for (p = acceptors_plus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      
      endlength = Substring_chimera_pos(acceptor);
      support = querylength - endlength;
#if 0
      /* Use same standard as sufficient_splice_prob_distant */
      if (sufficient_splice_prob_halfintron(support,Substring_nmismatches(acceptor),Substring_chimera_prob(acceptor))) {
#endif
	debug4h(printf("halfintron acceptor_plus: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(acceptor),Substring_chimera_pos(acceptor),Substring_nmismatches(acceptor)));
	hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),/*donor*/NULL,acceptor,/*distance*/0U,
							 /*shortdistancep*/false,splicing_penalty,support,
							 querylength,/*copyp*/true,/*sensep*/true));
#if 0
      } else {
	debug4h(printf("halfintron acceptor_plus, insufficient support: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(acceptor),Substring_chimera_pos(acceptor),Substring_nmismatches(acceptor)));
      }
#endif
    }
  
    /* End 3 */
    for (p = donors_minus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      
      support = Substring_chimera_pos(donor);
      /* endlength = querylength - support; */
#if 0
      /* Use same standard as sufficient_splice_prob_distant */
      if (sufficient_splice_prob_halfintron(support,Substring_nmismatches(donor),Substring_chimera_prob(donor))) {
#endif
	debug4h(printf("halfintron donor_minus: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(donor),Substring_chimera_pos(donor),Substring_nmismatches(donor)));
	hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),donor,/*acceptor*/NULL,/*distance*/0U,
							 /*shortdistancep*/false,splicing_penalty,support,
							 querylength,/*copyp*/true,/*sensep*/true));
#if 0
      } else {
	debug4h(printf("halfintron donor_minus, insufficient support: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(donor),Substring_chimera_pos(donor),Substring_nmismatches(donor)));
      }
#endif
    }
  
    /* End 4 */
    for (p = acceptors_minus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      
      endlength = Substring_chimera_pos(acceptor);
      support = querylength - endlength;
#if 0
      /* Use same standard as sufficient_splice_prob_distant */
      if (sufficient_splice_prob_halfintron(support,Substring_nmismatches(acceptor),Substring_chimera_prob(acceptor))) {
#endif
	debug4h(printf("halfintron acceptor_minus: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(acceptor),Substring_chimera_pos(acceptor),Substring_nmismatches(acceptor)));
	hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),/*donor*/NULL,acceptor,/*distance*/0U,
							 /*shortdistancep*/false,splicing_penalty,support,
							 querylength,/*copyp*/true,/*sensep*/true));
#if 0
      } else {
	debug4h(printf("halfintron acceptor_minus, insufficient support: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(acceptor),Substring_chimera_pos(acceptor),Substring_nmismatches(acceptor)));
      }
#endif
    }
  
    /* End 5 */
    for (p = antidonors_plus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      
      endlength = Substring_chimera_pos(donor);
      support = querylength - endlength;
#if 0
      /* Use same standard as sufficient_splice_prob_distant */
      if (sufficient_splice_prob_halfintron(support,Substring_nmismatches(donor),Substring_chimera_prob(donor))) {
#endif
	debug4h(printf("halfintron antidonor_plus: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(donor),Substring_chimera_pos(donor),Substring_nmismatches(donor)));
	hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),donor,/*acceptor*/NULL,/*distance*/0U,
							 /*shortdistancep*/false,splicing_penalty,support,
							 querylength,/*copyp*/true,/*sensep*/false));
#if 0
      } else {
	debug4h(printf("halfintron antidonor_plus, insufficient support: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(donor),Substring_chimera_pos(donor),Substring_nmismatches(donor)));
      }
#endif
    }
  
    /* End 6 */
    for (p = antiacceptors_plus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      
      support = Substring_chimera_pos(acceptor);
      /* endlength = querylength - support; */
#if 0
      /* Use same standard as sufficient_splice_prob_distant */
      if (sufficient_splice_prob_halfintron(support,Substring_nmismatches(acceptor),Substring_chimera_prob(acceptor))) {
#endif
	debug4h(printf("halfintron antiacceptor_plus: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(acceptor),Substring_chimera_pos(acceptor),Substring_nmismatches(acceptor)));
	hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),/*donor*/NULL,acceptor,/*distance*/0U,
							 /*shortdistancep*/false,splicing_penalty,support,
							 querylength,/*copyp*/true,/*sensep*/false));
#if 0
      } else {
	debug4h(printf("halfintron antiacceptor_plus, insufficient support: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(acceptor),Substring_chimera_pos(acceptor),Substring_nmismatches(acceptor)));
      }
#endif
    }
  
    /* End 7 */
    for (p = antidonors_minus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      
      endlength = Substring_chimera_pos(donor);
      support = querylength - endlength;
#if 0
      /* Use same standard as sufficient_splice_prob_distant */
      if (sufficient_splice_prob_halfintron(support,Substring_nmismatches(donor),Substring_chimera_prob(donor))) {
#endif
	debug4h(printf("halfintron antidonor_minus: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(donor),Substring_chimera_pos(donor),Substring_nmismatches(donor)));
	hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),donor,/*acceptor*/NULL,/*distance*/0U,
							 /*shortdistancep*/false,splicing_penalty,support,
							 querylength,/*copyp*/true,/*sensep*/false));
#if 0
      } else {
	debug4h(printf("halfintron antidonor_minus, insufficient support: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(donor),Substring_chimera_pos(donor),Substring_nmismatches(donor)));
      }
#endif
    }
  
    /* End 8 */
    for (p = antiacceptors_minus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      
      support = Substring_chimera_pos(acceptor);
      /* endlength = querylength - support; */
#if 0
      /* Use same standard as sufficient_splice_prob_distant */
      if (sufficient_splice_prob_halfintron(support,Substring_nmismatches(acceptor),Substring_chimera_prob(acceptor))) {
#endif
	debug4h(printf("halfintron antiacceptor_minus: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(acceptor),Substring_chimera_pos(acceptor),Substring_nmismatches(acceptor)));
	hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),/*donor*/NULL,acceptor,/*distance*/0U,
							 /*shortdistancep*/false,splicing_penalty,support,
							 querylength,/*copyp*/true,/*sensep*/false));
#if 0
      } else {
	debug4h(printf("halfintron antiacceptor_minus, insufficient support: %f at %d (%d mismatches)\n",
		      Substring_chimera_prob(acceptor),Substring_chimera_pos(acceptor),Substring_nmismatches(acceptor)));
      }
#endif
    }

    nmismatches++;
  }
  
  hits = Stage3_remove_duplicates(hits);

  return hits;
}


static bool
find_12mer_bounds (int *firstbound, int *lastbound, bool *omitted, int query_lastpos) {
  int nconsecutive;
  int querypos;

  querypos = 0;
  nconsecutive = 0;
  while (nconsecutive < 3 && querypos <= query_lastpos) {
    if (omitted[querypos]) {
      nconsecutive = 0;
    } else {
      nconsecutive++;
    }
    querypos++;
  }
  if (nconsecutive < 3) {
    *firstbound = 0;
    *lastbound = query_lastpos;
    debug(printf("nconsecutive from left < 3, so setting firstbound 0 and lastbound %d\n",query_lastpos));
    return false;
  } else {
    *firstbound = querypos+2 + INDEX1PART; /* From trial-and-error, this is the correct value */
    debug(printf("Assigning firstbound to be querypos %d + 2 + INDEX1PART = %d\n",
		 querypos,*firstbound));
  }

  querypos = query_lastpos;
  nconsecutive = 0;
  while (nconsecutive < 3 && querypos >= 0) {
    if (omitted[querypos]) {
      nconsecutive = 0;
    } else {
      nconsecutive++;
    }
    querypos--;
  }

  if (nconsecutive < 3) {
    *firstbound = 0;
    *lastbound = query_lastpos;
    debug(printf("nconsecutive from right < 3, so setting firstbound 0 and lastbound %d\n",query_lastpos));
    return false;
  } else {
    *lastbound = querypos-1; /* From trial-and-error, this is the correct value */
    debug(printf("Assigning lastbound to be querypos %d - 1 = %d\n",
		 querypos,*lastbound));
  }

  if (*firstbound > query_lastpos) {
    debug(printf("firstbound %d > query_lastpos %d, so setting firstbound 0 and lastbound %d\n",
		 *firstbound,query_lastpos,query_lastpos));
    *firstbound = 0;
    *lastbound = query_lastpos;
    return false;
#if 0
  } else if (*firstbound - INDEX1PART > *lastbound) {
    *firstbound = 0;
    *lastbound = query_lastpos;
    return false;
#endif
  } else if (*lastbound <= INDEX1PART) {
    debug(printf("lastbound %d <= %d, so setting firstbound 0 and lastbound %d\n",
		 *lastbound,INDEX1PART,query_lastpos));
    *firstbound = 0;
    *lastbound = query_lastpos;
    return false;
  } else {
    return true;
  }
}



static void
complete_set_mm_indels (int *found_score, bool *any_omitted_p, int *max_level, int *done_level, bool revise_levels_p,
			List_T *subs, List_T *indels, T this,
			char *queryuc_ptr, char *queryrc, int querylength, int query_lastpos,
			char *gbuffer, Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
			IIT_T chromosome_iit, UINT4 *genome_blocks, UINT4 *snp_blocks,
			Genome_T genome, Floors_T *floors_array,
			int subopt_levels, int indel_penalty, int max_middle_insertions, int max_middle_deletions,
			bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
			bool trim_ends_p, bool dibasep, bool cmetp,
			int fast_level, bool omit_frequent_p, bool omit_repetitive_p) {
  struct Indelsegment_T *plus_segments = NULL, *minus_segments = NULL;
  int plus_nsegments, minus_nsegments;
  Floors_T floors;
  int indel_level;
  int max_indel_sep;
  bool alloc_floors_p = false, all_omitted_p;
  int firstbound, lastbound;

  debug(printf("Starting complete_set_mm_indels\n"));

  /* 4 and 5. Mismatches and indels via complete set.  Requires compress and
     all positions fetched.  Omits oligos and creates segments for some diagonals. */
  if (this->query_compress_fwd == NULL) {
    this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
    this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
  }
  if (this->all_positions_fetched_p == true) {
    omit_oligos_clear(this,query_lastpos);
  } else {
    fetch_positions_for_all_12mers(this,indexdb,indexdb2,query_lastpos);
  }

  debug(printf("Omitting frequent/repetitive oligos\n"));
  omit_oligos(&all_omitted_p,&(*any_omitted_p),this,query_lastpos,indexdb_size_threshold,
	      omit_frequent_p,omit_repetitive_p);
  
  if (all_omitted_p == true) {
    debug(printf("Aborting because all oligos are omitted\n"));
    return;
  } else if (*any_omitted_p) {
    floors = Floors_new_omitted(querylength,max_end_insertions,this->omitted);
    alloc_floors_p = true;
  } else {
    if (floors_array[querylength] == NULL) {
      floors_array[querylength] = Floors_new_standard(querylength,max_end_insertions);
    }
    floors = floors_array[querylength];
  }
  
  if (find_12mer_bounds(&firstbound,&lastbound,this->omitted,query_lastpos) == false) {
    debug(printf("Cannot find 12_mer bounds\n"));
    allow_end_indels_p = false;
  } else {
    debug(printf("Found firstbound %d and lastbound %d\n",firstbound,lastbound));
  }

#if 0
  if (*subs != NULL) {
    indels_prior_penalty = 1;
  }

  if (*done_level < indel_penalty + indels_prior_penalty) {
    /* Prevents accumulation of segments for indels */
    max_indel_sep = 0;
    allow_end_indels_p = false;
  } else {
    max_indel_sep = (max_middle_insertions > max_middle_deletions) ? 
      max_middle_insertions + 1 : max_middle_deletions + 1;
  }
#else
  if (*done_level < indel_penalty) {
    /* Prevents accumulation of segments for indels */
    max_indel_sep = 0;
    allow_end_indels_p = false;
  } else {
    max_indel_sep = (max_middle_insertions > max_middle_deletions) ? 
      max_middle_insertions + 1 : max_middle_deletions + 1;
  }
#endif

  /* 4. Complete set mismatches */
  /* Done as a single batch */
  debug(printf("*** Stage 5.  Complete set mismatches up to %d ***\n",(*done_level < fast_level) ? -1 : *done_level));
  *subs = find_complete_mm_and_indelsegments(&(*found_score),*subs,&plus_segments,&plus_nsegments,
					     this->plus_positions,this->plus_npositions,
					     this->omitted,querylength,query_lastpos,
					     /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,gbuffer,/*query_compress*/this->query_compress_fwd,
					     genome_blocks,snp_blocks,genome,chromosome_iit,
					     floors,/*max_mismatches_allowed*/(*done_level < fast_level ? -1 : *done_level),
					     /*indel_mismatches_allowed*/*done_level - indel_penalty /* - indels_prior_penalty*/,
					     max_indel_sep,max_end_insertions,firstbound,lastbound,
					     trim_ends_p,dibasep,cmetp,/*plusp*/true);
  
  *subs = find_complete_mm_and_indelsegments(&(*found_score),*subs,&minus_segments,&minus_nsegments,
					     this->minus_positions,this->minus_npositions,
					     this->omitted,querylength,query_lastpos,
					     /*query*/queryuc_ptr,/*queryptr*/queryrc,gbuffer,/*query_compress*/this->query_compress_rev,
					     genome_blocks,snp_blocks,genome,chromosome_iit,
					     floors,/*max_mismatches_allowed*/(*done_level < fast_level ? -1 : *done_level),
					     /*indel_mismatches_allowed*/*done_level - indel_penalty /* - indels_prior_penalty*/,
					     max_indel_sep,max_end_insertions,firstbound,lastbound,
					     trim_ends_p,dibasep,cmetp,/*plusp*/false);
  if (revise_levels_p == true) {
    *max_level = (*found_score < *max_level) ? *found_score : *max_level;
    *done_level = *max_level + subopt_levels;
  }
  debug(printf("5> found_score = %d, max_level %d, done_level %d\n",*found_score,*max_level,*done_level));
  
  debug(printf("plus_nsegments = %d, minus_nsegments = %d\n",plus_nsegments,minus_nsegments));

#if 0
  if (*subs != NULL) {
    indels_prior_penalty = 1;
  }
#endif

  if (*done_level >= indel_penalty /* + indels_prior_penalty*/) {
    /* 5. Indels */
    /* Need to reverse, because middle indelsplicing procedure depends on ascending diagonal order */

    /* Done iteratively */
    indel_level = indel_penalty /* + indels_prior_penalty*/;
    while (indel_level <= *done_level) {
      debug(printf("*** Stage 6.  Middle indels with %d-%d mismatches allowed\n",indel_level,indel_penalty));
      *indels = find_middle_indels(&(*found_score),*indels,plus_segments,minus_segments,
				   plus_nsegments,minus_nsegments,genome,genome_blocks,snp_blocks,
				   queryuc_ptr,queryrc,floors,
				   querylength,query_lastpos,this->query_compress_fwd,this->query_compress_rev,
				   max_middle_insertions,max_middle_deletions,min_indel_end_matches,indel_penalty,
				   /*indel_mismatches_allowed*/indel_level - indel_penalty /* - indels_prior_penalty*/,
				   trim_ends_p,dibasep,cmetp);
      
      if (allow_end_indels_p == true) {
	debug(printf("*** Stage 6.  End indels with %d-%d mismatches allowed\n",indel_level,indel_penalty));
	*indels = find_end_indels(&(*found_score),*indels,plus_segments,minus_segments,
				  plus_nsegments,minus_nsegments,genome,genome_blocks,snp_blocks,
				  queryuc_ptr,queryrc,querylength,
				  firstbound,lastbound,this->query_compress_fwd,this->query_compress_rev,
				  max_end_insertions,max_end_deletions,min_indel_end_matches,indel_penalty,
				  /*indel_mismatches_allowed*/indel_level - indel_penalty /* - indels_prior_penalty*/,
				  trim_ends_p,dibasep,cmetp);
      }
      if (revise_levels_p == true) {
	*max_level = (*found_score < *max_level) ? *found_score : *max_level;
	*done_level = *max_level + subopt_levels;
      }
      indel_level++;
      debug(printf("6> found_score = %d, max_level %d, done_level %d\n",*found_score,*max_level,*done_level));
    }
  }

  FREE(minus_segments);
  FREE(plus_segments);

  if (alloc_floors_p) {
    Floors_free(&floors);
  }

  debug(printf("Finished with complete_set_mm_indels\n"));

  return;
}


static List_T
complete_set_localsplicing (int *found_score, bool *any_omitted_p, Floors_T *floors, bool *alloc_floors_p,
			    struct Splicesegment_T **plus_segments, int *plus_nsegments, 
			    struct Splicesegment_T **minus_segments, int *minus_nsegments,
			    List_T localsplicing, T this, char *queryuc_ptr, int querylength, int query_lastpos,
			    Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
			    IIT_T chromosome_iit, UINT4 *genome_blocks, UINT4 *snp_blocks,
			    Genome_T genome, Floors_T *floors_array, int max_end_insertions,
			    Genomicpos_T *splicesites, Splicetype_T *splicetypes, int nsplicesites, Genomicpos_T shortsplicedist,
			    bool novelsplicingp, int localsplicing_penalty,
			    int min_localsplicing_end_matches, int max_mismatches_allowed,
			    bool trim_ends_p, bool dibasep, bool cmetp, bool omit_frequent_p, bool omit_repetitive_p) {
  bool all_omitted_p;

  debug(printf("Starting complete_set_localsplicing with %d mismatches allowed\n",max_mismatches_allowed));

  if (this->query_compress_fwd == NULL) {
    this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
    this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
  }
  if (this->all_positions_fetched_p == true) {
    omit_oligos_clear(this,query_lastpos);
  } else {
    fetch_positions_for_all_12mers(this,indexdb,indexdb2,query_lastpos);
  }

  debug(printf("Omitting frequent/repetitive oligos\n"));
  omit_oligos(&all_omitted_p,&(*any_omitted_p),this,query_lastpos,indexdb_size_threshold,
	      omit_frequent_p,omit_repetitive_p);

  if (all_omitted_p == true) {
    debug(printf("Aborting because all oligos are omitted\n"));
    *alloc_floors_p = false;
    return (List_T) NULL;
  } else if (*any_omitted_p) {
    *alloc_floors_p = true;
    *floors = Floors_new_omitted(querylength,max_end_insertions,this->omitted);
  } else {
    *alloc_floors_p = false;
    if (floors_array[querylength] == NULL) {
      floors_array[querylength] = Floors_new_standard(querylength,max_end_insertions);
    }
    *floors = floors_array[querylength];
  }

   *plus_segments = identify_segments_all(&(*plus_nsegments),this->plus_positions,this->plus_npositions,
					  this->forward_oligos,this->omitted,querylength,query_lastpos,
					  chromosome_iit,/*plusp*/true);

   *minus_segments = identify_segments_all(&(*minus_nsegments),this->minus_positions,this->minus_npositions,
					   this->revcomp_oligos,this->omitted,querylength,query_lastpos,
					   chromosome_iit,/*plusp*/false);

   localsplicing = find_splicepairs_local_plus(&(*found_score),localsplicing,*plus_segments,*plus_nsegments,
					       genome_blocks,snp_blocks,genome,queryuc_ptr,
					       *floors,querylength,query_lastpos,/*query_compress*/this->query_compress_fwd,
					       splicesites,splicetypes,nsplicesites,novelsplicingp,
					       /*max_distance*/shortsplicedist,
					       /*splicing_penalty*/localsplicing_penalty,min_localsplicing_end_matches,
					       max_mismatches_allowed,trim_ends_p,dibasep,cmetp);

   localsplicing = find_splicepairs_local_minus(&(*found_score),localsplicing,*minus_segments,*minus_nsegments,
						genome_blocks,snp_blocks,genome,queryuc_ptr,
						*floors,querylength,query_lastpos,/*query_compress*/this->query_compress_rev,
						splicesites,splicetypes,nsplicesites,novelsplicingp,
						/*max_distance*/shortsplicedist,
						/*splicing_penalty*/localsplicing_penalty,min_localsplicing_end_matches,
						max_mismatches_allowed,trim_ends_p,dibasep,cmetp);

   debug(printf("Finished with complete_set_localsplicing\n"));

   return localsplicing;
}



/* done_level should probably be renamed final_level.  max_level
   should probably be renamed found_level or first_level. */
static List_T
align_end (int *cutoff_level, T this, char *queryuc_ptr, int querylength, int query_lastpos,
	   Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
	   IIT_T chromosome_iit, UINT4 *genome_blocks, UINT4 *snp_blocks,
	   Genome_T genome, Genome_T genomealt, Floors_T *floors_array,
	   Genomicpos_T *splicesites, Splicetype_T *splicetypes, int nsplicesites,
	   int usermax_level, int subopt_levels, Masktype_T masktype, int indel_penalty, 
	   int maxpaths, int maxchimerapaths, bool knownsplicingp, bool novelsplicingp, 
	   bool canonicalp, Genomicpos_T shortsplicedist,
	   int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
	   int min_distantsplicing_end_matches, double min_distantsplicing_identity,
	   int max_middle_insertions, int max_middle_deletions,
	   bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
	   bool allvalidp, bool trim_ends_p, bool dibasep, bool cmetp) {
  List_T subs = NULL, indels = NULL, localsplicing = NULL, halfintrons = NULL, distantsplicing = NULL;
  struct Splicesegment_T *plus_segments = NULL, *minus_segments = NULL;
  int plus_nsegments, minus_nsegments;
  char queryrc[MAX_QUERYLENGTH+1];
  int found_score, done_level, max_level, fast_level, mismatch_level, nmismatches;
  int nsplicepairs = 0;
  int max_splice_mismatches = -1, i;
  List_T *donors_plus, *antidonors_plus, *acceptors_plus, *antiacceptors_plus,
    *donors_minus, *antidonors_minus, *acceptors_minus, *antiacceptors_minus;
  bool any_omitted_p, alloc_floors_p = false;
  Floors_T floors;

  char gbuffer[MAX_QUERYLENGTH+1];
  gbuffer[querylength] = '\0';

  make_complement_buffered(queryrc,queryuc_ptr,querylength);
  found_score = querylength;
  fast_level = (querylength + 2)/INDEX1PART - NREQUIRED_FAST;
  if (fast_level < 1 && usermax_level < 0) {
    fast_level = 1;		/* Do at least 1 mismatch */
  }

  if (usermax_level < 0) {
    *cutoff_level = fast_level;
  } else {
    *cutoff_level = usermax_level;
  }

  if (dibasep) {
    max_level = querylength;	/* Allow extra because color errors may exceed nt errors */
  } else if (usermax_level < 0) {
    max_level = fast_level;
  } else {
    max_level = usermax_level;
  }
  done_level = max_level + subopt_levels;
  debug(printf("0> max_level %d, done_level %d\n",max_level,done_level));

  /* 1. Exact.  Requires compress if cmet or genomealt.  Creates and uses spanning set. */
  mismatch_level = 0;
  if (allvalidp == false) {
    debug(printf("Not all oligos are valid, so cannot perform spanning set\n"));
    fast_level = -1;
  } else {
    debug(printf("fast_level = %d\n",fast_level));
    debug(printf("*** Stage 1.  Exact ***\n"));
    if (cmetp || genomealt) {
      this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
      this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
    }
    subs = find_spanning_exact_matches(&found_score,this,queryuc_ptr,queryrc,querylength,query_lastpos,
				       indexdb,indexdb2,gbuffer,this->query_compress_fwd,this->query_compress_rev,
				       genome_blocks,snp_blocks,genome,genomealt,chromosome_iit,trim_ends_p,dibasep,cmetp);
    max_level = (found_score < max_level) ? found_score : max_level;
    done_level = max_level + subopt_levels;
    mismatch_level = 1;
    debug(printf("1> found_score = %d, max_level %d, done_level %d\n",found_score,max_level,done_level));
  }
    
  if (dibasep == false) {
    /* 2. One mismatch.  Requires spanning set and compress. */
    if (allvalidp && querylength >= 22 && done_level >= 1) {
      debug(printf("*** Stage 2.  One miss ***\n"));
      if (this->query_compress_fwd == NULL) {
	this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
	this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
      }
      subs = find_spanning_onemiss_matches(&found_score,subs,this,queryuc_ptr,queryrc,querylength,
					   gbuffer,this->query_compress_fwd,this->query_compress_rev,
					   genome_blocks,snp_blocks,genome,genomealt,chromosome_iit,trim_ends_p,dibasep,cmetp);
      max_level = (found_score < max_level) ? found_score : max_level;
      done_level = max_level + subopt_levels;
      mismatch_level = 2;
      debug(printf("2> found_score = %d, max_level %d, done_level %d\n",found_score,max_level,done_level));
    }
  
    /* 3. Mismatches via spanning set.  Requires spanning set and compress. */
    if (allvalidp && done_level >= 2) {
      if (this->query_compress_fwd == NULL) {
	this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
	this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
      }
      while (mismatch_level <= fast_level && mismatch_level <= done_level) {
	debug(printf("*** Stage 3 (level %d).  Spanning set mismatches ***\n",mismatch_level));
	subs = find_spanning_multimiss_matches(&found_score,subs,this,NREQUIRED_FAST,
					       queryuc_ptr,queryrc,querylength,
					       gbuffer,this->query_compress_fwd,this->query_compress_rev,
					       genome_blocks,snp_blocks,genome,genomealt,chromosome_iit,
					       /*nmisses_allowed*/mismatch_level,trim_ends_p,dibasep,cmetp);
	max_level = (found_score < max_level) ? found_score : max_level;
	done_level = max_level + subopt_levels;
	mismatch_level++;
	debug(printf("3> found_score = %d, max_level %d, done_level %d\n",found_score,max_level,done_level));
      }
    }
  }

  /* 4, 5.  Complete set mismatches and indels, omitting frequent oligos */
#if 0
  if (subs) {
    indels_prior_penalty = 1;
  }
#endif
  if (dibasep || done_level > fast_level || done_level >= indel_penalty /* + indels_prior_penalty*/) {
    if (masktype == MASK_NONE) {
      debug(printf("*** Stage 4,5.  Complete mm/indels with no masking with done_level %d ***\n",done_level));
      complete_set_mm_indels(&found_score,&any_omitted_p,&max_level,&done_level,/*revise_levels_p*/true,
			     &subs,&indels,this,queryuc_ptr,queryrc,
			     querylength,query_lastpos,gbuffer,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,genome,floors_array,
			     subopt_levels,indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     trim_ends_p,dibasep,cmetp,fast_level,
			     /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
    } else {
      debug(printf("*** Stage 4,5.  Complete mm/indels masking frequent oligos with done_level %d ***\n",done_level));
      complete_set_mm_indels(&found_score,&any_omitted_p,&max_level,&done_level,/*revise_levels_p*/true,
			     &subs,&indels,this,queryuc_ptr,queryrc,
			     querylength,query_lastpos,gbuffer,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,genome,floors_array,
			     subopt_levels,indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     trim_ends_p,dibasep,cmetp,fast_level,/*omit_frequent_p*/true,
			     /*omit_repetitive_p*/(masktype == MASK_REPETITIVE || masktype == MASK_GREEDY_REPETITIVE) ? true : false
			     );
      if ((masktype == MASK_GREEDY_FREQUENT || masktype == MASK_GREEDY_REPETITIVE) && subs == NULL && indels == NULL && any_omitted_p == true) {
	debug(printf("*** Stage 4,5.  Complete mm/indels with no masking with done_level %d ***\n",done_level));
	complete_set_mm_indels(&found_score,&any_omitted_p,&max_level,&done_level,/*revise_levels_p*/true,
			       &subs,&indels,this,queryuc_ptr,queryrc,
			       querylength,query_lastpos,gbuffer,indexdb,indexdb2,indexdb_size_threshold,
			       chromosome_iit,genome_blocks,snp_blocks,genome,floors_array,
			       subopt_levels,indel_penalty,max_middle_insertions,max_middle_deletions,
			       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			       trim_ends_p,dibasep,cmetp,fast_level,
			       /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
      }
    }
  }

  /* 6, 7, 8.  Splicing.  Requires compress and all positions fetched */
  if (knownsplicingp || novelsplicingp) {
    /* 6.  Local splicing */
#if 0
    if (subs || indels) {
      localsplicing_prior_penalty = 1;
    }
#endif

    debug(printf("Deciding whether to do localsplicing: done_level %d >=? localsplicing_penalty %d\n",
		 done_level,localsplicing_penalty));
    if (done_level >= localsplicing_penalty /* + localsplicing_prior_penalty*/) {
      debug(printf("*** Stage 6.  Local splicing masking frequent oligos with done_level %d ***\n",done_level));
      /* Always mask frequent oligos for splicing, which must be transcriptional */
      localsplicing = complete_set_localsplicing(&found_score,&any_omitted_p,&floors,&alloc_floors_p,
						 &plus_segments,&plus_nsegments,&minus_segments,&minus_nsegments,
						 localsplicing,this,queryuc_ptr,
						 querylength,query_lastpos,indexdb,indexdb2,indexdb_size_threshold,
						 chromosome_iit,genome_blocks,snp_blocks,genome,
						 floors_array,max_end_insertions,
						 splicesites,splicetypes,nsplicesites,shortsplicedist,novelsplicingp,
						 localsplicing_penalty,min_localsplicing_end_matches,
						 /*max_mismatches_allowed*/done_level - localsplicing_penalty,
						 trim_ends_p,dibasep,cmetp,/*omit_frequent_p*/true,/*omit_repetitive_p*/true);
      max_level = (found_score < max_level) ? found_score : max_level;
      done_level = max_level + subopt_levels;
    }

    /* 7, 8.  Distant and half splicing */
#if 0
    if (subs || indels || localsplicing) {
      distantsplicing_prior_penalty = 1;
    }
#endif

    max_splice_mismatches = done_level - distantsplicing_penalty /* - distantsplicing_prior_penalty*/;
    if (localsplicing == NULL /* && (knownsplicingp == true || querylength >= novel_half_intron_min_support) */) {
      if (done_level - localsplicing_penalty /* - localsplicing_prior_penalty*/ > max_splice_mismatches) {
	max_splice_mismatches = done_level - localsplicing_penalty /* - localsplicing_prior_penalty*/;
      }
    }

    debug(printf("max_splice_mismatches %d = done_level %d - distantsplicing_penalty %d\n",
		 max_splice_mismatches,done_level,distantsplicing_penalty));
    if (max_splice_mismatches >= 0) {
      donors_plus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      antidonors_plus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      acceptors_plus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      antiacceptors_plus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      donors_minus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      antidonors_minus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      acceptors_minus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      antiacceptors_minus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));

      if (plus_segments) {
	debug(printf("Starting identify_spliceends_plus\n"));
	identify_spliceends_plus(&donors_plus,&antidonors_plus,&acceptors_plus,&antiacceptors_plus,
				 plus_segments,plus_nsegments,genome_blocks,snp_blocks,genome,
				 queryuc_ptr,floors,querylength,query_lastpos,gbuffer,
				 /*query_compress*/this->query_compress_fwd,
				 splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
				 /*max_mismatches_allowed*/max_splice_mismatches,trim_ends_p,dibasep,cmetp);
	debug(printf("Finished identify_spliceends_plus\n"));
      }
      if (minus_segments) {
	debug(printf("Starting identify_spliceends_minus\n"));
	identify_spliceends_minus(&donors_minus,&antidonors_minus,&acceptors_minus,&antiacceptors_minus,
				  minus_segments,minus_nsegments,genome_blocks,snp_blocks,genome,
				  queryuc_ptr,queryrc,floors,querylength,query_lastpos,gbuffer,
				  /*query_compress*/this->query_compress_rev,
				  splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
				  /*max_mismatches_allowed*/max_splice_mismatches,trim_ends_p,dibasep,cmetp);
	debug(printf("Finished identify_spliceends_minus\n"));
      }

      /* 7.  Find distant splicing iteratively using both known and novel splice sites */
      nmismatches = 0;
      while (nmismatches <= done_level - distantsplicing_penalty /* - distantsplicing_prior_penalty*/ && 
	     nsplicepairs < maxchimerapaths) {
	debug(printf("*** Stage 7.  Distant splicing, allowing %d mismatches ***\n",nmismatches));

	debug4(printf("Sorting splice ends\n"));
	donors_plus[nmismatches] = Substring_sort_chimera_halves(donors_plus[nmismatches],/*ascendingp*/true);
	acceptors_plus[nmismatches] = Substring_sort_chimera_halves(acceptors_plus[nmismatches],/*ascendingp*/true);

	antidonors_plus[nmismatches] = Substring_sort_chimera_halves(antidonors_plus[nmismatches],/*ascendingp*/false);
	antiacceptors_plus[nmismatches] = Substring_sort_chimera_halves(antiacceptors_plus[nmismatches],/*ascendingp*/false);
	  
	donors_minus[nmismatches] = Substring_sort_chimera_halves(donors_minus[nmismatches],/*ascendingp*/false);
	acceptors_minus[nmismatches] = Substring_sort_chimera_halves(acceptors_minus[nmismatches],/*ascendingp*/false);

	antidonors_minus[nmismatches] = Substring_sort_chimera_halves(antidonors_minus[nmismatches],/*ascendingp*/true);
	antiacceptors_minus[nmismatches] = Substring_sort_chimera_halves(antiacceptors_minus[nmismatches],/*ascendingp*/true);

	debug4(printf("Splice ends at %d nmismatches: +donors/acceptors %d/%d, +antidonors/antiacceptors %d/%d, -donors/acceptors %d/%d, -antidonors/antiacceptors %d/%d\n",
		      nmismatches,
		      List_length(donors_plus[nmismatches]),List_length(acceptors_plus[nmismatches]),
		      List_length(antidonors_plus[nmismatches]),List_length(antiacceptors_plus[nmismatches]),
		      List_length(donors_minus[nmismatches]),List_length(acceptors_minus[nmismatches]),
		      List_length(antidonors_minus[nmismatches]),List_length(antiacceptors_minus[nmismatches])));

	distantsplicing = find_splicepairs_distant(&found_score,&nsplicepairs,distantsplicing,
						   donors_plus,antidonors_plus,acceptors_plus,antiacceptors_plus,
						   donors_minus,antidonors_minus,acceptors_minus,antiacceptors_minus,
						   shortsplicedist,localsplicing_penalty,distantsplicing_penalty,
						   min_distantsplicing_end_matches,min_distantsplicing_identity,
						   querylength,nmismatches,maxchimerapaths);
	if (distantsplicing) {
	  max_level = (found_score < max_level) ? found_score : max_level;
	  done_level = max_level + subopt_levels;
	  debug(printf("7> found_score = %d, max_level %d, done_level %d\n",found_score,max_level,done_level));
	}
	nmismatches++;
      }

      /* 8.  Half introns */
      debug(printf("%d local splices, %d distant splices\n",
		   List_length(localsplicing),List_length(distantsplicing)));
      if (localsplicing == NULL && distantsplicing == NULL) {
	debug(printf("*** Stage 8.  Half introns up to %d mismatches ***\n",max_splice_mismatches));
	halfintrons = find_splices_halfintrons(&found_score,donors_plus,antidonors_plus,acceptors_plus,antiacceptors_plus,
					       donors_minus,antidonors_minus,acceptors_minus,antiacceptors_minus,
					       /*max_mismatches*/max_splice_mismatches,/*splicing_penalty*/localsplicing_penalty,
					       querylength);
#if 0
	max_level = (found_score < max_level) ? found_score : max_level;
	done_level = max_level + subopt_levels;
	debug(printf("8> found_score = %d, max_level %d, done_level %d\n",found_score,max_level,done_level));
#endif
      }

      for (i = 0; i <= max_splice_mismatches; i++) {
	substringlist_gc(&(donors_plus[i]));
	substringlist_gc(&(antidonors_plus[i]));
	substringlist_gc(&(acceptors_plus[i]));
	substringlist_gc(&(antiacceptors_plus[i]));
	substringlist_gc(&(donors_minus[i]));
	substringlist_gc(&(antidonors_minus[i]));
	substringlist_gc(&(acceptors_minus[i]));
	substringlist_gc(&(antiacceptors_minus[i]));
      }
      FREE(donors_plus);
      FREE(antidonors_plus);
      FREE(acceptors_plus);
      FREE(antiacceptors_plus);
      FREE(donors_minus);
      FREE(antidonors_minus);
      FREE(acceptors_minus);
      FREE(antiacceptors_minus);

    }

    if (alloc_floors_p == true) {
      Floors_free(&floors);
    }
  }
    
  FREE(minus_segments);
  FREE(plus_segments);

  return List_append(subs,List_append(indels,List_append(localsplicing,List_append(halfintrons,distantsplicing))));
}


List_T
Stage1_single_read (Sequence_T queryseq, Indexdb_T indexdb, Indexdb_T indexdb2,
		    int indexdb_size_threshold, IIT_T geneprob_iit, IIT_T chromosome_iit, Genome_T genome,
		    Genome_T genomealt, Floors_T *floors_array,
		    bool knownsplicingp, bool novelsplicingp, bool canonicalp, bool trim_ends_p,
		    int maxpaths, int maxchimerapaths, double usermax_level_float, int subopt_levels, 
		    Masktype_T masktype, int indel_penalty, int max_middle_insertions, int max_middle_deletions,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    Genomicpos_T shortsplicedist,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
		    int min_distantsplicing_end_matches, double min_distantsplicing_identity,
		    Genomicpos_T *splicesites, Splicetype_T *splicetypes,
		    int nsplicesites, bool dibasep, bool cmetp) {
  List_T hits = NULL;
  T this = NULL;
  int usermax_level;
  int querylength, query_lastpos, cutoff_level;
  char *queryuc_ptr;
  bool allvalidp;
  UINT4 *genome_blocks, *snp_blocks;

  if ((querylength = Sequence_fulllength(queryseq)) <= INDEX1PART+2) {
    fprintf(stderr,"GSNAP cannot handle reads shorter than %d bp\n",INDEX1PART+2);
    return (List_T) NULL;

  } else if (querylength > MAX_QUERYLENGTH) {
    fprintf(stderr,"GSNAP cannot handle reads longer than %d bp.  Either recompile with a higher value of MAX_QUERYLENGTH, or consider using GMAP instead.\n",
	    MAX_QUERYLENGTH);
    return (List_T) NULL;

  } else {
    if (usermax_level_float < 0.0) {
      usermax_level = -1;
    } else if (usermax_level_float > 0.0 && usermax_level_float < 1.0) {
      usermax_level = (int) rint(usermax_level_float * (double) querylength);
    } else {
      usermax_level = (int) usermax_level_float;
    }

    this = Stage1_new(querylength);
    queryuc_ptr = Sequence_fullpointer_uc(queryseq);
    query_lastpos = querylength - INDEX1PART;
    if (read_oligos(&allvalidp,this,queryuc_ptr,querylength,query_lastpos,dibasep,cmetp) == 0) {
      debug(printf("Aborting because no hits found anywhere\n"));
      Stage1_free(&this,querylength);
      return (List_T) NULL;
    } else {
      genome_blocks = Genome_blocks(genome);
      snp_blocks = genomealt ? Genome_blocks(genomealt) : NULL;
      hits = align_end(&cutoff_level,this,queryuc_ptr,querylength,query_lastpos,
		       indexdb,indexdb2,indexdb_size_threshold,
		       chromosome_iit,genome_blocks,snp_blocks,
		       genome,genomealt,floors_array,
		       splicesites,splicetypes,nsplicesites,
		       usermax_level,subopt_levels,masktype,indel_penalty,
		       maxpaths,maxchimerapaths,knownsplicingp,novelsplicingp,
		       canonicalp,shortsplicedist,
		       localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
		       min_distantsplicing_end_matches,min_distantsplicing_identity,
		       max_middle_insertions,max_middle_deletions,
		       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
		       allvalidp,trim_ends_p,dibasep,cmetp);
      Stage1_free(&this,querylength);
      hits = Stage3_remove_duplicates(Stage3_optimal_score(hits,cutoff_level,subopt_levels));

      if (List_length(hits) > 1 && geneprob_iit != NULL) {
	Stage3_geneprob_eval(hits,geneprob_iit,chromosome_iit,queryseq);
      }

      return hits;
    }
  }
}


#define HITARRAY_SUBS 0
#define HITARRAY_INDELS 1
#define HITARRAY_LOCALSPLICING 2
#define HITARRAY_HALFINTRONS 3
#define HITARRAY_DISTANTSPLICING 4
#define HITARRAY_N 5

static List_T
align_pair (bool *abort_pairing_p, int *nconcordant, int *nsamechr, int *cutoff_level_5, int *cutoff_level_3,
	    List_T *hits5, List_T *hits3, T this5, T this3,
	    char *queryuc_ptr_5, char *queryuc_ptr_3,
	    int querylength5, int querylength3, int query5_lastpos, int query3_lastpos,
	    Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
	    IIT_T chromosome_iit, UINT4 *genome_blocks, UINT4 *snp_blocks,
	    Genome_T genome, Genome_T genomealt, Floors_T *floors_array,
	    Genomicpos_T *splicesites, Splicetype_T *splicetypes, int nsplicesites,
	    int usermax_level_5, int usermax_level_3, int subopt_levels, Masktype_T masktype, int indel_penalty,
	    int maxpaths, int maxchimerapaths, bool knownsplicingp, bool novelsplicingp, 
	    bool canonicalp, Genomicpos_T shortsplicedist,
	    int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
	    int min_distantsplicing_end_matches, double min_distantsplicing_identity,
	    int max_middle_insertions, int max_middle_deletions,
	    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
	    bool allvalidp5, bool allvalidp3, bool trim_ends_p, bool dibasep, bool cmetp,
	    Genomicpos_T pairmax, Genomicpos_T expected_pairlength, int maxpairedpaths) {

  List_T hitpairs = NULL;
  List_T hitarray5[HITARRAY_N], hitarray3[HITARRAY_N];
  List_T subs5 = NULL, indels5 = NULL, localsplicing5 = NULL, halfintrons5 = NULL, distantsplicing5 = NULL;
  List_T subs3 = NULL, indels3 = NULL, localsplicing3 = NULL, halfintrons3 = NULL, distantsplicing3 = NULL;
  struct Splicesegment_T *plus_segments_5 = NULL, *minus_segments_5 = NULL,
    *plus_segments_3 = NULL, *minus_segments_3 = NULL;
  int maxscore5, maxscore3;
  int plus_nsegments_5, minus_nsegments_5, plus_nsegments_3, minus_nsegments_3;
  char queryrc5[MAX_QUERYLENGTH+1], queryrc3[MAX_QUERYLENGTH+1];
  int found_score, ignore_score, done_level_5, done_level_3, max_level, fast_level_5, fast_level_3,
    mismatch_level_5, mismatch_level_3;
  int max_splice_mismatches_5 = -1, max_splice_mismatches_3 = -1, i;
  int nsplicepairs5 = 0, nsplicepairs3 = 0;
  List_T *donors_plus_5, *antidonors_plus_5, *acceptors_plus_5, *antiacceptors_plus_5,
    *donors_minus_5, *antidonors_minus_5, *acceptors_minus_5, *antiacceptors_minus_5;
  List_T *donors_plus_3, *antidonors_plus_3, *acceptors_plus_3, *antiacceptors_plus_3,
    *donors_minus_3, *antidonors_minus_3, *acceptors_minus_3, *antiacceptors_minus_3;
  bool any_omitted_p_5, any_omitted_p_3;
  Floors_T floors5, floors3;
  bool alloc_floors_p_5 = false, alloc_floors_p_3 = false, alloc5p = false, alloc3p = false;
  int nmismatches;

  char gbuffer5[MAX_QUERYLENGTH+1], gbuffer3[MAX_QUERYLENGTH+1];
  gbuffer5[querylength5] = '\0';
  gbuffer3[querylength3] = '\0';

  make_complement_buffered(queryrc5,queryuc_ptr_5,querylength5);
  make_complement_buffered(queryrc3,queryuc_ptr_3,querylength3);

  *nconcordant = *nsamechr = 0;

  found_score = querylength5 + querylength3;
  ignore_score = querylength5 + querylength3;

  fast_level_5 = (querylength5 + 2)/INDEX1PART - NREQUIRED_FAST;
  fast_level_3 = (querylength3 + 2)/INDEX1PART - NREQUIRED_FAST;
  if (fast_level_5 < 1 && usermax_level_5 < 0) {
    fast_level_5 = 1;		/* Do at least 1 mismatch */
  }
  if (fast_level_3 < 1 && usermax_level_3 < 0) {
    fast_level_3 = 1;		/* Do at least 1 mismatch */
  }

  if (usermax_level_5 < 0) {
    *cutoff_level_5 = fast_level_5;
  } else {
    *cutoff_level_5 = usermax_level_5;
  }

  if (usermax_level_3 < 0) {
    *cutoff_level_3 = fast_level_3;
  } else {
    *cutoff_level_3 = usermax_level_3;
  }

  if (dibasep) {
    max_level = querylength5 + querylength3;
    done_level_5 = querylength5;
    done_level_3 = querylength3;
  } else if (usermax_level_5 < 0) {
    max_level = fast_level_5 + fast_level_3;
    done_level_5 = fast_level_5 + subopt_levels;
    done_level_3 = fast_level_3 + subopt_levels;
  } else {
    max_level = usermax_level_5 + usermax_level_3;
    done_level_5 = usermax_level_5 + subopt_levels;
    done_level_3 = usermax_level_3 + subopt_levels;
  }
  debug(printf("0> max_level %d, done_level %d,%d\n",max_level,done_level_5,done_level_3));

  /* 1A. Exact.  Requires compress if cmet or genomealt.  Creates and uses spanning set. */
  mismatch_level_5 = 0;
  if (allvalidp5 == false) {
    debug(printf("Not all oligos in 5' end are valid, so cannot perform spanning set\n"));
    fast_level_5 = -1;
  } else {
    debug(printf("fast_level_5 = %d\n",fast_level_5));
    debug(printf("*** Stage 1.  Exact ***\n"));
    if (cmetp || genomealt) {
      this5->query_compress_fwd = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/true);
      this5->query_compress_rev = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/false);
    }
    subs5 = find_spanning_exact_matches(&ignore_score,this5,queryuc_ptr_5,queryrc5,querylength5,query5_lastpos,
					indexdb,indexdb2,gbuffer5,this5->query_compress_fwd,this5->query_compress_rev,
					genome_blocks,snp_blocks,genome,genomealt,chromosome_iit,trim_ends_p,dibasep,cmetp);
    mismatch_level_5 = 1;
  }

  /* 1B. Exact.  Requires compress if cmet or genomealt.  Creates and uses spanning set. */
  mismatch_level_3 = 0;
  if (allvalidp3 == false) {
    debug(printf("Not all oligos in 3' end are valid, so cannot perform spanning set\n"));
    fast_level_3 = -1;
  } else {
    debug(printf("fast_level_3 = %d\n",fast_level_3));
    debug(printf("*** Stage 1.  Exact ***\n"));
    if (cmetp || genomealt) {
      this3->query_compress_fwd = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/true);
      this3->query_compress_rev = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/false);
    }
    subs3 = find_spanning_exact_matches(&ignore_score,this3,queryuc_ptr_3,queryrc3,querylength3,query3_lastpos,
					indexdb,indexdb2,gbuffer3,this3->query_compress_fwd,this3->query_compress_rev,
					genome_blocks,snp_blocks,genome,genomealt,chromosome_iit,trim_ends_p,dibasep,cmetp);
    mismatch_level_3 = 1;
  }
    
  /* 1. Pairing */
  hitarray5[HITARRAY_SUBS] = subs5 = Stage3_remove_duplicates(subs5);
  hitarray3[HITARRAY_SUBS] = subs3 = Stage3_remove_duplicates(subs3);
  hitpairs = Stage3_pair_up(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*nsamechr),hitpairs,
			    hitarray5,/*narray5*/HITARRAY_SUBS+1,
			    hitarray3,/*narray3*/HITARRAY_SUBS+1,
			    *cutoff_level_5,*cutoff_level_3,subopt_levels,
			    pairmax,expected_pairlength,querylength5,querylength3,maxpairedpaths);
  debug(printf("After pairing exact, found %d concordant, found_score %d\n",*nconcordant,found_score));
  if (*abort_pairing_p == true) {
    *hits5 = subs5;
    *hits3 = subs3;
    return hitpairs;
  } else {
    max_level = (found_score < max_level) ? found_score : max_level;
    if (max_level + subopt_levels < done_level_5) {
      done_level_5 = max_level + subopt_levels;
    }
    if (max_level + subopt_levels < done_level_3) {
      done_level_3 = max_level + subopt_levels;
    }
    debug(printf("1> found_score = %d, max_level %d, done_level %d,%d\n",found_score,max_level,done_level_5,done_level_3));
  }


  if (dibasep == false) {
    /* 2A. One mismatch.  Requires spanning set and compress. */
    if (allvalidp5 && querylength5 >= 22 && done_level_5 >= 1) {
      debug(printf("*** Stage 2A.  One miss ***\n"));
      if (this5->query_compress_fwd == NULL) {
	this5->query_compress_fwd = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/true);
	this5->query_compress_rev = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/false);
      }
      subs5 = find_spanning_onemiss_matches(&ignore_score,subs5,this5,
					    queryuc_ptr_5,queryrc5,querylength5,
					    gbuffer5,this5->query_compress_fwd,this5->query_compress_rev,
					    genome_blocks,snp_blocks,genome,genomealt,chromosome_iit,trim_ends_p,dibasep,cmetp);
      mismatch_level_5 = 2;
    }

    /* 2B. One mismatch.  Requires spanning set and compress. */
    if (allvalidp3 && querylength3 >= 22 && done_level_3 >= 1) {
      debug(printf("*** Stage 2B.  One miss ***\n"));
      if (this3->query_compress_fwd == NULL) {
	this3->query_compress_fwd = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/true);
	this3->query_compress_rev = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/false);
      }
      subs3 = find_spanning_onemiss_matches(&ignore_score,subs3,this3,
					    queryuc_ptr_3,queryrc3,querylength3,
					    gbuffer3,this3->query_compress_fwd,this3->query_compress_rev,
					    genome_blocks,snp_blocks,genome,genomealt,chromosome_iit,trim_ends_p,dibasep,cmetp);
      mismatch_level_3 = 2;
    }
  
    /* 2. Pairing */
    hitarray5[HITARRAY_SUBS] = subs5 = Stage3_remove_duplicates(subs5);
    hitarray3[HITARRAY_SUBS] = subs3 = Stage3_remove_duplicates(subs3);
    hitpairs = Stage3_pair_up(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*nsamechr),hitpairs,
			      hitarray5,/*narray5*/HITARRAY_SUBS+1,
			      hitarray3,/*narray3*/HITARRAY_SUBS+1,
			      *cutoff_level_5,*cutoff_level_3,subopt_levels,
			      pairmax,expected_pairlength,querylength5,querylength3,maxpairedpaths);
    debug(printf("After pairing one mismatch, found %d concordant, found_score %d\n",*nconcordant,found_score));
    if (*abort_pairing_p == true) {
      *hits5 = subs5;
      *hits3 = subs3;
      return hitpairs;
    } else {
      max_level = (found_score < max_level) ? found_score : max_level;
      if (max_level + subopt_levels < done_level_5) {
	done_level_5 = max_level + subopt_levels;
      }
      if (max_level + subopt_levels < done_level_3) {
	done_level_3 = max_level + subopt_levels;
      }
      debug(printf("2> found_score = %d, max_level %d, done_level %d,%d\n",found_score,max_level,done_level_5,done_level_3));
    }
    

    /* 3A. Mismatches via spanning set.  Requires spanning set and compress. */
    if (allvalidp5 && done_level_5 >= 2) {
      if (this5->query_compress_fwd == NULL) {
	this5->query_compress_fwd = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/true);
	this5->query_compress_rev = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/false);
      }
      /* NOTE: Since done_level isn't updated, can do in one batch instead of iteratively */
      while (mismatch_level_5 <= fast_level_5 && mismatch_level_5 <= done_level_5) {
	debug(printf("*** Stage 3A (level %d).  Spanning set mismatches ***\n",mismatch_level_5));
	subs5 = find_spanning_multimiss_matches(&ignore_score,subs5,this5,NREQUIRED_FAST,
						queryuc_ptr_5,queryrc5,querylength5,
						gbuffer5,this5->query_compress_fwd,this5->query_compress_rev,
						genome_blocks,snp_blocks,genome,genomealt,chromosome_iit,
						/*nmisses_allowed*/mismatch_level_5,trim_ends_p,dibasep,cmetp);
	mismatch_level_5++;
      }
    }

    /* 3B. Mismatches via spanning set.  Requires spanning set and compress. */
    if (allvalidp3 && done_level_3 >= 2) {
      if (this3->query_compress_fwd == NULL) {
	this3->query_compress_fwd = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/true);
	this3->query_compress_rev = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/false);
      }
      /* NOTE: Since done_level isn't updated, can do in one batch instead of iteratively */
      while (mismatch_level_3 <= fast_level_3 && mismatch_level_3 <= done_level_3) {
	debug(printf("*** Stage 3B (level %d).  Spanning set mismatches ***\n",mismatch_level_3));
	subs3 = find_spanning_multimiss_matches(&ignore_score,subs3,this3,NREQUIRED_FAST,
						queryuc_ptr_3,queryrc3,querylength3,
						gbuffer3,this3->query_compress_fwd,this3->query_compress_rev,
						genome_blocks,snp_blocks,genome,genomealt,chromosome_iit,
						/*nmisses_allowed*/mismatch_level_3,trim_ends_p,dibasep,cmetp);
	mismatch_level_3++;
      }
    }
    
    /* 3. Pairing */
    hitarray5[HITARRAY_SUBS] = subs5 = Stage3_remove_duplicates(subs5);
    hitarray3[HITARRAY_SUBS] = subs3 = Stage3_remove_duplicates(subs3);
    hitpairs = Stage3_pair_up(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*nsamechr),hitpairs,
			      hitarray5,/*narray5*/HITARRAY_SUBS+1,
			      hitarray3,/*narray3*/HITARRAY_SUBS+1,
			      *cutoff_level_5,*cutoff_level_3,subopt_levels,
			      pairmax,expected_pairlength,querylength5,querylength3,maxpairedpaths);
    debug(printf("After pairing spanning set, found %d concordant, found_score %d\n",*nconcordant,found_score));
    if (*abort_pairing_p == true) {
      *hits5 = subs5;
      *hits3 = subs3;
      return hitpairs;
    } else {
      max_level = (found_score < max_level) ? found_score : max_level;
      if (max_level + subopt_levels < done_level_5) {
	done_level_5 = max_level + subopt_levels;
      }
      if (max_level + subopt_levels < done_level_3) {
	done_level_3 = max_level + subopt_levels;
      }
      debug(printf("3> found_score = %d, max_level %d, done_level %d,%d\n",found_score,max_level,done_level_5,done_level_3));
    }
  }


  /* 4/5A.  Complete set mismatches and indels, omitting frequent oligos */
#if 0
  if (subs5 != NULL) {
    indels_prior_penalty_5 = 1;
  }
#endif
  if (dibasep || done_level_5 > fast_level_5 || done_level_5 >= indel_penalty /* + indels_prior_penalty_5*/) {
    if (masktype == MASK_NONE) {
      debug(printf("*** Stage 4A,5A.  Complete mm/indels with no masking with done_level %d ***\n",done_level_5));
      complete_set_mm_indels(&ignore_score,&any_omitted_p_5,&max_level,&done_level_5,/*revise_levels_p*/false,
			     &subs5,&indels5,this5,queryuc_ptr_5,queryrc5,
			     querylength5,query5_lastpos,gbuffer5,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,genome,floors_array,
			     subopt_levels,indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     trim_ends_p,dibasep,cmetp,fast_level_5,
			     /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
    } else {
      debug(printf("*** Stage 4A,5A.  Complete mm/indels masking frequent oligos with done_level %d ***\n",done_level_5));
      complete_set_mm_indels(&ignore_score,&any_omitted_p_5,&max_level,&done_level_5,/*revise_levels_p*/false,
			     &subs5,&indels5,this5,queryuc_ptr_5,queryrc5,
			     querylength5,query5_lastpos,gbuffer5,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,genome,floors_array,
			     subopt_levels,indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     trim_ends_p,dibasep,cmetp,fast_level_5,/*omit_frequent_p*/true,
			     /*omit_repetitive_p*/(masktype == MASK_REPETITIVE || masktype == MASK_GREEDY_REPETITIVE) ? true : false
			     );
      if ((masktype == MASK_GREEDY_FREQUENT || masktype == MASK_GREEDY_REPETITIVE) && subs5 == NULL && indels5 == NULL && any_omitted_p_5 == true) {
	/* 4/5A.  Complete set mismatches and indels, with all oligos */
	debug(printf("*** Stage 4A,5A.  Complete mm/indels with no masking with done_level %d ***\n",done_level_5));
	complete_set_mm_indels(&ignore_score,&any_omitted_p_5,&max_level,&done_level_5,/*revise_levels_p*/false,
			       &subs5,&indels5,this5,queryuc_ptr_5,queryrc5,
			       querylength5,query5_lastpos,gbuffer5,indexdb,indexdb2,indexdb_size_threshold,
			       chromosome_iit,genome_blocks,snp_blocks,genome,floors_array,
			       subopt_levels,indel_penalty,max_middle_insertions,max_middle_deletions,
			       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			       trim_ends_p,dibasep,cmetp,fast_level_5,
			       /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
      }
    }
  }

  /* 4/5B.  Complete set mismatches and indels, omitting frequent oligos */
#if 0
  if (subs3 != NULL) {
    indels_prior_penalty_3 = 1;
  }
#endif
  if (dibasep || done_level_3 > fast_level_3 || done_level_3 >= indel_penalty /* + indels_prior_penalty_3*/) {
    if (masktype == MASK_NONE) {
      debug(printf("*** Stage 4B,5B.  Complete mm/indels with no masking with done_level %d ***\n",done_level_3));
      complete_set_mm_indels(&ignore_score,&any_omitted_p_3,&max_level,&done_level_3,/*revise_levels_p*/false,
			     &subs3,&indels3,this3,queryuc_ptr_3,queryrc3,
			     querylength3,query3_lastpos,gbuffer3,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,genome,floors_array,
			     subopt_levels,indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     trim_ends_p,dibasep,cmetp,fast_level_3,
			     /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
    } else {
      debug(printf("*** Stage 4B,5B.  Complete mm/indels masking frequent oligos with done_level %d ***\n",done_level_3));
      complete_set_mm_indels(&ignore_score,&any_omitted_p_3,&max_level,&done_level_3,/*revise_levels_p*/false,
			     &subs3,&indels3,this3,queryuc_ptr_3,queryrc3,
			     querylength3,query3_lastpos,gbuffer3,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,genome,floors_array,
			     subopt_levels,indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     trim_ends_p,dibasep,cmetp,fast_level_3,/*omit_frequent_p*/true,
			     /*omit_repetitive_p*/(masktype == MASK_REPETITIVE || masktype == MASK_GREEDY_REPETITIVE) ? true : false
			     );
      if ((masktype == MASK_GREEDY_FREQUENT || masktype == MASK_GREEDY_REPETITIVE) && subs3 == NULL && indels3 == NULL && any_omitted_p_3 == true) {
	/* 4/5B.  Complete set mismatches and indels, with all oligos */
	debug(printf("*** Stage 4B,5B.  Complete mm/indels with no masking with done_level %d ***\n",done_level_3));
	complete_set_mm_indels(&ignore_score,&any_omitted_p_3,&max_level,&done_level_3,/*revise_levels_p*/false,
			       &subs3,&indels3,this3,queryuc_ptr_3,queryrc3,
			       querylength3,query3_lastpos,gbuffer3,indexdb,indexdb2,indexdb_size_threshold,
			       chromosome_iit,genome_blocks,snp_blocks,genome,floors_array,
			       subopt_levels,indel_penalty,max_middle_insertions,max_middle_deletions,
			       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			       trim_ends_p,dibasep,cmetp,fast_level_3,
			       /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
      }
    }
  }

  /* 4/5. Pairing */
  debug(printf("Starting pairing of 4 and 5\n"));
  hitarray5[HITARRAY_SUBS] = subs5 = Stage3_remove_duplicates(subs5);
  hitarray5[HITARRAY_INDELS] = indels5 = Stage3_remove_duplicates(indels5);
  hitarray3[HITARRAY_SUBS] = subs3 = Stage3_remove_duplicates(subs3);
  hitarray3[HITARRAY_INDELS] = indels3 = Stage3_remove_duplicates(indels3);
  hitpairs = Stage3_pair_up(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*nsamechr),hitpairs,
			    hitarray5,/*narray5*/HITARRAY_INDELS+1,
			    hitarray3,/*narray3*/HITARRAY_INDELS+1,
			    *cutoff_level_5,*cutoff_level_3,subopt_levels,
			    pairmax,expected_pairlength,querylength5,querylength3,maxpairedpaths);
  debug(printf("After pairing complete set mismatches and indels, found %d concordant, found_score %d\n",*nconcordant,found_score));
  if (*abort_pairing_p == true) {
    *hits5 = List_append(subs5,indels5);
    *hits3 = List_append(subs3,indels3);
    return hitpairs;
  } else {
    max_level = (found_score < max_level) ? found_score : max_level;
    if (max_level + subopt_levels < done_level_5) {
      done_level_5 = max_level + subopt_levels;
    }
    if (max_level + subopt_levels < done_level_3) {
      done_level_3 = max_level + subopt_levels;
    }
    debug(printf("4/5> found_score = %d, max_level %d, done_level %d,%d\n",found_score,max_level,done_level_5,done_level_3));
  }


  /* 6/7/8A. Known splicing.  Requires compress and all positions fetched. */
  /* Subtract 1 from done_level for previous hits */
  if (knownsplicingp || novelsplicingp) {
#if 0
    if (subs5 || indels5) {
      localsplicing_prior_penalty_5 = 1;
    }
#endif
    debug(printf("Deciding whether to do localsplicing: done_level_5 %d >=? localsplicing_penalty %d\n",
		 done_level_5,localsplicing_penalty));

    if (*nconcordant == 0 && done_level_5 >= localsplicing_penalty /* + localsplicing_prior_penalty_5*/) {
      debug(printf("*** Stage 6A.  Local splicing masking frequent oligos with done_level %d ***\n",done_level_5));
      /* Always mask frequent oligos for splicing, which must be transcriptional */
      localsplicing5 = complete_set_localsplicing(&ignore_score,&any_omitted_p_5,&floors5,&alloc_floors_p_5,
						  &plus_segments_5,&plus_nsegments_5,&minus_segments_5,&minus_nsegments_5,
						  localsplicing5,this5,queryuc_ptr_5,
						  querylength5,query5_lastpos,indexdb,indexdb2,indexdb_size_threshold,
						  chromosome_iit,genome_blocks,snp_blocks,genome,
						  floors_array,max_end_insertions,
						  splicesites,splicetypes,nsplicesites,shortsplicedist,novelsplicingp,
						  localsplicing_penalty,min_localsplicing_end_matches,
						  /*max_mismatches_allowed*/done_level_5 - localsplicing_penalty,
						  trim_ends_p,dibasep,cmetp,/*omit_frequent_p*/true,/*omit_repetitive_p*/true);
    }


#if 0
    if (subs3 || indels3) {
      localsplicing_prior_penalty_3 = 1;
    }
#endif
    debug(printf("Deciding whether to do localsplicing: done_level_3 %d >=? localsplicing_penalty %d\n",
		 done_level_3,localsplicing_penalty));
    if (*nconcordant == 0 && done_level_3 >= localsplicing_penalty /* + localsplicing_prior_penalty_3 */) {
      debug(printf("*** Stage 6B.  Local splicing masking frequent oligos with done_level %d ***\n",done_level_3));
      /* Always mask frequent oligos for splicing, which must be transcriptional */
      localsplicing3 = complete_set_localsplicing(&ignore_score,&any_omitted_p_3,&floors3,&alloc_floors_p_3,
						  &plus_segments_3,&plus_nsegments_3,&minus_segments_3,&minus_nsegments_3,
						  localsplicing3,this3,queryuc_ptr_3,
						  querylength3,query3_lastpos,indexdb,indexdb2,indexdb_size_threshold,
						  chromosome_iit,genome_blocks,snp_blocks,genome,
						  floors_array,max_end_insertions,
						  splicesites,splicetypes,nsplicesites,shortsplicedist,novelsplicingp,
						  localsplicing_penalty,min_localsplicing_end_matches,
						  /*max_mismatches_allowed*/done_level_3 - localsplicing_penalty,
						  trim_ends_p,dibasep,cmetp,/*omit_frequent_p*/true,/*omit_repetitive_p*/true);
    }
    
    /* 6.  Pairing */
    hitarray5[HITARRAY_LOCALSPLICING] = localsplicing5; /* not removing duplicates */
    hitarray3[HITARRAY_LOCALSPLICING] = localsplicing3; /* not removing duplicates */
    if (*nconcordant == 0) {
      hitpairs = Stage3_pair_up(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*nsamechr),hitpairs,
				hitarray5,/*narray5*/HITARRAY_LOCALSPLICING+1,
				hitarray3,/*narray3*/HITARRAY_LOCALSPLICING+1,
				*cutoff_level_5,*cutoff_level_3,subopt_levels,
				pairmax,expected_pairlength,querylength5,querylength3,maxpairedpaths);
      debug(printf("After pairing local splicing, found %d concordant, found_score %d\n",*nconcordant,found_score));
      if (*abort_pairing_p == true) {
	*hits5 = List_append(subs5,List_append(indels5,localsplicing5));
	*hits3 = List_append(subs3,List_append(indels3,localsplicing3));
	return hitpairs;
      } else {
	max_level = (found_score < max_level) ? found_score : max_level;
	if (max_level + subopt_levels < done_level_5) {
	  done_level_5 = max_level + subopt_levels;
	}
	if (max_level + subopt_levels < done_level_3) {
	  done_level_3 = max_level + subopt_levels;
	}
	debug(printf("Pairing after 6A and 6B> found_score = %d, max_level %d, done_level %d,%d\n",
		     found_score,max_level,done_level_5,done_level_3));
      }
    }

    /* 7,8A.  Distant and half splicing */
#if 0
    if (subs5 || indels5 || localsplicing5) {
      distantsplicing_prior_penalty_5 = 1;
    }
#endif

    max_splice_mismatches_5 = done_level_5 - distantsplicing_penalty /* - distantsplicing_prior_penalty_5*/;
    if (localsplicing5 == NULL /* && (knownsplicingp == true || querylength5 >= novel_half_intron_min_support) */) {
      if (done_level_5 - localsplicing_penalty /* - localsplicing_prior_penalty_5*/ > max_splice_mismatches_5) {
	max_splice_mismatches_5 = done_level_5 - localsplicing_penalty /* - localsplicing_prior_penalty_5*/;
      }
    }

    if (*nconcordant == 0 && max_splice_mismatches_5 >= 0) {
      alloc5p = true;
      donors_plus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      antidonors_plus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      acceptors_plus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      antiacceptors_plus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      donors_minus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      antidonors_minus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      acceptors_minus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      antiacceptors_minus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));

      if (plus_segments_5) {
	identify_spliceends_plus(&donors_plus_5,&antidonors_plus_5,&acceptors_plus_5,&antiacceptors_plus_5,
				 plus_segments_5,plus_nsegments_5,genome_blocks,snp_blocks,genome,
				 queryuc_ptr_5,floors5,querylength5,query5_lastpos,gbuffer5,
				 /*query_compress*/this5->query_compress_fwd,
				 splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
				 /*max_mismatches_allowed*/max_splice_mismatches_5,trim_ends_p,dibasep,cmetp);
      }
      if (minus_segments_5) {
	identify_spliceends_minus(&donors_minus_5,&antidonors_minus_5,&acceptors_minus_5,&antiacceptors_minus_5,
				  minus_segments_5,minus_nsegments_5,genome_blocks,snp_blocks,genome,
				  queryuc_ptr_5,queryrc5,floors5,querylength5,query5_lastpos,gbuffer5,
				  /*query_compress*/this5->query_compress_rev,
				  splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
				  /*max_mismatches_allowed*/max_splice_mismatches_5,trim_ends_p,dibasep,cmetp);
      }

      /* 7A.  Distant splicing */
      nmismatches = 0;
      while (nmismatches <= done_level_5 - distantsplicing_penalty /* - distantsplicing_prior_penalty_5*/ && 
	     nsplicepairs5 < maxchimerapaths) {
	debug(printf("*** Stage 7A.  Distant splicing, allowing %d mismatches ***\n",nmismatches));

	debug4(printf("Sorting splice ends\n"));
	donors_plus_5[nmismatches] = Substring_sort_chimera_halves(donors_plus_5[nmismatches],/*ascendingp*/true);
	acceptors_plus_5[nmismatches] = Substring_sort_chimera_halves(acceptors_plus_5[nmismatches],/*ascendingp*/true);

	antidonors_plus_5[nmismatches] = Substring_sort_chimera_halves(antidonors_plus_5[nmismatches],/*ascendingp*/false);
	antiacceptors_plus_5[nmismatches] = Substring_sort_chimera_halves(antiacceptors_plus_5[nmismatches],/*ascendingp*/false);
	  
	donors_minus_5[nmismatches] = Substring_sort_chimera_halves(donors_minus_5[nmismatches],/*ascendingp*/false);
	acceptors_minus_5[nmismatches] = Substring_sort_chimera_halves(acceptors_minus_5[nmismatches],/*ascendingp*/false);

	antidonors_minus_5[nmismatches] = Substring_sort_chimera_halves(antidonors_minus_5[nmismatches],/*ascendingp*/true);
	antiacceptors_minus_5[nmismatches] = Substring_sort_chimera_halves(antiacceptors_minus_5[nmismatches],/*ascendingp*/true);

	debug4(printf("Splice ends at %d nmismatches: +donors/acceptors %d/%d, +antidonors/antiacceptors %d/%d, -donors/acceptors %d/%d, -antidonors/antiacceptors %d/%d\n",
		      nmismatches,
		      List_length(donors_plus_5[nmismatches]),List_length(acceptors_plus_5[nmismatches]),
		      List_length(antidonors_plus_5[nmismatches]),List_length(antiacceptors_plus_5[nmismatches]),
		      List_length(donors_minus_5[nmismatches]),List_length(acceptors_minus_5[nmismatches]),
		      List_length(antidonors_minus_5[nmismatches]),List_length(antiacceptors_minus_5[nmismatches])));

	distantsplicing5 = find_splicepairs_distant(&ignore_score,&nsplicepairs5,distantsplicing5,
						    donors_plus_5,antidonors_plus_5,acceptors_plus_5,antiacceptors_plus_5,
						    donors_minus_5,antidonors_minus_5,acceptors_minus_5,antiacceptors_minus_5,
						    shortsplicedist,localsplicing_penalty,distantsplicing_penalty,
						    min_distantsplicing_end_matches,min_distantsplicing_identity,
						    querylength5,nmismatches,maxchimerapaths);
	nmismatches++;
      }

      /* 8A.  Half introns */
      debug(printf("%d 5' local splices, %d 5' distant splices\n",
		   List_length(localsplicing5),List_length(distantsplicing5)));
      if (localsplicing5 == NULL && distantsplicing5 == NULL) {
	debug(printf("*** Stage 8A.  Half introns up to %d mismatches ***\n",max_splice_mismatches_5));
	halfintrons5 = find_splices_halfintrons(&ignore_score,donors_plus_5,antidonors_plus_5,acceptors_plus_5,antiacceptors_plus_5,
						donors_minus_5,antidonors_minus_5,acceptors_minus_5,antiacceptors_minus_5,
						/*max_mismatches*/max_splice_mismatches_5,/*splicing_penalty*/localsplicing_penalty,
						querylength5);
      }
    }


    /* 7,8B.  Distant and half splicing */
#if 0
    if (subs3 || indels3 || localsplicing3) {
      distantsplicing_prior_penalty_3 = 1;
    }
#endif

    max_splice_mismatches_3 = done_level_3 - distantsplicing_penalty /* - distantsplicing_prior_penalty_3*/;
    if (localsplicing3 == NULL /* && (knownsplicingp == true || querylength3 >= novel_half_intron_min_support) */) {
      if (done_level_3 - localsplicing_penalty /* - localsplicing_prior_penalty_3*/ > max_splice_mismatches_3) {
	max_splice_mismatches_3 = done_level_3 - localsplicing_penalty /* - localsplicing_prior_penalty_3*/;
      }
    }

    if (*nconcordant == 0 && max_splice_mismatches_3 >= 0) {
      alloc3p = true;
      donors_plus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      antidonors_plus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      acceptors_plus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      antiacceptors_plus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      donors_minus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      antidonors_minus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      acceptors_minus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      antiacceptors_minus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));

      if (plus_segments_3) {
	identify_spliceends_plus(&donors_plus_3,&antidonors_plus_3,&acceptors_plus_3,&antiacceptors_plus_3,
				 plus_segments_3,plus_nsegments_3,genome_blocks,snp_blocks,genome,
				 queryuc_ptr_3,floors3,querylength3,query3_lastpos,gbuffer3,
				 /*query_compress*/this3->query_compress_fwd,
				 splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
				 /*max_mismatches_allowed*/max_splice_mismatches_3,trim_ends_p,dibasep,cmetp);
      }
      if (minus_segments_3) {
	identify_spliceends_minus(&donors_minus_3,&antidonors_minus_3,&acceptors_minus_3,&antiacceptors_minus_3,
				  minus_segments_3,minus_nsegments_3,genome_blocks,snp_blocks,genome,
				  queryuc_ptr_3,queryrc3,floors3,querylength3,query3_lastpos,gbuffer3,
				  /*query_compress*/this3->query_compress_rev,
				  splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
				  /*max_mismatches_allowed*/max_splice_mismatches_3,trim_ends_p,dibasep,cmetp);
      }

      /* 7B.  Distant splicing */
      nmismatches = 0;
      while (nmismatches <= done_level_3 - distantsplicing_penalty /* - distantsplicing_prior_penalty_3*/ &&
	     nsplicepairs3 < maxchimerapaths) {
	debug(printf("*** Stage 7B.  Distant splicing, allowing %d mismatches ***\n",nmismatches));
	
	debug4(printf("Sorting splice ends\n"));
	donors_plus_3[nmismatches] = Substring_sort_chimera_halves(donors_plus_3[nmismatches],/*ascendingp*/true);
	acceptors_plus_3[nmismatches] = Substring_sort_chimera_halves(acceptors_plus_3[nmismatches],/*ascendingp*/true);
	
	antidonors_plus_3[nmismatches] = Substring_sort_chimera_halves(antidonors_plus_3[nmismatches],/*ascendingp*/false);
	antiacceptors_plus_3[nmismatches] = Substring_sort_chimera_halves(antiacceptors_plus_3[nmismatches],/*ascendingp*/false);
      
	donors_minus_3[nmismatches] = Substring_sort_chimera_halves(donors_minus_3[nmismatches],/*ascendingp*/false);
	acceptors_minus_3[nmismatches] = Substring_sort_chimera_halves(acceptors_minus_3[nmismatches],/*ascendingp*/false);
      
	antidonors_minus_3[nmismatches] = Substring_sort_chimera_halves(antidonors_minus_3[nmismatches],/*ascendingp*/true);
	antiacceptors_minus_3[nmismatches] = Substring_sort_chimera_halves(antiacceptors_minus_3[nmismatches],/*ascendingp*/true);

	debug4(printf("Splice ends at %d nmismatches: +donors/acceptors %d/%d, +antidonors/antiacceptors %d/%d, -donors/acceptors %d/%d, -antidonors/antiacceptors %d/%d\n",
		       nmismatches,
		       List_length(donors_plus_3[nmismatches]),List_length(acceptors_plus_3[nmismatches]),
		       List_length(antidonors_plus_3[nmismatches]),List_length(antiacceptors_plus_3[nmismatches]),
		       List_length(donors_minus_3[nmismatches]),List_length(acceptors_minus_3[nmismatches]),
		       List_length(antidonors_minus_3[nmismatches]),List_length(antiacceptors_minus_3[nmismatches])));
      
	distantsplicing3 = find_splicepairs_distant(&ignore_score,&nsplicepairs3,distantsplicing3,
						    donors_plus_3,antidonors_plus_3,acceptors_plus_3,antiacceptors_plus_3,
						    donors_minus_3,antidonors_minus_3,acceptors_minus_3,antiacceptors_minus_3,
						    shortsplicedist,localsplicing_penalty,distantsplicing_penalty,
						    min_distantsplicing_end_matches,min_distantsplicing_identity,
						    querylength3,nmismatches,maxchimerapaths);
#if 0
	/* May not be applicable to paired-end */
	if (distantsplicing3) {
	  max_level = (found_score < max_level) ? found_score : max_level;
	  done_level = max_level + subopt_levels;
	  debug(printf("7> found_score = %d, max_level %d, done_level %d\n",found_score,max_level,done_level));
	}
#endif

	nmismatches++;
      }

      debug(printf("%d 3' local splices, %d 3' distant splices\n",
		   List_length(localsplicing3),List_length(distantsplicing3)));
      /* 8B.  Half introns */
      if (localsplicing3 == NULL && distantsplicing3 == NULL) {
	debug(printf("*** Stage 8B.  Half introns up to %d mismatches ***\n",max_splice_mismatches_3));
	halfintrons3 = find_splices_halfintrons(&ignore_score,donors_plus_3,antidonors_plus_3,acceptors_plus_3,antiacceptors_plus_3,
						donors_minus_3,antidonors_minus_3,acceptors_minus_3,antiacceptors_minus_3,
						/*max_mismatches*/max_splice_mismatches_3,/*splicing_penalty*/localsplicing_penalty,
						querylength3);
      }
    }

    /* 7,8.  Pairing after distant and half splicing.  Try half splicing first. */
    hitarray5[HITARRAY_LOCALSPLICING] = localsplicing5; /* not removing duplicates */
    hitarray3[HITARRAY_LOCALSPLICING] = localsplicing3; /* not removing duplicates */
    hitarray5[HITARRAY_HALFINTRONS] = halfintrons5;	/* Duplicates already removed */
    hitarray3[HITARRAY_HALFINTRONS] = halfintrons3;	/* Duplicates already removed */
    if (*nconcordant == 0) {
      hitpairs = Stage3_pair_up(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*nsamechr),hitpairs,
				hitarray5,/*narray5*/HITARRAY_HALFINTRONS+1,
				hitarray3,/*narray3*/HITARRAY_HALFINTRONS+1,
				*cutoff_level_5,*cutoff_level_3,subopt_levels,
				pairmax,expected_pairlength,querylength5,querylength3,maxpairedpaths);
      debug(printf("After pairing half splicing, found %d concordant, found_score %d\n",*nconcordant,found_score));
      /* Don't return after abort_pairing, because we need to clean up memory */
      max_level = (found_score < max_level) ? found_score : max_level;
      if (max_level + subopt_levels < done_level_5) {
	done_level_5 = max_level + subopt_levels;
      }
      if (max_level + subopt_levels < done_level_3) {
	done_level_3 = max_level + subopt_levels;
      }
      debug(printf("8> found_score = %d, max_level %d, done_level %d,%d\n",found_score,max_level,done_level_5,done_level_3));
    }

    /* Clean up 5 */
    if (alloc5p == true) {
      for (i = 0; i <= max_splice_mismatches_5; i++) {
	substringlist_gc(&(donors_plus_5[i]));
	substringlist_gc(&(antidonors_plus_5[i]));
	substringlist_gc(&(acceptors_plus_5[i]));
	substringlist_gc(&(antiacceptors_plus_5[i]));
	substringlist_gc(&(donors_minus_5[i]));
	substringlist_gc(&(antidonors_minus_5[i]));
	substringlist_gc(&(acceptors_minus_5[i]));
	substringlist_gc(&(antiacceptors_minus_5[i]));
      }
      FREE(donors_plus_5);
      FREE(antidonors_plus_5);
      FREE(acceptors_plus_5);
      FREE(antiacceptors_plus_5);
      FREE(donors_minus_5);
      FREE(antidonors_minus_5);
      FREE(acceptors_minus_5);
      FREE(antiacceptors_minus_5);
    }

    /* Clean up 3 */
    if (alloc3p == true) {
      for (i = 0; i <= max_splice_mismatches_3; i++) {
	substringlist_gc(&(donors_plus_3[i]));
	substringlist_gc(&(antidonors_plus_3[i]));
	substringlist_gc(&(acceptors_plus_3[i]));
	substringlist_gc(&(antiacceptors_plus_3[i]));
	substringlist_gc(&(donors_minus_3[i]));
	substringlist_gc(&(antidonors_minus_3[i]));
	substringlist_gc(&(acceptors_minus_3[i]));
	substringlist_gc(&(antiacceptors_minus_3[i]));
      }
      FREE(donors_plus_3);
      FREE(antidonors_plus_3);
      FREE(acceptors_plus_3);
      FREE(antiacceptors_plus_3);
      FREE(donors_minus_3);
      FREE(antidonors_minus_3);
      FREE(acceptors_minus_3);
      FREE(antiacceptors_minus_3);
    }

    if (alloc_floors_p_5 == true) {
      Floors_free(&floors5);
    }
    if (alloc_floors_p_3 == true) {
      Floors_free(&floors3);
    }
  }

  FREE(minus_segments_5);
  FREE(plus_segments_5);

  FREE(minus_segments_3);
  FREE(plus_segments_3);

  /* Distant splicing */
  hitarray5[HITARRAY_LOCALSPLICING] = localsplicing5; /* not removing duplicates */
  hitarray5[HITARRAY_HALFINTRONS] = halfintrons5;     /* Duplicates already removed */
  hitarray5[HITARRAY_DISTANTSPLICING] = distantsplicing5; /* Don't want to remove duplicates */
  hitarray3[HITARRAY_LOCALSPLICING] = localsplicing3; /* not removing duplicates */
  hitarray3[HITARRAY_HALFINTRONS] = halfintrons3;     /* Duplicates already removed */
  hitarray3[HITARRAY_DISTANTSPLICING] = distantsplicing3; /* Don't want to remove duplicates */
  if (*nconcordant == 0) {
    hitpairs = Stage3_pair_up(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*nsamechr),hitpairs,
			      hitarray5,/*narray5*/HITARRAY_DISTANTSPLICING+1,
			      hitarray3,/*narray3*/HITARRAY_DISTANTSPLICING+1,
			      *cutoff_level_5,*cutoff_level_3,subopt_levels,
			      pairmax,expected_pairlength,querylength5,querylength3,maxpairedpaths);
    debug(printf("After pairing distant splicing, found %d concordant, found_score %d\n",*nconcordant,found_score));
  }

  *hits5 = List_append(subs5,List_append(indels5,List_append(localsplicing5,List_append(halfintrons5,distantsplicing5))));
  *hits3 = List_append(subs3,List_append(indels3,List_append(localsplicing3,List_append(halfintrons3,distantsplicing3))));

  return hitpairs;
}


List_T
Stage1_paired_read (List_T *singlehits5, List_T *singlehits3,
		    Sequence_T queryseq5, Sequence_T queryseq3,
		    Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
		    IIT_T geneprob_iit, IIT_T chromosome_iit, Genome_T genome, Genome_T genomealt, Floors_T *floors_array,
		    bool knownsplicingp, bool novelsplicingp, bool canonicalp, bool trim_ends_p,
		    int maxpaths, int maxchimerapaths, double usermax_level_float, int subopt_levels, 
		    Masktype_T masktype, int indel_penalty, int max_middle_insertions, int max_middle_deletions,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    Genomicpos_T shortsplicedist,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
		    int min_distantsplicing_end_matches, double min_distantsplicing_identity,
		    Genomicpos_T *splicesites, Splicetype_T *splicetypes,
		    int nsplicesites, bool dibasep, bool cmetp, int pairmax, int expected_pairlength) {
  
  List_T result = NULL, hitpairs = NULL, hits5 = NULL, hits3 = NULL, p;
  Stage3pair_T stage3pair;
  T this5, this3;
  char *queryuc_ptr_5, *queryuc_ptr_3;
  int usermax_level_5, usermax_level_3;
  int nconcordant, nsamechr;
  int cutoff_level_5, cutoff_level_3;
  int querylength5, querylength3, query5_lastpos, query3_lastpos;
  bool allvalidp5, allvalidp3;
  UINT4 *genome_blocks, *snp_blocks;
  int maxpairedpaths = 10*maxpaths; /* For computation, not for printing. */
  bool abort_pairing_p;

  querylength5 = Sequence_fulllength(queryseq5);
  querylength3 = Sequence_fulllength(queryseq3);

  *singlehits5 = *singlehits3 = (List_T) NULL;
  if (querylength5 <= INDEX1PART+2 || querylength3 <= INDEX1PART+2) {
    fprintf(stderr,"GSNAP cannot handle reads shorter than %d bp\n",INDEX1PART+2);
    return (List_T) NULL;
  } else if (querylength5 > MAX_QUERYLENGTH || querylength3 > MAX_QUERYLENGTH) {
    fprintf(stderr,"GSNAP cannot handle reads longer than %d bp.  Either recompile with a higher value of MAX_QUERYLENGTH, or consider using GMAP instead.\n",
	    MAX_QUERYLENGTH);
    return (List_T) NULL;

  } else {
    if (usermax_level_float < 0.0) {
      usermax_level_5 = usermax_level_3 = -1;
    } else if (usermax_level_float > 0.0 && usermax_level_float < 1.0) {
      usermax_level_5 = (int) rint(usermax_level_float * (double) querylength5);
      usermax_level_3 = (int) rint(usermax_level_float * (double) querylength3);
    } else {
      usermax_level_5 = usermax_level_3 = (int) usermax_level_float;
    }

    this5 = Stage1_new(querylength5);
    this3 = Stage1_new(querylength3);
    queryuc_ptr_5 = Sequence_fullpointer_uc(queryseq5);
    queryuc_ptr_3 = Sequence_fullpointer_uc(queryseq3);
    query5_lastpos = querylength5 - INDEX1PART;
    query3_lastpos = querylength3 - INDEX1PART;
  
    if (read_oligos(&allvalidp5,this5,queryuc_ptr_5,querylength5,query5_lastpos,dibasep,cmetp) == 0 || 
	read_oligos(&allvalidp3,this3,queryuc_ptr_3,querylength3,query3_lastpos,dibasep,cmetp) == 0) {
      debug(printf("Aborting because no hits found anywhere\n"));
      Stage1_free(&this3,querylength3);
      Stage1_free(&this5,querylength5);
      return (List_T) NULL;

    } else {
      genome_blocks = Genome_blocks(genome);
      snp_blocks = genomealt ? Genome_blocks(genomealt) : NULL;
      abort_pairing_p = false;
      hitpairs = align_pair(&abort_pairing_p,&nconcordant,&nsamechr,&cutoff_level_5,&cutoff_level_3,&hits5,&hits3,
			    this5,this3,queryuc_ptr_5,queryuc_ptr_3,
			    querylength5,querylength3,query5_lastpos,query3_lastpos,
			    indexdb,indexdb2,indexdb_size_threshold,
			    chromosome_iit,genome_blocks,snp_blocks,genome,genomealt,floors_array,
			    splicesites,splicetypes,nsplicesites,
			    usermax_level_5,usermax_level_3,subopt_levels,masktype,indel_penalty,
			    maxpaths,maxchimerapaths,knownsplicingp,novelsplicingp,
			    canonicalp,shortsplicedist,
			    localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
			    min_distantsplicing_end_matches,min_distantsplicing_identity,
			    max_middle_insertions,max_middle_deletions,
			    allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			    allvalidp5,allvalidp3,trim_ends_p,dibasep,cmetp,
			    pairmax,expected_pairlength,maxpairedpaths);
      Stage1_free(&this3,querylength3);
      Stage1_free(&this5,querylength5);
  
      if (abort_pairing_p == true) {
	/* Clean up all previous calculations */
	for (p = hitpairs; p != NULL; p = List_next(p)) {
	  stage3pair = (Stage3pair_T) List_head(p);
	  Stage3pair_free(&stage3pair);
	}
	List_free(&hitpairs);
	stage3list_gc(&hits3);
	stage3list_gc(&hits5);

	/* Treat 5' end as a single end */
	this5 = Stage1_new(querylength5);
	if (read_oligos(&allvalidp5,this5,queryuc_ptr_5,querylength5,query5_lastpos,dibasep,cmetp) == 0) {
	  debug(printf("Aborting because no hits found anywhere\n"));
	  *singlehits5 = (List_T) NULL;
	} else {
	  *singlehits5 = align_end(&cutoff_level_5,this5,queryuc_ptr_5,querylength5,query5_lastpos,
				   indexdb,indexdb2,indexdb_size_threshold,
				   chromosome_iit,genome_blocks,snp_blocks,
				   genome,genomealt,floors_array,
				   splicesites,splicetypes,nsplicesites,
				   usermax_level_5,subopt_levels,masktype,indel_penalty,
				   maxpaths,maxchimerapaths,knownsplicingp,novelsplicingp,
				   canonicalp,shortsplicedist,
				   localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
				   min_distantsplicing_end_matches,min_distantsplicing_identity,
				   max_middle_insertions,max_middle_deletions,
				   allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
				   allvalidp5,trim_ends_p,dibasep,cmetp);
	  *singlehits5 = Stage3_remove_duplicates(Stage3_optimal_score(*singlehits5,cutoff_level_5,subopt_levels));

	  if (List_length(*singlehits5) > 1 && geneprob_iit != NULL) {
	    Stage3_geneprob_eval(*singlehits5,geneprob_iit,chromosome_iit,queryseq5);
	  }
	}
	Stage1_free(&this5,querylength5);

	/* Treat 3' end as a single end */
	this3 = Stage1_new(querylength3);
	if (read_oligos(&allvalidp3,this3,queryuc_ptr_3,querylength3,query3_lastpos,dibasep,cmetp) == 0) {
	  debug(printf("Aborting because no hits found anywhere\n"));
	  *singlehits3 = (List_T) NULL;
	} else {
	  *singlehits3 = align_end(&cutoff_level_3,this3,queryuc_ptr_3,querylength3,query3_lastpos,
				   indexdb,indexdb2,indexdb_size_threshold,
				   chromosome_iit,genome_blocks,snp_blocks,
				   genome,genomealt,floors_array,
				   splicesites,splicetypes,nsplicesites,
				   usermax_level_3,subopt_levels,masktype,indel_penalty,
				   maxpaths,maxchimerapaths,knownsplicingp,novelsplicingp,
				   canonicalp,shortsplicedist,
				   localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
				   min_distantsplicing_end_matches,min_distantsplicing_identity,
				   max_middle_insertions,max_middle_deletions,
				   allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
				   allvalidp3,trim_ends_p,dibasep,cmetp);
	  *singlehits3 = Stage3_remove_duplicates(Stage3_optimal_score(*singlehits3,cutoff_level_3,subopt_levels));

	  if (List_length(*singlehits3) > 1 && geneprob_iit != NULL) {
	    Stage3_geneprob_eval(*singlehits3,geneprob_iit,chromosome_iit,queryseq3);
	  }
	}
	Stage1_free(&this3,querylength3);

	return (List_T) NULL;

      } else {
	if (hitpairs == NULL) {
	  result = (List_T) NULL;
	} else if (nconcordant > 0) {
	  result = Stage3pair_remove_duplicates(Stage3pair_optimal_score(Stage3pair_remove_samechr(hitpairs),
									 cutoff_level_5+cutoff_level_3,subopt_levels));
	  result = Stage3pair_sort_distance(result);
	} else {
	  result = Stage3pair_remove_duplicates(Stage3pair_optimal_score(hitpairs,cutoff_level_5+cutoff_level_3,subopt_levels));
	  result = Stage3pair_sort_distance(result);
	}

	if (result != NULL) {
	  stage3list_gc(&hits3);
	  stage3list_gc(&hits5);
	  *singlehits5 = *singlehits3 = (List_T) NULL;

	  if (List_length(result) > 1 && geneprob_iit != NULL) {
	    Stage3pair_geneprob_eval(result,geneprob_iit,chromosome_iit,queryseq5,queryseq3);
	  }
	  return result;

	} else {
	  /* printf("%d %d\n",List_length(hits5),List_length(hits3)); */
	  *singlehits5 = Stage3_remove_duplicates(Stage3_optimal_score(hits5,cutoff_level_5,subopt_levels));
	  *singlehits3 = Stage3_remove_duplicates(Stage3_optimal_score(hits3,cutoff_level_3,subopt_levels));
	  if (List_length(*singlehits5) > 1 && geneprob_iit != NULL) {
	    Stage3_geneprob_eval(*singlehits5,geneprob_iit,chromosome_iit,queryseq5);
	  }
	  if (List_length(*singlehits3) > 1 && geneprob_iit != NULL) {
	    Stage3_geneprob_eval(*singlehits3,geneprob_iit,chromosome_iit,queryseq3);
	  }

	  return (List_T) NULL;
	}
      }
    }
  }
}

