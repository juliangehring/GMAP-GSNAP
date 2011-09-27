static char rcsid[] = "$Id: stage1hr.c 36831 2011-03-19 00:30:36Z twu $";
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
#include "assert.h"
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
#include "maxent_hr.h"
#include "iitdef.h"
#include "interval.h"
#include "spanningelt.h"
#include "cmet.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#endif


#define A_CHAR 0x0
#define C_CHAR 0x1
#define G_CHAR 0x2
#define T_CHAR 0x3

#define NOT_APPLICABLE true

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
#define MAX_LOCALSPLICING_POTENTIAL 50


/* Overall flow */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* identify_segments */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
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

/* solve_doublesplice */
#ifdef DEBUG4D
#define debug4d(x) x
#else
#define debug4d(x)
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

/* Spliceends novel distant (long) details */
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

/* Terminals */
#ifdef DEBUG4T
#define debug4t(x) x
#else
#define debug4t(x)
#endif

/* Short overlaps */ 
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



typedef struct Segment_T *Segment_T;
struct Segment_T {
  int splicesites_i;		/* if no splicesites_iit, then splicesites_i is zero */
  Genomicpos_T diagonal;
  Genomicpos_T chroffset;
  Chrnum_T chrnum;

  int querypos5;
  int querypos3;

  int floor;
  int floor_xfirst;
  int floor_xlast;

  int floor_left;
  int floor_right;

  int leftmost;			/* For segmenti of local splice */
  int rightmost;		/* For segmentj of local splice */

#if 0
  int leftspan;			/* For segmentm of double splice */
  int rightspan;
#endif
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


#define KEEP_FLOORS 1

static Floors_T
Floors_new_standard (int querylength, int max_end_insertions) {
  Floors_T new;
  int query_lastpos, pos, from, to;
  int extra = 1 + max_end_insertions + max_end_insertions;

  query_lastpos = querylength - INDEX1PART;

#ifdef KEEP_FLOORS
#ifdef LEAKCHECK
  Mem_leak_check_deactivate();
#endif
#endif

  new = (Floors_T) MALLOC(sizeof(*new));
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

#ifdef KEEP_FLOORS
#ifdef LEAKCHECK
  Mem_leak_check_activate();
#endif
#endif

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
  Floors_T new;
  int query_lastpos, querypos, pos, from, to;
  int prev;
  int extra = 1 + max_end_insertions + max_end_insertions;

  query_lastpos = querylength - INDEX1PART;
  new = (Floors_T) MALLOC(sizeof(*new));
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


static bool free_positions_p;

void
Stage1_init_positions_free (bool positions_fileio_p) {
  if (positions_fileio_p == true) {
    free_positions_p = true;
  } else {
    free_positions_p = false;
  }
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

    if (free_positions_p == true) {
      for (i = -2; i < querylength; i++) {
	if ((*old)->plus_retrievedp[i] == true) {
	  FREE((*old)->plus_positions[i]);
	}
	if ((*old)->minus_retrievedp[i] == true) {
	  FREE((*old)->minus_positions[i]);
	}
      }
    } else {
      for (i = -2; i < querylength; i++) {
	if ((*old)->plus_allocp[i] == true) {
	  FREE((*old)->plus_positions[i]);
	}
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
  batch->diagonal = *diagonals;	/* Already in correct endianness */
  batch->npositions = ndiagonals;

  while (batch->npositions > 0 && batch->diagonal < querylength) {
    debug11(printf("Eliminating diagonal %u as straddling beginning of genome (Batch_init)\n",batch->diagonal));
    /* positions are really diagonals, already in correct endianness */
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
report_perfect_segment (int *found_score, List_T hits, Genomicpos_T left, Genomicpos_T diagonal,
			Chrnum_T chrnum, Genomicpos_T chroffset,
			char *query, char *queryptr, int querylength,
			Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			int nmisses_allowed, bool dibasep, bool cmetp, bool plusp) {
  Stage3_T hit;
  int nmismatches, ncolordiffs;

  if (cmetp) {
    /* Count actual number of mismatches.  May not be a perfect segment. */
    nmismatches = Genome_count_mismatches_limit(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
						left,/*pos5*/0,/*pos3*/querylength,
						/*max_mismatches_allowed*/nmisses_allowed,
						dibasep,/*cmetp*/true,plusp);
    if (nmismatches > nmisses_allowed) {
      return hits;
    } else {
      /* Don't use Stage3_new_exact, because need to mark mismatches */
      if ((hit = Stage3_new_substitution(&(*found_score),nmismatches,/*ncolordiffs*/0,
					 left,/*genomiclength*/querylength,
					 query_compress,genome_blocks,snp_blocks,
					 plusp,query,chrnum,chroffset,
					 dibasep,/*cmetp*/true)) == NULL) {
	return hits;
      } else {
	return List_push(hits,(void *) hit);
      }
    }

  } else if (snp_blocks) {
    if ((hit = Stage3_new_substitution(&(*found_score),/*nmismatches*/0,/*ncolordiffs*/0,
				       left,/*genomiclength*/querylength,
				       query_compress,genome_blocks,snp_blocks,
				       plusp,query,chrnum,chroffset,
				       dibasep,/*cmetp*/false)) == NULL) {
      return hits;
    } else {
      return List_push(hits,(void *) hit);
    }

  } else if (dibasep) {
    Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
				      left,/*pos5*/0,/*pos3*/querylength,dibasep,/*cmetp*/false,plusp);

    /* Need to fill buffer with nucleotide genome anyway */
    if ((hit = Stage3_new_substitution(&(*found_score),/*nmismatches*/0,ncolordiffs,
				       left,/*genomiclength*/querylength,
				       query_compress,genome_blocks,snp_blocks,
				       plusp,query,chrnum,chroffset,
				       dibasep,/*cmetp*/false)) == NULL) {
      return hits;
    } else {
      return List_push(hits,(void *) hit);
    }

  } else {
    if ((hit = Stage3_new_exact(&(*found_score),left,/*genomiclength*/querylength,
				query_compress,genome_blocks,snp_blocks,
				plusp,chrnum,chroffset)) == NULL) {
      return hits;
    } else {
      return List_push(hits,(void *) hit);
    }

  }
}


/* Called only by exact/sub:1 procedures, so need to do Bigendian conversion */
#ifdef WORDS_BIGENDIAN
static int
binary_search_bigendian (int lowi, int highi, Genomicpos_T *positions, Genomicpos_T goal) {
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
#endif


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



/* Generalization of identify_exact_iter and identify_onemiss_iter */
static List_T
identify_multimiss_iter (int *found_score, Chrnum_T *chrnum, Genomicpos_T *chroffset, Genomicpos_T *chrhigh,
			 List_T hits, Genomicpos_T goal, List_T prev, int *nempty,
			 int *global_miss_querypos5, int *global_miss_querypos3,
			 char *query, char *queryptr, int querylength,
			 Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			 IIT_T chromosome_iit, bool plusp,
			 int nmisses_allowed, int nmisses_seen, int miss_querypos5, int miss_querypos3,
			 bool dibasep, bool cmetp) {
  List_T spanningset;
  Stage3_T hit;
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
	debug7(printf("  (%d>>",elt->intersection_ndiagonals));
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
	    j = binary_search_bigendian(j >> 1,elt->partner_npositions,elt->partner_positions,local_goal);
	  } else {
	    j = binary_search_bigendian(j >> 1,j,elt->partner_positions,local_goal);
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
	    j = binary_search_bigendian(j >> 1,elt->npositions,elt->positions,local_goal);
	  } else {
	    j = binary_search_bigendian(j >> 1,j,elt->positions,local_goal);
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
    left = goal - querylength;
    if (goal > *chrhigh) {
      *chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,left,left);
      IIT_interval_bounds(&(*chroffset),&(*chrhigh),chromosome_iit,*chrnum);
      *chrhigh += 1U;
    }
    debug(printf("Reporting perfect segment at left %u and diagonal %u, with chroffset %u and chrhigh %u\n",
		 left,goal,*chroffset,*chrhigh));
    if (goal > *chrhigh) {
      /* Query goes over end of chromosome */
      debug(printf("  Ignore: goes over end of chromosome\n"));
      return hits;
    } else {
      return report_perfect_segment(&(*found_score),hits,left,/*diagonal*/goal,*chrnum,*chroffset,
				    query,queryptr,querylength,
				    query_compress,genome_blocks,snp_blocks,
				    nmisses_allowed,dibasep,cmetp,plusp);
    }
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
	nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
							left,/*pos5*/0,/*pos3*/querylength,dibasep,cmetp,plusp);
      } else {
	debug7(printf("  Testing in query bounds %d..%d\n",miss_querypos5,miss_querypos3));
	nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
							left,/*pos5*/miss_querypos5,/*pos3*/miss_querypos3,
							dibasep,cmetp,plusp);

      }
      debug7(printf("nmismatches = %d (vs %d misses allowed)\n",nmismatches,nmisses_allowed));

      if (nmismatches > nmisses_allowed) {
	debug7(printf("Result: too many mismatches\n"));
	return hits;
      } else {
	debug7(printf("Result: successful hit saved\n"));
	debug(printf("Reporting hit with %d mismatches\n",nmismatches));
	if ((hit = Stage3_new_substitution(&(*found_score),nmismatches,ncolordiffs,
					   left,/*genomiclength*/querylength,
					   query_compress,genome_blocks,snp_blocks,
					   plusp,query,*chrnum,*chroffset,
					   dibasep,cmetp)) == NULL) {
	  return hits;
	} else {
	  return List_push(hits,(void *) hit);
	}
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
find_spanning_exact_matches (int *found_score, T this, char *queryuc_ptr, char *queryrc,
			     int querylength, int query_lastpos, Indexdb_T indexdb, Indexdb_T indexdb2,
			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     UINT4 *genome_blocks, UINT4 *snp_blocks,
			     IIT_T chromosome_iit, bool dibasep, bool cmetp) {
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
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,Bigendian_convert_uint(*positions0),diagterm0));
	diagonal0 = Bigendian_convert_uint(*positions0++) + diagterm0;
#else
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,(*positions0),diagterm0));
	diagonal0 = (*positions0++) + diagterm0;
#endif
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/true,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3,
				       dibasep,cmetp);
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

      /* diagonals0 is now in correct endianness */
      diagonals0 = Spanningelt_diagonals(&ndiagonals0,(Spanningelt_T) sorted->first,&elt_miss_querypos5,&elt_miss_querypos3);
      spanningset = List_push(List_copy(sorted->rest),(void **) NULL); /* Add a dummy list elt to front */
      nempty = 0;
      global_miss_querypos5 = querylength;
      global_miss_querypos3 = 0;

      while (--ndiagonals0 >= 0 && nempty == 0) {
	debug7(printf("diag0 %d:%u advancing\n",ndiagonals0,(*diagonals0)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,*diagonals0++,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/true,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3,
				       dibasep,cmetp);
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
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,Bigendian_convert_uint(*positions0),diagterm0));
	diagonal0 = Bigendian_convert_uint(*positions0++) + diagterm0;
#else
	debug7(printf("diag0 %d:%u+%d advancing\n",npositions0,(*positions0),diagterm0));
	diagonal0 = (*positions0++) + diagterm0;
#endif
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/false,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3,
				       dibasep,cmetp);
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

      /* diagonals0 is now in correct endianness */
      diagonals0 = Spanningelt_diagonals(&ndiagonals0,(Spanningelt_T) sorted->first,&elt_miss_querypos5,&elt_miss_querypos3);
      spanningset = List_push(List_copy(sorted->rest),(void **) NULL); /* Add a dummy list elt to front */
      nempty = 0;
      global_miss_querypos5 = querylength;
      global_miss_querypos3 = 0;

      while (--ndiagonals0 >= 0 && nempty == 0) {
	debug7(printf("diag0 %d:%u advancing\n",ndiagonals0,(*diagonals0)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,*diagonals0++,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/false,/*nmisses_allowed*/0,
				       /*nmisses_seen*/0,global_miss_querypos5,global_miss_querypos3,
				       dibasep,cmetp);
      }
      List_free(&spanningset);
    }
  }

  return hits;
}


static List_T
find_spanning_onemiss_matches (int *found_score, List_T hits, T this,
			       char *queryuc_ptr, char *queryrc, int querylength,
			       Compress_T query_compress_fwd, Compress_T query_compress_rev,
			       UINT4 *genome_blocks, UINT4 *snp_blocks,
			       IIT_T chromosome_iit, bool dibasep, bool cmetp) {
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

    /* diagonals0 and diagonals1 are now in correct endianness */
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
	debug7(printf("diag0 %d:%u advancing\n",ndiagonals0,diagonal0));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/true,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3,
				       dibasep,cmetp);
	++diagonals0;
	--ndiagonals0;

      } else if (diagonal1 < diagonal0) {
	debug7(printf("diag1 %d:%u advancing\n",ndiagonals1,diagonal1));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal1,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/true,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3,
				       dibasep,cmetp);
	++diagonals1;
	--ndiagonals1;

      } else {
	debug7(printf("diag0&1 %d:%u == %d:%u advancing\n",ndiagonals0,diagonal0,ndiagonals1,diagonal1));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/true,/*nmisses_allowed*/1,
				       /*nmisses_seen*/nempty,global_miss_querypos5,global_miss_querypos3,
				       dibasep,cmetp);
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
				     /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
				     /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				     chromosome_iit,/*plusp*/true,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3,
				     dibasep,cmetp);
    }

    while (--ndiagonals1 >= 0 && nempty == 0) {
      debug7(printf("diag1 %d:%u advancing\n",ndiagonals1+1,(*diagonals1)));
      hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,*diagonals1++,
				     /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				     /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
				     /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				     chromosome_iit,/*plusp*/true,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3,
				     dibasep,cmetp);
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

    /* diagonals0 and diagonals1 are now in correct endianness */
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
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/false,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3,
				       dibasep,cmetp);
	++diagonals0;
	--ndiagonals0;

      } else if (diagonal1 < diagonal0) {
	debug7(printf("diag1 %d:%u advancing\n",ndiagonals1,(*diagonals1)));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal1,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/false,/*nmisses_allowed*/1,
				       /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3,
				       dibasep,cmetp);
	++diagonals1;
	--ndiagonals1;

      } else {
	debug7(printf("diag0&1 %d:%u == %d:%u advancing\n",ndiagonals0,diagonal0,ndiagonals1,diagonal1));
	hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,diagonal0,
				       /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/false,/*nmisses_allowed*/1,
				       /*nmisses_seen*/nempty,global_miss_querypos5,global_miss_querypos3,
				       dibasep,cmetp);
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
				     /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
				     /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				     chromosome_iit,/*plusp*/false,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss1_querypos5,miss1_querypos3,
				     dibasep,cmetp);
    }

    while (--ndiagonals1 >= 0 && nempty == 0) {
      debug7(printf("diag1 %d:%u advancing\n",ndiagonals1+1,(*diagonals1)));
      hits = identify_multimiss_iter(&(*found_score),&chrnum,&chroffset,&chrhigh,hits,*diagonals1++,
				     /*prev*/spanningset,&nempty,&global_miss_querypos5,&global_miss_querypos3,
				     /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
				     /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				     chromosome_iit,/*plusp*/false,/*nmisses_allowed*/1,
				     /*nmisses_seen*/1+nempty,miss0_querypos5,miss0_querypos3,
				     dibasep,cmetp);
    }

    List_free(&spanningset);
  }

  return hits;
}


static List_T
find_spanning_multimiss_matches (int *found_score, List_T hits, T this, int nrequired,
				 char *queryuc_ptr, char *queryrc, int querylength,
				 Compress_T query_compress_fwd, Compress_T query_compress_rev,
				 UINT4 *genome_blocks, UINT4 *snp_blocks,
				 IIT_T chromosome_iit, int nmisses_allowed, bool dibasep, bool cmetp) {
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
	/* These positions are diagonals, and already in correct endianness */
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
					   /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
					   /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
					   chromosome_iit,/*plusp*/true,nmisses_allowed,
					   /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3,
					   dibasep,cmetp);
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
	  /* These positions are diagonals, and already in correct endianness */
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
				       /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,querylength,
				       /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/true,nmisses_allowed,
				       /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3,
				       dibasep,cmetp);
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
	/* These positions are diagonals, and already in correct endianness */
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
					   /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
					   /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
					   chromosome_iit,/*plusp*/false,nmisses_allowed,
					   /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3,
					   dibasep,cmetp);
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
	  /* These positions are diagonals, and already in correct endianness */
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
				       /*query*/queryuc_ptr,/*queryptr*/queryrc,querylength,
				       /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
				       chromosome_iit,/*plusp*/false,nmisses_allowed,
				       /*nmisses_seen*/nunion-count+nempty,global_miss_querypos5,global_miss_querypos3,
				       dibasep,cmetp);
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
find_complete_mm (int *found_score, List_T hits, struct Segment_T *segments, int nsegments,
		  int querylength, char *query, Compress_T query_compress,
		  UINT4 *genome_blocks, UINT4 *snp_blocks,
		  IIT_T chromosome_iit, int max_mismatches_allowed,
		  bool dibasep, bool cmetp, bool plusp) {
  Stage3_T hit;
  int nmismatches, ncolordiffs;
  Genomicpos_T left;
  Segment_T segmenti;

  for (segmenti = segments; segmenti < &(segments[nsegments]); segmenti++) {
    if (segmenti->floor <= max_mismatches_allowed) {
      left = segmenti->diagonal - querylength;
      nmismatches = Genome_count_mismatches_limit(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
						  left,/*pos5*/0,/*pos3*/querylength,
						  max_mismatches_allowed,dibasep,cmetp,plusp);
      if (nmismatches <= max_mismatches_allowed) {
	if ((hit = Stage3_new_substitution(&(*found_score),nmismatches,ncolordiffs,
					   left,/*genomiclength*/querylength,
					   query_compress,genome_blocks,snp_blocks,
					   plusp,query,segmenti->chrnum,segmenti->chroffset,
					   dibasep,cmetp)) != NULL) {
	  hits = List_push(hits,(void *) hit);
	}
      }
    }
  }

  return hits;
}



static struct Segment_T *
identify_all_segments (int *nsegments, Genomicpos_T **positions, int *npositions,
		       bool *omitted, int querylength, int query_lastpos,
		       IIT_T chromosome_iit, Floors_T floors, 
		       Genomicpos_T *splicesites, int nsplicesites, bool plusp) {
  struct Segment_T *segments = NULL;
  Batch_T batch, sentinel;
  struct Batch_T sentinel_struct, *batchpool;
  Batch_T *heap;
  int heapsize = 0;
  int parenti, smallesti, righti, i;
  int querypos, first_querypos, last_querypos;
  int halfquerylength, halfquery_lastpos;
  int floor, floor_xfirst, floor_xlast;
  int floor_left, floor_right, floor_incr;
  int *floors_from_neg3, *floors_from_xfirst, *floors_to_xlast, *floors_to_pos3;
  /* int exclude_xfirst, exclude_xlast; */
  Genomicpos_T diagonal, segment_left, last_diagonal, chroffset = 0U, chrhigh = 0U;
  Chrnum_T chrnum;
#ifdef HAVE_64_BIT
  UINT8 diagonal_add_querypos;
#endif
  int total_npositions = 0;
  int joffset = 0, j;
  Segment_T ptr;

  debug(printf("*** Starting identify_all_segments ***\n"));

  halfquerylength = querylength / 2;
  halfquery_lastpos = halfquerylength - INDEX1PART;

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
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      } else if (npositions[querypos] > 0) {
	debug1(printf("Adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
	batch = &(batchpool[i]);
	Batch_init(batch,querypos,/*diagterm*/querylength - querypos,positions[querypos],npositions[querypos],querylength);
	total_npositions += npositions[querypos];
	if (batch->npositions > 0) {
	  min_heap_insert(heap,&heapsize,batch);
	  i++;
	}
      } else {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      }
    }
  } else {
    for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
      if (omitted[querypos] == true) {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      } else if (npositions[querypos] > 0) {
	debug1(printf("Adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
	batch = &(batchpool[i]);
	Batch_init(batch,querypos,/*diagterm*/querypos + INDEX1PART,positions[querypos],npositions[querypos],querylength);
	total_npositions += npositions[querypos];
	if (batch->npositions > 0) {
	  min_heap_insert(heap,&heapsize,batch);
	  i++;
	}
      } else {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      }
    }
  }


  if (i == 0) {
    FREE(heap);
    FREE(batchpool);
    *nsegments = 0;
    return (struct Segment_T *) NULL;
  }

  /* Set up rest of heap */
  for (i = heapsize+1; i <= 2*heapsize+1; i++) {
    heap[i] = sentinel;
  }

  segments = (struct Segment_T *) CALLOC(total_npositions,sizeof(struct Segment_T));
  ptr = &(segments[0]);

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
  if (INDEX1PART /* +max_end_insertions */ > query_lastpos + 3) {
    floors_from_xfirst = floors->scorefrom[query_lastpos+3];
  } else {
    floors_from_xfirst = floors->scorefrom[INDEX1PART /* +max_end_insertions */];
  }
  if (query_lastpos-INDEX1PART /* -max_end_insertions */ < -3) {
    floors_to_xlast = floors->scoreto[-3];
  } else {
    floors_to_xlast = floors->scoreto[query_lastpos-INDEX1PART /* -max_end_insertions */];
  }
#endif
  floors_from_neg3 = floors->scorefrom[-3];
  floors_to_pos3 = floors->scoreto[query_lastpos+3];


  /* Initialize loop */
  batch = heap[1];
  first_querypos = last_querypos = querypos = batch->querypos;
  last_diagonal = diagonal = batch->diagonal;

  floor_incr = floors_from_neg3[first_querypos];
  floor = floor_incr;
  floor_xlast = floor_incr;
  floor_xfirst = floors_from_xfirst[first_querypos] /* floors->scorefrom[xfirst_from][first_querypos] */;

#ifdef OLD_FLOOR_ENDS
  if (querypos < halfquery_lastpos) {
    floor_left = floor_incr;
  } else {
    floor_left = floors->scorefrom[-3][halfquery_lastpos];
  }
  if (querypos < halfquerylength) {
    floor_right = floors->scorefrom[halfquerylength-3][query_lastpos];
  } else {
    floor_right = floors->scorefrom[halfquerylength-3][first_querypos];
  }
#else
  floor_left = floor_incr;
#ifdef DEBUG1
  floor_right = -99;
#endif
#endif


  debug1(printf("multiple_mm_%s, diagonal %u, querypos %d\n",plusp ? "plus" : "minus",diagonal,querypos));
  debug1(printf("first_querypos = %d => initial values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
	        first_querypos,floor,floor_xfirst,floor_xlast,floor_left,floor_right));

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

#ifdef OLD_FLOOR_ENDS
      /* Why is this here?  Just set floor_left at start and floor_right at end. */
      if (querypos < halfquery_lastpos) {
	floor_left += floor_incr;
      } else if (last_querypos < halfquery_lastpos) {
	/* Finish floor_left */
	floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+3];
      }
      if (querypos >= halfquerylength) {
	if (last_querypos < halfquerylength) {
	  /* Start floor_right */
	  floor_right = floors->scorefrom[halfquerylength-3][querypos];
	} else {
	  floor_right += floor_incr;
	}
      }
#endif

      debug1(printf("diagonal %u unchanged: last_querypos = %d, querypos = %d => floor increments by %d\n",
		    diagonal,last_querypos,querypos,floor_incr));
      debug1(printf("*multiple_mm_%s, diagonal %u, querypos %d, floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		    plusp ? "plus" : "minus",diagonal,querypos,floor,floor_xfirst,floor_xlast,floor_left,floor_right));
    } else {
      /* End of diagonal */
      floor_incr = floors_to_pos3[last_querypos]  /* floors->score[last_querypos][query_lastpos+3] */;
      floor += floor_incr;
      floor_xfirst += floor_incr;
      floor_xlast += floors_to_xlast[last_querypos];  /* floors->score[last_querypos][xlast_to]; */

#ifdef OLD_FLOOR_ENDS
      if (last_querypos < halfquery_lastpos) {
	floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+3];
	floor_right = floors->scorefrom[halfquerylength-3][query_lastpos+3];
      }
      if (last_querypos >= halfquerylength) {
	floor_right += floor_incr;
      }
#else
      floor_right = floor_incr;
#endif

      debug1(printf("new diagonal %u > last diagonal %u: last_querypos = %d => final values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		    diagonal,last_diagonal,last_querypos,floor,floor_xfirst,floor_xlast,floor_left,floor_right));

      if (last_diagonal > chrhigh) {
	/* update chromosome bounds, based on low end */
	chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,last_diagonal-querylength,last_diagonal-querylength);
	IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
	chrhigh += 1U;
      }
      if (last_diagonal <= chrhigh) { /* FORMULA for high position */
	/* position of high end is within current chromosome */
	debug1(printf("  => multiple_mm, diagonal %u, query %d..%d, chrbounds %u..%u, floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		      last_diagonal,first_querypos,last_querypos,chroffset,chrhigh,floor,floor_xfirst,floor_xlast,floor_left,floor_right));

	/* Save segment, but first advance splicesites past segment_left */
	if (nsplicesites > 0) {
	  segment_left = last_diagonal - querylength;
	  if (*splicesites < segment_left) {
	    j = 1;
	    while (j < nsplicesites && splicesites[j] < segment_left) {
	      j <<= 1;		/* gallop by 2 */
	    }
	    if (j >= nsplicesites) {
	      j = binary_search(j >> 1,nsplicesites,splicesites,segment_left);
	    } else {
	      j = binary_search(j >> 1,j,splicesites,segment_left);
	    }
	    joffset += j;
	    splicesites += j;
	    nsplicesites -= j;
	  }
	}

	/* Save segment */
	ptr->splicesites_i = joffset;
	ptr->diagonal = last_diagonal;
	ptr->chrnum = chrnum;
	ptr->chroffset = chroffset;
	ptr->querypos5 = first_querypos;
	ptr->querypos3 = last_querypos;
	ptr->floor = floor;
	ptr->floor_xfirst = floor_xfirst;
	ptr->floor_xlast = floor_xlast;
	ptr->floor_left = floor_left;
	ptr->floor_right = floor_right;
	ptr->leftmost = ptr->rightmost = -1;
#if 0
	ptr->leftspan = ptr->rightspan = -1;
#endif
	ptr++;
      }

      /* Prepare next diagonal */
      first_querypos = querypos;
      last_diagonal = diagonal;
      floor_incr = floors_from_neg3[first_querypos] /* floors->score[-3][first_querypos] */;
      floor = floor_incr;
      floor_xlast = floor_incr;
      floor_xfirst = floors_from_xfirst[first_querypos];  /* floors->score[xfirst_from][first_querypos]; */

#ifdef OLD_FLOOR_ENDS
      if (querypos < halfquery_lastpos) {
	floor_left = floor_incr;
      } else {
	floor_left = floors->scorefrom[-3][halfquery_lastpos];
      }
      if (querypos < halfquerylength) {
	floor_right = floors->scorefrom[halfquerylength-3][query_lastpos];
      } else {
	floor_right = floors->scorefrom[halfquerylength-3][first_querypos];
      }
#else
      floor_left = floor_incr;
#ifdef DEBUG1
      floor_right = -99;	/* For debugging output */
#endif
#endif

      debug1(printf("*multiple_mm_%s, diagonal %u, querypos %d\n",plusp ? "plus" : "minus",diagonal,querypos));
      debug1(printf("start of diagonal %u, first_querypos = %d => initial values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		    diagonal,first_querypos,floor,floor_xfirst,floor_xlast,floor_left,floor_right));
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
  floor_incr = floors_to_pos3[last_querypos];   /* floors->score[last_querypos][query_lastpos+3]; */
  floor += floor_incr;
  floor_xfirst += floor_incr;
  floor_xlast += floors_to_xlast[last_querypos];  /* floors->score[last_querypos][xlast_to]; */

#ifdef OLD_FLOOR_ENDS
  if (last_querypos < halfquery_lastpos) {
    floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+3];
    floor_right = floors->scorefrom[halfquerylength-3][query_lastpos+3];
  }
  if (last_querypos >= halfquerylength) {
    floor_right += floor_incr;
  }
#else
  floor_right = floor_incr;
#endif
  
  debug1(printf("no more diagonals: last_querypos = %d => terminal values: floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		last_querypos,floor,floor_xfirst,floor_xlast,floor_left,floor_right));

  if (last_diagonal > chrhigh) {
    /* update chromosome bounds, based on low end */
    chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,last_diagonal-querylength,last_diagonal-querylength);
    IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
    chrhigh += 1U;
  }
  if (last_diagonal <= chrhigh) { /* FORMULA for high position */
    /* position of high end is within current chromosome */
    debug1(printf("  => multiple_mm, diagonal %u, query %d..%d, chrbounds %u..%u, floor %d, floor_xfirst %d, floor_xlast %d, floor_left %d, floor_right %d\n",
		  last_diagonal,first_querypos,last_querypos,chroffset,chrhigh,floor_xfirst,floor_xlast,floor_left,floor_right));

    /* Save segment, but first advance splicesites past segment_left */
    if (nsplicesites > 0) {
      segment_left = last_diagonal - querylength;
      if (*splicesites < segment_left) {
	j = 1;
	while (j < nsplicesites && splicesites[j] < segment_left) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= nsplicesites) {
	  j = binary_search(j >> 1,nsplicesites,splicesites,segment_left);
	} else {
	  j = binary_search(j >> 1,j,splicesites,segment_left);
	}
	joffset += j;
	splicesites += j;
	nsplicesites -= j;
      }
    }

    /* Save segment */
    ptr->splicesites_i = joffset;
    ptr->diagonal = last_diagonal;
    ptr->chrnum = chrnum;
    ptr->chroffset = chroffset;
    ptr->querypos5 = first_querypos;
    ptr->querypos3 = last_querypos;
    ptr->floor = floor;
    ptr->floor_xfirst = floor_xfirst;
    ptr->floor_xlast = floor_xlast;
    ptr->floor_left = floor_left;
    ptr->floor_right = floor_right;
    ptr->leftmost = ptr->rightmost = -1;
#if 0
    ptr->leftspan = ptr->rightspan = -1;
#endif
    ptr++;
  }

  FREE(heap);
  FREE(batchpool);

  /* Note: segments is in descending diagonal order.  Will need to
     reverse before solving middle deletions */

  *nsegments = ptr - segments;
  debug1(printf("nsegments = %d\n",*nsegments));
  debug(printf("nsegments = %d (total_npositions = %d)\n",*nsegments,total_npositions));

  return segments;
}


/* Specialized version of identify_all_segments that stores only floor_left and floor_right */
static struct Segment_T *
identify_all_segments_for_terminals (int *nsegments, Genomicpos_T **positions, int *npositions,
				     bool *omitted, int querylength, int query_lastpos,
				     IIT_T chromosome_iit, Floors_T floors, 
				     int max_mismatches_allowed, bool plusp) {
  struct Segment_T *segments = NULL;
  Batch_T batch, sentinel;
  struct Batch_T sentinel_struct, *batchpool;
  Batch_T *heap;
  int heapsize = 0;
  int parenti, smallesti, righti, i;
  int querypos, first_querypos, last_querypos;
  int halfquerylength, halfquery_lastpos;
  int floor_left, floor_right, floor_incr;
  int *floors_from_neg3, *floors_from_xfirst, *floors_to_xlast, *floors_to_pos3;
  /* int exclude_xfirst, exclude_xlast; */
  Genomicpos_T diagonal, last_diagonal, chroffset = 0U, chrhigh = 0U;
  Chrnum_T chrnum;
#ifdef HAVE_64_BIT
  UINT8 diagonal_add_querypos;
#endif
  int total_npositions = 0;
  Segment_T ptr;

  debug(printf("*** Starting identify_all_segments ***\n"));

  halfquerylength = querylength / 2;
  halfquery_lastpos = halfquerylength - INDEX1PART;

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
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      } else if (npositions[querypos] > 0) {
	debug1(printf("Adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
	batch = &(batchpool[i]);
	Batch_init(batch,querypos,/*diagterm*/querylength - querypos,positions[querypos],npositions[querypos],querylength);
	total_npositions += npositions[querypos];
	if (batch->npositions > 0) {
	  min_heap_insert(heap,&heapsize,batch);
	  i++;
	}
      } else {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      }
    }
  } else {
    for (querypos = 0, i = 0; querypos <= query_lastpos; querypos++) {
      if (omitted[querypos] == true) {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      } else if (npositions[querypos] > 0) {
	debug1(printf("Adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
	batch = &(batchpool[i]);
	Batch_init(batch,querypos,/*diagterm*/querypos + INDEX1PART,positions[querypos],npositions[querypos],querylength);
	total_npositions += npositions[querypos];
	if (batch->npositions > 0) {
	  min_heap_insert(heap,&heapsize,batch);
	  i++;
	}
      } else {
	debug1(printf("Not adding batch for querypos %d with %d positions, omitted %d\n",
		      querypos,npositions[querypos],omitted[querypos]));
      }
    }
  }


  if (i == 0) {
    FREE(heap);
    FREE(batchpool);
    *nsegments = 0;
    return (struct Segment_T *) NULL;
  }

  /* Set up rest of heap */
  for (i = heapsize+1; i <= 2*heapsize+1; i++) {
    heap[i] = sentinel;
  }

  segments = (struct Segment_T *) CALLOC(total_npositions,sizeof(struct Segment_T));
  ptr = &(segments[0]);

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
  if (INDEX1PART /* +max_end_insertions */ > query_lastpos + 3) {
    floors_from_xfirst = floors->scorefrom[query_lastpos+3];
  } else {
    floors_from_xfirst = floors->scorefrom[INDEX1PART /* +max_end_insertions */];
  }
  if (query_lastpos-INDEX1PART /* -max_end_insertions */ < -3) {
    floors_to_xlast = floors->scoreto[-3];
  } else {
    floors_to_xlast = floors->scoreto[query_lastpos-INDEX1PART /* -max_end_insertions */];
  }
#endif
  floors_from_neg3 = floors->scorefrom[-3];
  floors_to_pos3 = floors->scoreto[query_lastpos+3];


  /* Initialize loop */
  batch = heap[1];
  first_querypos = last_querypos = querypos = batch->querypos;
  last_diagonal = diagonal = batch->diagonal;

  floor_incr = floors_from_neg3[first_querypos];
#if 0
  floor = floor_incr;
  floor_xlast = floor_incr;
  floor_xfirst = floors_from_xfirst[first_querypos] /* floors->scorefrom[xfirst_from][first_querypos] */;
#endif

#ifdef OLD_FLOOR_ENDS
  if (querypos < halfquery_lastpos) {
    floor_left = floor_incr;
  } else {
    floor_left = floors->scorefrom[-3][halfquery_lastpos];
  }
  if (querypos < halfquerylength) {
    floor_right = floors->scorefrom[halfquerylength-3][query_lastpos];
  } else {
    floor_right = floors->scorefrom[halfquerylength-3][first_querypos];
  }
#else
  floor_left = floor_incr;
#ifdef DEBUG1
  floor_right = -99;
#endif
#endif


  debug1(printf("multiple_mm_%s, diagonal %u, querypos %d\n",plusp ? "plus" : "minus",diagonal,querypos));
  debug1(printf("first_querypos = %d => initial values: floor_left %d, floor_right %d\n",
	        first_querypos,floor_left,floor_right));

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
#if 0
      floor += floor_incr;
      floor_xfirst += floor_incr;
      floor_xlast += floor_incr;
#endif

#ifdef OLD_FLOOR_ENDS
      /* Why is this here?  Just set floor_left at start and floor_right at end. */
      if (querypos < halfquery_lastpos) {
	floor_left += floor_incr;
      } else if (last_querypos < halfquery_lastpos) {
	/* Finish floor_left */
	floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+3];
      }
      if (querypos >= halfquerylength) {
	if (last_querypos < halfquerylength) {
	  /* Start floor_right */
	  floor_right = floors->scorefrom[halfquerylength-3][querypos];
	} else {
	  floor_right += floor_incr;
	}
      }
#endif

      debug1(printf("diagonal %u unchanged: last_querypos = %d, querypos = %d => floor increments by %d\n",
		    diagonal,last_querypos,querypos,floor_incr));
      debug1(printf("*multiple_mm_%s, diagonal %u, querypos %d, floor_left %d, floor_right %d\n",
		    plusp ? "plus" : "minus",diagonal,querypos,floor_left,floor_right));
    } else {
      /* End of diagonal */
      floor_incr = floors_to_pos3[last_querypos]  /* floors->score[last_querypos][query_lastpos+3] */;
#if 0
      floor += floor_incr;
      floor_xfirst += floor_incr;
      floor_xlast += floors_to_xlast[last_querypos];  /* floors->score[last_querypos][xlast_to]; */
#endif

#ifdef OLD_FLOOR_ENDS
      if (last_querypos < halfquery_lastpos) {
	floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+3];
	floor_right = floors->scorefrom[halfquerylength-3][query_lastpos+3];
      }
      if (last_querypos >= halfquerylength) {
	floor_right += floor_incr;
      }
#else
      floor_right = floor_incr;
#endif

      debug1(printf("new diagonal %u > last diagonal %u: last_querypos = %d => final values: floor_left %d, floor_right %d\n",
		    diagonal,last_diagonal,last_querypos,floor_left,floor_right));

      if (last_diagonal > chrhigh) {
	/* update chromosome bounds, based on low end */
	chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,last_diagonal-querylength,last_diagonal-querylength);
	IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
	chrhigh += 1U;
      }
      if (last_diagonal <= chrhigh) { /* FORMULA for high position */
	/* position of high end is within current chromosome */
	debug1(printf("  => multiple_mm, diagonal %u, query %d..%d, chrbounds %u..%u, floor_left %d, floor_right %d\n",
		      last_diagonal,first_querypos,last_querypos,chroffset,chrhigh,floor_left,floor_right));
	if (floor_left <= max_mismatches_allowed || floor_right <= max_mismatches_allowed) {
	  /* Save segment */
	  ptr->diagonal = last_diagonal;
	  ptr->chrnum = chrnum;
	  ptr->chroffset = chroffset;
	  ptr->querypos5 = first_querypos;
	  ptr->querypos3 = last_querypos;
#if 0
	  ptr->floor = floor;
	  ptr->floor_xfirst = floor_xfirst;
	  ptr->floor_xlast = floor_xlast;
#endif
	  ptr->floor_left = floor_left;
	  ptr->floor_right = floor_right;
#if 0
	  ptr->leftmost = ptr->rightmost = -1;
	  ptr->leftspan = ptr->rightspan = -1;
#endif
	  ptr++;
	}
      }

      /* Prepare next diagonal */
      first_querypos = querypos;
      last_diagonal = diagonal;
      floor_incr = floors_from_neg3[first_querypos] /* floors->score[-3][first_querypos] */;
#if 0
      floor = floor_incr;
      floor_xlast = floor_incr;
      floor_xfirst = floors_from_xfirst[first_querypos];  /* floors->score[xfirst_from][first_querypos]; */
#endif

#ifdef OLD_FLOOR_ENDS
      if (querypos < halfquery_lastpos) {
	floor_left = floor_incr;
      } else {
	floor_left = floors->scorefrom[-3][halfquery_lastpos];
      }
      if (querypos < halfquerylength) {
	floor_right = floors->scorefrom[halfquerylength-3][query_lastpos];
      } else {
	floor_right = floors->scorefrom[halfquerylength-3][first_querypos];
      }
#else
      floor_left = floor_incr;
#ifdef DEBUG1
      floor_right = -99;
#endif
#endif

      debug1(printf("*multiple_mm_%s, diagonal %u, querypos %d\n",plusp ? "plus" : "minus",diagonal,querypos));
      debug1(printf("start of diagonal %u, first_querypos = %d => initial values: floor_left %d, floor_right %d\n",
		    diagonal,first_querypos,floor_left,floor_right));
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
  floor_incr = floors_to_pos3[last_querypos];   /* floors->score[last_querypos][query_lastpos+3]; */
#if 0
  floor += floor_incr;
  floor_xfirst += floor_incr;
  floor_xlast += floors_to_xlast[last_querypos];  /* floors->score[last_querypos][xlast_to]; */
#endif

#ifdef OLD_FLOOR_ENDS
  if (last_querypos < halfquery_lastpos) {
    floor_left += floors->scorefrom[last_querypos][halfquery_lastpos+3];
    floor_right = floors->scorefrom[halfquerylength-3][query_lastpos+3];
  }
  if (last_querypos >= halfquerylength) {
    floor_right += floor_incr;
  }
#else
  floor_right = floor_incr;
#endif

  
  debug1(printf("no more diagonals: last_querypos = %d => terminal values: floor_left %d, floor_right %d\n",
		last_querypos,floor_left,floor_right));

  if (last_diagonal > chrhigh) {
    /* update chromosome bounds, based on low end */
    chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,last_diagonal-querylength,last_diagonal-querylength);
    IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
    chrhigh += 1U;
  }
  if (last_diagonal <= chrhigh) { /* FORMULA for high position */
    /* position of high end is within current chromosome */
    debug1(printf("  => multiple_mm, diagonal %u, query %d..%d, chrbounds %u..%u, floor_left %d, floor_right %d\n",
		  last_diagonal,first_querypos,last_querypos,chroffset,chrhigh,floor_left,floor_right));
    if (floor_left <= max_mismatches_allowed || floor_right <= max_mismatches_allowed) {
      /* Save segment */
      ptr->diagonal = last_diagonal;
      ptr->chrnum = chrnum;
      ptr->chroffset = chroffset;
      ptr->querypos5 = first_querypos;
      ptr->querypos3 = last_querypos;
#if 0
      ptr->floor = floor;
      ptr->floor_xfirst = floor_xfirst;
      ptr->floor_xlast = floor_xlast;
#endif
      ptr->floor_left = floor_left;
      ptr->floor_right = floor_right;
#if 0
      ptr->leftmost = ptr->rightmost = -1;
      ptr->leftspan = ptr->rightspan = -1;
#endif
      ptr++;
    }
  }

  FREE(heap);
  FREE(batchpool);

  /* Note: segments is in descending diagonal order.  Will need to
     reverse before solving middle deletions */

  *nsegments = ptr - segments;
  debug1(printf("nsegments = %d\n",*nsegments));
  debug(printf("nsegments = %d (total_npositions = %d)\n",*nsegments,total_npositions));

  return segments;
}



/* indels is positive here */
static List_T
solve_middle_insertion (bool *foundp, int *found_score, List_T hits, Segment_T ptr, int indels,
			Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			char *query, char *queryptr, int querylength,
			int min_indel_end_matches, int indel_penalty, int max_mismatches_allowed,
			bool dibasep, bool cmetp, bool plusp) {
#ifdef DEBUG2
  int i;
  char gbuffer[MAX_QUERYLENGTH+1];
#endif
  Stage3_T hit;
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
					    query,query_compress,genome_blocks,snp_blocks,
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
					      query,query_compress,genome_blocks,snp_blocks,
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
	  debug2(printf("successful insertion with %d mismatches and indel_pos at %d\n",sum,indel_pos));
	  if (plusp == true) {
	    if ((hit = Stage3_new_insertion(&(*found_score),indels,indel_pos,
					    /*nmismatches1*/lefti,/*nmismatches2*/righti,
					    /*ncolordiffs1: colordiffs_left[lefti]*/0,
					    /*ncolordiffs2: colordiffs_right[righti]*/0,
					    /*left*/left+indels,/*genomiclength*/querylength-indels,
					    query_compress,genome_blocks,snp_blocks,
					    querylength,/*plusp*/true,query,
					    ptr->chrnum,ptr->chroffset,indel_penalty,
					    dibasep,cmetp)) == NULL) {
	      return hits;
	    } else {
	      return List_push(hits,(void *) hit);
	    }
	  } else {
	    if ((hit = Stage3_new_insertion(&(*found_score),indels,/*indel_pos*/querylength-indel_pos-indels,
					    /*nmismatches1*/righti,/*nmismatches2*/lefti,
					    /*ncolordiffs1: colordiffs_right[righti]*/0,
					    /*ncolordiffs2: colordiffs_left[lefti]*/0,
					    /*left*/left+indels,/*genomiclength*/querylength-indels,
					    query_compress,genome_blocks,snp_blocks,
					    querylength,/*plusp*/false,query,
					    ptr->chrnum,ptr->chroffset,indel_penalty,
					    dibasep,cmetp)) == NULL) {
	      return hits;
	    } else {
	      return List_push(hits,(void *) hit);
	    }
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
solve_middle_deletion (bool *foundp, int *found_score, List_T hits, Segment_T ptr, int indels,
		       Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		       char *query, char *queryptr, int querylength,
		       int min_indel_end_matches, int indel_penalty, int max_mismatches_allowed,
		       bool dibasep, bool cmetp, bool plusp) {
#ifdef DEBUG2
  int i;
  char gbuffer[MAX_QUERYLENGTH+1];
#endif
  Stage3_T hit;
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
					    query,query_compress,genome_blocks,snp_blocks,
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
					      query,query_compress,genome_blocks,snp_blocks,
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
	  debug2(printf("successful deletion with %d mismatches and indel_pos at %d\n",
			sum,indel_pos));
	  if (plusp == true) {
	    if ((hit = Stage3_new_deletion(&(*found_score),-indels,indel_pos,
					   /*nmismatches1*/lefti,/*nmismatches2*/righti,
					   /*ncolordiffs1: colordiffs_left[lefti]*/0,
					   /*ncolordiffs2: colordiffs_right[righti]*/0,
					   left,/*genomiclength*/querylength-indels,
					   query_compress,genome_blocks,snp_blocks,
					   querylength,/*plusp*/true,query,
					   ptr->chrnum,ptr->chroffset,indel_penalty,
					   dibasep,cmetp)) == NULL) {
	      return hits;
	    } else {
	      return List_push(hits,(void *) hit);
	    }
	  } else {
	    if ((hit = Stage3_new_deletion(&(*found_score),-indels,/*indel_pos*/querylength-indel_pos,
					   /*nmismatches1*/righti,/*nmismatches2*/lefti,
					   /*ncolordiffs1: colordiffs_right[righti]*/0,
					   /*ncolordiffs2: colordiffs_left[lefti]*/0,
					   left,/*genomiclength*/querylength-indels,
					   query_compress,genome_blocks,snp_blocks,
					   querylength,/*plusp*/false,query,
					   ptr->chrnum,ptr->chroffset,indel_penalty,
					   dibasep,cmetp)) == NULL) {
	      return hits;
	    } else {
	      return List_push(hits,(void *) hit);
	    }
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
		    struct Segment_T *plus_segments, struct Segment_T *minus_segments,
		    int plus_nsegments, int minus_nsegments,
		    UINT4 *genome_blocks, UINT4 *snp_blocks,
		    char *queryuc_ptr, char *queryrc, Floors_T floors, int querylength, int query_lastpos,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    int max_middle_insertions, int max_middle_deletions, int min_indel_end_matches,
		    int indel_penalty, int max_mismatches_allowed, bool dibasep, bool cmetp) {
  int indels, floor, pos, prev, middle;
  int *floors_from_neg3, *floors_to_pos3;
  int nfound = 0;
  Segment_T segmenti, segmentj;
  bool foundp;

  debug(printf("*** find_middle_indels with max_mismatches_allowed %d ***\n",
	       max_mismatches_allowed));

  if (plus_nsegments > 1) {
    floors_from_neg3 = floors->scorefrom[-3];
    floors_to_pos3 = floors->scoreto[query_lastpos+3];

    for (segmenti = plus_segments; segmenti < &(plus_segments[plus_nsegments]); segmenti++) {
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
					  /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
					  /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,
					  querylength,min_indel_end_matches,
					  indel_penalty,max_mismatches_allowed,
					  dibasep,cmetp,/*plusp*/true);
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
					 /*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
					 /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,
					 querylength,
					 min_indel_end_matches,indel_penalty,max_mismatches_allowed,
					 dibasep,cmetp,/*plusp*/true);
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
    floors_from_neg3 = floors->scorefrom[-3];
    floors_to_pos3 = floors->scoreto[query_lastpos+3];

    for (segmenti = minus_segments; segmenti < &(minus_segments[minus_nsegments]); segmenti++) {
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
					 /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
					 /*query*/queryuc_ptr,/*queryptr*/queryrc,
					 querylength,
					 min_indel_end_matches,indel_penalty,max_mismatches_allowed,
					 dibasep,cmetp,/*plusp*/false);
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
					  /*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
					  /*query*/queryuc_ptr,/*queryptr*/queryrc,
					  querylength,
					  min_indel_end_matches,indel_penalty,max_mismatches_allowed,
					  dibasep,cmetp,/*plusp*/false);
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

  return hits;
}


/************************************************************************/

static bool
compute_end_indels_right (int *indels, int *ncolordiffs, int *indel_pos, int length1, int querylength,
			  Genomicpos_T left, char *query, char *queryptr,
			  Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			  int min_indel_end_matches, int max_end_insertions, int max_end_deletions,
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
						     query,query_compress,genome_blocks,snp_blocks,
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
						     query,query_compress,genome_blocks,snp_blocks,
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
						     query,query_compress,genome_blocks,snp_blocks,
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
						     query,query_compress,genome_blocks,snp_blocks,
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
			 Genomicpos_T left, char *query, char *queryptr, Compress_T query_compress,
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
						     query,query_compress,genome_blocks,snp_blocks,
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
						   query,query_compress,genome_blocks,snp_blocks,
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
						     query,query_compress,genome_blocks,snp_blocks,
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
						   query,query_compress,genome_blocks,snp_blocks,
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
			Chrnum_T chrnum, Genomicpos_T chroffset, UINT4 *genome_blocks, UINT4 *snp_blocks,
			char *query, char *queryptr, int querylength,
			Compress_T query_compress /* expecting fwd */,
			int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
			int indel_penalty, int max_mismatches, bool dibasep, bool cmetp) {
#ifdef DEBUG2E
  char gbuffer[MAX_QUERYLENGTH];
#endif
  int i;
  Stage3_T hit;
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
					     query,query_compress,genome_blocks,snp_blocks,
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
				left,query,queryptr,query_compress,genome_blocks,snp_blocks,
				min_indel_end_matches,max_end_insertions,max_end_deletions,
#if 0
				/*max_mismatches_allowed*/max_mismatches-nmismatches_long,
#endif
				dibasep,cmetp,/*plusp*/true) == true) {
      debug2e(printf("Got indel_pos %d.\n",indel_pos));

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

	if ((hit = Stage3_new_insertion(&(*found_score),indels,/*indel_pos*/indel_pos-indels,
					/*nmismatches1*/0,/*nmismatches2*/nmismatches,
					/*ncolordiffs1*/ncolordiffs_short,/*ncolordiffs2*/ncolordiffs_long,
					left+indels,/*genomiclength*/querylength-indels,
					query_compress,genome_blocks,snp_blocks,
					querylength,/*plusp*/true,query,
					chrnum,chroffset,indel_penalty,dibasep,cmetp)) != NULL) {
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful insertion at %d with %d long + %d short mismatches\n",
			 indel_pos,nmismatches,0));
	}
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

	if ((hit = Stage3_new_deletion(&(*found_score),-indels,indel_pos,
				       /*nmismatches1*/0,/*nmismatches2*/nmismatches,
				       /*ncolordiffs1*/ncolordiffs_short,/*ncolordiffs2*/ncolordiffs_long,
				       left+indels,/*genomiclength*/querylength-indels,
				       query_compress,genome_blocks,snp_blocks,
				       querylength,/*plusp*/true,query,
				       chrnum,chroffset,indel_penalty,dibasep,cmetp)) != NULL) {
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful deletion at %d with %d long + %d short mismatches\n",
			 indel_pos,nmismatches,0));
	}
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
		       Chrnum_T chrnum, Genomicpos_T chroffset, UINT4 *genome_blocks, UINT4 *snp_blocks,
		       char *query, char *queryptr, int querylength,
		       Compress_T query_compress /* expecting fwd */,
		       int max_end_insertions, int max_end_deletions,
		       int min_indel_end_matches, int indel_penalty, int max_mismatches,
		       bool dibasep, bool cmetp) {
#ifdef DEBUG2E
  char gbuffer[MAX_QUERYLENGTH+1];
#endif
  int i;
  Stage3_T hit;
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
					    query,query_compress,genome_blocks,snp_blocks,
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
				 left,query,queryptr,query_compress,genome_blocks,snp_blocks,
				 min_indel_end_matches,max_end_insertions,max_end_deletions,
#if 0
				 /*max_mismatches_allowed*/max_mismatches-nmismatches_long,
#endif
				 dibasep,cmetp,/*plusp*/true) == true) {
      debug2e(printf("Got indel_pos %d\n",indel_pos));

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

	if ((hit = Stage3_new_insertion(&(*found_score),indels,indel_pos,
					/*nmismatches1*/nmismatches,/*nmismatches2*/0,
					/*ncolordiffs1*/ncolordiffs_long,
					/*ncolordiffs2*/ncolordiffs_short,
					left,/*genomiclength*/querylength-indels,
					query_compress,genome_blocks,snp_blocks,
					querylength,/*plusp*/true,query,chrnum,
					chroffset,indel_penalty,dibasep,cmetp)) != NULL) {
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful insertion at %d with %d long + %d short mismatches\n",
			 indel_pos,nmismatches,0));
	}
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

	if ((hit = Stage3_new_deletion(&(*found_score),-indels,indel_pos,
				       /*nmismatches1*/nmismatches,/*nmismatches2*/0,
				       /*ncolordiffs1*/ncolordiffs_long,
				       /*ncolordiffs2*/ncolordiffs_short,
				       left,/*genomiclength*/querylength-indels,
				       query_compress,genome_blocks,snp_blocks,
				       querylength,/*plusp*/true,query,chrnum,
				       chroffset,indel_penalty,dibasep,cmetp)) != NULL) {
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful deletion at %d with %d long + %d short mismatches\n",
			 indel_pos,nmismatches,0));
	}
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
			 Chrnum_T chrnum, Genomicpos_T chroffset, UINT4 *genome_blocks, UINT4 *snp_blocks,
			 char *query, char *queryptr, int querylength,
			 Compress_T query_compress /* expecting rev */,
			 int max_end_insertions, int max_end_deletions,
			 int min_indel_end_matches, int indel_penalty, int max_mismatches,
			 bool dibasep, bool cmetp) {
#ifdef DEBUG2E
  char gbuffer[MAX_QUERYLENGTH+1];
#endif
  int i;
  Stage3_T hit;
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
					    query,query_compress,genome_blocks,snp_blocks,
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
				 left,query,queryptr,query_compress,genome_blocks,snp_blocks,
				 min_indel_end_matches,max_end_insertions,max_end_deletions,
#if 0
				 /*max_mismatches_allowed*/max_mismatches-nmismatches_long,
#endif
				 dibasep,cmetp,/*plusp*/false) == true) {
      debug2e(printf("Got indel_pos %d\n",indel_pos));

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

	if ((hit = Stage3_new_insertion(&(*found_score),indels,/*indel_pos*/querylength-indel_pos-indels,
					/*nmismatches1*/0,/*nmismatches2*/nmismatches,
					/*ncolordiffs1*/ncolordiffs_short,
					/*ncolordiffs2*/ncolordiffs_long,
					left,/*genomiclength*/querylength-indels,
					query_compress,genome_blocks,snp_blocks,
					querylength,/*plusp*/false,query,chrnum,
					chroffset,indel_penalty,dibasep,cmetp)) != NULL) {
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful insertion at %d with %d long + %d short mismatches\n",
			 indel_pos,nmismatches,0));
	}
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

	if ((hit = Stage3_new_deletion(&(*found_score),-indels,querylength-indel_pos,
				       /*nmismatches1*/0,/*nmismatches2*/nmismatches,
				       /*ncolordiffs1*/ncolordiffs_short,
				       /*ncolordiffs2*/ncolordiffs_long,
				       left,/*genomiclength*/querylength-indels,
				       query_compress,genome_blocks,snp_blocks,
				       querylength,/*plusp*/false,query,chrnum,
				       chroffset,indel_penalty,dibasep,cmetp)) != NULL) {
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful deletion at %d with %d long + %d short mismatches\n",
			 indel_pos,nmismatches,0));
	}
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
			Chrnum_T chrnum, Genomicpos_T chroffset, UINT4 *genome_blocks, UINT4 *snp_blocks,
			char *query, char *queryptr, int querylength,
			Compress_T query_compress /* expecting rev */,
			int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
			int indel_penalty,int max_mismatches, bool dibasep, bool cmetp) {
#ifdef DEBUG2E
  char gbuffer[MAX_QUERYLENGTH];
#endif
  int i;
  Stage3_T hit;
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
					     query,query_compress,genome_blocks,snp_blocks,
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
				left,query,queryptr,query_compress,genome_blocks,snp_blocks,
				min_indel_end_matches,max_end_insertions,max_end_deletions,
#if 0
				/*max_mismatches_allowed*/max_mismatches-nmismatches_long,
#endif
				dibasep,cmetp,/*plusp*/false) == true) {
      debug2e(printf("Got indel_pos %d\n",indel_pos));

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

	if ((hit = Stage3_new_insertion(&(*found_score),indels,querylength-indel_pos,
					/*nmismatches1*/nmismatches,/*nmismatches2*/0,
					/*ncolordiffs1*/ncolordiffs_long,
					/*ncolordiffs2*/ncolordiffs_short,
					left+indels,/*genomiclength*/querylength-indels,
					query_compress,genome_blocks,snp_blocks,
					querylength,/*plusp*/false,query,
					chrnum,chroffset,indel_penalty,dibasep,cmetp)) != NULL) {
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful insertion at %d with %d long + %d short mismatches\n",
			 indel_pos,nmismatches,0));
	}
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

	if ((hit = Stage3_new_deletion(&(*found_score),-indels,querylength-indel_pos,
				       /*nmismatches1*/nmismatches,/*nmismatches2*/0,
				       /*ncolordiffs1*/ncolordiffs_long,
				       /*ncolordiffs2*/ncolordiffs_short,
				       left+indels,/*genomiclength*/querylength-indels,
				       query_compress,genome_blocks,snp_blocks,
				       querylength,/*plusp*/false,query,
				       chrnum,chroffset,indel_penalty,dibasep,cmetp)) != NULL) {
	  hits = List_push(hits,(void *) hit);
	  debug2e(printf("successful deletion at %d with %d long + %d short mismatches\n",
			 indel_pos,nmismatches,0));
	}
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
		 struct Segment_T *plus_segments, struct Segment_T *minus_segments,
		 int plus_nsegments, int minus_nsegments, UINT4 *genome_blocks, UINT4 *snp_blocks,
		 char *queryuc_ptr, char *queryrc, int querylength, int firstbound, int lastbound,
		 Compress_T query_compress_fwd, Compress_T query_compress_rev,
		 int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		 int indel_penalty, int max_mismatches_allowed, bool dibasep, bool cmetp) {
  Segment_T ptr;

  debug(printf("*** find_end_indels with max_mismatches_allowed %d ***\n",
	       max_mismatches_allowed));

  for (ptr = plus_segments; ptr < &(plus_segments[plus_nsegments]); ptr++) {

    if (ptr->floor_xfirst <= max_mismatches_allowed) {
      /* First indel, plus */
      debug2e(printf("floor_xfirst %d <= mismatches allowed %d\n",ptr->floor_xfirst,max_mismatches_allowed));
      hits = solve_first_indel_plus(&(*found_score),hits,ptr->diagonal,firstbound,
				    ptr->chrnum,ptr->chroffset,genome_blocks,snp_blocks,
				    /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,
				    querylength,/*query_compress*/query_compress_fwd,
				    max_end_insertions,max_end_deletions,min_indel_end_matches,
				    indel_penalty,max_mismatches_allowed,dibasep,cmetp);
    }

    if (ptr->floor_xlast <= max_mismatches_allowed) {
      /* Last indel, plus */
      debug2e(printf("floor_xlast %d <= mismatches allowed %d\n",ptr->floor_xlast,max_mismatches_allowed));
      hits = solve_last_indel_plus(&(*found_score),hits,ptr->diagonal,lastbound,
				   ptr->chrnum,ptr->chroffset,genome_blocks,snp_blocks,
				   /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,
				   querylength,/*query_compress*/query_compress_fwd,
				   max_end_insertions,max_end_deletions,min_indel_end_matches,
				   indel_penalty,max_mismatches_allowed,dibasep,cmetp);
    }
  }

  for (ptr = minus_segments; ptr < &(minus_segments[minus_nsegments]); ptr++) {

    if (ptr->floor_xfirst <= max_mismatches_allowed) {
      /* First indel, minus */
      debug2e(printf("floor_xfirst %d <= mismatches allowed %d\n",ptr->floor_xfirst,max_mismatches_allowed));
      hits = solve_first_indel_minus(&(*found_score),hits,ptr->diagonal,lastbound,
				     ptr->chrnum,ptr->chroffset,genome_blocks,snp_blocks,
				     /*query*/queryuc_ptr,/*queryptr*/queryrc,
				     querylength,/*query_compress*/query_compress_rev,
				     max_end_insertions,max_end_deletions,min_indel_end_matches,
				     indel_penalty,max_mismatches_allowed,dibasep,cmetp);
    }

    if (ptr->floor_xlast <= max_mismatches_allowed) {
      /* Last indel, minus */
      debug2e(printf("floor_xlast %d <= mismatches allowed %d\n",ptr->floor_xlast,max_mismatches_allowed));
      hits = solve_last_indel_minus(&(*found_score),hits,ptr->diagonal,firstbound,
				    ptr->chrnum,ptr->chroffset,genome_blocks,snp_blocks,
				    /*query*/queryuc_ptr,/*queryptr*/queryrc,
				    querylength,/*query_compress*/query_compress_rev,
				    max_end_insertions,max_end_deletions,min_indel_end_matches,
				    indel_penalty,max_mismatches_allowed,dibasep,cmetp);
    }
  }

  return hits;
}



/************************************************************************
 *   Splicing
 ************************************************************************/


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
  } else {
    return (spliceprob > 0.70);
  }
}

/* Generally need 25 positions to be unique in the human genome */
#define MIN_SPLICE_SUPPORT_DISTANT 25

/* Do not compare against true or false */
/* Moderate criterion */
static int
sufficient_splice_prob_distant (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support < MIN_SPLICE_SUPPORT_DISTANT) {
    return 0;
  } else if (support < 30) {
    return (spliceprob > 0.95);
  } else if (support < 35) {
    return (spliceprob > 0.90);
  } else if (support < 40) {
    return (spliceprob > 0.85);
  } else {
    return (spliceprob > 0.70);
  }
}

/* Do not compare against true or false */

#ifdef HALFINTRON
/* Strictest criterion */
static int
sufficient_splice_prob_halfintron (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support < 20) {
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
#endif



#if 0
static void
find_segmentm_span (Segment_T segmentm, int max_mismatches_allowed,
		    char *query, int querylength, Compress_T query_compress,
		    UINT4 *genome_blocks, UINT4 *snp_blocks, Genomicpos_T left,
		    bool dibasep, bool cmetp, bool plusp) {
  int mismatch_positions[MAX_QUERYLENGTH];
  int colordiffs[MAX_QUERYLENGTH];
  int nmismatches, i;
  int leftspan, rightspan, bestspan;

  /* Find all mismatches */
  nmismatches = Genome_mismatches_left(mismatch_positions,colordiffs,/*max_mismatches*/querylength,
				       query,query_compress,genome_blocks,snp_blocks,
				       left,/*pos5*/0,/*pos3*/querylength,dibasep,cmetp,plusp);

  if (nmismatches < max_mismatches_allowed) {
    segmentm->leftspan = 0;
    segmentm->rightspan = querylength;
  } else {
    segmentm->leftspan = 0;
    bestspan = segmentm->rightspan = mismatch_positions[max_mismatches_allowed] + /*slop*/ 1;
    for (i = 0; i < nmismatches - max_mismatches_allowed; i++) {
      leftspan = mismatch_positions[i];
      rightspan = mismatch_positions[i + max_mismatches_allowed + 1] + /*slop*/ 1;
      if (rightspan - leftspan > bestspan) {
	segmentm->leftspan = leftspan;
	segmentm->rightspan = rightspan;
	bestspan = rightspan - leftspan;
      } else if (rightspan - leftspan == bestspan) {
	segmentm->rightspan = rightspan;
      }
    }
  }
  return;
}
#endif


/* Note: knowni holds joffset + j + 1, so 0 represents no known site
   and values greater than 0 represent a known site.  Need to subtract
   1 to obtain joffset + j. */

static List_T
solve_singlesplice (int *found_score, int *nhits, List_T hits, Segment_T segmenti, Segment_T segmentj,
		    UINT4 *genome_blocks, UINT4 *snp_blocks,
		    char *query, int querylength, Compress_T query_compress,
		    int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
		    int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
		    int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
		    int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
		    int segmenti_donor_nknown, int segmentj_acceptor_nknown,
		    int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
		    bool novelsplicingp, int splicing_penalty,
		    int min_localsplicing_end_matches, int max_mismatches_allowed,
		    bool first_read_p, bool dibasep, bool cmetp, bool plusp) {
  Substring_T donor, acceptor;
  Genomicpos_T segmenti_left, segmentj_left;
  int best_splice_pos, splice_pos_start, splice_pos_end, splice_pos, i, j;
  int donor_positions_alloc[MAX_QUERYLENGTH+1], acceptor_positions_alloc[MAX_QUERYLENGTH+1];
  int donor_knowni_alloc[MAX_QUERYLENGTH+1], acceptor_knowni_alloc[MAX_QUERYLENGTH+1];

  int best_segmenti_nmismatches, best_segmentj_nmismatches, segmenti_nmismatches, segmentj_nmismatches, ncolordiffs;
  int donor_support, acceptor_support;
  int best_donor_knowni, best_acceptor_knowni;
  double best_prob, prob, best_donor_prob, best_acceptor_prob, probi, probj;
  bool orig_plusp, sensep;
  int sensedir;

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;


  segmenti_left = segmenti->diagonal - querylength;
  segmentj_left = segmentj->diagonal - querylength;
  debug4p(printf("solve_singlesplice: Getting genome at lefti %u and leftj %u\n",
		 segmenti_left,segmentj_left));

#if 0
  int sum, lefti, righti;
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
#else
  /* splice_pos_start = min_localsplicing_end_matches; */
  /* splice_pos_end = querylength - min_localsplicing_end_matches; */
  splice_pos_start = 2;
  splice_pos_end = querylength - 2;
#endif


  if (splice_pos_start <= splice_pos_end) {
    /* Originally from plus strand.  No complement.  */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(donor_positions_alloc,donor_knowni_alloc,
					     segmenti_donor_knownpos,segmenti_donor_knowni,
					     genome_blocks,segmenti_left,
					     splice_pos_start,splice_pos_end);
      donori_positions = donor_positions_alloc;
      donori_knowni = donor_knowni_alloc;
    } else {
      donori_nsites = segmenti_donor_nknown;
      donori_positions = segmenti_donor_knownpos;
      donori_knowni = segmenti_donor_knowni;
    }

#ifdef DEBUG4P
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
						   segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						   genome_blocks,segmentj_left,
						   splice_pos_start,splice_pos_end);
      acceptorj_positions = acceptor_positions_alloc;
      acceptorj_knowni = acceptor_knowni_alloc;
    } else {
      acceptorj_nsites = segmentj_acceptor_nknown;
      acceptorj_positions = segmentj_acceptor_knownpos;
      acceptorj_knowni = segmentj_acceptor_knowni;
    }

#ifdef DEBUG4P
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_prob = 0.0;
    orig_plusp = true;

    i = j = 0;
    while (i < donori_nsites && j < acceptorj_nsites) {
      if ((splice_pos = donori_positions[i]) < acceptorj_positions[j]) {
	i++;
      } else if (splice_pos > acceptorj_positions[j]) {
	j++;
      } else {
	if (donori_knowni[i] >= 0) {
	  probi = 1.0;
	} else {
	  probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,genome_blocks);
	}

	if (acceptorj_knowni[j] >= 0) {
	  probj = 1.0;
	} else {
	  probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,genome_blocks);
	}

	debug4p(
		if (plusp == true) {
		  printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		} else {
		  printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		});
	if ((prob = probi + probj) > best_prob) {
	  segmenti_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								   /*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   dibasep,cmetp,plusp);
	  segmentj_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								   /*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   dibasep,cmetp,plusp);
	  if (segmenti_nmismatches + segmentj_nmismatches <= max_mismatches_allowed) {
	    donor_support = splice_pos;
	    acceptor_support = querylength - splice_pos;
	    if (sufficient_splice_prob_local(donor_support,segmenti_nmismatches,probi) &&
		sufficient_splice_prob_local(acceptor_support,segmentj_nmismatches,probj)) {
	      /* Success */
	      best_donor_knowni = donori_knowni[i];
	      best_acceptor_knowni = acceptorj_knowni[j];
	      best_donor_prob = probi;
	      best_acceptor_prob = probj;
	      best_prob = prob;
	      best_splice_pos = splice_pos;
	      best_segmenti_nmismatches = segmenti_nmismatches;
	      best_segmentj_nmismatches = segmentj_nmismatches;
	    }
	  }
	}
	i++;
	j++;
      }
    }


    /* Originally from minus strand.  Complement. */
    /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
							   segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							   genome_blocks,segmenti_left,
							   splice_pos_start,splice_pos_end);
      antiacceptori_positions = acceptor_positions_alloc;
      antiacceptori_knowni = acceptor_knowni_alloc;
    } else {
      antiacceptori_nsites = segmenti_antiacceptor_nknown;
      antiacceptori_positions = segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG4P
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(donor_positions_alloc,donor_knowni_alloc,
						     segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						     genome_blocks,segmentj_left,
						     splice_pos_start,splice_pos_end);
      antidonorj_positions = donor_positions_alloc;
      antidonorj_knowni = donor_knowni_alloc;
    } else {
      antidonorj_nsites = segmentj_antidonor_nknown;
      antidonorj_positions = segmentj_antidonor_knownpos;
      antidonorj_knowni = segmentj_antidonor_knowni;
    }

#ifdef DEBUG4P
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    i = j = 0;
    while (i < antiacceptori_nsites && j < antidonorj_nsites) {
      if ((splice_pos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	i++;
      } else if (splice_pos > antidonorj_positions[j]) {
	j++;
      } else {
	if (antiacceptori_knowni[i] >= 0) {
	  probi = 1.0;
	} else {
	  probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,genome_blocks);
	}

	if (antidonorj_knowni[j] >= 0) {
	  probj = 1.0;
	} else {
	  probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,genome_blocks);
	}

	debug4p(
		if (plusp == true) {
		  printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		} else {
		  printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		});
	if ((prob = probi + probj) > best_prob) {
	  segmenti_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								   /*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								   dibasep,cmetp,plusp);
	  segmentj_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								   /*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								   dibasep,cmetp,plusp);
	  if (segmenti_nmismatches + segmentj_nmismatches <= max_mismatches_allowed) {
	    acceptor_support = splice_pos;
	    donor_support = querylength - splice_pos;
	    if (sufficient_splice_prob_local(acceptor_support,segmenti_nmismatches,probi) &&
		sufficient_splice_prob_local(donor_support,segmentj_nmismatches,probj)) {
	      /* Success */
	      best_donor_knowni = antidonorj_knowni[j];
	      best_acceptor_knowni = antiacceptori_knowni[i];
	      best_donor_prob = probj;
	      best_acceptor_prob = probi;
	      best_prob = prob;
	      best_splice_pos = splice_pos;
	      best_segmentj_nmismatches = segmentj_nmismatches;
	      best_segmenti_nmismatches = segmenti_nmismatches;
	      orig_plusp = false;
	    }
	  }
	}
	i++;
	j++;
      }
    }

    if (best_prob > 0.0) {
      debug4p(printf("best_prob = %f at splice_pos %d (%d,%d)\n",
		     best_prob,best_splice_pos,best_donor_knowni,best_acceptor_knowni));
      if (orig_plusp == true) {
	/* Originally from plus strand.  No complement. */
	sensep = (plusp == true) ? true : false;
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	donor = Substring_new_donor(best_donor_knowni,/*joffset*/0,best_splice_pos,best_segmenti_nmismatches,ncolordiffs,
				    best_donor_prob,/*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
				    querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

	acceptor = Substring_new_acceptor(best_acceptor_knowni,/*joffset*/0,best_splice_pos,best_segmentj_nmismatches,ncolordiffs,
					  best_acceptor_prob,/*left*/segmentj_left,query_compress,genome_blocks,snp_blocks,
					  querylength,plusp,sensep,query,segmentj->chrnum,segmentj->chroffset,dibasep,cmetp);
	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  *nhits += 1;
	  debug4p(printf("solve_singlesplice success\n"));
	  return List_push(hits,(void *) Stage3_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
							   donor,acceptor,/*distance*/segmentj_left - segmenti_left,
							   /*shortdistancep*/true,splicing_penalty,querylength,
							   /*ambi_left*/NULL,/*ambi_right*/NULL,
							   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
							   /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir));
	}

      } else {
	/* Originally from minus strand.  Complement. */
	sensep = (plusp == true) ? false : true;
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	donor = Substring_new_donor(best_donor_knowni,/*joffset*/0,best_splice_pos,best_segmentj_nmismatches,ncolordiffs,
				    best_donor_prob,/*left*/segmentj_left,query_compress,genome_blocks,snp_blocks,
				    querylength,plusp,sensep,query,segmentj->chrnum,segmentj->chroffset,dibasep,cmetp);

	acceptor = Substring_new_acceptor(best_acceptor_knowni,/*joffset*/0,best_splice_pos,best_segmenti_nmismatches,ncolordiffs,
					  best_acceptor_prob,/*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
					  querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);
	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  *nhits += 1;
	  debug4p(printf("solve_singlesplice success\n"));
	  return List_push(hits,(void *) Stage3_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
							   donor,acceptor,/*distance*/segmentj_left - segmenti_left,
							   /*shortdistancep*/true,splicing_penalty,querylength,
							   /*ambi_left*/NULL,/*ambi_right*/NULL,
							   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
							   /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir));
	}
      }
    }
  }

  debug4p(printf("solve_singlesplice fail\n"));
  return hits;
}


static List_T
solve_doublesplice (int *found_score, int *nhits, List_T hits,
		    Segment_T segmenti, Segment_T segmentm, Segment_T segmentj,
		    UINT4 *genome_blocks, UINT4 *snp_blocks,
		    char *query, int querylength, Compress_T query_compress,
		    int *segmenti_donor_knownpos, int *segmentm_acceptor_knownpos, int *segmentm_donor_knownpos, int *segmentj_acceptor_knownpos,
		    int *segmentj_antidonor_knownpos, int *segmentm_antiacceptor_knownpos, int *segmentm_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
		    int *segmenti_donor_knowni, int *segmentm_acceptor_knowni, int *segmentm_donor_knowni, int *segmentj_acceptor_knowni,
		    int *segmentj_antidonor_knowni, int *segmentm_antiacceptor_knowni, int *segmentm_antidonor_knowni, int *segmenti_antiacceptor_knowni,
		    int segmenti_donor_nknown, int segmentm_acceptor_nknown, int segmentm_donor_nknown, int segmentj_acceptor_nknown,
		    int segmentj_antidonor_nknown, int segmentm_antiacceptor_nknown, int segmentm_antidonor_nknown, int segmenti_antiacceptor_nknown,
		    bool novelsplicingp, int splicing_penalty, int min_localsplicing_end_matches, int max_mismatches_allowed,
		    bool first_read_p, bool dibasep, bool cmetp, bool plusp) {
  Substring_T donor, shortexon, acceptor;
  Genomicpos_T segmenti_left, segmentm_left, segmentj_left;
  int best_splice_pos_1, best_splice_pos_2, splice_pos_start, splice_pos_end, splice_pos_1, splice_pos_2;
  int i, a, b, j;
  int donor1_positions_alloc[MAX_QUERYLENGTH+1], acceptor1_positions_alloc[MAX_QUERYLENGTH+1],
    donor2_positions_alloc[MAX_QUERYLENGTH+1], acceptor2_positions_alloc[MAX_QUERYLENGTH+1];
  int donor1_knowni_alloc[MAX_QUERYLENGTH+1], acceptor1_knowni_alloc[MAX_QUERYLENGTH+1],
    donor2_knowni_alloc[MAX_QUERYLENGTH+1], acceptor2_knowni_alloc[MAX_QUERYLENGTH+1];

  int best_segmenti_nmismatches, best_segmentm_nmismatches, best_segmentj_nmismatches,
    segmenti_nmismatches, segmentm_nmismatches, segmentj_nmismatches, ncolordiffs;
  int donor_support, acceptor_support, middle_support;
  int best_donor1_knowni, best_acceptor1_knowni, best_donor2_knowni, best_acceptor2_knowni;
  double best_prob, best_donor1_prob, best_acceptor1_prob, best_donor2_prob, best_acceptor2_prob,
    prob, probi, proba, probb, probj;
  bool orig_plusp, sensep, matchp;
  int sensedir;

  int donori_nsites, acceptora_nsites, donorb_nsites, acceptorj_nsites,
    antiacceptori_nsites, antidonora_nsites, antiacceptorb_nsites, antidonorj_nsites;
  int *donori_positions, *acceptora_positions, *donorb_positions, *acceptorj_positions,
    *antiacceptori_positions, *antidonora_positions, *antiacceptorb_positions, *antidonorj_positions;
  int *donori_knowni, *acceptora_knowni, *donorb_knowni, *acceptorj_knowni,
    *antiacceptori_knowni, *antidonora_knowni, *antiacceptorb_knowni, *antidonorj_knowni;

  segmenti_left = segmenti->diagonal - querylength;
  segmentm_left = segmentm->diagonal - querylength;
  segmentj_left = segmentj->diagonal - querylength;
  debug4d(printf("solve_doublesplice: Getting genome at lefti %u, leftm %u, and leftj %u\n",
		 segmenti_left,segmentm_left,segmentj_left));

  splice_pos_start = 2;
  splice_pos_end = querylength - 2;

  if (splice_pos_start <= splice_pos_end) {
    /* Originally from plus strand.  No complement. */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */

    /* Segment i */
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(donor1_positions_alloc,donor1_knowni_alloc,
					     segmenti_donor_knownpos,segmenti_donor_knowni,
					     genome_blocks,segmenti_left,
					     splice_pos_start,splice_pos_end);
      donori_positions = donor1_positions_alloc;
      donori_knowni = donor1_knowni_alloc;
    } else {
      donori_nsites = segmenti_donor_nknown;
      donori_positions = segmenti_donor_knownpos;
      donori_knowni = segmenti_donor_knowni;
    }

#ifdef DEBUG4D
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m1 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptora_nsites = Genome_acceptor_positions(acceptor1_positions_alloc,acceptor1_knowni_alloc,
						   segmentm_acceptor_knownpos,segmentm_acceptor_knowni,
						   genome_blocks,segmentm_left,
						   splice_pos_start,splice_pos_end);
      acceptora_positions = acceptor1_positions_alloc;
      acceptora_knowni = acceptor1_knowni_alloc;
    } else {
      acceptora_nsites = segmentm_acceptor_nknown;
      acceptora_positions = segmentm_acceptor_knownpos;
      acceptora_knowni = segmentm_acceptor_knowni;
    }

#ifdef DEBUG4D
    printf("Found %d acceptora sites:",acceptora_nsites);
    for (i = 0; i < acceptora_nsites; i++) {
      printf(" %d",acceptora_positions[i]);
      if (acceptora_knowni[i] >= 0) {
	printf(" (%d)",acceptora_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m2 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donorb_nsites = Genome_donor_positions(donor2_positions_alloc,donor2_knowni_alloc,
					     segmentm_donor_knownpos,segmentm_donor_knowni,
					     genome_blocks,segmentm_left,
					     splice_pos_start,splice_pos_end);
      donorb_positions = donor2_positions_alloc;
      donorb_knowni = donor2_knowni_alloc;
    } else {
      donorb_nsites = segmentm_donor_nknown;
      donorb_positions = segmentm_donor_knownpos;
      donorb_knowni = segmentm_donor_knowni;
    }

#ifdef DEBUG4D
    printf("Found %d donorb sites:",donorb_nsites);
    for (i = 0; i < donorb_nsites; i++) {
      printf(" %d",donorb_positions[i]);
      if (donorb_knowni[i] >= 0) {
	printf(" (%d)",donorb_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment j */
    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(acceptor2_positions_alloc,acceptor2_knowni_alloc,
						   segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						   genome_blocks,segmentj_left,
						   splice_pos_start,splice_pos_end);
      acceptorj_positions = acceptor2_positions_alloc;
      acceptorj_knowni = acceptor2_knowni_alloc;
    } else {
      acceptorj_nsites = segmentj_acceptor_nknown;
      acceptorj_positions = segmentj_acceptor_knownpos;
      acceptorj_knowni = segmentj_acceptor_knowni;
    }

#ifdef DEBUG4D
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_prob = 0.0;
    orig_plusp = true;

    i = a = b = j = 0;
    while (i < donori_nsites && a < acceptora_nsites) {
      if ((splice_pos_1 = donori_positions[i]) < acceptora_positions[a]) {
	i++;
      } else if (splice_pos_1 > acceptora_positions[a]) {
	a++;
      } else {
	while (b < donorb_nsites && donorb_positions[b] <= splice_pos_1) {
	  b++;
	}
	while (j < acceptorj_nsites && acceptorj_positions[j] <= splice_pos_1) {
	  j++;
	}
	matchp = false;
	while (b < donorb_nsites && j < acceptorj_nsites && matchp == false) {
	  if ((splice_pos_2 = donorb_positions[b]) < acceptorj_positions[j]) {
	    b++;
	  } else if (splice_pos_2 > acceptorj_positions[j]) {
	    j++;
	  } else {
	    if (donori_knowni[i] >= 0) {
	      probi = 1.0;
	    } else {
	      probi = Maxent_hr_donor_prob(segmenti_left + splice_pos_1,genome_blocks);
	    }

	    if (acceptora_knowni[a] >= 0) {
	      proba = 1.0;
	    } else {
	      proba = Maxent_hr_acceptor_prob(segmentm_left + splice_pos_1,genome_blocks);
	    }

	    if (donorb_knowni[b] >= 0) {
	      probi = 1.0;
	    } else {
	      probb = Maxent_hr_donor_prob(segmentm_left + splice_pos_2,genome_blocks);
	    }

	    if (acceptorj_knowni[j] >= 0) {
	      probj = 1.0;
	    } else {
	      probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos_2,genome_blocks);
	    }

	    debug4d(
		    if (plusp == true) {
		      printf("plus sense splice_pos  %d, %d, i.donor %f, m.acceptor %f, m.donor %f, j.acceptor %f\n",
			     splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		    } else {
			printf("minus antisense splice_pos  %d %d, i.donor %f, m.acceptor %f, m.donor %f, j.acceptor %f\n",
			       splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		    });

	    if ((prob = probi + proba + probb + probj) > best_prob) {
	      segmenti_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								       /*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos_1,
								       dibasep,cmetp,plusp);
	      segmentm_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								       /*left*/segmentm_left,/*pos5*/splice_pos_1,/*pos3*/splice_pos_2,
								       dibasep,cmetp,plusp);
	      segmentj_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								       /*left*/segmentj_left,/*pos5*/splice_pos_2,/*pos3*/querylength,
								       dibasep,cmetp,plusp);
	      if (segmenti_nmismatches + segmentm_nmismatches + segmentj_nmismatches <= max_mismatches_allowed) {
		donor_support = splice_pos_1;
		middle_support = splice_pos_2 - splice_pos_1;
		acceptor_support = querylength - splice_pos_2;
		if (sufficient_splice_prob_local(donor_support,segmenti_nmismatches,probi) &&
		    sufficient_splice_prob_local(middle_support,segmentm_nmismatches,proba) &&
		    sufficient_splice_prob_local(middle_support,segmentm_nmismatches,probb) &&
		    sufficient_splice_prob_local(acceptor_support,segmentj_nmismatches,probj)) {
		  /* Success */
		  best_donor1_knowni = donori_knowni[i];
		  best_acceptor1_knowni = acceptora_knowni[a];
		  best_donor2_knowni = donorb_knowni[b];
		  best_acceptor2_knowni = acceptorj_knowni[j];
		  best_donor1_prob = probi;
		  best_acceptor1_prob = proba;
		  best_donor2_prob = probb;
		  best_acceptor2_prob = probj;
		  best_prob = prob;
		  best_splice_pos_1 = splice_pos_1;
		  best_splice_pos_2 = splice_pos_2;
		  best_segmenti_nmismatches = segmenti_nmismatches;
		  best_segmentm_nmismatches = segmentm_nmismatches;
		  best_segmentj_nmismatches = segmentj_nmismatches;
		}
	      }
	    }
	    /* b++; j++; Don't advance b or j, so next i/a can match */
	    matchp = true;
	  }
	}
	i++;
	a++;
      }
    }


    /* Originally from minus strand.  Complement. */
    /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */

    /* Segment i */
    if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(acceptor1_positions_alloc,acceptor1_knowni_alloc,
							   segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							   genome_blocks,segmenti_left,
							   splice_pos_start,splice_pos_end);
      antiacceptori_positions = acceptor1_positions_alloc;
      antiacceptori_knowni = acceptor1_knowni_alloc;
    } else {
      antiacceptori_nsites = segmenti_antiacceptor_nknown;
      antiacceptori_positions = segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG4D
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m1 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonora_nsites = Genome_antidonor_positions(donor1_positions_alloc,donor1_knowni_alloc,
						     segmentm_antidonor_knownpos,segmentm_antidonor_knowni,
						     genome_blocks,segmentm_left,
						     splice_pos_start,splice_pos_end);
      antidonora_positions = donor1_positions_alloc;
      antidonora_knowni = donor1_knowni_alloc;
    } else {
      antidonora_nsites = segmentm_antidonor_nknown;
      antidonora_positions = segmentm_antidonor_knownpos;
      antidonora_knowni = segmentm_antidonor_knowni;
    }

#ifdef DEBUG4D
    printf("Found %d antidonora sites:",antidonora_nsites);
    for (i = 0; i < antidonora_nsites; i++) {
      printf(" %d",antidonora_positions[i]);
      if (antidonora_knowni[i] >= 0) {
	printf(" (%d)",antidonora_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m2 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptorb_nsites = Genome_antiacceptor_positions(acceptor2_positions_alloc,acceptor2_knowni_alloc,
							   segmentm_antiacceptor_knownpos,segmentm_antiacceptor_knowni,
							   genome_blocks,segmentm_left,
							   splice_pos_start,splice_pos_end);
      antiacceptorb_positions = acceptor2_positions_alloc;
      antiacceptorb_knowni = acceptor2_knowni_alloc;
    } else {
      antiacceptorb_nsites = segmentm_antiacceptor_nknown;
      antiacceptorb_positions = segmentm_antiacceptor_knownpos;
      antiacceptorb_knowni = segmentm_antiacceptor_knowni;
    }

#ifdef DEBUG4D
    printf("Found %d antiacceptorb sites:",antiacceptorb_nsites);
    for (i = 0; i < antiacceptorb_nsites; i++) {
      printf(" %d",antiacceptorb_positions[i]);
      if (antiacceptorb_knowni[i] >= 0) {
	printf(" (%d)",antiacceptorb_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment j */
    if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(donor2_positions_alloc,donor2_knowni_alloc,
						     segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						     genome_blocks,segmentj_left,
						     splice_pos_start,splice_pos_end);
      antidonorj_positions = donor2_positions_alloc;
      antidonorj_knowni = donor2_knowni_alloc;
    } else {
      antidonorj_nsites = segmentj_antidonor_nknown;
      antidonorj_positions = segmentj_antidonor_knownpos;
      antidonorj_knowni = segmentj_antidonor_knowni;
    }

#ifdef DEBUG4D
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      }
    }
    printf("\n");
#endif


    i = a = b = j = 0;
    while (i < antiacceptori_nsites && a < antidonora_nsites) {
      if ((splice_pos_1 = antiacceptori_positions[i]) < antidonora_positions[a]) {
	i++;
      } else if (splice_pos_1 > antidonora_positions[a]) {
	a++;
      } else {
	while (b < antiacceptorb_nsites && antiacceptorb_positions[b] <= splice_pos_1) {
	  b++;
	}
	while (j < antidonorj_nsites && antidonorj_positions[j] <= splice_pos_1) {
	  j++;
	}
	matchp = false;
	while (b < antiacceptorb_nsites && j < antidonorj_nsites && matchp == false) {
	  if ((splice_pos_2 = antiacceptorb_positions[b]) < antidonorj_positions[j]) {
	    b++;
	  } else if (splice_pos_2 > antidonorj_positions[j]) {
	    j++;
	  } else {
	    if (antiacceptori_knowni[i] >= 0) {
	      probi = 1.0;
	    } else {
	      probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos_1,genome_blocks);
	    }
	    
	    if (antidonora_knowni[a] >= 0) {
	      proba = 1.0;
	    } else {
	      proba = Maxent_hr_antidonor_prob(segmentm_left + splice_pos_1,genome_blocks);
	    }

	    if (antiacceptorb_knowni[b] >= 0) {
	      probb = 1.0;
	    } else {
	      probb = Maxent_hr_antiacceptor_prob(segmentm_left + splice_pos_2,genome_blocks);
	    }

	    if (antidonorj_knowni[j] >= 0) {
	      probj = 1.0;
	    } else {
	      probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos_2,genome_blocks);
	    }

	    debug4d(
		    if (plusp == true) {
		      printf("plus antisense splice_pos  %d, %d, i.antiacceptor %f, m.antidonor %f, m.antiacceptor %f, j.antidonor %f\n",
			     splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		    } else {
		      printf("minus sense splice_pos  %d, %d, i.antiacceptor %f, m.antidonor %f, m.antiacceptor %f, j.antidonor %f\n",
			     splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		    });
	    if ((prob = probi + proba + probb + probj) > best_prob) {
	      segmenti_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								       /*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos_1,
								       dibasep,cmetp,plusp);
	      segmentm_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								       /*left*/segmentm_left,/*pos5*/splice_pos_1,/*pos3*/splice_pos_2,
								       dibasep,cmetp,plusp);
	      segmentj_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								       /*left*/segmentj_left,/*pos5*/splice_pos_2,/*pos3*/querylength,
								       dibasep,cmetp,plusp);
	      if (segmenti_nmismatches + segmentm_nmismatches + segmentj_nmismatches <= max_mismatches_allowed) {
		acceptor_support = splice_pos_1;
		middle_support = splice_pos_2 - splice_pos_1;
		donor_support = querylength - splice_pos_2;

		if (sufficient_splice_prob_local(acceptor_support,segmenti_nmismatches,probi) &&
		    sufficient_splice_prob_local(middle_support,segmentm_nmismatches,proba) &&
		    sufficient_splice_prob_local(middle_support,segmentm_nmismatches,probb) &&
		    sufficient_splice_prob_local(donor_support,segmentj_nmismatches,probj)) {
		  /* Success */
		  best_acceptor1_knowni = antiacceptori_knowni[i];
		  best_donor1_knowni = antidonora_knowni[a];
		  best_acceptor2_knowni = antiacceptorb_knowni[b];
		  best_donor2_knowni = antidonorj_knowni[j];
		  best_acceptor1_prob = probi;
		  best_donor1_prob = proba;
		  best_acceptor2_prob = probb;
		  best_donor2_prob = probj;
		  best_prob = prob;
		  best_splice_pos_1 = splice_pos_1;
		  best_splice_pos_2 = splice_pos_2;
		  best_segmenti_nmismatches = segmenti_nmismatches;
		  best_segmentm_nmismatches = segmentm_nmismatches;
		  best_segmentj_nmismatches = segmentj_nmismatches;
		  orig_plusp = false;
		}
	      }
	    }
	    /* b++; j++; Don't advance b or j, so next i/a can match */
	    matchp = true;
	  }
	}
	i++;
	a++;
      }
    }


    if (best_prob > 0.0) {
      debug4d(printf("best_prob = %f at splice_pos %d and %d\n",best_prob,best_splice_pos_1,best_splice_pos_2));
      if (orig_plusp == true) {
	/* Originally from plus strand.  No complement. */
	sensep = (plusp == true) ? true : false;
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	donor = Substring_new_donor(best_donor1_knowni,/*joffset*/0,best_splice_pos_1,best_segmenti_nmismatches,ncolordiffs,
				    best_donor1_prob,/*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
				    querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

	shortexon = Substring_new_shortexon(best_acceptor1_knowni,best_donor2_knowni,/*joffset*/0,
					    /*acceptor_pos*/best_splice_pos_1,/*donor_pos*/best_splice_pos_2,best_segmentm_nmismatches,
					    /*ncolordiffs*/0,/*acceptor_prob*/best_acceptor1_prob,/*donor_prob*/best_donor2_prob,
					    /*left*/segmentm_left,query_compress,genome_blocks,snp_blocks,
					    querylength,plusp,sensep,/*acceptor_ambp*/false,/*donor_ambp*/false,
					    query,segmentm->chrnum,segmentm->chroffset,dibasep,cmetp);

	acceptor = Substring_new_acceptor(best_acceptor2_knowni,/*joffset*/0,best_splice_pos_2,best_segmentj_nmismatches,ncolordiffs,
					  best_acceptor2_prob,/*left*/segmentj_left,query_compress,genome_blocks,snp_blocks,
					  querylength,plusp,sensep,query,segmentj->chrnum,segmentj->chroffset,dibasep,cmetp);

	if (donor == NULL || shortexon == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (shortexon != NULL) Substring_free(&shortexon);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  *nhits += 1;
	  hits = List_push(hits,(void *) Stage3_new_shortexon(&(*found_score),donor,acceptor,shortexon,
							      /*acceptor_distance*/segmentm_left - segmenti_left,
							      /*donor_distance*/segmentj_left - segmentm_left,
							      /*ambi_left*/NULL,/*ambi_right*/NULL,
							      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
							      /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
							      splicing_penalty,querylength,first_read_p,sensedir));

	}

      } else {
	/* Originally from minus strand.  Complement. */
	sensep = (plusp == true) ? false : true;
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	donor = Substring_new_donor(best_donor2_knowni,/*joffset*/0,best_splice_pos_2,best_segmentj_nmismatches,ncolordiffs,
				    best_donor2_prob,/*left*/segmentj_left,query_compress,genome_blocks,snp_blocks,
				    querylength,plusp,sensep,query,segmentj->chrnum,segmentj->chroffset,dibasep,cmetp);

	shortexon = Substring_new_shortexon(best_acceptor2_knowni,best_donor1_knowni,/*joffset*/0,
					    /*acceptor_pos*/best_splice_pos_2,/*donor_pos*/best_splice_pos_1,best_segmentm_nmismatches,
					    /*ncolordiffs*/0,/*acceptor_prob*/best_acceptor2_prob,/*donor_prob*/best_donor1_prob,
					    /*left*/segmentm_left,query_compress,genome_blocks,snp_blocks,
					    querylength,
					    plusp,sensep,/*acceptor_ambp*/false,/*donor_ambp*/false,
					    query,segmentm->chrnum,segmentm->chroffset,dibasep,cmetp);

	acceptor = Substring_new_acceptor(best_acceptor1_knowni,/*joffset*/0,best_splice_pos_1,best_segmenti_nmismatches,ncolordiffs,
					  best_acceptor1_prob,/*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
					  querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

	if (donor == NULL || shortexon == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (shortexon != NULL) Substring_free(&shortexon);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  *nhits += 1;
	  hits = List_push(hits,(void *) Stage3_new_shortexon(&(*found_score),donor,acceptor,shortexon,
							      /*acceptor_distance*/segmentj_left - segmentm_left,
							      /*donor_distance*/segmentm_left - segmenti_left,
							      /*ambi_left*/NULL,/*ambi_right*/NULL,
							      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
							      /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
							      splicing_penalty,querylength,first_read_p,sensedir));
	}
      }
    }
  }

  return hits;
}



static List_T
find_singlesplices_plus (int *found_score, List_T hits, struct Segment_T *plus_segments, int plus_nsegments,
			 UINT4 *genome_blocks, UINT4 *snp_blocks,
			 char *queryuc_ptr, Floors_T floors, int querylength, int query_lastpos,
			 Compress_T query_compress /* expecting fwd */,
			 Genomicpos_T *splicesites_orig, Splicetype_T *splicetypes_orig, Genomicpos_T *splicedists_orig,
			 int nsplicesites_orig, bool novelsplicingp, Genomicpos_T overall_max_distance,
			 int splicing_penalty, int min_localsplicing_end_matches,
			 int max_mismatches_allowed, bool first_read_p, bool dibasep, bool cmetp) {
#ifdef DEBUG4S
  int i;
#endif
  int j, joffset = 0;
  Segment_T segmenti, segmentj;
  Genomicpos_T segmenti_left, segmentj_left;
  int mismatch_positions_left[MAX_QUERYLENGTH], mismatch_positions_right[MAX_QUERYLENGTH];
  int colordiffs_left[MAX_QUERYLENGTH], colordiffs_right[MAX_QUERYLENGTH];
  int nmismatches_left, nmismatches_right;
  int segmenti_donor_knownpos[MAX_QUERYLENGTH+1], segmentj_acceptor_knownpos[MAX_QUERYLENGTH+1],
    segmentj_antidonor_knownpos[MAX_QUERYLENGTH+1], segmenti_antiacceptor_knownpos[MAX_QUERYLENGTH+1];
  int segmenti_donor_knowni[MAX_QUERYLENGTH+1], segmentj_acceptor_knowni[MAX_QUERYLENGTH+1],
    segmentj_antidonor_knowni[MAX_QUERYLENGTH+1], segmenti_antiacceptor_knowni[MAX_QUERYLENGTH+1];
  int segmenti_donor_nknown, segmentj_acceptor_nknown,
    segmentj_antidonor_nknown, segmenti_antiacceptor_nknown;
  
  Genomicpos_T *splicesites;
  Splicetype_T *splicetypes;
  Genomicpos_T *splicedists;
  int nsplicesites = nsplicesites_orig;
  Genomicpos_T max_distance;

  int floor_outer_i;
  int *floors_from_neg3, *floors_to_pos3;
  int nhits, npotential;

  debug(printf("*** Starting find_singlesplices_plus on %d segments with max_distance %u ***\n",
	       plus_nsegments,max_distance));
  debug(printf("Initially have %d hits\n",List_length(hits)));

  nhits = List_length(hits);

  if (plus_nsegments > 1) {
    floors_from_neg3 = floors->scorefrom[-3];
    floors_to_pos3 = floors->scoreto[query_lastpos+3];

    for (segmenti = plus_segments; segmenti < &(plus_segments[plus_nsegments]) && nhits < MAX_LOCALSPLICING_HITS; segmenti++) {
      segmenti_left = segmenti->diagonal - querylength;
      floor_outer_i = floors_from_neg3[segmenti->querypos5];

      /* Mark known splice sites first, so we know max_distance */
      joffset = segmenti->splicesites_i;
      splicesites = &(splicesites_orig[joffset]);
      splicetypes = &(splicetypes_orig[joffset]);
      splicedists = &(splicedists_orig[joffset]);
      nsplicesites = nsplicesites_orig - joffset;
      max_distance = overall_max_distance;

      /* Ends 1 (donor, plus) and 8 (antiacceptor, plus): mark known splice sites in segmenti */
      j = 0;
      segmenti_donor_nknown = 0;
      segmenti_antiacceptor_nknown = 0;
      while (j < nsplicesites && splicesites[j] < segmenti->diagonal) {
	if (splicetypes[j] == DONOR) {
	  debug4s(printf("Setting known donor %d for segmenti at %u\n",joffset+j,splicesites[j]));
	  segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - segmenti_left;
	  segmenti_donor_knowni[segmenti_donor_nknown++] = joffset + j;
	} else if (splicetypes[j] == ANTIACCEPTOR) {
	  debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",joffset+j,splicesites[j]));
	  segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - segmenti_left;
	  segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = joffset + j;
	}

	if (splicedists[j] > max_distance) {
	  debug4s(printf("Setting max_distance for known i %d to be %u\n",joffset+j,splicedists[j]));
	  max_distance = splicedists[j];
	}

	j++;
      }
      segmenti_donor_knownpos[segmenti_donor_nknown] = MAX_QUERYLENGTH;
      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = MAX_QUERYLENGTH;


      /* Identify potential segmentj for segmenti */
      npotential = 0;		/* Number of potential segmentj */
      for (segmentj = segmenti+1; segmentj < &(plus_segments[plus_nsegments]) &&
	     segmentj->diagonal <= segmenti->diagonal + max_distance &&
	     segmentj->chrnum == segmenti->chrnum &&
	     npotential++ < MAX_LOCALSPLICING_POTENTIAL; segmentj++) {
	debug4s(printf("plus local?  diagonal %u, querypos %d..%d => diagonal %u, querypos %d..%d => ",
		       segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
		       segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	/* i5 i3 j5 j3 */
	if (segmenti->querypos3 >= segmentj->querypos5) {
	  /* Fail querypos test */
	  debug4s(printf("Bad querypos\n"));

	} else if (floor_outer_i + floors_to_pos3[segmentj->querypos3] > max_mismatches_allowed) {
	  /* Fail outer floor test */
	  /* floors->score[-3][segmenti->querypos5] +floors->score[segmentj->querypos3][query_lastpos+3] */

	  debug4s(printf("too many mismatches, outer floor = %d+%d=%d > %d\n",
			 floors->scorefrom[-3][segmenti->querypos5],
			 floors->scorefrom[segmentj->querypos3][query_lastpos+3],
			 floors->scorefrom[-3][segmenti->querypos5] +
			 floors->scorefrom[segmentj->querypos3][query_lastpos+3],
			 max_mismatches_allowed));

	} else {
	  /* Apply leftmost/rightmost test */
	  if (segmenti->leftmost < 0) {
	    nmismatches_left = Genome_mismatches_left(mismatch_positions_left,colordiffs_left,max_mismatches_allowed,
						      /*query*/queryuc_ptr,query_compress,genome_blocks,snp_blocks,
						      /*left*/segmenti_left,/*pos5*/0,/*pos3*/querylength,
						      dibasep,cmetp,/*plusp*/true);
	    segmenti->leftmost = (nmismatches_left == 0) ? 0 : mismatch_positions_left[nmismatches_left-1];
	    debug4s(printf("%d mismatches on left at:",nmismatches_left);
		    for (i = 0; i <= nmismatches_left; i++) {
		      printf(" %d",mismatch_positions_left[i]);
		    }
		    printf("\n"));
	  }
	  
	  segmentj_left = segmentj->diagonal - querylength;
	  if (segmentj->rightmost < 0) {
	    nmismatches_right = Genome_mismatches_right(mismatch_positions_right,colordiffs_right,max_mismatches_allowed,
							/*query*/queryuc_ptr,query_compress,genome_blocks,snp_blocks,
							/*left*/segmentj_left,/*pos5*/0,/*pos3*/querylength,
							dibasep,cmetp,/*plusp*/true);
	    segmentj->rightmost = (nmismatches_right == 0) ? 0 : mismatch_positions_right[nmismatches_right-1];
	    debug4s(printf("%d mismatches on right at:",nmismatches_right);
		    for (i = 0; i <= nmismatches_right; i++) {
		      printf(" %d",mismatch_positions_right[i]);
		    }
		    printf("\n"));
	  }
	  
	  debug4s(printf("For a single splice, want leftmost %d > rightmost %d\n",segmenti->leftmost,segmentj->rightmost));
	    
	  if (segmenti->leftmost > segmentj->rightmost) {
	    /* Single splice is possible */

	    /* Ends 2 (acceptor, plus) and 7 (antidonor, plus): mark known splice sites in segmentj */
	    j = 0;
	    segmentj_acceptor_nknown = 0;
	    segmentj_antidonor_nknown = 0;
	    while (j < nsplicesites && splicesites[j] < segmentj_left) {
	      j++;			/* Advance to this segmentj */
	    }
	    while (j < nsplicesites && splicesites[j] < segmentj->diagonal) {
	      if (splicetypes[j] == ACCEPTOR) {
		debug4s(printf("Setting known acceptor %d for segmentj at %u\n",joffset+j,splicesites[j]));
		segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - segmentj_left;
		segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = joffset + j;
	      } else if (splicetypes[j] == ANTIDONOR) {
		debug4s(printf("Setting known antidonor %d for segmentj at %u\n",joffset+j,splicesites[j]));
		segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - segmentj_left;
		segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = joffset + j;
	      }
	      j++;
	    }
	    segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = MAX_QUERYLENGTH;
	    segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = MAX_QUERYLENGTH;


	    debug4s(printf("  => checking for single splice: solve_splicepair_local_plus\n"));
	    hits = solve_singlesplice(&(*found_score),&nhits,hits,segmenti,segmentj,
				      genome_blocks,snp_blocks,
				      /*query*/queryuc_ptr,querylength,query_compress,
				      segmenti_donor_knownpos,segmentj_acceptor_knownpos,
				      segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
				      segmenti_donor_knowni,segmentj_acceptor_knowni,
				      segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
				      segmenti_donor_nknown,segmentj_acceptor_nknown,
				      segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
				      novelsplicingp,
				      splicing_penalty,min_localsplicing_end_matches,max_mismatches_allowed,
				      first_read_p,dibasep,cmetp,/*plusp*/true);
	  }
	}
      }
    }
  }

  debug(printf("Finished find_singlesplices_plus with %d hits\n",List_length(hits)));

  return hits;
}


static List_T
find_singlesplices_minus (int *found_score, List_T hits, struct Segment_T *minus_segments, int minus_nsegments,
			  UINT4 *genome_blocks, UINT4 *snp_blocks,
			  char *query, char *queryptr /* expecting queryrc */,
			  Floors_T floors, int querylength, int query_lastpos, Compress_T query_compress /* expecting rev */,
			  Genomicpos_T *splicesites_orig, Splicetype_T *splicetypes_orig, Genomicpos_T *splicedists_orig,
			  int nsplicesites_orig, bool novelsplicingp, Genomicpos_T overall_max_distance,
			  int splicing_penalty, int min_localsplicing_end_matches,
			  int max_mismatches_allowed, bool first_read_p, bool dibasep, bool cmetp) {
#ifdef DEBUG4S
  int i;
#endif
  int j, joffset = 0;
  Segment_T segmenti, segmentj;
  Genomicpos_T segmenti_left, segmentj_left;
  int mismatch_positions_left[MAX_QUERYLENGTH], mismatch_positions_right[MAX_QUERYLENGTH];
  int colordiffs_left[MAX_QUERYLENGTH], colordiffs_right[MAX_QUERYLENGTH];
  int nmismatches_left, nmismatches_right;
  int segmenti_donor_knownpos[MAX_QUERYLENGTH+1], segmentj_acceptor_knownpos[MAX_QUERYLENGTH+1],
    segmentj_antidonor_knownpos[MAX_QUERYLENGTH+1], segmenti_antiacceptor_knownpos[MAX_QUERYLENGTH+1];
  int segmenti_donor_knowni[MAX_QUERYLENGTH+1], segmentj_acceptor_knowni[MAX_QUERYLENGTH+1],
    segmentj_antidonor_knowni[MAX_QUERYLENGTH+1], segmenti_antiacceptor_knowni[MAX_QUERYLENGTH+1];
  int segmenti_donor_nknown, segmentj_acceptor_nknown,
    segmentj_antidonor_nknown, segmenti_antiacceptor_nknown;

  Genomicpos_T *splicesites;
  Splicetype_T *splicetypes;
  Genomicpos_T *splicedists;
  int nsplicesites = nsplicesites_orig;
  Genomicpos_T max_distance;

  int floor_outer_i;
  int *floors_from_neg3, *floors_to_pos3;
  int nhits, npotential;


  debug(printf("*** Starting find_singlesplices_minus on %d segments with max_distance %u ***\n",
	       minus_nsegments,max_distance));
  debug(printf("Initially have %d hits\n",List_length(hits)));

  if (minus_nsegments > 1) {
    nhits = List_length(hits);

    floors_from_neg3 = floors->scorefrom[-3];
    floors_to_pos3 = floors->scoreto[query_lastpos+3];

    for (segmenti = minus_segments; segmenti < &(minus_segments[minus_nsegments]) && nhits < MAX_LOCALSPLICING_HITS; segmenti++) {
      segmenti_left = segmenti->diagonal - querylength;
      floor_outer_i = floors_to_pos3[segmenti->querypos3];

      /* Mark known splice sites first, so we know max_distance */
      joffset = segmenti->splicesites_i;
      splicesites = &(splicesites_orig[joffset]);
      splicetypes = &(splicetypes_orig[joffset]);
      splicedists = &(splicedists_orig[joffset]);
      nsplicesites = nsplicesites_orig - joffset;
      max_distance = overall_max_distance;
      debug4s(printf("Advanced to joffset = %d/%d (minus_segment)\n",joffset,nsplicesites_orig));

      /* Ends 4 and 5: mark known splice sites in segmenti */
      j = 0;
      segmenti_antiacceptor_nknown = 0;
      segmenti_donor_nknown = 0;
      while (j < nsplicesites && splicesites[j] < segmenti->diagonal) {
	if (splicetypes[j] == ANTIACCEPTOR) {
	  debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",joffset+j,splicesites[j]));
	  segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - segmenti_left;
	  segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = joffset + j;
	} else if (splicetypes[j] == DONOR) {
	  debug4s(printf("Setting known donor %d for segmenti at %u\n",joffset+j,splicesites[j]));
	  segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - segmenti_left;
	  segmenti_donor_knowni[segmenti_donor_nknown++] = joffset + j;
	}

	if (splicedists[j] > max_distance) {
	  debug4s(printf("Setting max_distance for known %d to be %u\n",joffset+j,splicedists[j]));
	  max_distance = splicedists[j];
	}

	j++;
      }
      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = MAX_QUERYLENGTH;
      segmenti_donor_knownpos[segmenti_donor_nknown] = MAX_QUERYLENGTH;


      /* Identify potential segmentj for segmenti */
      npotential = 0;		/* Number of potential segmentj and segmentm */
      for (segmentj = segmenti+1; segmentj < &(minus_segments[minus_nsegments]) &&
	     segmentj->diagonal <= segmenti->diagonal + max_distance &&
	     segmentj->chrnum == segmenti->chrnum &&
	     npotential++ < MAX_LOCALSPLICING_POTENTIAL; segmentj++) {
	debug4s(printf("minus local?  diagonal %u, querypos %d..%d => diagonal %u, querypos %d..%d => ",
		       segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
		       segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	/* j5 j3 i5 i3 */
	if (segmentj->querypos3 >= segmenti->querypos5) {
	  /* Fail querypos test */
	  debug4s(printf("Bad querypos\n"));

	} else if (floors_from_neg3[segmentj->querypos5] + floor_outer_i > max_mismatches_allowed) {
	  /* Fail outer floor test */
	  /* floors->score[-3][segmentj->querypos5] + floors->score[segmenti->querypos3][query_lastpos+3] */;
	  
	  debug4s(printf("too many mismatches, outer floor = %d+%d=%d > %d\n",
			 floors->scorefrom[-3][segmentj->querypos5],
			 floors->scorefrom[segmenti->querypos3][query_lastpos+3],
			 floors->scorefrom[-3][segmentj->querypos5] +
			 floors->scorefrom[segmenti->querypos3][query_lastpos+3],
			 max_mismatches_allowed));

	} else {
	  /* Apply leftmost/rightmost test */
	  if (segmenti->leftmost < 0) {
	    nmismatches_left = Genome_mismatches_left(mismatch_positions_left,colordiffs_left,max_mismatches_allowed,
						      query,query_compress,genome_blocks,snp_blocks,
						      /*left*/segmenti_left,/*pos5*/0,/*pos3*/querylength,
						      dibasep,cmetp,/*plusp*/false);
	    segmenti->leftmost = (nmismatches_left == 0) ? 0 : mismatch_positions_left[nmismatches_left-1];
	    debug4s(printf("%d mismatches on left at:",nmismatches_left);
		    for (i = 0; i <= nmismatches_left; i++) {
		      printf(" %d",mismatch_positions_left[i]);
		    }
		    printf("\n"));
	  }

	  segmentj_left = segmentj->diagonal - querylength;
	  if (segmentj->rightmost < 0) {
	    nmismatches_right = Genome_mismatches_right(mismatch_positions_right,colordiffs_right,max_mismatches_allowed,
							query,query_compress,genome_blocks,snp_blocks,
							/*left*/segmentj_left,/*pos5*/0,/*pos3*/querylength,
							dibasep,cmetp,/*plusp*/false);
	    segmentj->rightmost = (nmismatches_right == 0) ? 0 : mismatch_positions_right[nmismatches_right-1];
	    debug4s(printf("%d mismatches on right at:",nmismatches_right);
		    for (i = 0; i <= nmismatches_right; i++) {
		      printf(" %d",mismatch_positions_right[i]);
		    }
		    printf("\n"));
	  }

	  debug4s(printf("For a single splice, want leftmost %d > rightmost %d\n",segmenti->leftmost,segmentj->rightmost));

	  if (segmenti->leftmost > segmentj->rightmost) {
	    /* Single splice is possible */

	    /* Ends 3 and 6: mark known splice sites in segmentj */
	    j = 0;
	    segmentj_antidonor_nknown = 0;
	    segmentj_acceptor_nknown = 0;
	    while (j < nsplicesites && splicesites[j] < segmentj_left) {
	      j++;			/* Advance to this segmentj */
	    }
	    while (j < nsplicesites && splicesites[j] < segmentj->diagonal) {
	      if (splicetypes[j] == ANTIDONOR) {
		debug4s(printf("Setting known antidonor %d for segmentj at %u\n",joffset+j,splicesites[j]));
		segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - segmentj_left;
		segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = joffset + j;
	      } else if (splicetypes[j] == ACCEPTOR) {
		debug4s(printf("Setting known acceptor %d for segmentj at %u\n",joffset+j,splicesites[j]));
		segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - segmentj_left;
		segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = joffset + j;
	      }
	      j++;
	    }
	    segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = MAX_QUERYLENGTH;
	    segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = MAX_QUERYLENGTH;

	    debug4s(printf("  => checking for single splice: solve_singlesplice_minus\n"));
	    hits = solve_singlesplice(&(*found_score),&nhits,hits,segmenti,segmentj,
				      genome_blocks,snp_blocks,
				      query,querylength,query_compress,
				      segmenti_donor_knownpos,segmentj_acceptor_knownpos,
				      segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
				      segmenti_donor_knowni,segmentj_acceptor_knowni,
				      segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
				      segmenti_donor_nknown,segmentj_acceptor_nknown,
				      segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
				      novelsplicingp,splicing_penalty,min_localsplicing_end_matches,max_mismatches_allowed,
				      first_read_p,dibasep,cmetp,/*plusp*/false);
	  }
	}
      }
    }
  }

  debug(printf("Finished find_singlesplices_minus with %d hits\n",List_length(hits)));

  return hits;
}


static List_T
find_known_doublesplices (int *found_score, List_T hits, struct Segment_T *segments, int nsegments,
			  UINT4 *genome_blocks, UINT4 *snp_blocks,
			  char *query, char *queryptr, Floors_T floors,
			  int querylength, int query_lastpos, Compress_T query_compress,
			  Genomicpos_T *splicesites_orig, Splicetype_T *splicetypes_orig,
			  UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int nsplicesites_orig,
			  int *nsplicepartners_skip,
			  unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
			  unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
			  List_T *splicestrings,  Genomicpos_T max_distance, int splicing_penalty,
			  int min_localsplicing_end_matches, int min_shortend,
			  int max_mismatches_allowed, bool first_read_p, bool dibasep, bool cmetp, bool plusp) {
#ifdef DEBUG4S
  int i;
#endif
  Genomicpos_T *splicesites;
  Splicetype_T *splicetypes;
  int nsplicesites = nsplicesites_orig;

  int j, j1, j2, joffset = 0;
  Segment_T segmenti;
  Genomicpos_T segmenti_left;
  
  Intlist_T splicesites_i_left, splicesites_i_right;
  Intlist_T nmismatches_list_left, nmismatches_list_right;
  bool ambp_left, ambp_right;
  bool sensep;
  int sensedir;
  int *floors_from_neg3, *floors_to_pos3;
  int nhits;

  int nmismatches_shortexon_left, nmismatches_shortexon_middle, nmismatches_shortexon_right, ncolordiffs;
  int best_left_j, best_right_j;
  bool shortexon_orig_plusp, shortexon_orig_minusp, saw_antidonor_p, saw_acceptor_p;
  int leftpos, rightpos;
  Substring_T donor, acceptor, shortexon;


  debug(printf("*** Starting find_known_doublesplices on %d segments with max_distance %u ***\n",
	       nsegments,max_distance));
  debug(printf("Initially have %d hits\n",List_length(hits)));

  nhits = List_length(hits);

  if (nsegments > 0) {
    floors_from_neg3 = floors->scorefrom[-3];
    floors_to_pos3 = floors->scoreto[query_lastpos+3];

    for (segmenti = segments; segmenti < &(segments[nsegments]) && nhits < MAX_LOCALSPLICING_HITS; segmenti++) {
      segmenti_left = segmenti->diagonal - querylength;

      /* Advance through known splice sites */
      joffset = segmenti->splicesites_i;
      splicesites = &(splicesites_orig[joffset]);
      splicetypes = &(splicetypes_orig[joffset]);
      nsplicesites = nsplicesites_orig - joffset;
      debug4s(printf("Advanced to joffset = %d/%d\n",joffset,nsplicesites_orig));
      
      shortexon_orig_plusp = shortexon_orig_minusp = false;
      saw_acceptor_p = saw_antidonor_p = false;
      j = 0;
      while (j < nsplicesites && splicesites[j] < segmenti->diagonal) {
	if (splicetypes[j] == DONOR) {
	  debug4s(printf("Setting known donor %d for segmenti at %u\n",joffset+j,splicesites[j]));
	  if (saw_acceptor_p == true) {
	    /* acceptor...donor */
	    shortexon_orig_plusp = true;
	  }
	} else if (splicetypes[j] == ANTIACCEPTOR) {
	  debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",joffset+j,splicesites[j]));
	  if (saw_antidonor_p == true) {
	    /* antidonor...antiacceptor */
	    shortexon_orig_minusp = true;
	  }
	} else if (splicetypes[j] == ACCEPTOR) {
	  debug4s(printf("Saw known acceptor at %u\n",splicesites[j]));
	  saw_acceptor_p = true;
	} else if (splicetypes[j] == ANTIDONOR) {
	  debug4s(printf("Saw known antidonor at %u\n",splicesites[j]));
	  saw_antidonor_p = true;
	}
	j++;
      }


      /* Short exon, originally on plus strand */
      if (shortexon_orig_plusp == true) {
	debug4s(printf("Short exon candidate, orig_plusp.  Saw short exon acceptor...donor on segment i\n"));
	sensep = (plusp == true) ? true : false;
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;
	for (j1 = 0; j1 < j; j1++) {
	  if (splicetypes[j1] == ACCEPTOR) {
	    leftpos = splicesites[j1] - segmenti_left;
	    debug4s(printf("  Doing Splicetrie_find_left from leftpos %u\n",leftpos));
	    if ((splicesites_i_left =
		 Splicetrie_find_left(&nmismatches_shortexon_left,&nmismatches_list_left,joffset+j1,
				      /*origleft*/segmenti_left,/*pos5*/0,/*pos3*/leftpos,
				      splicesites_orig,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
				      trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				      trieoffsets_max,triecontents_max,nsplicepartners_max,
				      splicetypes_orig,splicestrings,
				      query_compress,genome_blocks,snp_blocks,query,queryptr,
				      max_mismatches_allowed,dibasep,cmetp,plusp)) != NULL) {
	      ambp_left = (leftpos < min_shortend || Intlist_length(splicesites_i_left) > 1) ? true : false;

	      for (j2 = j1 + 1; j2 < j; j2++) {
		if (splicetypes[j2] == DONOR && splicesites[j2] > splicesites[j1]) {
		  rightpos = splicesites[j2] - segmenti_left;
		  debug4s(printf("  Doing Splicetrie_find_right from rightpos %u\n",rightpos));
		  if ((nmismatches_shortexon_middle =
		       Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
							 segmenti_left,/*pos5*/leftpos,/*pos3*/rightpos,
							 dibasep,cmetp,plusp)) <= max_mismatches_allowed - nmismatches_shortexon_left &&
		      (splicesites_i_right =
		       Splicetrie_find_right(&nmismatches_shortexon_right,&nmismatches_list_right,joffset+j2,
					     /*origleft*/segmenti_left,/*pos5*/rightpos,/*pos3*/querylength,
					     splicesites_orig,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
					     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
					     trieoffsets_max,triecontents_max,nsplicepartners_max,
					     splicetypes_orig,splicestrings,
					     query_compress,genome_blocks,snp_blocks,query,queryptr,
					     max_mismatches_allowed - nmismatches_shortexon_left - nmismatches_shortexon_middle,
					     dibasep,cmetp,plusp)) != NULL) {
		    ambp_right = (querylength - rightpos < min_shortend || Intlist_length(splicesites_i_right) > 1) ? true : false;

		    debug4s(printf("  donor %s ... acceptor %d (%u) ... donor %d (%u) ... acceptor %s: %d + %d + %d mismatches\n",
				   Intlist_to_string(splicesites_i_left),j1,splicesites[j1],j2,splicesites[j2],Intlist_to_string(splicesites_i_right),
				   nmismatches_shortexon_left,nmismatches_shortexon_middle,nmismatches_shortexon_right));

		    if (ambp_left == true && ambp_right == true) {
		      shortexon = Substring_new_shortexon(j1,j2,joffset,/*acceptor_pos*/leftpos,/*donor_pos*/rightpos,
							  nmismatches_shortexon_middle,
							  /*ncolordiffs*/0,/*acceptor_prob*/2.0,/*donor_prob*/2.0,
							  /*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
							  querylength,plusp,sensep,/*acceptor_ambp*/true,/*donor_ambp*/true,
							  query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);
		      if (shortexon != NULL) {
			debug(printf("New one-third shortexon at left %u\n",segmenti_left));
			hits = List_push(hits,(void *) Stage3_new_shortexon(&(*found_score),/*donor*/NULL,/*acceptor*/NULL,shortexon,
									    /*acceptor_distance*/0U,/*donor_distance*/0U,
									    /*ambi_left*/splicesites_i_left,/*ambi_right*/splicesites_i_right,
									    nmismatches_list_left,nmismatches_list_right,
									    /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
									    splicing_penalty,querylength,first_read_p,sensedir));
		      }

		    } else if (ambp_left == true && ambp_right == false) {
		      debug(printf("ambp_left true, ambp_right false\n"));
		      best_right_j = Intlist_head(splicesites_i_right);

		      debug(printf("shortexon with amb_acceptor at %d+%d (%u) ... donor at %d+%d (%u)\n",
				   joffset,j1,splicesites[j1],joffset,j2,splicesites[j2]));
		      shortexon = Substring_new_shortexon(j1,j2,joffset,/*acceptor_pos*/leftpos,/*donor_pos*/rightpos,
							  nmismatches_shortexon_middle,
							  /*ncolordiffs*/0,/*acceptor_prob*/2.0,/*donor_prob*/2.0,
							  /*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
							  querylength,plusp,sensep,/*acceptor_ambp*/true,/*donor_ambp*/false,
							  query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      debug(printf("acceptor at %d (%u)\n",best_right_j,splicesites_orig[best_right_j]));
		      acceptor = Substring_new_acceptor(best_right_j,/*joffset*/0,/*splice_pos*/rightpos,nmismatches_shortexon_right,
							/*ncolordiffs*/0,/*prob*/2.0,/*left*/splicesites_orig[best_right_j]-rightpos,
							query_compress,genome_blocks,snp_blocks,
							querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      if (shortexon == NULL || acceptor == NULL) {
			if (shortexon != NULL) Substring_free(&shortexon);
			if (acceptor != NULL) Substring_free(&acceptor);
		      } else {
			debug(printf("ambp_left true, ambp_right false: New two-thirds shortexon at left %u\n",segmenti_left));
			hits = List_push(hits,(void *) Stage3_new_shortexon(&(*found_score),/*donor*/NULL,acceptor,shortexon,
									    /*acceptor_distance*/0U,
									    /*donor_distance*/splicesites_orig[best_right_j]-splicesites[j2],
									    /*ambi_left*/splicesites_i_left,/*ambi_right*/NULL,
									    nmismatches_list_left,/*amb_nmismatches_right*/NULL,
									    /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
									    splicing_penalty,querylength,first_read_p,sensedir));
		      }

		    } else if (ambp_left == false && ambp_right == true) {
		      debug(printf("ambp_left false, ambp_right true\n"));
		      best_left_j = Intlist_head(splicesites_i_left);

		      debug(printf("donor at %d (%u)\n",best_left_j,splicesites_orig[best_left_j]));
		      donor = Substring_new_donor(best_left_j,/*joffset*/0,/*splice_pos*/leftpos,nmismatches_shortexon_left,
						  /*ncolordiffs*/0,/*prob*/2.0,/*left*/splicesites_orig[best_left_j]-leftpos,
						  query_compress,genome_blocks,snp_blocks,
						  querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      debug(printf("shortexon with acceptor at %d+%d (%u) ... amb_donor %d+%d (%u)\n",
				   joffset,j1,splicesites[j1],joffset,j2,splicesites[j2]));
		      shortexon = Substring_new_shortexon(j1,j2,joffset,/*acceptor_pos*/leftpos,/*donor_pos*/rightpos,
							  nmismatches_shortexon_middle,
							  /*ncolordiffs*/0,/*acceptor_prob*/2.0,/*donor_prob*/2.0,
							  /*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
							  querylength,plusp,sensep,/*acceptor_ambp*/false,/*donor_ambp*/true,
							  query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      if (donor == NULL || shortexon == NULL) {
			if (donor != NULL) Substring_free(&donor);
			if (shortexon != NULL) Substring_free(&shortexon);
		      } else {
			hits = List_push(hits,(void *) Stage3_new_shortexon(&(*found_score),donor,/*acceptor*/NULL,shortexon,
									    /*acceptor_distance*/splicesites[j1]-splicesites_orig[best_left_j],
									    /*donor_distance*/0U,
									    /*ambi_left*/NULL,/*ambi_right*/splicesites_i_right,
									    /*amb_nmismatches_left*/NULL,nmismatches_list_right,
									    /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
									    splicing_penalty,querylength,first_read_p,sensedir));
		      }

		    } else { /* ambp_left == false && ambp_right == false */
		      debug(printf("ambp_left false, ambp_right false\n"));
		      best_left_j = Intlist_head(splicesites_i_left);
		      best_right_j = Intlist_head(splicesites_i_right);
		      donor = Substring_new_donor(best_left_j,/*joffset*/0,/*splice_pos*/leftpos,nmismatches_shortexon_left,
						  /*ncolordiffs*/0,/*prob*/2.0,/*left*/splicesites_orig[best_left_j]-leftpos,
						  query_compress,genome_blocks,snp_blocks,
						  querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      shortexon = Substring_new_shortexon(j1,j2,joffset,/*acceptor_pos*/leftpos,/*donor_pos*/rightpos,
							  nmismatches_shortexon_middle,/*ncolordiffs*/0,/*acceptor_prob*/2.0,/*donor_prob*/2.0,
							  /*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
							  querylength,plusp,sensep,/*acceptor_ambp*/false,/*donor_ambp*/false,
							  query,segmenti->chrnum,segmenti->chroffset,
							  dibasep,cmetp);
		      
		      acceptor = Substring_new_acceptor(best_right_j,/*joffset*/0,/*splice_pos*/rightpos,nmismatches_shortexon_right,
							/*ncolordiffs*/0,/*prob*/2.0,/*left*/splicesites_orig[best_right_j]-rightpos,
							query_compress,genome_blocks,snp_blocks,
							querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      if (donor == NULL || shortexon == NULL || acceptor == NULL) {
			if (donor != NULL) Substring_free(&donor);
			if (shortexon != NULL) Substring_free(&shortexon);
			if (acceptor != NULL) Substring_free(&acceptor);
		      } else {
			debug(printf("New shortexon at left %u\n",segmenti_left));
			hits = List_push(hits,(void *) Stage3_new_shortexon(&(*found_score),donor,acceptor,shortexon,
									    /*acceptor_distance*/splicesites[j1]-splicesites_orig[best_left_j],
									    /*donor_distance*/splicesites_orig[best_right_j]-splicesites[j2],
									    /*ambi_left*/NULL,/*ambi_right*/NULL,
									    /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
									    /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
									    splicing_penalty,querylength,first_read_p,sensedir));
		      }
		    }
		    Intlist_free(&nmismatches_list_right);
		    Intlist_free(&splicesites_i_right);
		  }
		}
	      }
	      Intlist_free(&nmismatches_list_left);
	      Intlist_free(&splicesites_i_left);
	    }
	  }
	}
	debug4s(printf("End of case 1\n"));
      }

      /* Short exon, originally on minus strand */
      if (shortexon_orig_minusp == true) {
	debug4s(printf("Short exon candidate, orig_minusp.  Saw short exon antidonor...antiacceptor on segment i\n"));
	sensep = (plusp == true) ? false : true;
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	for (j1 = 0; j1 < j; j1++) {
	  if (splicetypes[j1] == ANTIDONOR) {
	    leftpos = splicesites[j1] - segmenti_left;
	    debug4s(printf("  Doing Splicetrie_find_left from leftpos %u\n",leftpos));
	    if ((splicesites_i_left =
		 Splicetrie_find_left(&nmismatches_shortexon_left,&nmismatches_list_left,joffset+j1,
				      /*origleft*/segmenti_left,/*pos5*/0,/*pos3*/leftpos,
				      splicesites_orig,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
				      trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				      trieoffsets_max,triecontents_max,nsplicepartners_max,
				      splicetypes_orig,splicestrings,
				      query_compress,genome_blocks,snp_blocks,query,queryptr,
				      max_mismatches_allowed,dibasep,cmetp,plusp)) != NULL) {
	      ambp_left = (leftpos < min_shortend || Intlist_length(splicesites_i_left) > 1) ? true : false;
	      
	      for (j2 = j1 + 1; j2 < j; j2++) {
		if (splicetypes[j2] == ANTIACCEPTOR && splicesites[j2] > splicesites[j1]) {
		  rightpos = splicesites[j2] - segmenti_left;
		  debug4s(printf("  Doing Splicetrie_find_right from rightpos %u\n",rightpos));
		  if ((nmismatches_shortexon_middle =
		       Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
							 segmenti_left,/*pos5*/leftpos,/*pos3*/rightpos,
							 dibasep,cmetp,plusp)) <= max_mismatches_allowed - nmismatches_shortexon_left &&
		      (splicesites_i_right =
		       Splicetrie_find_right(&nmismatches_shortexon_right,&nmismatches_list_right,joffset+j2,
					     /*origleft*/segmenti_left,/*pos5*/rightpos,/*pos3*/querylength,
					     splicesites_orig,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
					     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
					     trieoffsets_max,triecontents_max,nsplicepartners_max,
					     splicetypes_orig,splicestrings,
					     query_compress,genome_blocks,snp_blocks,query,queryptr,
					     max_mismatches_allowed - nmismatches_shortexon_left - nmismatches_shortexon_middle,
					     dibasep,cmetp,plusp)) != NULL) {
		    ambp_right = (querylength - rightpos < min_shortend || Intlist_length(splicesites_i_right) > 1) ? true : false;

		    debug4s(printf("  antiacceptor %s ... antidonor %d (%u) ... antiacceptor %d (%u) ... antidonor %s: %d + %d + %d mismatches\n",
				   Intlist_to_string(splicesites_i_left),j1,splicesites[j1],j2,splicesites[j2],Intlist_to_string(splicesites_i_right),
				   nmismatches_shortexon_left,nmismatches_shortexon_middle,nmismatches_shortexon_right));

		    if (ambp_left == true && ambp_right == true) {
		      shortexon = Substring_new_shortexon(j2,j1,joffset,/*acceptor_pos*/rightpos,/*donor_pos*/leftpos,nmismatches_shortexon_middle,
							  /*ncolordiffs*/0,/*acceptor_prob*/2.0,/*donor_prob*/2.0,
							  /*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
							  querylength,plusp,sensep,/*acceptor_ambp*/true,/*donor_ambp*/true,
							  query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);
		      if (shortexon != NULL) {
			debug(printf("New one-third shortexon at left %u\n",segmenti_left));
			hits = List_push(hits,(void *) Stage3_new_shortexon(&(*found_score),/*donor*/NULL,/*acceptor*/NULL,shortexon,
									    /*acceptor_distance*/0U,/*donor_distance*/0U,
									    /*ambi_left*/splicesites_i_left,/*ambi_right*/splicesites_i_right,
									    nmismatches_list_left,nmismatches_list_right,
									    /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
									    splicing_penalty,querylength,first_read_p,sensedir));
		      }

		    } else if (ambp_left == true && ambp_right == false) {
		      debug(printf("ambp_left true, ambp_right false\n"));
		      best_right_j = Intlist_head(splicesites_i_right);

		      debug(printf("shortexon with amb_donor at %d+%d (%u) ... acceptor at %d+%d (%u)\n",
				   joffset,j1,splicesites[j1],joffset,j2,splicesites[j2]));
		      shortexon = Substring_new_shortexon(j2,j1,joffset,/*acceptor_pos*/rightpos,/*donor_pos*/leftpos,nmismatches_shortexon_middle,
							  /*ncolordiffs*/0,/*acceptor_prob*/2.0,/*donor_prob*/2.0,
							  /*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
							  querylength,plusp,sensep,/*acceptor_ambp*/false,/*donor_ambp*/true,
							  query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      debug(printf("donor at %d (%u)\n",best_right_j,splicesites_orig[best_right_j]));
		      donor = Substring_new_donor(best_right_j,/*joffset*/0,/*splice_pos*/rightpos,nmismatches_shortexon_right,
						  /*ncolordiffs*/0,/*prob*/2.0,/*left*/splicesites_orig[best_right_j]-rightpos,
						  query_compress,genome_blocks,snp_blocks,
						  querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      if (donor == NULL || shortexon == NULL) {
			if (donor != NULL) Substring_free(&donor);
			if (shortexon != NULL) Substring_free(&shortexon);
		      } else {
			hits = List_push(hits,(void *) Stage3_new_shortexon(&(*found_score),donor,/*acceptor*/NULL,shortexon,
									    /*acceptor_distance*/splicesites_orig[best_right_j]-splicesites[j2],
									    /*donor_distance*/0U,
									    /*ambi_left*/splicesites_i_left,/*ambi_right*/NULL,
									    nmismatches_list_left,/*amb_nmismatches_right*/NULL,
									    /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
									    splicing_penalty,querylength,first_read_p,sensedir));
		      }

		    } else if (ambp_left == false && ambp_right == true) {
		      debug(printf("ambp_left false, ambp_right true\n"));
		      best_left_j = Intlist_head(splicesites_i_left);

		      debug(printf("acceptor at %d (%u)\n",best_left_j,splicesites_orig[best_left_j]));
		      acceptor = Substring_new_acceptor(best_left_j,/*joffset*/0,/*splice_pos*/leftpos,nmismatches_shortexon_left,
							/*ncolordiffs*/0,/*prob*/2.0,/*left*/splicesites_orig[best_left_j]-leftpos,
							query_compress,genome_blocks,snp_blocks,
							querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      debug(printf("shortexon with donor at %d+%d (%u) ... amb_acceptor at %d+%d (%u)\n",
				   joffset,j2,splicesites[j2],joffset,j1,splicesites[j1]));
		      shortexon = Substring_new_shortexon(j2,j1,joffset,/*acceptor_pos*/rightpos,/*donor_pos*/leftpos,nmismatches_shortexon_middle,
							  /*ncolordiffs*/0,/*acceptor_prob*/2.0,/*donor_prob*/2.0,
							  /*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
							  querylength,plusp,sensep,/*acceptor_ambp*/true,/*donor_ambp*/false,
							  query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      if (shortexon == NULL || acceptor == NULL) {
			if (shortexon != NULL) Substring_free(&shortexon);
			if (acceptor != NULL) Substring_free(&acceptor);
		      } else {
			debug(printf("ambp_left false, ambp_right true: New splice at left %u\n",segmenti_left));
			hits = List_push(hits,(void *) Stage3_new_shortexon(&(*found_score),/*donor*/NULL,acceptor,shortexon,
									    /*acceptor_distance*/0U,
									    /*donor_distance*/splicesites[j1]-splicesites_orig[best_left_j],
									    /*ambi_left*/NULL,/*ambi_right*/splicesites_i_right,
									    /*amb_nmismatches_left*/NULL,nmismatches_list_right,
									    /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
									    splicing_penalty,querylength,first_read_p,sensedir));
		      }

		    } else {  /* ambp_left == false && ambp_right == false */
		      best_left_j = Intlist_head(splicesites_i_left);
		      best_right_j = Intlist_head(splicesites_i_right);
		      acceptor = Substring_new_acceptor(best_left_j,/*joffset*/0,/*splice_pos*/leftpos,nmismatches_shortexon_left,
							/*ncolordiffs*/0,/*prob*/2.0,/*left*/splicesites_orig[best_left_j]-leftpos,
							query_compress,genome_blocks,snp_blocks,
							querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      shortexon = Substring_new_shortexon(j2,j1,joffset,/*acceptor_pos*/rightpos,/*donor_pos*/leftpos,
							  nmismatches_shortexon_middle,/*ncolordiffs*/0,/*acceptor_prob*/2.0,/*donor_prob*/2.0,
							  /*left*/segmenti_left,query_compress,genome_blocks,snp_blocks,
							  querylength,plusp,sensep,/*acceptor_ambp*/false,/*donor_ambp*/false,
							  query,segmenti->chrnum,segmenti->chroffset,
							  dibasep,cmetp);

		      donor = Substring_new_donor(best_right_j,/*joffset*/0,/*splice_pos*/rightpos,nmismatches_shortexon_right,
						  /*ncolordiffs*/0,/*prob*/2.0,/*left*/splicesites_orig[best_right_j]-rightpos,
						  query_compress,genome_blocks,snp_blocks,
						  querylength,plusp,sensep,query,segmenti->chrnum,segmenti->chroffset,dibasep,cmetp);

		      if (acceptor == NULL || shortexon == NULL || donor == NULL) {
			if (acceptor != NULL) Substring_free(&acceptor);
			if (shortexon != NULL) Substring_free(&shortexon);
			if (donor != NULL) Substring_free(&donor);
		      } else {
			debug(printf("New shortexon at left %u\n",segmenti_left));
			hits = List_push(hits,(void *) Stage3_new_shortexon(&(*found_score),donor,acceptor,shortexon,
									    /*acceptor_distance*/splicesites_orig[best_right_j]-splicesites[j2],
									    /*donor_distance*/splicesites[j1]-splicesites_orig[best_left_j],
									    /*ambi_left*/NULL,/*ambi_right*/NULL,
									    /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
									    /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
									    splicing_penalty,querylength,first_read_p,sensedir));
		      }
		    }
		    Intlist_free(&nmismatches_list_right);
		    Intlist_free(&splicesites_i_right);
		  }
		}
	      }
	      Intlist_free(&nmismatches_list_left);
	      Intlist_free(&splicesites_i_left);
	    }
	  }
	}
	debug4s(printf("End of case 2\n"));
      }
      /* End of known splicesites, segment i */
    }
  }

  debug(printf("Finished find_known_doublesplices with %d hits\n",List_length(hits)));
  return hits;
}


/* For known and novel doublesplices */
static List_T
find_doublesplices (int *found_score, List_T hits, struct Segment_T *segments, int nsegments,
		    UINT4 *genome_blocks, UINT4 *snp_blocks, char *query,
		    Floors_T floors, int querylength, int query_lastpos, Compress_T query_compress,
		    Genomicpos_T *splicesites_orig, Splicetype_T *splicetypes_orig,
		    int nsplicesites_orig, Genomicpos_T max_distance, int splicing_penalty,
		    int min_localsplicing_end_matches, int max_mismatches_allowed,
		    bool novelsplicingp, bool first_read_p, bool dibasep, bool cmetp, bool plusp) {
#ifdef DEBUG4S
  int i;
#endif
  int j, joffset = 0;
  int k, l;
  Segment_T segmenti, segmentj, segmentm, potentialj[MAX_LOCALSPLICING_POTENTIAL];
  Genomicpos_T segmenti_left, segmentj_left, segmentm_left;
  int segmenti_donor_knownpos[MAX_QUERYLENGTH+1], segmentj_acceptor_knownpos[MAX_QUERYLENGTH+1],
    segmentj_antidonor_knownpos[MAX_QUERYLENGTH+1], segmenti_antiacceptor_knownpos[MAX_QUERYLENGTH+1],
    segmentm_donor_knownpos[MAX_QUERYLENGTH+1], segmentm_acceptor_knownpos[MAX_QUERYLENGTH+1],
    segmentm_antidonor_knownpos[MAX_QUERYLENGTH+1], segmentm_antiacceptor_knownpos[MAX_QUERYLENGTH+1];
  int segmenti_donor_knowni[MAX_QUERYLENGTH+1], segmentj_acceptor_knowni[MAX_QUERYLENGTH+1],
    segmentj_antidonor_knowni[MAX_QUERYLENGTH+1], segmenti_antiacceptor_knowni[MAX_QUERYLENGTH+1],
    segmentm_donor_knowni[MAX_QUERYLENGTH+1], segmentm_acceptor_knowni[MAX_QUERYLENGTH+1],
    segmentm_antidonor_knowni[MAX_QUERYLENGTH+1], segmentm_antiacceptor_knowni[MAX_QUERYLENGTH+1];
  int segmenti_donor_nknown, segmentj_acceptor_nknown,
    segmentj_antidonor_nknown, segmenti_antiacceptor_nknown,
    segmentm_donor_nknown, segmentm_acceptor_nknown,
    segmentm_antidonor_nknown, segmentm_antiacceptor_nknown;
  
  Genomicpos_T *splicesites;
  Splicetype_T *splicetypes;
  int nsplicesites = nsplicesites_orig;

  int ncolordiffs;

  int *floors_from_neg3, *floors_to_pos3;
  int nhits, npotential;


  debug(printf("*** Starting find_doublesplices on %d segments with max_distance %u ***\n",
	       nsegments,max_distance));
  debug(printf("Initially have %d hits\n",List_length(hits)));

  nhits = List_length(hits);

  if (nsegments > 1) {
    floors_from_neg3 = floors->scorefrom[-3];
    floors_to_pos3 = floors->scoreto[query_lastpos+3];

    for (segmenti = segments; segmenti < &(segments[nsegments]); segmenti++) {

      npotential = 0;		/* Number of potential segmentj and segmentm */
      for (segmentj = segmenti+1; segmentj < &(segments[nsegments]) &&
	     segmentj->diagonal <= segmenti->diagonal + max_distance &&
	     segmentj->chrnum == segmenti->chrnum &&
	     npotential < MAX_LOCALSPLICING_POTENTIAL; segmentj++) {
	debug4d(printf("plus local?  diagonal %u, querypos %d..%d => diagonal %u, querypos %d..%d => ",
		       segmenti->diagonal,segmenti->querypos5,segmenti->querypos3,
		       segmentj->diagonal,segmentj->querypos5,segmentj->querypos3));
	/* i5 i3 j5 j3 */
	if (plusp == true && segmenti->querypos3 >= segmentj->querypos5) {
	  debug4d(printf("Bad querypos\n"));
	} else if (plusp == false && segmentj->querypos3 >= segmenti->querypos5) {
	  debug4d(printf("Bad querypos\n"));
	} else {
	  debug4d(printf("Potential %d\n",npotential));
	  potentialj[npotential++] = segmentj;
	}
      }

      for (k = 0; k < npotential; k++) {
	segmentj = potentialj[k];
	
	if (segmenti->leftmost >= 0 && segmentj->rightmost >= 0) {
	  /* Presence of these values indicates that they passed the outer floor test in find_singlesplices */
	  
	  /* Determine if a double splice is possible.  Use querypos and leftspan/rightspan test. */
	  for (l = 0; l < k-1; l++) {
	    segmentm = potentialj[l];
	    if (plusp == true && segmenti->querypos3 >= segmentm->querypos5) {
	      /* Bad querypos from segmenti to middle */
	    } else if (plusp == true && segmentm->querypos3 >= segmentj->querypos5) {
	      /* Bad querypos from middle to segmentj */
	    } else if (plusp == false && segmentj->querypos3 >= segmentm->querypos5) {
	      /* Bad querypos from segmentj to middle */
	    } else if (plusp == false && segmentm->querypos3 >= segmenti->querypos5) {
	      /* Bad querypos from middle to segmenti */
	    } else {
	      segmentm_left = segmentm->diagonal - querylength;
#if 0
	      /* Apply span test */
	      if (segmentm->leftspan < 0) {
		find_segmentm_span(segmentm,max_mismatches_allowed,query,querylength,
				   query_compress,genome_blocks,snp_blocks,segmentm_left,
				   dibasep,cmetp,plusp);
	      }
	      debug4d(printf("Doublesplice span test (%d mismatches allowed): i.leftmost %d >= m.leftspan %d && m.rightspan %d >= j.rightmost %d\n",
			     max_mismatches_allowed,segmenti->leftmost,segmentm->leftspan,
			     segmentm->rightspan,segmentj->rightmost));
#else
	      debug4d(printf("Doublesplice span test (%d mismatches allowed): %d mismatches found from leftmost %d to j.rightmost %d\n",
			     max_mismatches_allowed,
			     Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,
							       genome_blocks,snp_blocks,segmentm_left,
							       /*pos5*/segmenti->leftmost,/*pos3*/segmentj->rightmost,
							       dibasep,cmetp,plusp),
			     segmenti->leftmost,segmentj->rightmost));
#endif

	      if (Genome_count_mismatches_limit(&ncolordiffs,query,query_compress,
						genome_blocks,snp_blocks,segmentm_left,
						/*pos5*/segmenti->leftmost,/*pos3*/segmentj->rightmost,
						max_mismatches_allowed,dibasep,cmetp,plusp) <= max_mismatches_allowed) {
		debug4d(printf("Double splice is possible\n"));
		segmenti_left = segmenti->diagonal - querylength;
		segmentj_left = segmentj->diagonal - querylength;

		joffset = segmenti->splicesites_i;
		splicesites = &(splicesites_orig[joffset]);
		splicetypes = &(splicetypes_orig[joffset]);
		nsplicesites = nsplicesites_orig - joffset;
		debug4d(printf("Advanced to joffset = %d/%d\n",joffset,nsplicesites_orig));

		/* Set known sites for segmenti */
		j = 0;
		segmenti_donor_nknown = 0;
		segmenti_antiacceptor_nknown = 0;
		while (j < nsplicesites && splicesites[j] < segmenti->diagonal) {
		  if (splicetypes[j] == DONOR) {
		    debug4d(printf("Setting known donor %d for segmenti at %u\n",joffset+j,splicesites[j]));
		    segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - segmenti_left;
		    segmenti_donor_knowni[segmenti_donor_nknown++] = joffset + j;
		  } else if (splicetypes[j] == ANTIACCEPTOR) {
		    debug4d(printf("Setting known antiacceptor %d for segmenti at %u\n",joffset+j,splicesites[j]));
		    segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - segmenti_left;
		    segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = joffset + j;
		  }
		  j++;
		}
		segmenti_donor_knownpos[segmenti_donor_nknown] = MAX_QUERYLENGTH;
		segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = MAX_QUERYLENGTH;
	      
		/* Set known sites for segmentm */
		j = 0;
		segmentm_donor_nknown = 0;
		segmentm_acceptor_nknown = 0;
		segmentm_antidonor_nknown = 0;
		segmentm_antiacceptor_nknown = 0;
		while (j < nsplicesites && splicesites[j] < segmentm_left) {
		  j++;			/* Advance to this segmentm */
		}
		while (j < nsplicesites && splicesites[j] < segmentm->diagonal) {
		  if (splicetypes[j] == DONOR) {
		    debug4d(printf("Setting known donor %d for segmentm at %u\n",joffset+j,splicesites[j]));
		    segmentm_donor_knownpos[segmentm_donor_nknown] = splicesites[j] - segmentm_left;
		    segmentm_donor_knowni[segmentm_donor_nknown++] = joffset + j;
		  } else if (splicetypes[j] == ACCEPTOR) {
		    debug4d(printf("Setting known acceptor %d for segmentm at %u\n",joffset+j,splicesites[j]));
		    segmentm_acceptor_knownpos[segmentm_acceptor_nknown] = splicesites[j] - segmentm_left;
		    segmentm_acceptor_knowni[segmentm_acceptor_nknown++] = joffset + j;
		  } else if (splicetypes[j] == ANTIDONOR) {
		    debug4d(printf("Setting known antidonor %d for segmentm at %u\n",joffset+j,splicesites[j]));
		    segmentm_antidonor_knownpos[segmentm_antidonor_nknown] = splicesites[j] - segmentm_left;
		    segmentm_antidonor_knowni[segmentm_antidonor_nknown++] = joffset + j;
		  } else if (splicetypes[j] == ANTIACCEPTOR) {
		    debug4d(printf("Setting known antiacceptor %d for segmentm at %u\n",joffset+j,splicesites[j]));
		    segmentm_antiacceptor_knownpos[segmentm_antiacceptor_nknown] = splicesites[j] - segmentm_left;
		    segmentm_antiacceptor_knowni[segmentm_antiacceptor_nknown++] = joffset + j;
		  }
		  j++;
		}
		segmentm_donor_knownpos[segmentm_donor_nknown] = MAX_QUERYLENGTH;
		segmentm_acceptor_knownpos[segmentm_acceptor_nknown] = MAX_QUERYLENGTH;
		segmentm_antidonor_knownpos[segmentm_antidonor_nknown] = MAX_QUERYLENGTH;
		segmentm_antiacceptor_knownpos[segmentm_antiacceptor_nknown] = MAX_QUERYLENGTH;

		/* Set known sites for segmentj */
		j = 0;
		segmentj_acceptor_nknown = 0;
		segmentj_antidonor_nknown = 0;
		while (j < nsplicesites && splicesites[j] < segmentj_left) {
		  j++;			/* Advance to this segmentj */
		}
		while (j < nsplicesites && splicesites[j] < segmentj->diagonal) {
		  if (splicetypes[j] == ACCEPTOR) {
		    debug4d(printf("Setting known acceptor %d for segmentj at %u\n",joffset+j,splicesites[j]));
		    segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - segmentj_left;
		    segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = joffset + j;
		  } else if (splicetypes[j] == ANTIDONOR) {
		    debug4d(printf("Setting known antidonor %d for segmentj at %u\n",joffset+j,splicesites[j]));
		    segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - segmentj_left;
		    segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = joffset + j;
		  }
		  j++;
		}
		segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = MAX_QUERYLENGTH;
		segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = MAX_QUERYLENGTH;

		debug4d(printf("  => checking for double splice: solve_doublesplice\n"));
		hits = solve_doublesplice(&(*found_score),&nhits,hits,segmenti,segmentm,segmentj,
					  genome_blocks,snp_blocks,query,querylength,query_compress,
					  segmenti_donor_knownpos,segmentm_acceptor_knownpos,segmentm_donor_knownpos,segmentj_acceptor_knownpos,
					  segmentj_antidonor_knownpos,segmentm_antiacceptor_knownpos,segmentm_antidonor_knownpos,segmenti_antiacceptor_knownpos,
					  segmenti_donor_knowni,segmentm_acceptor_knowni,segmentm_donor_knowni,segmentj_acceptor_knowni,
					  segmentj_antidonor_knowni,segmentm_antiacceptor_knowni,segmentm_antidonor_knowni,segmenti_antiacceptor_knowni,
					  segmenti_donor_nknown,segmentm_acceptor_nknown,segmentm_donor_nknown,segmentj_acceptor_nknown,
					  segmentj_antidonor_nknown,segmentm_antiacceptor_nknown,segmentm_antidonor_nknown,segmenti_antiacceptor_nknown,
					  novelsplicingp,splicing_penalty,min_localsplicing_end_matches,max_mismatches_allowed,
					  first_read_p,dibasep,cmetp,plusp);
	      }
	    }
	  }
	}
      }
    }
  }

  return hits;
}


static void
find_spliceends (List_T **donors, List_T **antidonors, List_T **acceptors, List_T **antiacceptors,
		 struct Segment_T *segments, int nsegments, UINT4 *genome_blocks, UINT4 *snp_blocks,
		 char *query, char *queryptr, Floors_T floors, int querylength,
		 int query_lastpos, Compress_T query_compress,
		 Genomicpos_T *splicesites_orig, Splicetype_T *splicetypes_orig, int nsplicesites_orig,
		 bool novelsplicingp, bool canonicalp,
		 int max_mismatches_allowed, bool dibasep, bool cmetp, bool plusp) {
#ifdef DEBUG4E
  char gbuffer[MAX_QUERYLENGTH];
#endif

  Segment_T segment;
  Substring_T hit;
  Genomicpos_T segment_left;
  int nmismatches, joffset = 0, j, i;
  int splice_pos;
  double prob;

  int mismatch_positions[MAX_QUERYLENGTH+1], colordiffs[MAX_QUERYLENGTH];
  int nmismatches_left, nmismatches_right;
  int ncolordiffs = 0;
  int *floors_from_neg3, *floors_to_pos3;
  bool sensep, advancep, filledp;

  Genomicpos_T *splicesites;
  Splicetype_T *splicetypes;
  int nsplicesites = nsplicesites_orig;

  int splice_pos_start, splice_pos_end;

  int segment_donor_knownpos[MAX_QUERYLENGTH+1], segment_acceptor_knownpos[MAX_QUERYLENGTH+1];
  int segment_antidonor_knownpos[MAX_QUERYLENGTH+1], segment_antiacceptor_knownpos[MAX_QUERYLENGTH+1];
  int segment_donor_knowni[MAX_QUERYLENGTH+1], segment_acceptor_knowni[MAX_QUERYLENGTH+1];
  int segment_antidonor_knowni[MAX_QUERYLENGTH+1], segment_antiacceptor_knowni[MAX_QUERYLENGTH+1];
  int segment_donor_nknown, segment_acceptor_nknown, segment_antidonor_nknown, segment_antiacceptor_nknown;

  int positions_alloc[MAX_QUERYLENGTH+1];
  int knowni_alloc[MAX_QUERYLENGTH+1];
  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;


  if (nsegments > 0) {
    floors_from_neg3 = floors->scorefrom[-3];
    floors_to_pos3 = floors->scoreto[query_lastpos+3];

    for (segment = segments; segment < &(segments[nsegments]); segment++) {
      segment_left = segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
      debug4e(printf("identify_spliceends: Checking up to %d mismatches at diagonal %u (querypos %d..%d) - querylength %d = %u, floors %d and %d\n",
		     max_mismatches_allowed,segment->diagonal,segment->querypos5,segment->querypos3,querylength,segment_left,
		     floors_from_neg3[segment->querypos5],floors_to_pos3[segment->querypos3]));

      debug4e(
	      Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	      printf("genome 0..: %s\n",gbuffer);
	      printf("query  0..: %s\n",queryptr);
	      );

      filledp = false;    	/* Signals that we haven't filled buffer yet */
      advancep = false;		/* Tells minus part whether the plus part advanced through known splicesites */

      /* Splice ends from left to splice site */
      if ((plusp == true && floors_from_neg3[segment->querypos5] <= max_mismatches_allowed) ||
	  (plusp == false && floors_to_pos3[segment->querypos3] <= max_mismatches_allowed)) {

	/* pos3 was trimpos */
	nmismatches_left = Genome_mismatches_left(mismatch_positions,colordiffs,max_mismatches_allowed,
						  /*query*/query,query_compress,genome_blocks,snp_blocks,
						  /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						  dibasep,cmetp,plusp);

	debug4e(
		printf("%d mismatches on left (%d allowed) at:",
		       nmismatches_left,max_mismatches_allowed);
		for (i = 0; i <= nmismatches_left; i++) {
		  printf(" %d",mismatch_positions[i]);
		}
		printf("\n");
		);

	splice_pos_start = INDEX1PART;
	if (nmismatches_left == 0) {
	  splice_pos_end = querylength - 1;
	} else if ((splice_pos_end = mismatch_positions[nmismatches_left-1]) > querylength - 1) {
	  splice_pos_end = querylength - 1;
	}

	if (splice_pos_start <= splice_pos_end) {
	  debug4e(printf("Search for splice sites from %d up to %d\n",splice_pos_start,splice_pos_end));
	  joffset = segment->splicesites_i;
	  splicesites = &(splicesites_orig[joffset]);
	  splicetypes = &(splicetypes_orig[joffset]);
	  nsplicesites = nsplicesites_orig - joffset;
	  advancep = true;
	
	  /* Known splicing */
	  /* Ends 1 (donor, plus) and 8 (antiacceptor, plus): mark known splice sites in segment */
	  j = 0;
	  segment_donor_nknown = 0;
	  segment_antiacceptor_nknown = 0;
	  while (j < nsplicesites && splicesites[j] <= segment_left + splice_pos_end) { /* Needs to be <= */
	    if (splicetypes[j] == DONOR) {
	      segment_donor_knownpos[segment_donor_nknown] = splicesites[j] - segment_left;
	      segment_donor_knowni[segment_donor_nknown++] = joffset + j;
	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      segment_antiacceptor_knownpos[segment_antiacceptor_nknown] = splicesites[j] - segment_left;
	      segment_antiacceptor_knowni[segment_antiacceptor_nknown++] = joffset + j;
	    }
	    j++;
	  }
	  segment_donor_knownpos[segment_donor_nknown] = MAX_QUERYLENGTH;
	  segment_antiacceptor_knownpos[segment_antiacceptor_nknown] = MAX_QUERYLENGTH;


	  /* Originally on plus strand.  No complement */
	  sensep = (plusp == true) ? true : false;
	  if (novelsplicingp && segment_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
	    donori_nsites = Genome_donor_positions(positions_alloc,knowni_alloc,
						   segment_donor_knownpos,segment_donor_knowni,
						   genome_blocks,segment_left,
						   splice_pos_start,splice_pos_end+1);
	    donori_positions = positions_alloc;
	    donori_knowni = knowni_alloc;
	  } else {
	    donori_nsites = segment_donor_nknown;
	    donori_positions = segment_donor_knownpos;
	    donori_knowni = segment_donor_knowni;
	  }

	  i = 0;
	  nmismatches = 0;
	  while (i < donori_nsites) {
	    splice_pos = donori_positions[i];
	    while (nmismatches < nmismatches_left && mismatch_positions[nmismatches] < splice_pos) { /* Changed from <= to < */
	      debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	      nmismatches++;
	    }
#if 0
	    assert(nmismatches == Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								    segment_left,/*pos5*/0,/*pos3*/splice_pos,dibasep,cmetp,plusp));
#endif
	    if (nmismatches <= max_mismatches_allowed) {
	      if (donori_knowni[i] >= 0) {
		prob = 1.0;
		debug4e(printf("Known donor for segment at %u, splice_pos %d (%d mismatches), stopi = %d\n",
			       segment_left,splice_pos,nmismatches,splice_pos_end));

	      } else {
		prob = Maxent_hr_donor_prob(segment_left + splice_pos,genome_blocks);
		debug4e(printf("Novel donor for segment at %u, splice_pos %d (%d mismatches), stopi = %d\n",
			       segment_left,splice_pos,nmismatches,splice_pos_end));

	      }

	      if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob)) {
		ncolordiffs = 0;
		if ((hit = Substring_new_donor(donori_knowni[i],/*joffset*/0,splice_pos,nmismatches,ncolordiffs,
					       prob,/*left*/segment_left,query_compress,genome_blocks,snp_blocks,
					       querylength,plusp,sensep,/*query*/query,segment->chrnum,segment->chroffset,
					       dibasep,cmetp)) != NULL) {
		  debug4e(printf("=> %s donor: %f at %d (%d mismatches)\n",
				 plusp == true ? "plus" : "minus",prob,Substring_chimera_pos(hit),nmismatches));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*donors)[nmismatches] = List_push((*donors)[nmismatches],(void *) hit);
		}
	      }
	    }

	    i++;
	  }


	  /* Splicing originally on minus strand.  Complement */
	  sensep = (plusp == true) ? false : true;
	  if (novelsplicingp && segment_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
	    antiacceptori_nsites = Genome_antiacceptor_positions(positions_alloc,knowni_alloc,
								 segment_antiacceptor_knownpos,segment_antiacceptor_knowni,
								 genome_blocks,segment_left,
								 splice_pos_start,splice_pos_end+1);
	    antiacceptori_positions = positions_alloc;
	    antiacceptori_knowni = knowni_alloc;
	  } else {
	    antiacceptori_nsites = segment_antiacceptor_nknown;
	    antiacceptori_positions = segment_antiacceptor_knownpos;
	    antiacceptori_knowni = segment_antiacceptor_knowni;
	  }

	  i = 0;
	  nmismatches = 0;
	  while (i < antiacceptori_nsites) {
	    splice_pos = antiacceptori_positions[i];
	    while (nmismatches < nmismatches_left && mismatch_positions[nmismatches] < splice_pos) { /* Changed from <= to < */
	      debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	      nmismatches++;
	    }
#if 0
	    assert(nmismatches == Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								    segment_left,/*pos5*/0,/*pos3*/splice_pos,dibasep,cmetp,plusp));
#endif
	    if (nmismatches <= max_mismatches_allowed) {
	      if (antiacceptori_knowni[i] >= 0) {
		prob = 1.0;
		debug4e(printf("Known antiacceptor for segment at %u, splice_pos %d (%d mismatches), stopi = %d\n",
			       segment_left,splice_pos,nmismatches,splice_pos_end));

	      } else {
		prob = Maxent_hr_antiacceptor_prob(segment_left + splice_pos,genome_blocks);
		debug4e(printf("Novel antiacceptor for segment at %u, splice_pos %d (%d mismatches), stopi = %d\n",
			       segment_left,splice_pos,nmismatches,splice_pos_end));

	      }

	      if (sufficient_splice_prob_distant(/*support*/splice_pos,nmismatches,prob)) {
		ncolordiffs = 0;
		if ((hit = Substring_new_acceptor(antiacceptori_knowni[i],/*joffset*/0,splice_pos,nmismatches,ncolordiffs,
						  prob,/*left*/segment_left,query_compress,genome_blocks,snp_blocks,
						  querylength,plusp,sensep,query,segment->chrnum,segment->chroffset,
						  dibasep,cmetp)) != NULL) {
		  debug4e(printf("=> %s antiacceptor : %f at %d (%d mismatches)\n",
				 plusp == true ? "plus" : "minus",prob,Substring_chimera_pos(hit),nmismatches));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*antiacceptors)[nmismatches] = List_push((*antiacceptors)[nmismatches],(void *) hit);
		}
	      }
	    }

	    i++;
	  }
	}

      }

      /* Splice ends from splice site to right end */
      if ((plusp == true && floors_to_pos3[segment->querypos3] <= max_mismatches_allowed) ||
	  (plusp == false && floors_from_neg3[segment->querypos5] <= max_mismatches_allowed)){

	/* pos5 was trimpos+1 */
	nmismatches_right = Genome_mismatches_right(mismatch_positions,colordiffs,max_mismatches_allowed,
						    query,query_compress,genome_blocks,snp_blocks,
						    /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						    dibasep,cmetp,plusp);

	debug4e(
		printf("%d mismatches on right (%d allowed) at:",nmismatches_right,max_mismatches_allowed);
		for (i = 0; i <= nmismatches_right; i++) {
		  printf(" %d",mismatch_positions[i]);
		}
		printf("\n");
		);

	splice_pos_end = query_lastpos;
	if (nmismatches_right == 0) {
	  splice_pos_start = 1;
	} else if ((splice_pos_start = mismatch_positions[nmismatches_right-1]) < 1) {
	  splice_pos_start = 1;
	}

	if (splice_pos_start <= splice_pos_end) {
	  debug4e(printf("Search for splice sites from %d down to %d\n",splice_pos_end,splice_pos_start));

	  /* Advance through known splice sites */
	  if (advancep == false) {
	    joffset = segment->splicesites_i;
	    splicesites = &(splicesites_orig[joffset]);
	    splicetypes = &(splicetypes_orig[joffset]);
	    nsplicesites = nsplicesites_orig - joffset;
	    /* advancep = true; */
	  }

	  /* Known splicing */
	  j = 0;
	  segment_acceptor_nknown = 0;
	  segment_antidonor_nknown = 0;
	  while (j < nsplicesites && splicesites[j] < segment_left + splice_pos_start) {
	    j++;			/* Advance to this segment */
	  }
	  while (j < nsplicesites && splicesites[j] <= segment_left + splice_pos_end) { /* Needs to be <= */
	    if (splicetypes[j] == ACCEPTOR) {
	      debug4s(printf("Setting known acceptor %d for segment at %u\n",joffset+j,splicesites[j]));
	      segment_acceptor_knownpos[segment_acceptor_nknown] = splicesites[j] - segment_left;
	      segment_acceptor_knowni[segment_acceptor_nknown++] = joffset + j;
	    } else if (splicetypes[j] == ANTIDONOR) {
	      debug4s(printf("Setting known antidonor %d for segment at %u\n",joffset+j,splicesites[j]));
	      segment_antidonor_knownpos[segment_antidonor_nknown] = splicesites[j] - segment_left;
	      segment_antidonor_knowni[segment_antidonor_nknown++] = joffset + j;
	    }
	    j++;
	  }
	  segment_acceptor_knownpos[segment_acceptor_nknown] = MAX_QUERYLENGTH;
	  segment_antidonor_knownpos[segment_antidonor_nknown] = MAX_QUERYLENGTH;


	  /* Splicing originally on plus strand.  No complement. */
	  sensep = (plusp == true) ? true : false;
	  if (novelsplicingp && segment_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
	    acceptorj_nsites = Genome_acceptor_positions(positions_alloc,knowni_alloc,
							 segment_acceptor_knownpos,segment_acceptor_knowni,
							 genome_blocks,segment_left,
							 splice_pos_start,splice_pos_end+1);
	    acceptorj_positions = positions_alloc;
	    acceptorj_knowni = knowni_alloc;
	  } else {
	    acceptorj_nsites = segment_acceptor_nknown;
	    acceptorj_positions = segment_acceptor_knownpos;
	    acceptorj_knowni = segment_acceptor_knowni;
	  }

	  i = acceptorj_nsites - 1;
	  nmismatches = 0;
	  while (i >= 0) {
	    splice_pos = acceptorj_positions[i];
	    while (nmismatches < nmismatches_right && mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	      debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	      nmismatches++;
	    }
#if 0
	    assert(nmismatches == Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								    segment_left,/*pos5*/splice_pos,/*pos3*/querylength,dibasep,cmetp,plusp));
#endif
	    if (nmismatches <= max_mismatches_allowed) {
	      if (acceptorj_knowni[i] >= 0) {
		prob = 1.0;
		debug4e(printf("Known acceptor for segment at %u, splice_pos %d (%d mismatches), stopi = %d\n",
			       segment_left,splice_pos,nmismatches,splice_pos_start));
	      } else {
		prob = Maxent_hr_acceptor_prob(segment_left + splice_pos,genome_blocks);
		debug4e(printf("Novel acceptor for segment at %u, splice_pos %d (%d mismatches), stopi = %d\n",
			       segment_left,splice_pos,nmismatches,splice_pos_start));
	      }

	      if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob)) {
		ncolordiffs = 0;
		if ((hit = Substring_new_acceptor(acceptorj_knowni[i],/*joffset*/0,splice_pos,nmismatches,ncolordiffs,prob,
						  /*left*/segment_left,query_compress,genome_blocks,snp_blocks,
						  querylength,plusp,sensep,query,segment->chrnum,segment->chroffset,
						  dibasep,cmetp)) != NULL) {
		  debug4e(printf("=> %s acceptor: %f at %d (%d mismatches)\n",
				 plusp == true ? "plus" : "minus",prob,Substring_chimera_pos(hit),nmismatches));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*acceptors)[nmismatches] = List_push((*acceptors)[nmismatches],(void *) hit);
		}
	      }
	    }
	  
	    i--;
	  }


	  /* Splicing originally on minus strand.  Complement.  */
	  sensep = (plusp == true) ? false : true;
	  if (novelsplicingp && segment_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
	    antidonorj_nsites = Genome_antidonor_positions(positions_alloc,knowni_alloc,
							   segment_antidonor_knownpos,segment_antidonor_knowni,
							   genome_blocks,segment_left,
							   splice_pos_start,splice_pos_end+1);
	    antidonorj_positions = positions_alloc;
	    antidonorj_knowni = knowni_alloc;
	  } else {
	    antidonorj_nsites = segment_antidonor_nknown;
	    antidonorj_positions = segment_antidonor_knownpos;
	    antidonorj_knowni = segment_antidonor_knowni;
	  }

	  i = antidonorj_nsites - 1;
	  nmismatches = 0;
	  while (i >= 0) {
	    splice_pos = antidonorj_positions[i];
	    while (nmismatches < nmismatches_right && mismatch_positions[nmismatches] >= splice_pos) { /* Must be >= */
	      debug4e(printf("  mismatch at %d\n",mismatch_positions[nmismatches]));
	      nmismatches++;
	    }
#if 0
	    assert(nmismatches == Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
								    segment_left,/*pos5*/splice_pos,/*pos3*/querylength,dibasep,cmetp,plusp));
#endif
	    if (nmismatches <= max_mismatches_allowed) {
	      if (antidonorj_knowni[i] >= 0) {
		prob = 1.0;
		debug4e(printf("Known antidonor for segmenti at %u, splice_pos %d (%d mismatches), stopi = %d\n",
			       segment_left,splice_pos,nmismatches,splice_pos_start));
	      } else {
		prob = Maxent_hr_antidonor_prob(segment_left + splice_pos,genome_blocks);
		debug4e(printf("Novel antidonor for segmenti at %u, splice_pos %d (%d mismatches), stopi = %d\n",
			       segment_left,splice_pos,nmismatches,splice_pos_start));
	      }

	      debug4e(printf("Checking sufficient_splice_prob for support %d, nmismatches %d, prob %f\n",
			     querylength - splice_pos,nmismatches,prob));
	      if (sufficient_splice_prob_distant(/*support*/querylength - splice_pos,nmismatches,prob)) {
		ncolordiffs = 0;
		if ((hit = Substring_new_donor(antidonorj_knowni[i],/*joffset*/0,splice_pos,nmismatches,ncolordiffs,prob,
					       /*left*/segment_left,query_compress,genome_blocks,snp_blocks,
					       querylength,plusp,sensep,query,segment->chrnum,segment->chroffset,
					       dibasep,cmetp)) != NULL) {
		  debug4e(printf("=> %s antidonor: %f at %d (%d mismatches)\n",
				 plusp == true ? "plus" : "minus",prob,Substring_chimera_pos(hit),nmismatches));
		  debug4e(printf("q: %s\ng: %s\n",queryptr,gbuffer));
		  (*antidonors)[nmismatches] = List_push((*antidonors)[nmismatches],(void *) hit);
		}
	      }
	    }

	    i--;
	  }
	}

      }
    }
  }

  return;
}




static List_T
find_terminals (struct Segment_T *plus_segments, int plus_nsegments,
		struct Segment_T *minus_segments, int minus_nsegments,
		UINT4 *genome_blocks, UINT4 *snp_blocks,
		char *queryuc_ptr, /* for debugging */ char *queryrc,
		Floors_T floors, int querylength, int query_lastpos,
		Compress_T query_compress_fwd, Compress_T query_compress_rev,
		int max_mismatches_allowed, int terminal_penalty, int max_terminal_length,
		bool dibasep, bool cmetp) {
#ifdef DEBUG4T
  char gbuffer[MAX_QUERYLENGTH];
#endif
  List_T terminals = (List_T) NULL;
  Segment_T segment;
  Stage3_T hit;
  Genomicpos_T segment_left;
  int nmismatches_left, nmismatches_right, nmismatches;

  int mismatch_positions[MAX_QUERYLENGTH+1], colordiffs[MAX_QUERYLENGTH];
  int ncolordiffs = 0;
  int *floors_from_neg3, *floors_to_pos3;

#ifdef DEBUG4T
  int i;
#endif

  debug(printf("identify_terminals: Checking up to %d mismatches\n",max_mismatches_allowed));

  if (floors == NULL) {
    return (List_T) NULL;

  } else {
    floors_from_neg3 = floors->scorefrom[-3];
    floors_to_pos3 = floors->scoreto[query_lastpos+3];

    if (max_terminal_length > querylength/2) {
      max_terminal_length = querylength/2;
    }
  }

  if (plus_nsegments > 0) {
    for (segment = plus_segments; segment < &(plus_segments[plus_nsegments]); segment++) {
      segment_left = segment->diagonal - querylength; /* FORMULA: Corresponds to querypos 0 */
      debug4t(printf("identify_terminals_plus: Checking up to %d mismatches at diagonal %u (querypos %d..%d) - querylength %d = %u\n",
		     max_mismatches_allowed,segment->diagonal,segment->querypos5,segment->querypos3,querylength,segment_left));

#if 0
      debug4t(
	      Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	      printf("genome 0..: %s\n",gbuffer);
	      printf("query  0..: %s\n",queryuc_ptr);
	      );
#endif

      if (segment->floor_left > max_mismatches_allowed) {
	debug4t(printf("Not checking left because floor_left %d > max_mismatches_allowed %d\n",
		       segment->floor_left,max_mismatches_allowed));
      } else {
	/* Check from left */
	debug4t(printf("Checking left because floor_left %d <= max_mismatches_allowed %d\n",
		       segment->floor_left,max_mismatches_allowed));

	nmismatches_left = Genome_mismatches_left(mismatch_positions,colordiffs,max_mismatches_allowed,
						  /*query*/queryuc_ptr,query_compress_fwd,genome_blocks,snp_blocks,
						  /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						  dibasep,cmetp,/*plusp*/true);

	debug4t(
		printf("%d mismatches on left at:",nmismatches_left);
		for (i = 0; i <= nmismatches_left; i++) {
		  printf(" %d",mismatch_positions[i]);
		}
		printf("\n");
		);


	if (nmismatches_left == 0 || nmismatches_left <= max_mismatches_allowed || 
	    mismatch_positions[nmismatches_left-1] > querylength - max_terminal_length) {  /* was querylength/2 */
	  debug4t(printf(" => Terminal at left"));
	  if ((hit = Stage3_new_terminal(/*querystart*/0,/*queryend*//*truncate_pos_left*/querylength,
					 nmismatches_left,/*ncolordiffs*/0,/*left*/segment_left,
					 query_compress_fwd,genome_blocks,snp_blocks,
					 querylength,/*plusp*/true,/*start_endtype*/END,/*end_endtype*/TERM,
					 /*query*/queryuc_ptr,segment->chrnum,segment->chroffset,
					 terminal_penalty,max_mismatches_allowed,
					 dibasep,cmetp)) != NULL) {
	    debug4t(printf(" => yes"));
	    terminals = List_push(terminals,(void *) hit);
	  }
	  debug4t(printf("\n"));
	}
      }

      if (segment->floor_right > max_mismatches_allowed) {
	debug4t(printf("Not checking right because floor_right %d > max_mismatches_allowed %d\n",
		       segment->floor_right,max_mismatches_allowed));
      } else {
	/* Check from right */
	debug4t(printf("Checking right because floor_right %d <= max_mismatches_allowed %d\n",
		       segment->floor_right,max_mismatches_allowed));
	nmismatches_right = Genome_mismatches_right(mismatch_positions,colordiffs,max_mismatches_allowed,
						    /*query*/queryuc_ptr,/*query_compress*/query_compress_fwd,genome_blocks,snp_blocks,
						    /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						    dibasep,cmetp,/*plusp*/true);

	debug4t(
		printf("%d mismatches on right at:",nmismatches_right);
		for (i = 0; i <= nmismatches_right; i++) {
		  printf(" %d",mismatch_positions[i]);
		}
		printf("\n");
		);
      
	if (nmismatches_right == 0 || nmismatches_right <= max_mismatches_allowed ||
	    mismatch_positions[nmismatches_right-1] < max_terminal_length) {   /* was querylength/2 */
	  debug4t(printf(" => Terminal at right"));
	  if ((hit = Stage3_new_terminal(/*querystart*//*truncate_pos_right*/0,/*queryend*/querylength,
					 nmismatches_right,/*ncolordiffs*/0,/*left*/segment_left,
					 query_compress_fwd,genome_blocks,snp_blocks,
					 querylength,/*plusp*/true,/*start_endtype*/TERM,/*end_endtype*/END,
					 /*query*/queryuc_ptr,segment->chrnum,segment->chroffset,
					 terminal_penalty,max_mismatches_allowed,
					 dibasep,cmetp)) != NULL) {
	    debug4t(printf(" => yes"));
	    terminals = List_push(terminals,(void *) hit);
	  }
	  debug4t(printf("\n"));
	}
      }
    }
  }

  if (minus_nsegments > 0) {
    for (segment = minus_segments; segment < &(minus_segments[minus_nsegments]); segment++) {
      segment_left = segment->diagonal - querylength;
      debug4t(printf("identify_terminals_minus: Getting genome at diagonal %u (querypos %d..%d) + 12 - querylength %d = %u\n",
		     segment->diagonal,segment->querypos5,segment->querypos3,querylength,segment_left));

#if 0
      debug4t(
	      Genome_fill_buffer_simple(genome,segment_left,querylength,gbuffer);
	      printf("genome   0..: %s\n",gbuffer);
	      printf("query.rc 0..: %s\n",queryrc);
	      );
#endif

      if (segment->floor_left > max_mismatches_allowed) {
	debug4t(printf("Not checking left because floor_left %d > max_mismatches_allowed %d\n",
		       segment->floor_left,max_mismatches_allowed));
      } else {
	/* Check from left */
	debug4t(printf("Not checking left because floor_left %d <= max_mismatches_allowed %d\n",
		       segment->floor_left,max_mismatches_allowed));
	nmismatches_left = Genome_mismatches_left(mismatch_positions,colordiffs,max_mismatches_allowed,
						  /*query*/queryuc_ptr,/*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
						  /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						  dibasep,cmetp,/*plusp*/false);

	debug4t(
		printf("%d mismatches on left at:",nmismatches_left);
		for (i = 0; i <= nmismatches_left; i++) {
		  printf(" %d",mismatch_positions[i]);
		}
		printf("\n");
		);

	if (nmismatches_left == 0 || nmismatches_left <= max_mismatches_allowed || 
	    mismatch_positions[nmismatches_left-1] > querylength - max_terminal_length) {  /* was querylength/2 */
	  debug4t(printf(" => Terminal at left"));
	  if ((hit = Stage3_new_terminal(/*querystart*//*querylength-truncate_pos_left*/0,/*queryend*/querylength,
					 nmismatches_left,/*ncolordiffs*/0,/*left*/segment_left,
					 query_compress_rev,genome_blocks,snp_blocks,
					 querylength,/*plusp*/false,/*start_endtype*/TERM,/*end_endtype*/END,
					 /*query*/queryuc_ptr,segment->chrnum,segment->chroffset,
					 terminal_penalty,max_mismatches_allowed,
					 dibasep,cmetp)) != NULL) {
	    debug4t(printf(" => yes"));
	    terminals = List_push(terminals,(void *) hit);
	  }
	  debug4t(printf("\n"));
	}
      }

      if (segment->floor_right > max_mismatches_allowed) {
	debug4t(printf("Not checking right because floor_right %d > max_mismatches_allowed %d\n",
		       segment->floor_right,max_mismatches_allowed));
      } else {
	/* Check from right */
	debug4t(printf("Checking right because floor_right %d <= max_mismatches_allowed %d\n",
		       segment->floor_right,max_mismatches_allowed));
	nmismatches_right = Genome_mismatches_right(mismatch_positions,colordiffs,max_mismatches_allowed,
						    /*query*/queryuc_ptr,/*query_compress*/query_compress_rev,genome_blocks,snp_blocks,
						    /*left*/segment_left,/*pos5*/0,/*pos3*/querylength,
						    dibasep,cmetp,/*plusp*/false);

	debug4t(
		printf("%d mismatches on right at:",nmismatches_right);
		for (i = 0; i <= nmismatches_right; i++) {
		  printf(" %d",mismatch_positions[i]);
		}
		printf("\n");
		);

	if (nmismatches_right == 0 || nmismatches_right <= max_mismatches_allowed ||
	    mismatch_positions[nmismatches_right-1] < max_terminal_length) {   /* was querylength/2 */
	  debug4t(printf(" => Terminal at right"));
	  if ((hit = Stage3_new_terminal(/*querystart*/0,/*queryend*//*querylength-truncate_pos_right*/querylength,
					 nmismatches_right,/*ncolordiffs*/0,/*left*/segment_left,
					 query_compress_rev,genome_blocks,snp_blocks,
					 querylength,/*plusp*/false,/*start_endtype*/END,/*end_endtype*/TERM,
					 /*query*/queryuc_ptr,segment->chrnum,segment->chroffset,
					 terminal_penalty,max_mismatches_allowed,
					 dibasep,cmetp)) != NULL) {
	    debug4t(printf(" => yes"));
	    terminals = List_push(terminals,(void *) hit);
	  }
	  debug4t(printf("\n"));
	}
      }
    }
  }

  debug4t(printf("Total number of terminals: %d\n",List_length(terminals)));
  return terminals;
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
find_splicepairs_distant (int *found_score, int *nsplicepairs, List_T distantsplicing_orig,
			  List_T *donors_plus, List_T *antidonors_plus,
			  List_T *acceptors_plus, List_T *antiacceptors_plus,
			  List_T *donors_minus, List_T *antidonors_minus,
			  List_T *acceptors_minus, List_T *antiacceptors_minus,
			  Genomicpos_T shortsplicedist, int localsplicing_penalty, int distantsplicing_penalty,
			  int min_distantsplicing_end_matches, double min_distantsplicing_identity,
			  int querylength, int nmismatches_allowed, int maxchimerapaths, bool first_read_p) {
  List_T distantsplicing = NULL, p, q, qsave;
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
	/* Generate all pairs at this splice_pos */
	qsave = q;
	while (p != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really a continuous match. */
	    } else {
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
	      debug4ld(printf("1-2. Pushing a candidate at splice_pos %d (%d..%d), donor %u to acceptor %u.  shortdistancep = %d\n",
			      pos,min_endlength_1,querylength-min_endlength_2,
			      Substring_genomicstart(donor),Substring_genomicstart(acceptor),shortdistancep));
#if 0
	      if (shortdistancep) {
		*localsplicing = List_push(*localsplicing,
					   (void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								      donor,acceptor,distance,
								      shortdistancep,splicing_penalty,querylength,
								      /*ambi_left*/NULL,/*ambi_right*/NULL,
								      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								      /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								      /*sensedir*/SENSE_FORWARD));
	      } else {
#endif
		distantsplicing = List_push(distantsplicing,
					    (void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								       donor,acceptor,distance,
								       shortdistancep,splicing_penalty,querylength,
								       /*ambi_left*/NULL,/*ambi_right*/NULL,
								       /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								       /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								       /*sensedir*/SENSE_FORWARD));
#if 0
	      }
#endif
	      (*nsplicepairs)++;

	    }
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
	while (p != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really a continuous match. */
	    } else {
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
	      debug4ld(printf("3-4. Pushing a candidate at splice_pos %d (%d..%d), donor %u to acceptor %u.  shortdistancep = %d.\n",
			      pos,min_endlength_1,querylength-min_endlength_2,
			      Substring_genomicstart(donor),Substring_genomicstart(acceptor),shortdistancep));
#if 0
	      if (shortdistancep) {
		*localsplicing = List_push(*localsplicing,
					   (void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								      donor,acceptor,distance,
								      shortdistancep,distantsplicing_penalty,querylength,
								      /*ambi_left*/NULL,/*ambi_right*/NULL,
								      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								      /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								      /*sensedir*/SENSE_FORWARD));
	      } else {
#endif
		distantsplicing = List_push(distantsplicing,
					    (void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								       donor,acceptor,distance,
								       shortdistancep,splicing_penalty,querylength,
								       /*ambi_left*/NULL,/*ambi_right*/NULL,
								       /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								       /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								       /*sensedir*/SENSE_FORWARD));
#if 0
	      }
#endif
	      (*nsplicepairs)++;
	    }
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
	while (p != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really an continuous match. */
	    } else {
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

	      debug4ld(printf("5-6. Pushing a candidate at splice_pos %d (%d..%d), donor %u to acceptor %u.  shortdistancep = %d\n",
			      pos,min_endlength_2,querylength-min_endlength_1,
			      Substring_genomicstart(donor),Substring_genomicstart(acceptor),shortdistancep));
#if 0
	      if (shortdistancep) {
		*localsplicing = List_push(*localsplicing,
					   (void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								      donor,acceptor,distance,
								      shortdistancep,distantsplicing_penalty,querylength,
								      /*ambi_left*/NULL,/*ambi_right*/NULL,
								      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								      /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								      /*sensedir*/SENSE_ANTI));
	      } else {
#endif
		distantsplicing = List_push(distantsplicing,
					    (void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								       donor,acceptor,distance,
								       shortdistancep,splicing_penalty,querylength,
								       /*ambi_left*/NULL,/*ambi_right*/NULL,
								       /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								       /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								       /*sensedir*/SENSE_ANTI));
#if 0
	      }
#endif
	      (*nsplicepairs)++;
	    }
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

	while (p != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_genomicstart(acceptor) == Substring_genomicstart(donor)) {
	      /* Skip.  Really a continuous match. */
	    } else {
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
	      debug4ld(printf("7-8. Pushing a candidate at splice_pos %d (%d..%d), donor %u to acceptor %u.  shortdistancep = %d.\n",
			      pos,min_endlength_2,querylength-min_endlength_1,
			      Substring_genomicstart(donor),Substring_genomicstart(acceptor),shortdistancep));
#if 0
	      if (shortdistancep) {
		*localsplicing = List_push(*localsplicing,
					   (void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								      donor,acceptor,distance,
								      shortdistancep,distantsplicing_penalty,querylength,
								      /*ambi_left*/NULL,/*ambi_right*/NULL,
								      /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								      /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								      /*sensedir*/SENSE_ANTI));
	      } else {
#endif
		distantsplicing = List_push(distantsplicing,
					    (void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								       donor,acceptor,distance,
								       shortdistancep,distantsplicing_penalty,querylength,
								       /*ambi_left*/NULL,/*ambi_right*/NULL,
								       /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								       /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								       /*sensedir*/SENSE_ANTI));
#if 0
	      }
#endif
	      (*nsplicepairs)++;
	    }
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
	while (p != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) - pos) > (Substring_genomicstart(donor) + pos)) {
	      distance = (Substring_genomicstart(acceptor) - pos) - (Substring_genomicstart(donor) + pos);
	    } else {
	      distance = (Substring_genomicstart(donor) + pos) - (Substring_genomicstart(acceptor) - pos);
	    }
	    debug4ld(printf("1-4. Pushing a candidate at splice_pos %d (%d..%d), donor %u to acceptor %u.  Different strands, so not shortdistance.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								   donor,acceptor,distance,
								   /*shortdistancep*/false,distantsplicing_penalty,querylength,
								   /*ambi_left*/NULL,/*ambi_right*/NULL,
								   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								   /*sensedir*/SENSE_FORWARD));
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
	while (p != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if (Substring_genomicstart(acceptor) > Substring_genomicstart(donor)) {
	      distance = (Substring_genomicstart(acceptor) + pos) - (Substring_genomicstart(donor) - pos);
	    } else {
	      distance = (Substring_genomicstart(donor) - pos) - (Substring_genomicstart(acceptor) + pos);
	    }
	    debug4ld(printf("3-2. Pushing a candidate at splice_pos %d (%d..%d), donor %u to acceptor %u.  Different strands so not shortdistance.\n",
			    pos,min_endlength_1,querylength-min_endlength_2,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								   donor,acceptor,distance,
								   /*shortdistancep*/false,distantsplicing_penalty,querylength,
								   /*ambi_left*/NULL,/*ambi_right*/NULL,
								   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								   /*sensedir*/SENSE_FORWARD));
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
	while (p != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) - pos) > (Substring_genomicstart(donor) + pos)) {
	      distance = (Substring_genomicstart(acceptor) - pos) - (Substring_genomicstart(donor) + pos);
	    } else {
	      distance = (Substring_genomicstart(donor) + pos) - (Substring_genomicstart(acceptor) - pos);
	    }
	    debug4ld(printf("5-8. Pushing a candidate at splice_pos %d (%d..%d), donor %u to acceptor %u.  Different strands so not shortdistance.\n",
			    pos,min_endlength_2,querylength-min_endlength_1,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								   donor,acceptor,distance,
								   /*shortdistancep*/false,distantsplicing_penalty,querylength,
								   /*ambi_left*/NULL,/*ambi_right*/NULL,
								   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								   /*sensedir*/SENSE_ANTI));
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
	while (p != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(p->first) == pos) {
	  donor = (Substring_T) p->first;
	  debug4ld(printf("donor at %u, pos %d\n",Substring_genomicstart(donor),pos));
	  q = qsave;
	  while (q != NULL && *nsplicepairs < maxchimerapaths && Substring_chimera_pos(q->first) == pos) {
	    acceptor = (Substring_T) q->first;
	    debug4ld(printf("acceptor at %u, pos %d\n",Substring_genomicstart(acceptor),pos));
	    if (Substring_chrnum(donor) != Substring_chrnum(acceptor)) {
	      distance = 0U;
	    } else if ((Substring_genomicstart(acceptor) + pos) > (Substring_genomicstart(donor) - pos)) {
	      distance = (Substring_genomicstart(acceptor) + pos) - (Substring_genomicstart(donor) - pos);
	    } else {
	      distance = (Substring_genomicstart(donor) - pos) - (Substring_genomicstart(acceptor) + pos);
	    }
	    debug4ld(printf("7-6. Pushing a candidate at splice_pos %d (%d..%d), donor %u to acceptor %u.  Different strands so not shortdistance.\n",
			    pos,min_endlength_2,querylength-min_endlength_1,
			    Substring_genomicstart(donor),Substring_genomicstart(acceptor)));
	    distantsplicing = List_push(distantsplicing,
					(void *) Stage3_new_splice(&(*found_score),nmismatches1,nmismatches2,
								   donor,acceptor,distance,
								   /*shortdistancep*/false,distantsplicing_penalty,querylength,
								   /*ambi_left*/NULL,/*ambi_right*/NULL,
								   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								   /*copy_donor_p*/true,/*copy_acceptor_p*/true,first_read_p,
								   /*sensedir*/SENSE_ANTI));
	    (*nsplicepairs)++;
	    q = q->rest;
	  }
	  p = p->rest;
	}
      }
    }
  }

  if (*nsplicepairs >= maxchimerapaths) {
    stage3list_gc(&distantsplicing);
    return distantsplicing_orig;
  } else {
    return List_append(distantsplicing_orig,distantsplicing);
  }
}


static List_T
find_splicepairs_short_overlaps (int *found_score, List_T hits,
				 List_T *donors_plus, List_T *antidonors_plus,
				 List_T *acceptors_plus, List_T *antiacceptors_plus,
				 List_T *donors_minus, List_T *antidonors_minus,
				 List_T *acceptors_minus, List_T *antiacceptors_minus,
				 Compress_T query_compress_fwd, Compress_T query_compress_rev,
				 UINT4 *genome_blocks, UINT4 *snp_blocks,
				 char *queryuc_ptr, char *queryrc, int min_shortend,
				 Genomicpos_T *splicesites, Splicetype_T *splicetypes,
				 UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int nsplicesites, int *nsplicepartners_skip,
				 unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
				 unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
				 List_T *splicestrings, Genomicpos_T shortsplicedist, int localsplicing_penalty,
				 int min_localsplicing_end_matches, int max_mismatches_allowed,
				 int querylength, bool first_read_p, bool dibasep, bool cmetp) {
  List_T p;
  Substring_T donor, acceptor;
  Intlist_T splicesites_i;
  Intlist_T nmismatches_list;
  int nmismatches, nmismatches_shortend, nmisses_allowed, ncolordiffs, support, endlength;
#ifdef DEBUG4H
  Genomicpos_T chroffset, leftbound, rightbound;
#endif
  Genomicpos_T bestleft, origleft;
  int i;
  int bestj = 0;


  debug(printf("Starting find_splices_short_overlaps\n"));
  debug(
	for (nmismatches = 0; nmismatches <= max_mismatches_allowed; nmismatches++) {
	  printf("At %d nmismatches: +donors/acceptors %d/%d, +antidonors/antiacceptors %d/%d, -donors/acceptors %d/%d, -antidonors/antiacceptors %d/%d\n",
		 nmismatches,
		 List_length(donors_plus[nmismatches]),
		 List_length(acceptors_plus[nmismatches]),
		 List_length(antidonors_plus[nmismatches]),
		 List_length(antiacceptors_plus[nmismatches]),
		 List_length(donors_minus[nmismatches]),
		 List_length(acceptors_minus[nmismatches]),
		 List_length(antidonors_minus[nmismatches]),
		 List_length(antiacceptors_minus[nmismatches]));
	});

  /* Donors and antiacceptors => Want chimera_pos to be at end */
  /* Acceptors and antidonors => Want chimera_pos to be at beginning */

  /* Don't want to end when first set of hits found */
  for (nmismatches = 0; /* hits == NULL && */ nmismatches <= max_mismatches_allowed;
       nmismatches++) {
    nmisses_allowed = max_mismatches_allowed - nmismatches;

    /* End 1 */
    for (p = donors_plus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      support = Substring_chimera_pos(donor);
      endlength = querylength - support;
      
#ifdef DEBUG4H
      chroffset = Substring_chroffset(donor);
      leftbound = Substring_alignend_trim(donor) + 1U;
#endif
      debug4h(printf("End 1: short-overlap donor_plus: #%d:%u, endlength %d\n",
		     Substring_chrnum(donor),leftbound-1-chroffset,endlength));

      if (endlength <= min_localsplicing_end_matches) {
	debug4h(printf("End 1: short-overlap donor_plus: #%d:%u (%d mismatches) => searching right\n",
		       Substring_chrnum(donor),leftbound-1-chroffset,Substring_nmismatches_whole(donor)));

	if ((i = Substring_splicesites_i(donor)) >= 0) {
	  origleft = Substring_genomicstart(donor);
	  if ((splicesites_i = 
	       Splicetrie_find_right(&nmismatches_shortend,&nmismatches_list,i,
				     origleft,/*pos5*/support,/*pos3*/querylength,
				     splicesites,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
				     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				     trieoffsets_max,triecontents_max,nsplicepartners_max,
				     splicetypes,splicestrings,query_compress_fwd,
				     genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,
				     nmisses_allowed,dibasep,cmetp,/*plusp*/true)) != NULL) {
	    
	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      debug4h(printf("End 1: short-overlap donor_plus: Successful ambiguous from donor #%d\n",
			     Substring_splicesites_i(donor)));
	      hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
							       donor,/*acceptor*/NULL,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,
							       /*ambi_left*/NULL,/*ambi_right*/splicesites_i,
							       /*amb_nmismatches_left*/NULL,nmismatches_list,
							       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
							       /*sensedir*/SENSE_FORWARD));
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - support;
	      if ((acceptor = Substring_new_acceptor(/*splicesites_i*/bestj,/*joffset*/0,Substring_chimera_pos(donor),nmismatches_shortend,/*ncolordiffs*/0,
						     /*prob*/2.0,/*left*/bestleft,query_compress_fwd,genome_blocks,snp_blocks,
						     querylength,/*plusp*/true,/*sensep*/true,
						     /*query*/queryuc_ptr,Substring_chrnum(donor),Substring_chroffset(donor),
						     dibasep,cmetp)) != NULL) {
		debug4h(printf("End 1: short-overlap donor_plus: Successful splice from donor #%d to acceptor #%d\n",
			       Substring_splicesites_i(donor),Substring_splicesites_i(acceptor)));
		hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								 donor,acceptor,/*distance*/bestleft-origleft,
								 /*shortdistancep*/true,localsplicing_penalty,querylength,
								 /*ambi_left*/NULL,/*ambi_right*/NULL,
								 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								 /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								 /*sensedir*/SENSE_FORWARD));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 2 */
    for (p = acceptors_plus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      endlength = Substring_chimera_pos(acceptor);
      support = querylength - endlength;

#ifdef DEBUG4H
      chroffset = Substring_chroffset(acceptor);
      rightbound = Substring_alignstart_trim(acceptor);
#endif

      debug4h(printf("End 2: short-overlap acceptor_plus: #%d:%u, endlength %d\n",
		     Substring_chrnum(acceptor),rightbound+1-chroffset,endlength));

      if (endlength <= min_localsplicing_end_matches) {
	debug4h(printf("End 2: short-overlap acceptor_plus: #%d:%u (%d mismatches) => searching left\n",
		       Substring_chrnum(acceptor),rightbound+1-chroffset,Substring_nmismatches_whole(acceptor)));

	if ((i = Substring_splicesites_i(acceptor)) >= 0) {
	  origleft = Substring_genomicstart(acceptor);
	  if ((splicesites_i =
	       Splicetrie_find_left(&nmismatches_shortend,&nmismatches_list,i,
				    origleft,/*pos5*/0,/*pos3*/endlength,
				    splicesites,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
				    trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				    trieoffsets_max,triecontents_max,nsplicepartners_max,
				    splicetypes,splicestrings,query_compress_fwd,
				    genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,
				    nmisses_allowed,dibasep,cmetp,/*plusp*/true)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      debug4h(printf("End 2: short-overlap acceptor_plus: Successful ambiguous from acceptor #%d\n",
			     Substring_splicesites_i(acceptor)));
	      hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
							       /*donor*/NULL,acceptor,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,
							       /*ambi_left*/splicesites_i,/*ambi_right*/NULL,
							       nmismatches_list,/*amb_nmismatches_right*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
							       /*sensedir*/SENSE_FORWARD));

	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - endlength;
	      if ((donor = Substring_new_donor(/*splicesites_i*/bestj,/*joffset*/0,Substring_chimera_pos(acceptor),nmismatches_shortend,/*ncolordiffs*/0,
					       /*prob*/2.0,/*left*/bestleft,query_compress_fwd,genome_blocks,snp_blocks,
					       querylength,/*plusp*/true,/*sensep*/true,
					       queryuc_ptr,Substring_chrnum(acceptor),Substring_chroffset(acceptor),
					       dibasep,cmetp)) != NULL) {
		debug4h(printf("End 2: short-overlap acceptor_plus: Successful splice from acceptor #%d to donor #%d\n",
			       Substring_splicesites_i(acceptor),Substring_splicesites_i(donor)));

		hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								 donor,acceptor,/*distance*/origleft-bestleft,
								 /*shortdistancep*/true,localsplicing_penalty,querylength,
								 /*ambi_left*/NULL,/*ambi_right*/NULL,
								 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								 /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								 /*sensedir*/SENSE_FORWARD));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 3 */
    for (p = donors_minus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      support = Substring_chimera_pos(donor);
      endlength = querylength - support;

#ifdef DEBUG4H
      chroffset = Substring_chroffset(donor);
      rightbound = Substring_alignend_trim(donor);
#endif

      debug4h(printf("End 3: short-overlap donor_minus: #%d:%u, endlength %d\n",
		     Substring_chrnum(donor),rightbound+1U-chroffset,endlength));

      if (endlength <= min_localsplicing_end_matches) {
	debug4h(printf("End 3: short-overlap donor_minus: #%d:%u (%d mismatches) => searching left\n",
		       Substring_chrnum(donor),rightbound+1U-chroffset,Substring_nmismatches_whole(donor)));

	if ((i = Substring_splicesites_i(donor)) >= 0) {
	  origleft = Substring_genomicend(donor);
	  if ((splicesites_i =
	       Splicetrie_find_left(&nmismatches_shortend,&nmismatches_list,i,
				    origleft,/*pos5*/0,/*pos3*/endlength,
				    splicesites,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
				    trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				    trieoffsets_max,triecontents_max,nsplicepartners_max,
				    splicetypes,splicestrings,query_compress_rev,
				    genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*queryptr*/queryrc,
				    nmisses_allowed,dibasep,cmetp,/*plusp*/false)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      debug4h(printf("End 3: short-overlap donor_minus: Successful ambiguous from donor #%d\n",
			     Substring_splicesites_i(donor)));
	      hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
							       donor,/*acceptor*/NULL,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,
							       /*ambi_left*/splicesites_i,/*ambi_right*/NULL,
							       nmismatches_list,/*amb_nmismatches_right*/NULL,
							       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
							       /*sensedir*/SENSE_FORWARD));
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - endlength;
	      if ((acceptor = Substring_new_acceptor(/*splicesites_i*/bestj,/*joffset*/0,
						     querylength-Substring_chimera_pos(donor),nmismatches_shortend,/*ncolordiffs*/0,
						     /*prob*/2.0,/*left*/bestleft,query_compress_rev,genome_blocks,snp_blocks,
						     querylength,/*plusp*/false,/*sensep*/true,
						     queryuc_ptr,Substring_chrnum(donor),Substring_chroffset(donor),
						     dibasep,cmetp)) != NULL) {
		debug4h(printf("End 3: short-overlap donor_minus: Successful splice from donor #%d to acceptor #%d\n",
			       Substring_splicesites_i(donor),Substring_splicesites_i(acceptor)));
		hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								 donor,acceptor,/*distance*/origleft-bestleft,
								 /*shortdistancep*/true,localsplicing_penalty,querylength,
								 /*ambi_left*/NULL,/*ambi_right*/NULL,
								 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								 /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								 /*sensedir*/SENSE_FORWARD));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 4 */
    for (p = acceptors_minus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      endlength = Substring_chimera_pos(acceptor);
      support = querylength - endlength;

#ifdef DEBUG4H
      chroffset = Substring_chroffset(acceptor);
      leftbound = Substring_alignstart_trim(acceptor) + 1U;
#endif

      debug4h(printf("End 4: short-overlap acceptor_minus: #%d:%u, endlength %d\n",
		     Substring_chrnum(acceptor),leftbound-1U-chroffset,endlength));

      if (endlength <= min_localsplicing_end_matches) {
	debug4h(printf("End 4: short-overlap acceptor_minus: #%d:%u (%d mismatches) => searching right\n",
		       Substring_chrnum(acceptor),leftbound-1U-chroffset,Substring_nmismatches_whole(acceptor)));

	if ((i = Substring_splicesites_i(acceptor)) >= 0) {
	  origleft = Substring_genomicend(acceptor);
	  if ((splicesites_i =
	       Splicetrie_find_right(&nmismatches_shortend,&nmismatches_list,i,
				     origleft,/*pos5*/support,/*pos3*/querylength,
				     splicesites,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
				     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				     trieoffsets_max,triecontents_max,nsplicepartners_max,
				     splicetypes,splicestrings,query_compress_rev,
				     genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*queryptr*/queryrc,
				     nmisses_allowed,dibasep,cmetp,/*plusp*/false)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      debug4h(printf("End 4: short-overlap acceptor_minus: Successful ambiguous from acceptor #%d\n",
			     Substring_splicesites_i(acceptor)));
	      hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
							       /*donor*/NULL,acceptor,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,
							       /*ambi_left*/NULL,/*ambi_right*/splicesites_i,
							       /*amb_nmismatches_left*/NULL,nmismatches_list,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
							       /*sensedir*/SENSE_FORWARD));
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - support;
	      if ((donor = Substring_new_donor(/*splicesites_i*/bestj,/*joffset*/0,
					       querylength-Substring_chimera_pos(acceptor),nmismatches_shortend,/*ncolordiffs*/0,
					       /*prob*/2.0,/*left*/bestleft,query_compress_rev,genome_blocks,snp_blocks,
					       querylength,/*plusp*/false,/*sensep*/true,
					       queryuc_ptr,Substring_chrnum(acceptor),Substring_chroffset(acceptor),
					       dibasep,cmetp)) != NULL) {
		debug4h(printf("End 4: short-overlap acceptor_minus: Successful splice from acceptor #%d to #%d\n",
			       Substring_splicesites_i(acceptor),Substring_splicesites_i(donor)));
		hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								 donor,acceptor,/*distance*/bestleft-origleft,
								 /*shortdistancep*/true,localsplicing_penalty,querylength,
								 /*ambi_left*/NULL,/*ambi_right*/NULL,
								 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								 /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								 /*sensedir*/SENSE_FORWARD));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 5 */
    for (p = antidonors_plus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      endlength = Substring_chimera_pos(donor);
      support = querylength - endlength;

#ifdef DEBUG4H
      chroffset = Substring_chroffset(donor);
      rightbound = Substring_alignstart_trim(donor);
#endif

      debug4h(printf("End 5: short-overlap antidonor_plus: #%d:%u, endlength %d\n",
		     Substring_chrnum(donor),rightbound+1U-chroffset,endlength));

      if (endlength <= min_localsplicing_end_matches) {
	debug4h(printf("End 5: short-overlap antidonor_plus: #%d:%u (%d mismatches) => searching left\n",
		       Substring_chrnum(donor),rightbound+1U-chroffset,Substring_nmismatches_whole(donor)));

	if ((i = Substring_splicesites_i(donor)) >= 0) {
	  origleft = Substring_genomicstart(donor);
	  if ((splicesites_i =
	       Splicetrie_find_left(&nmismatches_shortend,&nmismatches_list,i,
				    origleft,/*pos5*/0,/*pos3*/endlength,
				    splicesites,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
				    trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				    trieoffsets_max,triecontents_max,nsplicepartners_max,
				    splicetypes,splicestrings,query_compress_fwd,
				    genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,
				    nmisses_allowed,dibasep,cmetp,/*plusp*/true)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      debug4h(printf("End 5: short-overlap antidonor_plus: Successful ambiguous from antidonor #%d\n",
			     Substring_splicesites_i(donor)));
	      hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
							       donor,/*acceptor*/NULL,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,
							       /*ambi_left*/splicesites_i,/*ambi_right*/NULL,
							       nmismatches_list,/*amb_nmismatches_right*/NULL,
							       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
							       /*sensedir*/SENSE_ANTI));
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - endlength;
	      if ((acceptor = Substring_new_acceptor(/*splicesites_i*/bestj,/*joffset*/0,Substring_chimera_pos(donor),nmismatches_shortend,/*ncolordiffs*/0,
						     /*prob*/2.0,/*left*/bestleft,query_compress_fwd,genome_blocks,snp_blocks,
						     querylength,/*plusp*/true,/*sensep*/false,
						     queryuc_ptr,Substring_chrnum(donor),Substring_chroffset(donor),
						     dibasep,cmetp)) != NULL) {
		debug4h(printf("End 5: short-overlap antidonor_plus: Successful splice from antidonor #%d to antiacceptor #%d\n",
			       Substring_splicesites_i(donor),Substring_splicesites_i(acceptor)));
		hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								 donor,acceptor,/*distance*/origleft-bestleft,
								 /*shortdistancep*/true,localsplicing_penalty,querylength,
								 /*ambi_left*/NULL,/*ambi_right*/NULL,
								 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								 /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								 /*sensedir*/SENSE_ANTI));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 6 */
    for (p = antiacceptors_plus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      support = Substring_chimera_pos(acceptor);
      endlength = querylength - support;

#ifdef DEBUG4H
      chroffset = Substring_chroffset(acceptor);
      leftbound = Substring_alignend_trim(acceptor) + 1U;
#endif

      debug4h(printf("End 6: short-overlap antiacceptor_plus: #%d:%u, endlength %d\n",
		     Substring_chrnum(acceptor),leftbound-1U-chroffset,endlength));

      if (endlength <= min_localsplicing_end_matches) {
	debug4h(printf("End 6: short-overlap antiacceptor_plus: #%d:%u (%d mismatches) => searching right\n",
		       Substring_chrnum(acceptor),leftbound-1U-chroffset,Substring_nmismatches_whole(acceptor)));

	if ((i = Substring_splicesites_i(acceptor)) >= 0) {
	  origleft = Substring_genomicstart(acceptor);
	  if ((splicesites_i =
	       Splicetrie_find_right(&nmismatches_shortend,&nmismatches_list,i,
				     origleft,/*pos5*/support,/*pos3*/querylength,
				     splicesites,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
				     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				     trieoffsets_max,triecontents_max,nsplicepartners_max,
				     splicetypes,splicestrings,query_compress_fwd,
				     genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*query*/queryuc_ptr,
				     nmisses_allowed,dibasep,cmetp,/*plusp*/true)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      debug4h(printf("End 6: short-overlap antiacceptor_plus: Successful ambiguous from antiacceptor #%d\n",
			     Substring_splicesites_i(acceptor)));
	      hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
							       /*donor*/NULL,acceptor,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,
							       /*ambi_left*/NULL,/*ambi_right*/splicesites_i,
							       /*amb_nmismatches_left*/NULL,nmismatches_list,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
							       /*sensedir*/SENSE_ANTI));
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - support;
	      if ((donor = Substring_new_donor(/*splicesites_i*/bestj,/*joffset*/0,Substring_chimera_pos(acceptor),nmismatches_shortend,/*ncolordiffs*/0,
					       /*prob*/2.0,/*left*/bestleft,query_compress_fwd,genome_blocks,snp_blocks,
					       querylength,/*plusp*/true,/*sensep*/false,
					       queryuc_ptr,Substring_chrnum(acceptor),Substring_chroffset(acceptor),
					       dibasep,cmetp)) != NULL) {
		debug4h(printf("End 6: short-overlap antiacceptor_plus: Successful splice from antiacceptor #%d to antidonor #%d\n",
			       Substring_splicesites_i(acceptor),Substring_splicesites_i(donor)));
		hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								 donor,acceptor,/*distance*/bestleft-origleft,
								 /*shortdistancep*/true,localsplicing_penalty,querylength,
								 /*ambi_left*/NULL,/*ambi_right*/NULL,
								 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								 /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								 /*sensedir*/SENSE_ANTI));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 7 */
    for (p = antidonors_minus[nmismatches]; p != NULL; p = p->rest) {
      donor = (Substring_T) p->first;
      endlength = Substring_chimera_pos(donor);
      support = querylength - endlength;

#ifdef DEBUG4H
      chroffset = Substring_chroffset(donor);
      leftbound = Substring_alignstart_trim(donor) + 1U;
#endif

      debug4h(printf("End 7: short-overlap antidonor_minus: #%d:%u, endlength %d\n",
		     Substring_chrnum(donor),leftbound-1U-chroffset,endlength));

      if (endlength <= min_localsplicing_end_matches) {
	debug4h(printf("End 7: short-overlap antidonor_minus: #%d:%u (%d mismatches) => searching right\n",
		       Substring_chrnum(donor),leftbound-1U-chroffset,Substring_nmismatches_whole(donor)));

	if ((i = Substring_splicesites_i(donor)) >= 0) {
	  origleft = Substring_genomicend(donor);
	  if ((splicesites_i =
	       Splicetrie_find_right(&nmismatches_shortend,&nmismatches_list,i,
				     origleft,/*pos5*/support,/*pos3*/querylength,
				     splicesites,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
				     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				     trieoffsets_max,triecontents_max,nsplicepartners_max,
				     splicetypes,splicestrings,query_compress_rev,
				     genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*queryptr*/queryrc,
				     nmisses_allowed,dibasep,cmetp,/*plusp*/false)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      debug4h(printf("End 7: short-overlap antidonor_minus: Successful ambiguous from antidonor #%d\n",
			     Substring_splicesites_i(donor)));
	      hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
							       donor,/*acceptor*/NULL,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,
							       /*ambi_left*/NULL,/*ambi_right*/splicesites_i,
							       /*amb_nmismatches_left*/NULL,nmismatches_list,
							       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
							       /*sensedir*/SENSE_ANTI));
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - support;
	      if ((acceptor = Substring_new_acceptor(/*splicesites_i*/bestj,/*joffset*/0,
						     querylength-Substring_chimera_pos(donor),nmismatches_shortend,/*ncolordiffs*/0,
						     /*prob*/2.0,/*left*/bestleft,query_compress_rev,genome_blocks,snp_blocks,
						     querylength,/*plusp*/false,/*sensep*/false,
						     queryuc_ptr,Substring_chrnum(donor),Substring_chroffset(donor),
						     dibasep,cmetp)) != NULL) {
		debug4h(printf("End 7: short-overlap antidonor_minus: Successful splice from antidonor #%d to antiacceptor #%d\n",
			       Substring_splicesites_i(donor),Substring_splicesites_i(acceptor)));
		hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches,nmismatches_shortend,
								 donor,acceptor,/*distance*/bestleft-origleft,
								 /*shortdistancep*/true,localsplicing_penalty,querylength,
								 /*ambi_left*/NULL,/*ambi_right*/NULL,
								 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								 /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
								 /*sensedir*/SENSE_ANTI));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

    /* End 8 */
    for (p = antiacceptors_minus[nmismatches]; p != NULL; p = p->rest) {
      acceptor = (Substring_T) p->first;
      support = Substring_chimera_pos(acceptor);
      endlength = querylength - support;

#ifdef DEBUG4H
      chroffset = Substring_chroffset(acceptor);
      rightbound = Substring_alignend_trim(acceptor);
#endif

      debug4h(printf("End 8: short-overlap antiacceptor_minus: #%d:%u, endlength %d\n",
		     Substring_chrnum(acceptor),rightbound+1-chroffset,endlength));

      if (endlength <= min_localsplicing_end_matches) {
	debug4h(printf("End 8: short-overlap antiacceptor_minus: #%d:%u (%d mismatches) => searching left\n",
		       Substring_chrnum(acceptor),rightbound+1-chroffset,Substring_nmismatches_whole(acceptor)));

	if ((i = Substring_splicesites_i(acceptor)) >= 0) {
	  origleft = Substring_genomicend(acceptor);
	  if ((splicesites_i =
	       Splicetrie_find_left(&nmismatches_shortend,&nmismatches_list,i,
				    origleft,/*pos5*/0,/*pos3*/endlength,
				    splicesites,splicefrags_ref,splicefrags_alt,nsplicepartners_skip,
				    trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				    trieoffsets_max,triecontents_max,nsplicepartners_max,
				    splicetypes,splicestrings,query_compress_rev,
				    genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*queryptr*/queryrc,
				    nmisses_allowed,dibasep,cmetp,/*plusp*/false)) != NULL) {

	    if (endlength < min_shortend || Intlist_length(splicesites_i) > 1) {
	      debug4h(printf("End 8: short-overlap antiacceptor_minus: Successful ambiguous from antiacceptor #%d\n",
			     Substring_splicesites_i(acceptor)));
	      hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
							       /*donor*/NULL,acceptor,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,
							       /*ambi_left*/splicesites_i,/*ambi_right*/NULL,
							       nmismatches_list,/*amb_nmismatches_right*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
							       /*sensedir*/SENSE_ANTI));
	    } else {
	      bestj = Intlist_head(splicesites_i);
	      bestleft = splicesites[bestj] - endlength;
	      if ((donor = Substring_new_donor(/*splicesites_i*/bestj,/*joffset*/0,
					       querylength-Substring_chimera_pos(acceptor),nmismatches_shortend,/*ncolordiffs*/0,
					       /*prob*/2.0,/*left*/bestleft,query_compress_rev,genome_blocks,snp_blocks,
					       querylength,/*plusp*/false,/*sensep*/false,
					       queryuc_ptr,Substring_chrnum(acceptor),Substring_chroffset(acceptor),
					       dibasep,cmetp)) != NULL) {
		debug4h(printf("End 8: short-overlap antiacceptor_minus: Successful splice from antiacceptor #%d to antidonor #%d\n",
			       Substring_splicesites_i(acceptor),Substring_splicesites_i(donor)));
		hits = List_push(hits,(void *) Stage3_new_splice(&(*found_score),nmismatches_shortend,nmismatches,
								 donor,acceptor,/*distance*/origleft-bestleft,
								 /*shortdistancep*/true,localsplicing_penalty,querylength,
								 /*ambi_left*/NULL,/*ambi_right*/NULL,
								 /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								 /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
								 /*sensedir*/SENSE_ANTI));
	      }
	    }
	    Intlist_free(&nmismatches_list);
	    Intlist_free(&splicesites_i);
	  }
	}
      }
    }

  }
  debug(printf("Ending find_splices_short_overlaps\n"));

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



static Floors_T
compute_floors (bool *any_omitted_p, bool *alloc_floors_p, Floors_T *floors_array,
		T this, int querylength, int query_lastpos, Indexdb_T indexdb, Indexdb_T indexdb2,
		int indexdb_size_threshold, int max_end_insertions,
		bool omit_frequent_p, bool omit_repetitive_p) {
  Floors_T floors;
  bool all_omitted_p;

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
    return (Floors_T) NULL;
  } else if (*any_omitted_p) {
    floors = Floors_new_omitted(querylength,max_end_insertions,this->omitted);
    *alloc_floors_p = true;
  } else {
#ifdef KEEP_FLOORS
    if (floors_array[querylength] == NULL) {
      floors_array[querylength] = Floors_new_standard(querylength,max_end_insertions);
    }
    floors = floors_array[querylength];
    *alloc_floors_p = false;
#else
    floors = Floors_new_standard(querylength,max_end_insertions);
    *alloc_floors_p = true;
#endif
  }

  return floors;
}


static void
complete_set_mm_indels (int *found_score, bool *segments_computed_p,
			struct Segment_T **plus_segments, int *plus_nsegments,
			struct Segment_T **minus_segments, int *minus_nsegments,
			bool *any_omitted_p, int *opt_level, int *done_level, int user_maxlevel,
			bool revise_levels_p, List_T *subs, List_T *indels, T this,
			char *queryuc_ptr, char *queryrc, int querylength, int query_lastpos,
			Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
			IIT_T chromosome_iit, UINT4 *genome_blocks, UINT4 *snp_blocks,
			Floors_T *floors_array, Genomicpos_T *splicesites, int nsplicesites,
			int subopt_levels, int indel_penalty, int max_middle_insertions, int max_middle_deletions,
			bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
			bool dibasep, bool cmetp,
			int fast_level, bool omit_frequent_p, bool omit_repetitive_p) {
  Floors_T floors;
  int indel_level;
  int max_indel_sep;
  bool alloc_floors_p = false;
  int firstbound, lastbound;
  int max_mismatches_allowed;

  debug(printf("Starting complete_set_mm_indels\n"));
  *plus_segments = NULL;
  *minus_segments = NULL;

  /* 4 and 5. Mismatches and indels via complete set.  Requires compress and
     all positions fetched.  Omits oligos and creates segments for some diagonals. */
  if (this->query_compress_fwd == NULL) {
    this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
    this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
  }

  if ((floors = compute_floors(&(*any_omitted_p),&alloc_floors_p,floors_array,
			       this,querylength,query_lastpos,indexdb,indexdb2,
			       indexdb_size_threshold,max_end_insertions,
			       omit_frequent_p,omit_repetitive_p)) == NULL) {
    debug(printf("Aborting because all oligos are omitted\n"));
    return;
  }

  if (find_12mer_bounds(&firstbound,&lastbound,this->omitted,query_lastpos) == false) {
    debug(printf("Cannot find 12_mer bounds\n"));
    allow_end_indels_p = false;
  } else {
    debug(printf("Found firstbound %d and lastbound %d\n",firstbound,lastbound));
  }

  if (*done_level < indel_penalty) {
    /* Prevents accumulation of segments for indels */
    max_indel_sep = 0;
    allow_end_indels_p = false;
  } else {
    max_indel_sep = (max_middle_insertions > max_middle_deletions) ?
      max_middle_insertions + 1 : max_middle_deletions + 1;
  }

  /* 4. Complete set mismatches */
  /* Done as a single batch */
  max_mismatches_allowed = (*done_level <= fast_level) ? -1 : *done_level;
  debug(printf("*** Stage 5.  Complete set mismatches up to %d (done_level %d, fast_level %d) ***\n",
	       max_mismatches_allowed,*done_level,fast_level));

  if (max_mismatches_allowed >= 0) {
    *plus_segments = identify_all_segments(&(*plus_nsegments),this->plus_positions,this->plus_npositions,
					   this->omitted,querylength,query_lastpos,
					   chromosome_iit,floors,splicesites,nsplicesites,/*plusp*/true);
    *minus_segments = identify_all_segments(&(*minus_nsegments),this->minus_positions,this->minus_npositions,
					    this->omitted,querylength,query_lastpos,
					    chromosome_iit,floors,splicesites,nsplicesites,/*plusp*/false);

    *subs = find_complete_mm(&(*found_score),*subs,*plus_segments,*plus_nsegments,
			     querylength,/*query*/queryuc_ptr,/*queryptr:queryuc_ptr,*/
			     /*query_compress*/this->query_compress_fwd,
			     genome_blocks,snp_blocks,chromosome_iit,
			     max_mismatches_allowed,
			     dibasep,cmetp,/*plusp*/true);

    *subs = find_complete_mm(&(*found_score),*subs,*minus_segments,*minus_nsegments,
			     querylength,/*query*/queryuc_ptr,/*queryptr:queryrc,*/
			     /*query_compress*/this->query_compress_rev,
			     genome_blocks,snp_blocks,chromosome_iit,
			     max_mismatches_allowed,
			     dibasep,cmetp,/*plusp*/false);
    *segments_computed_p = true;

    debug(printf("5> found_score = %d, opt_level %d, done_level %d\n",*found_score,*opt_level,*done_level));
    debug(printf("plus_nsegments = %d, minus_nsegments = %d\n",*plus_nsegments,*minus_nsegments));
  }

  if (revise_levels_p == true) {
    *opt_level = (*found_score < *opt_level) ? *found_score : *opt_level;
    if ((*done_level = *opt_level + subopt_levels) > user_maxlevel) {
      *done_level = user_maxlevel;
    }
  }

  if (*done_level >= indel_penalty) {
    /* 6. Indels */
    /* Need to reverse, because middle indelsplicing procedure depends on ascending diagonal order */

    if (*segments_computed_p == false) {
      *plus_segments = identify_all_segments(&(*plus_nsegments),this->plus_positions,this->plus_npositions,
					     this->omitted,querylength,query_lastpos,
					     chromosome_iit,floors,splicesites,nsplicesites,/*plusp*/true);
      *minus_segments = identify_all_segments(&(*minus_nsegments),this->minus_positions,this->minus_npositions,
					      this->omitted,querylength,query_lastpos,
					      chromosome_iit,floors,splicesites,nsplicesites,/*plusp*/false);
      *segments_computed_p = true;
    }

    /* Done iteratively */
    indel_level = indel_penalty;
    while (indel_level <= *done_level) {
      debug(printf("*** Stage 6.  Middle indels with %d-%d mismatches allowed\n",indel_level,indel_penalty));
      *indels = find_middle_indels(&(*found_score),*indels,*plus_segments,*minus_segments,
				   *plus_nsegments,*minus_nsegments,genome_blocks,snp_blocks,
				   queryuc_ptr,queryrc,floors,
				   querylength,query_lastpos,this->query_compress_fwd,this->query_compress_rev,
				   max_middle_insertions,max_middle_deletions,min_indel_end_matches,indel_penalty,
				   /*indel_mismatches_allowed*/indel_level - indel_penalty,
				   dibasep,cmetp);

      if (allow_end_indels_p == true) {
	debug(printf("*** Stage 6.  End indels with %d-%d mismatches allowed\n",indel_level,indel_penalty));
	*indels = find_end_indels(&(*found_score),*indels,*plus_segments,*minus_segments,
				  *plus_nsegments,*minus_nsegments,genome_blocks,snp_blocks,
				  queryuc_ptr,queryrc,querylength,
				  firstbound,lastbound,this->query_compress_fwd,this->query_compress_rev,
				  max_end_insertions,max_end_deletions,min_indel_end_matches,indel_penalty,
				  /*indel_mismatches_allowed*/indel_level - indel_penalty,
				  dibasep,cmetp);
      }
      if (revise_levels_p == true) {
	*opt_level = (*found_score < *opt_level) ? *found_score : *opt_level;
	if ((*done_level = *opt_level + subopt_levels) > user_maxlevel) {
	  *done_level = user_maxlevel;
	}
      }
      indel_level++;
      debug(printf("6> found_score = %d, opt_level %d, done_level %d\n",*found_score,*opt_level,*done_level));
    }
  }

  if (alloc_floors_p) {
    Floors_free(&floors);
  }

  debug(printf("Finished with complete_set_mm_indels\n"));

  return;
}



static List_T
complete_set_singlesplicing (int *found_score, Floors_T floors,
			     struct Segment_T *plus_segments, int plus_nsegments,
			     struct Segment_T *minus_segments, int minus_nsegments,
			     List_T localsplicing, T this, char *queryuc_ptr, char *queryrc, 
			     int querylength, int query_lastpos,
			     IIT_T chromosome_iit, UINT4 *genome_blocks, UINT4 *snp_blocks,
			     Genomicpos_T *splicesites, Splicetype_T *splicetypes, Genomicpos_T *splicedists,
			     int nsplicesites, Genomicpos_T shortsplicedist, bool novelsplicingp, int localsplicing_penalty,
			     int min_localsplicing_end_matches, int max_mismatches_allowed,
			     bool first_read_p, bool dibasep, bool cmetp) {

  debug(printf("Starting complete_set_singlesplicing with %d mismatches allowed\n",max_mismatches_allowed));

  if (floors == NULL) {
    return (List_T) NULL;
  }

  if (this->query_compress_fwd == NULL) {
    this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
    this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
  }

  localsplicing = find_singlesplices_plus(&(*found_score),localsplicing,plus_segments,plus_nsegments,
					  genome_blocks,snp_blocks,/*queryptr*/queryuc_ptr,
					  floors,querylength,query_lastpos,/*query_compress*/this->query_compress_fwd,
					  splicesites,splicetypes,splicedists,nsplicesites,
					  novelsplicingp,/*max_distance*/shortsplicedist,
					  /*splicing_penalty*/localsplicing_penalty,min_localsplicing_end_matches,
					  max_mismatches_allowed,first_read_p,dibasep,cmetp);

  localsplicing = find_singlesplices_minus(&(*found_score),localsplicing,minus_segments,minus_nsegments,
					   genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*queryptr*/queryrc,
					   floors,querylength,query_lastpos,/*query_compress*/this->query_compress_rev,
					   splicesites,splicetypes,splicedists,nsplicesites,
					   novelsplicingp,/*max_distance*/shortsplicedist,
					   /*splicing_penalty*/localsplicing_penalty,min_localsplicing_end_matches,
					   max_mismatches_allowed,first_read_p,dibasep,cmetp);

  debug(printf("Finished with complete_set_singlesplicing\n"));

  return localsplicing;
}


static List_T
complete_set_doublesplicing (int *found_score, Floors_T floors,
			     struct Segment_T *plus_segments, int plus_nsegments,
			     struct Segment_T *minus_segments, int minus_nsegments,
			     List_T localsplicing, T this, char *queryuc_ptr, char *queryrc, 
			     int querylength, int query_lastpos,
			     IIT_T chromosome_iit, UINT4 *genome_blocks, UINT4 *snp_blocks,
			     Genomicpos_T *splicesites, Splicetype_T *splicetypes,
			     UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int nsplicesites, int *nsplicepartners_skip,
			     unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
			     unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
			     List_T *splicestrings, Genomicpos_T shortsplicedist,
			     bool knownsplicingp, bool novelsplicingp, int localsplicing_penalty,
			     int min_localsplicing_end_matches, int min_shortend, bool find_novel_doublesplices_p,
			     int max_mismatches_allowed, bool first_read_p, bool dibasep, bool cmetp) {
  
  debug(printf("Starting complete_set_doublesplicing with %d mismatches allowed\n",max_mismatches_allowed));

  if (floors == NULL) {
    return (List_T) NULL;
  }

  if (this->query_compress_fwd == NULL) {
    this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
    this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
  }

  if (find_novel_doublesplices_p == true) {
    localsplicing = find_doublesplices(&(*found_score),localsplicing,plus_segments,plus_nsegments,
				       genome_blocks,snp_blocks,/*query*/queryuc_ptr,
				       floors,querylength,query_lastpos,/*query_compress*/this->query_compress_fwd,
				       splicesites,splicetypes,nsplicesites,/*max_distance*/shortsplicedist,
				       /*splicing_penalty*/localsplicing_penalty,min_localsplicing_end_matches,
				       max_mismatches_allowed,novelsplicingp,first_read_p,dibasep,cmetp,
				       /*plusp*/true);

    localsplicing = find_doublesplices(&(*found_score),localsplicing,minus_segments,minus_nsegments,
				       genome_blocks,snp_blocks,/*query*/queryuc_ptr,
				       floors,querylength,query_lastpos,/*query_compress*/this->query_compress_rev,
				       splicesites,splicetypes,nsplicesites,/*max_distance*/shortsplicedist,
				       /*splicing_penalty*/localsplicing_penalty,min_localsplicing_end_matches,
				       max_mismatches_allowed,novelsplicingp,first_read_p,dibasep,cmetp,
				       /*plusp*/false);
  }

  if (knownsplicingp) {
    localsplicing = find_known_doublesplices(&(*found_score),localsplicing,plus_segments,plus_nsegments,
					     genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,
					     floors,querylength,query_lastpos,/*query_compress*/this->query_compress_fwd,
					     splicesites,splicetypes,splicefrags_ref,splicefrags_alt,
					     nsplicesites,nsplicepartners_skip,
					     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
					     trieoffsets_max,triecontents_max,nsplicepartners_max,
					     splicestrings,/*max_distance*/shortsplicedist,
					     /*splicing_penalty*/localsplicing_penalty,min_localsplicing_end_matches,
					     min_shortend,max_mismatches_allowed,first_read_p,dibasep,cmetp,/*plusp*/true);

    localsplicing = find_known_doublesplices(&(*found_score),localsplicing,minus_segments,minus_nsegments,
					     genome_blocks,snp_blocks,/*query*/queryuc_ptr,/*queryptr*/queryrc,
					     floors,querylength,query_lastpos,/*query_compress*/this->query_compress_rev,
					     splicesites,splicetypes,splicefrags_ref,splicefrags_alt,
					     nsplicesites,nsplicepartners_skip,
					     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
					     trieoffsets_max,triecontents_max,nsplicepartners_max,
					     splicestrings,/*max_distance*/shortsplicedist,
					     /*splicing_penalty*/localsplicing_penalty,min_localsplicing_end_matches,
					     min_shortend,max_mismatches_allowed,first_read_p,dibasep,cmetp,/*plusp*/false);
  }


  debug(printf("Finished with complete_set_doublesplicing\n"));

  return localsplicing;
}






/* done_level should probably be renamed final_level.  opt_level
   should probably be renamed found_level or opt_level. */
static List_T
align_end (int *cutoff_level, T this, char *queryuc_ptr, int querylength, int query_lastpos,
	   Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
	   IIT_T chromosome_iit, UINT4 *genome_blocks, UINT4 *snp_blocks, Floors_T *floors_array,
	   Genomicpos_T *splicesites, Splicetype_T *splicetypes, Genomicpos_T *splicedists,
	   UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int nsplicesites, int *nsplicepartners_skip,
	   unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
	   unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
	   List_T *splicestrings, int user_maxlevel, int subopt_levels,
	   Masktype_T masktype, int terminal_penalty, int max_terminal_length,
	   int indel_penalty, int maxchimerapaths, bool knownsplicingp, bool novelsplicingp,
	   bool canonicalp, Genomicpos_T shortsplicedist,
	   int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
	   int min_distantsplicing_end_matches, double min_distantsplicing_identity, int min_shortend,
	   bool find_novel_doublesplices_p, int max_middle_insertions, int max_middle_deletions,
	   bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
	   bool allvalidp, bool dibasep, bool cmetp) {
  List_T subs = NULL, indels = NULL, singlesplicing = NULL, doublesplicing = NULL, shortendsplicing = NULL,
    distantsplicing = NULL, terminals = NULL;
  struct Segment_T *plus_segments = NULL, *minus_segments = NULL;
  int plus_nsegments, minus_nsegments;
  char queryrc[MAX_QUERYLENGTH+1];
  int found_score, done_level, opt_level, fast_level, mismatch_level, nmismatches;
  int max_splice_mismatches = -1, i;
  int nsplicepairs = 0;
  List_T *donors_plus, *antidonors_plus, *acceptors_plus, *antiacceptors_plus,
    *donors_minus, *antidonors_minus, *acceptors_minus, *antiacceptors_minus;
  bool any_omitted_p, ambiguousp, alloc_floors_p = false, floors_computed_p = false;
  Floors_T floors;
  bool segments_computed_p = false;

  make_complement_buffered(queryrc,queryuc_ptr,querylength);
  found_score = querylength;
  fast_level = (querylength + 2)/INDEX1PART - NREQUIRED_FAST;
  if (fast_level < 1 && user_maxlevel < 0) {
    fast_level = 1;		/* Do at least 1 mismatch */
  }

  if (user_maxlevel < 0) {
    *cutoff_level = fast_level;
  } else {
    *cutoff_level = user_maxlevel;
  }

  if (user_maxlevel < 0) {
    user_maxlevel = fast_level;
  }

  if (dibasep) {
    opt_level = querylength;	/* Allow extra because color errors may exceed nt errors */
  } else {
    opt_level = user_maxlevel;
  }
  done_level = opt_level /* + subopt_levels.  -- Initially the same */;
  debug(printf("0> opt_level %d, done_level %d\n",opt_level,done_level));

  /* 1. Exact.  Requires compress if cmet or genomealt.  Creates and uses spanning set. */
  mismatch_level = 0;
  if (allvalidp == false) {
    debug(printf("Not all oligos are valid, so cannot perform spanning set\n"));
    fast_level = -1;
  } else {
    debug(printf("fast_level = %d\n",fast_level));
    debug(printf("*** Stage 1.  Exact ***\n"));
    if (cmetp || snp_blocks) {
      this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
      this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
    }
    subs = find_spanning_exact_matches(&found_score,this,queryuc_ptr,queryrc,
				       querylength,query_lastpos,indexdb,indexdb2,
				       this->query_compress_fwd,this->query_compress_rev,
				       genome_blocks,snp_blocks,chromosome_iit,
				       dibasep,cmetp);
    opt_level = (found_score < opt_level) ? found_score : opt_level;
    if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
      done_level = user_maxlevel;
    }
    mismatch_level = 1;
    debug(printf("1> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
  }

  if (dibasep == false) {
    /* 2. One mismatch.  Requires spanning set and compress. */
    if (allvalidp && querylength >= 22 && done_level >= 1) {
      debug(printf("*** Stage 2.  One miss ***\n"));
      if (this->query_compress_fwd == NULL) {
	this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
	this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
      }
      subs = find_spanning_onemiss_matches(&found_score,subs,this,queryuc_ptr,queryrc,
					   querylength,this->query_compress_fwd,this->query_compress_rev,
					   genome_blocks,snp_blocks,chromosome_iit,
					   dibasep,cmetp);
      opt_level = (found_score < opt_level) ? found_score : opt_level;
      if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	done_level = user_maxlevel;
      }
      mismatch_level = 2;
      debug(printf("2> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
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
					       this->query_compress_fwd,this->query_compress_rev,
					       genome_blocks,snp_blocks,chromosome_iit,
					       /*nmisses_allowed*/mismatch_level,dibasep,cmetp);
	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	  done_level = user_maxlevel;
	}
	mismatch_level++;
	debug(printf("3> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
      }
    }
  }

  /* 4, 5.  Complete set mismatches and indels, omitting frequent oligos */
  if (dibasep || done_level > fast_level || done_level >= indel_penalty) {
    if (masktype == MASK_NONE) {
      debug(printf("*** Stage 4,5.  Complete mm/indels with no masking with done_level %d ***\n",done_level));
      complete_set_mm_indels(&found_score,&segments_computed_p,&plus_segments,&plus_nsegments,&minus_segments,&minus_nsegments,
			     &any_omitted_p,&opt_level,&done_level,user_maxlevel,/*revise_levels_p*/true,
			     &subs,&indels,this,queryuc_ptr,queryrc,
			     querylength,query_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,floors_array,
			     splicesites,nsplicesites,subopt_levels,
			     indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     dibasep,cmetp,fast_level,
			     /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
    } else {
      debug(printf("*** Stage 4,5.  Complete mm/indels masking frequent oligos with done_level %d ***\n",done_level));
      complete_set_mm_indels(&found_score,&segments_computed_p,&plus_segments,&plus_nsegments,&minus_segments,&minus_nsegments,
			     &any_omitted_p,&opt_level,&done_level,user_maxlevel,/*revise_levels_p*/true,
			     &subs,&indels,this,queryuc_ptr,queryrc,
			     querylength,query_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,floors_array,
			     splicesites,nsplicesites,subopt_levels,
			     indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     dibasep,cmetp,fast_level,/*omit_frequent_p*/true,
			     /*omit_repetitive_p*/(masktype == MASK_REPETITIVE || masktype == MASK_GREEDY_REPETITIVE) ? true : false
			     );
      if ((masktype == MASK_GREEDY_FREQUENT || masktype == MASK_GREEDY_REPETITIVE) && subs == NULL && indels == NULL && any_omitted_p == true) {
	FREE(minus_segments);
	FREE(plus_segments);

	debug(printf("*** Stage 4,5.  Complete mm/indels with no masking with done_level %d ***\n",done_level));
	complete_set_mm_indels(&found_score,&segments_computed_p,&plus_segments,&plus_nsegments,&minus_segments,&minus_nsegments,
			       &any_omitted_p,&opt_level,&done_level,user_maxlevel,/*revise_levels_p*/true,
			       &subs,&indels,this,queryuc_ptr,queryrc,
			       querylength,query_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			       chromosome_iit,genome_blocks,snp_blocks,floors_array,
			       splicesites,nsplicesites,subopt_levels,
			       indel_penalty,max_middle_insertions,max_middle_deletions,
			       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			       dibasep,cmetp,fast_level,
			       /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
      }
    }
  }

  /* 6, 7, 8, 9.  Splicing.  Requires compress and all positions fetched */
  if (knownsplicingp || novelsplicingp) {
    /* 6.  Single splicing */
    debug(printf("Deciding whether to do singlesplicing: done_level %d >=? localsplicing_penalty %d\n",
		 done_level,localsplicing_penalty));
    if (done_level >= localsplicing_penalty) {
      debug(printf("*** Stage 6.  Single splicing masking frequent oligos with done_level %d ***\n",done_level));
      /* Always mask frequent oligos for splicing, which must be transcriptional */
      floors = compute_floors(&any_omitted_p,&alloc_floors_p,floors_array,this,querylength,query_lastpos,
			      indexdb,indexdb2,indexdb_size_threshold,max_end_insertions,
			      /*omit_frequent_p*/true,/*omit_repetitive_p*/true);
      floors_computed_p = true;

      if (segments_computed_p == false) {
	plus_segments = identify_all_segments(&plus_nsegments,this->plus_positions,this->plus_npositions,
					      this->omitted,querylength,query_lastpos,
					      chromosome_iit,floors,splicesites,nsplicesites,/*plusp*/true);
	minus_segments = identify_all_segments(&minus_nsegments,this->minus_positions,this->minus_npositions,
					       this->omitted,querylength,query_lastpos,
					       chromosome_iit,floors,splicesites,nsplicesites,/*plusp*/false);
	segments_computed_p = true;
      }

      singlesplicing = complete_set_singlesplicing(&found_score,floors,
						   plus_segments,plus_nsegments,
						   minus_segments,minus_nsegments,
						   singlesplicing,this,queryuc_ptr,queryrc,
						   querylength,query_lastpos,chromosome_iit,
						   genome_blocks,snp_blocks,
						   splicesites,splicetypes,splicedists,nsplicesites,
						   shortsplicedist,novelsplicingp,
						   localsplicing_penalty,min_localsplicing_end_matches,
						   /*max_mismatches_allowed*/done_level - localsplicing_penalty,
						   /*first_read_p*/NOT_APPLICABLE,dibasep,cmetp);

#if 0
      /* Mark ambiguous splices only for single-end reads */
      singlesplicing = Stage3_mark_ambiguous_splices(&ambiguousp,singlesplicing);
#endif
      singlesplicing = Stage3_optimal_score(singlesplicing,/*cutoff_level*/opt_level,subopt_levels);

      if (singlesplicing) {
	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	  done_level = user_maxlevel;
	}
      }
    }

    /* 7.  Double splicing */
    debug(printf("Deciding whether to do doublesplicing: done_level %d >=? localsplicing_penalty %d\n",
		 done_level,localsplicing_penalty));
    if (done_level >= localsplicing_penalty) {
      debug(printf("*** Stage 7.  Double splicing masking frequent oligos with done_level %d ***\n",done_level));
      doublesplicing = complete_set_doublesplicing(&found_score,floors,plus_segments,plus_nsegments,
						   minus_segments,minus_nsegments,
						   doublesplicing,this,queryuc_ptr,queryrc,
						   querylength,query_lastpos,chromosome_iit,
						   genome_blocks,snp_blocks,
						   splicesites,splicetypes,splicefrags_ref,splicefrags_alt,
						   nsplicesites,nsplicepartners_skip,
						   trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
						   trieoffsets_max,triecontents_max,nsplicepartners_max,
						   splicestrings,shortsplicedist,knownsplicingp,novelsplicingp,
						   localsplicing_penalty,min_localsplicing_end_matches,min_shortend,
						   find_novel_doublesplices_p,/*max_mismatches_allowed*/done_level - localsplicing_penalty,
						   /*first_read_p*/NOT_APPLICABLE,dibasep,cmetp);
      
#if 0
      /* Mark ambiguous splices only for single-end reads */
      doublesplicing = Stage3_mark_ambiguous_splices(&ambiguousp,doublesplicing);
#endif
      doublesplicing = Stage3_optimal_score(doublesplicing,/*cutoff_level*/opt_level,subopt_levels);

      if (doublesplicing) {
	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	  done_level = user_maxlevel;
	}
      }
    }

    /* 8, 9.  Shortends and distant splicing */
    if (knownsplicingp == true && localsplicing_penalty < distantsplicing_penalty) {
      max_splice_mismatches = done_level - localsplicing_penalty;
      debug(printf("max_splice_mismatches %d = done_level %d - localsplicing_penalty %d\n",
		   max_splice_mismatches,done_level,localsplicing_penalty));
    } else {
      max_splice_mismatches = done_level - distantsplicing_penalty;
      debug(printf("max_splice_mismatches %d = done_level %d - distantsplicing_penalty %d\n",
		   max_splice_mismatches,done_level,distantsplicing_penalty));
    }

    if (max_splice_mismatches >= 0) {
      donors_plus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      antidonors_plus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      acceptors_plus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      antiacceptors_plus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      donors_minus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      antidonors_minus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      acceptors_minus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));
      antiacceptors_minus = (List_T *) CALLOC(max_splice_mismatches+1,sizeof(List_T));

      debug(printf("Starting find_spliceends (plus)\n"));
      find_spliceends(&donors_plus,&antidonors_plus,&acceptors_plus,&antiacceptors_plus,
		      plus_segments,plus_nsegments,genome_blocks,snp_blocks,
		      /*query*/queryuc_ptr,/*queryptr*/queryuc_ptr,
		      floors,querylength,query_lastpos,
		      /*query_compress*/this->query_compress_fwd,
		      splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
		      /*max_mismatches_allowed*/max_splice_mismatches,
		      dibasep,cmetp,/*plusp*/true);
      debug(printf("Finished find_spliceends (plus)\n"));

      debug(printf("Starting find_spliceends (minus)\n"));
      find_spliceends(&antidonors_minus,&donors_minus,&antiacceptors_minus,&acceptors_minus,
		      minus_segments,minus_nsegments,genome_blocks,snp_blocks,
		      /*query*/queryuc_ptr,/*queryptr*/queryrc,
		      floors,querylength,query_lastpos,
		      /*query_compress*/this->query_compress_rev,
		      splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
		      /*max_mismatches_allowed*/max_splice_mismatches,
		      dibasep,cmetp,/*plusp*/false);
      debug(printf("Finished find_spliceends (minus)\n"));

      if (knownsplicingp == true) {
	/* 8.  Find short-overlap splicing using known splice sites */
	shortendsplicing = find_splicepairs_short_overlaps(&found_score,shortendsplicing,
							   donors_plus,antidonors_plus,acceptors_plus,antiacceptors_plus,
							   donors_minus,antidonors_minus,acceptors_minus,antiacceptors_minus,
							   this->query_compress_fwd,this->query_compress_rev,
							   genome_blocks,snp_blocks,queryuc_ptr,queryrc,
							   min_shortend,
							   splicesites,splicetypes,splicefrags_ref,splicefrags_alt,
							   nsplicesites,nsplicepartners_skip,
							   trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
							   trieoffsets_max,triecontents_max,nsplicepartners_max,
							   splicestrings,shortsplicedist,localsplicing_penalty,
							   min_localsplicing_end_matches,
							   /*max_mismatches_allowed*/max_splice_mismatches,querylength,
							   /*first_read_p*/NOT_APPLICABLE,dibasep,cmetp);
	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	  done_level = user_maxlevel;
	}
	debug(printf("8> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
      }


      /* 9.  Find distant splicing iteratively using both known and novel splice sites */
      nmismatches = 0;
      ambiguousp = false;
      while (nmismatches <= done_level - distantsplicing_penalty &&
	     nsplicepairs < maxchimerapaths && ambiguousp == false) {
	debug(printf("*** Stage 8.  Distant splicing, allowing %d mismatches ***\n",nmismatches));

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
						   querylength,nmismatches,maxchimerapaths,/*first_read_p*/NOT_APPLICABLE);

#if 0
	/* Mark ambiguous splices only for single-end reads */
	distantsplicing = Stage3_mark_ambiguous_splices(&ambiguousp,distantsplicing);
#endif
	distantsplicing = Stage3_optimal_score(distantsplicing,opt_level,subopt_levels);
	if (distantsplicing) {
	  opt_level = (found_score < opt_level) ? found_score : opt_level;
	  if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
	    done_level = user_maxlevel;
	  }
	  debug(printf("8> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
	}
	nmismatches++;
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
    debug(printf("%d single splices, %d double splices, %d short-end splices, %d distant splices\n",
		 List_length(singlesplicing),List_length(doublesplicing),
		 List_length(shortendsplicing),List_length(distantsplicing)));
  }

  if (done_level >= terminal_penalty) {
    /* 10.  Terminals */
    debug(printf("*** Stage 10.  Terminals up to %d mismatches ***\n",10));
    if (floors_computed_p == false) {
      floors = compute_floors(&any_omitted_p,&alloc_floors_p,floors_array,this,querylength,query_lastpos,
			      indexdb,indexdb2,indexdb_size_threshold,max_end_insertions,
			      /*omit_frequent_p*/true,/*omit_repetitive_p*/true);
    }

    if (segments_computed_p == false) {
      plus_segments = identify_all_segments_for_terminals(&plus_nsegments,this->plus_positions,this->plus_npositions,
							  this->omitted,querylength,query_lastpos,
							  chromosome_iit,floors,
							  /*max_mismatches_allowed*/done_level - terminal_penalty,
							  /*plusp*/true);
      minus_segments = identify_all_segments_for_terminals(&minus_nsegments,this->minus_positions,this->minus_npositions,
							   this->omitted,querylength,query_lastpos,
							   chromosome_iit,floors,
							  /*max_mismatches_allowed*/done_level - terminal_penalty,
							   /*plusp*/false);
    }

    if (this->query_compress_fwd == NULL) {
      this->query_compress_fwd = Compress_new(queryuc_ptr,querylength,/*plusp*/true);
      this->query_compress_rev = Compress_new(queryuc_ptr,querylength,/*plusp*/false);
    }
    terminals = find_terminals(plus_segments,plus_nsegments,minus_segments,minus_nsegments,
			       genome_blocks,snp_blocks,
			       queryuc_ptr,queryrc,floors,querylength,query_lastpos,
			       this->query_compress_fwd,this->query_compress_rev,
			       /*max_mismatches_allowed*/done_level - terminal_penalty,
			       /*terminal_penalty*/terminal_penalty,max_terminal_length,
			       dibasep,cmetp);
#if 0
    opt_level = (found_score < opt_level) ? found_score : opt_level;
    if ((done_level = opt_level + subopt_levels) > user_maxlevel) {
      done_level = user_maxlevel;
    }
    debug(printf("10> found_score = %d, opt_level %d, done_level %d\n",found_score,opt_level,done_level));
#endif
  }

  if (alloc_floors_p == true) {
    Floors_free(&floors);
  }

  FREE(minus_segments);
  FREE(plus_segments);

  return List_append(subs,
		     List_append(indels,
				 List_append(singlesplicing,
					     List_append(doublesplicing,
							 List_append(shortendsplicing,
								     List_append(distantsplicing,terminals))))));
  }


Stage3_T *
Stage1_single_read (int *npaths, Shortread_T queryseq, Indexdb_T indexdb, Indexdb_T indexdb2,
		    int indexdb_size_threshold, IIT_T chromosome_iit, Genome_T genome,
		    Genome_T genomealt, Floors_T *floors_array,
		    bool knownsplicingp, bool novelsplicingp, bool canonicalp,
		    int maxpaths, int maxchimerapaths, double user_maxlevel_float, int subopt_levels,
		    Masktype_T masktype, int terminal_penalty, int max_terminal_length,
		    int indel_penalty, int max_middle_insertions, int max_middle_deletions,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    Genomicpos_T shortsplicedist,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
		    int min_distantsplicing_end_matches, double min_distantsplicing_identity, int min_shortend,
		    bool find_novel_doublesplices_p, Genomicpos_T *splicesites, Splicetype_T *splicetypes,
		    Genomicpos_T *splicedists, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int nsplicesites,
		    int *nsplicepartners_skip,
		    unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
		    unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
		    List_T *splicestrings, bool dibasep, bool cmetp) {
  Stage3_T *stage3array;
  List_T hits = NULL;
  T this = NULL;
  int user_maxlevel;
  int querylength, query_lastpos, cutoff_level;
  char *queryuc_ptr, *quality_string;
  bool allvalidp;
  UINT4 *genome_blocks, *snp_blocks;

  if ((querylength = Shortread_fulllength(queryseq)) <= INDEX1PART+2) {
    fprintf(stderr,"GSNAP cannot handle reads shorter than %d bp\n",INDEX1PART+2);
    *npaths = 0;
    return (Stage3_T *) NULL;

  } else if (querylength > MAX_QUERYLENGTH) {
    fprintf(stderr,"GSNAP cannot handle reads longer than %d bp.  Either recompile with a higher value of MAX_QUERYLENGTH, or consider using GMAP instead.\n",
	    MAX_QUERYLENGTH);
    *npaths = 0;
    return (Stage3_T *) NULL;

  } else {
    if (user_maxlevel_float < 0.0) {
      user_maxlevel = -1;
    } else if (user_maxlevel_float > 0.0 && user_maxlevel_float < 1.0) {
      user_maxlevel = (int) rint(user_maxlevel_float * (double) querylength);
    } else {
      user_maxlevel = (int) user_maxlevel_float;
    }

    this = Stage1_new(querylength);
    queryuc_ptr = Shortread_fullpointer_uc(queryseq);
    quality_string = Shortread_quality_string(queryseq);
    query_lastpos = querylength - INDEX1PART;
    if (read_oligos(&allvalidp,this,queryuc_ptr,querylength,query_lastpos,dibasep,cmetp) == 0) {
      debug(printf("Aborting because no hits found anywhere\n"));
      *npaths = 0;
      Stage1_free(&this,querylength);
      return (Stage3_T *) NULL;

    } else {
      genome_blocks = Genome_blocks(genome);
      snp_blocks = genomealt ? Genome_blocks(genomealt) : NULL;
      hits = align_end(&cutoff_level,this,queryuc_ptr,querylength,query_lastpos,
		       indexdb,indexdb2,indexdb_size_threshold,
		       chromosome_iit,genome_blocks,snp_blocks,floors_array,
		       splicesites,splicetypes,splicedists,splicefrags_ref,splicefrags_alt,
		       nsplicesites,nsplicepartners_skip,
		       trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
		       trieoffsets_max,triecontents_max,nsplicepartners_max,
		       splicestrings,user_maxlevel,subopt_levels,masktype,
		       terminal_penalty,max_terminal_length,indel_penalty,
		       maxchimerapaths,knownsplicingp,novelsplicingp,
		       canonicalp,shortsplicedist,
		       localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
		       min_distantsplicing_end_matches,min_distantsplicing_identity,min_shortend,
		       find_novel_doublesplices_p,max_middle_insertions,max_middle_deletions,
		       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
		       allvalidp,dibasep,cmetp);
      debug(printf("Before remove_duplicates at cutoff level %d: %d\n",cutoff_level,List_length(hits)));
      hits = Stage3_remove_duplicates(Stage3_optimal_score(hits,cutoff_level,subopt_levels));
      debug(printf("After remove_duplicates: %d\n",List_length(hits)));

      if ((*npaths = List_length(hits)) == 0) {
	stage3array = (Stage3_T *) NULL;
      } else {
	stage3array = (Stage3_T *) List_to_array(hits,NULL);
	stage3array = Stage3_eval_and_sort(stage3array,*npaths,maxpaths,queryseq,
					   this->query_compress_fwd,this->query_compress_rev,
					   genome_blocks,snp_blocks,genome,quality_string,
					   dibasep,cmetp);
	List_free(&hits);
      }

      Stage1_free(&this,querylength); 
      return stage3array;
    }
  }
}


#define HITARRAY_SUBS 0
#define HITARRAY_INDELS 1
#define HITARRAY_SINGLESPLICING 2
#define HITARRAY_DOUBLESPLICING 3
#define HITARRAY_SHORTENDSPLICING 4
#define HITARRAY_DISTANTSPLICING 5
#define HITARRAY_TERMINALS 6
#define HITARRAY_N 7

static List_T
align_pair (bool *abort_pairing_p, int *nconcordant, int *cutoff_level_5, int *cutoff_level_3,
	    List_T *samechr, List_T *hits5, List_T *hits3, T this5, T this3,
	    char *queryuc_ptr_5, char *queryuc_ptr_3,
	    int querylength5, int querylength3, int query5_lastpos, int query3_lastpos,
	    Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
	    IIT_T chromosome_iit, UINT4 *genome_blocks, UINT4 *snp_blocks, Floors_T *floors_array,
	    Genomicpos_T *splicesites, Splicetype_T *splicetypes, Genomicpos_T *splicedists,
	    UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int nsplicesites, int *nsplicepartners_skip,
	    unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
	    unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
	    List_T *splicestrings, int user_maxlevel_5, int user_maxlevel_3, int subopt_levels,
	    Masktype_T masktype, int terminal_penalty, int max_terminal_length, int indel_penalty,
	    int maxchimerapaths, bool knownsplicingp, bool novelsplicingp,
	    bool canonicalp, Genomicpos_T shortsplicedist,
	    int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
	    int min_distantsplicing_end_matches, double min_distantsplicing_identity, int min_shortend,
	    bool find_novel_doublesplices_p, int max_middle_insertions, int max_middle_deletions,
	    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
	    bool allvalidp5, bool allvalidp3, bool dibasep, bool cmetp,
	    Genomicpos_T pairmax, Genomicpos_T expected_pairlength, Genomicpos_T pairlength_deviation,
	    int maxpairedpaths) {

  List_T hitpairs = NULL;
  List_T hitarray5[HITARRAY_N], hitarray3[HITARRAY_N];
  List_T subs5 = NULL, indels5 = NULL, singlesplicing5 = NULL, doublesplicing5 = NULL,
    shortends5 = NULL, distantsplicing5 = NULL, terminals5 = NULL;
  List_T subs3 = NULL, indels3 = NULL, singlesplicing3 = NULL, doublesplicing3 = NULL,
    shortends3 = NULL, distantsplicing3 = NULL, terminals3 = NULL;
  struct Segment_T *plus_segments_5 = NULL, *minus_segments_5 = NULL,
    *plus_segments_3 = NULL, *minus_segments_3 = NULL;
  int plus_nsegments_5, minus_nsegments_5, plus_nsegments_3, minus_nsegments_3;
  char queryrc5[MAX_QUERYLENGTH+1], queryrc3[MAX_QUERYLENGTH+1];
  int found_score, ignore_found_score, done_level_5, done_level_3, opt_level, fast_level_5, fast_level_3,
    mismatch_level_5, mismatch_level_3, nmismatches;
  int max_splice_mismatches_5 = -1, max_splice_mismatches_3 = -1, i;
  int nsplicepairs5 = 0, nsplicepairs3 = 0;
  List_T *donors_plus_5, *antidonors_plus_5, *acceptors_plus_5, *antiacceptors_plus_5,
    *donors_minus_5, *antidonors_minus_5, *acceptors_minus_5, *antiacceptors_minus_5;
  List_T *donors_plus_3, *antidonors_plus_3, *acceptors_plus_3, *antiacceptors_plus_3,
    *donors_minus_3, *antidonors_minus_3, *acceptors_minus_3, *antiacceptors_minus_3;

  bool any_omitted_p_5, any_omitted_p_3;
  Floors_T floors5, floors3;
  bool alloc_floors_p_5 = false, alloc_floors_p_3 = false, floors5_computed_p = false, floors3_computed_p = false,
    segments5_computed_p = false, segments3_computed_p = false, alloc5p = false, alloc3p = false;


  make_complement_buffered(queryrc5,queryuc_ptr_5,querylength5);
  make_complement_buffered(queryrc3,queryuc_ptr_3,querylength3);

  *nconcordant = 0;
  *samechr = (List_T) NULL;

  /* For paired-end alignment, ignore found_scores from single-end
     alignments.  Use only the found_score from
     Stage3_pair_up_concordant. */
  found_score = querylength5 + querylength3;
  ignore_found_score = querylength5 + querylength3;

  fast_level_5 = (querylength5 + 2)/INDEX1PART - NREQUIRED_FAST;
  fast_level_3 = (querylength3 + 2)/INDEX1PART - NREQUIRED_FAST;
  if (fast_level_5 < 1 && user_maxlevel_5 < 0) {
    fast_level_5 = 1;		/* Do at least 1 mismatch */
  }
  if (fast_level_3 < 1 && user_maxlevel_3 < 0) {
    fast_level_3 = 1;		/* Do at least 1 mismatch */
  }

  if (user_maxlevel_5 < 0) {
    user_maxlevel_5 = fast_level_5;
  }
  if (user_maxlevel_3 < 0) {
    user_maxlevel_3 = fast_level_3;
  }

  *cutoff_level_5 = user_maxlevel_5;
  *cutoff_level_3 = user_maxlevel_3;

  if (dibasep) {
    opt_level = querylength5 + querylength3;
    done_level_5 = querylength5;
    done_level_3 = querylength3;
  } else {
    opt_level = user_maxlevel_5 + user_maxlevel_3;
    done_level_5 = user_maxlevel_5 /* + subopt_levels */;
    done_level_3 = user_maxlevel_3 /* + subopt_levels */;
  }
  debug(printf("0> opt_level %d, done_level %d,%d\n",opt_level,done_level_5,done_level_3));

  /* 1A. Exact.  Requires compress if cmet or genomealt.  Creates and uses spanning set. */
  mismatch_level_5 = 0;
  if (allvalidp5 == false) {
    debug(printf("Not all oligos in 5' end are valid, so cannot perform spanning set\n"));
    fast_level_5 = -1;
  } else {
    debug(printf("fast_level_5 = %d\n",fast_level_5));
    debug(printf("*** Stage 1.  Exact ***\n"));
    if (cmetp || snp_blocks) {
      this5->query_compress_fwd = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/true);
      this5->query_compress_rev = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/false);
    }
    subs5 = find_spanning_exact_matches(&ignore_found_score,this5,queryuc_ptr_5,queryrc5,
					querylength5,query5_lastpos,indexdb,indexdb2,
					this5->query_compress_fwd,this5->query_compress_rev,
					genome_blocks,snp_blocks,chromosome_iit,
					dibasep,cmetp);
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
    if (cmetp || snp_blocks) {
      this3->query_compress_fwd = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/true);
      this3->query_compress_rev = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/false);
    }
    subs3 = find_spanning_exact_matches(&ignore_found_score,this3,queryuc_ptr_3,queryrc3,
					querylength3,query3_lastpos,indexdb,indexdb2,
					this3->query_compress_fwd,this3->query_compress_rev,
					genome_blocks,snp_blocks,chromosome_iit,
					dibasep,cmetp);
    mismatch_level_3 = 1;
  }

  for (i = 0; i < HITARRAY_N; i++) {
    hitarray5[i] = hitarray3[i] = (List_T) NULL;
  }

  /* 1. Pairing after exact */
  hitarray5[HITARRAY_SUBS] = subs5 = Stage3_remove_duplicates(subs5);
  hitarray3[HITARRAY_SUBS] = subs3 = Stage3_remove_duplicates(subs3);
  hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*samechr),hitpairs,
				       hitarray5,/*narray5*/HITARRAY_SUBS+1,
				       hitarray3,/*narray3*/HITARRAY_SUBS+1,
				       *cutoff_level_5,*cutoff_level_3,subopt_levels,
				       splicesites,queryuc_ptr_5,queryuc_ptr_3,this5->query_compress_fwd,this5->query_compress_rev,
				       this3->query_compress_fwd,this3->query_compress_rev,genome_blocks,snp_blocks,
				       pairmax,expected_pairlength,pairlength_deviation,
				       querylength5,querylength3,maxpairedpaths,
				       /*allow_concordant_translocations_p*/true,localsplicing_penalty,dibasep,cmetp);
  debug(printf("After pairing exact, found %d concordant, found_score %d\n",*nconcordant,found_score));
  if (*abort_pairing_p == true) {
    *hits5 = subs5;
    *hits3 = subs3;
    return hitpairs;
  } else {
    opt_level = (found_score < opt_level) ? found_score : opt_level;
    if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
      done_level_5 = user_maxlevel_5;
    }
    if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
      done_level_3 = user_maxlevel_3;
    }
    debug(printf("1> found_score = %d, opt_level %d, done_level %d,%d\n",found_score,opt_level,done_level_5,done_level_3));
  }


  if (dibasep == false) {
    /* 2A. One mismatch.  Requires spanning set and compress. */
    if (allvalidp5 && querylength5 >= 22 && done_level_5 >= 1) {
      debug(printf("*** Stage 2A.  One miss ***\n"));
      if (this5->query_compress_fwd == NULL) {
	this5->query_compress_fwd = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/true);
	this5->query_compress_rev = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/false);
      }
      subs5 = find_spanning_onemiss_matches(&ignore_found_score,subs5,this5,
					    queryuc_ptr_5,queryrc5,querylength5,
					    this5->query_compress_fwd,this5->query_compress_rev,
					    genome_blocks,snp_blocks,chromosome_iit,
					    dibasep,cmetp);
      mismatch_level_5 = 2;
    }

    /* 2B. One mismatch.  Requires spanning set and compress. */
    if (allvalidp3 && querylength3 >= 22 && done_level_3 >= 1) {
      debug(printf("*** Stage 2B.  One miss ***\n"));
      if (this3->query_compress_fwd == NULL) {
	this3->query_compress_fwd = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/true);
	this3->query_compress_rev = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/false);
      }
      subs3 = find_spanning_onemiss_matches(&ignore_found_score,subs3,this3,
					    queryuc_ptr_3,queryrc3,querylength3,
					    this3->query_compress_fwd,this3->query_compress_rev,
					    genome_blocks,snp_blocks,chromosome_iit,
					    dibasep,cmetp);
      mismatch_level_3 = 2;
    }

    /* 2. Pairing after one mismatch */
    hitarray5[HITARRAY_SUBS] = subs5 = Stage3_remove_duplicates(subs5);
    hitarray3[HITARRAY_SUBS] = subs3 = Stage3_remove_duplicates(subs3);
    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*samechr),hitpairs,
					 hitarray5,/*narray5*/HITARRAY_SUBS+1,
					 hitarray3,/*narray3*/HITARRAY_SUBS+1,
					 *cutoff_level_5,*cutoff_level_3,subopt_levels,
					 splicesites,queryuc_ptr_5,queryuc_ptr_3,this5->query_compress_fwd,this5->query_compress_rev,
					 this3->query_compress_fwd,this3->query_compress_rev,genome_blocks,snp_blocks,
					 pairmax,expected_pairlength,pairlength_deviation,
					 querylength5,querylength3,maxpairedpaths,
					 /*allow_concordant_translocations_p*/true,localsplicing_penalty,dibasep,cmetp);
    debug(printf("After pairing one mismatch, found %d concordant, found_score %d\n",*nconcordant,found_score));
    if (*abort_pairing_p == true) {
      *hits5 = subs5;
      *hits3 = subs3;
      return hitpairs;
    } else {
      opt_level = (found_score < opt_level) ? found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("2> found_score = %d, opt_level %d, done_level %d,%d\n",found_score,opt_level,done_level_5,done_level_3));
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
	subs5 = find_spanning_multimiss_matches(&ignore_found_score,subs5,this5,NREQUIRED_FAST,
						queryuc_ptr_5,queryrc5,querylength5,
						this5->query_compress_fwd,this5->query_compress_rev,
						genome_blocks,snp_blocks,chromosome_iit,
						/*nmisses_allowed*/mismatch_level_5,dibasep,cmetp);
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
	subs3 = find_spanning_multimiss_matches(&ignore_found_score,subs3,this3,NREQUIRED_FAST,
						queryuc_ptr_3,queryrc3,querylength3,
						this3->query_compress_fwd,this3->query_compress_rev,
						genome_blocks,snp_blocks,chromosome_iit,
						/*nmisses_allowed*/mismatch_level_3,dibasep,cmetp);
	mismatch_level_3++;
      }
    }

    /* 3. Pairing after spanning set subs */
    hitarray5[HITARRAY_SUBS] = subs5 = Stage3_remove_duplicates(subs5);
    hitarray3[HITARRAY_SUBS] = subs3 = Stage3_remove_duplicates(subs3);
    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*samechr),hitpairs,
					 hitarray5,/*narray5*/HITARRAY_SUBS+1,
					 hitarray3,/*narray3*/HITARRAY_SUBS+1,
					 *cutoff_level_5,*cutoff_level_3,subopt_levels,
					 splicesites,queryuc_ptr_5,queryuc_ptr_3,this5->query_compress_fwd,this5->query_compress_rev,
					 this3->query_compress_fwd,this3->query_compress_rev,genome_blocks,snp_blocks,
					 pairmax,expected_pairlength,pairlength_deviation,
					 querylength5,querylength3,maxpairedpaths,
					 /*allow_concordant_translocations_p*/true,localsplicing_penalty,dibasep,cmetp);
    debug(printf("After pairing spanning set, found %d concordant, found_score %d\n",*nconcordant,found_score));
    if (*abort_pairing_p == true) {
      *hits5 = subs5;
      *hits3 = subs3;
      return hitpairs;
    } else {
      opt_level = (found_score < opt_level) ? found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("3> found_score = %d, opt_level %d, done_level %d,%d\n",found_score,opt_level,done_level_5,done_level_3));
    }
  }


  /* 4/5A.  Complete set mismatches and indels, omitting frequent oligos */
  if (dibasep || done_level_5 > fast_level_5 || done_level_5 >= indel_penalty) {
    if (masktype == MASK_NONE) {
      debug(printf("*** Stage 4A,5A.  Complete mm/indels with no masking with done_level %d ***\n",done_level_5));
      complete_set_mm_indels(&ignore_found_score,&segments5_computed_p,&plus_segments_5,&plus_nsegments_5,&minus_segments_5,&minus_nsegments_5,
			     &any_omitted_p_5,&opt_level,&done_level_5,user_maxlevel_5,/*revise_levels_p*/false,
			     &subs5,&indels5,this5,queryuc_ptr_5,queryrc5,
			     querylength5,query5_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,floors_array,
			     splicesites,nsplicesites,subopt_levels,
			     indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     dibasep,cmetp,fast_level_5,
			     /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
    } else {
      debug(printf("*** Stage 4A,5A.  Complete mm/indels masking frequent oligos with done_level %d ***\n",done_level_5));
      complete_set_mm_indels(&ignore_found_score,&segments5_computed_p,&plus_segments_5,&plus_nsegments_5,&minus_segments_5,&minus_nsegments_5,
			     &any_omitted_p_5,&opt_level,&done_level_5,user_maxlevel_5,/*revise_levels_p*/false,
			     &subs5,&indels5,this5,queryuc_ptr_5,queryrc5,
			     querylength5,query5_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,floors_array,
			     splicesites,nsplicesites,subopt_levels,
			     indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     dibasep,cmetp,fast_level_5,/*omit_frequent_p*/true,
			     /*omit_repetitive_p*/(masktype == MASK_REPETITIVE || masktype == MASK_GREEDY_REPETITIVE) ? true : false
			     );
      if ((masktype == MASK_GREEDY_FREQUENT || masktype == MASK_GREEDY_REPETITIVE) && subs5 == NULL && indels5 == NULL && any_omitted_p_5 == true) {
	FREE(minus_segments_5);
	FREE(plus_segments_5);

	/* 4/5A.  Complete set mismatches and indels, with all oligos */
	debug(printf("*** Stage 4A,5A.  Complete mm/indels with no masking with done_level %d ***\n",done_level_5));
	complete_set_mm_indels(&ignore_found_score,&segments5_computed_p,&plus_segments_5,&plus_nsegments_5,&minus_segments_5,&minus_nsegments_5,
			       &any_omitted_p_5,&opt_level,&done_level_5,user_maxlevel_5,/*revise_levels_p*/false,
			       &subs5,&indels5,this5,queryuc_ptr_5,queryrc5,
			       querylength5,query5_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			       chromosome_iit,genome_blocks,snp_blocks,floors_array,
			       splicesites,nsplicesites,subopt_levels,
			       indel_penalty,max_middle_insertions,max_middle_deletions,
			       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			       dibasep,cmetp,fast_level_5,
			       /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
      }
    }
  }

  /* 4/5B.  Complete set mismatches and indels, omitting frequent oligos */
  if (dibasep || done_level_3 > fast_level_3 || done_level_3 >= indel_penalty) {
    if (masktype == MASK_NONE) {
      debug(printf("*** Stage 4B,5B.  Complete mm/indels with no masking with done_level %d ***\n",done_level_3));
      complete_set_mm_indels(&ignore_found_score,&segments3_computed_p,&plus_segments_3,&plus_nsegments_3,&minus_segments_3,&minus_nsegments_3,
			     &any_omitted_p_3,&opt_level,&done_level_3,user_maxlevel_3,/*revise_levels_p*/false,
			     &subs3,&indels3,this3,queryuc_ptr_3,queryrc3,
			     querylength3,query3_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,floors_array,
			     splicesites,nsplicesites,subopt_levels,
			     indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     dibasep,cmetp,fast_level_3,
			     /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
    } else {
      debug(printf("*** Stage 4B,5B.  Complete mm/indels masking frequent oligos with done_level %d ***\n",done_level_3));
      complete_set_mm_indels(&ignore_found_score,&segments3_computed_p,&plus_segments_3,&plus_nsegments_3,&minus_segments_3,&minus_nsegments_3,
			     &any_omitted_p_3,&opt_level,&done_level_3,user_maxlevel_3,/*revise_levels_p*/false,
			     &subs3,&indels3,this3,queryuc_ptr_3,queryrc3,
			     querylength3,query3_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			     chromosome_iit,genome_blocks,snp_blocks,floors_array,
			     splicesites,nsplicesites,subopt_levels,
			     indel_penalty,max_middle_insertions,max_middle_deletions,
			     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			     dibasep,cmetp,fast_level_3,/*omit_frequent_p*/true,
			     /*omit_repetitive_p*/(masktype == MASK_REPETITIVE || masktype == MASK_GREEDY_REPETITIVE) ? true : false
			     );
      if ((masktype == MASK_GREEDY_FREQUENT || masktype == MASK_GREEDY_REPETITIVE) && subs3 == NULL && indels3 == NULL && any_omitted_p_3 == true) {
	FREE(minus_segments_3);
	FREE(plus_segments_3);

	/* 4/5B.  Complete set mismatches and indels, with all oligos */
	debug(printf("*** Stage 4B,5B.  Complete mm/indels with no masking with done_level %d ***\n",done_level_3));
	complete_set_mm_indels(&ignore_found_score,&segments3_computed_p,&plus_segments_3,&plus_nsegments_3,&minus_segments_3,&minus_nsegments_3,
			       &any_omitted_p_3,&opt_level,&done_level_3,user_maxlevel_3,/*revise_levels_p*/false,
			       &subs3,&indels3,this3,queryuc_ptr_3,queryrc3,
			       querylength3,query3_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			       chromosome_iit,genome_blocks,snp_blocks,floors_array,
			       splicesites,nsplicesites,subopt_levels,
			       indel_penalty,max_middle_insertions,max_middle_deletions,
			       allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			       dibasep,cmetp,fast_level_3,
			       /*omit_frequent_p*/false,/*omit_repetitive_p*/false);
      }
    }
  }


  /* 4/5. Pairing after complete set subs and indels */
  debug(printf("Starting pairing of 4 and 5\n"));
  hitarray5[HITARRAY_SUBS] = subs5 = Stage3_remove_duplicates(subs5);
  hitarray5[HITARRAY_INDELS] = indels5 = Stage3_remove_duplicates(indels5);
  hitarray3[HITARRAY_SUBS] = subs3 = Stage3_remove_duplicates(subs3);
  hitarray3[HITARRAY_INDELS] = indels3 = Stage3_remove_duplicates(indels3);
  hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*samechr),hitpairs,
				       hitarray5,/*narray5*/HITARRAY_INDELS+1,
				       hitarray3,/*narray3*/HITARRAY_INDELS+1,
				       *cutoff_level_5,*cutoff_level_3,subopt_levels,
				       splicesites,queryuc_ptr_5,queryuc_ptr_3,this5->query_compress_fwd,this5->query_compress_rev,
				       this3->query_compress_fwd,this3->query_compress_rev,genome_blocks,snp_blocks,
				       pairmax,expected_pairlength,pairlength_deviation,
				       querylength5,querylength3,maxpairedpaths,
				       /*allow_concordant_translocations_p*/true,localsplicing_penalty,dibasep,cmetp);
  debug(printf("After pairing complete set mismatches and indels, found %d concordant, found_score %d\n",*nconcordant,found_score));
  if (*abort_pairing_p == true) {
    *hits5 = List_append(subs5,indels5);
    *hits3 = List_append(subs3,indels3);
    return hitpairs;
  } else {
    opt_level = (found_score < opt_level) ? found_score : opt_level;
    if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
      done_level_5 = user_maxlevel_5;
    }
    if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
      done_level_3 = user_maxlevel_3;
    }
    debug(printf("4/5> found_score = %d, opt_level %d, done_level %d,%d\n",found_score,opt_level,done_level_5,done_level_3));
  }


  /* 6/7/8/9. Known splicing.  Requires compress and all positions fetched. */
  /* Subtract 1 from done_level for previous hits */
  if (knownsplicingp || novelsplicingp) {
    debug(printf("Deciding whether to do singlesplicing: done_level_5 %d >=? localsplicing_penalty %d\n",
		 done_level_5,localsplicing_penalty));

    if (done_level_5 >= localsplicing_penalty) {
      debug(printf("*** Stage 6A.  Single splicing masking frequent oligos with done_level %d ***\n",done_level_5));
      /* Always mask frequent oligos for splicing, which must be transcriptional */
      floors5 = compute_floors(&any_omitted_p_5,&alloc_floors_p_5,floors_array,this5,
			       querylength5,query5_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			       max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true);
      floors5_computed_p = true;

      if (segments5_computed_p == false) {
	plus_segments_5 = identify_all_segments(&plus_nsegments_5,this5->plus_positions,this5->plus_npositions,
						this5->omitted,querylength5,query5_lastpos,
						chromosome_iit,floors5,splicesites,nsplicesites,/*plusp*/true);
	minus_segments_5 = identify_all_segments(&minus_nsegments_5,this5->minus_positions,this5->minus_npositions,
						 this5->omitted,querylength5,query5_lastpos,
						 chromosome_iit,floors5,splicesites,nsplicesites,/*plusp*/false);
	segments5_computed_p = true;
      }

      singlesplicing5 = complete_set_singlesplicing(&ignore_found_score,floors5,
						    plus_segments_5,plus_nsegments_5,
						    minus_segments_5,minus_nsegments_5,
						    singlesplicing5,this5,queryuc_ptr_5,queryrc5,
						    querylength5,query5_lastpos,chromosome_iit,
						    genome_blocks,snp_blocks,
						    splicesites,splicetypes,splicedists,nsplicesites,
						    shortsplicedist,novelsplicingp,
						    localsplicing_penalty,min_localsplicing_end_matches,
						    /*max_mismatches_allowed*/done_level_5 - localsplicing_penalty,
						    /*first_read_p*/true,dibasep,cmetp);
    }

    debug(printf("Deciding whether to do singlesplicing: done_level_3 %d >=? localsplicing_penalty %d\n",
		 done_level_3,localsplicing_penalty));
    if (done_level_3 >= localsplicing_penalty) {
      debug(printf("*** Stage 6B.  Single splicing masking frequent oligos with done_level %d ***\n",done_level_3));
      /* Always mask frequent oligos for splicing, which must be transcriptional */
      floors3 = compute_floors(&any_omitted_p_3,&alloc_floors_p_3,floors_array,this3,
			       querylength3,query3_lastpos,indexdb,indexdb2,indexdb_size_threshold,
			       max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true);
      floors3_computed_p = true;

      if (segments3_computed_p == false) {
	plus_segments_3 = identify_all_segments(&plus_nsegments_3,this3->plus_positions,this3->plus_npositions,
						this3->omitted,querylength3,query3_lastpos,
						chromosome_iit,floors3,splicesites,nsplicesites,/*plusp*/true);
	minus_segments_3 = identify_all_segments(&minus_nsegments_3,this3->minus_positions,this3->minus_npositions,
						 this3->omitted,querylength3,query3_lastpos,
						 chromosome_iit,floors3,splicesites,nsplicesites,/*plusp*/false);
	segments3_computed_p = true;
      }

      singlesplicing3 = complete_set_singlesplicing(&ignore_found_score,floors3,
						    plus_segments_3,plus_nsegments_3,
						    minus_segments_3,minus_nsegments_3,
						    singlesplicing3,this3,queryuc_ptr_3,queryrc3,
						    querylength3,query3_lastpos,chromosome_iit,
						    genome_blocks,snp_blocks,
						    splicesites,splicetypes,splicedists,nsplicesites,
						    shortsplicedist,novelsplicingp,
						    localsplicing_penalty,min_localsplicing_end_matches,
						    /*max_mismatches_allowed*/done_level_3 - localsplicing_penalty,
						    /*first_read_p*/false,dibasep,cmetp);
    }

    /* 6.  Pairing after single splicing */
    /* Mark ambiguous splices only for single-end reads */
    hitarray5[HITARRAY_SINGLESPLICING] = singlesplicing5;
    hitarray3[HITARRAY_SINGLESPLICING] = singlesplicing3;

    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*samechr),hitpairs,
					 hitarray5,/*narray5*/HITARRAY_SINGLESPLICING+1,
					 hitarray3,/*narray3*/HITARRAY_SINGLESPLICING+1,
					 *cutoff_level_5,*cutoff_level_3,subopt_levels,
					 splicesites,queryuc_ptr_5,queryuc_ptr_3,this5->query_compress_fwd,this5->query_compress_rev,
					 this3->query_compress_fwd,this3->query_compress_rev,genome_blocks,snp_blocks,
					 pairmax,expected_pairlength,pairlength_deviation,
					 querylength5,querylength3,maxpairedpaths,
					 /*allow_concordant_translocations_p*/true,localsplicing_penalty,dibasep,cmetp);
    debug(printf("After pairing single splicing, found %d concordant, found_score %d\n",*nconcordant,found_score));
    if (*abort_pairing_p == true) {
      FREE(minus_segments_3);
      FREE(plus_segments_3);
      FREE(minus_segments_5);
      FREE(plus_segments_5);
      if (alloc_floors_p_5 == true) {
	Floors_free(&floors5);
      }
      if (alloc_floors_p_3 == true) {
	Floors_free(&floors3);
      }
      *hits5 = List_append(subs5,List_append(indels5,singlesplicing5));
      *hits3 = List_append(subs3,List_append(indels3,singlesplicing3));
      return hitpairs;

    } else {
      opt_level = (found_score < opt_level) ? found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("Pairing after 6A and 6B> found_score = %d, opt_level %d, done_level %d,%d\n",
		   found_score,opt_level,done_level_5,done_level_3));
    }


    /* 7.  Double splicing */
    if (done_level_5 >= localsplicing_penalty) {
      debug(printf("*** Stage 7A.  Double splicing masking frequent oligos with done_level %d ***\n",done_level_5));
      doublesplicing5 = complete_set_doublesplicing(&ignore_found_score,floors5,
						    plus_segments_5,plus_nsegments_5,
						    minus_segments_5,minus_nsegments_5,
						    doublesplicing5,this5,queryuc_ptr_5,queryrc5,
						    querylength5,query5_lastpos,chromosome_iit,
						    genome_blocks,snp_blocks,
						    splicesites,splicetypes,splicefrags_ref,splicefrags_alt,
						    nsplicesites,nsplicepartners_skip,
						    trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
						    trieoffsets_max,triecontents_max,nsplicepartners_max,
						    splicestrings,shortsplicedist,knownsplicingp,novelsplicingp,
						    localsplicing_penalty,min_localsplicing_end_matches,min_shortend,
						    find_novel_doublesplices_p,/*max_mismatches_allowed*/done_level_5 - localsplicing_penalty,
						    /*first_read_p*/true,dibasep,cmetp);
    }

    if (done_level_3 >= localsplicing_penalty) {
      debug(printf("*** Stage 7B.  Double splicing masking frequent oligos with done_level %d ***\n",done_level_3));
      doublesplicing3 = complete_set_doublesplicing(&ignore_found_score,floors3,
						    plus_segments_3,plus_nsegments_3,
						    minus_segments_3,minus_nsegments_3,
						    doublesplicing3,this3,queryuc_ptr_3,queryrc3,
						    querylength3,query3_lastpos,chromosome_iit,
						    genome_blocks,snp_blocks,
						    splicesites,splicetypes,splicefrags_ref,splicefrags_alt,
						    nsplicesites,nsplicepartners_skip,
						    trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
						    trieoffsets_max,triecontents_max,nsplicepartners_max,
						    splicestrings,shortsplicedist,knownsplicingp,novelsplicingp,
						    localsplicing_penalty,min_localsplicing_end_matches,min_shortend,
						    find_novel_doublesplices_p,/*max_mismatches_allowed*/done_level_3 - localsplicing_penalty,
						    /*first_read_p*/false,dibasep,cmetp);
    }
    
    /* 7.  Pairing after double splicing */
    /* Mark ambiguous splices only for single-end reads */
    hitarray5[HITARRAY_DOUBLESPLICING] = doublesplicing5;
    hitarray3[HITARRAY_DOUBLESPLICING] = doublesplicing3;

    hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*samechr),hitpairs,
					 hitarray5,/*narray5*/HITARRAY_DOUBLESPLICING+1,
					 hitarray3,/*narray3*/HITARRAY_DOUBLESPLICING+1,
					 *cutoff_level_5,*cutoff_level_3,subopt_levels,
					 splicesites,queryuc_ptr_5,queryuc_ptr_3,this5->query_compress_fwd,this5->query_compress_rev,
					 this3->query_compress_fwd,this3->query_compress_rev,genome_blocks,snp_blocks,
					 pairmax,expected_pairlength,pairlength_deviation,
					 querylength5,querylength3,maxpairedpaths,
					 /*allow_concordant_translocations_p*/true,localsplicing_penalty,dibasep,cmetp);
    debug(printf("After pairing double splicing, found %d concordant, found_score %d\n",*nconcordant,found_score));
    if (*abort_pairing_p == true) {
      FREE(minus_segments_3);
      FREE(plus_segments_3);
      FREE(minus_segments_5);
      FREE(plus_segments_5);
      if (alloc_floors_p_5 == true) {
	Floors_free(&floors5);
      }
      if (alloc_floors_p_3 == true) {
	Floors_free(&floors3);
      }
      *hits5 = List_append(subs5,List_append(indels5,List_append(singlesplicing5,doublesplicing5)));
      *hits3 = List_append(subs3,List_append(indels3,List_append(singlesplicing3,doublesplicing3)));
      return hitpairs;

    } else {
      opt_level = (found_score < opt_level) ? found_score : opt_level;
      if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	done_level_5 = user_maxlevel_5;
      }
      if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	done_level_3 = user_maxlevel_3;
      }
      debug(printf("Pairing after 7A and 7B> found_score = %d, opt_level %d, done_level %d,%d\n",
		   found_score,opt_level,done_level_5,done_level_3));
    }


    /* 8A.  Short-Overlap splicing */
    if (knownsplicingp == true && localsplicing_penalty < distantsplicing_penalty) {
      max_splice_mismatches_5 = done_level_5 - localsplicing_penalty;
    } else {
      max_splice_mismatches_5 = done_level_5 - distantsplicing_penalty;
    }

    if (max_splice_mismatches_5 >= 0) {
      alloc5p = true;
      donors_plus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      antidonors_plus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      acceptors_plus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      antiacceptors_plus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      donors_minus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      antidonors_minus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      acceptors_minus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));
      antiacceptors_minus_5 = (List_T *) CALLOC(max_splice_mismatches_5+1,sizeof(List_T));

      find_spliceends(&donors_plus_5,&antidonors_plus_5,&acceptors_plus_5,&antiacceptors_plus_5,
		      plus_segments_5,plus_nsegments_5,genome_blocks,snp_blocks,
		      queryuc_ptr_5,queryuc_ptr_5,floors5,querylength5,query5_lastpos,
		      /*query_compress*/this5->query_compress_fwd,
		      splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
		      /*max_mismatches_allowed*/max_splice_mismatches_5,
		      dibasep,cmetp,/*plusp*/true);

      find_spliceends(&antidonors_minus_5,&donors_minus_5,&antiacceptors_minus_5,&acceptors_minus_5,
		      minus_segments_5,minus_nsegments_5,genome_blocks,snp_blocks,
		      /*query*/queryuc_ptr_5,/*queryptr*/queryrc5,
		      floors5,querylength5,query5_lastpos,
		      /*query_compress*/this5->query_compress_rev,
		      splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
		      /*max_mismatches_allowed*/max_splice_mismatches_5,
		      dibasep,cmetp,/*plusp*/false);

      if (knownsplicingp == true) {
	/* 8A.  Find short-overlap splicing using known splice sites */
	shortends5 = find_splicepairs_short_overlaps(&ignore_found_score,/*hits*/(List_T) NULL,
						     donors_plus_5,antidonors_plus_5,acceptors_plus_5,antiacceptors_plus_5,
						     donors_minus_5,antidonors_minus_5,acceptors_minus_5,antiacceptors_minus_5,
						     this5->query_compress_fwd,this5->query_compress_rev,
						     genome_blocks,snp_blocks,queryuc_ptr_5,queryrc5,
						     min_shortend,
						     splicesites,splicetypes,splicefrags_ref,splicefrags_alt,
						     nsplicesites,nsplicepartners_skip,
						     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
						     trieoffsets_max,triecontents_max,nsplicepartners_max,
						     splicestrings,shortsplicedist,localsplicing_penalty,min_localsplicing_end_matches,
						     /*max_mismatches_allowed*/max_splice_mismatches_5,querylength5,
						     /*first_read_p*/true,dibasep,cmetp);
      }
    }


    /* 8B.  Short-Overlap splicing */
    if (knownsplicingp == true && localsplicing_penalty < distantsplicing_penalty) {
      max_splice_mismatches_3 = done_level_3 - localsplicing_penalty;
    } else {
      max_splice_mismatches_3 = done_level_3 - distantsplicing_penalty;
    }

    if (max_splice_mismatches_3 >= 0) {
      alloc3p = true;
      donors_plus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      antidonors_plus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      acceptors_plus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      antiacceptors_plus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      donors_minus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      antidonors_minus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      acceptors_minus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));
      antiacceptors_minus_3 = (List_T *) CALLOC(max_splice_mismatches_3+1,sizeof(List_T));

      find_spliceends(&donors_plus_3,&antidonors_plus_3,&acceptors_plus_3,&antiacceptors_plus_3,
		      plus_segments_3,plus_nsegments_3,genome_blocks,snp_blocks,
		      queryuc_ptr_3,queryuc_ptr_3,floors3,querylength3,query3_lastpos,
		      /*query_compress*/this3->query_compress_fwd,
		      splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
		      /*max_mismatches_allowed*/max_splice_mismatches_3,
		      dibasep,cmetp,/*plusp*/true);

      find_spliceends(&antidonors_minus_3,&donors_minus_3,&antiacceptors_minus_3,&acceptors_minus_3,
		      minus_segments_3,minus_nsegments_3,genome_blocks,snp_blocks,
		      /*query*/queryuc_ptr_3,/*queryptr*/queryrc3,
		      floors3,querylength3,query3_lastpos,
		      /*query_compress*/this3->query_compress_rev,
		      splicesites,splicetypes,nsplicesites,novelsplicingp,canonicalp,
		      /*max_mismatches_allowed*/max_splice_mismatches_3,
		      dibasep,cmetp,/*plusp*/false);

      if (knownsplicingp == true) {
	/* 8B.  Find short-overlap splicing using known splice sites */
	shortends3 = find_splicepairs_short_overlaps(&ignore_found_score,/*hits*/(List_T) NULL,
						     donors_plus_3,antidonors_plus_3,acceptors_plus_3,antiacceptors_plus_3,
						     donors_minus_3,antidonors_minus_3,acceptors_minus_3,antiacceptors_minus_3,
						     this3->query_compress_fwd,this3->query_compress_rev,
						     genome_blocks,snp_blocks,queryuc_ptr_3,queryrc3,min_shortend,
						     splicesites,splicetypes,splicefrags_ref,splicefrags_alt,
						     nsplicesites,nsplicepartners_skip,
						     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
						     trieoffsets_max,triecontents_max,nsplicepartners_max,
						     splicestrings,shortsplicedist,localsplicing_penalty,min_localsplicing_end_matches,
						     /*max_mismatches_allowed*/max_splice_mismatches_3,querylength3,
						     /*first_read_p*/false,dibasep,cmetp);
      }
    }

    if (knownsplicingp == true) {
      /* 8.  Pairing after short-overlaps */
      hitarray5[HITARRAY_SHORTENDSPLICING] = shortends5 = Stage3_remove_duplicates(shortends5);
      hitarray3[HITARRAY_SHORTENDSPLICING] = shortends3 = Stage3_remove_duplicates(shortends3);
      hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*samechr),hitpairs,
					   hitarray5,/*narray5*/HITARRAY_SHORTENDSPLICING+1,
					   hitarray3,/*narray3*/HITARRAY_SHORTENDSPLICING+1,
					   *cutoff_level_5,*cutoff_level_3,subopt_levels,
					   splicesites,queryuc_ptr_5,queryuc_ptr_3,this5->query_compress_fwd,this5->query_compress_rev,
					   this3->query_compress_fwd,this3->query_compress_rev,genome_blocks,snp_blocks,
					   pairmax,expected_pairlength,pairlength_deviation,
					   querylength5,querylength3,maxpairedpaths,
					   /*allow_concordant_translocations_p*/true,localsplicing_penalty,dibasep,cmetp);
      debug(printf("After pairing short-overlap splicing, found %d concordant, found_score %d\n",*nconcordant,found_score));
      if (*abort_pairing_p == false) {
	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	  done_level_5 = user_maxlevel_5;
	}
	if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	  done_level_3 = user_maxlevel_3;
	}
	debug(printf("Pairing after 8A and 8B> found_score = %d, opt_level %d, done_level %d,%d\n",
		     found_score,opt_level,done_level_5,done_level_3));
      }
    }

    if (alloc5p == true) {
      if (*abort_pairing_p == false) {
	/* 9A.  Distant splicing */
	nmismatches = 0;
	while (nmismatches <= done_level_5 - distantsplicing_penalty && nsplicepairs5 < maxchimerapaths) {
	  debug(printf("*** Stage 9A.  Distant splicing, allowing %d mismatches ***\n",nmismatches));

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

	  distantsplicing5 = find_splicepairs_distant(&ignore_found_score,&nsplicepairs5,distantsplicing5,
						      donors_plus_5,antidonors_plus_5,acceptors_plus_5,antiacceptors_plus_5,
						      donors_minus_5,antidonors_minus_5,acceptors_minus_5,antiacceptors_minus_5,
						      shortsplicedist,localsplicing_penalty,distantsplicing_penalty,
						      min_distantsplicing_end_matches,min_distantsplicing_identity,
						      querylength5,nmismatches,maxchimerapaths,/*first_read_p*/true);
	  nmismatches++;
	}
      }

      /* Clean up 5 */
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

    if (alloc3p == true) {
      if (*abort_pairing_p == false) {
	/* 9B.  Distant splicing */
	nmismatches = 0;
	while (nmismatches <= done_level_3 - distantsplicing_penalty && nsplicepairs3 < maxchimerapaths) {
	  debug(printf("*** Stage 9B.  Distant splicing, allowing %d mismatches ***\n",nmismatches));

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

	  distantsplicing3 = find_splicepairs_distant(&ignore_found_score,&nsplicepairs3,distantsplicing3,
						      donors_plus_3,antidonors_plus_3,acceptors_plus_3,antiacceptors_plus_3,
						      donors_minus_3,antidonors_minus_3,acceptors_minus_3,antiacceptors_minus_3,
						      shortsplicedist,localsplicing_penalty,distantsplicing_penalty,
						      min_distantsplicing_end_matches,min_distantsplicing_identity,
						      querylength3,nmismatches,maxchimerapaths,/*first_read_p*/false);
	  nmismatches++;
	}
      }

      /* Clean up 3 */
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

    /* 9.  Pairing after distant splicing */
#if 0
    /* Mark ambiguous splices only for single-end reads */
    hitarray5[HITARRAY_DISTANTSPLICING] = Stage3_mark_ambiguous_splices(&ambiguousp,distantsplicing5);
    hitarray3[HITARRAY_DISTANTSPLICING] = Stage3_mark_ambiguous_splices(&ambiguousp,distantsplicing3);
#else
    hitarray5[HITARRAY_DISTANTSPLICING] = distantsplicing5;
    hitarray3[HITARRAY_DISTANTSPLICING] = distantsplicing3;
#endif

    if (*nconcordant == 0) {
      hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*samechr),hitpairs,
					   hitarray5,/*narray5*/HITARRAY_DISTANTSPLICING+1,
					   hitarray3,/*narray3*/HITARRAY_DISTANTSPLICING+1,
					   *cutoff_level_5,*cutoff_level_3,subopt_levels,
					   splicesites,queryuc_ptr_5,queryuc_ptr_3,this5->query_compress_fwd,this5->query_compress_rev,
					   this3->query_compress_fwd,this3->query_compress_rev,genome_blocks,snp_blocks,
					   pairmax,expected_pairlength,pairlength_deviation,
					   querylength5,querylength3,maxpairedpaths,
					   /*allow_concordant_translocations_p*/true,localsplicing_penalty,dibasep,cmetp);
      debug(printf("After pairing distant splicing, found %d concordant, found_score %d\n",*nconcordant,found_score));

      if (*abort_pairing_p == false) {
	opt_level = (found_score < opt_level) ? found_score : opt_level;
	if ((done_level_5 = opt_level + subopt_levels) > user_maxlevel_5) {
	  done_level_5 = user_maxlevel_5;
	}
	if ((done_level_3 = opt_level + subopt_levels) > user_maxlevel_3) {
	  done_level_3 = user_maxlevel_3;
	}
	debug(printf("9> found_score = %d, opt_level %d, done_level %d,%d\n",found_score,opt_level,done_level_5,done_level_3));
      }
    }

  }

  /* 10A,B.  Terminals */
  if (/* *nconcordant == 0 && */ *abort_pairing_p == false) {
    debug(printf("Stage 10A.  Finding terminals5, done_level_5 = %d, terminal_penalty = %d\n",
		 done_level_5,terminal_penalty));
    if (done_level_5 >= terminal_penalty) {
      if (floors5_computed_p == false) {
	floors5 = compute_floors(&any_omitted_p_5,&alloc_floors_p_5,floors_array,this5,
				 querylength5,query5_lastpos,indexdb,indexdb2,indexdb_size_threshold,
				 max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true);
      }

      if (segments5_computed_p == false) {
	plus_segments_5 = identify_all_segments_for_terminals(&plus_nsegments_5,this5->plus_positions,this5->plus_npositions,
							      this5->omitted,querylength5,query5_lastpos,
							      chromosome_iit,floors5,
							      /*max_mismatches_allowed*/done_level_5 - terminal_penalty,
							      /*plusp*/true);
	minus_segments_5 = identify_all_segments_for_terminals(&minus_nsegments_5,this5->minus_positions,this5->minus_npositions,
							       this5->omitted,querylength5,query5_lastpos,
							       chromosome_iit,floors5,
							       /*max_mismatches_allowed*/done_level_5 - terminal_penalty,
							       /*plusp*/false);
      }

      /* Don't run Stage3_remove_duplicates until after concordant pairs are found */
      if (this5->query_compress_fwd == NULL) {
	this5->query_compress_fwd = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/true);
	this5->query_compress_rev = Compress_new(queryuc_ptr_5,querylength5,/*plusp*/false);
      }
      hitarray5[HITARRAY_TERMINALS] = terminals5 =
	find_terminals(plus_segments_5,plus_nsegments_5,minus_segments_5,minus_nsegments_5,
		       genome_blocks,snp_blocks,queryuc_ptr_5,queryrc5,
		       floors5,querylength5,query5_lastpos,
		       this5->query_compress_fwd,this5->query_compress_rev,
		       /*max_mismatches_allowed*/done_level_5 - terminal_penalty,
		       terminal_penalty,max_terminal_length,dibasep,cmetp);
    }

    debug(printf("Stage 10B.  Finding terminals3, done_level_3 = %d, terminal_penalty = %d\n",
		 done_level_3,terminal_penalty));
    if (done_level_3 >= terminal_penalty) {
      if (floors3_computed_p == false) {
	floors3 = compute_floors(&any_omitted_p_3,&alloc_floors_p_3,floors_array,this3,
				 querylength3,query3_lastpos,indexdb,indexdb2,indexdb_size_threshold,
				 max_end_insertions,/*omit_frequent_p*/true,/*omit_repetitive_p*/true);
      }
      if (segments3_computed_p == false) {
	plus_segments_3 = identify_all_segments_for_terminals(&plus_nsegments_3,this3->plus_positions,this3->plus_npositions,
							      this3->omitted,querylength3,query3_lastpos,
							      chromosome_iit,floors3,
							      /*max_mismatches_allowed*/done_level_3 - terminal_penalty,
							      /*plusp*/true);
	minus_segments_3 = identify_all_segments_for_terminals(&minus_nsegments_3,this3->minus_positions,this3->minus_npositions,
							       this3->omitted,querylength3,query3_lastpos,
							       chromosome_iit,floors3,
							       /*max_mismatches_allowed*/done_level_3 - terminal_penalty,
							       /*plusp*/false);
      }

      /* Don't run Stage3_remove_duplicates until after concordant pairs are found */
      if (this3->query_compress_fwd == NULL) {
	this3->query_compress_fwd = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/true);
	this3->query_compress_rev = Compress_new(queryuc_ptr_3,querylength3,/*plusp*/false);
      }
      hitarray3[HITARRAY_TERMINALS] = terminals3 =
	find_terminals(plus_segments_3,plus_nsegments_3,minus_segments_3,minus_nsegments_3,
		       genome_blocks,snp_blocks,queryuc_ptr_3,queryrc3,
		       floors3,querylength3,query3_lastpos,
		       this3->query_compress_fwd,this3->query_compress_rev,
		       /*max_mismatches_allowed*/done_level_3 - terminal_penalty,
		       terminal_penalty,max_terminal_length,dibasep,cmetp);
    }

    if (done_level_5 >= terminal_penalty || done_level_3 >= terminal_penalty) {
      hitpairs = Stage3_pair_up_concordant(&(*abort_pairing_p),&found_score,&(*nconcordant),&(*samechr),hitpairs,
					   hitarray5,/*narray5*/HITARRAY_TERMINALS+1,
					   hitarray3,/*narray3*/HITARRAY_TERMINALS+1,
					   *cutoff_level_5,*cutoff_level_3,subopt_levels,
					   splicesites,queryuc_ptr_5,queryuc_ptr_3,this5->query_compress_fwd,this5->query_compress_rev,
					   this3->query_compress_fwd,this3->query_compress_rev,genome_blocks,snp_blocks,
					   pairmax,expected_pairlength,pairlength_deviation,
					   querylength5,querylength3,maxpairedpaths,
					   /*allow_concordant_translocations_p*/false,localsplicing_penalty,dibasep,cmetp);
      debug(printf("10> After pairing terminals, found %d concordant, found_score %d\n",*nconcordant,found_score));
    }
  }

  if (alloc_floors_p_5 == true) {
    Floors_free(&floors5);
  }
  if (alloc_floors_p_3 == true) {
    Floors_free(&floors3);
  }

  FREE(minus_segments_5);
  FREE(plus_segments_5);

  FREE(minus_segments_3);
  FREE(plus_segments_3);


  *hits5 = List_append(subs5,List_append(indels5,List_append(singlesplicing5,List_append(doublesplicing5,List_append(shortends5,List_append(distantsplicing5,terminals5))))));
  *hits3 = List_append(subs3,List_append(indels3,List_append(singlesplicing3,List_append(doublesplicing3,List_append(shortends3,List_append(distantsplicing3,terminals3))))));

  debug(printf("Ending with %d hitpairs\n",List_length(hitpairs)));

  return hitpairs;
}


Stage3pair_T *
Stage1_paired_read (int *npaths, bool *concordantp, Stage3_T **stage3array5, int *nhits5, Stage3_T **stage3array3, int *nhits3,
		    Shortread_T queryseq5, Shortread_T queryseq3,
		    Indexdb_T indexdb, Indexdb_T indexdb2, int indexdb_size_threshold,
		    IIT_T chromosome_iit, Genome_T genome, Genome_T genomealt, Floors_T *floors_array,
		    bool knownsplicingp, bool novelsplicingp, bool canonicalp,
		    int maxpaths, int maxchimerapaths, double user_maxlevel_float, int subopt_levels,
		    Masktype_T masktype, int terminal_penalty, int max_terminal_length,
		    int indel_penalty, int max_middle_insertions, int max_middle_deletions,
		    bool allow_end_indels_p, int max_end_insertions, int max_end_deletions, int min_indel_end_matches,
		    Genomicpos_T shortsplicedist,
		    int localsplicing_penalty, int distantsplicing_penalty, int min_localsplicing_end_matches,
		    int min_distantsplicing_end_matches, double min_distantsplicing_identity, int min_shortend,
		    bool find_novel_doublesplices_p, Genomicpos_T *splicesites, Splicetype_T *splicetypes,
		    Genomicpos_T *splicedists, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int nsplicesites,
		    int *nsplicepartners_skip,
		    unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
		    unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
		    List_T *splicestrings, Genomicpos_T pairmax, Genomicpos_T expected_pairlength,
		    Genomicpos_T pairlength_deviation, bool dibasep, bool cmetp) {
  Stage3pair_T *stage3pairarray, stage3pair;
  List_T result = NULL, hitpairs = NULL, samechr = NULL, hits5 = NULL, hits3 = NULL, p;
  List_T singlehits5 = NULL, singlehits3 = NULL;
  Stage3_T hit5, hit3;
  Pairtype_T pairtype;
  T this5, this3;
  char *queryuc_ptr_5, *queryuc_ptr_3, *quality_string_5, *quality_string_3;
  int user_maxlevel_5, user_maxlevel_3;
  int nconcordant, nsamechr;
  int cutoff_level_5, cutoff_level_3;
  int querylength5, querylength3, query5_lastpos, query3_lastpos;
  bool allvalidp5, allvalidp3;
  UINT4 *genome_blocks, *snp_blocks;
#if 0
  int maxpairedpaths = 10*maxpaths; /* For computation, not for printing. */
#else
  int maxpairedpaths = 100000;
#endif
  bool abort_pairing_p;

  querylength5 = Shortread_fulllength(queryseq5);
  querylength3 = Shortread_fulllength(queryseq3);

  if (querylength5 <= INDEX1PART+2 || querylength3 <= INDEX1PART+2) {
    fprintf(stderr,"GSNAP cannot handle reads shorter than %d bp\n",INDEX1PART+2);
    *npaths = *nhits5 = *nhits3 = 0;
    *stage3array5 = *stage3array3 = (Stage3_T *) NULL;
    return (Stage3pair_T *) NULL;

  } else if (querylength5 > MAX_QUERYLENGTH || querylength3 > MAX_QUERYLENGTH) {
    fprintf(stderr,"GSNAP cannot handle reads longer than %d bp.  Either recompile with a higher value of MAX_QUERYLENGTH, or consider using GMAP instead.\n",
	    MAX_QUERYLENGTH);
    *npaths = *nhits5 = *nhits3 = 0;
    *stage3array5 = *stage3array3 = (Stage3_T *) NULL;
    return (Stage3pair_T *) NULL;

  } else {
    if (user_maxlevel_float < 0.0) {
      user_maxlevel_5 = user_maxlevel_3 = -1;
    } else if (user_maxlevel_float > 0.0 && user_maxlevel_float < 1.0) {
      user_maxlevel_5 = (int) rint(user_maxlevel_float * (double) querylength5);
      user_maxlevel_3 = (int) rint(user_maxlevel_float * (double) querylength3);
    } else {
      user_maxlevel_5 = user_maxlevel_3 = (int) user_maxlevel_float;
    }

    this5 = Stage1_new(querylength5);
    this3 = Stage1_new(querylength3);
    queryuc_ptr_5 = Shortread_fullpointer_uc(queryseq5);
    queryuc_ptr_3 = Shortread_fullpointer_uc(queryseq3);
    quality_string_5 = Shortread_quality_string(queryseq5);
    quality_string_3 = Shortread_quality_string(queryseq3);
    query5_lastpos = querylength5 - INDEX1PART;
    query3_lastpos = querylength3 - INDEX1PART;

    if (read_oligos(&allvalidp5,this5,queryuc_ptr_5,querylength5,query5_lastpos,dibasep,cmetp) == 0 ||
	read_oligos(&allvalidp3,this3,queryuc_ptr_3,querylength3,query3_lastpos,dibasep,cmetp) == 0) {
      debug(printf("Aborting because no hits found anywhere\n"));
      Stage1_free(&this3,querylength3);
      Stage1_free(&this5,querylength5);

      *npaths = *nhits5 = *nhits3 = 0;
      *stage3array5 = *stage3array3 = (Stage3_T *) NULL;
      return (Stage3pair_T *) NULL;

    } else {
      genome_blocks = Genome_blocks(genome);
      snp_blocks = genomealt ? Genome_blocks(genomealt) : NULL;
      abort_pairing_p = false;
      hitpairs = align_pair(&abort_pairing_p,&nconcordant,&cutoff_level_5,&cutoff_level_3,&samechr,&hits5,&hits3,
			    this5,this3,queryuc_ptr_5,queryuc_ptr_3,
			    querylength5,querylength3,query5_lastpos,query3_lastpos,
			    indexdb,indexdb2,indexdb_size_threshold,
			    chromosome_iit,genome_blocks,snp_blocks,floors_array,
			    splicesites,splicetypes,splicedists,splicefrags_ref,splicefrags_alt,
			    nsplicesites,nsplicepartners_skip,
			    trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
			    trieoffsets_max,triecontents_max,nsplicepartners_max,
			    splicestrings,user_maxlevel_5,user_maxlevel_3,subopt_levels,masktype,
			    terminal_penalty,max_terminal_length,indel_penalty,
			    maxchimerapaths,knownsplicingp,novelsplicingp,
			    canonicalp,shortsplicedist,
			    localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
			    min_distantsplicing_end_matches,min_distantsplicing_identity,min_shortend,
			    find_novel_doublesplices_p,max_middle_insertions,max_middle_deletions,
			    allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			    allvalidp5,allvalidp3,dibasep,cmetp,
			    pairmax,expected_pairlength,pairlength_deviation,maxpairedpaths);

      if (abort_pairing_p == true) {
	/* Clean up all previous calculations */
	for (p = hitpairs; p != NULL; p = List_next(p)) {
	  stage3pair = (Stage3pair_T) List_head(p);
	  Stage3pair_free(&stage3pair);
	}
	List_free(&hitpairs);

	for (p = samechr; p != NULL; p = List_next(p)) {
	  stage3pair = (Stage3pair_T) List_head(p);
	  Stage3pair_free(&stage3pair);
	}
	List_free(&samechr);

	stage3list_gc(&hits3);
	stage3list_gc(&hits5);
	Stage1_free(&this3,querylength3);
	Stage1_free(&this5,querylength5);


	/* Re-align 5' end as a single end */
	this5 = Stage1_new(querylength5);
	if (read_oligos(&allvalidp5,this5,queryuc_ptr_5,querylength5,query5_lastpos,dibasep,cmetp) == 0) {
	  debug(printf("Aborting because no hits found anywhere\n"));
	  singlehits5 = (List_T) NULL;
	} else {
	  singlehits5 = align_end(&cutoff_level_5,this5,queryuc_ptr_5,querylength5,query5_lastpos,
				  indexdb,indexdb2,indexdb_size_threshold,
				  chromosome_iit,genome_blocks,snp_blocks,floors_array,
				  splicesites,splicetypes,splicedists,splicefrags_ref,splicefrags_alt,
				  nsplicesites,nsplicepartners_skip,
				  trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				  trieoffsets_max,triecontents_max,nsplicepartners_max,
				  splicestrings,user_maxlevel_5,subopt_levels,masktype,
				  terminal_penalty,max_terminal_length,indel_penalty,
				  maxchimerapaths,knownsplicingp,novelsplicingp,
				  canonicalp,shortsplicedist,
				  localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
				  min_distantsplicing_end_matches,min_distantsplicing_identity,min_shortend,
				  find_novel_doublesplices_p,max_middle_insertions,max_middle_deletions,
				  allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
				  allvalidp5,dibasep,cmetp);
	  singlehits5 = Stage3_remove_duplicates(Stage3_optimal_score(singlehits5,cutoff_level_5,subopt_levels));
	}

	if ((*nhits5 = List_length(singlehits5)) == 0) {
	  *stage3array5 = (Stage3_T *) NULL;
	} else {
	  *stage3array5 = (Stage3_T *) List_to_array(singlehits5,NULL);
	  *stage3array5 = Stage3_eval_and_sort(*stage3array5,*nhits5,maxpaths,queryseq5,
					       this5->query_compress_fwd,this5->query_compress_rev,
					       genome_blocks,snp_blocks,genome,quality_string_5,
					       dibasep,cmetp);
	  List_free(&singlehits5);
	}

	/* Re-align 3' end as a single end */
	this3 = Stage1_new(querylength3);
	if (read_oligos(&allvalidp3,this3,queryuc_ptr_3,querylength3,query3_lastpos,dibasep,cmetp) == 0) {
	  debug(printf("Aborting because no hits found anywhere\n"));
	  singlehits3 = (List_T) NULL;
	} else {
	  singlehits3 = align_end(&cutoff_level_3,this3,queryuc_ptr_3,querylength3,query3_lastpos,
				  indexdb,indexdb2,indexdb_size_threshold,
				  chromosome_iit,genome_blocks,snp_blocks,floors_array,
				  splicesites,splicetypes,splicedists,splicefrags_ref,splicefrags_alt,
				  nsplicesites,nsplicepartners_skip,
				  trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				  trieoffsets_max,triecontents_max,nsplicepartners_max,
				  splicestrings,user_maxlevel_3,subopt_levels,masktype,
				  terminal_penalty,max_terminal_length,indel_penalty,
				  maxchimerapaths,knownsplicingp,novelsplicingp,
				  canonicalp,shortsplicedist,
				  localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
				  min_distantsplicing_end_matches,min_distantsplicing_identity,min_shortend,
				  find_novel_doublesplices_p,max_middle_insertions,max_middle_deletions,
				  allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
				  allvalidp3,dibasep,cmetp);
	  singlehits3 = Stage3_remove_duplicates(Stage3_optimal_score(singlehits3,cutoff_level_3,subopt_levels));
	}

	if ((*nhits3 = List_length(singlehits3)) == 0) {
	  *stage3array3 = (Stage3_T *) NULL;
	} else {
	  *stage3array3 = (Stage3_T *) List_to_array(singlehits3,NULL);
	  *stage3array3 = Stage3_eval_and_sort(*stage3array3,*nhits3,maxpaths,queryseq3,
					       this3->query_compress_fwd,this3->query_compress_rev,
					       genome_blocks,snp_blocks,genome,quality_string_3,
					       dibasep,cmetp);
	  List_free(&singlehits3);
	}

	*npaths = 0;
	*concordantp = false;
	Stage1_free(&this5,querylength5);
	Stage1_free(&this3,querylength3);
	return (Stage3pair_T *) NULL;

      } else {
	if (hitpairs == NULL) {
	  *concordantp = false;
	  if (samechr == NULL) {
	    /* No concordant or paired results */
	    result = (List_T) NULL;
	  } else {
	    /* No concordant results, but have paired results */
	    result = Stage3pair_remove_duplicates(Stage3pair_optimal_score(samechr,cutoff_level_5+cutoff_level_3,subopt_levels));
	  }

	} else {
	  /* Have concordant results */
	  *concordantp = true;
	  for (p = samechr; p != NULL; p = List_next(p)) {
	    stage3pair = (Stage3pair_T) List_head(p);
	    Stage3pair_free(&stage3pair);
	  }
	  List_free(&samechr);
	  
	  if (nconcordant > 0) {
	    if (novelsplicingp || knownsplicingp) {
	      hitpairs = Stage3pair_remove_excess_terminals(hitpairs);
	    }
	    result = Stage3pair_remove_duplicates(Stage3pair_optimal_score(hitpairs,cutoff_level_5+cutoff_level_3,subopt_levels));
	    /* result = Stage3pair_sort_distance(result); */
	    
	  } else {
	    /* nconcordant == 0, but hitpairs is not NULL, e.g., concordant translocations that don't increase nconcordant */
	    result = Stage3pair_remove_duplicates(Stage3pair_optimal_score(hitpairs,cutoff_level_5+cutoff_level_3,subopt_levels));
	    /* result = Stage3pair_sort_distance(result); */
	  }
	}

	if (result != NULL) {
	  /* Paired or concordant pairs found.  Remove single hits.
	     Value of *concordantp assigned above. */
	  stage3list_gc(&hits3);
	  stage3list_gc(&hits5);

	  *nhits5 = *nhits3 = 0;
	  *stage3array5 = *stage3array3 = (Stage3_T *) NULL;

	  if ((*npaths = List_length(result)) == 0) {
	    stage3pairarray = (Stage3pair_T *) NULL;
	  } else {
	    stage3pairarray = (Stage3pair_T *) List_to_array(result,NULL);
	    Stage3pair_eval(stage3pairarray,*npaths,maxpaths,queryseq5,queryseq3,
			    this5->query_compress_fwd,this5->query_compress_rev,
			    this3->query_compress_fwd,this3->query_compress_rev,
			    genome_blocks,snp_blocks,genome,quality_string_5,
			    quality_string_3,dibasep,cmetp);
	    List_free(&result);
	  }

	  Stage1_free(&this3,querylength3);
	  Stage1_free(&this5,querylength5);
	  return stage3pairarray;

	} else {
	  /* Unpaired */
	  singlehits5 = Stage3_remove_duplicates(Stage3_optimal_score(hits5,cutoff_level_5,subopt_levels));
	  singlehits3 = Stage3_remove_duplicates(Stage3_optimal_score(hits3,cutoff_level_3,subopt_levels));

	  if (List_length(singlehits5) == 1 && List_length(singlehits3) == 1) {
	    hit5 = (Stage3_T) List_head(singlehits5);
	    hit3 = (Stage3_T) List_head(singlehits3);
	    if ((pairtype = Stage3_determine_pairtype(hit5,hit3)) != UNPAIRED) {
	      /* Can convert unpaired uniq to a paired uniq */
	      stage3pairarray = (Stage3pair_T *) CALLOC(1,sizeof(Stage3pair_T));
	      stage3pairarray[0] = Stage3pair_new(hit5,hit3,queryuc_ptr_5,queryuc_ptr_3,
						  splicesites,this5->query_compress_fwd,this5->query_compress_rev,
						  this3->query_compress_fwd,this3->query_compress_rev,
						  genome_blocks,snp_blocks,expected_pairlength,pairlength_deviation,
						  pairtype,localsplicing_penalty,dibasep,cmetp);
	      stage3list_gc(&singlehits3);
	      stage3list_gc(&singlehits5);

	      *nhits5 = *nhits3 = 0;
	      *stage3array5 = *stage3array3 = (Stage3_T *) NULL;

	      *npaths = 1;
	      *concordantp = false;
	      Stage3pair_eval(stage3pairarray,/*npaths*/1,maxpaths,queryseq5,queryseq3,
			      this5->query_compress_fwd,this5->query_compress_rev,
			      this3->query_compress_fwd,this3->query_compress_rev,
			      genome_blocks,snp_blocks,genome,quality_string_5,
			      quality_string_3,dibasep,cmetp);

	      Stage1_free(&this3,querylength3);
	      Stage1_free(&this5,querylength5);
	      return stage3pairarray;
	    }
	  }

	  /* Unpaired hits */
	  *npaths = 0;
	  *concordantp = false;

	  if ((*nhits5 = List_length(singlehits5)) == 0) {
	    *stage3array5 = (Stage3_T *) NULL;
	  } else {
	    *stage3array5 = (Stage3_T *) List_to_array(singlehits5,NULL);
	    *stage3array5 = Stage3_eval_and_sort(*stage3array5,*nhits5,maxpaths,queryseq5,
						 this5->query_compress_fwd,this5->query_compress_rev,
						 genome_blocks,snp_blocks,genome,quality_string_5,
						 dibasep,cmetp);
	    List_free(&singlehits5);
	  }

	  if ((*nhits3 = List_length(singlehits3)) == 0) {
	    *stage3array3 = (Stage3_T *) NULL;
	  } else {
	    *stage3array3 = (Stage3_T *) List_to_array(singlehits3,NULL);
	    *stage3array3 = Stage3_eval_and_sort(*stage3array3,*nhits3,maxpaths,queryseq3,
						 this3->query_compress_fwd,this3->query_compress_rev,
						 genome_blocks,snp_blocks,genome,quality_string_3,
						 dibasep,cmetp);
	    List_free(&singlehits3);
	  }

	  Stage1_free(&this3,querylength3);
	  Stage1_free(&this5,querylength5);
	  return (Stage3pair_T *) NULL;
	}
      }
    }
  }
}

