static char rcsid[] = "$Id: oligoindex.c,v 1.42 2005/05/20 17:40:36 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "oligoindex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy and memset */
#include "bool.h"
#include "mem.h"
#include "types.h"
#include "genomicpos.h"
#include "list.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* For querycounts, used for trimming ends of queryseq */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif



/* A short oligomer, typically 8 nt, requiring 16 bits or 2 bytes, but
   could be 9 or 10 nt */
typedef UINT4 Shortoligomer_T;

#define T Oligoindex_T
struct T {
  bool *overabundant;
  bool *inquery;
  int *counts;
  Genomicpos_T **positions;
  Genomicpos_T *data;
};

static int indexsize;
static int oligospace;
static Shortoligomer_T mask;

#define ALLOCSIZE 200
#define MAXHITS 200		/* Certain code below depends on whether this is greater than ALLOCSIZE */
#define MAXHITS_GT_ALLOCSIZE 0

static int
power (int base, int exponent) {
  int result = 1, i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}

void
Oligoindex_init (int indexsize0) {
  indexsize = indexsize0;
  mask = ~(~0U << 2*indexsize);
  oligospace = power(4,indexsize);
  return;
}

T
Oligoindex_new () {
  T new = (T) MALLOC(sizeof(*new));

  new->overabundant = (bool *) CALLOC(oligospace,sizeof(bool));
  new->inquery = (bool *) CALLOC(oligospace,sizeof(bool));
  new->counts = (int *) CALLOC(oligospace,sizeof(int));
  new->positions = (Genomicpos_T **) CALLOC(oligospace,sizeof(Genomicpos_T *));
  new->data = (Genomicpos_T *) CALLOC(oligospace*ALLOCSIZE,sizeof(Genomicpos_T));
  return new;
}


/*                87654321 */
#define RIGHT_A 0x00000000
#define RIGHT_C 0x00000001
#define RIGHT_G 0x00000002
#define RIGHT_T 0x00000003

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

static char *
shortoligo_nt (Shortoligomer_T oligo, int oligosize) {
  char *nt;
  int i, j;
  Shortoligomer_T lowbits;

  nt = (char *) CALLOC(oligosize+1,sizeof(char));
  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  return nt;
}


#define POLY_A 0x0000
#define POLY_C 0x5555
#define POLY_G 0xAAAA
#define POLY_T 0xFFFF

double
Oligoindex_set_inquery (int *badoligos, int *trim_start, int *trim_end, T this, Sequence_T queryseq) {
  double oligodepth;
  int nunique = 0, querylength, trimlength;
  int i = 0, noligos = 0;
  int in_counter = 0, querypos;
  Shortoligomer_T oligo = 0U, masked;
  char *p, *nt;
  bool prevunique = false;

  memset((void *) this->inquery,0,oligospace*sizeof(bool));
  memset((void *) this->counts,0,oligospace*sizeof(int));

  querylength = Sequence_fulllength(queryseq);
  for (i = 0, p = Sequence_fullpointer(queryseq); i < querylength; i++, p++) {
    in_counter++;

    switch (toupper(*p)) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }

    if (in_counter == indexsize) {
      masked = oligo&mask;
      noligos++;
      debug(nt = shortoligo_nt(oligo,indexsize);
	    printf("At querypos %d, oligo %s seen\n",i,nt);
	    FREE(nt));
      this->counts[masked] += 1;
      if (this->inquery[masked] == false) {
	nunique += 1;
	this->inquery[masked] = true;
      }
      in_counter--;
    }
  }

  /*
  printf("Saw %d oligos out of %d expected\n",noligos,seqlength-indexsize+1);
  */

  *badoligos = querylength-indexsize+1 - noligos;

  /* Determine where to trim */
  *trim_start = -1;
  *trim_end = 0;

  in_counter = 0;
  querypos = -indexsize;
  for (i = 0, p = Sequence_fullpointer(queryseq); i < querylength; i++, p++) {
    in_counter++;
    querypos++;

    switch (toupper(*p)) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }

    debug1(if (querypos < 0) {
	     /* Print nothing */
	   } else if (in_counter == indexsize) {
	     printf("%d\n",this->counts[oligo&mask]);
	   } else {
	     printf("0\n");
	   });

    if (in_counter == indexsize) {
      if (this->counts[oligo&mask] == 1) {
	if (prevunique == true) {
	  if (*trim_start < 0) {
	    *trim_start = querypos; /* trim 5' at first unique 8-mer */
	  }
	  *trim_end = querypos;	/* trim 3' at last unique 8-mer */
	}
	prevunique = true;
      } else {
	prevunique = false;
      }
      in_counter--;
    } else {
      prevunique = false;
    }
  }

  if (*trim_start < 0) {
    *trim_start = 0;
  } else {
    *trim_start -= 1;		/* To compensate for prevunique */
  }
  *trim_end += indexsize;

  if (nunique == 0) {
    return 1000000.0;
  } else {
    trimlength = *trim_end - *trim_start;
    oligodepth = (double) noligos/(double) nunique;
    return oligodepth;
  }
}



/* Run genomicseg through this procedure */
static void
store_positions (Genomicpos_T **positions, Genomicpos_T *data, 
		 bool *overabundant, bool *inquery,
		 int *counts, char *sequence, int seqlength) {
  int i = 0, sequencepos = -indexsize;
  int in_counter = 0, count;
  Genomicpos_T oldposition, *ptr;
  Shortoligomer_T oligo = 0U, masked;
  char *p, *nt;
  int availslot = 0;

  for (i = 0, p = sequence; i < seqlength; i++, p++) {
    in_counter++;
    sequencepos++;

    debug(printf("At genomicpos %u, char is %c\n",sequencepos,*p));

    switch (toupper(*p)) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }

    if (in_counter == indexsize) {
      masked = oligo&mask;
      if (overabundant[masked] == true) {
	/* Don't bother */
      } else if (inquery[masked] == false) {
	/* Don't bother, because it's not in the queryseq */
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s wasn't seen in querypos\n",sequencepos,nt,counts[masked]);
	      FREE(nt));

      } else if ((count = counts[masked]) == 0) {
	/* Initial allocation; use pre-allocated space */
	positions[masked] = &(data[availslot]);
	availslot += ALLOCSIZE;
	positions[masked][0] = (Genomicpos_T) sequencepos;
	counts[masked] = 1;
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s seen, counts is now %d\n",sequencepos,nt,counts[masked]);
	      FREE(nt));

      } else if (count >= MAXHITS) {
	/* Overabundant: Turn this oligo off */
	overabundant[masked] = true;
	counts[masked] = 0;	/* This is okay as long as check for overabundant is first */
	/* This also prevents allocated memory from being freed later */
	/* Free previously allocated new space */
	/* if (count > ALLOCSIZE) { -- Should always be true if MAXHITS > ALLOCSIZE */
#if MAXHITS_GT_ALLOCSIZE
	FREE(positions[masked]);
#endif
	/* } */

#if MAXHITS_GT_ALLOCSIZE
      } else if ((count % ALLOCSIZE) == 0) {
	/* Change to new space */
	/*
	  positions[masked] = 
	  (Genomicpos_T *) RESIZE(positions[masked],(count+ALLOCSIZE)*sizeof(Genomicpos_T)); 
	*/
	ptr = (Genomicpos_T *) CALLOC(count+ALLOCSIZE,sizeof(Genomicpos_T));
	memcpy((void *) ptr,positions[masked],count*sizeof(Genomicpos_T));
	if (count > ALLOCSIZE) {
	  /* Free previously allocated new space */
	  FREE(positions[masked]);
	}

	positions[masked] = ptr;
	positions[masked][count] = (Genomicpos_T) sequencepos;
	counts[masked] += 1;
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s seen, counts is now %d\n",sequencepos,nt,counts[masked]);
	      FREE(nt));
#endif

      } else {
	positions[masked][count] = (Genomicpos_T) sequencepos;
	counts[masked] += 1;
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s seen, counts is now %d\n",sequencepos,nt,counts[masked]);
	      FREE(nt));
      }
      in_counter--;
    }
  }

  debug(
	for (i = 0; i < oligospace; i++) {
	  nt = shortoligo_nt(i,indexsize);
	  if (counts[i] >= 1) {
	    printf("Oligo %s => %d entries: %u...%u\n",
		   nt,counts[i],positions[i][0],positions[i][counts[i]-1]);
	  }
	  FREE(nt);
	}
	);

  return;
}

/* Run queryseq through this procedure */
static void
clear_positions (Genomicpos_T **positions, int *counts, char *sequence, int seqlength) {
  List_T list;
  int i = 0;
  int in_counter = 0;
  Shortoligomer_T oligo = 0U, masked;
  char *p;

  for (i = 0, p = sequence; i < seqlength; i++, p++) {
    in_counter++;

    switch (toupper(*p)) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }

    if (in_counter == indexsize) {
      masked = oligo&mask;
      if (counts[masked] > ALLOCSIZE) {
	FREE(positions[masked]);
	counts[masked] = 0;	/* To avoid freeing list twice */
      }
      in_counter--;
    }
  }

  return;
}

/* First, we count occurrences of each oligo in queryseq.  This allows
 * us to scan genomicseg intelligently, because then we know whether
 * to store positions for that oligo. */

/* Returns true if the sequence is repetitive */
void
Oligoindex_tally (T this, Sequence_T queryseq, Sequence_T genomicseg) {
  int nunique;
  int trimlength;

  memset((void *) this->counts,0,oligospace*sizeof(int));
  memset((void *) this->overabundant,0,oligospace*sizeof(bool));

  /* These values will prevent oligoindex from getting mappings later */
  this->overabundant[POLY_A] = true;
  this->overabundant[POLY_C] = true;
  this->overabundant[POLY_G] = true;
  this->overabundant[POLY_T] = true;

  trimlength = Sequence_trimlength(queryseq);
  store_positions(this->positions,this->data,this->overabundant,this->inquery,
		  this->counts,Sequence_trimpointer(genomicseg),
		  Sequence_trimlength(genomicseg));
  return;
}
  
void
Oligoindex_cleanup (T this, Sequence_T queryseq) {

#if MAXHITS_GT_ALLOCSIZE
  clear_positions(this->positions,this->counts,Sequence_trimpointer(queryseq),Sequence_trimlength(queryseq));
#endif

  return;
}


void
Oligoindex_free (T *old) {
  if (*old) {
    FREE((*old)->data);
    FREE((*old)->positions);
    FREE((*old)->counts);
    FREE((*old)->inquery);
    FREE((*old)->overabundant);
    FREE(*old);
  }
  return;
}

static Genomicpos_T *
lookup (int *nhits, T this, Shortoligomer_T oligo) {
  List_T p;
  int i;

  oligo = oligo & mask;

  if ((*nhits = this->counts[oligo]) >= 1) {
    debug(printf("Oligo %04X => %d entries: %u...%u\n",
		 oligo,*nhits,this->positions[oligo][0],this->positions[oligo][*nhits-1]));
    return this->positions[oligo];
  } else {
    debug(printf("Oligo %04X not found\n",oligo));
    /* Warning: *nhits might be -1 here, but shouldn't affect anything */
    return NULL;
  }
}


/* Retrieves appropriate oligo information for a given position and
   copies it to that position */
/* Note: Be careful on totalpositions, because nhits may be < 0 */
unsigned int **
Oligoindex_get_mappings (int **npositions, int *totalpositions, T this, 
			 Sequence_T queryseq) {
  unsigned int **mappings, *positions;
  int nhits, k, i;
  int querylength, trimlength, trim_start, trim_end;

  char *p;
  int in_counter = 0, querypos, left_querypos, right_querypos;
  Shortoligomer_T oligo = 0U;

  querylength = Sequence_fulllength(queryseq);
  trimlength = Sequence_trimlength(queryseq);
  trim_start = Sequence_trim_start(queryseq);
  trim_end = Sequence_trim_end(queryseq);
  mappings = (unsigned int **) CALLOC(querylength,sizeof(unsigned int *));
  *npositions = (int *) CALLOC(querylength,sizeof(int));
  *totalpositions = 0;

  querypos = trim_start - indexsize;
  for (i = 0, p = Sequence_trimpointer(queryseq); i < trimlength; i++, p++) {
    in_counter++;
    querypos++;
    
    switch (toupper(*p)) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }

    if (in_counter == indexsize) {
      mappings[querypos] = lookup(&nhits,this,oligo);
      (*npositions)[querypos] = nhits;
      if (nhits > 0) {
	*totalpositions += nhits;
      }
      in_counter--;
    }
  }

  return mappings;
}

