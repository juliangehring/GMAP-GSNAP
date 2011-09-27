static char rcsid[] = "$Id: oligoindex.c,v 1.35 2005/02/15 01:52:34 twu Exp $";
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

#define MAXHITS 1000

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* A short oligomer, typically 8 nt, requiring 16 bits or 2 bytes, but
   could be 9 or 10 nt */
typedef UINT4 Shortoligomer_T;

typedef union {
  Genomicpos_T *multiple;
  Genomicpos_T single;
} Positionarray_T;

#define T Oligoindex_T
struct T {
  int *counts;
  Positionarray_T *positions;
};

static int indexsize;
static int oligospace;
static Shortoligomer_T mask;


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

  new->counts = (int *) CALLOC(oligospace,sizeof(int));
  new->positions = (Positionarray_T *) CALLOC(oligospace,sizeof(Positionarray_T));
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



/* Run queryseq through this procedure.  Returns true if there are
   enough bad characters. */
static bool
set_counts (int *counts, char *sequence, int seqlength) {
  int i = 0, noligos = 0;
  int in_counter = 0;
  Shortoligomer_T oligo = 0U;
  char *p, *nt;

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
      noligos++;
      debug(nt = shortoligo_nt(oligo,indexsize);
	    printf("At querypos %d, oligo %s seen\n",i,nt);
	    FREE(nt));
      counts[oligo&mask] = -1;	/* Signifies that the oligo is present */
      in_counter--;
    }
  }

  counts[0] = 0;		/* Ignore poly-A */
  counts[~0U & mask] = 0;	/* Ignore poly-T */

  /*
  printf("Saw %d oligos out of %d expected\n",noligos,seqlength-indexsize+1);
  */

  if (seqlength-indexsize+1 - noligos > 16) {
    return true;
  } else {
    return false;
  }
}

#define ALLOCSIZE 100

/* Run genomicseg through this procedure */
static void
store_positions (Positionarray_T *positions, int *counts, char *sequence, int seqlength) {
  int i = 0, sequencepos = -indexsize;
  int in_counter = 0, count;
  Genomicpos_T oldposition, *ptr;
  Shortoligomer_T oligo = 0U, masked;
  char *p, *nt;

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
      if ((count = counts[masked]) == 0) {
	/* Don't bother, because it's not in the queryseq */
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s wasn't seen in querypos\n",sequencepos,nt,counts[masked]);
	      FREE(nt));
      } else if (count == -1) {
	/* Really means that it's in the queryseq with a genomic count of zero */
	/* We are going to cheat by putting the atomic information
           where the pointer should be, thereby saving us the cost of
           allocating an array. */
	positions[masked].single = (Genomicpos_T) sequencepos;
	counts[masked] = 1;
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s seen, counts is now %d\n",sequencepos,nt,counts[masked]);
	      FREE(nt));
      } else if (count == 1) {
	/* Convert atomic information to an array */
	oldposition = positions[masked].single;
	positions[masked].multiple = (Genomicpos_T *) CALLOC(ALLOCSIZE,sizeof(Genomicpos_T));
	(positions[masked].multiple)[0] = oldposition;
	(positions[masked].multiple)[1] = (Genomicpos_T) sequencepos;
	counts[masked] = 2;
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s seen, counts is now %d\n",sequencepos,nt,counts[masked]);
	      FREE(nt));
      } else if ((count % ALLOCSIZE) == 0) {
	/*
	  positions[masked] = 
	  (Genomicpos_T *) RESIZE(positions[masked],(count+ALLOCSIZE)*sizeof(Genomicpos_T)); 
	*/
	ptr = (Genomicpos_T *) CALLOC(count+ALLOCSIZE,sizeof(Genomicpos_T));
	memcpy((void *) ptr,positions[masked].multiple,count*sizeof(Genomicpos_T));
	FREE(positions[masked].multiple);
	positions[masked].multiple = ptr;

	(positions[masked].multiple)[count] = (Genomicpos_T) sequencepos;
	counts[masked] += 1;
	debug(nt = shortoligo_nt(masked,indexsize);
	      printf("At genomicpos %u, oligo %s seen, counts is now %d\n",sequencepos,nt,counts[masked]);
	      FREE(nt));
      } else {
	(positions[masked].multiple)[count] = (Genomicpos_T) sequencepos;
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
	  if (counts[i] == 1) {
	    printf("Oligo %s => %d entries: %u...%u\n",
		   nt,counts[i],positions[i].single,positions[i].single);
	  } else if (counts[i] > 1) {
	    printf("Oligo %s => %d entries: %u...%u\n",
		   nt,counts[i],(positions[i].multiple)[0],(positions[i].multiple)[counts[i]-1]);
	  }
	  FREE(nt);
	}
	);

  return;
}

/* Run queryseq through this procedure */
static void
clear_positions (Positionarray_T *positions, int *counts, char *sequence, int seqlength) {
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
      if (counts[masked] > 1) {
	FREE(positions[masked].multiple);
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
bool
Oligoindex_tally (T this, Sequence_T queryseq, Sequence_T genomicseg) {
  bool badsequencep;

  memset((void *) this->counts,0,oligospace*sizeof(int));
  badsequencep = set_counts(this->counts,Sequence_pointer(queryseq),Sequence_length(queryseq));
  store_positions(this->positions,this->counts,Sequence_pointer(genomicseg),Sequence_length(genomicseg));
  return badsequencep;
}
  
void
Oligoindex_cleanup (T this, Sequence_T queryseq) {
  clear_positions(this->positions,this->counts,Sequence_pointer(queryseq),Sequence_length(queryseq));
  return;
}


void
Oligoindex_free (T *old) {
  if (*old) {
    FREE((*old)->positions);
    FREE((*old)->counts);
    FREE(*old);
  }
  return;
}

static Genomicpos_T *
lookup (int *nhits, T this, Shortoligomer_T oligo) {
  List_T p;
  int i;

  oligo = oligo & mask;

  if ((*nhits = this->counts[oligo]) == 1) {
    debug(printf("Oligo %04X => 1 entry: %u\n",oligo,this->positions[oligo].single));
    /* Take the atomic information out of the pointer */
    return (Genomicpos_T *) (&(this->positions[oligo].single));
  } else if (*nhits > 1) {
    debug(printf("Oligo %04X => %d entries: %u...%u\n",
		 oligo,*nhits,(this->positions[oligo].multiple)[0],(this->positions[oligo].multiple)[*nhits-1]));
    return this->positions[oligo].multiple;
  } else {
    debug(printf("Oligo %04X not found\n",oligo));
    /* Warning: *nhits might be -1 here, but shouldn't affect anything */
    return NULL;
  }
}


/* Note: Be careful on totalpositions, because nhits may be < 0 */
unsigned int **
Oligoindex_get_mappings (int **npositions, int *totalpositions, T this,
			 Sequence_T queryseq) {
  unsigned int **mappings, *positions;
  int nhits, k, i;
  int querylength;

  char *p;
  int in_counter = 0, querypos = -indexsize;
  Shortoligomer_T oligo = 0U;

  querylength = Sequence_length(queryseq);
  mappings = (unsigned int **) CALLOC(querylength,sizeof(unsigned int *));
  *npositions = (int *) CALLOC(querylength,sizeof(int));
  *totalpositions = 0;

  for (i = 0, p = Sequence_pointer(queryseq); i < querylength; i++, p++) {
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
      if (nhits > MAXHITS) {
	(*npositions)[querypos] = 0;
      } else {
	(*npositions)[querypos] = nhits;
	if (nhits > 0) {
	  *totalpositions += nhits;
	}
      }
      in_counter--;
    }
  }

  return mappings;
}

