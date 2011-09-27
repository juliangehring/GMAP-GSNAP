static char rcsid[] = "$Id: indexdb_hr.c 45940 2011-08-29 21:09:29Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include "indexdb_hr.h"
#include "indexdbdef.h"
#include "genome_hr.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
/* #define CONVERT_TO_LITTLEENDIAN 1 */
/* Because only call is to generate intersection_diagonals, which are assumed to be in native form */
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "mem.h"
#include "listdef.h"


/* ALLOW_DUPLICATES is possible only if we permit alternative strains */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* results of merge_batches */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* merge_batches */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* binary_search, identify_doubles */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* heapify */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Compoundpos_intersect */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* Compoundpos_find */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif

/* straddling at beginning of genome.  May want to turn on DEBUG11 in stage1hr.c */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif


#define T Indexdb_T

static int index1part;		/* For debugging */
static unsigned int kmer_mask;	/* Was LOW12MER     0x00FFFFFF */
#define right_subst  0x00000001
static unsigned int left_subst; /* Was LEFT_SUBST   0x00100000 */
static unsigned int top_subst;  /* Was TOP_SUBST    0x00400000 */

void
Indexdb_hr_setup (int index1part_in) {
  index1part = index1part_in;
  kmer_mask = ~(~0U << 2*index1part);

  top_subst = (1 << 2*(index1part-1));
  left_subst = (1 << 2*(index1part-2));

  return;
}



typedef struct Batch_T *Batch_T;
struct Batch_T {
  int nentries;
  Genomicpos_T position;
  Genomicpos_T *positionptr;
};

typedef struct Header_T *Header_T;
struct Header_T {
  int heapsize;
  int delta;
};


/* We want to handle 16 nodes.  In the typical heap structure, we need
a node 8 with one left child, node 16, and then sentinels to handle
nodes 17-31.  Instead, we use a different heap structure, where node 1
has only one child, node 2.  Then, each parent node i has a right node
(2*i) and a left node (2*i-1). */

#define PARENT2(i) ((i+1) >> 1)
#define LEFTSIBLING2(i) (i-1)

static const unsigned int PARENT_EVEN[17] =
/* 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 */
  {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8};

static void
heap_insert_even (Batch_T *heap, int *heapsize, Batch_T elt, Genomicpos_T key) {
  unsigned int i, parenti;

  i = ++(*heapsize);
  while (i > 1 && heap[parenti = PARENT_EVEN[i]]->position > key) {
    heap[i] = heap[parenti];
    i = parenti;
  }
  heap[i] = elt;

  return;
}


#define NCASES 48


#if 0
static void
check_heap_even (Batch_T *heap, int heapsize) {
  int i, j;

  for (i = 1; i <= heapsize; i++) {
    if (heap[i]->position > heap[2*i-1]->position) {
      fprintf(stderr,"Failed because position %u at heap %d is > position %u at heap %d\n",
	      heap[i]->position,i,heap[2*i-1]->position,2*i-1);
      for (j = 1; j <= heapsize*2; j++) {
	fprintf(stderr,"%02d %u\n",j,heap[j]->position);
      }
      abort();
    }
    if (heap[i]->position > heap[2*i]->position) {
      fprintf(stderr,"Failed because position %u at heap %d is > position %u at heap %d\n",
	      heap[i]->position,i,heap[2*i]->position,2*i);
      for (j = 1; j <= heapsize*2; j++) {
	fprintf(stderr,"%02d %u\n",j,heap[j]->position);
      }
      abort();
    }
  }
}
#endif

#define READ_THEN_WRITE 1

static Genomicpos_T *
merge_batches_one_heap_16_existing (int *nmerged, struct Batch_T *batchpool, int nentries, int diagterm) {
  Genomicpos_T *positions, *ptr, position, last_position, this_position;
  struct Batch_T sentinel_struct;
  Batch_T batch, sentinel, heap[17];
  int heapsize;
  unsigned int i;
#ifdef READ_THEN_WRITE
  unsigned int smallesti_1, smallesti_2, smallesti;
#else
  unsigned int parenti, smallesti;
#endif

  debug3(printf("starting merge_batches_one_heap_16_existing\n"));

  debug0(int nentries_save = nentries);

  ptr = positions = (Genomicpos_T *) CALLOC(nentries,sizeof(Genomicpos_T));

  /* Set up heap */
  heapsize = 0;
  for (i = 0; i < 16; i++) {
    batch = &(batchpool[i]);
    if (batch->nentries > 0) {
#ifdef WORDS_BIGENDIAN
      batch->position = Bigendian_convert_uint(*batch->positionptr++);
#else
      batch->position = *batch->positionptr++;
#endif
      heap_insert_even(heap,&heapsize,batch,batch->position);
    }
  }

  sentinel_struct.position = (Genomicpos_T) -1; /* infinity */
  sentinel_struct.positionptr = &(sentinel_struct.position);
  sentinel = &sentinel_struct;

  for (i = heapsize+1; i <= 16; i++) {
    heap[i] = sentinel;
  }

  last_position = 0U;
  while (--nentries >= 1) {
    debug3(printf("nentries = %d, top of heap is %u (%d)\n",
		  nentries+1,heap[1]->position,heapsize));

    /* Get minimum */
    batch = heap[1];
#ifdef CONVERT_TO_LITTLEENDIAN
    this_position = Bigendian_convert_uint(batch->position) + diagterm;
#else
    this_position = batch->position + diagterm;
#endif
    if (this_position != last_position) {
      *ptr++ = this_position;
    }
    last_position = this_position;

    if (--batch->nentries <= 0) {
      /* Use last batch (or sentinel) in heap for insertion */
      heap[1] = batch = (heapsize == 1) ? sentinel : heap[heapsize];
      heap[heapsize--] = sentinel;

    } else {
      /* Advance heap, and use this batch for insertion */
#ifdef WORDS_BIGENDIAN
      batch->position = Bigendian_convert_uint(*batch->positionptr++);
#else
      batch->position = *batch->positionptr++;
#endif
    }

    position = batch->position;
    debug3(printf("starting heapify with %u\n",position));

#ifdef READ_THEN_WRITE
    /* Comparison 0/3 */
    debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
    if (position <= heap[2]->position) {
      debug3(printf("Inserting at 1\n"));
      /* heap[1] = batch; -- not necessary because batch is already at heap[1] */
    } else {
      /* Comparison 1/3 */
      debug3(printf("Comparing left %d/right %d: %u and %u\n",
		    3,4,heap[3]->position,heap[4]->position));
      smallesti = 4 - (heap[3]->position < heap[4]->position);
      if (position <= heap[smallesti]->position) {
	debug3(printf("Inserting at 2\n"));
	heap[1] = heap[2];
	heap[2] = batch;
      } else {
	smallesti_1 = smallesti;
	smallesti <<= 1;
	/* Comparison 2/3 */
	debug3(printf("Comparing left %d/right %d: %u and %u\n",
		      smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	if (position <= heap[smallesti]->position) {
	  debug3(printf("Inserting at %d\n",smallesti_1));
	  heap[1] = heap[2];
	  heap[2] = heap[smallesti_1];
	  heap[smallesti_1] = batch;
	} else {
	  smallesti_2 = smallesti;
	  smallesti <<= 1;
	  /* Comparison 3/3 */
	  debug3(printf("Comparing left %d/right %d: %u and %u\n",
			smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	  smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	  if (position <= heap[smallesti]->position) {
	    debug3(printf("Inserting at %d\n",smallesti_2));
	    heap[1] = heap[2];
	    heap[2] = heap[smallesti_1];
	    heap[smallesti_1] = heap[smallesti_2];
	    heap[smallesti_2] = batch;
	  } else {
	    debug3(printf("Inserting at %d\n",smallesti));
	    heap[1] = heap[2];
	    heap[2] = heap[smallesti_1];
	    heap[smallesti_1] = heap[smallesti_2];
	    heap[smallesti_2] = heap[smallesti];
	    heap[smallesti] = batch;
	  }
	}
      }
    }
#else
    /* Comparison 0/3 */
    debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
    if (position <= heap[2]->position) {
      debug3(printf("Inserting at 1\n"));
      /* heap[1] = batch; -- not necessary because batch is already at heap[1] */
    } else {
      heap[1] = heap[2];
      /* Comparison 1/3 */
      debug3(printf("Comparing left %d/right %d: %u and %u\n",
		    3,4,heap[3]->position,heap[4]->position));
      smallesti = 4 - (heap[3]->position < heap[4]->position);
      if (position <= heap[smallesti]->position) {
	debug3(printf("Inserting at 2\n"));
	heap[2] = batch;
      } else {
	heap[2] = heap[smallesti];
	parenti = smallesti;
	smallesti <<= 1;
	/* Comparison 2/3 */
	debug3(printf("Comparing left %d/right %d: %u and %u\n",
		      smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	if (position <= heap[smallesti]->position) {
	  debug3(printf("Inserting at %d\n",parenti));
	  heap[parenti] = batch;
	} else {
	  heap[parenti] = heap[smallesti];
	  parenti = smallesti;
	  smallesti <<= 1;
	  /* Comparison 3/3 */
	  debug3(printf("Comparing left %d/right %d: %u and %u\n",
			smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	  smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	  if (position <= heap[smallesti]->position) {
	    debug3(printf("Inserting at %d\n",parenti));
	    heap[parenti] = batch;
	  } else {
	    heap[parenti] = heap[smallesti];
	    debug3(printf("Inserting at %d\n",smallesti));
	    heap[smallesti] = batch;
	  }
	}
      }
    }
#endif
  }

#ifdef CONVERT_TO_LITTLEENDIAN
  this_position = Bigendian_convert_uint(heap[1]->position) + diagterm;
#else
  this_position = heap[1]->position + diagterm;
#endif
  if (this_position != last_position) {
    *ptr++ = this_position;
  }

  *nmerged = (ptr - positions);

#if 0
  position = positions[0];
  for (i = 1; i < nentries_save; i++) {
    if (positions[i] <= position) {
      abort();
    }
    position = positions[i];
  }
#endif

  debug0(
	 for (i = 0; i < nentries_save; i++) {
	   printf("%u\n",positions[i]);
	 }
	 printf("\n");
	 )

  return positions;
}


static Genomicpos_T *
merge_batches_one_heap_4_existing (int *nmerged, struct Batch_T *batchpool, int nentries, int diagterm) {
  Genomicpos_T *positions, *ptr, position, last_position, this_position;
  struct Batch_T sentinel_struct;
  Batch_T batch, sentinel, heap[5];
  int heapsize;
  unsigned int i;
#ifdef READ_THEN_WRITE
  unsigned int smallesti;
#else
  unsigned int parenti, smallesti;
#endif

  debug3(printf("starting merge_batches_one_heap_4_existing\n"));

  debug0(int nentries_save = nentries);

  ptr = positions = (Genomicpos_T *) CALLOC(nentries,sizeof(Genomicpos_T));

  /* Set up heap */
  heapsize = 0;
  for (i = 0; i < 4; i++) {
    batch = &(batchpool[i]);
    if (batch->nentries > 0) {
#ifdef WORDS_BIGENDIAN
      batch->position = Bigendian_convert_uint(*batch->positionptr++);
#else
      batch->position = *batch->positionptr++;
#endif
      heap_insert_even(heap,&heapsize,batch,batch->position);
    }
  }

  sentinel_struct.position = (Genomicpos_T) -1; /* infinity */
  sentinel_struct.positionptr = &(sentinel_struct.position);
  sentinel = &sentinel_struct;

  for (i = heapsize+1; i <= 4; i++) {
    heap[i] = sentinel;
  }

  last_position = 0U;
  while (--nentries >= 1) {
    debug3(printf("nentries = %d, top of heap is %u (%d)\n",
		  nentries+1,heap[1]->position,heapsize));

    /* Get minimum */
    batch = heap[1];
#ifdef CONVERT_TO_LITTLEENDIAN
    this_position = Bigendian_convert_uint(batch->position) + diagterm;
#else
    this_position = batch->position + diagterm;
#endif
    if (this_position != last_position) {
      *ptr++ = this_position;
    }
    last_position = this_position;


    if (--batch->nentries <= 0) {
      /* Use last batch (or sentinel) in heap for insertion */
      heap[1] = batch = (heapsize == 1) ? sentinel : heap[heapsize];
      heap[heapsize--] = sentinel;

    } else {
      /* Advance heap, and use this batch for insertion */
#ifdef WORDS_BIGENDIAN
      batch->position = Bigendian_convert_uint(*batch->positionptr++);
#else
      batch->position = *batch->positionptr++;
#endif
    }

    position = batch->position;
    debug3(printf("starting heapify with %u\n",position));

#ifdef READ_THEN_WRITE
    /* Comparison 0/3 */
    debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
    if (position <= heap[2]->position) {
      debug3(printf("Inserting at 1\n"));
      /* heap[1] = batch; -- not necessary because batch is already at heap[1] */
    } else {
      /* Comparison 1/3 */
      debug3(printf("Comparing left %d/right %d: %u and %u\n",
		    3,4,heap[3]->position,heap[4]->position));
      smallesti = 4 - (heap[3]->position < heap[4]->position);
      if (position <= heap[smallesti]->position) {
	debug3(printf("Inserting at 2\n"));
	heap[1] = heap[2];
	heap[2] = batch;
      } else {
	debug3(printf("Inserting at %d\n",smallesti));
	heap[1] = heap[2];
	heap[2] = heap[smallesti];
	heap[smallesti] = batch;
      }
    }

#else
    /* Comparison 0/3 */
    debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
    if (position <= heap[2]->position) {
      debug3(printf("Inserting at 1\n"));
      /* heap[1] = batch; -- not necessary because batch is already at heap[1] */
    } else {
      heap[1] = heap[2];
      /* Comparison 1/3 */
      debug3(printf("Comparing left %d/right %d: %u and %u\n",
		    3,4,heap[3]->position,heap[4]->position));
      smallesti = 4 - (heap[3]->position < heap[4]->position);
      if (position <= heap[smallesti]->position) {
	debug3(printf("Inserting at 2\n"));
	heap[2] = batch;
      } else {
	heap[2] = heap[smallesti];
	heap[smallesti] = batch;
      }
    }

#endif
  }

#ifdef CONVERT_TO_LITTLEENDIAN
  this_position = Bigendian_convert_uint(heap[1]->position) + diagterm;
#else
  this_position = heap[1]->position + diagterm;
#endif
  if (this_position != last_position) {
    *ptr++ = this_position;
  }

  *nmerged = (ptr - positions);

#if 0
  position = positions[0];
  for (i = 1; i < nentries_save; i++) {
    if (positions[i] <= position) {
      abort();
    }
    position = positions[i];
  }
#endif

  debug0(
	 for (i = 0; i < nentries_save; i++) {
	   printf("%u\n",positions[i]);
	 }
	 printf("\n");
	 )


  return positions;
}


/************************************************************************
 *  The following positions functions are take from indexdb.c
 ************************************************************************/

static void
positions_move_absolute (int positions_fd, Positionsptr_T ptr) {
  off_t offset = ptr*((off_t) sizeof(Genomicpos_T));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %u*%lu=%lu\n",
	    ptr,sizeof(Genomicpos_T),(long unsigned int) offset);
    perror("Error in indexdb.c, positions_move_absolute");
    exit(9);
  }
  return;
}

static void
positions_read_multiple (int positions_fd, Genomicpos_T *values, int n) {
  int i;
  Genomicpos_T value;
  unsigned char buffer[4];

#ifdef WORDS_BIGENDIAN
  /* Need to keep in bigendian format */
  for (i = 0; i < n; i++) {
    read(positions_fd,buffer,4);

    value = (buffer[0] & 0xff);
    value <<= 8;
    value |= (buffer[1] & 0xff);
    value <<= 8;
    value |= (buffer[2] & 0xff);
    value <<= 8;
    value |= (buffer[3] & 0xff);

    values[i] = value;
  }
#else
  for (i = 0; i < n; i++) {
    read(positions_fd,buffer,4);

    value = (buffer[3] & 0xff);
    value <<= 8;
    value |= (buffer[2] & 0xff);
    value <<= 8;
    value |= (buffer[1] & 0xff);
    value <<= 8;
    value |= (buffer[0] & 0xff);

    values[i] = value;
  }
#endif

  return;
}


static Genomicpos_T *
point_one_shift (int *nentries, T this, Storedoligomer_T subst) {
  Genomicpos_T *positions;
  Positionsptr_T ptr0, end0;
#ifdef DEBUG
  int i;
#endif

  if (this->offsets) {
#ifdef WORDS_BIGENDIAN
    if (this->offsets_access == ALLOCATED) {
      ptr0 = this->offsets[subst];
      end0 = this->offsets[subst+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsets[subst]);
      end0 = Bigendian_convert_uint(this->offsets[subst+1]);
    }
#else
    ptr0 = this->offsets[subst];
    end0 = this->offsets[subst+1];
#endif

  } else {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,subst);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,subst);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,subst);
#endif
  }

  debug(printf("point_one_shift: %08X %u %u\n",subst,ptr0,end0));

  if ((*nentries = end0 - ptr0) == 0) {
    return (Genomicpos_T *) NULL;
  } else {
    if (this->positions_access == FILEIO) {
      positions = (Genomicpos_T *) CALLOC(*nentries,sizeof(Genomicpos_T));
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
      positions_move_absolute(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else {
      /* ALLOCATED or MMAPPED */
      positions = &(this->positions[ptr0]);
    }
  }
      
#ifdef WORDS_BIGENDIAN
  debug(
	printf("%d entries:",*nentries);
	for (i = 0; i < *nentries; i++) {
	  printf(" %u",Bigendian_convert_uint(positions[i]));
	}
	printf("\n");
	);
#else
  debug(
	printf("%d entries:",*nentries);
	for (i = 0; i < *nentries; i++) {
	  printf(" %u",positions[i]);
	}
	printf("\n");
	);
#endif
  
  return positions;
}





/*                87654321 */
#define RIGHT_A 0x00000000
#define RIGHT_C 0x00000001
#define RIGHT_G 0x00000002
#define RIGHT_T 0x00000003

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

#ifdef DEBUG
static char *
shortoligo_nt (Storedoligomer_T oligo, int oligosize) {
  char *nt;
  int i, j;
  Storedoligomer_T lowbits;

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
#endif

static int
count_one_shift (T this, Storedoligomer_T subst, int nadjacent) {
  Positionsptr_T ptr0, end0;

  if (this->offsets) {
#ifdef WORDS_BIGENDIAN
    if (this->offsets_access == ALLOCATED) {
      ptr0 = this->offsets[subst];
      end0 = this->offsets[subst+nadjacent];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsets[subst]);
      end0 = Bigendian_convert_uint(this->offsets[subst+nadjacent]);
    }
#else
    ptr0 = this->offsets[subst];
    end0 = this->offsets[subst+nadjacent];
#endif

  } else {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_only_from_gammas(this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,subst);
      end0 = Genome_offsetptr_only_from_gammas(this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,subst+nadjacent);
    } else {
      ptr0 = Genome_offsetptr_only_from_gammas_bigendian(this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,subst);
      end0 = Genome_offsetptr_only_from_gammas_bigendian(this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,subst+nadjacent);
    }
#else
    ptr0 = Genome_offsetptr_only_from_gammas(this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,subst);
    end0 = Genome_offsetptr_only_from_gammas(this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,subst+nadjacent);
#endif
  }

  debug(printf("count_one_shift: oligo = %06X (%s), %u - %u = %u\n",
	       subst,shortoligo_nt(subst,index1part),end0,ptr0,end0-ptr0));
  return (end0 - ptr0);

}



/************************************************************************
 *   Counting procedures
 ************************************************************************/

/* Don't mask out leftmost nucleotides with LOWXXMER */
int
Indexdb_count_left_subst_2 (T this, Storedoligomer_T oligo) {
  int nentries = 0;
  Storedoligomer_T base;
  int i;

  debug(printf("count_left_subst_2: oligo = %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

#ifdef ALLOW_DUPLICATES
  /* Right shift */
  base = (oligo >> 4);
  for (i = 0; i < 16; i++, base += left_subst) {
    nentries += count_one_shift(this,base);
  }
#else
  /* Right shift */
  base = (oligo >> 4);
  debug(printf("shift right => %06X (%s)\n",base,shortoligo_nt(base,index1part)));
  for (i = 0; i < 16; i++, base += left_subst) {
    nentries += count_one_shift(this,base,/*nadjacent*/1);
  }
#endif
      
  return nentries;
}


/* Don't mask out leftmost nucleotides with LOWXXMER */
int
Indexdb_count_left_subst_1 (T this, Storedoligomer_T oligo) {
  int nentries = 0;
  Storedoligomer_T base;
  int i;

  debug(printf("count_left_subst_1: oligo = %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

#ifdef ALLOW_DUPLICATES
  /* Zero shift. */
  base = (oligo >> 2);
  for (i = 0; i < 4; i++, base += top_subst) {
    nentries += count_one_shift(this,base);
  }
#else
  /* Zero shift. */
  base = (oligo >> 2);
  for (i = 0; i < 4; i++, base += top_subst) {
    nentries += count_one_shift(this,base,/*nadjacent*/1);
  }
#endif
      
  return nentries;
}


int
Indexdb_count_right_subst_2 (T this, Storedoligomer_T oligo) {
  int nentries;
  Storedoligomer_T base;
#ifdef ALLOW_DUPLICATES
  int i;
#endif
#ifdef DEBUG
  int i;
#endif

  debug(printf("count_right_subst_2: oligo = %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

#ifdef ALLOW_DUPLICATES
  /* Left shift */
  base = (oligo << 4) & kmer_mask;
  nentries = 0;
  for (i = 0; i < 16; i++, base += right_subst) {
    nentries += count_one_shift(this,base);
  }
#else
  /* Left shift */
  base = (oligo << 4) & kmer_mask;
  nentries = count_one_shift(this,base,/*nadjacent*/16);

  debug(
	printf("Details\n");
	nentries = 0;
	for (i = 0; i < 16; i++, base += right_subst) {
	  nentries += count_one_shift(this,base,/*nadjacent*/1);
	}
	);
#endif
      
  return nentries;
}


int
Indexdb_count_right_subst_1 (T this, Storedoligomer_T oligo) {
  int nentries;
  Storedoligomer_T base;
#ifdef ALLOW_DUPLICATES
  int i;
#endif
#ifdef DEBUG
  int i;
#endif

  debug(printf("count_right_subst_1: oligo = %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

#ifdef ALLOW_DUPLICATES
  /* Zero shift */
  base = (oligo << 2) & kmer_mask;
  nentries = 0;
  for (i = 0; i < 4; i++, base += right_subst) {
    nentries += count_one_shift(this,base);
  }
#else
  /* Zero shift */
  base = (oligo << 2) & kmer_mask;
  nentries = count_one_shift(this,base,/*nadjacent*/4);

  debug(
	printf("Details\n");
	nentries = 0;
	for (i = 0; i < 4; i++, base += right_subst) {
	  nentries += count_one_shift(this,base,/*nadjacent*/1);
	}
	);
#endif
      
  return nentries;
}


/************************************************************************/


static bool free_positions_p;	/* Needs to be true if Indexdb positions are FILEIO */

void
Compoundpos_init_positions_free (bool positions_fileio_p) {
  if (positions_fileio_p == true) {
    free_positions_p = true;
  } else {
    free_positions_p = false;
  }
  return;
}



struct Compoundpos_T {
  int n;

  Genomicpos_T *positions[16];
  int npositions[16];

  struct Batch_T batchpool[16];
  Batch_T heap[17];
  int heapsize;
  struct Batch_T sentinel_struct;
  Batch_T sentinel;

  Genomicpos_T *positions_reset[16]; /* altered by find_nomiss_aux and find_onemiss_aux */
  int npositions_reset[16]; /* altered by find_nomiss_aux and find_onemiss_aux */
};


void
Compoundpos_set (Compoundpos_T compoundpos) {
  int i;

  for (i = 0; i < compoundpos->n; i++) {
    compoundpos->positions_reset[i] = compoundpos->positions[i];
    compoundpos->npositions_reset[i] = compoundpos->npositions[i];
  }
  return;
}

void
Compoundpos_reset (Compoundpos_T compoundpos) {
  int i;

  for (i = 0; i < compoundpos->n; i++) {
    compoundpos->positions[i] = compoundpos->positions_reset[i];
    compoundpos->npositions[i] = compoundpos->npositions_reset[i];
  }
  return;
}


void
Compoundpos_print_sizes (Compoundpos_T compoundpos) {
  int i;

  for (i = 0; i < compoundpos->n; i++) {
    printf(" %d",compoundpos->npositions[i]);
  }

  return;
}


void
Compoundpos_dump (Compoundpos_T compoundpos, int diagterm) {
  int i, j;

  printf("%d diagonals: ",compoundpos->n);
  for (i = 0; i < compoundpos->n; i++) {
    printf(" %d",compoundpos->npositions[i]);
  }
  printf("\n");

  for (i = 0; i < compoundpos->n; i++) {
    for (j = 0; j < compoundpos->npositions[i]; j++) {
#ifdef WORDS_BIGENDIAN
      printf(" compound%d.%d:%u+%d\n",
	     i,j,Bigendian_convert_uint(compoundpos->positions[i][j]),diagterm);
#else
      printf(" compound%d.%d:%u+%d\n",i,j,compoundpos->positions[i][j],diagterm);
#endif
    }
  }
  return;
}


void
Compoundpos_free (Compoundpos_T *old) {
  int i;

  if (*old) {
    if (free_positions_p == true) {
      for (i = 0; i < (*old)->n; i++) {
	FREE((*old)->positions[i]);
      }
    }

    /* No need, since allocated statically.  FREE((*old)->npositions); */
    /* No need, since allocated statically.  FREE((*old)->positions); */
  
    FREE(*old);
  }
  return;
}


Compoundpos_T
Indexdb_compoundpos_left_subst_2 (T this, Storedoligomer_T oligo) {
  Compoundpos_T compoundpos = (Compoundpos_T) MALLOC(sizeof(*compoundpos));
  Storedoligomer_T base;
  int i;

  debug(printf("compoundpos_left_subst_2: %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

  compoundpos->n = 16;
  /* compoundpos->npositions = (int *) CALLOC(16,sizeof(int)); */
  /* compoundpos->positions = (Genomicpos_T **) CALLOC(16,sizeof(Genomicpos_T *)); */

  /* Right shift */
  base = (oligo >> 4);
  for (i = 0; i < 16; i++, base += left_subst) {
    compoundpos->positions[i] = point_one_shift(&(compoundpos->npositions[i]),this,base);
  }

  return compoundpos;
}

Compoundpos_T
Indexdb_compoundpos_left_subst_1 (T this, Storedoligomer_T oligo) {
  Compoundpos_T compoundpos = (Compoundpos_T) MALLOC(sizeof(*compoundpos));
  Storedoligomer_T base;
  int i;

  debug(printf("compoundpos_left_subst_1: %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

  compoundpos->n = 4;
  /* compoundpos->npositions = (int *) CALLOC(4,sizeof(int)); */
  /* compoundpos->positions = (Genomicpos_T **) CALLOC(4,sizeof(Genomicpos_T *)); */

  /* Zero shift */
  base = (oligo >> 2);
  for (i = 0; i < 4; i++, base += top_subst) {
    compoundpos->positions[i] = point_one_shift(&(compoundpos->npositions[i]),this,base);
  }

  return compoundpos;
}

Compoundpos_T
Indexdb_compoundpos_right_subst_2 (T this, Storedoligomer_T oligo) {
  Compoundpos_T compoundpos = (Compoundpos_T) MALLOC(sizeof(*compoundpos));
  Storedoligomer_T base;
  int i;

  debug(printf("compoundpos_right_subst_2: %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

  compoundpos->n = 16;
  /* compoundpos->npositions = (int *) CALLOC(16,sizeof(int)); */
  /* compoundpos->positions = (Genomicpos_T **) CALLOC(16,sizeof(Genomicpos_T *)); */

  /* Left shift */
  base = (oligo << 4) & kmer_mask;
  for (i = 0; i < 16; i++, base += right_subst) {
    compoundpos->positions[i] = point_one_shift(&(compoundpos->npositions[i]),this,base);
  }

  return compoundpos;
}

Compoundpos_T
Indexdb_compoundpos_right_subst_1 (T this, Storedoligomer_T oligo) {
  Compoundpos_T compoundpos = (Compoundpos_T) MALLOC(sizeof(*compoundpos));
  Storedoligomer_T base;
  int i;

  debug(printf("compoundpos_right_subst_1: %06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));

  compoundpos->n = 4;
  /* compoundpos->npositions = (int *) CALLOC(4,sizeof(int)); */
  /* compoundpos->positions = (Genomicpos_T **) CALLOC(4,sizeof(Genomicpos_T *)); */

  /* Zero shift */
  base = (oligo << 2) & kmer_mask;
  for (i = 0; i < 4; i++, base += right_subst) {
    compoundpos->positions[i] = point_one_shift(&(compoundpos->npositions[i]),this,base);
  }

  return compoundpos;
}



/************************************************************************/

static int
binary_search (int lowi, int highi, Genomicpos_T *positions, Genomicpos_T goal) {
  bool foundp = false;
  int middlei;

#ifdef NOBINARY
  return lowi;
#endif

  if (goal == 0U) {
    return lowi;
  }

  while (!foundp && lowi < highi) {
#ifdef WORDS_BIGENDIAN
    middlei = (lowi+highi)/2;
    debug2(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		  lowi,Bigendian_convert_uint(positions[lowi]),
		  middlei,Bigendian_convert_uint(positions[middlei]),
		  highi,Bigendian_convert_uint(positions[highi]),goal));
    if (goal < Bigendian_convert_uint(positions[middlei])) {
      highi = middlei;
    } else if (goal > Bigendian_convert_uint(positions[middlei])) {
      lowi = middlei + 1;
    } else {
      foundp = true;
    }
#else
    middlei = (lowi+highi)/2;
    debug2(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		  lowi,positions[lowi],middlei,positions[middlei],
		  highi,positions[highi],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      foundp = true;
    }
#endif
  }

  if (foundp == true) {
    return middlei;
  } else {
    return highi;
  }
}


void
Compoundpos_heap_init (Compoundpos_T compoundpos, int querylength, int diagterm) {
  Batch_T batch;
  int startbound, i;

  compoundpos->heapsize = 0;
  for (i = 0; i < compoundpos->n; i++) {
    batch = &(compoundpos->batchpool[i]);
    batch->positionptr = compoundpos->positions[i];
    batch->nentries = compoundpos->npositions[i];
    if (diagterm < querylength) {
      startbound = querylength - diagterm;
#ifdef WORDS_BIGENDIAN
      while (batch->nentries > 0 && Bigendian_convert_uint(*batch->positionptr) < (unsigned int) startbound) {
	debug11(printf("Eliminating diagonal %u as straddling beginning of genome (Compoundpos_heap_init)\n",
		       Bigendian_convert_uint(*batch->positionptr)));
	++batch->positionptr;
	--batch->nentries;
      }
#else
      while (batch->nentries > 0 && *batch->positionptr < (unsigned int) startbound) {
	debug11(printf("Eliminating diagonal %u as straddling beginning of genome (Compoundpos_heap_init)\n",
		       *batch->positionptr));
	++batch->positionptr;
	--batch->nentries;
      }
#endif
    }
    if (batch->nentries > 0) {
#ifdef WORDS_BIGENDIAN
      batch->position = Bigendian_convert_uint(*batch->positionptr);
#else
      batch->position = *batch->positionptr;
#endif
      heap_insert_even(compoundpos->heap,&compoundpos->heapsize,batch,batch->position);
    }
  }

  compoundpos->sentinel_struct.position = (Genomicpos_T) -1U; /* infinity */
  compoundpos->sentinel_struct.positionptr = &(compoundpos->sentinel_struct.position);
  compoundpos->sentinel = &compoundpos->sentinel_struct;

  for (i = compoundpos->heapsize+1; i <= compoundpos->n; i++) {
    compoundpos->heap[i] = compoundpos->sentinel;
  }

  return;
}


/* Used by DEBUG3 and DEBUG6 */
static void
heap_even_dump (Batch_T *heap, int heapsize) {
  int i;
  Batch_T batch;

  for (i = 1; i <= heapsize; i++) {
    batch = heap[i];
    printf("#%d--%d:%u  ",i,batch->nentries,batch->position);
  }
  printf("\n");
}



/* Returns true if found.  emptyp is true only if every batch is
   empty.  If procedure returns true, empty is guaranteed to be
   false. */
bool
Compoundpos_find (bool *emptyp, Compoundpos_T compoundpos, Genomicpos_T local_goal) {
  Batch_T *heap = compoundpos->heap, batch;
  int i, j;

  debug6(printf("\nEntering Compoundpos_find with local_goal %u\n",local_goal));

  *emptyp = true;
  i = 1;
  while (i <= compoundpos->heapsize) {
    debug6(printf("Compoundpos_find iteration, heapsize %d:\n",compoundpos->heapsize));
    debug6(heap_even_dump(heap,compoundpos->heapsize));

    batch = heap[i];
#ifdef WORDS_BIGENDIAN
    if (batch->nentries > 0 && Bigendian_convert_uint(*batch->positionptr) < local_goal) {
      j = 1;
      while (j < batch->nentries && Bigendian_convert_uint(batch->positionptr[j]) < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= batch->nentries) {
	j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
      } else {
	j = binary_search(j >> 1,j,batch->positionptr,local_goal);
      }
      batch->positionptr += j;
      batch->nentries -= j;
      debug6(printf("binary search jump %d positions to %d:%u\n",
		    j,batch->nentries,Bigendian_convert_uint(*batch->positionptr)));
    }
#else
    if (batch->nentries > 0 && *batch->positionptr < local_goal) {
      j = 1;
      while (j < batch->nentries && batch->positionptr[j] < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= batch->nentries) {
	j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
      } else {
	j = binary_search(j >> 1,j,batch->positionptr,local_goal);
      }
      batch->positionptr += j;
      batch->nentries -= j;
      debug6(printf("binary search jump %d positions to %d:%u\n",
		    j,batch->nentries,*batch->positionptr));
    }
#endif

    if (batch->nentries <= 0) {
      /* Empty, so continue with loop */
      /* Move last heap to this one, and reduce heapsize */
      compoundpos->heap[i] = compoundpos->heap[compoundpos->heapsize];
      --compoundpos->heapsize;

#ifdef WORDS_BIGENDIAN
    } else if (Bigendian_convert_uint(*batch->positionptr) > local_goal) {
      /* Already advanced past goal, so continue with loop */
      debug6(printf("Setting emptyp to be false\n"));
      *emptyp = false;
      i++;
#else
    } else if (*batch->positionptr > local_goal) {
      /* Already advanced past goal, so continue with loop */
      debug6(printf("Setting emptyp to be false\n"));
      *emptyp = false;
      i++;
#endif
    } else {
      /* Found goal, so return */
      debug6(printf("Setting emptyp to be false\n"));
      *emptyp = false;
#ifdef WORDS_BIGENDIAN
      debug6(printf("Found! Returning position %u\n",Bigendian_convert_uint(*batch->positionptr)));
#else
      debug6(printf("Found! Returning position %u\n",*batch->positionptr));
#endif
      ++batch->positionptr;
      --batch->nentries;
      return true;
    }
  }

  /* Done with loop: Fail. */
  debug6(printf("Returning emptyp %d\n",*emptyp));
  return false;
}



/* Returns 0 if heapsize is 0, else 1, and returns smallest value >= local_goal */
int
Compoundpos_search (Genomicpos_T *value, Compoundpos_T compoundpos, Genomicpos_T local_goal) {
  int parenti, smallesti, j;
  Batch_T batch, *heap = compoundpos->heap;
  Genomicpos_T position;

  debug3(printf("\nEntering Compoundpos_search with local_goal %u\n",local_goal));
  if (compoundpos->heapsize <= 0) {
    debug3(printf("Returning because heapsize is %d\n",compoundpos->heapsize));
    return 0;
  }

  if (compoundpos->n == 4) {
    while (compoundpos->heapsize > 0 && (batch = heap[1])->position < local_goal) {
      debug3(printf("Compoundpos_search iteration, heapsize %d:\n",compoundpos->heapsize));
      debug3(heap_even_dump(heap,compoundpos->heapsize));
#ifdef WORDS_BIGENDIAN
      if (batch->nentries > 0 && Bigendian_convert_uint(*batch->positionptr) < local_goal) {
	j = 1;
	while (j < batch->nentries && Bigendian_convert_uint(batch->positionptr[j]) < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= batch->nentries) {
	  j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
	} else {
	  j = binary_search(j >> 1,j,batch->positionptr,local_goal);
	}
	batch->positionptr += j;
	batch->nentries -= j;
	debug3(printf("binary search jump %d positions to %d:%u\n",
		      j,batch->nentries,Bigendian_convert_uint(*batch->positionptr)));
      }
      batch->position = Bigendian_convert_uint(*batch->positionptr);
#else
      if (batch->nentries > 0 && *batch->positionptr < local_goal) {
	j = 1;
	while (j < batch->nentries && batch->positionptr[j] < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= batch->nentries) {
	  j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
	} else {
	  j = binary_search(j >> 1,j,batch->positionptr,local_goal);
	}
	batch->positionptr += j;
	batch->nentries -= j;
	debug3(printf("binary search jump %d positions to %d:%u\n",
		      j,batch->nentries,*batch->positionptr));
      }
      batch->position = *batch->positionptr;
#endif

      if (batch->nentries <= 0) {
	debug3(printf("top of heap found to be empty\n"));
	heap[1] = batch = (compoundpos->heapsize == 1) ? 
	  compoundpos->sentinel : heap[compoundpos->heapsize];
	heap[compoundpos->heapsize--] = compoundpos->sentinel;
      }
      
      position = batch->position;
      debug3(printf("heapify downward on %u\n",position));
      debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
      if (position <= heap[2]->position) {
	debug3(printf("Inserting at 1\n"));
	/* heap[1] = batch; -- not necessary because batch is already at heap[1] */
      } else {
	heap[1] = heap[2];
	debug3(printf("Comparing left %d/right %d: %u and %u\n",
		      3,4,heap[3]->position,heap[4]->position));
	smallesti = 4 - (heap[3]->position < heap[4]->position);
	if (position <= heap[smallesti]->position) {
	  debug3(printf("Inserting at 2\n"));
	  heap[2] = batch;
	} else {
	  debug3(printf("Inserting at %d\n",smallesti));
	  heap[2] = heap[smallesti];
	  heap[smallesti] = batch;
	}
      }
    }
    if (batch->position == local_goal) {
      *value = batch->position;
      debug3(printf("Found! Returning position %u\n",*value));
      return 1;
    }

  } else {
    /* 16 batches */
    while (compoundpos->heapsize > 0 && (batch = heap[1])->position < local_goal) {
      debug3(printf("Compoundpos_search iteration, heapsize %d:\n",compoundpos->heapsize));
      debug3(heap_even_dump(heap,compoundpos->heapsize));
#ifdef WORDS_BIGENDIAN
      if (batch->nentries > 0 && Bigendian_convert_uint(*batch->positionptr) < local_goal) {
	j = 1;
	while (j < batch->nentries && Bigendian_convert_uint(batch->positionptr[j]) < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= batch->nentries) {
	  j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
	} else {
	  j = binary_search(j >> 1,j,batch->positionptr,local_goal);
	}
	batch->positionptr += j;
	batch->nentries -= j;
	debug3(printf("binary search jump %d positions to %d:%u\n",
		      j,batch->nentries,Bigendian_convert_uint(*batch->positionptr)));
      }
      batch->position = Bigendian_convert_uint(*batch->positionptr);
#else
      if (batch->nentries > 0 && *batch->positionptr < local_goal) {
	j = 1;
	while (j < batch->nentries && batch->positionptr[j] < local_goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= batch->nentries) {
	  j = binary_search(j >> 1,batch->nentries,batch->positionptr,local_goal);
	} else {
	  j = binary_search(j >> 1,j,batch->positionptr,local_goal);
	}
	batch->positionptr += j;
	batch->nentries -= j;
	debug3(printf("binary search jump %d positions to %d:%u\n",
		      j,batch->nentries,*batch->positionptr));
      }
      batch->position = *batch->positionptr;
#endif

      if (batch->nentries <= 0) {
	debug3(printf("top of heap found to be empty\n"));
	heap[1] = batch = (compoundpos->heapsize == 1) ? 
	  compoundpos->sentinel : heap[compoundpos->heapsize];
	heap[compoundpos->heapsize--] = compoundpos->sentinel;
      }
      
      position = batch->position;
      debug3(printf("heapify downward on %u\n",position));
      /* Comparison 0/3 */
      debug3(printf("Comparing right %d: %u\n",2,heap[2]->position));
      if (position <= heap[2]->position) {
	debug3(printf("Inserting at 1\n"));
	/* heap[1] = batch; -- not necessary because batch is already at heap[1] */
      } else {
	heap[1] = heap[2];
	/* Comparison 1/3 */
	debug3(printf("Comparing left %d/right %d: %u and %u\n",
		      3,4,heap[3]->position,heap[4]->position));
	smallesti = 4 - (heap[3]->position < heap[4]->position);
	if (position <= heap[smallesti]->position) {
	  debug3(printf("Inserting at 2\n"));
	  heap[2] = batch;
	} else {
	  heap[2] = heap[smallesti];
	  parenti = smallesti;
	  smallesti <<= 1;
	  /* Comparison 2/3 */
	  debug3(printf("Comparing left %d/right %d: %u and %u\n",
			smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	  smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	  if (position <= heap[smallesti]->position) {
	    debug3(printf("Inserting at %d\n",parenti));
	    heap[parenti] = batch;
	  } else {
	    heap[parenti] = heap[smallesti];
	    parenti = smallesti;
	    smallesti <<= 1;
	    /* Comparison 3/3 */
	    debug3(printf("Comparing left %d/right %d: %u and %u\n",
			  smallesti-1,smallesti,heap[smallesti-1]->position,heap[smallesti]->position));
	    smallesti -= (heap[LEFTSIBLING2(smallesti)]->position < heap[smallesti]->position);
	    if (position <= heap[smallesti]->position) {
	      debug3(printf("Inserting at %d\n",parenti));
	      heap[parenti] = batch;
	    } else {
	      heap[parenti] = heap[smallesti];
	      debug3(printf("Inserting at %d\n",smallesti));
	      heap[smallesti] = batch;
	    }
	  }
	}
      }
    }
    if (batch->position == local_goal) {
      *value = batch->position;
      debug3(printf("Found! Returning position %u\n",*value));
      return 1;
    }
  }

  *value = batch->position;
  debug3(printf("Returning position %u\n",*value));
  return 1;
}


Genomicpos_T *
Indexdb_merge_compoundpos (int *nmerged, Compoundpos_T compoundpos, int diagterm) {
  int i;
  Batch_T batch;
  struct Batch_T batchpool[16];
  int nentries = 0;

  debug(printf("merge_compoundpos, sizes:"));

  for (i = 0; i < compoundpos->n; i++) {
    batch = &(batchpool[i]);
    batch->positionptr = compoundpos->positions[i];
    batch->nentries = compoundpos->npositions[i];
    debug(printf(" %d",batch->nentries));
    nentries += batch->nentries;
  }
  debug(printf("\n"));

  if (nentries == 0) {
    *nmerged = 0;
    return (Genomicpos_T *) NULL;
  } else if (compoundpos->n == 4) {
    return merge_batches_one_heap_4_existing(&(*nmerged),batchpool,nentries,diagterm);
  } else {
    return merge_batches_one_heap_16_existing(&(*nmerged),batchpool,nentries,diagterm);
  }
}



int
Indexdb_count_no_subst (T this, Storedoligomer_T oligo) {
  if (this->index1interval == 1) {
    abort();
  } else if (this->index1interval == 3) {
#ifdef ALLOW_DUPLICATES
    return count_one_shift(this,oligo);
#else
    return count_one_shift(this,oligo,/*nadjacent*/1);
#endif
  } else {
    abort();
  }
}

#if 0
int
Indexdb_gsnapbase (T this) {
  if (this->index1interval == 1) {
    return index1part;
  } else if (this->index1interval == 3) {
    return index1part - 2;
  } else {
    abort();
  }
}
#endif
