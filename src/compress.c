static char rcsid[] = "$Id: compress.c 110675 2013-10-10 02:33:10Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "compress.h"
#include "compress-write.h"


#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>		/* For isalpha, toupper */
#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
#include "complement.h"
#include "assert.h"
#include "mem.h"		/* For Compress_new */

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif


#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* SIMD */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif




#define T Compress_T
struct T {
  Genomecomp_T *blocks;
  int nblocks;
  Genomecomp_T **shift_array;
  bool availp[32];
};


void
Compress_free (T *old) {
  if (*old) {
    FREE((*old)->shift_array[0]);
    FREE((*old)->shift_array);
#ifdef HAVE_SSE2
    _mm_free((*old)->blocks);
#else
    FREE((*old)->blocks);
#endif
    FREE(*old);
  }
  return;
}
void
Compress_print (T this) {
  int ptr = 0;

  while (ptr < this->nblocks*COMPRESS_BLOCKSIZE) {
    printf("high: %08X  low: %08X  flags: %08X\n",
	   this->blocks[ptr],this->blocks[ptr+1],this->blocks[ptr+2]);
    ptr += COMPRESS_BLOCKSIZE;
  }
  printf("\n");
  return;
}


int
Compress_nblocks (T this) {
  return this->nblocks;
}


static void
write_chars (Genomecomp_T high, Genomecomp_T low, Genomecomp_T flags) {
  char Buffer[33];
  int i;

  Buffer[32] = '\0';
  /* printf("%08X %08X %08X => ",high,low,flags); */

  for (i = 0; i < 32; i++) {
    switch (((high & 0x01) << 1) | (low & 0x01)) {
    case 0U: Buffer[i] = 'A'; break;
    case 1U: Buffer[i] = 'C'; break;
    case 2U: Buffer[i] = 'G'; break;
    case 3U: Buffer[i] = 'T'; break;
    default: abort();
    }
    high >>= 1;
    low >>= 1;
  }

  if (flags != 0U) {
    for (i = 0; i < 32; i++) {
      if (flags & 0x01) {
	Buffer[i] = 'N';
      }
      flags >>= 1;
    }
  }

  printf("%s",Buffer);
  return;
}


void
Compress_print_blocks (Genomecomp_T *blocks, int nblocks) {
  int ptr = 0;

  while (ptr < nblocks*COMPRESS_BLOCKSIZE) {
    printf("high: %08X  low: %08X  flags: %08X\t",
	   blocks[ptr],blocks[ptr+1],blocks[ptr+2]);
    write_chars(blocks[ptr],blocks[ptr+1],blocks[ptr+2]);
    printf("\n");
    ptr += COMPRESS_BLOCKSIZE;
  }
  printf("\n");
  return;
}


/*                   87654321 */
#define LEFT_SET   0x80000000
#define LEFT_CLEAR 0x00000000


T
Compress_new_fwd (char *gbuffer, Chrpos_T length) {
  T new = (T) MALLOC(sizeof(*new));
  Genomecomp_T high, low, flags;
  Chrpos_T ptr;
  Chrpos_T position;
  int c, i;
  int in_counter = 0;

  new->nblocks = (length+31)/32U;

#ifdef HAVE_SSE2
  new->blocks = (Genomecomp_T *) _mm_malloc((new->nblocks+1)*COMPRESS_BLOCKSIZE * sizeof(Genomecomp_T),16);
#else
  new->blocks = (Genomecomp_T *) MALLOC((new->nblocks+1)*COMPRESS_BLOCKSIZE * sizeof(Genomecomp_T));
#endif

  /* Note that elements of shift_array do not have extra block at beginning */
  new->shift_array = (Genomecomp_T **) MALLOC(32 * sizeof(Genomecomp_T *));
  new->shift_array[0] = (Genomecomp_T *) MALLOC(32*(new->nblocks+1)*COMPRESS_BLOCKSIZE * sizeof(Genomecomp_T));
  new->availp[0] = false;
  for (i = 1; i < 32; i++) {
    new->shift_array[i] = &(new->shift_array[i-1][(new->nblocks+1)*COMPRESS_BLOCKSIZE]);
    new->availp[i] = false;
  }

  ptr = 0;

  position = 0U;
  while (position < length) {

    high = low = flags = 0U;
    in_counter = 0;
    while (position < length && in_counter < 32) {
      c = gbuffer[position++];
      high >>= 1;
      low >>= 1;
      flags >>= 1;

      /* Assume that gbuffer is upper case */
      switch /*(uppercaseCode[c])*/ (c) {
      case 'A': /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'C': /* high |= LEFT_CLEAR; */    low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      case 'G':    high |= LEFT_SET;      /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'T':    high |= LEFT_SET;         low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      default:  /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */    flags |= LEFT_SET;
      }
      in_counter++;
    }
      
    while (in_counter < 32) {
      high >>= 1;
      low >>= 1;
      flags >>= 1;
      in_counter++;
    }

    new->blocks[ptr] = high;
    new->blocks[ptr+1] = low;
    new->blocks[ptr+2] = flags;
    /* new->blocks[ptr+3] = 0U; */

    ptr += COMPRESS_BLOCKSIZE;
  }

  /* Compress_shift will access these values */
  new->blocks[ptr] = 0U;
  new->blocks[ptr+1] = 0U;
  new->blocks[ptr+2] = 0U;
  /* new->blocks[ptr+3] = 0U; */

  assert(ptr+3 < (new->nblocks+1)*COMPRESS_BLOCKSIZE);

  debug1(printf("Compress_new_fwd\n"));
  debug1(Compress_print_blocks(new->blocks,new->nblocks));
  debug1(printf("\n"));

  return new;
}

T
Compress_new_rev (char *gbuffer, Chrpos_T length) {
  T new = (T) MALLOC(sizeof(*new));
  Genomecomp_T high, low, flags;
  Chrpos_T ptr;
  Chrpos_T position;
  int c, i;
  int in_counter = 0;

  new->nblocks = (length+31)/32U;

#ifdef HAVE_SSE2
  new->blocks = (Genomecomp_T *) _mm_malloc((new->nblocks+1)*COMPRESS_BLOCKSIZE * sizeof(Genomecomp_T),16);
#else
  new->blocks = (Genomecomp_T *) MALLOC((new->nblocks+1)*COMPRESS_BLOCKSIZE * sizeof(Genomecomp_T));
#endif

  /* Note that elements of shift_array do not have extra block at beginning */
  new->shift_array = (Genomecomp_T **) MALLOC(32 * sizeof(Genomecomp_T *));
  new->shift_array[0] = (Genomecomp_T *) MALLOC(32*(new->nblocks+1)*COMPRESS_BLOCKSIZE * sizeof(Genomecomp_T));
  new->availp[0] = false;
  for (i = 1; i < 32; i++) {
    new->shift_array[i] = &(new->shift_array[i-1][(new->nblocks+1)*COMPRESS_BLOCKSIZE]);
    new->availp[i] = false;
  }

  ptr = 0;

  position = length;
  while (position > 0) {

    high = low = flags = 0U;
    in_counter = 0;
    while (position > 0 && in_counter < 32) {
      c = gbuffer[--position];
      high >>= 1;
      low >>= 1;
      flags >>= 1;

      /* Assume that gbuffer is upper case */
      switch /*(uppercaseCode[c])*/ (c) {
      case 'T': /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'G': /* high |= LEFT_CLEAR; */    low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      case 'C':    high |= LEFT_SET;      /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'A':    high |= LEFT_SET;         low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      default:  /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */    flags |= LEFT_SET;
      }
      in_counter++;
    }

    while (in_counter < 32) {
      high >>= 1;
      low >>= 1;
      flags >>= 1;
      in_counter++;
    }

    new->blocks[ptr] = high;
    new->blocks[ptr+1] = low;
    new->blocks[ptr+2] = flags;
    /* new->blocks[ptr+3] = 0U; */
    
    ptr += COMPRESS_BLOCKSIZE;
  }

  /* Compress_shift will access these values */
  new->blocks[ptr] = 0U;
  new->blocks[ptr+1] = 0U;
  new->blocks[ptr+2] = 0U;
  /* new->blocks[ptr+3] = 0U; */

  assert(ptr+3 < (new->nblocks+1)*COMPRESS_BLOCKSIZE);

  debug1(printf("Compress_new_rev\n"));
  debug1(Compress_print_blocks(new->blocks,new->nblocks));
  debug1(printf("\n"));

  return new;
}



Genomecomp_T *
Compress_shift (T this, int nshift) {
  Genomecomp_T *shifted;
  int rightshift;
  int ptr;
#ifdef HAVE_SSE2
  __m128i out, current, next;
#endif
#ifdef DEBUG9
  Genomecomp_T high, low, flags;
#endif

  if (this->availp[nshift] == true) {
    return this->shift_array[nshift];

  } else {
    shifted = this->shift_array[nshift];

    /* Shift */
    ptr = this->nblocks*COMPRESS_BLOCKSIZE;
    if (nshift == 0) {
      while (ptr >= 0) {
	memcpy(&(shifted[ptr]),&(this->blocks[ptr]),COMPRESS_BLOCKSIZE*sizeof(Genomecomp_T));
	ptr -= COMPRESS_BLOCKSIZE;
      }

    } else {
      rightshift = 32 - nshift;

#ifdef DEFECTIVE_SSE2_COMPILER
      while (ptr > 0) {
	shifted[ptr+2] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr-2] >> rightshift);
	shifted[ptr+1] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr-3] >> rightshift);
	shifted[ptr] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-4] >> rightshift);
	ptr -= COMPRESS_BLOCKSIZE;
      }

      shifted[2] = this->blocks[2] << nshift;
      shifted[1] = this->blocks[1] << nshift;
      shifted[0] = this->blocks[0] << nshift;

#elif defined(HAVE_SSE2)
      next = _mm_load_si128((__m128i *) &(this->blocks[ptr]));
      while (ptr > 0) {
	current = next;
	next = _mm_load_si128((__m128i *) &(this->blocks[ptr-4]));
	out = _mm_or_si128(_mm_slli_epi32(current,nshift),_mm_srli_epi32(next,rightshift));
	_mm_store_si128((__m128i *) &(shifted[ptr]),out);

#ifdef DEBUG9
	flags = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr-2] >> rightshift);
	low = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr-3] >> rightshift);
	high = (this->blocks[ptr] << nshift) | (this->blocks[ptr-4] >> rightshift);
	/* printf("Comparing %08X with %08X, %08X with %08X, %08X with %08X\n",
	   shifted[ptr],high,shifted[ptr+1],low,shifted[ptr+2],flags); */
	if (shifted[ptr+2] != flags) abort();
	if (shifted[ptr+1] != low) abort();
	if (shifted[ptr] != high) abort();
#endif

	ptr -= COMPRESS_BLOCKSIZE;
      }

      out = _mm_slli_epi32(next,nshift);
      _mm_store_si128((__m128i *) &(shifted[0]),out);

#ifdef DEBUG9
      flags = this->blocks[2] << nshift;
      low = this->blocks[1] << nshift;
      high = this->blocks[0] << nshift;
      /* printf("Comparing %08X with %08X, %08X with %08X, %08X with %08X\n",
	 shifted[0],high,shifted[1],low,shifted[2],flags); */
      if (shifted[2] != flags) abort();
      if (shifted[1] != low) abort();
      if (shifted[0] != high) abort();
#endif

#else
      while (ptr > 0) {
	shifted[ptr+2] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr-2] >> rightshift);
	shifted[ptr+1] = (this->blocks[ptr+1] << nshift) | (this->blocks[ptr-3] >> rightshift);
	shifted[ptr] = (this->blocks[ptr] << nshift) | (this->blocks[ptr-4] >> rightshift);
	ptr -= COMPRESS_BLOCKSIZE;
      }

      shifted[2] = this->blocks[2] << nshift;
      shifted[1] = this->blocks[1] << nshift;
      shifted[0] = this->blocks[0] << nshift;
#endif

    }

    this->availp[nshift] = true;

    debug1(Compress_print_blocks(shifted,this->nblocks+1));
    return shifted;
  }
}


