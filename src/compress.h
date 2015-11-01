/* $Id: compress.h 102728 2013-07-24 22:48:20Z twu $ */
#ifndef COMPRESS_INCLUDED
#define COMPRESS_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "types.h"
#include "genomicpos.h"

/* Store blocks every 4 uints, so they are on 128-bit boundaries.
   This will waste one out of every 4 uints, but allows for use of
   SIMD in Compress_shift */

#define COMPRESS_BLOCKSIZE 4	/* 4 unsigned ints per block */


#define T Compress_T
typedef struct T *T;

extern void
Compress_free (T *old);
extern void
Compress_print (T this);
extern int
Compress_nblocks (T this);
extern void
Compress_print_blocks (Genomecomp_T *blocks, int nblocks);
extern T
Compress_new_fwd (char *gbuffer, Chrpos_T length);
extern T
Compress_new_rev (char *gbuffer, Chrpos_T length);
extern Genomecomp_T *
Compress_shift (T this, int nshift);

#undef T
#endif

