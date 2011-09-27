/* $Id: compress.h 46784 2011-09-07 23:06:49Z twu $ */
#ifndef COMPRESS_INCLUDED
#define COMPRESS_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "types.h"
#include "genomicpos.h"

#define T Compress_T
typedef struct T *T;

extern int
Compress_get_char (FILE *sequence_fp, Genomicpos_T position, bool uncompressedp);
extern void
Compress_compress (FILE *fp);
extern void
Compress_uncompress (FILE *fp, int wraplength);
extern int
Compress_update_file (int nbadchars, FILE *fp, char *gbuffer, Genomicpos_T startpos,
		      Genomicpos_T endpos, int index1part);
extern int
Compress_update_memory (int nbadchars, UINT4 *genomeseq, char *gbuffer, Genomicpos_T startpos,
			Genomicpos_T endpos);
extern void
Compress_free (T *old);
extern void
Compress_print (T this);
extern int
Compress_nblocks (T this);
extern T
Compress_new (char *gbuffer, int length, bool plusp);
extern UINT4 *
Compress_shift (T this, int nshift);

#undef T
#endif

