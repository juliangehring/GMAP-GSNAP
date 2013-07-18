/* $Id: compress.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef COMPRESS_INCLUDED
#define COMPRESS_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "types.h"
#include "genomicpos.h"

#define T Compress_T
typedef struct T *T;

extern int
Compress_get_char (FILE *sequence_fp, Univcoord_T position, bool uncompressedp);
extern void
Compress_compress (FILE *fp);
extern void
Compress_uncompress (FILE *fp, int wraplength);
extern int
Compress_update_file (int nbadchars, FILE *fp, char *gbuffer, Univcoord_T startpos,
		      Univcoord_T endpos, int index1part);
extern int
Compress_update_memory (int nbadchars, Genomecomp_T *genomecomp, char *gbuffer, Univcoord_T startpos,
			Univcoord_T endpos);
extern void
Compress_free (T *old);
extern void
Compress_print (T this);
extern int
Compress_nblocks (T this);
extern T
Compress_new (char *gbuffer, Chrpos_T length, bool plusp);
extern Genomecomp_T *
Compress_shift (T this, int nshift);

#undef T
#endif

