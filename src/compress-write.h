/* $Id: compress-write.h 102137 2013-07-19 23:05:33Z twu $ */
#ifndef COMPRESS_WRITE_INCLUDED
#define COMPRESS_WRITE_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "types.h"		/* Needed also for HAVE_64_BIT */
#include "genomicpos.h"

extern int
Compress_get_char (FILE *sequence_fp, Univcoord_T position, bool uncompressedp);
extern int
Compress_update_file (int nbadchars, FILE *fp, char *gbuffer, Univcoord_T startpos,
		      Univcoord_T endpos, int index1part);
extern int
Compress_update_memory (int nbadchars, Genomecomp_T *genomecomp, char *gbuffer, Univcoord_T startpos,
			Univcoord_T endpos);
extern void
Compress_unshuffle (FILE *out, FILE *in);
extern Genomecomp_T *
Compress_create_blocks_comp (char *genomicseg, Univcoord_T genomelength);
extern Genomecomp_T *
Compress_create_blocks_bits (Genomecomp_T *genomecomp, Univcoord_T genomelength);


#undef T
#endif

