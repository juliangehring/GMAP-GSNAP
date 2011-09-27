/* $Id: compress.h,v 1.2 2005/02/14 12:51:56 twu Exp $ */
#ifndef COMPRESS_INCLUDED
#define COMPRESS_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "types.h"
#include "genomicpos.h"

extern int
Compress_get_char (FILE *sequence_fp, Genomicpos_T position, bool uncompressedp);
extern void
Compress_compress (FILE *fp);
extern void
Compress_uncompress (FILE *fp);
extern void
Compress_update_file (FILE *fp, char *gbuffer, Genomicpos_T startpos,
		      Genomicpos_T endpos);
extern void
Compress_update_memory (UINT4 *genomeseq, char *gbuffer, Genomicpos_T startpos,
			Genomicpos_T endpos);

#endif
