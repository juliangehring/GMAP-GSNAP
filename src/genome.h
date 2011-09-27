/* $Id: genome.h,v 1.23 2007/02/18 20:53:08 twu Exp $ */
#ifndef GENOME_INCLUDED
#define GENOME_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "sequence.h"
#include "iit-read.h"

#define T Genome_T
typedef struct T *T;

extern void
Genome_free (T *old);
extern T
Genome_new (char *genomesubdir, char *fileroot, bool uncompressedp, bool batchp);
extern void
Genome_replace_x (void);
extern Sequence_T
Genome_get_segment (T this, Genomicpos_T left, Genomicpos_T length, 
		    bool revcomp, char *gbuffer1, char *gbuffer2, int gbufferlen);
extern void
Genome_fill_buffer (T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1);
extern int
Genome_next_char (T this);
extern Sequence_T
Genome_patch_strain (int *indices, int nindices, IIT_T altstrain_iit, 
		     Genomicpos_T refL, Genomicpos_T reflen,
		     bool revcomp, char *gbuffer1, char *gbuffer2, char *gbuffer3,
		     int gbuffer3len);

#undef T
#endif
