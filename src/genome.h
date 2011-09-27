/* $Id: genome.h,v 1.35 2009/08/29 00:31:30 twu Exp $ */
#ifndef GENOME_INCLUDED
#define GENOME_INCLUDED
#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "sequence.h"
#include "iit-read.h"
#include "chrnum.h"

#define OUTOFBOUNDS '*'

#define T Genome_T
typedef struct T *T;

extern void
Genome_free (T *old);
extern UINT4 *
Genome_blocks (T this);
extern Genomicpos_T
Genome_totallength (T this);
extern T
Genome_new (char *genomesubdir, char *fileroot, char *snps_root, bool genome_lc_p, bool batchp);
extern void
Genome_uncompress_mmap (char *gbuffer1, UINT4 *blocks, Genomicpos_T startpos, 
			Genomicpos_T endpos, const char defaultchars[],
			const char flagchars[]);
extern bool
Genome_fill_buffer (Chrnum_T *chrnum, int *nunknowns, T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1,
		    IIT_T chromosome_iit);
extern void
Genome_fill_buffer_simple (T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1);
extern Sequence_T
Genome_get_segment (T this, Genomicpos_T left, Genomicpos_T length, IIT_T chromosome_iit,
		    bool revcomp, char *gbuffer1, char *gbuffer2, int gbufferlen);
extern Sequence_T
Genome_get_segment_alt (T this, Genomicpos_T left, Genomicpos_T length, IIT_T chromosome_iit,
			bool revcomp, char *gbuffer1, char *gbuffer2, int gbufferlen);
extern Sequence_T
Genome_get_segment_snp (T this, Genomicpos_T left, Genomicpos_T length, IIT_T chromosome_iit,
			bool revcomp, char *gbuffer1, char *gbuffer2, int gbufferlen);
extern int
Genome_next_char (T this);
extern Sequence_T
Genome_patch_strain (int *indices, int nindices, IIT_T altstrain_iit, 
		     Genomicpos_T refL, Genomicpos_T reflen,
		     bool revcomp, char *gbuffer1, char *gbuffer2, char *gbuffer3,
		     int gbuffer3len);

#undef T
#endif
