/* $Id: genome.h 79247 2012-11-15 19:09:06Z twu $ */
#ifndef GENOME_INCLUDED
#define GENOME_INCLUDED

#include "bool.h"
#include "access.h"
#include "types.h"
#include "genomicpos.h"
#include "sequence.h"
#include "iit-read.h"
#include "chrnum.h"
#include "mode.h"


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
Genome_new (char *genomesubdir, char *fileroot, char *snps_root, bool genome_lc_p, Access_mode_T access);

extern void
Genome_setup (T genome_in, T genomealt_in, Mode_T mode_in, int circular_typeint_in);
extern void
Genome_user_setup (UINT4 *genome_blocks_in);

extern void
Genome_uncompress_mmap (char *gbuffer1, UINT4 *blocks, Genomicpos_T startpos, 
			Genomicpos_T endpos);
extern bool
Genome_fill_buffer (Chrnum_T *chrnum, int *nunknowns, T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1,
		    IIT_T chromosome_iit);
extern void
Genome_fill_buffer_simple (T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1);
extern void
Genome_fill_buffer_blocks (Genomicpos_T left, Genomicpos_T length, char *gbuffer1);
extern void
Genome_fill_buffer_blocks_noterm (Genomicpos_T left, Genomicpos_T length, char *gbuffer1, char *gbuffer2);
extern void
Genome_fill_buffer_simple_alt (T genome, T genomealt, Genomicpos_T left, Genomicpos_T length, char *gbuffer1);
extern void
Genome_fill_buffer_nucleotides (T this, Genomicpos_T left, Genomicpos_T length, unsigned char *gbuffer);
extern char
Genome_get_char (T this, Genomicpos_T left);
extern char
Genome_get_char_blocks (char *charalt, Genomicpos_T left);
extern Sequence_T
Genome_get_segment (T this, Genomicpos_T left, Genomicpos_T length, IIT_T chromosome_iit,
		    bool revcomp);
extern Sequence_T
Genome_get_segment_alt (T this, Genomicpos_T left, Genomicpos_T length, IIT_T chromosome_iit,
			bool revcomp);
extern Sequence_T
Genome_get_segment_snp (T this, Genomicpos_T left, Genomicpos_T length, IIT_T chromosome_iit,
			bool revcomp);
extern int
Genome_next_char (T this);
extern int
Genome_ntcounts (int *na, int *nc, int *ng, int *nt,
		 T this, Genomicpos_T left, Genomicpos_T length);
extern Sequence_T
Genome_patch_strain (int *indices, int nindices, IIT_T altstrain_iit, 
		     Genomicpos_T refL, Genomicpos_T reflen,
		     bool revcomp, char *gbuffer1, char *gbuffer2, char *gbuffer3,
		     int gbuffer3len);

#undef T
#endif
