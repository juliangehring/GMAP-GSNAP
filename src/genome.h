/* $Id: genome.h 101462 2013-07-15 15:07:35Z twu $ */
#ifndef GENOME_INCLUDED
#define GENOME_INCLUDED

#include "bool.h"
#include "access.h"
#include "types.h"
#include "genomicpos.h"
#include "sequence.h"
#include "iit-read-univ.h"
#include "chrnum.h"
#include "mode.h"


#define OUTOFBOUNDS '*'

typedef enum {GENOME_OLIGOS, GENOME_BITS} Genometype_T;

#define T Genome_T
typedef struct T *T;

extern void
Genome_free (T *old);
extern Genomecomp_T *
Genome_blocks (T this);
extern Univcoord_T
Genome_totallength (T this);
extern T
Genome_new (char *genomesubdir, char *fileroot, char *snps_root,
	    Genometype_T genometype, bool genome_lc_p, Access_mode_T access);

extern void
Genome_setup (T genome_in, T genomealt_in, Mode_T mode_in, int circular_typeint_in);
extern void
Genome_user_setup (Genomecomp_T *genome_blocks_in);

extern void
Genome_uncompress_mmap (char *gbuffer1, Genomecomp_T *blocks, Univcoord_T startpos, 
			Univcoord_T endpos);
extern bool
Genome_fill_buffer (Chrnum_T *chrnum, int *nunknowns, T this, Univcoord_T left, Chrpos_T length, char *gbuffer1,
		    Univ_IIT_T chromosome_iit);
extern void
Genome_fill_buffer_simple (T this, Univcoord_T left, Chrpos_T length, char *gbuffer1);
extern void
Genome_fill_buffer_blocks (Univcoord_T left, Chrpos_T length, char *gbuffer1);
extern void
Genome_fill_buffer_blocks_noterm (Univcoord_T left, Chrpos_T length, char *gbuffer1, char *gbuffer2);
extern void
Genome_fill_buffer_simple_alt (T genome, T genomealt, Univcoord_T left, Chrpos_T length, char *gbuffer1);
extern void
Genome_fill_buffer_nucleotides (T this, Univcoord_T left, Chrpos_T length, unsigned char *gbuffer);
extern void
Genome_fill_buffer_int_string (T this, Univcoord_T left, Chrpos_T length, unsigned char *gbuffer);
extern char
Genome_get_char (T this, Univcoord_T left);
extern char
Genome_get_char_blocks (char *charalt, Univcoord_T left);
extern char *
Genome_get_segment_blocks_right (char **segmentalt, Univcoord_T left, Chrpos_T length, Univcoord_T chrhigh,
				 bool revcomp);
extern char *
Genome_get_segment_blocks_left (char **segmentalt, Univcoord_T left, Chrpos_T length, Univcoord_T chroffset,
				bool revcomp);
extern Sequence_T
Genome_get_segment (T this, Univcoord_T left, Chrpos_T length, Univ_IIT_T chromosome_iit,
		    bool revcomp);
extern Sequence_T
Genome_get_segment_alt (T this, Univcoord_T left, Chrpos_T length, Univ_IIT_T chromosome_iit,
			bool revcomp);
extern Sequence_T
Genome_get_segment_snp (T this, Univcoord_T left, Chrpos_T length, Univ_IIT_T chromosome_iit,
			bool revcomp);
extern int
Genome_next_char (T this);
extern int
Genome_ntcounts (int *na, int *nc, int *ng, int *nt,
		 T this, Univcoord_T left, Chrpos_T length);

#undef T
#endif
