/* $Id: genome_hr.h 49447 2011-10-09 17:38:10Z twu $ */
#ifndef GENOME_HR_INCLUDED
#define GENOME_HR_INCLUDED
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "compress.h"
#include "chrnum.h"
#include "indexdbdef.h"		/* For Storedoligomer_T */

extern void
Genome_hr_setup (UINT4 *ref_blocks_in, UINT4 *snp_blocks_in,
		 bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
		 Mode_T mode_in);

/* Procedures for indexdb */
extern int
Genome_read_gamma (unsigned int **ptr, int ctr, unsigned int *cum);
extern Positionsptr_T
Genome_offsetptr_from_gammas (Positionsptr_T *end0, UINT4 *gammaptrs, Positionsptr_T *offsetscomp,
			      unsigned int offsets_blocksize, Storedoligomer_T oligo);
extern Positionsptr_T
Genome_offsetptr_only_from_gammas (UINT4 *gammaptrs, Positionsptr_T *offsetscomp,
				   unsigned int offsets_blocksize, Storedoligomer_T oligo);

#ifdef WORDS_BIGENDIAN
extern int
Genome_read_gamma_bigendian (unsigned int **ptr, int ctr, unsigned int *cum);
extern Positionsptr_T
Genome_offsetptr_from_gammas_bigendian (Positionsptr_T *end0, UINT4 *gammaptrs, Positionsptr_T *offsetscomp,
					unsigned int offsets_blocksize, Storedoligomer_T oligo);
extern Positionsptr_T
Genome_offsetptr_only_from_gammas_bigendian (UINT4 *gammaptrs, Positionsptr_T *offsetscomp,
					     unsigned int offsets_blocksize, Storedoligomer_T oligo);
#endif


extern int
Genome_count_mismatches (Compress_T query_compress, Genomicpos_T left, Genomicpos_T left_plus_length);
extern int
Genome_count_mismatches_limit (Compress_T query_compress, Genomicpos_T left, int pos5, int pos3,
			       int max_mismatches, bool plusp, int genestrand);
extern int
Genome_count_mismatches_substring_ref (Compress_T query_compress, Genomicpos_T left, int pos5, int pos3,
				       bool plusp, int genestrand);
extern int
Genome_count_mismatches_substring (Compress_T query_compress, Genomicpos_T left, int pos5, int pos3,
				   bool plusp, int genestrand);
extern UINT4
Genome_query_shift_fragment_right (UINT4 *flags, UINT4 *mask, Compress_T query_compress, int pos5, int pos3);
extern UINT4
Genome_query_shift_fragment_left (UINT4 *flags, UINT4 *mask, Compress_T query_compress, int pos5, int pos3);
extern int
Genome_count_mismatches_fragment (UINT4 query_shifted, UINT4 flags, UINT4 mask,
				  UINT4 ref_fragment, UINT4 alt_fragment);

extern int
Genome_mismatches_left (int *mismatch_positions, int max_mismatches, Compress_T query_compress,
			Genomicpos_T left, int pos5, int pos3, bool plusp, int genestrand);
extern int
Genome_mismatches_right (int *mismatch_positions, int max_mismatches, Compress_T query_compress,
			 Genomicpos_T left, int pos5, int pos3, bool plusp, int genestrand);

extern int
Genome_mark_mismatches_ref (char *genomic, int querylength, Compress_T query_compress,
			    Genomicpos_T left, int pos5, int pos3, int mismatch_offset,
			    bool plusp, int genestrand);
extern int
Genome_mark_mismatches (char *genomic, int querylength, Compress_T query_compress,
			Genomicpos_T left, int pos5, int pos3, int mismatch_offset, bool plusp, int genestrand);

extern int
Genome_trim_left (Compress_T query_compress, Genomicpos_T left, int pos5, int pos3,
		  bool plusp, int genestrand);

extern int
Genome_trim_right (Compress_T query_compress, Genomicpos_T left, int pos5, int pos3,
		   bool plusp, int genestrand);

extern char
Genome_get_dinucleotide (char *altdinucl, Genomicpos_T pos);

extern int
Genome_donor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			Genomicpos_T left, int pos5, int pos3);

extern int
Genome_acceptor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			   Genomicpos_T left, int pos5, int pos3);

extern int
Genome_antidonor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			    Genomicpos_T left, int pos5, int pos3);

extern int
Genome_antiacceptor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			       Genomicpos_T left, int pos5, int pos3);


extern int
Genome_prev_donor_position (int pos, Genomicpos_T genomicstart, Genomicpos_T genomicend, int pos5, bool plusp);
extern int
Genome_prev_acceptor_position (int pos, Genomicpos_T genomicstart, Genomicpos_T genomicend, int pos5, bool plusp);
extern int
Genome_prev_antidonor_position (int pos, Genomicpos_T genomicstart, Genomicpos_T genomicend, int pos5, bool plusp);
extern int
Genome_prev_antiacceptor_position (int pos, Genomicpos_T genomicstart, Genomicpos_T genomicend, int pos5, bool plusp);



extern void
Genome_last_donor_positions (int *last_position, Genomicpos_T genomicstart, int margin5, int margin3, int genomiclength,
			     bool plusp);

extern void
Genome_last_acceptor_positions (int *last_position, Genomicpos_T genomicstart, int margin5, int margin3, int genomiclength,
				bool plusp);

extern void
Genome_last_antidonor_positions (int *last_position, Genomicpos_T genomicstart, int margin5, int margin3, int genomiclength,
				 bool plusp);

extern void
Genome_last_antiacceptor_positions (int *last_position, Genomicpos_T genomicstart, int margin5, int margin3, int genomiclength,
				    bool plusp);


#endif

