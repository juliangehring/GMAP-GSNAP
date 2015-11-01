/* $Id: genome_hr.h 101271 2013-07-12 02:44:39Z twu $ */
#ifndef GENOME_HR_INCLUDED
#define GENOME_HR_INCLUDED
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "compress.h"


extern void
Genome_hr_setup (Genomecomp_T *ref_blocks_in, Genomecomp_T *snp_blocks_in,
		 bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
		 Mode_T mode, bool genomebits_avail_p);

extern void
Genome_hr_user_setup (UINT4 *ref_blocks_in,
		      bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
		      Mode_T mode_in);

/* Procedures for indexdb */
extern int
Genome_read_gamma (Positionsptr_T **ptr, int ctr, Positionsptr_T *cum);
extern Positionsptr_T
Genome_offsetptr_from_gammas (Positionsptr_T *end0, Gammaptr_T *gammaptrs, Offsetscomp_T *offsetscomp,
			      Blocksize_T offsetscomp_blocksize, Storedoligomer_T oligo);
extern Positionsptr_T
Genome_offsetptr_only_from_gammas (Gammaptr_T *gammaptrs, Offsetscomp_T *offsetscomp,
				   Blocksize_T offsetscomp_blocksize, Storedoligomer_T oligo);

#ifdef WORDS_BIGENDIAN
extern int
Genome_read_gamma_bigendian (Positionsptr_T **ptr, int ctr, Positionsptr_T *cum);
extern Positionsptr_T
Genome_offsetptr_from_gammas_bigendian (Positionsptr_T *end0, Gammaptr_T *gammaptrs, Offsetscomp_T *offsetscomp,
					Blocksize_T offsetscomp_blocksize, Storedoligomer_T oligo);
extern Positionsptr_T
Genome_offsetptr_only_from_gammas_bigendian (Gammaptr_T *gammaptrs, Offsetscomp_T *offsetscomp,
					     Blocksize_T offsetscomp_blocksize, Storedoligomer_T oligo);
#endif


extern int
Genome_consecutive_matches_rightward (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				      bool plusp, int genestrand);
extern int
Genome_consecutive_matches_leftward (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				     bool plusp, int genestrand);
extern int
Genome_count_mismatches (Compress_T query_compress, Univcoord_T left, Univcoord_T left_plus_length);
extern int
Genome_count_mismatches_limit (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
			       int max_mismatches, bool plusp, int genestrand);
extern int
Genome_count_mismatches_substring_ref (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				       bool plusp, int genestrand);
extern int
Genome_count_mismatches_substring (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
				   bool plusp, int genestrand);
extern int
Genome_count_mismatches_fragment_left (Compress_T query_compress, int pos5, int pos3,
				       Genomecomp_T ref_fragment, Genomecomp_T alt_fragment);
extern int
Genome_count_mismatches_fragment_right (Compress_T query_compress, int pos5, int pos3,
					Genomecomp_T ref_fragment, Genomecomp_T alt_fragment);

extern int
Genome_mismatches_left (int *mismatch_positions, int max_mismatches, Compress_T query_compress,
			Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand);
extern int
Genome_mismatches_right (int *mismatch_positions, int max_mismatches, Compress_T query_compress,
			 Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand);

extern int
Genome_mismatches_left_trim (int *mismatch_positions, int max_mismatches, Compress_T query_compress,
			     Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand);
extern int
Genome_mismatches_right_trim (int *mismatch_positions, int max_mismatches, Compress_T query_compress,
			      Univcoord_T left, int pos5, int pos3, bool plusp, int genestrand);

extern int
Genome_mark_mismatches_ref (char *genomic, int querylength, Compress_T query_compress,
			    Univcoord_T left, int pos5, int pos3, int mismatch_offset,
			    bool plusp, int genestrand);
extern int
Genome_mark_mismatches (char *genomic, int querylength, Compress_T query_compress,
			Univcoord_T left, int pos5, int pos3, int mismatch_offset, bool plusp, int genestrand);

extern int
Genome_trim_left (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
		  bool plusp, int genestrand);

extern int
Genome_trim_right (Compress_T query_compress, Univcoord_T left, int pos5, int pos3,
		   bool plusp, int genestrand);

extern char
Genome_get_dinucleotide (char *altdinucl, Univcoord_T pos);


#endif

