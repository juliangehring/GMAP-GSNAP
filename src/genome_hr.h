/* $Id: genome_hr.h 99749 2013-06-27 21:03:51Z twu $ */
#ifndef GENOME_HR_INCLUDED
#define GENOME_HR_INCLUDED
#include "types.h"
#include "mode.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "compress.h"
#include "chrnum.h"

extern void
Genome_hr_setup (Genomecomp_T *ref_blocks_in, Genomecomp_T *snp_blocks_in,
		 bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
		 Mode_T mode_in);

extern void
Genome_hr_user_setup (UINT4 *ref_blocks_in,
		      bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
		      Mode_T mode_in);

/* Procedures for indexdb */
extern int
Genome_read_gamma (Positionsptr_T **ptr, int ctr, Positionsptr_T *cum);
extern Positionsptr_T
Genome_offsetptr_from_gammas (Positionsptr_T *end0, Gammaptr_T *gammaptrs, Positionsptr_T *offsetscomp,
			      Blocksize_T offsetscomp_blocksize, Storedoligomer_T oligo);
extern Positionsptr_T
Genome_offsetptr_only_from_gammas (Gammaptr_T *gammaptrs, Positionsptr_T *offsetscomp,
				   Blocksize_T offsetscomp_blocksize, Storedoligomer_T oligo);

#ifdef WORDS_BIGENDIAN
extern int
Genome_read_gamma_bigendian (Positionsptr_T **ptr, int ctr, Positionsptr_T *cum);
extern Positionsptr_T
Genome_offsetptr_from_gammas_bigendian (Positionsptr_T *end0, Gammaptr_T *gammaptrs, Positionsptr_T *offsetscomp,
					Blocksize_T offsetscomp_blocksize, Storedoligomer_T oligo);
extern Positionsptr_T
Genome_offsetptr_only_from_gammas_bigendian (Gammaptr_T *gammaptrs, Positionsptr_T *offsetscomp,
					     Blocksize_T offsetscomp_blocksize, Storedoligomer_T oligo);
#endif


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
extern Genomecomp_T
Genome_query_shift_fragment_right (Genomecomp_T *flags, Genomecomp_T *mask, Compress_T query_compress, int pos5, int pos3);
extern Genomecomp_T
Genome_query_shift_fragment_left (Genomecomp_T *flags, Genomecomp_T *mask, Compress_T query_compress, int pos5, int pos3);
extern int
Genome_count_mismatches_fragment (Genomecomp_T query_shifted, Genomecomp_T flags, Genomecomp_T mask,
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

extern int
Genome_donor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			Univcoord_T left, int pos5, int pos3);

extern int
Genome_acceptor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			   Univcoord_T left, int pos5, int pos3);

extern int
Genome_antidonor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			    Univcoord_T left, int pos5, int pos3);

extern int
Genome_antiacceptor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			       Univcoord_T left, int pos5, int pos3);


Chrpos_T
Genome_prev_donor_position (Chrpos_T pos, Chrpos_T prevpos,
			    Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp);
Chrpos_T
Genome_prev_acceptor_position (Chrpos_T pos, Chrpos_T prevpos,
			       Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp);
Chrpos_T
Genome_prev_antidonor_position (Chrpos_T pos, Chrpos_T prevpos,
				Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp);
Chrpos_T
Genome_prev_antiacceptor_position (Chrpos_T pos, Chrpos_T prevpos,
				   Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp);


#if 0
extern void
Genome_last_donor_positions (int *last_position, Univcoord_T genomicstart, int margin5, int margin3, int genomiclength,
			     bool plusp);

extern void
Genome_last_acceptor_positions (int *last_position, Univcoord_T genomicstart, int margin5, int margin3, int genomiclength,
				bool plusp);

extern void
Genome_last_antidonor_positions (int *last_position, Univcoord_T genomicstart, int margin5, int margin3, int genomiclength,
				 bool plusp);

extern void
Genome_last_antiacceptor_positions (int *last_position, Univcoord_T genomicstart, int margin5, int margin3, int genomiclength,
				    bool plusp);
#endif


#endif

