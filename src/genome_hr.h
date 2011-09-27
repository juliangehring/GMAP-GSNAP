/* $Id: genome_hr.h 34786 2011-02-06 22:20:41Z twu $ */
#ifndef GENOME_HR_INCLUDED
#define GENOME_HR_INCLUDED
#include "types.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "compress.h"
#include "chrnum.h"


extern int
Genome_count_mismatches (Compress_T query_compress, UINT4 *blocks, Genomicpos_T left, Genomicpos_T left_plus_length);
extern int
Genome_count_mismatches_limit (int *ncolordiffs, char *query, Compress_T query_compress,
			       UINT4 *blocks, UINT4 *snp_blocks, Genomicpos_T left, int pos5, int pos3,
			       int max_mismatches, bool dibasep, bool cmetp, bool plusp);
extern int
Genome_count_mismatches_substring (int *ncolordiffs, char *query, Compress_T query_compress,
				   UINT4 *blocks, UINT4 *snp_blocks, Genomicpos_T left, int pos5, int pos3,
				   bool dibasep, bool cmetp, bool plusp);
extern UINT4
Genome_query_shift_fragment_right (UINT4 *flags, UINT4 *mask, Compress_T query_compress, int pos5, int pos3);
extern UINT4
Genome_query_shift_fragment_left (UINT4 *flags, UINT4 *mask, Compress_T query_compress, int pos5, int pos3);
extern int
Genome_count_mismatches_fragment (UINT4 query_shifted, UINT4 flags, UINT4 mask,
				  UINT4 ref_fragment, UINT4 alt_fragment);

extern int
Genome_mismatches_left (int *mismatch_positions, int *colordiffs, int max_mismatches,
			char *query, Compress_T query_compress, UINT4 *blocks, UINT4 *snp_blocks,
			Genomicpos_T left, int pos5, int pos3, bool dibasep, bool cmetp, bool plusp);
extern int
Genome_mismatches_right (int *mismatch_positions, int *colordiffs, int max_mismatches,
			 char *query, Compress_T query_compress, UINT4 *blocks, UINT4 *snp_blocks,
			 Genomicpos_T left, int pos5, int pos3, bool dibasep, bool cmetp, bool plusp);

extern int
Genome_mark_mismatches (char *genomic, int querylength, Compress_T query_compress, UINT4 *blocks, UINT4 *snp_blocks,
			Genomicpos_T left, int pos5, int pos3, int mismatch_offset, bool dibasep, bool cmetp, bool plusp);

extern int
Genome_trim_left (Compress_T query_compress,
		  UINT4 *blocks, UINT4 *snp_blocks, Genomicpos_T left, int pos5, int pos3,
		  bool dibasep, bool cmetp, bool plusp);

extern int
Genome_trim_right (Compress_T query_compress,
		   UINT4 *blocks, UINT4 *snp_blocks, Genomicpos_T left, int pos5, int pos3,
		   bool dibasep, bool cmetp, bool plusp);

extern char
Genome_get_dinucleotide (char *altdinucl, UINT4 *ref_blocks, UINT4 *alt_blocks,
			 Genomicpos_T pos);

extern int
Genome_donor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			UINT4 *blocks, Genomicpos_T left, int pos5, int pos3);

extern int
Genome_acceptor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			   UINT4 *blocks, Genomicpos_T left, int pos5, int pos3);

extern int
Genome_antidonor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			    UINT4 *blocks, Genomicpos_T left, int pos5, int pos3);

extern int
Genome_antiacceptor_positions (int *site_positions, int *site_knowni, int *knownpos, int *knowni,
			       UINT4 *blocks, Genomicpos_T left, int pos5, int pos3);


#endif

