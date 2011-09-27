/* $Id: genome_hr.h,v 1.18 2010-07-19 20:27:49 twu Exp $ */
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
Genome_count_mismatches_limit (int *ncolordiffs, char *queryuc_ptr, Compress_T query_compress,
			       UINT4 *blocks, UINT4 *snp_blocks, Genomicpos_T left, int pos5, int pos3,
			       int max_mismatches, bool dibasep, bool cmetp, bool plusp);
extern int
Genome_count_mismatches_substring (int *ncolordiffs, char *queryuc_ptr, Compress_T query_compress,
				   UINT4 *blocks, UINT4 *snp_blocks, Genomicpos_T left, int pos5, int pos3,
				   bool dibasep, bool cmetp, bool plusp);
extern int
Genome_mismatches_left (int *mismatch_positions, int *colordiffs, int max_mismatches,
			char *queryuc_ptr, Compress_T query_compress, UINT4 *blocks, UINT4 *snp_blocks,
			Genomicpos_T left, int pos5, int pos3, bool dibasep, bool cmetp, bool plusp);
extern int
Genome_mismatches_right (int *mismatch_positions, int *colordiffs, int max_mismatches,
			 char *queryuc_ptr, Compress_T query_compress, UINT4 *blocks, UINT4 *snp_blocks,
			 Genomicpos_T left, int pos5, int pos3, bool dibasep, bool cmetp, bool plusp);

extern int
Genome_trim_left (char *queryuc_ptr, Compress_T query_compress,
		  UINT4 *blocks, UINT4 *snp_blocks, Genomicpos_T left, int pos5, int pos3,
		  bool dibasep, bool cmetp, bool plusp);

extern int
Genome_trim_right (char *queryuc_ptr, Compress_T query_compress,
		   UINT4 *blocks, UINT4 *snp_blocks, Genomicpos_T left, int pos5, int pos3,
		   bool dibasep, bool cmetp, bool plusp);

extern char
Genome_get_dinucleotide (char *altdinucl, UINT4 *ref_blocks, UINT4 *alt_blocks,
			 Genomicpos_T pos);

#endif

