#ifndef DIBASE_INCLUDED
#define DIBASE_INCLUDED

#include "types.h"
#include "genomicpos.h"

extern int
Dibase_count_mismatches_limit (int *color_nmismatches, char *queryuc_ptr, int pos5, int pos3, UINT4 *genome_blocks,
			       Genomicpos_T startpos, Genomicpos_T endpos, int max_mismatches);
extern int
Dibase_count_mismatches_substring (int *color_nmismatches, char *queryuc_ptr, int pos5, int pos3, UINT4 *genome_blocks,
				   Genomicpos_T startpos, Genomicpos_T endpos);
extern int
Dibase_mismatches_left (int *mismatch_positions, int *colordiffs, int max_mismatches, char *queryuc_ptr,
			int pos5, int pos3, UINT4 *genome_blocks, Genomicpos_T startpos, Genomicpos_T endpos);
extern int
Dibase_mismatches_right (int *mismatch_positions, int *colordiffs, int max_mismatches, char *queryuc_ptr,
			 int pos5, int pos3, UINT4 *genome_blocks, Genomicpos_T startpos, Genomicpos_T endpos);
#endif

