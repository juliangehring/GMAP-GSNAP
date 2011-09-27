/* $Id: mapq.h 35469 2011-02-21 16:59:57Z twu $ */
#ifndef MAPQ_INCLUDED
#define MAPQ_INCLUDED

#include "types.h"
#include "compress.h"
#include "genomicpos.h"


extern void
MAPQ_init (int quality_score_adj_input);
extern int
MAPQ_max_quality_score (char *quality_string, int querylength);
extern double
MAPQ_loglik_exact (char *quality_string, int querystart, int queryend);
extern double
MAPQ_loglik (char *query, Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
	     Genomicpos_T left, int querystart, int queryend, int querylength,
	     char *quality_string, bool dibasep, bool cmetp, bool plusp);

#endif

