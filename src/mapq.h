/* $Id: mapq.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef MAPQ_INCLUDED
#define MAPQ_INCLUDED

#include "types.h"
#include "compress.h"
#include "genomicpos.h"


#define MAX_QUALITY_SCORE_INPUT 96	/* Was 40 */
#define MAX_QUALITY_SCORE 40

extern void
MAPQ_init (int quality_score_adj_in);
extern int
MAPQ_max_quality_score (char *quality_string, int querylength);
extern double
MAPQ_loglik_exact (char *quality_string, int querystart, int queryend);
extern double
MAPQ_loglik (Compress_T query_compress, Univcoord_T left, int querystart, int queryend,
	     int querylength, char *quality_string, bool plusp, int genestrand);

#endif

