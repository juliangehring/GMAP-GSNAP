/* $Id: stage2.h,v 1.68 2010-07-10 19:21:13 twu Exp $ */
#ifndef STAGE2_INCLUDED
#define STAGE2_INCLUDED
#include "bool.h"
#include "list.h"
#include "oligoindex.h"
#include "pairpool.h"
#include "diagpool.h"
#include "stopwatch.h"

#define T Stage2_T
typedef struct T *T;

extern List_T
Stage2_path (T this);
extern int
Stage2_ncanonical (T this);
extern int
Stage2_nnoncanonical (T this);
extern int
Stage2_cdna_direction (T this);
extern double
Stage2_diag_runtime (T this);
extern double
Stage2_align_runtime (T this);
extern int
Stage2_indexsize (T this);
extern int
Stage2_pathlength (T this);
	    
extern void
Stage2_free (T *old);

extern int
Stage2_scan (int *stage2_source,
	     char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,
	     char *genomicseg_ptr, char *genomicuc_ptr, int genomiclength, int genomic_offset,
	     Oligoindex_T *oligoindices, int noligoindices,
	     Diagpool_T diagpool, bool debug_graphic_p, bool diagnosticp);

extern List_T
Stage2_compute (int *stage2_source, int *stage2_indexsize,
		char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,
		char *genomicseg_ptr, char *genomicuc_ptr, int genomiclength, int genomic_offset,
		Oligoindex_T *oligoindices, int noligoindices, double proceed_pctcoverage,
		Pairpool_T pairpool, Diagpool_T diagpool, int sufflookback, int nsufflookback,
		int maxintronlen, bool localp, bool skip_repetitive_p, bool debug_graphic_p,
		bool diagnosticp, Stopwatch_T stopwatch, bool diag_debug);

#undef T
#endif

