/* $Id: stage2.h 79302 2012-11-15 23:55:58Z twu $ */
#ifndef STAGE2_INCLUDED
#define STAGE2_INCLUDED
#include "bool.h"
#include "list.h"
#include "oligoindex.h"
#include "pairpool.h"
#include "diagpool.h"
#include "stopwatch.h"
#include "mode.h"

#define T Stage2_T
typedef struct T *T;

extern void
Stage2_setup (bool splicingp_in, int suboptimal_score_start_in, int suboptimal_score_end_in,
	      Mode_T mode_in, bool snps_p_in);
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
Stage2_scan (int *stage2_source, char *queryuc_ptr, int querylength,
#ifdef PMAP
	     char *genomicuc_ptr,
#endif
	     Genomicpos_T chrstart, Genomicpos_T chrend,
	     Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp,
	     int genestrand, Oligoindex_T *oligoindices, int noligoindices,
	     Diagpool_T diagpool, bool debug_graphic_p, bool diagnosticp);

extern List_T
Stage2_compute (int *stage2_source, int *stage2_indexsize,
		char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,

#ifdef PMAP
		char *genomicuc_ptr,
#endif
		Genomicpos_T chrstart, Genomicpos_T chrend,
		Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp, int genestrand,
		Oligoindex_T *oligoindices, int noligoindices, double proceed_pctcoverage,
		Pairpool_T pairpool, Diagpool_T diagpool, int sufflookback, int nsufflookback,
		int maxintronlen, bool localp, bool skip_repetitive_p, bool use_shifted_canonical_p,
		bool favor_right_p, int max_nalignments, bool debug_graphic_p, bool diagnosticp,
		Stopwatch_T stopwatch, bool diag_debug);

extern List_T
Stage2_compute_one (int *stage2_source, int *stage2_indexsize,
		    char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,	
#ifdef PMAP
		    char *genomicuc_ptr,
#endif
		    Genomicpos_T chrstart, Genomicpos_T chrend,
		    Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp, int genestrand,

		    Oligoindex_T *oligoindices, int noligoindices, double proceed_pctcoverage,
		    Pairpool_T pairpool, Diagpool_T diagpool, int sufflookback, int nsufflookback,
		    int maxintronlen, bool localp, bool skip_repetitive_p, bool use_shifted_canonical_p,
		    bool favor_right_p, bool debug_graphic_p, bool diagnosticp,
		    Stopwatch_T stopwatch, bool diag_debug);
#undef T
#endif

