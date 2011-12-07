/* $Id: stage2.h 51826 2011-11-07 17:39:39Z twu $ */
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

extern void
Stage2_setup (bool splicingp_in, int suboptimal_score_start_in, int suboptimal_score_end_in);
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
	     char *genomicuc_ptr, int genomiclength,
	     Oligoindex_T *oligoindices, int noligoindices,
	     Diagpool_T diagpool, bool debug_graphic_p, bool diagnosticp);

extern List_T
Stage2_compute (int *stage2_source, int *stage2_indexsize,
		char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,

		char *genomicseg_ptr, char *genomicuc_ptr,
		Genomicpos_T genomicstart, Genomicpos_T genomicend,
		Genomicpos_T mappingstart, Genomicpos_T mappingend,
		bool plusp, int genomiclength, int genomic_offset,

		Oligoindex_T *oligoindices, int noligoindices, double proceed_pctcoverage,
		Pairpool_T pairpool, Diagpool_T diagpool, int sufflookback, int nsufflookback,
		int maxintronlen, bool localp, bool skip_repetitive_p, bool use_shifted_canonical_p,
		bool favor_right_p, bool just_one_p, bool debug_graphic_p, bool diagnosticp,
		Stopwatch_T stopwatch, bool diag_debug);

extern List_T
Stage2_compute_one (int *stage2_source, int *stage2_indexsize,
		    char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,	

		    char *genomicseg_ptr, char *genomicuc_ptr,
		    Genomicpos_T genomicstart, Genomicpos_T genomicend,
		    Genomicpos_T mappingstart, Genomicpos_T mappingend,
		    bool plusp, int genomiclength, int genomic_offset,

		    Oligoindex_T *oligoindices, int noligoindices, double proceed_pctcoverage,
		    Pairpool_T pairpool, Diagpool_T diagpool, int sufflookback, int nsufflookback,
		    int maxintronlen, bool localp, bool skip_repetitive_p, bool use_shifted_canonical_p,
		    bool favor_right_p, bool debug_graphic_p, bool diagnosticp,
		    Stopwatch_T stopwatch, bool diag_debug);
#undef T
#endif

