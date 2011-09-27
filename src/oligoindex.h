/* $Id: oligoindex.h,v 1.44 2007/04/16 23:46:20 twu Exp $ */
#ifndef OLIGOINDEX_INCLUDED
#define OLIGOINDEX_INCLUDED
#include "bool.h"
#include "sequence.h"
#include "intlist.h"
#include "intpool.h"
#include "diagpool.h"

/* This number should be determined by the probability of two
   oligomers occurring within local_lookback */
#ifdef PMAP
#define SUFFNCONSECUTIVE 10
#else
#define SUFFNCONSECUTIVE 20
#endif

#define SUFF_PCTCOVERAGE 0.20

#define T Oligoindex_T
typedef struct T *T;

extern T
Oligoindex_new (int indexsize);
extern double
Oligoindex_set_inquery (int *badoligos, int *repoligos, int *trimoligos, int *trim_start, int *trim_end,
			T this, Sequence_T queryuc, bool trimp);
extern void
Oligoindex_tally (T this, Sequence_T genomicuc, Sequence_T queryuc, int maxoligohits);
extern void
Oligoindex_cleanup (T this, Sequence_T queryuc);
extern void
Oligoindex_free (T *old);

extern unsigned int **
Oligoindex_get_mappings (int *ndiagonals, int *maxnconsecutive, double *pct_coverage, double *pct_clear_coverage,
			 int **npositions, int *totalpositions,
			 unsigned int *minactive, unsigned int *maxactive,
			 T this, Sequence_T queryuc, Sequence_T genomicuc, Diagpool_T diagpool, int local_lookback,
			 bool debug_graphic_p);
extern void
Oligoindex_set_active (int *maxlength, int **active, unsigned int *minactive, unsigned int *maxactive,
		       Intlist_T *genomicdiag_querypos, Intlist_T *genomicdiag_hit, 
		       Diagpool_T diagpool, int genomiclength, int querylength, int local_lookback);
extern int *
Oligoindex_active_guide (int **active, int *npositions, unsigned int *minactive, unsigned int *maxactive,
			 unsigned int **mappings, int querylength, int minindexsize);

#undef T
#endif

