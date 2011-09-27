/* $Id: stage1.h,v 1.61 2006/11/28 00:50:01 twu Exp $ */
#ifndef STAGE1_INCLUDED
#define STAGE1_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "sequence.h"
#include "list.h"
#include "match.h"
#include "matchpool.h"
#include "iit-read.h"
#include "chrsubset.h"
#include "genome.h"
#include "stopwatch.h"

#define T Stage1_T
typedef struct T *T;

extern List_T
Stage1_matchlist (T this);

extern double
Stage1_runtime (T this);
extern int
Stage1_trial (T this);
extern int
Stage1_firstpair_searched_5 (T this);
extern int
Stage1_firstpair_searched_3 (T this);
extern int
Stage1_stutter_searched_5 (T this);
extern int
Stage1_stutter_searched_3 (T this);
extern int
Stage1_stutter_matchpairs (T this);
extern int
Stage1_stutter_matches_5 (T this);
extern int
Stage1_stutter_matches_3 (T this);
extern int
Stage1_dangling_5 (T this);
extern int
Stage1_dangling_3 (T this);
extern int
Stage1_dangling_matchpairs (T this);
extern int
Stage1_sampling_rounds (T this);
extern int
Stage1_sampling_nskip (T this);
extern int
Stage1_nentries_total (T this);
extern int
Stage1_nentries_used (T this);

extern void
Stage1_find_extensions (Genomicpos_T *extension5, Genomicpos_T *extension3, T this, 
			Match_T match5, Match_T match3, int maxintronlen_bound, bool maponlyp);

extern void
Stage1_free (T *old, bool completep);

extern T
Stage1_compute (Sequence_T queryuc, 
#ifdef PMAP
		Indexdb_T indexdb_fwd,
		Indexdb_T indexdb_rev,
#else
		Indexdb_T indexdb,
#endif
		IIT_T chromosome_iit, Chrsubset_T chrsubset, Matchpool_T matchpool,
		int maxintronlen_bound, int maxtotallen_bound, int stuttercycles, int stutterhits, bool crossspeciesp,
		bool subsequencep, Genome_T genome, bool diagnosticp, Stopwatch_T stopwatch);

#undef T
#endif

