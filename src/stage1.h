/* $Id: stage1.h,v 1.73 2009/08/04 20:37:34 twu Exp $ */
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
#include "diagnostic.h"

#define T Stage1_T
typedef struct T *T;

extern void
Stage1_find_extensions (Genomicpos_T *extension5, Genomicpos_T *extension3, T this, 
			Match_T match5, Match_T match3, int maxintronlen_bound, bool maponlyp);

extern List_T
Stage1_compute (bool *samplingp, Sequence_T queryuc,
#ifdef PMAP
		Indexdb_T indexdb_fwd,
		Indexdb_T indexdb_rev,
#else
		Indexdb_T indexdb,
#endif
		int indexdb_size_threshold, IIT_T chromosome_iit, Chrsubset_T chrsubset,
		Matchpool_T matchpool, int maxintronlen_bound, int maxtotallen_bound,
		int stuttercycles, int stutterhits, bool subsequencep, Genome_T genome,
		Genomicpos_T genome_totallength, Diagnostic_T diagnostic, Stopwatch_T stopwatch);

#undef T
#endif

