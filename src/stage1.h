/* $Id: stage1.h 82071 2012-12-19 21:43:45Z twu $ */
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

#ifdef PMAP
extern void
Stage1_setup (int index1part_aa_in, Genomicpos_T maxextension_in, Genomicpos_T maxtotallen_bound_in,
	      int min_extra_end_in, int circular_typeint_in);
#else
extern void
Stage1_setup (int index1part_in, Genomicpos_T maxextension_in, Genomicpos_T maxtotallen_bound_in,
	      int min_extra_end_in, int circular_typeint_in);
#endif


extern List_T
Stage1_compute (bool *samplingp, Sequence_T queryuc,
		Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev,
		int indexdb_size_threshold, IIT_T chromosome_iit, Chrsubset_T chrsubset,
		Matchpool_T matchpool, int stutterhits, Diagnostic_T diagnostic, Stopwatch_T stopwatch,
		int nbest);

extern List_T
Stage1_compute_nonstranded (bool *lowidentityp, Sequence_T queryuc,
			    Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev,
			    int indexdb_size_threshold, IIT_T chromosome_iit, Chrsubset_T chrsubset,
			    Matchpool_T matchpool, int stutterhits, Diagnostic_T diagnostic, Stopwatch_T stopwatch,
			    int nbest);
#undef T
#endif

