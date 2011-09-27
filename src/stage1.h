/* $Id: stage1.h,v 1.56 2006/04/07 01:21:09 twu Exp $ */
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

#define T Stage1_T
typedef struct T *T;

extern List_T
Stage1_matchlist (T this);

extern void
Stage1_find_extensions (Genomicpos_T *extension5, Genomicpos_T *extension3, T this, 
			Match_T match5, Match_T match3, int maxextension, bool maponlyp);

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
		int maxtotallen_bound, int stuttercycles, int stutterhits, bool crossspeciesp,
		bool subsequencep, Genome_T genome);

#undef T
#endif

