/* $Id: stage1.h,v 1.49 2005/10/01 05:06:42 twu Exp $ */
#ifndef STAGE1_INCLUDED
#define STAGE1_INCLUDED
#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "sequence.h"
#include "list.h"
#include "match.h"
#include "iit-read.h"
#include "chrsubset.h"

#define T Stage1_T
typedef struct T *T;

extern List_T
Stage1_matchlist (T this,
#ifdef PMAP
		  Indexdb_T indexdb_fwd,
		  Indexdb_T indexdb_rev,
#else
		  Indexdb_T indexdb,
#endif
		  IIT_T chromosome_iit, Chrsubset_T chrsubset);

extern void
Stage1_find_extensions (Genomicpos_T *extension5, Genomicpos_T *extension3, T this, 
			Match_T match5, Match_T match3, int maxextension);

extern void
Stage1_free (T *old);

extern T
Stage1_compute (Sequence_T queryuc, 
#ifdef PMAP
		Indexdb_T indexdb_fwd,
		Indexdb_T indexdb_rev,
#else
		Indexdb_T indexdb,
#endif
		IIT_T chromosome_iit, Chrsubset_T chrsubset, int maxintronlen_bound, int stuttercycles, 
		int stutterhits, bool crossspeciesp);

#undef T
#endif

