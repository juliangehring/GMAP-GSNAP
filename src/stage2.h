/* $Id: stage2.h,v 1.51 2006/10/04 19:30:22 twu Exp $ */
#ifndef STAGE2_INCLUDED
#define STAGE2_INCLUDED
#include "bool.h"
#include "sequence.h"
#include "list.h"
#include "oligoindex.h"
#include "pairpool.h"
#include "gbuffer.h"
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
Stage2_runtime (T this);
extern int
Stage2_indexsize (T this);
extern double
Stage2_mapfraction (T this);
extern int
Stage2_maxconsecutive (T this);
extern double
Stage2_defectrate (T this);
extern int
Stage2_pathlength (T this);
	    
extern void
Stage2_free (T *old);

extern T
Stage2_compute (Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		Sequence_T queryntseq,
#endif
		Sequence_T genomicseg, Sequence_T genomicuc, Oligoindex_T *oligoindices,
		Gbuffer_T gbuffer, int maxoligohits, int minindexsize, int maxindexsize, Pairpool_T pairpool, 
		int sufflookback, int nsufflookback, int maxintronlen, bool crossspeciesp,
		Stopwatch_T stopwatch);

#undef T
#endif

