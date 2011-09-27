/* $Id: stage2.h,v 1.58 2007/02/12 17:36:18 twu Exp $ */
#ifndef STAGE2_INCLUDED
#define STAGE2_INCLUDED
#include "bool.h"
#include "sequence.h"
#include "list.h"
#include "oligoindex.h"
#include "pairpool.h"
#include "intpool.h"
#include "diagpool.h"
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
Stage2_defectrate (T this);
extern int
Stage2_pathlength (T this);
	    
extern void
Stage2_free (T *old);

extern T
Stage2_compute (Sequence_T queryseq, Sequence_T queryuc,
		Sequence_T genomicseg, Sequence_T genomicuc, Oligoindex_T oligoindex,
		Gbuffer_T gbuffer, int maxoligohits, int minindexsize, int maxindexsize, Pairpool_T pairpool, 
		Intpool_T intpool, Diagpool_T diagpool, int sufflookback, int nsufflookback, int maxintronlen,
		bool crossspeciesp, bool debug_graphic_p, Stopwatch_T stopwatch);

#undef T
#endif

