/* $Id: stage2.h,v 1.44 2005/10/01 05:13:00 twu Exp $ */
#ifndef STAGE2_INCLUDED
#define STAGE2_INCLUDED
#include "bool.h"
#include "sequence.h"
#include "list.h"
#include "oligoindex.h"
#include "pairpool.h"

#define T Stage2_T
typedef struct T *T;

extern List_T
Stage2_path (T this);
extern double
Stage2_defect_rate (T this);
extern int
Stage2_nfwdintrons (T this);
extern int
Stage2_nrevintrons (T this);
extern int
Stage2_nunkintrons (T this);
extern int
Stage2_pathlength (T this);
	    
extern void
Stage2_free (T *old);

extern T
Stage2_compute (Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		Sequence_T queryntseq,
#endif
		Sequence_T genomicseg, Sequence_T genomicuc,
		Oligoindex_T oligoindex, int indexsize, Pairpool_T pairpool, 
		int sufflookback, int nsufflookback, int badoligos, bool crossspeciesp);

#undef T
#endif

