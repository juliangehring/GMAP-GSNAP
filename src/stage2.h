/* $Id: stage2.h,v 1.39 2005/02/15 01:52:34 twu Exp $ */
#ifndef STAGE2_INCLUDED
#define STAGE2_INCLUDED
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
Stage2_compute (Sequence_T queryseq, Sequence_T genomicseg,
		Oligoindex_T oligoindex, int indexsize, Pairpool_T pairpool, 
		int sufflookback, int nsufflookback);

#undef T
#endif

