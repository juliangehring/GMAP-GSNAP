/* $Id: oligoindex.h,v 1.31 2006/10/09 16:58:55 twu Exp $ */
#ifndef OLIGOINDEX_INCLUDED
#define OLIGOINDEX_INCLUDED
#include "bool.h"
#include "sequence.h"

#define T Oligoindex_T
typedef struct T *T;

extern T
Oligoindex_new (int indexsize);
extern double
Oligoindex_set_inquery (int *badoligos, int *repoligos, int *trimoligos, int *trim_start, int *trim_end,
			T this, Sequence_T queryuc, bool trimp);
extern void
Oligoindex_tally (T this, Sequence_T genomicuc, int maxoligohits);
extern void
Oligoindex_cleanup (T this, Sequence_T queryuc);
extern void
Oligoindex_free (T *old);

extern unsigned int **
Oligoindex_get_mappings (double *mapfraction, int *maxconsecutive, int **npositions, int *totalpositions, 
			 T this, Sequence_T queryuc);


#undef T
#endif

