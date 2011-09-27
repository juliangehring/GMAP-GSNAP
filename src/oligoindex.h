/* $Id: oligoindex.h,v 1.59 2010/02/03 18:12:10 twu Exp $ */
#ifndef OLIGOINDEX_INCLUDED
#define OLIGOINDEX_INCLUDED
#include "bool.h"
#include "list.h"
#include "intlist.h"
#include "intpool.h"
#include "diagpool.h"

#define T Oligoindex_T
typedef struct T *T;

extern int
Oligoindex_indexsize (T this);

extern T *
Oligoindex_new_major (int *noligoindices);

extern T *
Oligoindex_new_minor (int *noligoindices);

extern int
Oligoindex_indexsize (T this);

extern double
Oligoindex_set_inquery (int *badoligos, int *repoligos, int *trimoligos, int *trim_start, int *trim_end,
			T this, char *queryuc_ptr, int querylength, bool trimp);
extern void
Oligoindex_tally (T this, char *genomicuc_trimptr, int genomicuc_trimlength,
		  char *queryuc_ptr, int querylength);
extern void
Oligoindex_untally (T this);
extern void
Oligoindex_clear_inquery (T this);
extern void
Oligoindex_free_array (T **oligoindices, int noligoindices);

extern List_T
Oligoindex_get_mappings (List_T diagonals, bool *coveredp, unsigned int **mappings, int *npositions,
			 int *totalpositions, bool *oned_matrix_p, int *maxnconsecutive, 
			 unsigned int *minactive, unsigned int *maxactive,
			 T this, char *queryuc_ptr, int querylength, int genomiclength, Diagpool_T diagpool);

#undef T
#endif

