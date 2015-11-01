/* $Id: diag.h 106198 2013-08-28 23:07:34Z twu $ */
#ifndef DIAG_INCLUDED
#define DIAG_INCLUDED
#include "bool.h"
#include "list.h"
#include "genomicpos.h"
#include "types.h"

#define T Diag_T
typedef struct T *T;

#ifndef USE_DIAGPOOL
extern T
Diag_new (Chrpos_T diagonal, int querystart, int queryend, int nconsecutive);
extern void
Diag_free (T *old);
#endif

extern Chrpos_T
Diag_diagonal (T this);
extern int
Diag_querystart (T this);
extern int
Diag_queryend (T this);
extern int
Diag_nconsecutive (T this);
extern bool
Diag_dominatedp (T this);
extern void
Diag_set_dominatedp (T this);
extern int
Diag_compare_nconsecutive (const void *x, const void *y);
extern int
Diag_compare_diagonal (const void *x, const void *y);
extern double
Diag_update_coverage (bool *coveredp, int *ncovered, List_T diagonals, int querylength);
extern int
Diag_compare_querystart (const void *x, const void *y);
extern void
Diag_print_segments (List_T diagonals, char *queryseq_ptr, char *genomicseg_ptr);
extern void
Diag_range (int *start, int *end, List_T diagonals, int querylength);
extern int
Diag_compute_bounds (int *diag_querystart, int *diag_queryend,
		     Chrpos_T *minactive, Chrpos_T *maxactive, List_T diagonals,
		     int querylength, bool debug_graphic_p,
		     Chrpos_T chrstart, Chrpos_T chrend,
		     Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp);
extern void
Diag_max_bounds (Chrpos_T *minactive, Chrpos_T *maxactive,
		 int querylength, Chrpos_T chrstart, Chrpos_T chrend,
		 Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp);


#undef T
#endif

