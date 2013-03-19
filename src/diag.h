/* $Id: diag.h 79302 2012-11-15 23:55:58Z twu $ */
#ifndef DIAG_INCLUDED
#define DIAG_INCLUDED
#include "bool.h"
#include "list.h"
#include "genomicpos.h"

#define T Diag_T
typedef struct T *T;

#ifndef USE_DIAGPOOL
extern T
Diag_new (Genomicpos_T diagonal, int querystart, int queryend, int nconsecutive);
extern void
Diag_free (T *old);
#endif

extern Genomicpos_T
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
Diag_compute_bounds (Genomicpos_T *minactive, Genomicpos_T *maxactive, List_T diagonals,
		     int querylength, bool debug_graphic_p,
		     Genomicpos_T chrstart, Genomicpos_T chrend,
		     Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp);

#undef T
#endif

